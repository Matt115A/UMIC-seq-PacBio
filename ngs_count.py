#!/usr/bin/env python3
"""
NGS pool counting via UMI-to-consensus exact matching (circular, both strands).

- Streams paired-end reads per pool (R1/R2), extracts UMIs using probe and umi_loc
- Builds an index of all possible UMI-length substrings from circular consensus
  sequences on both strands for fast lookup
- On perfect UMI match, increments counts for that consensus' VCF variant rows
- Outputs a wide CSV: rows = unique VCF entries (CHROM, POS, REF, ALT),
  columns = pools

Designed for low memory and speed: consensus index built once; reads streamed.
"""

import os
import sys
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Iterable


def open_maybe_gzip(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def read_fastq_sequences(path: str) -> Iterable[str]:
    """Stream sequences from FASTQ (text or gz)."""
    with open_maybe_gzip(path) as fh:
        i = 0
        seq = None
        for line in fh:
            i += 1
            m = i % 4
            if m == 2:
                seq = line.strip()
                yield seq


def reverse_complement(seq: str) -> str:
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]


def load_probe_sequences(probe_fasta: str) -> Tuple[str, str]:
    """Load probe FASTA (single record)."""
    with open(probe_fasta, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    seq_lines = [l for l in lines if not l.startswith('>')]
    probe = ''.join(seq_lines)
    return probe.upper(), reverse_complement(probe.upper())


def extract_umi_from_read(seq: str, probe_fwd: str, probe_rev: str, umi_len: int, umi_loc: str) -> str:
    """Extract UMI using exact probe match (fast). Returns '' if not found/invalid.
    - umi_loc 'up': UMI upstream of probe
    - umi_loc 'down': UMI downstream of probe
    """
    s = seq.upper()
    # Try forward probe
    pos = s.find(probe_fwd)
    if pos != -1:
        if umi_loc == 'up':
            start = pos - umi_len
            end = pos
            if start >= 0:
                return s[start:end]
        else:
            start = pos + len(probe_fwd)
            end = start + umi_len
            if end <= len(s):
                return s[start:end]
        return ''
    # Try reverse probe
    pos = s.find(probe_rev)
    if pos != -1:
        if umi_loc == 'up':
            # UMI is upstream relative to probe orientation; on reverse this is to the right
            start = pos + len(probe_rev)
            end = start + umi_len
            if end <= len(s):
                return reverse_complement(s[start:end])
        else:
            # downstream relative to probe orientation; on reverse this is to the left
            start = pos - umi_len
            end = pos
            if start >= 0:
                return reverse_complement(s[start:end])
    return ''


def load_consensus_fasta(consensus_dir: str) -> Dict[str, str]:
    cons: Dict[str, str] = {}
    for fn in os.listdir(consensus_dir):
        if not fn.endswith('.fasta') and not fn.endswith('.fa'):
            continue
        path = os.path.join(consensus_dir, fn)
        name = Path(fn).stem
        # read first record only
        with open(path, 'r') as f:
            seq_lines = []
            for line in f:
                if line.startswith('>'):
                    if seq_lines:
                        break
                    continue
                seq_lines.append(line.strip())
        cons[name] = ''.join(seq_lines).upper()
    return cons


def build_umi_index(consensus_map: Dict[str, str], umi_len: int) -> Dict[str, str]:
    """Build mapping from UMI (exact) to consensus name.
    If collisions occur, keep the first seen (deterministic by os.listdir order).
    Index is built from circular sequences on both strands.
    """
    idx: Dict[str, str] = {}
    for name, seq in consensus_map.items():
        if not seq:
            continue
        seq_circ = (seq + seq)
        seq_rc = reverse_complement(seq)
        seq_rc_circ = (seq_rc + seq_rc)
        L = len(seq)
        if L < umi_len:
            continue
        # forward
        for i in range(L):
            umi = seq_circ[i:i+umi_len]
            if len(umi) == umi_len and umi not in idx:
                idx[umi] = name
        # reverse
        for i in range(L):
            umi = seq_rc_circ[i:i+umi_len]
            if len(umi) == umi_len and umi not in idx:
                idx[umi] = name
    return idx


def read_vcf_variants(vcf_path: str) -> List[Tuple[str, int, str, str]]:
    rows: List[Tuple[str, int, str, str]] = []
    if not os.path.exists(vcf_path):
        return rows
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[2] if parts[0] == '' else parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4].split(',')[0]
            rows.append((chrom, pos, ref, alt))
    return rows


def collect_all_variant_rows(variants_dir: str) -> Tuple[List[Tuple[str,int,str,str]], Dict[str, List[int]]]:
    """Return list of unique variant keys and mapping from consensus name to indices in that list."""
    unique: Dict[Tuple[str,int,str,str], int] = {}
    per_consensus: Dict[str, List[int]] = {}
    for fn in os.listdir(variants_dir):
        if not fn.endswith('.vcf'):
            continue
        name = Path(fn).stem  # matches consensus stem
        vpath = os.path.join(variants_dir, fn)
        idxs: List[int] = []
        for key in read_vcf_variants(vpath):
            if key not in unique:
                unique[key] = len(unique)
            idxs.append(unique[key])
        per_consensus[name] = idxs
    # invert to list
    rows: List[Tuple[str,int,str,str]] = [None] * len(unique)
    for key, i in unique.items():
        rows[i] = key
    return rows, per_consensus


def find_pool_folders(pools_dir: str) -> List[str]:
    folders = []
    for entry in os.scandir(pools_dir):
        if entry.is_dir():
            folders.append(entry.path)
    folders.sort()
    return folders


def find_r1_r2(folder: str) -> Tuple[str, str]:
    r1 = r2 = ''
    for fn in os.listdir(folder):
        if fn.endswith('.fastq') or fn.endswith('.fastq.gz'):
            if '_R1' in fn or 'R1.' in fn:
                r1 = os.path.join(folder, fn)
            elif '_R2' in fn or 'R2.' in fn:
                r2 = os.path.join(folder, fn)
    return r1, r2


def run_ngs_count(pools_dir: str, consensus_dir: str, variants_dir: str, probe_fasta: str,
                  umi_len: int, umi_loc: str, output_csv: str) -> bool:
    # Load probe
    probe_fwd, probe_rev = load_probe_sequences(probe_fasta)

    # Load consensus and build UMI index
    consensus = load_consensus_fasta(consensus_dir)
    umi_index = build_umi_index(consensus, umi_len)

    # Collect variant rows and per-consensus mapping
    var_rows, per_consensus = collect_all_variant_rows(variants_dir)

    # Prepare counts table: rows x pools
    pool_folders = find_pool_folders(pools_dir)
    pool_names = [Path(p).name for p in pool_folders]
    counts = [[0 for _ in pool_folders] for _ in range(len(var_rows))]

    for j, pool in enumerate(pool_folders):
        r1, r2 = find_r1_r2(pool)
        if not r1 or not r2:
            continue
        # Stream paired reads in lockstep
        for s1, s2 in zip(read_fastq_sequences(r1), read_fastq_sequences(r2)):
            # Prefer R1 for UMI extraction; fallback to R2
            umi = extract_umi_from_read(s1, probe_fwd, probe_rev, umi_len, umi_loc)
            if not umi:
                umi = extract_umi_from_read(s2, probe_fwd, probe_rev, umi_len, umi_loc)
            if not umi:
                continue
            cons_name = umi_index.get(umi)
            if not cons_name:
                continue
            idxs = per_consensus.get(cons_name, [])
            for i in idxs:
                counts[i][j] += 1

    # Write CSV
    with open(output_csv, 'w') as out:
        header = ['CHROM', 'POS', 'REF', 'ALT'] + pool_names
        out.write(','.join(header) + '\n')
        for i, key in enumerate(var_rows):
            chrom, pos, ref, alt = key
            row = [str(chrom), str(pos), ref, alt] + [str(counts[i][j]) for j in range(len(pool_folders))]
            out.write(','.join(row) + '\n')

    return True


if __name__ == '__main__':
    # Minimal CLI for ad-hoc execution
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--pools_dir', required=True)
    ap.add_argument('--consensus_dir', required=True)
    ap.add_argument('--variants_dir', required=True)
    ap.add_argument('--probe', required=True)
    ap.add_argument('--umi_len', type=int, required=True)
    ap.add_argument('--umi_loc', choices=['up','down'], required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    ok = run_ngs_count(a.pools_dir, a.consensus_dir, a.variants_dir, a.probe, a.umi_len, a.umi_loc, a.output)
    sys.exit(0 if ok else 1)


