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
import traceback
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# Force immediate output - write directly to stderr if needed
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(line_buffering=True)


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


def find_approx(query: str, target: str, max_mismatch: int) -> int:
    """Return start index of first occurrence of query in target allowing up to max_mismatch mismatches; -1 if not found."""
    q = query.upper()
    t = target.upper()
    qlen = len(q)
    tlen = len(t)
    if qlen == 0 or qlen > tlen:
        return -1
    # Fast exact pass first
    pos = t.find(q)
    if pos != -1:
        return pos
    # Sliding Hamming window
    for i in range(0, tlen - qlen + 1):
        window = t[i:i+qlen]
        mism = 0
        for a, b in zip(q, window):
            if a != b:
                mism += 1
                if mism > max_mismatch:
                    break
        if mism <= max_mismatch:
            return i
    return -1


def load_probe_sequences(probe_fasta: str) -> Tuple[str, str]:
    """Load probe FASTA (single record)."""
    with open(probe_fasta, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    seq_lines = [l for l in lines if not l.startswith('>')]
    probe = ''.join(seq_lines)
    return probe.upper(), reverse_complement(probe.upper())


def merge_trimmed_internal(s1: str, s2: str, left_ignore: int, right_ignore: int, min_overlap: int = 12) -> str:
    """Create merged internal fragment from paired-end reads.
    - Trim first 22 and last 24 bases from each
    - Reverse-complement R2 internal
    - Merge by maximal suffix/prefix overlap; fallback to concatenation
    """
    s1f = s1.upper()
    s2f = s2.upper()
    if len(s1f) <= left_ignore + right_ignore or len(s2f) <= left_ignore + right_ignore:
        return ''
    r1 = s1f[left_ignore:len(s1f) - right_ignore]
    r2 = s2f[left_ignore:len(s2f) - right_ignore]
    r2rc = reverse_complement(r2)
    # find maximal overlap of r1 suffix with r2rc prefix
    max_ol = 0
    max_possible = min(len(r1), len(r2rc))
    for ol in range(max_possible, min_overlap - 1, -1):
        if r1.endswith(r2rc[:ol]):
            max_ol = ol
            break
    if max_ol >= min_overlap:
        return r1 + r2rc[max_ol:]
    # also try the other orientation (if library orientation unexpected)
    # overlap of r2rc suffix with r1 prefix
    for ol in range(max_possible, min_overlap - 1, -1):
        if r2rc.endswith(r1[:ol]):
            max_ol = ol
            return r2rc + r1[max_ol:]
    # fallback concatenate
    return r1 + r2rc


def extract_umi_from_assembled(seq: str, umi_len: int, left_ignore: int, right_ignore: int) -> str:
    """Extract UMI by removing first 22 and last 24 bases and returning the internal UMI-length substring.
    For Illumina, UMI lives in the internal window. If trimmed internal is too short, return empty.
    """
    s = seq.upper()
    if len(s) < left_ignore + right_ignore + umi_len:
        return ''
    internal = s[left_ignore:len(s)-right_ignore]
    if len(internal) < umi_len:
        return ''
    # Return first umi_len bases from internal window
    return internal[:umi_len]


def load_consensus_fasta(consensus_dir: str) -> Dict[str, str]:
    cons: Dict[str, str] = {}
    files = [f for f in os.listdir(consensus_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    total = len(files)
    for idx, fn in enumerate(files):
        if idx % 1000 == 0 and idx > 0:
            print(f"  Loaded {idx}/{total} consensus files...", flush=True)
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
    total = len(consensus_map)
    for count, (name, seq) in enumerate(consensus_map.items()):
        if count > 0 and count % 1000 == 0:
            print(f"  Indexed {count}/{total} consensus sequences...", flush=True)
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


def read_vcf_haplotype(vcf_path: str) -> Tuple[str, List[Tuple[str, int, str, str]]]:
    """Return a semicolon-joined, position-sorted mutation list and list of variant tuples.
    Format per mutation: CHROM:POS:REF>ALT. Returns empty string if no variants.
    """
    muts = []
    variants = []
    if not os.path.exists(vcf_path):
        return '', variants
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[0] if parts[0] else parts[2]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4].split(',')[0]
            variants.append((chrom, pos, ref, alt))
            muts.append((pos, f"{chrom}:{pos}:{ref}>{alt}"))
    muts.sort(key=lambda x: x[0])
    return ';'.join(m for _, m in muts), variants


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


def collect_haplotype_rows(variants_dir: str, ref_seq: str) -> Tuple[List[Tuple[str, str, str]], Dict[str, int]]:
    """Return list of (consensus, mutations_str, aa_mutations_str) and mapping from consensus name to row index."""
    rows: List[Tuple[str, str, str]] = []
    index: Dict[str, int] = {}
    ref = Seq(ref_seq)
    
    for fn in os.listdir(variants_dir):
        if not fn.endswith('.vcf'):
            continue
        name = Path(fn).stem
        vpath = os.path.join(variants_dir, fn)
        muts_str, variants = read_vcf_haplotype(vpath)
        
        # Convert to amino acid mutations - prevariant neess fall within codon
        codons: Dict[int, List[Tuple[int, str, str]]] = {}
        for chrom, pos, ref_allele, alt_allele in variants:
            pos_0 = pos - 1  # 0-based
            codon_pos = pos_0 // 3
            if codon_pos not in codons:
                codons[codon_pos] = []
            codons[codon_pos].append((pos_0, ref_allele, alt_allele))
        
        aa_muts = []
        for codon_pos in sorted(codons.keys()):
            codon_start = codon_pos * 3
            if codon_start + 3 > len(ref):
                # Codon extends beyond reference, skip this codon
                continue
            
            wt_codon = str(ref[codon_start:codon_start + 3])
            mut_codon = list(wt_codon)
            
            # Apply all mutations in this codon
            valid = True
            for pos_0, ref_allele, alt_allele in codons[codon_pos]:
                if len(ref_allele) != 1 or len(alt_allele) != 1:
                    valid = False
                    break
                rel_pos = pos_0 - codon_start
                if rel_pos < 3 and pos_0 < len(ref) and ref[pos_0] == ref_allele:
                    mut_codon[rel_pos] = alt_allele
                else:
                    valid = False
                    break
            
            if not valid:
                # Fallback: record as codon number (amino acid position)
                wt_aa = str(Seq(wt_codon).translate())
                aa_muts.append(f"{wt_aa}{codon_pos + 1}X")
                continue
            
            mut_codon_str = ''.join(mut_codon)
            wt_aa = str(Seq(wt_codon).translate())
            mut_aa = str(Seq(mut_codon_str).translate())
            
            if wt_aa != mut_aa:
                # Non-synonymous - record as amino acid
                aa_muts.append(f"{wt_aa}{codon_pos + 1}{mut_aa}")
            # Skip synonymous mutations (silent changes)
        
        aa_muts_str = '+'.join(aa_muts) if aa_muts else 'WT'
        rows.append((name, muts_str, aa_muts_str))
        index[name] = len(rows) - 1
    return rows, index


def find_pool_folders(pools_dir: str) -> List[str]:
    folders = []
    for entry in os.scandir(pools_dir):
        if entry.is_dir():
            folders.append(entry.path)
    folders.sort()
    return folders


def find_r1_r2_pairs(folder: str) -> List[Tuple[str, str]]:
    """Find all R1/R2 pairs in a folder. Pairs are identified by matching basenames
    that differ only by R1/R2 designation.
    Returns list of (r1_path, r2_path) tuples.
    """
    r1_files = {}
    r2_files = {}
    
    for fn in os.listdir(folder):
        if not (fn.endswith('.fastq') or fn.endswith('.fastq.gz')):
            continue
        
        full_path = os.path.join(folder, fn)
        
        # Normalize filename to extract base key (remove R1/R2 and extension)
        # Handle patterns: sample_R1.fastq.gz, sample_R1_001.fastq.gz, sampleR1.fastq.gz
        key = None
        
        if '_R1' in fn:
            # Replace _R1 with placeholder, then remove extension
            key = fn.replace('_R1', '_R*')
            # Remove extension
            if key.endswith('.fastq.gz'):
                key = key[:-9]  # Remove .fastq.gz
            elif key.endswith('.fastq'):
                key = key[:-6]  # Remove .fastq
            r1_files[key] = full_path
        elif 'R1.' in fn:
            # Handle R1.fastq.gz pattern
            key = fn.replace('R1.', 'R*.')
            if key.endswith('.fastq.gz'):
                key = key[:-9]
            elif key.endswith('.fastq'):
                key = key[:-6]
            r1_files[key] = full_path
        
        if '_R2' in fn:
            key = fn.replace('_R2', '_R*')
            if key.endswith('.fastq.gz'):
                key = key[:-9]
            elif key.endswith('.fastq'):
                key = key[:-6]
            r2_files[key] = full_path
        elif 'R2.' in fn:
            key = fn.replace('R2.', 'R*.')
            if key.endswith('.fastq.gz'):
                key = key[:-9]
            elif key.endswith('.fastq'):
                key = key[:-6]
            r2_files[key] = full_path
    
    # Match pairs by key
    pairs = []
    for key in r1_files:
        if key in r2_files:
            pairs.append((r1_files[key], r2_files[key]))
    
    return pairs


def run_ngs_count(pools_dir: str, consensus_dir: str, variants_dir: str, probe_fasta: str,
                  umi_len: int, umi_loc: str, output_csv: str, reference_fasta: str,
                  left_ignore: int = 22, right_ignore: int = 24) -> bool:
    try:
        print("Starting...", flush=True)
        print("Loading probe sequences...", flush=True)
        probe_fwd, probe_rev = load_probe_sequences(probe_fasta)
        print(f"Probe length: {len(probe_fwd)} bp", flush=True)

        print("Loading consensus sequences...", flush=True)
        consensus = load_consensus_fasta(consensus_dir)
        print(f"Loaded {len(consensus)} consensus sequences", flush=True)

        print("Building UMI index (circular, both strands)...", flush=True)
        umi_index = build_umi_index(consensus, umi_len)
        print(f"Index contains {len(umi_index)} unique UMIs", flush=True)

        print("Collecting variant rows...", flush=True)
        var_rows, per_consensus = collect_all_variant_rows(variants_dir)
        print(f"Found {len(var_rows)} unique variant rows across {len(per_consensus)} consensus", flush=True)
        
        print("Loading reference sequence...", flush=True)
        ref_record = SeqIO.read(reference_fasta, "fasta")
        ref_seq = str(ref_record.seq)
        print(f"Reference length: {len(ref_seq)} bp", flush=True)

        print(f"Finding pool folders in: {pools_dir}", flush=True)
        pool_folders = find_pool_folders(pools_dir)
        pool_names = [Path(p).name for p in pool_folders]
        print(f"Found {len(pool_folders)} pools: {pool_names}", flush=True)
    except Exception as e:
        print(f"ERROR in initialization: {e}", flush=True)
        traceback.print_exc()
        return False

    counts = [[0 for _ in pool_folders] for _ in range(len(var_rows))]
    # Haplotype rows (consensus-level, preserve multi-mutations)
    hap_rows, hap_index = collect_haplotype_rows(variants_dir, ref_seq)
    hap_counts = [[0 for _ in pool_folders] for _ in range(len(hap_rows))]

    for j, pool in enumerate(pool_folders):
        pool_name = pool_names[j]
        print(f"\nProcessing pool {j+1}/{len(pool_folders)}: {pool_name}")
        
        # Find all R1/R2 pairs in this pool
        pairs = find_r1_r2_pairs(pool)
        if not pairs:
            print(f"  WARNING: No R1/R2 pairs found, skipping")
            continue
        
        print(f"  Found {len(pairs)} R1/R2 pair(s)")
        
        # Accumulate counts across all pairs in this pool
        pool_read_count = 0
        pool_umi_extracted = 0
        pool_umi_matched = 0
        
        for pair_idx, (r1, r2) in enumerate(pairs):
            pair_name = Path(r1).stem.replace('_R1', '').replace('R1', '')
            if len(pairs) > 1:
                print(f"  Processing pair {pair_idx+1}/{len(pairs)}: {Path(r1).name} / {Path(r2).name}")
            else:
                print(f"  R1: {Path(r1).name}")
                print(f"  R2: {Path(r2).name}")
            
            # Prefer PEAR merging for robust overlap handling
            # Use unique prefix for each pair to avoid conflicts
            prefix = os.path.join(pool, f"{pool_name}_pair{pair_idx+1}")
            assembled_fastq = prefix + '.assembled.fastq'
            
            # Check if PEAR output already exists
            if os.path.exists(assembled_fastq):
                print("    Using existing PEAR output...", flush=True)
            else:
                # Run PEAR if output doesn't exist
                try:
                    print("    Merging with PEAR...", flush=True)
                    cmd = [
                        'pear',
                        '-f', r1,
                        '-r', r2,
                        '-o', prefix,
                        '-q', '20',
                        '-t', '100'
                    ]
                    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    if not os.path.exists(assembled_fastq):
                        print("    WARNING: PEAR did not produce assembled reads; falling back to on-the-fly merge")
                        assembled_fastq = ''
                except Exception as e:
                    print(f"    WARNING: PEAR merge failed ({e}); falling back to on-the-fly merge")
                    assembled_fastq = ''

            read_count = 0
            umi_extracted = 0
            umi_matched = 0

            if assembled_fastq:
                print("    Scanning assembled reads for UMI matches...", flush=True)
                for seq in read_fastq_sequences(assembled_fastq):
                    read_count += 1
                    if read_count % 100000 == 0:
                        print(f"    Processed {read_count:,} merged reads, extracted {umi_extracted}, matched {umi_matched}", end='\r')
                    umi = extract_umi_from_assembled(seq, umi_len, left_ignore, right_ignore)
                    if not umi:
                        continue
                    umi_extracted += 1
                    cons_name = umi_index.get(umi)
                    if not cons_name:
                        continue
                    umi_matched += 1
                    idxs = per_consensus.get(cons_name, [])
                    for i in idxs:
                        counts[i][j] += 1
                    if cons_name in hap_index:
                        hap_counts[hap_index[cons_name]][j] += 1
            else:
                # Fallback: stream paired reads and merge on the fly
                print("    Merging reads on-the-fly...", flush=True)
                for s1, s2 in zip(read_fastq_sequences(r1), read_fastq_sequences(r2)):
                    read_count += 1
                    if read_count % 100000 == 0:
                        print(f"    Processed {read_count:,} reads, extracted {umi_extracted}, matched {umi_matched}", end='\r')
                    # Merge, trim and extract from merged
                    merged = merge_trimmed_internal(s1, s2, left_ignore, right_ignore)
                    if not merged:
                        continue
                    umi = extract_umi_from_assembled(merged, umi_len, left_ignore, right_ignore)
                    if not umi:
                        continue
                    umi_extracted += 1
                    cons_name = umi_index.get(umi)
                    if not cons_name:
                        continue
                    umi_matched += 1
                    idxs = per_consensus.get(cons_name, [])
                    for i in idxs:
                        counts[i][j] += 1
                    if cons_name in hap_index:
                        hap_counts[hap_index[cons_name]][j] += 1
            
            pool_read_count += read_count
            pool_umi_extracted += umi_extracted
            pool_umi_matched += umi_matched
            
            if len(pairs) > 1:
                print(f"\n    Pair {pair_idx+1} summary: {read_count:,} reads, {umi_extracted:,} UMIs extracted, {umi_matched:,} matched")

        print(f"\n  Pool {pool_name} total: {pool_read_count:,} reads, {pool_umi_extracted:,} UMIs extracted, {pool_umi_matched:,} matched to consensus")

    print(f"\nWriting output to: {output_csv}")
    with open(output_csv, 'w') as out:
        header = ['CHROM', 'POS', 'REF', 'ALT'] + pool_names
        out.write(','.join(header) + '\n')
        for i, key in enumerate(var_rows):
            chrom, pos, ref, alt = key
            row = [str(chrom), str(pos), ref, alt] + [str(counts[i][j]) for j in range(len(pool_folders))]
            out.write(','.join(row) + '\n')
    # Write haplotype counts
    hap_csv = str(Path(output_csv).parent / 'pool_haplotype_counts.csv')
    print(f"Writing haplotype counts to: {hap_csv}")
    with open(hap_csv, 'w') as out:
        header = ['CONSENSUS', 'MUTATIONS', 'AA_MUTATIONS'] + pool_names
        out.write(','.join(header) + '\n')
        for i, (cons_name, muts, aa_muts) in enumerate(hap_rows):
            row = [cons_name, muts, aa_muts] + [str(hap_counts[i][j]) for j in range(len(pool_folders))]
            out.write(','.join(row) + '\n')

    # Write counts merged on non-synonymous amino acid mutations
    merged_csv = str(Path(output_csv).parent / 'merged_on_nonsyn_counts.csv')
    print(f"Writing merged non-synonymous counts to: {merged_csv}")
    merged: Dict[str, Dict[str, object]] = {}
    for i, (cons_name, muts, aa_muts) in enumerate(hap_rows):
        entry = merged.setdefault(aa_muts, {
            'consensus': [],
            'mutations': set(),
            'counts': [0] * len(pool_folders)
        })
        entry['consensus'].append(cons_name)
        if muts:
            entry['mutations'].add(muts)
        counts_row = entry['counts']
        for j in range(len(pool_folders)):
            counts_row[j] += hap_counts[i][j]

    with open(merged_csv, 'w') as out:
        header = ['AA_MUTATIONS', 'CONSENSUS_IDS', 'NUC_MUTATIONS'] + pool_names
        out.write(','.join(header) + '\n')
        for aa_muts in sorted(merged.keys()):
            entry = merged[aa_muts]
            consensus_ids = ';'.join(sorted(entry['consensus']))
            nuc_muts = ';'.join(sorted(entry['mutations'])) if entry['mutations'] else ''
            row = [aa_muts, consensus_ids, nuc_muts] + [str(c) for c in entry['counts']]
            out.write(','.join(row) + '\n')
    
    print("Done!")
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
    ap.add_argument('--reference', required=True, help='Reference FASTA file for amino acid mapping')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    ok = run_ngs_count(a.pools_dir, a.consensus_dir, a.variants_dir, a.probe, a.umi_len, a.umi_loc, a.output, a.reference)
    sys.exit(0 if ok else 1)


