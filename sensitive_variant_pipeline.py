#!/usr/bin/env python3
"""
Sensitive variant calling pipeline that calls ALL variants, including single mismatches.
Uses custom SAM parsing to bypass bcftools filtering.
"""

import os
import subprocess
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import shutil
import re

def extract_variants_from_sam(sam_file, reference_file, output_vcf):
    """Extract all variants from SAM alignment and create VCF."""
    
    # Read reference sequence
    with open(reference_file, 'r') as f:
        ref_lines = [line.strip() for line in f if not line.startswith('>')]
        ref_seq = ''.join(ref_lines)
    
    # Read SAM alignment
    # Use dictionary to deduplicate variants by position (chr:pos:ref:alt)
    # This ensures only one variant per position is kept
    variants_dict = {}
    
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            
            # Check SAM flag - skip secondary alignments (flag 0x100 = 256)
            # Only process primary alignments
            try:
                flag = int(parts[1])
                if flag & 0x100:  # Secondary alignment flag
                    continue
            except (ValueError, IndexError):
                continue
                
            cigar = parts[5]
            seq = parts[9]
            
            # Skip if CIGAR is * (unmapped)
            if cigar == '*':
                continue
            
            # Parse CIGAR string to find mismatches
            pos = int(parts[3]) - 1  # Convert to 0-based
            seq_pos = 0
            ref_pos = pos
            
            # Simple CIGAR parsing (only handles M, D, I, S)
            cigar_parts = re.findall(r'(\d+)([MDIS])', cigar)
            
            for length, op in cigar_parts:
                length = int(length)
                
                if op == 'M':  # Match or mismatch
                    for i in range(length):
                        if ref_pos < len(ref_seq) and seq_pos < len(seq):
                            if ref_seq[ref_pos] != seq[seq_pos]:
                                # Found a mismatch!
                                # Use position as key to deduplicate
                                # If same position has different ALT, keep the first one encountered
                                var_key = (parts[2], ref_pos + 1)  # chr, pos (1-based)
                                if var_key not in variants_dict:
                                    variants_dict[var_key] = {
                                        'chr': parts[2],
                                        'pos': ref_pos + 1,  # Convert back to 1-based
                                        'ref': ref_seq[ref_pos],
                                        'alt': seq[seq_pos],
                                        'qual': 30  # Default quality
                                    }
                            ref_pos += 1
                            seq_pos += 1
                
                elif op == 'D':  # Deletion in read (insertion in reference)
                    # Skip reference positions
                    ref_pos += length
                
                elif op == 'I':  # Insertion in read (deletion in reference)
                    # Skip read positions
                    seq_pos += length
                
                elif op == 'S':  # Soft clip
                    # Skip read positions
                    seq_pos += length
    
    # Convert dictionary to sorted list
    variants = sorted(variants_dict.values(), key=lambda x: (x['chr'], x['pos']))
    
    # Write VCF
    with open(output_vcf, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##reference={reference_file}\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        
        # Write variants
        for var in variants:
            f.write(f"{var['chr']}\t{var['pos']}\t.\t{var['ref']}\t{var['alt']}\t"
                   f"{var['qual']}\tPASS\tDP=1\tGT:AD:DP\t1/1:0,1:1\n")
    
    return len(variants)

def align_and_call_variants_sensitive(consensus_file, reference_file, output_dir):
    """Align consensus to reference and call ALL variants using sensitive method."""
    try:
        consensus_name = Path(consensus_file).stem
        vcf_file = os.path.join(output_dir, f"{consensus_name}.vcf")
        
        # Skip if VCF already exists
        if os.path.exists(vcf_file):
            return vcf_file
        
        # Step 1: Align consensus to reference using minimap2
        sam_file = os.path.join(output_dir, f"{consensus_name}.sam")
        
        # Run minimap2 alignment
        minimap2_cmd = [
            "minimap2",
            "-ax", "map-ont",
            "-t", "1",
            reference_file,
            consensus_file
        ]
        
        with open(sam_file, 'w') as sam_out:
            subprocess.run(minimap2_cmd, check=True, stdout=sam_out, stderr=subprocess.PIPE, text=True)
        
        # Step 2: Extract variants directly from SAM
        variant_count = extract_variants_from_sam(sam_file, reference_file, vcf_file)
        
        # Clean up intermediate files
        os.remove(sam_file)
        
        if variant_count == 0:
            return f"No variants found for {consensus_name}"
        
        return vcf_file
        
    except Exception as e:
        return f"Error processing {consensus_name}: {str(e)}"

def process_consensus_file_sensitive(consensus_file, reference_file, output_dir):
    """Process a single consensus file to generate VCF using sensitive method."""
    consensus_name = Path(consensus_file).stem
    
    try:
        # Align and call variants
        vcf_result = align_and_call_variants_sensitive(consensus_file, reference_file, output_dir)
        
        if isinstance(vcf_result, str) and vcf_result.endswith('.vcf'):
            # Count variants in the VCF
            with open(vcf_result, 'r') as f:
                variant_count = sum(1 for line in f if not line.startswith('#'))
            return f"SUCCESS: {consensus_name} ({variant_count} variants)"
        else:
            return f"FAILED: {consensus_name} - {vcf_result}"
            
    except Exception as e:
        return f"ERROR: {consensus_name} - {str(e)}"

def combine_vcf_files(vcf_dir, output_vcf):
    """Combine all VCF files into a single VCF with one entry per consensus."""
    try:
        # Get list of VCF files
        vcf_files = []
        for file in os.listdir(vcf_dir):
            if file.endswith('.vcf'):
                vcf_files.append(os.path.join(vcf_dir, file))
        
        if not vcf_files:
            print("No VCF files found to combine")
            return False
        
        print(f"Combining {len(vcf_files)} VCF files...")
        
        # Create combined VCF header
        with open(output_vcf, 'w') as outfile:
            # Write VCF header
            outfile.write("##fileformat=VCFv4.2\n")
            outfile.write("##INFO=<ID=CLUSTER,Number=1,Type=String,Description=\"Cluster ID\">\n")
            outfile.write("##INFO=<ID=VARIANT_COUNT,Number=1,Type=Integer,Description=\"Number of variants in this consensus\">\n")
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Process each VCF file
            for vcf_file in vcf_files:
                cluster_name = Path(vcf_file).stem
                
                try:
                    # Read VCF file and extract variants
                    with open(vcf_file, 'r') as infile:
                        variants = []
                        for line in infile:
                            if line.startswith('#') or line.strip() == '':
                                continue
                            
                            parts = line.strip().split('\t')
                            if len(parts) >= 8:
                                # Add cluster information to INFO field
                                info_parts = parts[7].split(';')
                                info_parts.append(f"CLUSTER={cluster_name}")
                                parts[7] = ';'.join(info_parts)
                                variants.append('\t'.join(parts))
                        
                        # Write variants to combined VCF
                        for variant in variants:
                            outfile.write(variant + '\n')
                            
                except Exception as e:
                    print(f"Error processing {vcf_file}: {e}")
                    continue
        
        print(f"Combined VCF written to: {output_vcf}")
        return True
        
    except Exception as e:
        print(f"Error combining VCF files: {e}")
        return False

def main():
    # Configuration
    consensus_dir = '/Volumes/Matt115A1TB_1/temp-dict-epgh11/simple_consensus_results'
    reference_file = 'GH11_epPCR/reference.fasta'
    output_dir = '/Volumes/Matt115A1TB_1/temp-dict-epgh11/sensitive_variant_results'
    combined_vcf = '/Volumes/Matt115A1TB_1/temp-dict-epgh11/combined_variants_sensitive.vcf'
    
    # Use 4 threads for parallel processing
    max_workers = 4
    
    print(f"Starting SENSITIVE variant calling pipeline...")
    print(f"This will call ALL variants, including single mismatches!")
    print(f"Consensus directory: {consensus_dir}")
    print(f"Reference file: {reference_file}")
    print(f"Output directory: {output_dir}")
    print(f"Combined VCF: {combined_vcf}")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of consensus files
    consensus_files = []
    for file in os.listdir(consensus_dir):
        if file.endswith('_consensus.fasta'):
            consensus_files.append(os.path.join(consensus_dir, file))
    
    total_consensus = len(consensus_files)
    print(f"Found {total_consensus:,} consensus files to process")
    
    # Track progress
    success_count = 0
    failed_count = 0
    start_time = time.time()
    
    # Thread-safe progress tracking
    progress_lock = threading.Lock()
    
    def update_progress():
        nonlocal success_count, failed_count
        with progress_lock:
            processed = success_count + failed_count
            elapsed_time = time.time() - start_time
            rate = processed / elapsed_time if elapsed_time > 0 else 0
            
            # Calculate ETA
            if rate > 0:
                remaining = total_consensus - processed
                eta_seconds = remaining / rate
                eta_minutes = eta_seconds / 60
                eta_hours = eta_minutes / 60
                
                if eta_hours >= 1:
                    eta_str = f"{eta_hours:.1f}h"
                elif eta_minutes >= 1:
                    eta_str = f"{eta_minutes:.1f}m"
                else:
                    eta_str = f"{eta_seconds:.0f}s"
            else:
                eta_str = "unknown"
            
            # Format elapsed time
            if elapsed_time >= 3600:
                elapsed_str = f"{elapsed_time/3600:.1f}h"
            elif elapsed_time >= 60:
                elapsed_str = f"{elapsed_time/60:.1f}m"
            else:
                elapsed_str = f"{elapsed_time:.0f}s"
            
            # Progress bar
            progress_percent = (processed / total_consensus) * 100
            bar_length = 50
            filled_length = int(bar_length * processed // total_consensus)
            bar = '█' * filled_length + '-' * (bar_length - filled_length)
            
            # Print progress line
            progress_line = (f"\rProgress: |{bar}| {progress_percent:.1f}% "
                           f"({processed:,}/{total_consensus:,}) | "
                           f"Success: {success_count:,} | Failed: {failed_count:,} | "
                           f"Rate: {rate:.1f} consensus/s | ETA: {eta_str} | "
                           f"Elapsed: {elapsed_str}")
            
            print(progress_line, end='', flush=True)
    
    # Process consensus files in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_consensus_file_sensitive, consensus_file, reference_file, output_dir) for consensus_file in consensus_files]
        
        # Process completed jobs
        for future in futures:
            result = future.result()
            
            with progress_lock:
                if "SUCCESS" in result:
                    success_count += 1
                else:
                    failed_count += 1
                    if failed_count <= 10:  # Only show first 10 failures
                        print(f"\nFailed: {result}")
            
            update_progress()
    
    # Combine VCF files
    print(f"\n\nCombining VCF files...")
    combine_success = combine_vcf_files(output_dir, combined_vcf)
    
    # Final summary
    total_time = time.time() - start_time
    print(f"\nSensitive variant calling pipeline completed!")
    print(f"Total consensus processed: {total_consensus:,}")
    print(f"Successful: {success_count:,}")
    print(f"Failed: {failed_count:,}")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average rate: {total_consensus/total_time:.1f} consensus/second")
    print(f"Individual VCF files: {output_dir}")
    if combine_success:
        print(f"Combined VCF: {combined_vcf}")
    else:
        print("Failed to create combined VCF")

if __name__ == "__main__":
    main()
