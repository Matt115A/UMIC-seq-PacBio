#!/usr/bin/env python3
"""
Simple consensus pipeline for PacBio data.
Generates consensus sequences directly from cluster files using abpoa.
No VCF generation - just consensus sequences.
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

def run_abpoa_consensus(cluster_file, output_dir, max_reads=20):
    """Run abpoa to generate consensus sequence from first N reads."""
    temp_file = None
    try:
        cluster_name = Path(cluster_file).stem
        consensus_file = os.path.join(output_dir, f"{cluster_name}_consensus.fasta")
        
        # Skip if already exists
        if os.path.exists(consensus_file):
            return consensus_file
        
        # Create temporary file with first N sequences
        temp_file = os.path.join(output_dir, f"{cluster_name}_temp_{max_reads}.fasta")
        
        # Extract first N sequences (fixed to include full sequences)
        with open(cluster_file, 'r') as infile, open(temp_file, 'w') as outfile:
            count = 0
            in_sequence = False
            for line in infile:
                if line.startswith('>'):
                    if count >= max_reads:
                        break
                    count += 1
                    in_sequence = True
                if in_sequence:
                    outfile.write(line)
        
        # Run abpoa consensus
        cmd = [
            "abpoa", 
            "-r", "0",    # output consensus in FASTA format
            "-a", "0",    # heaviest bundling path algorithm
            "-o", consensus_file,
            temp_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.returncode != 0:
            return f"abpoa failed for {cluster_name}: {result.stderr}"
        
        return consensus_file
        
    except subprocess.TimeoutExpired:
        return f"abpoa timeout for {cluster_name}"
    except Exception as e:
        return f"Error in abpoa for {cluster_name}: {str(e)}"
    finally:
        # Always clean up temp file
        if temp_file and os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except:
                pass

def process_cluster_simple(cluster_file, output_dir, max_reads=20):
    """Process a single cluster to generate consensus sequence."""
    cluster_name = Path(cluster_file).stem
    
    try:
        # Generate consensus using abpoa
        consensus_result = run_abpoa_consensus(cluster_file, output_dir, max_reads)
        
        if isinstance(consensus_result, str) and consensus_result.endswith('.fasta'):
            return f"SUCCESS: {cluster_name}"
        else:
            return f"FAILED: {cluster_name} - {consensus_result}"
            
    except Exception as e:
        return f"ERROR: {cluster_name} - {str(e)}"

def main():
    # Configuration
    cluster_files_dir = '/Volumes/Matt115A1TB_1/temp-dict-epgh11/UMIclusterfull_fast'
    output_dir = '/Volumes/Matt115A1TB_1/temp-dict-epgh11/simple_consensus_results'
    max_reads = 20  # Use first 20 reads for consensus
    
    # Use 4 threads for parallel processing
    max_workers = 4
    
    print(f"Starting simple consensus pipeline...")
    print(f"Cluster files directory: {cluster_files_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Max reads per consensus: {max_reads}")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of cluster files
    cluster_files = []
    for file in os.listdir(cluster_files_dir):
        if file.endswith('.fasta') and file.startswith('cluster_'):
            cluster_files.append(os.path.join(cluster_files_dir, file))
    
    total_clusters = len(cluster_files)
    print(f"Found {total_clusters:,} cluster files to process")
    
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
                remaining = total_clusters - processed
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
            progress_percent = (processed / total_clusters) * 100
            bar_length = 50
            filled_length = int(bar_length * processed // total_clusters)
            bar = '█' * filled_length + '-' * (bar_length - filled_length)
            
            # Print progress line
            progress_line = (f"\rProgress: |{bar}| {progress_percent:.1f}% "
                           f"({processed:,}/{total_clusters:,}) | "
                           f"Success: {success_count:,} | Failed: {failed_count:,} | "
                           f"Rate: {rate:.1f} clusters/s | ETA: {eta_str} | "
                           f"Elapsed: {elapsed_str}")
            
            print(progress_line, end='', flush=True)
    
    # Process clusters in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_cluster_simple, cluster_file, output_dir, max_reads) for cluster_file in cluster_files]
        
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
    
    # Final summary
    total_time = time.time() - start_time
    print(f"\nSimple consensus pipeline completed!")
    print(f"Total clusters processed: {total_clusters:,}")
    print(f"Successful: {success_count:,}")
    print(f"Failed: {failed_count:,}")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average rate: {total_clusters/total_time:.1f} clusters/second")
    print(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()
