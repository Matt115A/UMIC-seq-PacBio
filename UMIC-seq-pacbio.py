#!/usr/bin/env python3
"""
UMIC-seq PacBio Pipeline - Main Entry Point
Complete pipeline for processing PacBio data from raw FASTQ to detailed mutation analysis.

This script orchestrates the entire UMIC-seq pipeline:
1. UMI extraction from raw PacBio reads
2. Clustering of similar UMIs
3. Consensus generation using abpoa
4. Variant calling with sensitive parameters
5. Detailed mutation analysis and CSV output

Usage:
    python UMIC-seq-pacbio.py --help
    python UMIC-seq-pacbio.py all --input raw_reads.fastq.gz --probe probe.fasta --reference reference.fasta --output_dir /path/to/output
    python UMIC-seq-pacbio.py extract --input raw_reads.fastq.gz --probe probe.fasta --output umis.fasta
    python UMIC-seq-pacbio.py cluster --input_umi umis.fasta --input_reads raw_reads.fastq.gz --output_dir clusters/
    python UMIC-seq-pacbio.py consensus --input_dir clusters/ --output_dir consensus/
    python UMIC-seq-pacbio.py variants --input_dir consensus/ --reference reference.fasta --output_dir variants/
    python UMIC-seq-pacbio.py analyze --input_vcf combined.vcf --reference reference.fasta --output detailed.csv
"""

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path

def run_command(cmd, description="", check=True):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    print(f"COMMAND: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    start_time = time.time()
    try:
        result = subprocess.run(cmd, check=check, capture_output=False, text=True)
        elapsed = time.time() - start_time
        print(f"\n✓ COMPLETED: {description} ({elapsed:.1f}s)")
        return True
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        print(f"\n✗ FAILED: {description} ({elapsed:.1f}s)")
        print(f"Exit code: {e.returncode}")
        if e.stdout:
            print(f"STDOUT: {e.stdout}")
        if e.stderr:
            print(f"STDERR: {e.stderr}")
        return False

def run_umi_extraction(args):
    """Run UMI extraction step."""
    cmd = [
        "python", "UMIC-seq.py", "UMIextract",
        "-i", args.input,
        "-o", args.output,
        "--probe", args.probe,
        "--umi_loc", args.umi_loc,  # Use the actual umi_loc parameter
        "--umi_len", str(args.umi_len),
        "--min_probe_score", "15"  # Use the working threshold
    ]
    
    return run_command(cmd, "UMI Extraction")

def run_clustering(args):
    """Run clustering step."""
    # Determine which clustering method to use
    use_fast = args.fast and not args.slow
    
    if use_fast:
        # Use fast CD-HIT clustering
        cmd = [
            "python", "UMIC-seq.py", "clusterfull_fast",
            "-i", args.input_umi,
            "-o", args.output_dir,
            "--reads", args.input_reads,
            "--identity", str(args.identity),
            "--size_thresh", str(args.size_thresh)
        ]
        description = "Fast UMI Clustering (CD-HIT)"
    else:
        # Use slow alignment-based clustering
        cmd = [
            "python", "UMIC-seq.py", "clusterfull",
            "-i", args.input_umi,
            "-o", args.output_dir,
            "--reads", args.input_reads,
            "--aln_thresh", str(int(args.aln_thresh * 100)),  # Convert to integer percentage
            "--size_thresh", str(args.size_thresh)
        ]
        description = "Slow UMI Clustering (alignment-based)"
    
    return run_command(cmd, description)

def run_consensus_generation(args):
    """Run consensus generation step."""
    # Import the consensus functions directly
    import simple_consensus_pipeline
    
    cluster_files_dir = args.input_dir
    output_dir = args.output_dir
    max_reads = args.max_reads
    max_workers = args.max_workers
    
    print(f"Starting simple consensus pipeline...")
    print(f"Cluster files directory: {cluster_files_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Max reads per consensus: {max_reads}")
    print(f"Max workers: {max_workers}")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of cluster files
    cluster_files = []
    for file in os.listdir(cluster_files_dir):
        if file.endswith('.fasta') and file.startswith('cluster_'):
            cluster_files.append(os.path.join(cluster_files_dir, file))
    
    total_clusters = len(cluster_files)
    print(f"Found {total_clusters:,} cluster files to process")
    
    if total_clusters == 0:
        print("❌ No cluster files found!")
        return False
    
    # Track progress
    success_count = 0
    failed_count = 0
    start_time = time.time()
    
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import gc
    
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
                
                if eta_minutes >= 1:
                    eta_str = f"{eta_minutes:.1f}m"
                else:
                    eta_str = f"{eta_seconds:.0f}s"
            else:
                eta_str = "unknown"
            
            # Progress percentage
            progress_percent = (processed / total_clusters) * 100
            
            # Print progress
            progress_line = (f"\rProgress: {progress_percent:.1f}% "
                           f"({processed:,}/{total_clusters:,}) | "
                           f"Success: {success_count:,} | Failed: {failed_count:,} | "
                           f"Rate: {rate:.1f} clusters/s | ETA: {eta_str}")
            
            print(progress_line, end='', flush=True)
            
            # Periodic garbage collection every 5000 clusters
            if processed > 0 and processed % 5000 == 0:
                gc.collect()
                print(f"\n[GC at {processed:,} clusters]", flush=True)
    
    # Process clusters in parallel with batched submission
    BATCH_SIZE = 1000  # Submit 1000 at a time to avoid memory buildup
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for batch_start in range(0, total_clusters, BATCH_SIZE):
            batch_end = min(batch_start + BATCH_SIZE, total_clusters)
            batch_files = cluster_files[batch_start:batch_end]
            
            # Submit batch
            futures = {executor.submit(simple_consensus_pipeline.process_cluster_simple, cf, output_dir, max_reads): cf 
                      for cf in batch_files}
            
            # Process as they complete
            for future in as_completed(futures):
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
    print(f"\n\nConsensus generation completed!")
    print(f"Total clusters processed: {total_clusters:,}")
    print(f"Successful: {success_count:,}")
    print(f"Failed: {failed_count:,}")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Output directory: {output_dir}")
    
    return success_count > 0

def run_variant_calling(args):
    """Run variant calling step."""
    import sensitive_variant_pipeline
    import threading
    import time
    from pathlib import Path
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    consensus_dir = args.input_dir
    reference_file = args.reference
    output_dir = args.output_dir
    combined_vcf = args.combined_vcf
    max_workers = args.max_workers
    
    print(f"\n{'='*60}")
    print(f"RUNNING: Variant Calling")
    print(f"{'='*60}")
    
    start_time = time.time()
    
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
    print(f"Found {total_consensus:,} consensus files to process\n")
    
    # Track progress
    success_count = 0
    failed_count = 0
    
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
    
    # Process consensus files in parallel with batched submission
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(sensitive_variant_pipeline.process_consensus_file_sensitive, 
                                  cf, reference_file, output_dir): cf 
                  for cf in consensus_files}
        
        # Process as they complete
        for future in as_completed(futures):
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
    combine_success = sensitive_variant_pipeline.combine_vcf_files(output_dir, combined_vcf)
    
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
    
    print(f"\n✓ COMPLETED: Variant Calling ({total_time:.1f}s)")
    
    return success_count > 0

def run_analysis(args):
    """Run detailed analysis step."""
    cmd = [
        "python", "vcf2csv_detailed.py",
        "--input", args.input_vcf,
        "--reference", args.reference,
        "--output", args.output
    ]
    
    return run_command(cmd, "Detailed Analysis")

def run_full_pipeline(args):
    """Run the complete pipeline."""
    print(f"\n{'='*80}")
    print("STARTING COMPLETE UMIC-seq PACBIO PIPELINE")
    print(f"{'='*80}")
    print(f"Input FASTQ: {args.input}")
    print(f"Probe file: {args.probe}")
    print(f"Reference file: {args.reference}")
    print(f"Output directory: {args.output_dir}")
    print(f"UMI length: {args.umi_len}")
    print(f"Alignment threshold: {args.aln_thresh}")
    print(f"Size threshold: {args.size_thresh}")
    print(f"Max reads per consensus: {args.max_reads}")
    print(f"Max workers: {args.max_workers}")
    print(f"{'='*80}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: UMI Extraction
    umi_file = os.path.join(args.output_dir, "ExtractedUMIs.fasta")
    extract_args = argparse.Namespace(
        input=args.input,
        probe=args.probe,
        umi_len=args.umi_len,
        umi_loc=args.umi_loc,
        output=umi_file
    )
    
    if not run_umi_extraction(extract_args):
        print("❌ Pipeline failed at UMI extraction step")
        return False
    
    # Step 2: Clustering
    cluster_dir = os.path.join(args.output_dir, "clusters")
    cluster_args = argparse.Namespace(
        input_umi=umi_file,
        input_reads=args.input,
        aln_thresh=args.aln_thresh,
        identity=getattr(args, 'identity', 0.90),
        size_thresh=args.size_thresh,
        output_dir=cluster_dir,
        fast=getattr(args, 'fast', True),
        slow=getattr(args, 'slow', False)
    )
    
    if not run_clustering(cluster_args):
        print("❌ Pipeline failed at clustering step")
        return False
    
    # Step 3: Consensus Generation
    consensus_dir = os.path.join(args.output_dir, "consensus")
    consensus_args = argparse.Namespace(
        input_dir=cluster_dir,
        output_dir=consensus_dir,
        max_reads=args.max_reads,
        max_workers=args.max_workers
    )
    
    if not run_consensus_generation(consensus_args):
        print("❌ Pipeline failed at consensus generation step")
        return False
    
    # Step 4: Variant Calling
    variant_dir = os.path.join(args.output_dir, "variants")
    combined_vcf = os.path.join(args.output_dir, "combined_variants.vcf")
    variant_args = argparse.Namespace(
        input_dir=consensus_dir,
        reference=args.reference,
        output_dir=variant_dir,
        combined_vcf=combined_vcf,
        max_workers=args.max_workers
    )
    
    if not run_variant_calling(variant_args):
        print("❌ Pipeline failed at variant calling step")
        return False
    
    # Step 5: Detailed Analysis
    analysis_file = os.path.join(args.output_dir, "detailed_mutations.csv")
    analysis_args = argparse.Namespace(
        input_vcf=combined_vcf,
        reference=args.reference,
        output=analysis_file
    )
    
    if not run_analysis(analysis_args):
        print("❌ Pipeline failed at analysis step")
        return False
    
    print(f"\n{'='*80}")
    print("✅ PIPELINE COMPLETED SUCCESSFULLY!")
    print(f"{'='*80}")
    print(f"Final output: {analysis_file}")
    print(f"Cluster directory: {cluster_dir}")
    print(f"Consensus directory: {consensus_dir}")
    print(f"Variant directory: {variant_dir}")
    print(f"Combined VCF: {combined_vcf}")
    print(f"{'='*80}")
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description="UMIC-seq PacBio Pipeline - Complete pipeline for processing PacBio data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python UMIC-seq-pacbio.py all --input raw_reads.fastq.gz --probe probe.fasta --reference reference.fasta --output_dir /path/to/output
  
  # Run individual steps
  python UMIC-seq-pacbio.py extract --input raw_reads.fastq.gz --probe probe.fasta --output umis.fasta
  python UMIC-seq-pacbio.py cluster --input_umi umis.fasta --input_reads raw_reads.fastq.gz --output_dir clusters/
  python UMIC-seq-pacbio.py consensus --input_dir clusters/ --output_dir consensus/
  python UMIC-seq-pacbio.py variants --input_dir consensus/ --reference reference.fasta --output_dir variants/
  python UMIC-seq-pacbio.py analyze --input_vcf combined.vcf --reference reference.fasta --output detailed.csv
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline step to run')
    
    # Full pipeline command
    all_parser = subparsers.add_parser('all', help='Run complete pipeline')
    all_parser.add_argument('--input', required=True, help='Input FASTQ file (can be .gz)')
    all_parser.add_argument('--probe', required=True, help='Probe FASTA file')
    all_parser.add_argument('--reference', required=True, help='Reference FASTA file')
    all_parser.add_argument('--output_dir', required=True, help='Output directory')
    all_parser.add_argument('--umi_len', type=int, default=52, help='UMI length (default: 52)')
    all_parser.add_argument('--aln_thresh', type=float, default=0.47, help='Alignment threshold for slow clustering (default: 0.47)')
    all_parser.add_argument('--identity', type=float, default=0.90, help='Sequence identity for fast clustering (default: 0.90)')
    all_parser.add_argument('--size_thresh', type=int, default=10, help='Size threshold (default: 10)')
    all_parser.add_argument('--max_reads', type=int, default=20, help='Max reads per consensus (default: 20)')
    all_parser.add_argument('--max_workers', type=int, default=4, help='Max parallel workers (default: 4)')
    all_parser.add_argument('--fast', action='store_true', default=True, help='Use fast CD-HIT clustering (default: True)')
    all_parser.add_argument('--slow', action='store_true', help='Use slow alignment-based clustering')
    all_parser.add_argument('--umi_loc', type=str, default='up', choices=['up', 'down'], help='UMI location relative to probe (up or down, default: up)')
    all_parser.add_argument('--min_probe_score', type=int, default=15, help='Minimal alignment score of probe for processing (default: 15)')
    
    # Individual step commands
    extract_parser = subparsers.add_parser('extract', help='Extract UMIs from raw reads')
    extract_parser.add_argument('--input', required=True, help='Input FASTQ file (can be .gz)')
    extract_parser.add_argument('--probe', required=True, help='Probe FASTA file')
    extract_parser.add_argument('--umi_len', type=int, default=52, help='UMI length (default: 52)')
    extract_parser.add_argument('--output', required=True, help='Output FASTA file')
    
    cluster_parser = subparsers.add_parser('cluster', help='Cluster UMIs')
    cluster_parser.add_argument('--input_umi', required=True, help='Input UMI FASTA file')
    cluster_parser.add_argument('--input_reads', required=True, help='Input reads FASTQ file (can be .gz)')
    cluster_parser.add_argument('--aln_thresh', type=float, default=0.47, help='Alignment threshold for slow method (default: 0.47)')
    cluster_parser.add_argument('--identity', type=float, default=0.90, help='Sequence identity for fast method (default: 0.90)')
    cluster_parser.add_argument('--size_thresh', type=int, default=10, help='Size threshold (default: 10)')
    cluster_parser.add_argument('--output_dir', required=True, help='Output directory for clusters')
    cluster_parser.add_argument('--fast', action='store_true', default=True, help='Use fast CD-HIT clustering (default: True)')
    cluster_parser.add_argument('--slow', action='store_true', help='Use slow alignment-based clustering')
    
    consensus_parser = subparsers.add_parser('consensus', help='Generate consensus sequences')
    consensus_parser.add_argument('--input_dir', required=True, help='Input cluster directory')
    consensus_parser.add_argument('--output_dir', required=True, help='Output consensus directory')
    consensus_parser.add_argument('--max_reads', type=int, default=20, help='Max reads per consensus (default: 20)')
    consensus_parser.add_argument('--max_workers', type=int, default=4, help='Max parallel workers (default: 4)')
    
    variant_parser = subparsers.add_parser('variants', help='Call variants')
    variant_parser.add_argument('--input_dir', required=True, help='Input consensus directory')
    variant_parser.add_argument('--reference', required=True, help='Reference FASTA file')
    variant_parser.add_argument('--output_dir', required=True, help='Output variant directory')
    variant_parser.add_argument('--combined_vcf', required=True, help='Combined VCF output file')
    variant_parser.add_argument('--max_workers', type=int, default=4, help='Max parallel workers (default: 4)')
    
    analyze_parser = subparsers.add_parser('analyze', help='Analyze variants and generate detailed CSV')
    analyze_parser.add_argument('--input_vcf', required=True, help='Input combined VCF file')
    analyze_parser.add_argument('--reference', required=True, help='Reference FASTA file')
    analyze_parser.add_argument('--output', required=True, help='Output CSV file')
    
    # NGS counting command
    ngs_parser = subparsers.add_parser('ngs_count', help='Count pool reads per variant via UMI matching')
    ngs_parser.add_argument('--pools_dir', required=True, help='Directory containing per-pool folders with R1/R2 fastqs')
    ngs_parser.add_argument('--consensus_dir', required=True, help='Consensus directory (from consensus step)')
    ngs_parser.add_argument('--variants_dir', required=True, help='Variants directory with per-consensus VCFs')
    ngs_parser.add_argument('--probe', required=True, help='Probe FASTA file (same used for UMI extraction)')
    ngs_parser.add_argument('--umi_len', type=int, default=52, help='UMI length (default: 52)')
    ngs_parser.add_argument('--umi_loc', type=str, default='up', choices=['up','down'], help='UMI location relative to probe (default: up)')
    ngs_parser.add_argument('--output', required=True, help='Output counts CSV file')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Check that required scripts exist
    required_scripts = ['UMIC-seq.py', 'simple_consensus_pipeline.py', 'sensitive_variant_pipeline.py', 'vcf2csv_detailed.py']
    missing_scripts = [script for script in required_scripts if not os.path.exists(script)]
    
    if missing_scripts:
        print(f"❌ Missing required scripts: {', '.join(missing_scripts)}")
        return 1
    
    # Run the appropriate command
    if args.command == 'all':
        success = run_full_pipeline(args)
    elif args.command == 'extract':
        success = run_umi_extraction(args)
    elif args.command == 'cluster':
        success = run_clustering(args)
    elif args.command == 'consensus':
        success = run_consensus_generation(args)
    elif args.command == 'variants':
        success = run_variant_calling(args)
    elif args.command == 'analyze':
        success = run_analysis(args)
    elif args.command == 'ngs_count':
        print("=" * 60)
        print("NGS POOL COUNTING")
        print("=" * 60)
        from ngs_count import run_ngs_count
        success = run_ngs_count(args.pools_dir, args.consensus_dir, args.variants_dir, args.probe, args.umi_len, args.umi_loc, args.output)
    else:
        print(f"❌ Unknown command: {args.command}")
        return 1
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
