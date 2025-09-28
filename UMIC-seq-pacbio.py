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
        "--umi_loc", "down",  # UMI is downstream of probe
        "--umi_len", str(args.umi_len),
        "--min_probe_score", "15"  # Use the working threshold
    ]
    
    return run_command(cmd, "UMI Extraction")

def run_clustering(args):
    """Run clustering step."""
    cmd = [
        "python", "UMIC-seq.py", "clusterfull",
        "-i", args.input_umi,
        "-o", args.output_dir,
        "--reads", args.input_reads,
        "--aln_thresh", str(int(args.aln_thresh * 100)),  # Convert to integer percentage
        "--size_thresh", str(args.size_thresh)
    ]
    
    return run_command(cmd, "UMI Clustering")

def run_consensus_generation(args):
    """Run consensus generation step."""
    # Create a temporary script with the correct paths
    script_content = f'''#!/usr/bin/env python3
import os
import sys
sys.path.append('.')

# Override the configuration in simple_consensus_pipeline.py
import simple_consensus_pipeline

# Update the configuration
simple_consensus_pipeline.main = lambda: None

# Set the configuration
cluster_files_dir = '{args.input_dir}'
output_dir = '{args.output_dir}'
max_reads = {args.max_reads}
max_workers = {args.max_workers}

# Run the consensus pipeline
from simple_consensus_pipeline import run_abpoa_consensus, main as original_main
import tempfile
import shutil
from concurrent.futures import ThreadPoolExecutor
import time
import threading

def main():
    print(f"Starting simple consensus pipeline...")
    print(f"Cluster files directory: {{cluster_files_dir}}")
    print(f"Output directory: {{output_dir}}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open('clusterfiles.txt', 'w') as f:
        for cluster_file in os.listdir(cluster_files_dir):
            if cluster_file.endswith('.fasta'):
                f.write(os.path.join(cluster_files_dir, cluster_file) + '\\n')
    
    # Run the original main function
    original_main()

if __name__ == "__main__":
    main()
'''
    
    with open('temp_consensus_script.py', 'w') as f:
        f.write(script_content)
    
    cmd = ["python", "temp_consensus_script.py"]
    success = run_command(cmd, "Consensus Generation")
    
    # Clean up
    if os.path.exists('temp_consensus_script.py'):
        os.remove('temp_consensus_script.py')
    if os.path.exists('clusterfiles.txt'):
        os.remove('clusterfiles.txt')
    
    return success

def run_variant_calling(args):
    """Run variant calling step."""
    # Create a temporary script with the correct paths
    script_content = f'''#!/usr/bin/env python3
import os
import sys
sys.path.append('.')

# Override the configuration in sensitive_variant_pipeline.py
import sensitive_variant_pipeline

# Update the configuration
consensus_dir = '{args.input_dir}'
reference_file = '{args.reference}'
output_dir = '{args.output_dir}'
combined_vcf = '{args.combined_vcf}'
max_workers = {args.max_workers}

# Run the variant calling pipeline
from sensitive_variant_pipeline import main as original_main

def main():
    print(f"Starting SENSITIVE variant calling pipeline...")
    print(f"This will call ALL variants, including single mismatches!")
    print(f"Consensus directory: {{consensus_dir}}")
    print(f"Reference file: {{reference_file}}")
    print(f"Output directory: {{output_dir}}")
    print(f"Combined VCF: {{combined_vcf}}")
    
    # Run the original main function
    original_main()

if __name__ == "__main__":
    main()
'''
    
    with open('temp_variant_script.py', 'w') as f:
        f.write(script_content)
    
    cmd = ["python", "temp_variant_script.py"]
    success = run_command(cmd, "Variant Calling")
    
    # Clean up
    if os.path.exists('temp_variant_script.py'):
        os.remove('temp_variant_script.py')
    
    return success

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
        size_thresh=args.size_thresh,
        output_dir=cluster_dir
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
    all_parser.add_argument('--aln_thresh', type=float, default=0.47, help='Alignment threshold (default: 0.47)')
    all_parser.add_argument('--size_thresh', type=int, default=10, help='Size threshold (default: 10)')
    all_parser.add_argument('--max_reads', type=int, default=20, help='Max reads per consensus (default: 20)')
    all_parser.add_argument('--max_workers', type=int, default=4, help='Max parallel workers (default: 4)')
    
    # Individual step commands
    extract_parser = subparsers.add_parser('extract', help='Extract UMIs from raw reads')
    extract_parser.add_argument('--input', required=True, help='Input FASTQ file (can be .gz)')
    extract_parser.add_argument('--probe', required=True, help='Probe FASTA file')
    extract_parser.add_argument('--umi_len', type=int, default=52, help='UMI length (default: 52)')
    extract_parser.add_argument('--output', required=True, help='Output FASTA file')
    
    cluster_parser = subparsers.add_parser('cluster', help='Cluster UMIs')
    cluster_parser.add_argument('--input_umi', required=True, help='Input UMI FASTA file')
    cluster_parser.add_argument('--input_reads', required=True, help='Input reads FASTQ file (can be .gz)')
    cluster_parser.add_argument('--aln_thresh', type=float, default=0.47, help='Alignment threshold (default: 0.47)')
    cluster_parser.add_argument('--size_thresh', type=int, default=10, help='Size threshold (default: 10)')
    cluster_parser.add_argument('--output_dir', required=True, help='Output directory for clusters')
    
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
    else:
        print(f"❌ Unknown command: {args.command}")
        return 1
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
