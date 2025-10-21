# UMIC-seq-PacBio


## Performance Optimizations (2025)

## UMIC-seq-pacbio: Complete Pipeline

For PacBio data analysis, use the complete pipeline entry-point that handles the entire workflow from UMI extraction to final variant analysis:

```bash
python UMIC-seq-pacbio.py all \
  --input reads.fastq.gz \
  --probe probe.fasta \
  --reference reference.fasta \
  --output_dir /path/to/output
```

probe.fasta should contain an approximately 50 bp sequence
reference.fasta should contain the reference gene
output_dir should lead to where you want to place the output files

### Pipeline Steps

The `all` command runs the complete pipeline:

1. **UMI Extraction**: Extract UMIs from PacBio reads
2. **Clustering**: Cluster similar UMIs using ultra-fast hash-based algorithm
3. **Consensus Generation**: Generate consensus sequences using abpoa
4. **Variant Calling**: Call variants using minimap2 and bcftools
5. **Analysis**: Generate detailed CSV with mutation analysis

### Individual Commands

You can also run individual steps:

```bash
# Extract UMIs
python UMIC-seq-pacbio.py extract \
  --input reads.fastq.gz \
  --probe probe.fasta \
  --umi_len 52 \
  --output ExtractedUMIs.fasta \
  --umi_loc up \
  --min_probe_score 15

# Cluster UMIs
python UMIC-seq-pacbio.py cluster \
  --input_umi ExtractedUMIs.fasta \
  --input_reads reads.fastq.gz \
  --output_dir UMIclusterfull_fast \
  --aln_thresh 0.47 \
  --size_thresh 10

# Generate consensus sequences
python UMIC-seq-pacbio.py consensus \
  --input_dir UMIclusterfull_fast \
  --output_dir consensus_results \
  --max_reads 20 \
  --max_workers 4

# Call variants
python UMIC-seq-pacbio.py variants \
  --input_dir consensus_results \
  --reference reference.fasta \
  --output_dir variant_results \
  --combined_vcf combined_variants.vcf \
  --max_workers 4

# Analyze variants
python UMIC-seq-pacbio.py analyze \
  --input_vcf combined_variants.vcf \
  --reference reference.fasta \
  --output final_results.csv
```

### PacBio-Specific Parameters

- `--umi_loc up`: UMI is upstream of probe (default: up)
- `--min_probe_score 15`: Minimum probe alignment score (default: 15)
- `--aln_thresh 0.47`
- `--max_reads 20`: Use first 20 reads per cluster for consensus (faster processing)

### Output Files

The pipeline generates:
- `ExtractedUMIs.fasta`: Extracted UMI sequences
- `UMIclusterfull_fast/`: Cluster files (cluster_1.fasta, cluster_2.fasta, ...)
- `consensus_results/`: Consensus sequences per cluster
- `variant_results/`: Individual VCF files per cluster
- `combined_variants.vcf`: Combined variant calls
- `final_results.csv`: Detailed analysis with amino acid mutations, Hamming distance, stop codons, and indels