#!/usr/bin/env python3
"""
Fitness analysis module for merged non-synonymous counts.

Analyzes fitness from input/output pool comparisons, generates mutability plots,
epistasis analysis, and fitness distributions.
"""

import os
import sys
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


def parse_aa_mutations(aa_mut_str: str) -> List[Tuple[str, int, str]]:
    """Parse AA mutation string like 'S45F+Y76P' into list of (wt, pos, mut) tuples.
    Returns empty list for 'WT'.
    """
    if aa_mut_str == 'WT' or not aa_mut_str:
        return []
    mutations = []
    # Pattern: single letter, number, single letter (e.g., S45F)
    pattern = r'([A-Z])(\d+)([A-Z*])'
    for match in re.finditer(pattern, aa_mut_str):
        wt, pos, mut = match.groups()
        mutations.append((wt, int(pos), mut))
    return mutations


def get_hamming_distance(mut1: List[Tuple[str, int, str]], mut2: List[Tuple[str, int, str]]) -> int:
    """Calculate Hamming distance between two mutation sets (number of different positions)."""
    pos1 = {pos: mut for _, pos, mut in mut1}
    pos2 = {pos: mut for _, pos, mut in mut2}
    all_pos = set(pos1.keys()) | set(pos2.keys())
    distance = 0
    for pos in all_pos:
        if pos1.get(pos) != pos2.get(pos):
            distance += 1
    return distance


def is_stop_codon(mut: Tuple[str, int, str]) -> bool:
    """Check if mutation introduces a stop codon."""
    return mut[2] == '*'


def is_proline(mut: Tuple[str, int, str]) -> bool:
    """Check if mutation is to proline."""
    return mut[2] == 'P'


def run_fitness_analysis(
    input_csv: str,
    output_dir: str,
    input_pools: List[str],
    output_pools: List[str],
    min_input: int = 10,
    aa_filter: Optional[str] = None
) -> bool:
    """
    Run fitness analysis on merged counts CSV.
    
    Args:
        input_csv: Path to merged_on_nonsyn_counts.csv
        output_dir: Directory to save plots and results
        input_pools: List of pool names that are inputs
        output_pools: List of pool names that are outputs (paired with inputs)
        min_input: Minimum count threshold in input pools
        aa_filter: Optional amino acid to filter mutability plot (e.g., 'S')
    """
    try:
        print("=" * 60)
        print("FITNESS ANALYSIS")
        print("=" * 60)
        
        # Validate input/output pairing
        if len(input_pools) != len(output_pools):
            print(f"ERROR: Number of input pools ({len(input_pools)}) must match output pools ({len(output_pools)})")
            return False
        
        print(f"Input pools: {input_pools}")
        print(f"Output pools: {output_pools}")
        print(f"Min input threshold: {min_input}")
        
        # Load data
        print(f"\nLoading data from: {input_csv}")
        df = pd.read_csv(input_csv)
        print(f"Loaded {len(df)} rows")
        
        # Get all pool columns (exclude metadata columns)
        metadata_cols = ['AA_MUTATIONS', 'CONSENSUS_IDS', 'NUC_MUTATIONS']
        all_pool_cols = [col for col in df.columns if col not in metadata_cols]
        
        # Validate pool names exist
        for pool in input_pools + output_pools:
            if pool not in all_pool_cols:
                print(f"ERROR: Pool '{pool}' not found in CSV. Available pools: {all_pool_cols}")
                return False
        
        # Filter by min_input threshold
        print(f"\nFiltering by min_input threshold ({min_input})...")
        input_mask = df[input_pools].min(axis=1) >= min_input
        df_filtered = df[input_mask].copy()
        print(f"Kept {len(df_filtered)} rows after filtering (removed {len(df) - len(df_filtered)})")
        
        if len(df_filtered) == 0:
            print("ERROR: No rows remain after filtering!")
            return False
        
        # Calculate relative frequencies (column normalize)
        print("\nCalculating relative frequencies...")
        for pool in input_pools + output_pools:
            total = df_filtered[pool].sum()
            if total > 0:
                df_filtered[f'rel_{pool}'] = df_filtered[pool] / total
            else:
                df_filtered[f'rel_{pool}'] = 0.0
        
        # Calculate log fitness ratios for each pair
        print("\nCalculating fitness (log ratios)...")
        fitness_cols = []
        for i, (in_pool, out_pool) in enumerate(zip(input_pools, output_pools)):
            rel_in = df_filtered[f'rel_{in_pool}']
            rel_out = df_filtered[f'rel_{out_pool}']
            # Avoid log(0) by adding small epsilon
            epsilon = 1e-10
            ratio = (rel_out + epsilon) / (rel_in + epsilon)
            fitness_col = f'fitness_{i+1}'
            df_filtered[fitness_col] = np.log(ratio)
            fitness_cols.append(fitness_col)
            print(f"  {in_pool} -> {out_pool}: {fitness_col}")
        
        # Calculate average fitness across all pairs
        df_filtered['fitness_avg'] = df_filtered[fitness_cols].mean(axis=1)
        
        # Calculate bootstrap confidence intervals for fitness (if multiple replicates)
        if len(fitness_cols) >= 2:
            print("\nCalculating bootstrap confidence intervals...")
            bootstrap_results = calculate_bootstrap_ci(df_filtered, fitness_cols, n_bootstrap=1000, ci_level=0.95)
            df_filtered['fitness_ci_lower'] = bootstrap_results['ci_lower']
            df_filtered['fitness_ci_upper'] = bootstrap_results['ci_upper']
            df_filtered['fitness_std'] = bootstrap_results['std']
            print(f"  Calculated 95% CI for {len(df_filtered)} variants")
        
        # Parse mutations
        print("\nParsing mutations...")
        df_filtered['mutations_parsed'] = df_filtered['AA_MUTATIONS'].apply(parse_aa_mutations)
        df_filtered['n_mutations'] = df_filtered['mutations_parsed'].apply(len)
        df_filtered['hamming'] = df_filtered['n_mutations']
        
        # Identify mutation types
        df_filtered['has_stop'] = df_filtered['mutations_parsed'].apply(
            lambda muts: any(is_stop_codon(m) for m in muts)
        )
        df_filtered['has_proline'] = df_filtered['mutations_parsed'].apply(
            lambda muts: any(is_proline(m) for m in muts) and not any(is_stop_codon(m) for m in muts)
        )
        df_filtered['other'] = ~(df_filtered['has_stop'] | df_filtered['has_proline'])
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Save processed data
        output_csv = os.path.join(output_dir, 'fitness_analysis_results.csv')
        df_filtered.to_csv(output_csv, index=False)
        print(f"\nSaved processed data to: {output_csv}")
        
        # Generate plots for all mutants
        print("\nGenerating plots for all mutants...")
        
        # 1. Mutability plot
        plot_mutability(df_filtered, output_dir, aa_filter)
        
        # 2. Epistasis plot (single + single = double)
        plot_epistasis_single_double(df_filtered, output_dir)
        
        # 3. Fitness distributions
        plot_fitness_distributions(df_filtered, output_dir)
        
        # 4. Hamming distance distributions
        plot_hamming_distributions(df_filtered, output_dir)
        
        # 5. Reproducibility plot (if multiple replicate pairs)
        if len(fitness_cols) >= 2:
            plot_reproducibility(df_filtered, output_dir, fitness_cols)
        
        # 6. Substitution matrix heatmap
        plot_substitution_matrix(df_filtered, output_dir)
        
        print(f"\n{'='*60}")
        print("FITNESS ANALYSIS COMPLETE")
        print(f"{'='*60}")
        print(f"Results saved to: {output_dir}")
        
        return True
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def calculate_bootstrap_ci(df: pd.DataFrame, fitness_cols: List[str], n_bootstrap: int = 1000, ci_level: float = 0.95) -> Dict[str, np.ndarray]:
    """
    Calculate bootstrap confidence intervals for average fitness.
    
    Method: Replicate-level bootstrap
    - Resample replicate pairs (with replacement) n_bootstrap times
    - For each bootstrap iteration, recalculate average fitness
    - Use percentiles of bootstrap distribution as confidence intervals
    
    Args:
        df: DataFrame with fitness columns
        fitness_cols: List of fitness column names (one per replicate pair)
        n_bootstrap: Number of bootstrap iterations (default: 1000)
        ci_level: Confidence level (default: 0.95 for 95% CI)
    
    Returns:
        Dictionary with 'ci_lower', 'ci_upper', and 'std' arrays
    """
    n_variants = len(df)
    n_replicates = len(fitness_cols)
    
    # Arrays to store bootstrap results
    bootstrap_fitness = np.zeros((n_bootstrap, n_variants))
    
    # Perform bootstrap resampling
    for b in range(n_bootstrap):
        # Resample replicate indices with replacement
        # This simulates having different sets of replicates
        resampled_indices = np.random.choice(n_replicates, size=n_replicates, replace=True)
        resampled_cols = [fitness_cols[i] for i in resampled_indices]
        
        # Calculate average fitness from resampled replicates
        bootstrap_fitness[b, :] = df[resampled_cols].mean(axis=1).values
    
    # Calculate percentiles for confidence intervals
    alpha = 1 - ci_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100
    
    ci_lower = np.percentile(bootstrap_fitness, lower_percentile, axis=0)
    ci_upper = np.percentile(bootstrap_fitness, upper_percentile, axis=0)
    std = np.std(bootstrap_fitness, axis=0)
    
    return {
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'std': std
    }


def plot_mutability(df: pd.DataFrame, output_dir: str, aa_filter: Optional[str] = None):
    """Plot average fitness at each position for Hamming 1 mutants, relative to WT."""
    print("  Generating mutability plot...")
    
    # Get WT fitness
    wt_rows = df[df['AA_MUTATIONS'] == 'WT']
    if len(wt_rows) == 0:
        print("    WARNING: No WT found in dataset, cannot calculate relative fitness")
        return
    wt_fitness = wt_rows['fitness_avg'].iloc[0]
    
    # Filter to Hamming 1 (single mutants) and Hamming 2 (double mutants) if filtering by AA
    if aa_filter is not None:
        # Include singles and doubles (doubles will be included if they contain the filtered AA)
        hamming1 = df[df['hamming'] == 1].copy()
        hamming2 = df[df['hamming'] == 2].copy()
        combined = pd.concat([hamming1, hamming2], ignore_index=True)
    else:
        # Only singles if no filter
        combined = df[df['hamming'] == 1].copy()
    
    if len(combined) == 0:
        print("    WARNING: No mutants found, skipping mutability plot")
        return
    
    # Extract position and fitness for each mutation
    positions = []
    fitnesses = []
    wt_aas = []
    mut_aas = []
    
    for _, row in combined.iterrows():
        muts = row['mutations_parsed']
        # For singles, include directly
        if len(muts) == 1:
            wt, pos, mut = muts[0]
            # Filter by mutant amino acid (not wild-type)
            if aa_filter is None or mut == aa_filter:
                positions.append(pos)
                fitnesses.append(row['fitness_avg'])
                wt_aas.append(wt)
                mut_aas.append(mut)
        # For doubles, include if any mutation matches the filter
        elif len(muts) == 2 and aa_filter is not None:
            for wt, pos, mut in muts:
                if mut == aa_filter:
                    positions.append(pos)
                    fitnesses.append(row['fitness_avg'])
                    wt_aas.append(wt)
                    mut_aas.append(mut)
    
    if len(positions) == 0:
        print(f"    WARNING: No mutations found for filter '{aa_filter}'")
        return
    
    # Group by position and calculate average fitness
    pos_fitness = defaultdict(list)
    for pos, fit in zip(positions, fitnesses):
        pos_fitness[pos].append(fit)
    
    pos_avg = {pos: np.mean(fits) for pos, fits in pos_fitness.items()}
    
    # Calculate relative fitness (position average - WT fitness)
    sorted_pos = sorted(pos_avg.keys())
    avg_fits = [pos_avg[p] - wt_fitness for p in sorted_pos]
    
    # Plot as bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Color bars: positive = beneficial (green), negative = deleterious (red), neutral = gray
    colors = ['green' if f > 0 else 'red' if f < 0 else 'gray' for f in avg_fits]
    ax.bar(sorted_pos, avg_fits, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Amino Acid Position')
    ax.set_ylabel('Average Fitness Relative to WT (log ratio)')
    title = 'Mutability: Average Fitness by Position (Hamming 1, Relative to WT)'
    if aa_filter:
        title += f' - {aa_filter} only'
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'mutability_plot.png')
    if aa_filter:
        output_file = os.path.join(output_dir, f'mutability_plot_{aa_filter}.png')
    plt.savefig(output_file)
    plt.close()
    print(f"    Saved: {output_file}")


def plot_epistasis_single_double(df: pd.DataFrame, output_dir: str):
    """Plot epistasis: sum of single mutant fitnesses vs double mutant fitness.
    Excludes stop codons from analysis.
    """
    print("  Generating epistasis plot (single + single = double)...")
    
    # Get all single mutants, excluding stop codons
    singles = df[df['hamming'] == 1].copy()
    singles = singles[~singles['has_stop']].copy()  # Remove stop codons
    
    # Get all double mutants, excluding those with stop codons
    doubles = df[df['hamming'] == 2].copy()
    doubles = doubles[~doubles['has_stop']].copy()  # Remove stop codons
    
    if len(singles) == 0 or len(doubles) == 0:
        print("    WARNING: Need both single and double mutants for epistasis plot")
        return
    
    # Build index of single mutants by their full mutation string
    single_index = {}
    for _, row in singles.iterrows():
        aa_mut_str = row['AA_MUTATIONS']
        if aa_mut_str != 'WT':
            single_index[aa_mut_str] = row['fitness_avg']
    
    # For each double mutant, check if both constituent singles exist
    x_vals = []
    y_vals = []
    labels = []
    
    for _, row in doubles.iterrows():
        muts = row['mutations_parsed']
        if len(muts) == 2:
            (wt1, pos1, mut1), (wt2, pos2, mut2) = muts
            # Skip if either mutation is a stop codon
            if mut1 == '*' or mut2 == '*':
                continue
            # Create single mutation strings
            single1 = f"{wt1}{pos1}{mut1}"
            single2 = f"{wt2}{pos2}{mut2}"
            
            if single1 in single_index and single2 in single_index:
                sum_singles = single_index[single1] + single_index[single2]
                double_fit = row['fitness_avg']
                x_vals.append(sum_singles)
                y_vals.append(double_fit)
                labels.append(f"{single1}+{single2}")
    
    if len(x_vals) == 0:
        print("    WARNING: No double mutants with both singles present")
        return
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(x_vals, y_vals, alpha=0.6, s=50)
    
    # Add diagonal line (additivity)
    min_val = min(min(x_vals), min(y_vals))
    max_val = max(max(x_vals), max(y_vals))
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='Additivity')
    
    ax.set_xlabel('Sum of Single Mutant Fitnesses')
    ax.set_ylabel('Double Mutant Fitness')
    ax.set_title('Epistasis: Single + Single = Double Mutants')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Add correlation coefficient
    corr = np.corrcoef(x_vals, y_vals)[0, 1]
    ax.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'epistasis_plot.png')
    plt.savefig(output_file)
    plt.close()
    print(f"    Saved: {output_file}")


def plot_fitness_distributions(df: pd.DataFrame, output_dir: str):
    """Plot overlaid KDE plots for single mutants: stop codons, proline, others."""
    print("  Generating fitness distribution plot...")
    
    # Filter to single mutants only
    singles = df[df['hamming'] == 1].copy()
    
    if len(singles) == 0:
        print("    WARNING: No single mutants found")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot KDEs for each category
    if singles['has_stop'].any():
        stops = singles[singles['has_stop']]['fitness_avg']
        sns.kdeplot(stops, label='Stop Codons', ax=ax, linewidth=2)
    
    if singles['has_proline'].any():
        prolines = singles[singles['has_proline']]['fitness_avg']
        sns.kdeplot(prolines, label='Proline Mutations', ax=ax, linewidth=2)
    
    if singles['other'].any():
        others = singles[singles['other']]['fitness_avg']
        sns.kdeplot(others, label='Other Mutations', ax=ax, linewidth=2)
    
    ax.set_xlabel('Fitness (log ratio)')
    ax.set_ylabel('Density')
    ax.set_title('Fitness Distribution: Single Mutants by Type')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'fitness_distributions.png')
    plt.savefig(output_file)
    plt.close()
    print(f"    Saved: {output_file}")


def plot_hamming_distributions(df: pd.DataFrame, output_dir: str):
    """Plot overlaid KDE plots for fitness at different Hamming distances (1-5)."""
    print("  Generating Hamming distance distribution plot...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, 5))
    
    for hamming in range(1, 6):
        hamming_data = df[df['hamming'] == hamming]['fitness_avg']
        if len(hamming_data) > 0:
            sns.kdeplot(hamming_data, label=f'Hamming {hamming} (n={len(hamming_data)})', 
                       ax=ax, linewidth=2, color=colors[hamming-1])
    
    ax.set_xlabel('Fitness (log ratio)')
    ax.set_ylabel('Density')
    ax.set_title('Fitness Distribution by Hamming Distance')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'hamming_distributions.png')
    plt.savefig(output_file)
    plt.close()
    print(f"    Saved: {output_file}")


def plot_reproducibility(df: pd.DataFrame, output_dir: str, fitness_cols: List[str]):
    """Plot reproducibility heatmaps comparing fitness across replicate pairs.
    Creates a multi-panel plot with heatmaps showing fitness correlation between replicates.
    Uses blue->red color scheme.
    """
    print("  Generating reproducibility plot...")
    
    n_replicates = len(fitness_cols)
    if n_replicates < 2:
        print("    WARNING: Need at least 2 replicate pairs for reproducibility plot")
        return
    
    # Determine number of panels (all pairwise comparisons)
    # For 3 replicates: 3 panels (1vs2, 1vs3, 2vs3)
    # For 4 replicates: 6 panels (1vs2, 1vs3, 1vs4, 2vs3, 2vs4, 3vs4)
    panels = []
    for i in range(n_replicates):
        for j in range(i + 1, n_replicates):
            panels.append((i, j, fitness_cols[i], fitness_cols[j]))
    
    n_panels = len(panels)
    if n_panels == 0:
        return
    
    # Calculate grid dimensions
    if n_panels <= 3:
        n_rows, n_cols = 1, n_panels
    elif n_panels <= 6:
        n_rows, n_cols = 2, 3
    else:
        n_rows = int(np.ceil(np.sqrt(n_panels)))
        n_cols = int(np.ceil(n_panels / n_rows))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    if n_panels == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Get global fitness range for consistent color scale
    all_fitness = []
    for col in fitness_cols:
        all_fitness.extend(df[col].dropna().tolist())
    vmin = np.percentile(all_fitness, 1)
    vmax = np.percentile(all_fitness, 99)
    
    for idx, (i, j, col1, col2) in enumerate(panels):
        ax = axes[idx]
        
        # Get data for this pair
        x_data = df[col1].dropna()
        y_data = df[col2].dropna()
        
        # Align data (only use rows where both replicates have data)
        common_idx = x_data.index.intersection(y_data.index)
        x_vals = df.loc[common_idx, col1].values
        y_vals = df.loc[common_idx, col2].values
        
        if len(x_vals) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Replicate {i+1} vs {j+1}')
            continue
        
        # Create 2D histogram (heatmap)
        hist, xedges, yedges = np.histogram2d(x_vals, y_vals, bins=50)
        
        # Create meshgrid for pcolormesh
        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
        
        # Plot heatmap with blue->red colormap
        im = ax.pcolormesh(X, Y, hist.T, cmap='coolwarm', shading='auto', vmin=0, vmax=np.percentile(hist[hist > 0], 95))
        
        # Add scatter overlay colored by fitness (optional - can be commented out)
        # avg_fitness = (x_vals + y_vals) / 2
        # scatter = ax.scatter(x_vals, y_vals, c=avg_fitness, cmap='coolwarm', 
        #                     alpha=0.3, s=10, vmin=vmin, vmax=vmax, edgecolors='none')
        
        # Add diagonal line (perfect correlation)
        min_val = min(np.min(x_vals), np.min(y_vals))
        max_val = max(np.max(x_vals), np.max(y_vals))
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, linewidth=1, label='y=x')
        
        # Calculate correlation
        corr = np.corrcoef(x_vals, y_vals)[0, 1]
        
        ax.set_xlabel(f'Fitness (Replicate {i+1})')
        ax.set_ylabel(f'Fitness (Replicate {j+1})')
        ax.set_title(f'Replicate {i+1} vs {j+1}\n(r={corr:.3f})')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=8)
        
        # Add colorbar for heatmap
        plt.colorbar(im, ax=ax, label='Density')
    
    # Hide unused subplots
    for idx in range(n_panels, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'reproducibility_plot.png')
    plt.savefig(output_file)
    plt.close()
    print(f"    Saved: {output_file}")


def plot_substitution_matrix(df: pd.DataFrame, output_dir: str):
    """Plot substitution matrix heatmap showing average fitness for each amino acid mutation type.
    Amino acids are arranged by similarity (hydrophobic, polar, charged, etc.).
    """
    print("  Generating substitution matrix heatmap...")
    
    # Define amino acid order based on physicochemical properties
    # This groups similar amino acids together for better visualization
    aa_order = [
        # Hydrophobic aliphatic
        'A', 'V', 'L', 'I', 'M',
        # Aromatic
        'F', 'W', 'Y',
        # Polar uncharged
        'S', 'T', 'N', 'Q',
        # Positively charged
        'K', 'R', 'H',
        # Negatively charged
        'D', 'E',
        # Special
        'C', 'G', 'P',
        # Stop codon
        '*'
    ]
    
    # Collect all (wt, mut) pairs and their fitnesses
    substitution_data = defaultdict(list)
    
    for _, row in df.iterrows():
        muts = row['mutations_parsed']
        fitness = row['fitness_avg']
        
        for wt, pos, mut in muts:
            # Skip if either is not a standard amino acid or stop
            if wt not in aa_order or mut not in aa_order:
                continue
            substitution_data[(wt, mut)].append(fitness)
    
    if len(substitution_data) == 0:
        print("    WARNING: No valid amino acid substitutions found")
        return
    
    # Calculate average fitness for each substitution
    substitution_matrix = {}
    substitution_counts = {}
    for (wt, mut), fitnesses in substitution_data.items():
        substitution_matrix[(wt, mut)] = np.mean(fitnesses)
        substitution_counts[(wt, mut)] = len(fitnesses)
    
    # Create matrix (rows = wild-type, cols = mutant)
    matrix = np.full((len(aa_order), len(aa_order)), np.nan)
    count_matrix = np.zeros((len(aa_order), len(aa_order)), dtype=int)
    
    for i, wt in enumerate(aa_order):
        for j, mut in enumerate(aa_order):
            if (wt, mut) in substitution_matrix:
                matrix[i, j] = substitution_matrix[(wt, mut)]
                count_matrix[i, j] = substitution_counts[(wt, mut)]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Determine color scale (centered around 0, with symmetric range)
    vmax = np.nanpercentile(np.abs(matrix), 95)
    vmin = -vmax
    
    # Create heatmap
    im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', vmin=vmin, vmax=vmax, interpolation='nearest')
    
    # Add text annotations for substitutions with data
    for i, wt in enumerate(aa_order):
        for j, mut in enumerate(aa_order):
            if not np.isnan(matrix[i, j]):
                # Format fitness value
                fitness_val = matrix[i, j]
                count = count_matrix[i, j]
                # Show fitness value and count
                text = f'{fitness_val:.2f}\n(n={count})'
                ax.text(j, i, text, ha='center', va='center', 
                       fontsize=7, color='white' if abs(fitness_val) > vmax * 0.5 else 'black',
                       weight='bold' if count >= 10 else 'normal')
    
    # Set ticks and labels
    ax.set_xticks(range(len(aa_order)))
    ax.set_yticks(range(len(aa_order)))
    ax.set_xticklabels(aa_order, fontsize=10)
    ax.set_yticklabels(aa_order, fontsize=10)
    
    # Labels
    ax.set_xlabel('Mutant Amino Acid', fontsize=12, fontweight='bold')
    ax.set_ylabel('Wild-Type Amino Acid', fontsize=12, fontweight='bold')
    ax.set_title('Substitution Matrix: Average Fitness by Mutation Type\n' +
                '(Red = Beneficial, Blue = Deleterious, White = No Data)', 
                fontsize=13, fontweight='bold', pad=20)
    
    # Add grid
    ax.set_xticks(np.arange(len(aa_order)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(aa_order)) - 0.5, minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Average Fitness (log ratio)', rotation=270, labelpad=20, fontsize=11)
    
    # Highlight diagonal (no change mutations - should be rare/zero)
    for i in range(len(aa_order)):
        if not np.isnan(matrix[i, i]):
            # Draw a box around diagonal elements
            rect = plt.Rectangle((i-0.5, i-0.5), 1, 1, fill=False, 
                               edgecolor='yellow', linewidth=2, linestyle='--')
            ax.add_patch(rect)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'substitution_matrix.png')
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"    Saved: {output_file}")
    
    # Also save the matrix as CSV for further analysis
    matrix_df = pd.DataFrame(matrix, index=aa_order, columns=aa_order)
    matrix_csv = os.path.join(output_dir, 'substitution_matrix.csv')
    matrix_df.to_csv(matrix_csv)
    print(f"    Saved matrix data to: {matrix_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Fitness analysis on merged non-synonymous counts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single input/output pair
  python fitness_analysis.py \\
    --input merged_on_nonsyn_counts.csv \\
    --output_dir results/ \\
    --input_pools pool1 \\
    --output_pools pool2 \\
    --min_input 10

  # Multiple pairs
  python fitness_analysis.py \\
    --input merged_on_nonsyn_counts.csv \\
    --output_dir results/ \\
    --input_pools pool1 pool2 \\
    --output_pools pool3 pool4 \\
    --min_input 10 \\
    --aa_filter S
        """
    )
    
    parser.add_argument('--input', required=True,
                       help='Input CSV file (merged_on_nonsyn_counts.csv)')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for plots and results')
    parser.add_argument('--input_pools', required=True, nargs='+',
                       help='Input pool names (space-separated)')
    parser.add_argument('--output_pools', required=True, nargs='+',
                       help='Output pool names (space-separated, paired with inputs)')
    parser.add_argument('--min_input', type=int, default=10,
                       help='Minimum count threshold in input pools (default: 10)')
    parser.add_argument('--aa_filter', type=str, default=None,
                       help='Filter mutability plot to specific mutant amino acid (e.g., S for serine, P for proline, * for stop codons)')
    
    args = parser.parse_args()
    
    success = run_fitness_analysis(
        args.input,
        args.output_dir,
        args.input_pools,
        args.output_pools,
        args.min_input,
        args.aa_filter
    )
    
    sys.exit(0 if success else 1)

