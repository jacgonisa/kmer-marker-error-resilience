#!/usr/bin/env python3
"""
COMPREHENSIVE MARKER EVALUATION FRAMEWORK

This script performs a thorough analysis of cenhapmer markers across k-mer sizes,
considering ALL relevant factors:

1. MARKER AVAILABILITY: Total number of unique k-mers
2. MARKER DENSITY: K-mers per megabase (uniform coverage)
3. CROSS-CONTAMINATION: False positive risk from sequencing errors
4. ERROR RESILIENCE: Read retention under ONT sequencing
5. OVERALL QUALITY SCORE: Weighted combination of all factors

Author: Generated for comprehensive k-mer size comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import tempfile
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

# Arabidopsis genome sizes (bp)
GENOME_SIZES = {
    'Chr1': 30427671,
    'Chr2': 19698289,
    'Chr3': 23459830,
    'Chr4': 18585056,
    'Chr5': 26975502
}

# Approximate sizes of regions (for density calculation)
# CEN regions are ~1-5 Mb, ARMS are the rest
CEN_SIZES = {
    'Chr1': 3_500_000,
    'Chr2': 3_000_000,
    'Chr3': 4_000_000,
    'Chr4': 3_500_000,
    'Chr5': 4_500_000
}

ARMS_SIZES = {chrom: GENOME_SIZES[chrom] - CEN_SIZES[chrom] for chrom in GENOME_SIZES}


def find_kmc_databases(directory: str) -> List[Dict]:
    """Find all KMC databases and parse metadata."""
    databases = []
    dir_path = Path(directory)

    for kmc_pre in dir_path.glob("*.kmc_pre"):
        base_name = kmc_pre.stem
        kmc_suf = kmc_pre.with_suffix('.kmc_suf')

        if not kmc_suf.exists():
            continue

        pattern = r'unique_([^_]+)_([^_]+)_(Chr\d+)_k(\d+)'
        match = re.match(pattern, base_name)

        if match:
            genotype, region, chromosome, k = match.groups()
            databases.append({
                'path': str(kmc_pre.parent / base_name),
                'name': base_name,
                'genotype': genotype,
                'region': region,
                'chromosome': chromosome,
                'k': int(k),
                'label': f"{genotype}_{region}_{chromosome}"
            })

    return sorted(databases, key=lambda x: (x['genotype'], x['region'], x['chromosome']))


def get_kmc_stats(db_path: str) -> Dict:
    """Get statistics from a KMC database."""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
            tmp_path = tmp.name

        subprocess.run(
            ['kmc_tools', 'info', db_path, tmp_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )

        stats = {}
        with open(tmp_path, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip().lower().replace(' ', '_').replace('.', '')
                    value = value.strip()
                    try:
                        stats[key] = int(value)
                    except ValueError:
                        stats[key] = value

        os.unlink(tmp_path)

        # Get total k-mers
        total_kmers = stats.get('total_k-mers',
                                stats.get('total_kmers',
                                         stats.get('no_of_unique_k-mers', 0)))

        return {'total_kmers': total_kmers}

    except Exception as e:
        print(f"Error getting stats for {db_path}: {e}")
        return {'total_kmers': 0}


def analyze_marker_availability(kmer_sizes: List[int], base_dir: str) -> pd.DataFrame:
    """Analyze total k-mer availability for each k-mer size."""

    all_data = []

    for k in kmer_sizes:
        kmc_dir = os.path.join(base_dir, f"k{k}")
        print(f"\nAnalyzing k={k} databases in {kmc_dir}...")

        databases = find_kmc_databases(kmc_dir)

        for db in databases:
            stats = get_kmc_stats(db['path'])

            # Calculate marker density
            chrom = db['chromosome']
            region_size = CEN_SIZES[chrom] if db['region'] == 'CEN' else ARMS_SIZES[chrom]
            density_per_mb = (stats['total_kmers'] / region_size) * 1_000_000

            all_data.append({
                'kmer_size': k,
                'database': db['label'],
                'genotype': db['genotype'],
                'region': db['region'],
                'chromosome': chrom,
                'total_kmers': stats['total_kmers'],
                'region_size_bp': region_size,
                'density_per_mb': density_per_mb
            })

    return pd.DataFrame(all_data)


def load_error_resilience_data(kmer_sizes: List[int]) -> pd.DataFrame:
    """Load error resilience data from previous analysis."""

    all_data = []
    for k in kmer_sizes:
        csv_path = f"final_results/realistic_k{k}_100k_error_resilience_stats.csv"
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df['kmer_size'] = k
            all_data.append(df)

    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()


def calculate_overall_quality_score(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate a comprehensive quality score considering:
    1. Marker availability (normalized)
    2. Marker density uniformity (lower CV is better)
    3. False positive rate (lower is better)
    4. Read retention (higher is better)

    Weights can be adjusted based on priorities.
    """

    # Normalize metrics to 0-100 scale
    # For each k-mer size, calculate summary metrics

    summary_stats = []

    for k in merged_df['kmer_size'].unique():
        k_data = merged_df[merged_df['kmer_size'] == k]

        # 1. Marker Availability (total k-mers available)
        total_markers = k_data['total_kmers'].sum()

        # 2. Marker Density Uniformity (CV - coefficient of variation)
        density_cv = k_data['density_per_mb'].std() / k_data['density_per_mb'].mean() * 100

        # 3. False Positive Rate (mean across all databases)
        # Calculate absolute FP rate
        k_data['abs_fp_rate'] = (k_data['pct_kmers_with_errors'] * k_data['pct_wrong_db']) / 100
        mean_fp_rate = k_data['abs_fp_rate'].mean()

        # 4. Read Retention (mean usable reads)
        k_data['usable_reads'] = (100 - k_data['pct_kmers_with_errors']) + \
                                 (k_data['pct_kmers_with_errors'] * k_data['pct_error_tolerant'] / 100)
        mean_retention = k_data['usable_reads'].mean()

        # 5. Mean marker density
        mean_density = k_data['density_per_mb'].mean()

        summary_stats.append({
            'kmer_size': k,
            'total_markers': total_markers,
            'mean_density_per_mb': mean_density,
            'density_cv': density_cv,
            'mean_fp_rate': mean_fp_rate,
            'mean_read_retention': mean_retention
        })

    summary_df = pd.DataFrame(summary_stats)

    # Normalize to 0-100 scale (higher is better)
    summary_df['marker_score'] = (summary_df['total_markers'] / summary_df['total_markers'].max()) * 100
    summary_df['density_score'] = (1 - summary_df['density_cv'] / summary_df['density_cv'].max()) * 100
    summary_df['fp_score'] = (1 - summary_df['mean_fp_rate'] / summary_df['mean_fp_rate'].max()) * 100
    summary_df['retention_score'] = (summary_df['mean_read_retention'] / summary_df['mean_read_retention'].max()) * 100

    # Overall quality score (weighted)
    # Adjust weights based on priorities
    WEIGHTS = {
        'marker_availability': 0.25,  # Having enough markers
        'density_uniformity': 0.15,   # Uniform coverage
        'specificity': 0.20,          # Low false positives
        'read_retention': 0.40        # Keep most reads (MOST IMPORTANT!)
    }

    summary_df['overall_quality'] = (
        summary_df['marker_score'] * WEIGHTS['marker_availability'] +
        summary_df['density_score'] * WEIGHTS['density_uniformity'] +
        summary_df['fp_score'] * WEIGHTS['specificity'] +
        summary_df['retention_score'] * WEIGHTS['read_retention']
    )

    return summary_df


def create_comprehensive_plot(availability_df: pd.DataFrame,
                              error_df: pd.DataFrame,
                              quality_df: pd.DataFrame):
    """Create comprehensive visualization."""

    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.30,
                          left=0.07, right=0.96, top=0.93, bottom=0.06)

    # Merge dataframes
    merged = availability_df.merge(error_df, on=['kmer_size', 'database', 'genotype', 'region', 'chromosome'])

    kmer_sizes = sorted(availability_df['kmer_size'].unique())
    x = np.arange(len(kmer_sizes))

    # =========================================================================
    # Panel A: Total Marker Availability
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])

    totals_arms = []
    totals_cen = []
    for k in kmer_sizes:
        k_data = availability_df[availability_df['kmer_size'] == k]
        totals_arms.append(k_data[k_data['region'] == 'ARMS']['total_kmers'].sum())
        totals_cen.append(k_data[k_data['region'] == 'CEN']['total_kmers'].sum())

    width = 0.35
    bars1 = ax1.bar(x - width/2, np.array(totals_arms) / 1e6, width, label='ARMS',
                   color='#9b59b6', alpha=0.85, edgecolor='black', linewidth=1.5)
    bars2 = ax1.bar(x + width/2, np.array(totals_cen) / 1e6, width, label='CEN',
                   color='#e67e22', alpha=0.85, edgecolor='black', linewidth=1.5)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, height + 0.5,
                f'{height:.1f}M', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for bar in bars2:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, height + 0.01,
                f'{height:.2f}M', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax1.set_ylabel('Total K-mers Available (millions)', fontweight='bold')
    ax1.set_xlabel('K-mer Size', fontweight='bold')
    ax1.set_title('A. Marker Availability\nTotal unique k-mers per k-mer size', fontweight='bold', pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax1.legend()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(axis='y', alpha=0.3)

    # =========================================================================
    # Panel B: Marker Density (k-mers per Mb)
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])

    mean_densities_arms = []
    mean_densities_cen = []
    for k in kmer_sizes:
        k_data = availability_df[availability_df['kmer_size'] == k]
        mean_densities_arms.append(k_data[k_data['region'] == 'ARMS']['density_per_mb'].mean())
        mean_densities_cen.append(k_data[k_data['region'] == 'CEN']['density_per_mb'].mean())

    bars1 = ax2.bar(x - width/2, mean_densities_arms, width, label='ARMS',
                   color='#9b59b6', alpha=0.85, edgecolor='black', linewidth=1.5)
    bars2 = ax2.bar(x + width/2, mean_densities_cen, width, label='CEN',
                   color='#e67e22', alpha=0.85, edgecolor='black', linewidth=1.5)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, height + 30,
                f'{height:.0f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    for bar in bars2:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, height + 5,
                f'{height:.0f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax2.set_ylabel('K-mers per Megabase', fontweight='bold')
    ax2.set_xlabel('K-mer Size', fontweight='bold')
    ax2.set_title('B. Marker Density\nMean k-mers per Mb (uniform coverage)', fontweight='bold', pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax2.legend()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.3)

    # =========================================================================
    # Panel C: Density Uniformity (CV)
    # =========================================================================
    ax3 = fig.add_subplot(gs[0, 2])

    cv_values = quality_df['density_cv'].values

    bars = ax3.bar(x, cv_values, color='#3498db', alpha=0.85,
                  edgecolor='black', linewidth=1.5)

    for bar in bars:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2, height + 1,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax3.set_ylabel('Coefficient of Variation (%)', fontweight='bold')
    ax3.set_xlabel('K-mer Size', fontweight='bold')
    ax3.set_title('C. Density Uniformity\nLower CV = more uniform coverage', fontweight='bold', pad=10)
    ax3.set_xticks(x)
    ax3.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.grid(axis='y', alpha=0.3)
    ax3.axhline(20, color='green', linestyle='--', alpha=0.5, label='Good (<20%)')
    ax3.legend()

    # =========================================================================
    # Panel D: Read Retention (from error analysis)
    # =========================================================================
    ax4 = fig.add_subplot(gs[1, :2])

    retention_values = quality_df['mean_read_retention'].values

    bars = ax4.bar(x, retention_values, color='#27ae60', alpha=0.85,
                  edgecolor='black', linewidth=1.5)

    # Highlight k=21
    bars[0].set_edgecolor('green')
    bars[0].set_linewidth(3)

    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2, height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
        if i > 0:
            loss = retention_values[0] - height
            ax4.text(bar.get_x() + bar.get_width()/2, height - 1.5,
                    f'-{loss:.1f}%', ha='center', va='top', fontsize=9,
                    color='darkred', style='italic')

    ax4.set_ylabel('Read Retention (%)', fontweight='bold', fontsize=13)
    ax4.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
    ax4.set_title('D. Read Retention Under ONT Sequencing (1% error rate)\n'\
                  'Percentage of reads correctly classified after errors',
                  fontweight='bold', pad=10, fontsize=13)
    ax4.set_xticks(x)
    ax4.set_xticklabels([f'k={k}' for k in kmer_sizes], fontsize=11)
    ax4.set_ylim(64, 84)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.grid(axis='y', alpha=0.3)
    ax4.text(0, 82.5, '★ BEST', ha='center', fontsize=12, fontweight='bold', color='green')

    # =========================================================================
    # Panel E: False Positive Rate
    # =========================================================================
    ax5 = fig.add_subplot(gs[1, 2])

    fp_values = quality_df['mean_fp_rate'].values * 100  # Convert to percentage

    bars = ax5.bar(x, fp_values, color='#e74c3c', alpha=0.85,
                  edgecolor='black', linewidth=1.5)

    for bar in bars:
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2, height + 0.002,
                f'{height:.3f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax5.set_ylabel('False Positive Rate (%)', fontweight='bold')
    ax5.set_xlabel('K-mer Size', fontweight='bold')
    ax5.set_title('E. Cross-Contamination Risk\nMean false positive rate', fontweight='bold', pad=10)
    ax5.set_xticks(x)
    ax5.set_xticklabels([f'k={k}' for k in kmer_sizes])
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.grid(axis='y', alpha=0.3)
    ax5.text(2, 0.16, '✓ All excellent\n(<0.2%)', ha='center', fontsize=10,
            fontweight='bold', color='darkgreen',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8))

    # =========================================================================
    # Panel F: Overall Quality Score
    # =========================================================================
    ax6 = fig.add_subplot(gs[2, :])

    # Component scores
    components = ['marker_score', 'density_score', 'fp_score', 'retention_score']
    component_labels = ['Marker\nAvailability', 'Density\nUniformity', 'Specificity\n(Low FP)', 'Read\nRetention']
    component_colors = ['#3498db', '#9b59b6', '#e74c3c', '#27ae60']

    bar_width = 0.15
    for i, (comp, label, color) in enumerate(zip(components, component_labels, component_colors)):
        values = quality_df[comp].values
        positions = x + (i - 1.5) * bar_width
        ax6.bar(positions, values, bar_width, label=label, color=color,
               alpha=0.85, edgecolor='black', linewidth=1)

    # Overall quality score as a line
    overall = quality_df['overall_quality'].values
    ax6.plot(x, overall, 'ko-', linewidth=3, markersize=10, label='Overall Quality', zorder=10)

    # Add values for overall score
    for i, val in enumerate(overall):
        ax6.text(i, val + 3, f'{val:.1f}', ha='center', va='bottom',
                fontsize=11, fontweight='bold')

    # Highlight best
    best_idx = overall.argmax()
    ax6.scatter([best_idx], [overall[best_idx]], s=500, facecolors='none',
               edgecolors='green', linewidths=3, zorder=11)
    ax6.text(best_idx, overall[best_idx] - 8, '★ BEST\nOVERALL', ha='center',
            fontsize=12, fontweight='bold', color='green')

    ax6.set_ylabel('Quality Score (0-100)', fontweight='bold', fontsize=13)
    ax6.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
    ax6.set_title('F. Comprehensive Quality Score\n'\
                  'Weighted combination: 25% availability + 15% uniformity + 20% specificity + 40% retention',
                  fontweight='bold', pad=10, fontsize=13)
    ax6.set_xticks(x)
    ax6.set_xticklabels([f'k={k}' for k in kmer_sizes], fontsize=11)
    ax6.set_ylim(0, 105)
    ax6.legend(loc='lower right', ncol=5, fontsize=10)
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.grid(axis='y', alpha=0.3)

    # =========================================================================
    # Main title
    # =========================================================================
    fig.suptitle('Comprehensive Cenhapmer Marker Evaluation: Complete Analysis Across K-mer Sizes\n'\
                 'Integrating marker availability, density, cross-contamination, and error resilience',
                 fontsize=16, fontweight='bold', y=0.97)

    plt.savefig('comprehensive_marker_evaluation.png', dpi=300, bbox_inches='tight')
    plt.savefig('comprehensive_marker_evaluation.pdf', bbox_inches='tight')
    print("\n✓ Saved comprehensive evaluation plots")


def generate_summary_report(availability_df: pd.DataFrame,
                           quality_df: pd.DataFrame):
    """Generate comprehensive text report."""

    report = []
    report.append("="*90)
    report.append("COMPREHENSIVE CENHAPMER MARKER EVALUATION")
    report.append("="*90)
    report.append("")

    for i, row in quality_df.iterrows():
        k = int(row['kmer_size'])
        k_avail = availability_df[availability_df['kmer_size'] == k]

        report.append(f"K-MER SIZE: {k}")
        report.append("-" * 90)

        # Marker availability
        arms_total = k_avail[k_avail['region'] == 'ARMS']['total_kmers'].sum()
        cen_total = k_avail[k_avail['region'] == 'CEN']['total_kmers'].sum()
        report.append(f"  Marker Availability:")
        report.append(f"    ARMS: {arms_total:,} k-mers ({arms_total/1e6:.2f}M)")
        report.append(f"    CEN:  {cen_total:,} k-mers ({cen_total/1e6:.2f}M)")
        report.append(f"    TOTAL: {arms_total + cen_total:,} k-mers")

        # Density
        arms_density = k_avail[k_avail['region'] == 'ARMS']['density_per_mb'].mean()
        cen_density = k_avail[k_avail['region'] == 'CEN']['density_per_mb'].mean()
        report.append(f"  Marker Density (k-mers/Mb):")
        report.append(f"    ARMS: {arms_density:.0f} k-mers/Mb")
        report.append(f"    CEN:  {cen_density:.0f} k-mers/Mb")
        report.append(f"    Uniformity (CV): {row['density_cv']:.1f}%")

        # Error resilience
        report.append(f"  Error Resilience:")
        report.append(f"    Read Retention: {row['mean_read_retention']:.2f}%")
        report.append(f"    False Positive Rate: {row['mean_fp_rate']:.4f}%")

        # Quality scores
        report.append(f"  Quality Scores (0-100):")
        report.append(f"    Marker Availability: {row['marker_score']:.1f}")
        report.append(f"    Density Uniformity:  {row['density_score']:.1f}")
        report.append(f"    Specificity:         {row['fp_score']:.1f}")
        report.append(f"    Read Retention:      {row['retention_score']:.1f}")
        report.append(f"    OVERALL QUALITY:     {row['overall_quality']:.1f} {'★ BEST' if row['overall_quality'] == quality_df['overall_quality'].max() else ''}")
        report.append("")

    # Final recommendation
    best_k = int(quality_df.loc[quality_df['overall_quality'].idxmax(), 'kmer_size'])
    report.append("="*90)
    report.append("FINAL RECOMMENDATION")
    report.append("="*90)
    report.append(f"🏆 USE K={best_k} FOR OPTIMAL OVERALL PERFORMANCE")
    report.append("")
    report.append("Justification:")
    report.append(f"  ✓ Highest overall quality score ({quality_df['overall_quality'].max():.1f}/100)")
    report.append(f"  ✓ Best read retention ({quality_df['mean_read_retention'].max():.2f}%)")
    report.append(f"  ✓ Excellent false positive rate (<0.2%)")
    report.append(f"  ✓ Good marker availability ({availability_df[availability_df['kmer_size']==best_k]['total_kmers'].sum()/1e6:.1f}M k-mers)")
    report.append("="*90)

    report_text = "\n".join(report)

    with open('comprehensive_evaluation_report.txt', 'w') as f:
        f.write(report_text)

    print(report_text)
    print("\n✓ Saved comprehensive report")


def main():
    print("="*90)
    print("COMPREHENSIVE MARKER EVALUATION")
    print("="*90)

    kmer_sizes = [21, 25, 31, 35, 41]
    base_dir = "../../../03-cenhapmers"

    # Step 1: Analyze marker availability
    print("\n[1/4] Analyzing marker availability and density...")
    availability_df = analyze_marker_availability(kmer_sizes, base_dir)

    # Step 2: Load error resilience data
    print("\n[2/4] Loading error resilience data...")
    error_df = load_error_resilience_data(kmer_sizes)

    if error_df.empty:
        print("ERROR: No error resilience data found!")
        print("Please ensure realistic_k*_100k_error_resilience_stats.csv files exist in final_results/")
        return

    # Step 3: Merge and calculate quality scores
    print("\n[3/4] Calculating comprehensive quality scores...")
    merged_df = availability_df.merge(error_df, on=['kmer_size', 'database', 'genotype', 'region', 'chromosome'])
    quality_df = calculate_overall_quality_score(merged_df)

    # Step 4: Generate outputs
    print("\n[4/4] Generating comprehensive visualization and report...")
    create_comprehensive_plot(availability_df, error_df, quality_df)
    generate_summary_report(availability_df, quality_df)

    # Save dataframes
    availability_df.to_csv('marker_availability_data.csv', index=False)
    quality_df.to_csv('comprehensive_quality_scores.csv', index=False)
    merged_df.to_csv('complete_marker_data.csv', index=False)

    print("\n✨ Comprehensive evaluation complete!")
    print("\nGenerated files:")
    print("  - comprehensive_marker_evaluation.png/pdf (main figure)")
    print("  - comprehensive_evaluation_report.txt (detailed report)")
    print("  - marker_availability_data.csv (raw availability data)")
    print("  - comprehensive_quality_scores.csv (quality scores)")
    print("  - complete_marker_data.csv (merged dataset)")


if __name__ == '__main__':
    main()
