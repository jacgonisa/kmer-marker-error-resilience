#!/usr/bin/env python3
"""
Error Resilience Analysis - Realistic ONT Sequencing Error Simulation

This script simulates realistic sequencing errors at a per-base error rate
(e.g., 1% per base, similar to ONT R9 chemistry).

Key differences from previous version:
- Each BASE has error_rate chance of being mutated (not 1 guaranteed error)
- This means some k-mers have 0 errors, some have 1, some have 2+
- More realistic for actual sequencing conditions
"""

import argparse
import subprocess
import tempfile
import os
import re
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Set, List, Tuple
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


# --- KMC FUNCTIONS ---

def find_kmc_databases(directory: str) -> List[Dict[str, str]]:
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


def sample_kmers_from_db(db_path: str, n_sample: int, seed: int = 42) -> List[str]:
    """Randomly sample k-mers from a KMC database."""
    random.seed(seed)

    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
            tmp_path = tmp.name

        subprocess.run(
            ['kmc_tools', 'transform', db_path, 'dump', tmp_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )

        all_kmers = []
        with open(tmp_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    all_kmers.append(parts[0])

        os.unlink(tmp_path)

        if len(all_kmers) <= n_sample:
            return all_kmers
        else:
            return random.sample(all_kmers, n_sample)

    except subprocess.CalledProcessError as e:
        print(f"Error sampling from {db_path}: {e}")
        return []


# --- REALISTIC ERROR SIMULATION ---

def introduce_realistic_sequencing_errors(kmer: str, error_rate: float = 0.01) -> Tuple[str, int]:
    """
    Introduce realistic sequencing errors with per-base error rate.

    Args:
        kmer: Original k-mer sequence
        error_rate: Per-base error probability (default 0.01 = 1%)

    Returns:
        (mutated_kmer, n_errors_introduced)
    """
    bases = ['A', 'C', 'G', 'T']
    mutated = list(kmer)
    n_errors = 0

    for i in range(len(kmer)):
        if random.random() < error_rate:
            # This base gets an error
            original_base = mutated[i]
            alternative_bases = [b for b in bases if b != original_base]
            mutated[i] = random.choice(alternative_bases)
            n_errors += 1

    return ''.join(mutated), n_errors


def check_kmer_matches(kmer: str, all_db_kmers: Dict[str, Set[str]]) -> List[str]:
    """
    Check which databases contain this k-mer.

    Returns list of database labels that contain this k-mer.
    """
    matches = []
    for db_label, kmer_set in all_db_kmers.items():
        if kmer in kmer_set:
            matches.append(db_label)
    return matches


# --- ANALYSIS ---

def analyze_error_resilience(databases: List[Dict],
                             n_sample_per_db: int = 100000,
                             error_rate: float = 0.01,
                             seed: int = 42) -> Tuple[Dict, pd.DataFrame]:
    """
    Analyze marker resilience under realistic sequencing errors.

    For each sampled k-mer:
    1. Apply per-base error rate to simulate realistic sequencing
    2. Count how many errors were introduced
    3. Check if mutated k-mer still maps uniquely to original database

    Returns:
        error_resilience_data: Per-database statistics
        events_df: DataFrame with all events for detailed analysis
    """

    random.seed(seed)

    # Step 1: Sample k-mers from all databases
    print(f"\nStep 1: Sampling {n_sample_per_db} k-mers from each database...")

    db_kmers = {}  # db_label -> list of sampled k-mer sequences
    all_db_kmers = {}  # db_label -> set of ALL k-mers (for matching)

    for db in databases:
        label = db['label']
        print(f"  Sampling from {label}...")

        sampled = sample_kmers_from_db(db['path'], n_sample_per_db, seed=seed + hash(label) % 10000)
        db_kmers[label] = sampled

        # Also load ALL k-mers for matching (for fast lookup)
        all_kmers = sample_kmers_from_db(db['path'], n_sample=999999999, seed=seed)
        all_db_kmers[label] = set(all_kmers)

        print(f"    Sampled: {len(sampled):,} k-mers, Total in DB: {len(all_kmers):,}")

    # Step 2: Simulate sequencing errors and check resilience
    print(f"\nStep 2: Simulating sequencing errors (per-base error rate: {error_rate*100:.1f}%)...")

    error_resilience_data = {}
    all_events = []

    for db in databases:
        label = db['label']
        sampled_kmers = db_kmers[label]

        print(f"\n  Processing {label} ({len(sampled_kmers):,} k-mers)...")

        # Track outcomes
        n_tested = 0
        n_had_errors = 0  # How many k-mers had at least 1 error
        n_error_tolerant = 0  # Still matches ONLY original DB after errors
        n_becomes_novel = 0  # Doesn't match any DB
        n_becomes_wrong = 0  # Matches a different DB
        n_becomes_ambiguous = 0  # Matches multiple DBs
        n_no_errors = 0  # No errors introduced

        error_count_dist = defaultdict(int)  # Distribution of error counts

        for i, original_kmer in enumerate(sampled_kmers):
            if i % 10000 == 0 and i > 0:
                print(f"    Progress: {i:,}/{len(sampled_kmers):,} k-mers...")

            # Introduce realistic errors
            mutated_kmer, n_errors = introduce_realistic_sequencing_errors(original_kmer, error_rate)

            n_tested += 1
            error_count_dist[n_errors] += 1

            if n_errors == 0:
                n_no_errors += 1
                # If no errors, it still matches correctly - don't count in results
                continue

            n_had_errors += 1

            # Check where the mutated k-mer matches
            matches = check_kmer_matches(mutated_kmer, all_db_kmers)

            # Classify outcome
            if len(matches) == 0:
                outcome = 'novel'
                n_becomes_novel += 1
            elif len(matches) == 1 and matches[0] == label:
                outcome = 'error_tolerant'
                n_error_tolerant += 1
            elif len(matches) == 1 and matches[0] != label:
                outcome = 'wrong_db'
                n_becomes_wrong += 1
            else:  # len(matches) > 1
                outcome = 'ambiguous'
                n_becomes_ambiguous += 1

            # Record event
            all_events.append({
                'database': label,
                'original_kmer': original_kmer,
                'mutated_kmer': mutated_kmer,
                'n_errors': n_errors,
                'outcome': outcome,
                'matches': ','.join(matches) if matches else 'none'
            })

        # Calculate percentages (out of k-mers that HAD errors)
        if n_had_errors > 0:
            pct_error_tolerant = (n_error_tolerant / n_had_errors) * 100
            pct_novel = (n_becomes_novel / n_had_errors) * 100
            pct_wrong = (n_becomes_wrong / n_had_errors) * 100
            pct_ambiguous = (n_becomes_ambiguous / n_had_errors) * 100
        else:
            pct_error_tolerant = pct_novel = pct_wrong = pct_ambiguous = 0.0

        # Store results
        error_resilience_data[label] = {
            'n_tested': n_tested,
            'n_had_errors': n_had_errors,
            'n_no_errors': n_no_errors,
            'pct_kmers_with_errors': (n_had_errors / n_tested * 100) if n_tested > 0 else 0,
            'n_error_tolerant': n_error_tolerant,
            'n_novel': n_becomes_novel,
            'n_wrong': n_becomes_wrong,
            'n_ambiguous': n_becomes_ambiguous,
            'pct_error_tolerant': pct_error_tolerant,
            'pct_novel': pct_novel,
            'pct_wrong': pct_wrong,
            'pct_ambiguous': pct_ambiguous,
            'error_count_dist': dict(error_count_dist),
            'mean_errors_per_kmer': sum(k*v for k,v in error_count_dist.items()) / n_tested
        }

        print(f"    Results:")
        print(f"      K-mers with errors: {n_had_errors:,}/{n_tested:,} ({n_had_errors/n_tested*100:.2f}%)")
        print(f"      Error-tolerant: {pct_error_tolerant:.2f}%")
        print(f"      Becomes novel: {pct_novel:.2f}%")
        print(f"      Mean errors/k-mer: {error_resilience_data[label]['mean_errors_per_kmer']:.3f}")

    # Convert events to DataFrame
    events_df = pd.DataFrame(all_events)

    return error_resilience_data, events_df


def generate_report(resilience_data: Dict, databases: List[Dict], output_prefix: str, error_rate: float):
    """Generate text report."""

    # Sort by error tolerance (descending)
    sorted_dbs = sorted(resilience_data.items(), key=lambda x: x[1]['pct_error_tolerant'], reverse=True)

    report_path = f"{output_prefix}_error_resilience_report.txt"
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("REALISTIC ERROR RESILIENCE ANALYSIS REPORT\n")
        f.write(f"Per-base sequencing error rate: {error_rate*100:.1f}% (ONT-like)\n")
        f.write("=" * 70 + "\n\n")

        f.write("ERROR RESILIENCE RANKING (best first)\n")
        f.write("-" * 40 + "\n")
        for i, (db_label, stats) in enumerate(sorted_dbs, 1):
            f.write(f"{i}. {db_label}: {stats['pct_error_tolerant']:.2f}% error-tolerant\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("DETAILED STATISTICS\n")
        f.write("=" * 70 + "\n\n")

        for db_label, stats in sorted_dbs:
            f.write(f"{db_label}\n")
            f.write(f"  Tested: {stats['n_tested']:,} k-mers\n")
            f.write(f"  Had sequencing errors: {stats['n_had_errors']:,} ({stats['pct_kmers_with_errors']:.2f}%)\n")
            f.write(f"  Mean errors per k-mer: {stats['mean_errors_per_kmer']:.3f}\n")
            f.write(f"  \n")
            f.write(f"  Of k-mers WITH errors:\n")
            f.write(f"    Error-tolerant: {stats['pct_error_tolerant']:.2f}%\n")
            f.write(f"    Becomes novel: {stats['pct_novel']:.2f}%\n")
            f.write(f"    Matches wrong DB: {stats['pct_wrong']:.2f}%\n")
            f.write(f"    Becomes ambiguous: {stats['pct_ambiguous']:.2f}%\n")
            f.write(f"\n")

    print(f"\n✓ Saved report: {report_path}")


def create_summary_csv(resilience_data: Dict, databases: List[Dict], output_prefix: str):
    """Create CSV summary."""

    rows = []
    for db in databases:
        label = db['label']
        stats = resilience_data[label]

        rows.append({
            'database': label,
            'genotype': db['genotype'],
            'region': db['region'],
            'chromosome': db['chromosome'],
            'n_tested': stats['n_tested'],
            'n_had_errors': stats['n_had_errors'],
            'pct_kmers_with_errors': stats['pct_kmers_with_errors'],
            'mean_errors_per_kmer': stats['mean_errors_per_kmer'],
            'pct_error_tolerant': stats['pct_error_tolerant'],
            'pct_becomes_novel': stats['pct_novel'],
            'pct_wrong_db': stats['pct_wrong'],
            'pct_ambiguous': stats['pct_ambiguous']
        })

    df = pd.DataFrame(rows)
    csv_path = f"{output_prefix}_error_resilience_stats.csv"
    df.to_csv(csv_path, index=False)
    print(f"✓ Saved CSV: {csv_path}")

    return df


def create_plots(resilience_data: Dict, databases: List[Dict], output_prefix: str, error_rate: float):
    """Create visualization plots."""

    # Prepare data
    plot_data = []
    for db in databases:
        label = db['label']
        stats = resilience_data[label]
        plot_data.append({
            'Database': label,
            'Genotype': db['genotype'],
            'Region': db['region'],
            'Error Tolerant (%)': stats['pct_error_tolerant'],
            'Becomes Novel (%)': stats['pct_novel'],
            'Wrong DB (%)': stats['pct_wrong'],
            'K-mers with Errors (%)': stats['pct_kmers_with_errors'],
            'Mean Errors': stats['mean_errors_per_kmer']
        })

    df = pd.DataFrame(plot_data)

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Error resilience by region
    ax1 = axes[0, 0]
    regions = ['ARMS', 'CEN']
    region_data = [df[df['Region'] == r]['Error Tolerant (%)'].values for r in regions]

    bp = ax1.boxplot(region_data, labels=regions, patch_artist=True)
    for patch, color in zip(bp['boxes'], ['#FF6B35', '#004E89']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax1.set_ylabel('% Error-Tolerant K-mers', fontweight='bold')
    ax1.set_title(f'Error Resilience by Region\n({error_rate*100:.1f}% per-base error rate)', fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)

    # Plot 2: K-mers affected by errors
    ax2 = axes[0, 1]
    arms_pct = df[df['Region'] == 'ARMS']['K-mers with Errors (%)'].mean()
    cen_pct = df[df['Region'] == 'CEN']['K-mers with Errors (%)'].mean()

    bars = ax2.bar(['ARMS', 'CEN'], [arms_pct, cen_pct], color=['#FF6B35', '#004E89'], alpha=0.7)
    ax2.set_ylabel('% K-mers with ≥1 Error', fontweight='bold')
    ax2.set_title('Sequencing Error Impact', fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)

    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2f}%', ha='center', va='bottom', fontweight='bold')

    # Plot 3: Outcome distribution
    ax3 = axes[1, 0]
    arms_df = df[df['Region'] == 'ARMS']
    cen_df = df[df['Region'] == 'CEN']

    x = np.arange(2)
    width = 0.2

    outcomes = ['Error Tolerant (%)', 'Becomes Novel (%)', 'Wrong DB (%)']
    colors_outcomes = ['#2E8B57', '#FFD700', '#DC143C']

    for i, outcome in enumerate(outcomes):
        arms_vals = [arms_df[outcome].mean()]
        cen_vals = [cen_df[outcome].mean()]
        ax3.bar(x[0] + i*width, arms_vals, width, label=outcome.replace(' (%)', ''),
                color=colors_outcomes[i], alpha=0.7)
        ax3.bar(x[1] + i*width, cen_vals, width, color=colors_outcomes[i], alpha=0.7)

    ax3.set_ylabel('Percentage', fontweight='bold')
    ax3.set_title('Outcome Distribution by Region', fontweight='bold')
    ax3.set_xticks(x + width)
    ax3.set_xticklabels(['ARMS', 'CEN'])
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)

    # Plot 4: Per-database heatmap
    ax4 = axes[1, 1]

    # Create matrix
    databases_sorted = sorted(df['Database'].values,
                              key=lambda x: (x.split('_')[0], x.split('_')[1], x.split('_')[2]))

    matrix_data = df.set_index('Database').loc[databases_sorted, 'Error Tolerant (%)'].values.reshape(-1, 1)

    im = ax4.imshow(matrix_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=max(1, matrix_data.max()))
    ax4.set_yticks(np.arange(len(databases_sorted)))
    ax4.set_yticklabels(databases_sorted, fontsize=8)
    ax4.set_xticks([0])
    ax4.set_xticklabels(['Error Tolerance'])
    ax4.set_title('Per-Database Error Tolerance', fontweight='bold')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax4)
    cbar.set_label('% Error-Tolerant', fontweight='bold')

    # Add text annotations
    for i in range(len(databases_sorted)):
        val = matrix_data[i, 0]
        text_color = 'white' if val > matrix_data.max()/2 else 'black'
        ax4.text(0, i, f'{val:.2f}', ha='center', va='center',
                color=text_color, fontsize=7, fontweight='bold')

    plt.tight_layout()

    png_path = f"{output_prefix}_error_comparison.png"
    pdf_path = f"{output_prefix}_error_comparison.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()

    print(f"✓ Saved plots: {png_path}, {pdf_path}")


# --- MAIN ---

def main():
    parser = argparse.ArgumentParser(
        description='Analyze marker resilience under realistic sequencing errors (per-base error rate)'
    )
    parser.add_argument('kmc_dir', help='Directory containing KMC databases')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('--sample', type=int, default=100000,
                       help='Number of k-mers to sample per database (default: 100000)')
    parser.add_argument('--error-rate', type=float, default=0.01,
                       help='Per-base sequencing error rate (default: 0.01 = 1%%)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')

    args = parser.parse_args()

    print("=" * 70)
    print("REALISTIC ERROR RESILIENCE ANALYSIS")
    print(f"Per-base error rate: {args.error_rate*100:.1f}% (similar to ONT sequencing)")
    print("=" * 70)

    # Find databases
    print("\nFinding KMC databases...")
    databases = find_kmc_databases(args.kmc_dir)

    if not databases:
        print("ERROR: No KMC databases found!")
        return

    print(f"Found {len(databases)} databases")

    # Run analysis
    resilience_data, events_df = analyze_error_resilience(
        databases,
        n_sample_per_db=args.sample,
        error_rate=args.error_rate,
        seed=args.seed
    )

    # Generate outputs
    print("\nGenerating outputs...")
    generate_report(resilience_data, databases, args.output, args.error_rate)
    create_summary_csv(resilience_data, databases, args.output)
    create_plots(resilience_data, databases, args.output, args.error_rate)

    # Save detailed events
    events_path = f"{args.output}_events.csv.gz"
    events_df.to_csv(events_path, index=False, compression='gzip')
    print(f"✓ Saved detailed events: {events_path}")

    print("\n✨ Analysis complete!")


if __name__ == '__main__':
    main()
