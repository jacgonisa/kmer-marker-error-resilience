#!/usr/bin/env python3
"""
Error Resilience Analysis - Simulating Sequencing Errors

This script tests marker resilience under realistic sequencing error conditions.

Key Question: If a k-mer gets 1 sequencing error, does it still map uniquely
to the correct database, or does it become ambiguous?

Strategy:
1. Sample k-mers from each database
2. For each k-mer, introduce 1 random sequencing error (substitution)
3. Check if the mutated k-mer:
   - Still only matches the original database (GOOD - error-tolerant)
   - Matches a different database (BAD - becomes ambiguous)
   - Matches multiple databases (WORST - complete loss of specificity)

This directly answers: "How do k=21 vs k=31 markers perform with real sequencing errors?"
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
import matplotlib.patches as mpatches
from pathlib import Path
from typing import Dict, Set, List, Tuple, Optional
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')


# --- KMC FUNCTIONS (reuse from previous scripts) ---

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


# --- ERROR SIMULATION ---

def introduce_sequencing_error(kmer: str, error_rate: float = 1.0, seed: int = None) -> List[str]:
    """
    Introduce sequencing errors (substitutions) into a k-mer.

    Args:
        kmer: Original k-mer sequence
        error_rate: Number of errors to introduce (1.0 = 1 error, 2.0 = 2 errors)
        seed: Random seed for reproducibility

    Returns:
        List of mutated k-mers (one per error rate position)
    """
    if seed is not None:
        random.seed(seed)

    bases = ['A', 'C', 'G', 'T']
    n_errors = int(error_rate)

    mutated_kmers = []

    # For each error position
    for _ in range(n_errors):
        # Choose a random position
        pos = random.randint(0, len(kmer) - 1)
        original_base = kmer[pos]

        # Choose a different base
        alternative_bases = [b for b in bases if b != original_base]
        new_base = random.choice(alternative_bases)

        # Create mutated k-mer
        mutated = kmer[:pos] + new_base + kmer[pos+1:]
        mutated_kmers.append(mutated)

    return mutated_kmers


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
                             n_sample_per_db: int = 1000,
                             n_errors: int = 1,
                             seed: int = 42,
                             show_progress: bool = True) -> Tuple[Dict, List[Dict]]:
    """
    Analyze marker resilience under sequencing errors.

    For each sampled k-mer:
    1. Introduce n_errors random substitutions
    2. Check if mutated k-mer still maps uniquely to original database

    Returns:
        error_resilience_data: Per-database error resilience statistics
        error_events: List of all error simulation events
    """

    # Step 1: Sample k-mers from all databases
    print(f"\nStep 1: Sampling {n_sample_per_db} k-mers from each database...")

    db_kmers = {}  # db_label -> list of sampled k-mer sequences
    all_db_kmers = {}  # db_label -> set of ALL k-mers (for matching)

    for db in databases:
        sampled = sample_kmers_from_db(db['path'], n_sample_per_db, seed=seed)
        db_kmers[db['label']] = sampled
        all_db_kmers[db['label']] = set(sampled)  # For fast lookup

        if show_progress:
            print(f"  {db['label']}: sampled {len(sampled)} k-mers")

    # Step 2: Simulate errors and test resilience
    print(f"\nStep 2: Simulating {n_errors} sequencing error(s) per k-mer...")

    error_resilience_data = {}
    error_events = []

    for db in databases:
        db_label = db['label']
        test_kmers = db_kmers[db_label]

        n_error_tolerant = 0  # Mutated k-mer still unique to this database
        n_becomes_ambiguous = 0  # Mutated k-mer matches multiple databases
        n_becomes_wrong = 0  # Mutated k-mer matches only OTHER database
        n_becomes_novel = 0  # Mutated k-mer doesn't match any database

        for i, original_kmer in enumerate(test_kmers):
            # Introduce error
            mutated_kmers = introduce_sequencing_error(original_kmer, n_errors, seed=seed+i)

            for mutated_kmer in mutated_kmers:
                # Check which databases match
                matches = check_kmer_matches(mutated_kmer, all_db_kmers)

                # Classify outcome
                if len(matches) == 0:
                    n_becomes_novel += 1
                    outcome = 'novel'
                elif len(matches) == 1 and matches[0] == db_label:
                    n_error_tolerant += 1
                    outcome = 'error_tolerant'
                elif len(matches) == 1 and matches[0] != db_label:
                    n_becomes_wrong += 1
                    outcome = 'wrong_db'
                else:  # len(matches) > 1
                    n_becomes_ambiguous += 1
                    outcome = 'ambiguous'

                # Record event
                error_events.append({
                    'original_db': db_label,
                    'original_kmer': original_kmer,
                    'mutated_kmer': mutated_kmer,
                    'outcome': outcome,
                    'matches': ','.join(matches) if matches else 'none'
                })

        # Calculate statistics
        n_tested = len(test_kmers) * n_errors

        error_resilience_data[db_label] = {
            'database': db_label,
            'genotype': db['genotype'],
            'region': db['region'],
            'chromosome': db['chromosome'],
            'n_tested': n_tested,
            'n_error_tolerant': n_error_tolerant,
            'n_becomes_ambiguous': n_becomes_ambiguous,
            'n_becomes_wrong': n_becomes_wrong,
            'n_becomes_novel': n_becomes_novel,
            'pct_error_tolerant': 100 * n_error_tolerant / n_tested,
            'pct_becomes_ambiguous': 100 * n_becomes_ambiguous / n_tested,
            'pct_becomes_wrong': 100 * n_becomes_wrong / n_tested,
            'pct_becomes_novel': 100 * n_becomes_novel / n_tested,
        }

        if show_progress:
            print(f"  {db_label}: {n_error_tolerant}/{n_tested} error-tolerant "
                  f"({100*n_error_tolerant/n_tested:.1f}%)")

    return error_resilience_data, error_events


# --- VISUALIZATION ---

def plot_error_resilience(error_resilience_data: Dict, databases: List[Dict],
                         n_errors: int, output_prefix: str):
    """Plot error resilience statistics."""

    labels = [db['label'] for db in databases]
    n = len(labels)

    # Extract data
    error_tolerant = [error_resilience_data[l]['pct_error_tolerant'] for l in labels]
    becomes_ambiguous = [error_resilience_data[l]['pct_becomes_ambiguous'] for l in labels]
    becomes_wrong = [error_resilience_data[l]['pct_becomes_wrong'] for l in labels]
    becomes_novel = [error_resilience_data[l]['pct_becomes_novel'] for l in labels]

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(16, 8))

    x = np.arange(n)
    width = 0.8

    colors_geno = ['#e41a1c' if 'Col-0' in l else '#377eb8' for l in labels]

    # Stacked bars
    p1 = ax.bar(x, error_tolerant, width, label='Error-tolerant (still unique)',
                color='#2ca02c', alpha=0.8)
    p2 = ax.bar(x, becomes_novel, width, bottom=error_tolerant,
                label='Becomes novel (no match)', color='#d62728', alpha=0.8)
    p3 = ax.bar(x, becomes_wrong, width,
                bottom=np.array(error_tolerant) + np.array(becomes_novel),
                label='Matches wrong database', color='#ff7f0e', alpha=0.8)
    p4 = ax.bar(x, becomes_ambiguous, width,
                bottom=np.array(error_tolerant) + np.array(becomes_novel) + np.array(becomes_wrong),
                label='Becomes ambiguous (multiple matches)', color='#9467bd', alpha=0.8)

    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title(f'Marker Error Resilience (with {n_errors} sequencing error(s))',
                fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=9)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_ylim(0, 100)
    ax.grid(axis='y', alpha=0.3)

    # Add horizontal line at 90%
    ax.axhline(90, color='green', linestyle='--', alpha=0.5, linewidth=1)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_error_resilience.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved error resilience plot to {output_prefix}_error_resilience.png")


def plot_error_comparison(error_resilience_data: Dict, databases: List[Dict],
                         output_prefix: str):
    """Compare error resilience across genotypes and regions."""

    # Group by region
    arms_data = [error_resilience_data[db['label']] for db in databases if db['region'] == 'ARMS']
    cen_data = [error_resilience_data[db['label']] for db in databases if db['region'] == 'CEN']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Error tolerance by region
    ax1 = axes[0]
    arms_tolerant = [d['pct_error_tolerant'] for d in arms_data]
    cen_tolerant = [d['pct_error_tolerant'] for d in cen_data]

    bp = ax1.boxplot([arms_tolerant, cen_tolerant], labels=['ARMS', 'CEN'], patch_artist=True)
    bp['boxes'][0].set_facecolor('#ff7f00')
    bp['boxes'][1].set_facecolor('#4daf4a')

    ax1.set_ylabel('% Error-Tolerant K-mers', fontsize=11)
    ax1.set_title('Error Resilience by Region', fontsize=12, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)

    # Plot 2: Outcome distribution
    ax2 = axes[1]

    outcomes = ['Error-tolerant', 'Novel', 'Wrong DB', 'Ambiguous']
    arms_means = [
        np.mean([d['pct_error_tolerant'] for d in arms_data]),
        np.mean([d['pct_becomes_novel'] for d in arms_data]),
        np.mean([d['pct_becomes_wrong'] for d in arms_data]),
        np.mean([d['pct_becomes_ambiguous'] for d in arms_data])
    ]
    cen_means = [
        np.mean([d['pct_error_tolerant'] for d in cen_data]),
        np.mean([d['pct_becomes_novel'] for d in cen_data]),
        np.mean([d['pct_becomes_wrong'] for d in cen_data]),
        np.mean([d['pct_becomes_ambiguous'] for d in cen_data])
    ]

    x = np.arange(len(outcomes))
    width = 0.35

    ax2.bar(x - width/2, arms_means, width, label='ARMS', color='#ff7f00', alpha=0.8)
    ax2.bar(x + width/2, cen_means, width, label='CEN', color='#4daf4a', alpha=0.8)

    ax2.set_ylabel('Percentage (%)', fontsize=11)
    ax2.set_title('Mean Outcome Distribution', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(outcomes, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_error_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved error comparison to {output_prefix}_error_comparison.png")


def generate_report(error_resilience_data: Dict, databases: List[Dict],
                   n_errors: int, output_prefix: str):
    """Generate comprehensive error resilience report."""

    report_path = f'{output_prefix}_error_resilience_report.txt'

    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("ERROR RESILIENCE ANALYSIS REPORT\n")
        f.write(f"Simulating {n_errors} sequencing error(s) per k-mer\n")
        f.write("=" * 70 + "\n\n")

        # Ranking
        ranked = sorted(error_resilience_data.items(),
                       key=lambda x: x[1]['pct_error_tolerant'], reverse=True)

        f.write("ERROR RESILIENCE RANKING (best first)\n")
        f.write("-" * 40 + "\n")
        for i, (label, data) in enumerate(ranked, 1):
            f.write(f"{i}. {label}: {data['pct_error_tolerant']:.1f}% error-tolerant\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("DETAILED STATISTICS\n")
        f.write("=" * 70 + "\n\n")

        for db in databases:
            label = db['label']
            data = error_resilience_data[label]

            f.write(f"{label}\n")
            f.write(f"  Tested: {data['n_tested']} k-mers with errors\n")
            f.write(f"  Error-tolerant: {data['pct_error_tolerant']:.1f}%\n")
            f.write(f"  Becomes novel: {data['pct_becomes_novel']:.1f}%\n")
            f.write(f"  Matches wrong DB: {data['pct_becomes_wrong']:.1f}%\n")
            f.write(f"  Becomes ambiguous: {data['pct_becomes_ambiguous']:.1f}%\n")
            f.write("\n")

    print(f"Saved report to {report_path}")

    # Export CSV
    df = pd.DataFrame(error_resilience_data.values())
    df.to_csv(f'{output_prefix}_error_resilience_stats.csv', index=False)


# --- MAIN ---

def main():
    parser = argparse.ArgumentParser(
        description='Error resilience analysis with sequencing error simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script simulates sequencing errors and tests marker resilience.

Critical for comparing k=21 vs k=31 under realistic conditions:
- Longer k-mers have higher error probability
- Tests if errors cause markers to become ambiguous

Example:
  %(prog)s k21/ -o error_k21 --sample 1000 --errors 1
  %(prog)s k31/ -o error_k31 --sample 1000 --errors 1

Then compare error_k21 vs error_k31 results!
        """
    )

    parser.add_argument('input_dir', help='Directory with KMC databases')
    parser.add_argument('-o', '--output', default='error_resilience',
                        help='Output prefix')
    parser.add_argument('--sample', type=int, default=1000,
                        help='K-mers per database (default: 1000)')
    parser.add_argument('--errors', type=int, default=1,
                        help='Number of errors to simulate (default: 1)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')

    args = parser.parse_args()

    print("=" * 70)
    print("ERROR RESILIENCE ANALYSIS")
    print("=" * 70)

    # Find databases
    databases = find_kmc_databases(args.input_dir)

    if not databases:
        print("ERROR: No KMC databases found")
        return 1

    print(f"\nFound {len(databases)} databases")
    print(f"Will simulate {args.errors} sequencing error(s) per k-mer")

    # Analyze
    error_resilience_data, error_events = analyze_error_resilience(
        databases,
        n_sample_per_db=args.sample,
        n_errors=args.errors,
        seed=args.seed
    )

    # Visualize and report
    print("\nGenerating reports and visualizations...")
    plot_error_resilience(error_resilience_data, databases, args.errors, args.output)
    plot_error_comparison(error_resilience_data, databases, args.output)
    generate_report(error_resilience_data, databases, args.errors, args.output)

    # Summary
    mean_tolerant = np.mean([d['pct_error_tolerant'] for d in error_resilience_data.values()])
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Mean error tolerance: {mean_tolerant:.1f}%")

    return 0


if __name__ == '__main__':
    exit(main())
