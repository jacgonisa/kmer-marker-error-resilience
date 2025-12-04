#!/usr/bin/env python3
"""
Plot 1: Marker Availability and Density Analysis
Shows total k-mer counts and marker density for each database.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import subprocess

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def count_kmers_in_database(db_path):
    """Count unique k-mers in a database file."""
    try:
        result = subprocess.run(
            f"kmc_tools info {db_path} 2>/dev/null | grep -i 'total k-mers' | awk '{{print $4}}'",
            shell=True,
            capture_output=True,
            text=True
        )
        if result.stdout.strip():
            return int(result.stdout.strip())
    except Exception as e:
        print(f"Error counting k-mers for {db_path}: {e}")
    return None

def get_region_size(genotype, chrom, region_type):
    """
    Estimate region size based on typical sizes.
    For ARMS: typically 10-50kb
    For CEN: typically 100-500kb
    """
    if region_type == "ARMS":
        return 30_000  # 30kb average
    else:  # CEN
        return 300_000  # 300kb average

# K-mer sizes to analyze
k_sizes = [21, 25, 31, 35, 41]
base_path = Path("../../../03-cenhapmers")

# Collect data
data = []

for k in k_sizes:
    k_dir = base_path / f"k{k}"
    if not k_dir.exists():
        print(f"Warning: {k_dir} does not exist, skipping k={k}")
        continue

    # Find all .kmc_pre files
    for kmc_file in k_dir.glob("*.kmc_pre"):
        db_name = kmc_file.stem
        db_path = str(kmc_file.parent / db_name)

        # Parse database name: e.g., "unique_Col-0_ARMS_Chr1_k21"
        # Remove "unique_" prefix and "_k##" suffix
        clean_name = db_name.replace('unique_', '')
        clean_name = clean_name.rsplit('_k', 1)[0]  # Remove _k## suffix

        parts = clean_name.split('_')
        if len(parts) >= 3:
            genotype = parts[0]
            region_type = parts[1]
            chrom = '_'.join(parts[2:])  # Handle cases like "Chr1" or "Chr_1"

            # Count k-mers
            kmer_count = count_kmers_in_database(db_path)
            if kmer_count is not None:
                # Estimate region size
                region_size = get_region_size(genotype, chrom, region_type)

                # Calculate density (k-mers per Mb)
                density = (kmer_count / region_size) * 1_000_000

                data.append({
                    'k_size': k,
                    'database': clean_name,
                    'genotype': genotype,
                    'region': region_type,
                    'chromosome': chrom,
                    'total_kmers': kmer_count,
                    'estimated_size_kb': region_size / 1000,
                    'density_per_Mb': density
                })
                print(f"k={k:2d} | {clean_name:25s} | {kmer_count:10,} k-mers | {density:10,.0f} k-mers/Mb")

# Create DataFrame
df = pd.DataFrame(data)

if df.empty:
    print("ERROR: No data collected! Check that KMC databases exist.")
    exit(1)

# Save summary
df.to_csv("final_results/marker_availability_summary.csv", index=False)

print(f"\n✓ Analyzed {len(df)} databases across {len(k_sizes)} k-mer sizes")
print(f"✓ Total unique combinations: {df.groupby(['k_size', 'region']).ngroups}")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Marker Availability and Density Analysis', fontsize=16, fontweight='bold', y=0.995)

# Panel A: Total k-mers by k-size and region
ax = axes[0, 0]
summary = df.groupby(['k_size', 'region'])['total_kmers'].agg(['mean', 'std']).reset_index()
x = np.arange(len(k_sizes))
width = 0.35

arms_data = summary[summary['region'] == 'ARMS'].sort_values('k_size')
cen_data = summary[summary['region'] == 'CEN'].sort_values('k_size')

bars1 = ax.bar(x - width/2, arms_data['mean'], width, label='ARMS',
               color='#3498db', yerr=arms_data['std'], capsize=5)
bars2 = ax.bar(x + width/2, cen_data['mean'], width, label='CEN',
               color='#e74c3c', yerr=cen_data['std'], capsize=5)

ax.set_xlabel('K-mer Size', fontweight='bold')
ax.set_ylabel('Average K-mer Count', fontweight='bold')
ax.set_title('A. Average Marker Count per Database', fontweight='bold', loc='left')
ax.set_xticks(x)
ax.set_xticklabels([f'k={k}' for k in k_sizes])
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height):,}',
                   ha='center', va='bottom', fontsize=8)

# Panel B: Marker density (k-mers per Mb)
ax = axes[0, 1]
density_summary = df.groupby(['k_size', 'region'])['density_per_Mb'].agg(['mean', 'std']).reset_index()

arms_dens = density_summary[density_summary['region'] == 'ARMS'].sort_values('k_size')
cen_dens = density_summary[density_summary['region'] == 'CEN'].sort_values('k_size')

bars1 = ax.bar(x - width/2, arms_dens['mean'], width, label='ARMS',
               color='#3498db', yerr=arms_dens['std'], capsize=5)
bars2 = ax.bar(x + width/2, cen_dens['mean'], width, label='CEN',
               color='#e74c3c', yerr=cen_dens['std'], capsize=5)

ax.set_xlabel('K-mer Size', fontweight='bold')
ax.set_ylabel('K-mers per Megabase', fontweight='bold')
ax.set_title('B. Marker Density (Coverage)', fontweight='bold', loc='left')
ax.set_xticks(x)
ax.set_xticklabels([f'k={k}' for k in k_sizes])
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height):,}',
                   ha='center', va='bottom', fontsize=8)

# Panel C: Distribution of k-mer counts across databases
ax = axes[1, 0]
for region, color in [('ARMS', '#3498db'), ('CEN', '#e74c3c')]:
    for k in k_sizes:
        subset = df[(df['k_size'] == k) & (df['region'] == region)]
        if not subset.empty:
            ax.scatter([k]*len(subset), subset['total_kmers'],
                      alpha=0.6, s=80, color=color,
                      label=region if k == k_sizes[0] else "")

ax.set_xlabel('K-mer Size', fontweight='bold')
ax.set_ylabel('Total K-mers', fontweight='bold')
ax.set_title('C. K-mer Count Distribution per Database', fontweight='bold', loc='left')
ax.set_xticks(k_sizes)
ax.set_xticklabels([f'k={k}' for k in k_sizes])
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Panel D: Coefficient of Variation (uniformity)
ax = axes[1, 1]
cv_data = []
for k in k_sizes:
    for region in ['ARMS', 'CEN']:
        subset = df[(df['k_size'] == k) & (df['region'] == region)]
        if not subset.empty and len(subset) > 1:
            mean_count = subset['total_kmers'].mean()
            std_count = subset['total_kmers'].std()
            cv = (std_count / mean_count) * 100 if mean_count > 0 else 0
            cv_data.append({'k_size': k, 'region': region, 'cv': cv})

cv_df = pd.DataFrame(cv_data)
if not cv_df.empty:
    arms_cv = cv_df[cv_df['region'] == 'ARMS'].sort_values('k_size')
    cen_cv = cv_df[cv_df['region'] == 'CEN'].sort_values('k_size')

    bars1 = ax.bar(x - width/2, arms_cv['cv'], width, label='ARMS', color='#3498db')
    bars2 = ax.bar(x + width/2, cen_cv['cv'], width, label='CEN', color='#e74c3c')

    ax.set_xlabel('K-mer Size', fontweight='bold')
    ax.set_ylabel('Coefficient of Variation (%)', fontweight='bold')
    ax.set_title('D. Marker Uniformity Across Databases', fontweight='bold', loc='left')
    ax.set_xticks(x)
    ax.set_xticklabels([f'k={k}' for k in k_sizes])
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    ax.axhline(y=20, color='green', linestyle='--', alpha=0.5, label='Good (<20%)')
    ax.text(len(k_sizes)-0.5, 22, '✓ Lower is better\n(more uniform)',
            ha='right', va='bottom', fontsize=9, color='green', fontweight='bold')

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.1f}%',
                       ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig('final_results/01_marker_availability.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/01_marker_availability.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/01_marker_availability.png")
print(f"✓ Saved: final_results/01_marker_availability.pdf")

# Print summary statistics
print("\n" + "="*80)
print("MARKER AVAILABILITY SUMMARY")
print("="*80)
for k in k_sizes:
    k_data = df[df['k_size'] == k]
    arms_mean = k_data[k_data['region'] == 'ARMS']['total_kmers'].mean()
    cen_mean = k_data[k_data['region'] == 'CEN']['total_kmers'].mean()
    print(f"k={k:2d} | ARMS: {arms_mean:10,.0f} k-mers | CEN: {cen_mean:10,.0f} k-mers")
print("="*80)
