#!/usr/bin/env python3
"""
Plot 6: Coverage Loss Analysis - Error-Affected K-mers per Megabase
This is the CORRECT metric: How much classification ability do you lose per Mb?
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("Set2")

# Load data
k_sizes = [21, 25, 31, 35, 41]
error_data = []

for k in k_sizes:
    csv_file = Path(f"final_results/realistic_k{k}_100k_error_resilience_stats.csv")
    if csv_file.exists():
        df = pd.read_csv(csv_file)
        df['k_size'] = k
        error_data.append(df)

if not error_data:
    print("ERROR: No error resilience data found!")
    exit(1)

error_df = pd.concat(error_data, ignore_index=True)

# Load marker availability (k-mer density)
marker_df = pd.read_csv("final_results/marker_availability_summary.csv")
error_df = error_df.merge(marker_df[['k_size', 'database', 'density_per_Mb']],
                           on=['k_size', 'database'], how='left')

print(f"✓ Loaded data for {len(error_df)} databases across {len(k_sizes)} k-mer sizes")

# Calculate coverage loss metrics
error_df['error_affected_kmers_per_Mb'] = (error_df['pct_kmers_with_errors'] / 100) * error_df['density_per_Mb']
error_df['usable_kmers_per_Mb'] = ((100 - error_df['pct_kmers_with_errors']) / 100) * error_df['density_per_Mb']

# Calculate summary statistics
summary_stats = []
for k in k_sizes:
    k_data = error_df[error_df['k_size'] == k]

    for region in ['ARMS', 'CEN']:
        region_data = k_data[k_data['region'] == region]

        summary_stats.append({
            'k_size': k,
            'region': region,
            'avg_density_per_Mb': region_data['density_per_Mb'].mean(),
            'avg_error_rate': region_data['pct_kmers_with_errors'].mean(),
            'avg_lost_kmers_per_Mb': region_data['error_affected_kmers_per_Mb'].mean(),
            'avg_usable_kmers_per_Mb': region_data['usable_kmers_per_Mb'].mean(),
            'coverage_retention': ((100 - region_data['pct_kmers_with_errors'].mean()) / 100) * 100
        })

summary_df = pd.DataFrame(summary_stats)

print("\n" + "="*100)
print("COVERAGE LOSS ANALYSIS: Error-Affected K-mers per Megabase")
print("="*100)
print(f"{'K-mer':<8} {'Region':<8} {'Density/Mb':<15} {'Error %':<12} {'Lost/Mb':<20} {'Usable/Mb':<20}")
print("-"*100)
for _, row in summary_df.iterrows():
    print(f"k={row['k_size']:<5} {row['region']:<8} "
          f"{row['avg_density_per_Mb']:>12,.0f}    "
          f"{row['avg_error_rate']:>6.1f}%     "
          f"{row['avg_lost_kmers_per_Mb']:>16,.0f}    "
          f"{row['avg_usable_kmers_per_Mb']:>16,.0f}")
print("="*100)

# Create visualization
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

fig.suptitle('Coverage Loss Analysis: Error-Affected K-mers per Megabase\n'
             'The Correct Metric for ONT Sequencing Performance',
             fontsize=16, fontweight='bold', y=0.98)

# Panel A: Lost k-mers per Mb (THE KEY METRIC!)
ax1 = fig.add_subplot(gs[0, :])

x = np.arange(len(k_sizes))
width = 0.35

arms_lost = summary_df[summary_df['region'] == 'ARMS']['avg_lost_kmers_per_Mb'].values
cen_lost = summary_df[summary_df['region'] == 'CEN']['avg_lost_kmers_per_Mb'].values

bars1 = ax1.bar(x - width/2, arms_lost / 1_000_000, width, label='ARMS',
                color='#66c2a5', edgecolor='black', linewidth=1.5)
bars2 = ax1.bar(x + width/2, cen_lost / 1_000_000, width, label='CEN',
                color='#fc8d62', edgecolor='black', linewidth=1.5)

ax1.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
ax1.set_ylabel('Lost K-mers per Mb (millions)', fontweight='bold', fontsize=13)
ax1.set_title('A. Coverage Loss: Error-Affected K-mers per Megabase ⭐ KEY METRIC',
              fontweight='bold', loc='left', fontsize=14, color='darkred')
ax1.set_xticks(x)
ax1.set_xticklabels([f'k={k}' for k in k_sizes])
ax1.legend(fontsize=12)
ax1.grid(axis='y', alpha=0.3)

# Highlight k=21 as best
bars1[0].set_edgecolor('green')
bars1[0].set_linewidth(3)
bars2[0].set_edgecolor('green')
bars2[0].set_linewidth(3)

# Add value labels
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}M',
                ha='center', va='bottom', fontsize=10, fontweight='bold')

# Add annotation
ax1.text(0.02, 0.95, '← k=21: LEAST coverage loss (highlighted in green)',
         transform=ax1.transAxes, fontsize=12, color='green', fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.4))

ax1.text(0.98, 0.95, 'Higher = MORE coverage lost →',
         transform=ax1.transAxes, fontsize=12, color='red', fontweight='bold',
         ha='right', bbox=dict(boxstyle='round', facecolor='#ffcccc', alpha=0.4))

# Panel B: Total k-mer density per Mb
ax2 = fig.add_subplot(gs[1, 0])

arms_density = summary_df[summary_df['region'] == 'ARMS']['avg_density_per_Mb'].values
cen_density = summary_df[summary_df['region'] == 'CEN']['avg_density_per_Mb'].values

ax2.plot(k_sizes, arms_density / 1_000_000, 'o-', linewidth=3, markersize=10,
         label='ARMS', color='#66c2a5')
ax2.plot(k_sizes, cen_density / 1_000_000, 'o-', linewidth=3, markersize=10,
         label='CEN', color='#fc8d62')

ax2.set_xlabel('K-mer Size', fontweight='bold')
ax2.set_ylabel('K-mer Density (millions/Mb)', fontweight='bold')
ax2.set_title('B. Total K-mer Density\n(before errors)', fontweight='bold', loc='left')
ax2.legend()
ax2.grid(alpha=0.3)

# Add annotation
ax2.text(0.5, 0.95, 'More k-mers = More opportunities to classify',
         transform=ax2.transAxes, ha='center', fontsize=9,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.6))

# Panel C: Error rate (% k-mers affected)
ax3 = fig.add_subplot(gs[1, 1])

arms_error = summary_df[summary_df['region'] == 'ARMS']['avg_error_rate'].values
cen_error = summary_df[summary_df['region'] == 'CEN']['avg_error_rate'].values

ax3.plot(k_sizes, arms_error, 'o-', linewidth=3, markersize=10,
         label='ARMS', color='#66c2a5')
ax3.plot(k_sizes, cen_error, 'o-', linewidth=3, markersize=10,
         label='CEN', color='#fc8d62')

ax3.set_xlabel('K-mer Size', fontweight='bold')
ax3.set_ylabel('K-mers Affected by Errors (%)', fontweight='bold')
ax3.set_title('C. Error Impact Rate', fontweight='bold', loc='left')
ax3.legend()
ax3.grid(alpha=0.3)

# Add annotation
ax3.text(0.5, 0.05, 'Longer k-mers = More likely to contain errors',
         transform=ax3.transAxes, ha='center', fontsize=9,
         bbox=dict(boxstyle='round', facecolor='#ffcccc', alpha=0.6))

# Panel D: Usable k-mers per Mb
ax4 = fig.add_subplot(gs[1, 2])

arms_usable = summary_df[summary_df['region'] == 'ARMS']['avg_usable_kmers_per_Mb'].values
cen_usable = summary_df[summary_df['region'] == 'CEN']['avg_usable_kmers_per_Mb'].values

ax4.plot(k_sizes, arms_usable / 1_000_000, 'o-', linewidth=3, markersize=10,
         label='ARMS', color='#66c2a5')
ax4.plot(k_sizes, cen_usable / 1_000_000, 'o-', linewidth=3, markersize=10,
         label='CEN', color='#fc8d62')

ax4.set_xlabel('K-mer Size', fontweight='bold')
ax4.set_ylabel('Usable K-mers per Mb (millions)', fontweight='bold')
ax4.set_title('D. Remaining Classification Ability', fontweight='bold', loc='left')
ax4.legend()
ax4.grid(alpha=0.3)

# Highlight best
ax4.axvline(x=21, color='green', linestyle='--', alpha=0.5, linewidth=2)
ax4.text(21, ax4.get_ylim()[1] * 0.95, 'k=21\nBEST', ha='center', fontsize=10,
         color='green', fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

# Panel E: Stacked bar showing composition
ax5 = fig.add_subplot(gs[2, :2])

# Calculate for overall (averaging ARMS and CEN)
overall_usable = []
overall_lost = []
for k in k_sizes:
    k_summary = summary_df[summary_df['k_size'] == k]
    avg_usable = k_summary['avg_usable_kmers_per_Mb'].mean()
    avg_lost = k_summary['avg_lost_kmers_per_Mb'].mean()
    overall_usable.append(avg_usable / 1_000_000)
    overall_lost.append(avg_lost / 1_000_000)

bars_usable = ax5.bar(x, overall_usable, label='Usable k-mers',
                       color='#2ecc71', edgecolor='black')
bars_lost = ax5.bar(x, overall_lost, bottom=overall_usable,
                     label='Lost to errors', color='#e74c3c', edgecolor='black', alpha=0.7)

ax5.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax5.set_ylabel('K-mers per Mb (millions)', fontweight='bold', fontsize=12)
ax5.set_title('E. Coverage Composition: Usable vs Lost K-mers per Megabase',
              fontweight='bold', loc='left', fontsize=13)
ax5.set_xticks(x)
ax5.set_xticklabels([f'k={k}' for k in k_sizes])
ax5.legend(loc='upper left', fontsize=11)
ax5.grid(axis='y', alpha=0.3)

# Add percentage labels
for i, (usable, lost) in enumerate(zip(overall_usable, overall_lost)):
    total = usable + lost
    pct_usable = (usable / total) * 100
    pct_lost = (lost / total) * 100

    # Usable percentage
    ax5.text(i, usable/2, f'{pct_usable:.1f}%',
            ha='center', va='center', fontsize=10, fontweight='bold', color='white')

    # Lost percentage
    ax5.text(i, usable + lost/2, f'{pct_lost:.1f}%',
            ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Panel F: Key insights table
ax6 = fig.add_subplot(gs[2, 2])
ax6.axis('tight')
ax6.axis('off')

insight_text = """
KEY INSIGHTS

What "Coverage Loss" Means:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Each sequencing error affects a
LOCAL window of k-mers:
  • k=21: error affects 21 k-mers
  • k=41: error affects 41 k-mers

You don't lose entire reads!
You lose CLASSIFICATION ABILITY
in that local genomic region.

Why k=21 is Better:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ARMS markers, k=21:
  ~22M lost k-mers/Mb

ARMS markers, k=41:
  ~48M lost k-mers/Mb

k=41 loses 2× MORE coverage
per megabase!

Practical Impact:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
With k=21: ~81% coverage retained
With k=41: ~66% coverage retained

For a 10 Mb region:
  k=21: 8.1 Mb classifiable ✓
  k=41: 6.6 Mb classifiable ✗

You lose 1.5 Mb more with k=41!
"""

ax6.text(0.05, 0.5, insight_text, fontsize=9, family='monospace',
        verticalalignment='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9,
                  edgecolor='orange', linewidth=2))

plt.savefig('final_results/06_coverage_loss_analysis.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/06_coverage_loss_analysis.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/06_coverage_loss_analysis.png")
print(f"✓ Saved: final_results/06_coverage_loss_analysis.pdf")

print("\n" + "="*100)
print("💡 CORRECTED INTERPRETATION:")
print("="*100)
print("You DON'T lose entire reads with longer k-mers!")
print("You lose COVERAGE per megabase because:")
print("  1. Longer k-mers have MORE bases")
print("  2. More bases = higher probability of containing an error")
print("  3. Each error affects k consecutive overlapping k-mers")
print("  4. Those k-mers can't be classified → local coverage loss")
print("")
print("ARMS markers example:")
print(f"  k=21: Lose {arms_lost[0]/1_000_000:.1f}M k-mers per Mb to errors")
print(f"  k=41: Lose {arms_lost[-1]/1_000_000:.1f}M k-mers per Mb to errors")
print(f"  → k=41 loses {(arms_lost[-1]/arms_lost[0]):.1f}× MORE coverage!")
print("="*100)
