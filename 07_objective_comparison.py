#!/usr/bin/env python3
"""
Plot 7: Objective K-mer Comparison - No Arbitrary Scoring
Present the trade-offs clearly and let the user decide.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

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
error_df['absolute_false_positive_rate'] = (error_df['pct_kmers_with_errors'] / 100) * (error_df['pct_wrong_db'] / 100) * 100

# Load marker availability
marker_df = pd.read_csv("final_results/marker_availability_summary.csv")
error_df = error_df.merge(marker_df[['k_size', 'database', 'total_kmers', 'density_per_Mb']],
                           on=['k_size', 'database'], how='left')

print(f"✓ Loaded data for {len(error_df)} databases across {len(k_sizes)} k-mer sizes")

# Calculate key metrics
summary = []
for k in k_sizes:
    k_data = error_df[error_df['k_size'] == k]

    # ARMS markers (more important for analysis)
    arms_data = k_data[k_data['region'] == 'ARMS']

    avg_density = arms_data['density_per_Mb'].mean()
    avg_error_rate = arms_data['pct_kmers_with_errors'].mean()
    usable_kmers_per_mb = avg_density * (100 - avg_error_rate) / 100
    avg_fp_rate = arms_data['absolute_false_positive_rate'].mean()

    summary.append({
        'k_size': k,
        'total_density': avg_density / 1e6,  # millions
        'error_rate': avg_error_rate,
        'usable_per_mb': usable_kmers_per_mb / 1e6,  # millions
        'fp_rate': avg_fp_rate,
        'total_markers': arms_data['total_kmers'].mean() / 1e6  # millions
    })

summary_df = pd.DataFrame(summary)

print("\n" + "="*100)
print("OBJECTIVE K-MER COMPARISON (ARMS markers)")
print("="*100)
print(summary_df.to_string(index=False))
print("="*100)

# Create visualization
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.35)

fig.suptitle('Objective K-mer Size Comparison: Understanding the Trade-offs\n'
             'No arbitrary scoring - Make your own informed decision',
             fontsize=17, fontweight='bold', y=0.98)

# Define colors - all neutral, no highlighting
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# Panel A: Usable k-mers per Mb (THE KEY METRIC)
ax1 = fig.add_subplot(gs[0, :2])

bars = ax1.bar(range(len(k_sizes)), summary_df['usable_per_mb'],
               color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

ax1.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
ax1.set_ylabel('Usable K-mers per Mb (millions)', fontweight='bold', fontsize=13)
ax1.set_title('A. Classification Ability: Usable K-mers per Megabase',
              fontweight='bold', loc='left', fontsize=14)
ax1.set_xticks(range(len(k_sizes)))
ax1.set_xticklabels([f'k={k}' for k in k_sizes])
ax1.grid(axis='y', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, summary_df['usable_per_mb'])):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
            f'{val:.1f}M',
            ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Show difference from k=21
    if i > 0:
        diff = ((val - summary_df['usable_per_mb'].iloc[0]) / summary_df['usable_per_mb'].iloc[0]) * 100
        ax1.text(bar.get_x() + bar.get_width()/2., height - 8,
                f'+{diff:.1f}%',
                ha='center', va='top', fontsize=10, color='darkgreen', fontweight='bold')

ax1.text(0.02, 0.97, 'Higher = More classification attempts possible',
         transform=ax1.transAxes, fontsize=11, style='italic',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

# Panel B: False positive rate
ax2 = fig.add_subplot(gs[0, 2:])

bars = ax2.bar(range(len(k_sizes)), summary_df['fp_rate'],
               color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

ax2.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
ax2.set_ylabel('False Positive Rate (%)', fontweight='bold', fontsize=13)
ax2.set_title('B. Specificity: False Positive Rate',
              fontweight='bold', loc='left', fontsize=14)
ax2.set_xticks(range(len(k_sizes)))
ax2.set_xticklabels([f'k={k}' for k in k_sizes])
ax2.grid(axis='y', alpha=0.3)

# Add value labels
for bar, val in zip(bars, summary_df['fp_rate']):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.005,
            f'{val:.3f}%',
            ha='center', va='bottom', fontsize=11, fontweight='bold')

ax2.axhline(y=0.2, color='orange', linestyle='--', alpha=0.5, linewidth=2)
ax2.text(len(k_sizes)-0.5, 0.205, 'All excellent (<0.2%)',
        ha='right', fontsize=10, color='orange', fontweight='bold')

ax2.text(0.02, 0.97, 'Lower = Better specificity (fewer false positives)',
         transform=ax2.transAxes, fontsize=11, style='italic',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

# Panel C: Error impact
ax3 = fig.add_subplot(gs[1, 0])

line = ax3.plot(k_sizes, summary_df['error_rate'], 'o-',
                linewidth=3, markersize=12, color='#e74c3c')

ax3.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax3.set_ylabel('K-mers Affected by Errors (%)', fontweight='bold', fontsize=12)
ax3.set_title('C. Error Resilience', fontweight='bold', loc='left', fontsize=13)
ax3.grid(alpha=0.3)

# Add value labels
for k, val in zip(k_sizes, summary_df['error_rate']):
    ax3.text(k, val + 1, f'{val:.1f}%', ha='center', fontsize=10, fontweight='bold')

ax3.text(0.5, 0.05, 'Longer k-mers more affected',
         transform=ax3.transAxes, ha='center', fontsize=10, style='italic',
         bbox=dict(boxstyle='round', facecolor='#ffcccc', alpha=0.5))

# Panel D: Total marker count
ax4 = fig.add_subplot(gs[1, 1])

line = ax4.plot(k_sizes, summary_df['total_markers'], 'o-',
                linewidth=3, markersize=12, color='#2ca02c')

ax4.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax4.set_ylabel('Average Markers per Database (M)', fontweight='bold', fontsize=12)
ax4.set_title('D. Marker Availability', fontweight='bold', loc='left', fontsize=13)
ax4.grid(alpha=0.3)

# Add value labels
for k, val in zip(k_sizes, summary_df['total_markers']):
    ax4.text(k, val + 0.1, f'{val:.1f}M', ha='center', fontsize=10, fontweight='bold')

ax4.text(0.5, 0.05, 'All have sufficient markers',
         transform=ax4.transAxes, ha='center', fontsize=10, style='italic',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

# Panel E: Density breakdown
ax5 = fig.add_subplot(gs[1, 2:])

x = np.arange(len(k_sizes))
usable = summary_df['usable_per_mb'].values
lost = summary_df['total_density'].values - usable

bars1 = ax5.bar(x, usable, label='Usable', color='#2ecc71', edgecolor='black')
bars2 = ax5.bar(x, lost, bottom=usable, label='Lost to errors',
                color='#e74c3c', edgecolor='black', alpha=0.7)

ax5.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax5.set_ylabel('K-mers per Mb (millions)', fontweight='bold', fontsize=12)
ax5.set_title('E. Density Breakdown: Usable vs Lost', fontweight='bold', loc='left', fontsize=13)
ax5.set_xticks(x)
ax5.set_xticklabels([f'k={k}' for k in k_sizes])
ax5.legend(fontsize=11)
ax5.grid(axis='y', alpha=0.3)

# Add percentage labels
for i, (u, l) in enumerate(zip(usable, lost)):
    total = u + l
    pct_usable = (u / total) * 100
    ax5.text(i, u/2, f'{pct_usable:.0f}%', ha='center', va='center',
            fontsize=10, fontweight='bold', color='white')

# Panel F: Trade-off summary table
ax6 = fig.add_subplot(gs[2, :])
ax6.axis('tight')
ax6.axis('off')

# Create comprehensive comparison table
table_data = [
    ['K-mer', 'Usable K-mers/Mb', 'False Positive', 'Error Impact', 'Markers', 'Best For'],
    ['', '(Higher = Better)', '(Lower = Better)', '(Lower = Better)', '', ''],
]

use_cases = [
    'Very noisy data (>2% errors)',
    'Balanced general use',
    'Good balance + coverage',
    'High coverage needs',
    'Maximum resolution/specificity'
]

for i, (_, row) in enumerate(summary_df.iterrows()):
    table_data.append([
        f"k={int(row['k_size'])}",
        f"{row['usable_per_mb']:.1f}M",
        f"{row['fp_rate']:.3f}%",
        f"{row['error_rate']:.1f}%",
        f"{row['total_markers']:.1f}M",
        use_cases[i]
    ])

table = ax6.table(cellText=table_data, cellLoc='center', loc='center',
                 bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.8)

# Style header rows
for i in range(6):
    cell = table[(0, i)]
    cell.set_facecolor('#34495e')
    cell.set_text_props(weight='bold', color='white', size=11)

    cell2 = table[(1, i)]
    cell2.set_facecolor('#7f8c8d')
    cell2.set_text_props(style='italic', color='white', size=9)

# Style data rows with alternating colors
for i in range(2, len(table_data)):
    row_color = '#f8f9fa' if i % 2 == 0 else 'white'
    for j in range(6):
        cell = table[(i, j)]
        cell.set_facecolor(row_color)

        # Bold the k-mer column
        if j == 0:
            cell.set_text_props(weight='bold', size=11)

# Highlight metric winners in their respective columns
# Usable k-mers - highest is k=41
table[(6, 1)].set_text_props(color='darkgreen', weight='bold')
# False positive - lowest is k=41
table[(6, 2)].set_text_props(color='darkgreen', weight='bold')
# Error impact - lowest is k=21
table[(2, 3)].set_text_props(color='darkgreen', weight='bold')

ax6.set_title('F. Comprehensive Comparison: Understanding Your Options',
              fontweight='bold', fontsize=14, pad=15)

# Add decision guide
ax7 = fig.add_subplot(gs[2, :])
ax7.axis('off')

decision_text = """
DECISION GUIDE - No "one size fits all" answer!

Your choice depends on your priorities:

Prioritize ERROR RESILIENCE → k=21 or k=25
  • Noisy ONT data (>2% per-base errors)
  • Minimal false positives acceptable
  • Can accept lower coverage

Prioritize BALANCE → k=31 ⭐ GENERAL PURPOSE
  • Standard ONT quality (~1% errors)
  • Good coverage + reasonable specificity
  • Most versatile choice

Prioritize MAXIMUM COVERAGE → k=35 or k=41
  • High-quality ONT data (<1% errors)
  • Need maximum resolution
  • Best specificity possible
  • Can tolerate more error impact
"""

ax7.text(0.5, -1.3, decision_text, fontsize=10, family='monospace',
        ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9,
                  edgecolor='orange', linewidth=2))

plt.savefig('final_results/07_objective_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/07_objective_comparison.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/07_objective_comparison.png")
print(f"✓ Saved: final_results/07_objective_comparison.pdf")

print("\n" + "="*100)
print("OBJECTIVE SUMMARY - THE TRADE-OFFS")
print("="*100)
print("There is NO universally 'best' k-mer size!")
print("")
print("The choice depends on YOUR priorities:")
print("  • k=21: Best error resilience, lowest coverage")
print("  • k=31: Balanced - good coverage, good specificity, moderate errors")
print("  • k=41: Maximum coverage & specificity, most affected by errors")
print("")
print("All choices are valid depending on your data quality and needs.")
print("="*100)
