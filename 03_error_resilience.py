#!/usr/bin/env python3
"""
Plot 3: Error Resilience and Coverage Retention (OBJECTIVE VERSION)
Shows coverage retention across k-mer sizes without bias.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')

# Define neutral colors for all k-mer sizes
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# Load error resilience data
k_sizes = [21, 25, 31, 35, 41]
all_data = []

for k in k_sizes:
    csv_file = Path(f"final_results/realistic_k{k}_100k_error_resilience_stats.csv")
    if csv_file.exists():
        df = pd.read_csv(csv_file)
        df['k_size'] = k
        all_data.append(df)

if not all_data:
    print("ERROR: No error resilience data found!")
    exit(1)

df = pd.concat(all_data, ignore_index=True)

# Calculate overall statistics per k-size
overall_stats = []
for k in k_sizes:
    k_data = df[df['k_size'] == k]

    pct_with_errors = k_data['pct_kmers_with_errors'].mean()
    pct_error_tolerant = k_data['pct_error_tolerant'].mean()
    usable_kmers_pct = (100 - pct_with_errors) + (pct_with_errors * pct_error_tolerant / 100)

    overall_stats.append({
        'k_size': k,
        'pct_kmers_with_errors': pct_with_errors,
        'pct_becomes_novel': k_data['pct_becomes_novel'].mean(),
        'pct_error_tolerant': pct_error_tolerant,
        'usable_kmers_pct': usable_kmers_pct
    })

overall_df = pd.DataFrame(overall_stats)

print(f"✓ Loaded data for {len(df)} databases across {len(k_sizes)} k-mer sizes")

# Create visualization
fig = plt.figure(figsize=(18, 11))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

fig.suptitle('Error Resilience and Coverage Retention Analysis\n(1% per-base error rate, ONT-like sequencing)',
             fontsize=16, fontweight='bold', y=0.98)

# Panel A: Coverage retention (THE KEY METRIC)
ax1 = fig.add_subplot(gs[0, :])

bars = ax1.bar(range(len(k_sizes)), overall_df['usable_kmers_pct'],
               color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

ax1.set_xlabel('K-mer Size', fontweight='bold', fontsize=13)
ax1.set_ylabel('Usable K-mers (%)', fontweight='bold', fontsize=13)
ax1.set_title('A. Coverage Retention: Percentage of K-mers Remaining Classifiable',
              fontweight='bold', loc='left', fontsize=14)
ax1.set_xticks(range(len(k_sizes)))
ax1.set_xticklabels([f'k={k}' for k in k_sizes])
ax1.set_ylim(60, 85)
ax1.grid(axis='y', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, overall_df['usable_kmers_pct'])):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{val:.2f}%',
            ha='center', va='bottom', fontsize=11, fontweight='bold')

# Add trend annotation
ax1.text(0.5, 0.05, 'Trade-off: Longer k-mers have lower retention but higher total density',
         transform=ax1.transAxes, ha='center', fontsize=11, style='italic',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.6))

# Panel B: K-mers affected by errors
ax2 = fig.add_subplot(gs[1, 0])

bars = ax2.bar(range(len(k_sizes)), overall_df['pct_kmers_with_errors'],
               color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)

ax2.set_xlabel('K-mer Size', fontweight='bold')
ax2.set_ylabel('K-mers With Errors (%)', fontweight='bold')
ax2.set_title('B. Error Impact', fontweight='bold', loc='left')
ax2.set_xticks(range(len(k_sizes)))
ax2.set_xticklabels([f'k={k}' for k in k_sizes])
ax2.grid(axis='y', alpha=0.3)

for bar, val in zip(bars, overall_df['pct_kmers_with_errors']):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{val:.1f}%',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# Panel C: Error tolerance by region
ax3 = fig.add_subplot(gs[1, 1])

summary = df.groupby(['k_size', 'region'])['pct_error_tolerant'].agg(['mean', 'std']).reset_index()
x = np.arange(len(k_sizes))
width = 0.35

arms_data = summary[summary['region'] == 'ARMS'].sort_values('k_size')
cen_data = summary[summary['region'] == 'CEN'].sort_values('k_size')

ax3.bar(x - width/2, arms_data['mean'], width, label='ARMS',
        color='#66c2a5', yerr=arms_data['std'], capsize=3, edgecolor='black')
ax3.bar(x + width/2, cen_data['mean'], width, label='CEN',
        color='#fc8d62', yerr=cen_data['std'], capsize=3, edgecolor='black')

ax3.set_xlabel('K-mer Size', fontweight='bold')
ax3.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax3.set_title('C. Error Tolerance by Region', fontweight='bold', loc='left')
ax3.set_xticks(x)
ax3.set_xticklabels([f'k={k}' for k in k_sizes])
ax3.legend()
ax3.grid(axis='y', alpha=0.3)

# Panel D: Violin plot - error tolerance distribution
ax4 = fig.add_subplot(gs[1, 2])

plot_data = df[df['k_size'].isin([21, 31, 41])].copy()

import matplotlib.patches as mpatches
arms_color = '#66c2a5'
cen_color = '#fc8d62'

for i, k in enumerate([21, 31, 41]):
    for j, region in enumerate(['ARMS', 'CEN']):
        subset = plot_data[(plot_data['k_size'] == k) & (plot_data['region'] == region)]
        if not subset.empty:
            pos = i * 2.5 + j * 1
            color = arms_color if region == 'ARMS' else cen_color
            parts = ax4.violinplot([subset['pct_error_tolerant'].values],
                                   positions=[pos],
                                   widths=0.7,
                                   showmeans=True,
                                   showmedians=True)
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)

ax4.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax4.set_title('D. Error Tolerance Distribution', fontweight='bold', loc='left')
ax4.set_xticks([0.5, 2.5*1 + 0.5, 2.5*2 + 0.5])
ax4.set_xticklabels(['k=21', 'k=31', 'k=41'])
ax4.grid(axis='y', alpha=0.3)

arms_patch = mpatches.Patch(color=arms_color, label='ARMS', alpha=0.7)
cen_patch = mpatches.Patch(color=cen_color, label='CEN', alpha=0.7)
ax4.legend(handles=[arms_patch, cen_patch], loc='upper right')

# Panel E: Error outcomes comparison
ax5 = fig.add_subplot(gs[2, :2])

x_pos = np.arange(len(k_sizes))
width = 0.25

novel_pct = overall_df['pct_becomes_novel'].values
tolerant_pct = overall_df['pct_error_tolerant'].values
# Calculate wrong_db as remainder (simplification)
wrong_pct = 100 - novel_pct - tolerant_pct

bars1 = ax5.bar(x_pos - width, novel_pct, width, label='Becomes Novel (lost)',
                color='#95a5a6', edgecolor='black')
bars2 = ax5.bar(x_pos, tolerant_pct, width, label='Stays Correct',
                color='#2ecc71', edgecolor='black')
bars3 = ax5.bar(x_pos + width, wrong_pct, width, label='Wrong DB (false pos)',
                color='#e74c3c', edgecolor='black', alpha=0.8)

ax5.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax5.set_ylabel('Percentage of Errors (%)', fontweight='bold', fontsize=12)
ax5.set_title('E. Error Outcome Distribution', fontweight='bold', loc='left', fontsize=13)
ax5.set_xticks(x_pos)
ax5.set_xticklabels([f'k={k}' for k in k_sizes])
ax5.legend(loc='upper right', fontsize=10)
ax5.grid(axis='y', alpha=0.3)
ax5.set_ylim(0, 105)

# Panel F: Summary table
ax6 = fig.add_subplot(gs[2, 2])
ax6.axis('tight')
ax6.axis('off')

table_data = [['K-mer', 'Retention', 'Error Impact', 'Interpretation']]

interpretations = [
    'Best resilience',
    'Good balance',
    'Balanced',
    'More affected',
    'Most affected'
]

for i, (_, row) in enumerate(overall_df.iterrows()):
    table_data.append([
        f"k={int(row['k_size'])}",
        f"{row['usable_kmers_pct']:.1f}%",
        f"{row['pct_kmers_with_errors']:.1f}%",
        interpretations[i]
    ])

table = ax6.table(cellText=table_data, cellLoc='center', loc='center',
                 bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 3)

# Style header
for i in range(4):
    cell = table[(0, i)]
    cell.set_facecolor('#34495e')
    cell.set_text_props(weight='bold', color='white')

# Style data rows
for i in range(1, len(table_data)):
    row_color = '#f8f9fa' if i % 2 == 0 else 'white'
    for j in range(4):
        table[(i, j)].set_facecolor(row_color)

ax6.set_title('F. Summary Table', fontweight='bold', fontsize=12, pad=10)

plt.savefig('final_results/03_error_resilience.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/03_error_resilience.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/03_error_resilience.png")
print(f"✓ Saved: final_results/03_error_resilience.pdf")

# Print summary
print("\n" + "="*80)
print("ERROR RESILIENCE SUMMARY")
print("="*80)
print(f"{'K-mer':<8} {'Retention':<15} {'Error Impact':<15} {'Balance':<20}")
print("-"*80)
for _, row in overall_df.iterrows():
    usable = row['usable_kmers_pct']
    errors = row['pct_kmers_with_errors']
    print(f"k={row['k_size']:<5} {usable:>6.2f}%         {errors:>6.2f}%")
print("="*80)
print("\n💡 All k-mer sizes are viable - choice depends on your priorities")
print("   Higher retention = better for noisy data")
print("   But remember: longer k-mers have higher total density!")
print("="*80)
