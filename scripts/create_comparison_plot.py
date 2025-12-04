#!/usr/bin/env python3
"""
Create a beautiful publication-quality comparison plot for k-mer error resilience.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set publication-quality style
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

# Load data from all k-mer sizes
kmer_sizes = [21, 25, 31, 35, 41]
all_data = []

for k in kmer_sizes:
    df = pd.read_csv(f'error_k{k}_100k_error_resilience_stats.csv')
    df['kmer_size'] = k
    all_data.append(df)

combined_df = pd.concat(all_data, ignore_index=True)

# Calculate summary statistics
summary_stats = []
for k in kmer_sizes:
    k_data = combined_df[combined_df['kmer_size'] == k]

    # Overall
    mean_tolerant = k_data['pct_error_tolerant'].mean()

    # By region
    arms_data = k_data[k_data['region'] == 'ARMS']
    cen_data = k_data[k_data['region'] == 'CEN']

    mean_arms = arms_data['pct_error_tolerant'].mean()
    mean_cen = cen_data['pct_error_tolerant'].mean()

    summary_stats.append({
        'kmer_size': k,
        'overall': mean_tolerant,
        'arms': mean_arms,
        'cen': mean_cen
    })

summary_df = pd.DataFrame(summary_stats)

# Create the figure
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# Color palette
colors = {
    'overall': '#2E86AB',
    'arms': '#A23B72',
    'cen': '#F18F01',
    'Col-0': '#e63946',
    'Ler-0': '#457b9d'
}

# ============================================================================
# Panel A: Error Tolerance by K-mer Size (Bar chart)
# ============================================================================
ax1 = fig.add_subplot(gs[0, :2])

x = np.arange(len(kmer_sizes))
width = 0.25

bars1 = ax1.bar(x - width, summary_df['overall'], width,
                label='Overall', color=colors['overall'], alpha=0.8, edgecolor='black', linewidth=1)
bars2 = ax1.bar(x, summary_df['arms'], width,
                label='ARMS', color=colors['arms'], alpha=0.8, edgecolor='black', linewidth=1)
bars3 = ax1.bar(x + width, summary_df['cen'], width,
                label='CEN', color=colors['cen'], alpha=0.8, edgecolor='black', linewidth=1)

# Add value labels on bars
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        if height > 0.01:
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}%', ha='center', va='bottom', fontsize=8)

ax1.set_xlabel('K-mer Size', fontweight='bold')
ax1.set_ylabel('Error-Tolerant K-mers (%)', fontweight='bold')
ax1.set_title('A. Error Tolerance Comparison Across K-mer Sizes',
              fontweight='bold', pad=15)
ax1.set_xticks(x)
ax1.set_xticklabels([f'k={k}' for k in kmer_sizes])
ax1.legend(frameon=True, fancybox=True, shadow=True)
ax1.set_ylim(0, max(summary_df['cen']) * 1.2)
ax1.grid(axis='y', alpha=0.3, linestyle='--')

# ============================================================================
# Panel B: Expected Read Retention (showing practical impact)
# ============================================================================
ax2 = fig.add_subplot(gs[0, 2])

# Calculate expected read retention (assuming 0.1% sequencing error rate)
error_rate = 0.001
read_retention = []
for k in kmer_sizes:
    # Probability of NO error in k-mer of length k
    p_no_error = (1 - error_rate) ** k
    read_retention.append(p_no_error * 100)

bars = ax2.barh(range(len(kmer_sizes)), read_retention,
                color=sns.color_palette("RdYlGn", len(kmer_sizes))[::-1],
                edgecolor='black', linewidth=1)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, read_retention)):
    reads_lost = 100 - val
    ax2.text(val - 0.1, bar.get_y() + bar.get_height()/2,
            f'{val:.2f}%', ha='right', va='center', fontweight='bold', fontsize=9)
    ax2.text(val + 0.05, bar.get_y() + bar.get_height()/2,
            f'(-{reads_lost:.2f}%)', ha='left', va='center', fontsize=8, style='italic')

ax2.set_yticks(range(len(kmer_sizes)))
ax2.set_yticklabels([f'k={k}' for k in kmer_sizes])
ax2.set_xlabel('Read Retention (%)', fontweight='bold')
ax2.set_title('B. Expected Read Retention\n(0.1% sequencing error)',
              fontweight='bold', pad=15)
ax2.set_xlim(95, 100)
ax2.axvline(98, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax2.grid(axis='x', alpha=0.3, linestyle='--')

# ============================================================================
# Panel C: Distribution by Region (Violin plot)
# ============================================================================
ax3 = fig.add_subplot(gs[1, :])

# Prepare data for violin plot
plot_data = []
for k in kmer_sizes:
    k_data = combined_df[combined_df['kmer_size'] == k]

    for _, row in k_data.iterrows():
        plot_data.append({
            'K-mer Size': f'k={k}',
            'Region': row['region'],
            'Error Tolerance (%)': row['pct_error_tolerant']
        })

plot_df = pd.DataFrame(plot_data)

# Create violin plot
positions = []
for i, k in enumerate(kmer_sizes):
    positions.extend([i*3, i*3 + 1])

parts = ax3.violinplot(
    [plot_df[(plot_df['K-mer Size'] == f'k={k}') & (plot_df['Region'] == 'ARMS')]['Error Tolerance (%)'].values
     for k in kmer_sizes] +
    [plot_df[(plot_df['K-mer Size'] == f'k={k}') & (plot_df['Region'] == 'CEN')]['Error Tolerance (%)'].values
     for k in kmer_sizes],
    positions=[i*3 for i in range(len(kmer_sizes))] + [i*3+1 for i in range(len(kmer_sizes))],
    widths=0.7,
    showmeans=True,
    showmedians=True
)

# Color the violins
for i, pc in enumerate(parts['bodies']):
    if i < len(kmer_sizes):
        pc.set_facecolor(colors['arms'])
        pc.set_alpha(0.7)
    else:
        pc.set_facecolor(colors['cen'])
        pc.set_alpha(0.7)

# Set x-axis
ax3.set_xticks([i*3 + 0.5 for i in range(len(kmer_sizes))])
ax3.set_xticklabels([f'k={k}' for k in kmer_sizes])
ax3.set_xlabel('K-mer Size', fontweight='bold')
ax3.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax3.set_title('C. Error Tolerance Distribution by Region',
              fontweight='bold', pad=15)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=colors['arms'], alpha=0.7, label='ARMS'),
    Patch(facecolor=colors['cen'], alpha=0.7, label='CEN')
]
ax3.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=True, shadow=True)
ax3.grid(axis='y', alpha=0.3, linestyle='--')
ax3.set_ylim(-0.1, max(plot_df['Error Tolerance (%)']) * 1.1)

# ============================================================================
# Add main title and save
# ============================================================================
fig.suptitle('Sequencing Error Resilience Analysis Across K-mer Sizes\n'
             'Impact of Single-Base Sequencing Errors on Marker Specificity',
             fontsize=16, fontweight='bold', y=0.98)

plt.savefig('kmer_error_resilience_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('kmer_error_resilience_comparison.pdf', bbox_inches='tight')
print("✓ Saved publication-quality figures:")
print("  - kmer_error_resilience_comparison.png (300 dpi)")
print("  - kmer_error_resilience_comparison.pdf (vector)")

# ============================================================================
# Create a second figure: Detailed heatmap
# ============================================================================
fig2, ax = plt.subplots(figsize=(12, 8))

# Pivot data for heatmap
heatmap_data = combined_df.pivot_table(
    values='pct_error_tolerant',
    index='database',
    columns='kmer_size'
)

# Sort by genotype and region
sort_order = combined_df.groupby('database').first().sort_values(['genotype', 'region', 'chromosome']).index
heatmap_data = heatmap_data.reindex(sort_order)

# Create heatmap
im = ax.imshow(heatmap_data.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1.5)

# Set ticks and labels
ax.set_xticks(np.arange(len(kmer_sizes)))
ax.set_yticks(np.arange(len(heatmap_data)))
ax.set_xticklabels([f'k={k}' for k in kmer_sizes], fontsize=11)
ax.set_yticklabels(heatmap_data.index, fontsize=9)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Error Tolerance (%)', fontweight='bold', fontsize=11)

# Add text annotations
for i in range(len(heatmap_data)):
    for j in range(len(kmer_sizes)):
        val = heatmap_data.iloc[i, j]
        if not np.isnan(val):
            text_color = 'white' if val > 0.75 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                   color=text_color, fontsize=7, fontweight='bold')

# Add region separators
genotypes = [heatmap_data.index[i].split('_')[0] for i in range(len(heatmap_data))]
regions = [heatmap_data.index[i].split('_')[1] for i in range(len(heatmap_data))]

for i in range(1, len(heatmap_data)):
    if genotypes[i] != genotypes[i-1] or regions[i] != regions[i-1]:
        ax.axhline(i - 0.5, color='white', linewidth=2)

ax.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax.set_ylabel('Marker Database', fontweight='bold', fontsize=12)
ax.set_title('Error Tolerance Heatmap: All Databases × All K-mer Sizes\n'
             'Percentage of k-mers remaining specific after 1 sequencing error',
             fontweight='bold', fontsize=13, pad=15)

plt.tight_layout()
plt.savefig('kmer_error_resilience_heatmap.png', dpi=300, bbox_inches='tight')
plt.savefig('kmer_error_resilience_heatmap.pdf', bbox_inches='tight')
print("✓ Saved heatmap figures:")
print("  - kmer_error_resilience_heatmap.png (300 dpi)")
print("  - kmer_error_resilience_heatmap.pdf (vector)")

print("\n✨ All beautiful plots created successfully!")
