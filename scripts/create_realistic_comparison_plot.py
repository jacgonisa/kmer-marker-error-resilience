#!/usr/bin/env python3
"""
Create publication-quality comparison plot for REALISTIC k-mer error resilience.
Uses per-base 1% error rate (ONT-like sequencing).
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
    df = pd.read_csv(f'realistic_k{k}_100k_error_resilience_stats.csv')
    df['kmer_size'] = k
    all_data.append(df)

combined_df = pd.concat(all_data, ignore_index=True)

# Calculate summary statistics
summary_stats = []
for k in kmer_sizes:
    k_data = combined_df[combined_df['kmer_size'] == k]

    # Overall
    mean_tolerant = k_data['pct_error_tolerant'].mean()
    mean_affected = k_data['pct_kmers_with_errors'].mean()

    # By region
    arms_data = k_data[k_data['region'] == 'ARMS']
    cen_data = k_data[k_data['region'] == 'CEN']

    mean_arms_tolerant = arms_data['pct_error_tolerant'].mean()
    mean_cen_tolerant = cen_data['pct_error_tolerant'].mean()

    mean_arms_affected = arms_data['pct_kmers_with_errors'].mean()
    mean_cen_affected = cen_data['pct_kmers_with_errors'].mean()

    summary_stats.append({
        'kmer_size': k,
        'overall_tolerant': mean_tolerant,
        'arms_tolerant': mean_arms_tolerant,
        'cen_tolerant': mean_cen_tolerant,
        'overall_affected': mean_affected,
        'arms_affected': mean_arms_affected,
        'cen_affected': mean_cen_affected
    })

summary_df = pd.DataFrame(summary_stats)

# Create the figure
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

# Color palette
colors = {
    'overall': '#2E86AB',
    'arms': '#A23B72',
    'cen': '#F18F01',
    'Col-0': '#e63946',
    'Ler-0': '#457b9d'
}

# ============================================================================
# Panel A: K-mers Affected by Errors (THE KEY DIFFERENCE!)
# ============================================================================
ax1 = fig.add_subplot(gs[0, :2])

x = np.arange(len(kmer_sizes))
width = 0.35

bars1 = ax1.bar(x - width/2, summary_df['arms_affected'], width,
                label='ARMS', color=colors['arms'], alpha=0.8, edgecolor='black', linewidth=1.5)
bars2 = ax1.bar(x + width/2, summary_df['cen_affected'], width,
                label='CEN', color=colors['cen'], alpha=0.8, edgecolor='black', linewidth=1.5)

# Add value labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

ax1.set_xlabel('K-mer Size', fontweight='bold')
ax1.set_ylabel('K-mers with ≥1 Sequencing Error (%)', fontweight='bold')
ax1.set_title('A. Sequencing Error Impact: K-mers Affected\n(1% per-base error rate, ONT-like)',
              fontweight='bold', pad=15)
ax1.set_xticks(x)
ax1.set_xticklabels([f'k={k}' for k in kmer_sizes])
ax1.legend(frameon=True, fancybox=True, shadow=True, loc='upper left')
ax1.set_ylim(0, 40)
ax1.grid(axis='y', alpha=0.3, linestyle='--')

# Add annotation showing the problem
ax1.annotate('', xy=(4, 33.7), xytext=(0, 19),
            arrowprops=dict(arrowstyle='->', lw=2, color='red'))
ax1.text(2, 28, '+76% more\nk-mers lost!', fontsize=11, color='red',
        fontweight='bold', ha='center',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))

# ============================================================================
# Panel B: Error Tolerance of Affected K-mers
# ============================================================================
ax2 = fig.add_subplot(gs[0, 2])

bars1 = ax2.barh(range(len(kmer_sizes)), summary_df['arms_tolerant'],
                color=colors['arms'], alpha=0.8, edgecolor='black', linewidth=1.5)
bars2 = ax2.barh([i+0.4 for i in range(len(kmer_sizes))], summary_df['cen_tolerant'],
                color=colors['cen'], alpha=0.8, edgecolor='black', linewidth=1.5, height=0.35)

# Add value labels
for i, val in enumerate(summary_df['arms_tolerant']):
    ax2.text(val + 0.01, i, f'{val:.2f}%', va='center', fontsize=8, fontweight='bold')
for i, val in enumerate(summary_df['cen_tolerant']):
    ax2.text(val + 0.01, i+0.4, f'{val:.2f}%', va='center', fontsize=8, fontweight='bold')

ax2.set_yticks([i+0.2 for i in range(len(kmer_sizes))])
ax2.set_yticklabels([f'k={k}' for k in kmer_sizes])
ax2.set_xlabel('Error-Tolerant (%)', fontweight='bold')
ax2.set_title('B. Error Tolerance\n(of k-mers WITH errors)',
              fontweight='bold', pad=15)
ax2.set_xlim(0, 1.2)
ax2.grid(axis='x', alpha=0.3, linestyle='--')
ax2.legend([bars1, bars2], ['ARMS', 'CEN'], loc='upper right',
          frameon=True, fancybox=True, shadow=True)

# ============================================================================
# Panel C: Read Retention Comparison (PRACTICAL IMPACT)
# ============================================================================
ax3 = fig.add_subplot(gs[1, :])

# Calculate expected usable reads
# (100% - % with errors) + (% with errors × % error-tolerant)
usable_reads = []
for i, row in summary_df.iterrows():
    # Overall usable = no errors + (has errors but still maps correctly)
    pct_no_errors = 100 - row['overall_affected']
    pct_has_errors_but_ok = row['overall_affected'] * (row['overall_tolerant'] / 100)
    pct_usable = pct_no_errors + pct_has_errors_but_ok
    usable_reads.append(pct_usable)

bars = ax3.bar(range(len(kmer_sizes)), usable_reads,
              color=sns.color_palette("RdYlGn_r", len(kmer_sizes)),
              edgecolor='black', linewidth=1.5, alpha=0.8)

# Add value labels and loss indicators
baseline = usable_reads[0]  # k=21
for i, (bar, val) in enumerate(zip(bars, usable_reads)):
    loss = baseline - val
    ax3.text(bar.get_x() + bar.get_width()/2, val + 0.2,
            f'{val:.2f}%', ha='center', va='bottom', fontweight='bold', fontsize=11)
    if i > 0:
        ax3.text(bar.get_x() + bar.get_width()/2, val - 1.5,
                f'-{loss:.2f}%', ha='center', va='top', fontsize=9,
                style='italic', color='darkred', fontweight='bold')

ax3.set_xticks(range(len(kmer_sizes)))
ax3.set_xticklabels([f'k={k}' for k in kmer_sizes])
ax3.set_xlabel('K-mer Size', fontweight='bold')
ax3.set_ylabel('Usable Reads (%)', fontweight='bold')
ax3.set_title('C. Expected Read Retention After Sequencing Errors\n'\
              'Percentage of reads that remain correctly classified',
              fontweight='bold', pad=15)
ax3.set_ylim(78, 82)
ax3.axhline(80, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax3.grid(axis='y', alpha=0.3, linestyle='--')

# ============================================================================
# Panel D: Error Tolerance Distribution by Region (Violin plot)
# ============================================================================
ax4 = fig.add_subplot(gs[2, :])

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
positions_arms = [i*3 for i in range(len(kmer_sizes))]
positions_cen = [i*3 + 1 for i in range(len(kmer_sizes))]

arms_data = [plot_df[(plot_df['K-mer Size'] == f'k={k}') & (plot_df['Region'] == 'ARMS')]['Error Tolerance (%)'].values
             for k in kmer_sizes]
cen_data = [plot_df[(plot_df['K-mer Size'] == f'k={k}') & (plot_df['Region'] == 'CEN')]['Error Tolerance (%)'].values
            for k in kmer_sizes]

parts1 = ax4.violinplot(arms_data, positions=positions_arms, widths=0.7,
                       showmeans=True, showmedians=True)
parts2 = ax4.violinplot(cen_data, positions=positions_cen, widths=0.7,
                       showmeans=True, showmedians=True)

# Color the violins
for pc in parts1['bodies']:
    pc.set_facecolor(colors['arms'])
    pc.set_alpha(0.7)
for pc in parts2['bodies']:
    pc.set_facecolor(colors['cen'])
    pc.set_alpha(0.7)

# Set x-axis
ax4.set_xticks([i*3 + 0.5 for i in range(len(kmer_sizes))])
ax4.set_xticklabels([f'k={k}' for k in kmer_sizes])
ax4.set_xlabel('K-mer Size', fontweight='bold')
ax4.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax4.set_title('D. Error Tolerance Distribution by Region\n'\
              '(Only k-mers with sequencing errors are counted)',
              fontweight='bold', pad=15)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=colors['arms'], alpha=0.7, label='ARMS'),
    Patch(facecolor=colors['cen'], alpha=0.7, label='CEN')
]
ax4.legend(handles=legend_elements, loc='upper right', frameon=True, fancybox=True, shadow=True)
ax4.grid(axis='y', alpha=0.3, linestyle='--')
ax4.set_ylim(-0.05, max(plot_df['Error Tolerance (%)']) * 1.1)

# ============================================================================
# Add main title and save
# ============================================================================
fig.suptitle('Realistic Sequencing Error Resilience Analysis (1% per-base error rate)\n'\
             'Impact of ONT-like Sequencing Errors on Cenhapmer Marker Performance',
             fontsize=16, fontweight='bold', y=0.99)

plt.savefig('realistic_kmer_error_resilience_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('realistic_kmer_error_resilience_comparison.pdf', bbox_inches='tight')
print("✓ Saved publication-quality figures:")
print("  - realistic_kmer_error_resilience_comparison.png (300 dpi)")
print("  - realistic_kmer_error_resilience_comparison.pdf (vector)")

# ============================================================================
# Create a summary table
# ============================================================================
print("\n" + "="*80)
print("REALISTIC ERROR RESILIENCE SUMMARY")
print("="*80)
print("\nK-mer Size | K-mers Affected | ARMS Tolerant | CEN Tolerant | Usable Reads")
print("-----------|-----------------|--------------|--------------|--------------")
for i, row in summary_df.iterrows():
    print(f"K={int(row['kmer_size']):2d}       | {row['overall_affected']:5.2f}%         | "
          f"{row['arms_tolerant']:6.3f}%      | {row['cen_tolerant']:6.3f}%     | {usable_reads[i]:5.2f}%")

print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)
print(f"1. Longer k-mers lose MORE reads to errors:")
print(f"   - K=21: {summary_df.loc[0, 'overall_affected']:.1f}% of k-mers affected")
print(f"   - K=41: {summary_df.loc[4, 'overall_affected']:.1f}% of k-mers affected "
      f"(+{summary_df.loc[4, 'overall_affected']-summary_df.loc[0, 'overall_affected']:.1f}% more!)")
print(f"\n2. Expected usable reads:")
print(f"   - K=21: {usable_reads[0]:.2f}% of reads remain correctly classified")
print(f"   - K=41: {usable_reads[4]:.2f}% of reads remain correctly classified "
      f"({usable_reads[0]-usable_reads[4]:.2f}% loss vs K=21)")
print(f"\n3. ARMS markers show low but NON-ZERO error tolerance:")
print(f"   - This is realistic! Unique sequences occasionally have error-tolerant positions")
print(f"   - CEN markers are 10-20× more error-tolerant due to repetitive nature")
print(f"\n4. RECOMMENDATION: Use K=21 for optimal read retention under ONT sequencing")
print("="*80)

print("\n✨ Analysis complete!")
