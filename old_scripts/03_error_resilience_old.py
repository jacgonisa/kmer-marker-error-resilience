#!/usr/bin/env python3
"""
Plot 3: Error Resilience and Read Retention
Shows how many reads remain usable after ONT sequencing errors.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("Set2")

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

# Parse database names
df['genotype'] = df['database'].str.split('_').str[0]
df['region'] = df['database'].str.split('_').str[1]

# Calculate overall statistics per k-size
overall_stats = []
for k in k_sizes:
    k_data = df[df['k_size'] == k]

    # Overall metrics
    pct_with_errors = k_data['pct_kmers_with_errors'].mean()
    pct_becomes_novel = k_data['pct_becomes_novel'].mean()

    # Calculate usable reads: k-mers without errors + k-mers with errors that stay correct
    # For simplicity: usable = 100% - pct_with_errors + (pct_with_errors * pct_error_tolerant/100)
    pct_error_tolerant = k_data['pct_error_tolerant'].mean()
    usable_reads = (100 - pct_with_errors) + (pct_with_errors * pct_error_tolerant / 100)

    overall_stats.append({
        'k_size': k,
        'pct_kmers_with_errors': pct_with_errors,
        'pct_becomes_novel': pct_becomes_novel,
        'pct_error_tolerant': pct_error_tolerant,
        'usable_reads': usable_reads
    })

overall_df = pd.DataFrame(overall_stats)

print(f"✓ Loaded data for {len(df)} databases across {len(k_sizes)} k-mer sizes")

# Create visualization
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

fig.suptitle('Error Resilience and Read Retention Analysis\n(1% per-base error rate, ONT-like sequencing)',
             fontsize=16, fontweight='bold', y=0.98)

# Panel A: Read retention (MOST IMPORTANT - make it big!)
ax1 = fig.add_subplot(gs[0, :2])

bars = ax1.bar(range(len(k_sizes)), overall_df['usable_reads'],
               color=['#2ecc71' if i == 0 else '#95a5a6' for i in range(len(k_sizes))],
               edgecolor='black', linewidth=1.5)

ax1.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax1.set_ylabel('Usable Reads (%)', fontweight='bold', fontsize=12)
ax1.set_title('A. Read Retention After Sequencing Errors ★ KEY METRIC',
              fontweight='bold', loc='left', fontsize=13, color='darkred')
ax1.set_xticks(range(len(k_sizes)))
ax1.set_xticklabels([f'k={k}' for k in k_sizes])
ax1.set_ylim(60, 85)
ax1.grid(axis='y', alpha=0.3)

# Highlight k=21 as best
bars[0].set_edgecolor('green')
bars[0].set_linewidth(3)

# Add value labels and loss percentages
baseline = overall_df.iloc[0]['usable_reads']
for i, (bar, val) in enumerate(zip(bars, overall_df['usable_reads'])):
    height = bar.get_height()
    # Value label
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{val:.2f}%',
            ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Loss percentage
    if i > 0:
        loss = val - baseline
        ax1.text(bar.get_x() + bar.get_width()/2., height - 2,
                f'{loss:.2f}%',
                ha='center', va='top', fontsize=10, color='red',
                fontweight='bold', style='italic')

# Add annotation
ax1.annotate('← BEST: Highest read retention',
            xy=(0, baseline), xytext=(0.5, baseline - 3),
            arrowprops=dict(arrowstyle='->', color='green', lw=2),
            fontsize=11, color='green', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

# Panel B: K-mers affected by errors
ax2 = fig.add_subplot(gs[0, 2])

bars = ax2.bar(range(len(k_sizes)), overall_df['pct_kmers_with_errors'],
               color='#e74c3c', edgecolor='black', linewidth=1)

ax2.set_xlabel('K-mer Size', fontweight='bold')
ax2.set_ylabel('K-mers With Errors (%)', fontweight='bold')
ax2.set_title('B. Error Impact', fontweight='bold', loc='left')
ax2.set_xticks(range(len(k_sizes)))
ax2.set_xticklabels([f'k={k}' for k in k_sizes], rotation=45)
ax2.grid(axis='y', alpha=0.3)

for bar in bars:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.1f}%',
            ha='center', va='bottom', fontsize=9)

# Add trend arrow
ax2.annotate('', xy=(len(k_sizes)-0.5, overall_df.iloc[-1]['pct_kmers_with_errors']),
            xytext=(0.5, overall_df.iloc[0]['pct_kmers_with_errors']),
            arrowprops=dict(arrowstyle='->', color='red', lw=2, linestyle='--'))
ax2.text(len(k_sizes)/2, 25, 'More errors →',
        ha='center', fontsize=9, color='red', fontweight='bold')

# Panel C: Error tolerance by region
ax3 = fig.add_subplot(gs[1, 0])

summary = df.groupby(['k_size', 'region'])['pct_error_tolerant'].agg(['mean', 'std']).reset_index()
x = np.arange(len(k_sizes))
width = 0.35

arms_data = summary[summary['region'] == 'ARMS'].sort_values('k_size')
cen_data = summary[summary['region'] == 'CEN'].sort_values('k_size')

bars1 = ax3.bar(x - width/2, arms_data['mean'], width, label='ARMS',
                color='#3498db', yerr=arms_data['std'], capsize=3)
bars2 = ax3.bar(x + width/2, cen_data['mean'], width, label='CEN',
                color='#e74c3c', yerr=cen_data['std'], capsize=3)

ax3.set_xlabel('K-mer Size', fontweight='bold')
ax3.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax3.set_title('C. Error Tolerance by Region', fontweight='bold', loc='left')
ax3.set_xticks(x)
ax3.set_xticklabels([f'k={k}' for k in k_sizes])
ax3.legend()
ax3.grid(axis='y', alpha=0.3)

# Panel D: Violin plot - error tolerance distribution
ax4 = fig.add_subplot(gs[1, 1])

plot_data = df[df['k_size'].isin([21, 31, 41])].copy()  # Show subset for clarity
plot_data['k_region'] = plot_data['k_size'].astype(str) + '\n' + plot_data['region']

import matplotlib.patches as mpatches
arms_color = '#3498db'
cen_color = '#e74c3c'

positions = []
labels = []
for i, k in enumerate([21, 31, 41]):
    for j, region in enumerate(['ARMS', 'CEN']):
        subset = plot_data[(plot_data['k_size'] == k) & (plot_data['region'] == region)]
        if not subset.empty:
            pos = i * 2.5 + j * 1
            positions.append(pos)
            color = arms_color if region == 'ARMS' else cen_color
            parts = ax4.violinplot([subset['pct_error_tolerant'].values],
                                   positions=[pos],
                                   widths=0.7,
                                   showmeans=True,
                                   showmedians=True)
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.7)

            if i == 0:  # Only label first occurrence
                labels.append(f'k={k}\n{region}')

ax4.set_ylabel('Error Tolerance (%)', fontweight='bold')
ax4.set_title('D. Error Tolerance Distribution', fontweight='bold', loc='left')
ax4.set_xticks([0.5, 2.5*1 + 0.5, 2.5*2 + 0.5])
ax4.set_xticklabels(['k=21', 'k=31', 'k=41'])
ax4.grid(axis='y', alpha=0.3)

# Add legend
arms_patch = mpatches.Patch(color=arms_color, label='ARMS', alpha=0.7)
cen_patch = mpatches.Patch(color=cen_color, label='CEN', alpha=0.7)
ax4.legend(handles=[arms_patch, cen_patch], loc='upper right')

# Panel E: What happens to errors?
ax5 = fig.add_subplot(gs[1, 2])

k21_data = df[df['k_size'] == 21]
outcomes = {
    'Becomes Novel\n(Lost)': k21_data['pct_becomes_novel'].mean(),
    'Stays Correct': k21_data['pct_error_tolerant'].mean(),
    'Wrong DB\n(False Pos)': k21_data['pct_wrong_db'].mean()
}

colors_pie = ['#95a5a6', '#2ecc71', '#e74c3c']
wedges, texts, autotexts = ax5.pie(outcomes.values(),
                                     labels=outcomes.keys(),
                                     colors=colors_pie,
                                     autopct='%1.1f%%',
                                     startangle=90,
                                     textprops={'fontweight': 'bold'})

for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontsize(10)

ax5.set_title('E. Error Outcomes (k=21)', fontweight='bold', loc='left')

# Panel F: Practical impact (with 1M reads)
ax6 = fig.add_subplot(gs[2, :])

read_counts = 1_000_000
usable_counts = (overall_df['usable_reads'] / 100) * read_counts
lost_counts = read_counts - usable_counts

x_pos = np.arange(len(k_sizes))
width = 0.7

bars1 = ax6.bar(x_pos, usable_counts/1000, width, label='Usable Reads',
                color='#2ecc71', edgecolor='black')
bars2 = ax6.bar(x_pos, lost_counts/1000, width, bottom=usable_counts/1000,
                label='Lost to Errors', color='#e74c3c', edgecolor='black', alpha=0.7)

ax6.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax6.set_ylabel('Read Count (thousands)', fontweight='bold', fontsize=12)
ax6.set_title('F. Practical Impact with 1 Million ONT Reads', fontweight='bold', loc='left', fontsize=13)
ax6.set_xticks(x_pos)
ax6.set_xticklabels([f'k={k}' for k in k_sizes])
ax6.legend(loc='upper right', fontsize=11)
ax6.set_ylim(0, 1100)
ax6.grid(axis='y', alpha=0.3)

# Add value labels
for i, (bar, usable, lost) in enumerate(zip(bars1, usable_counts, lost_counts)):
    # Usable count
    ax6.text(bar.get_x() + bar.get_width()/2., usable/1000 - 15,
            f'{int(usable/1000)}k',
            ha='center', va='top', fontsize=10, fontweight='bold', color='white')
    # Lost count
    ax6.text(bar.get_x() + bar.get_width()/2., (usable + lost/2)/1000,
            f'-{int(lost/1000)}k',
            ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Highlight difference
if len(usable_counts) > 1:
    diff = usable_counts.iloc[0] - usable_counts.iloc[-1]
    ax6.annotate(f'{int(diff/1000)}k more reads\nwith k=21!',
                xy=(len(k_sizes)-1, usable_counts.iloc[-1]/1000),
                xytext=(len(k_sizes)/2, 950),
                arrowprops=dict(arrowstyle='->', color='green', lw=2),
                fontsize=11, color='green', fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

plt.savefig('final_results/03_error_resilience.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/03_error_resilience.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/03_error_resilience.png")
print(f"✓ Saved: final_results/03_error_resilience.pdf")

# Print summary
print("\n" + "="*80)
print("ERROR RESILIENCE SUMMARY")
print("="*80)
print(f"{'K-mer':<8} {'Usable Reads':<15} {'Lost Reads':<15} {'Loss vs k=21':<15}")
print("-"*80)
baseline = overall_df.iloc[0]['usable_reads']
for i, row in overall_df.iterrows():
    usable = row['usable_reads']
    lost = 100 - usable
    diff = usable - baseline
    diff_str = f"{diff:+.2f}%" if i > 0 else "baseline"
    print(f"k={row['k_size']:<5} {usable:>6.2f}%         {lost:>6.2f}%         {diff_str:>12}")
print("="*80)
print(f"\n💡 With 1M reads:")
print(f"   k=21: {int(overall_df.iloc[0]['usable_reads']*10000):,} usable reads")
print(f"   k=41: {int(overall_df.iloc[-1]['usable_reads']*10000):,} usable reads")
print(f"   Difference: {int((overall_df.iloc[0]['usable_reads'] - overall_df.iloc[-1]['usable_reads'])*10000):,} more reads with k=21!")
print("="*80)
