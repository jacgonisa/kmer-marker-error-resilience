#!/usr/bin/env python3
"""
Plot 4: Final Recommendation - Integrated Decision Matrix
Combines all factors to recommend optimal k-mer size.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import subprocess

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# Load error resilience data
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

# Calculate absolute false positive rate
error_df['absolute_false_positive_rate'] = (error_df['pct_kmers_with_errors'] / 100) * (error_df['pct_wrong_db'] / 100) * 100

# Load marker availability data
marker_avail_file = Path("final_results/marker_availability_summary.csv")
if marker_avail_file.exists():
    marker_df = pd.read_csv(marker_avail_file)
    # Merge marker counts with error data
    error_df = error_df.merge(marker_df[['k_size', 'database', 'total_kmers']],
                               on=['k_size', 'database'], how='left')

# Calculate comprehensive scores for each k-mer size
scores = []

for k in k_sizes:
    k_data = error_df[error_df['k_size'] == k]

    # 1. Read retention (higher is better)
    pct_with_errors = k_data['pct_kmers_with_errors'].mean()
    pct_error_tolerant = k_data['pct_error_tolerant'].mean()
    usable_reads = (100 - pct_with_errors) + (pct_with_errors * pct_error_tolerant / 100)

    # 2. Specificity (lower false positive is better)
    false_positive_rate = k_data['absolute_false_positive_rate'].mean()

    # 3. Uniformity (lower CV is better)
    arms_counts = k_data[k_data['region'] == 'ARMS']['total_kmers']
    cen_counts = k_data[k_data['region'] == 'CEN']['total_kmers']

    arms_cv = (arms_counts.std() / arms_counts.mean()) * 100 if len(arms_counts) > 1 else 0
    cen_cv = (cen_counts.std() / cen_counts.mean()) * 100 if len(cen_counts) > 1 else 0
    avg_cv = (arms_cv + cen_cv) / 2

    # 4. Average marker count
    avg_marker_count = k_data['total_kmers'].mean()

    scores.append({
        'k_size': k,
        'usable_reads': usable_reads,
        'false_positive_rate': false_positive_rate,
        'uniformity_cv': avg_cv,
        'avg_marker_count': avg_marker_count,
        'specificity_score': 100 - (false_positive_rate * 100),  # Higher is better
    })

scores_df = pd.DataFrame(scores)

# Normalize scores to 0-100 scale
def normalize_score(values, higher_is_better=True):
    """Normalize to 0-100 scale."""
    vmin, vmax = values.min(), values.max()
    if vmax == vmin:
        return pd.Series([50] * len(values))
    if higher_is_better:
        return ((values - vmin) / (vmax - vmin)) * 100
    else:
        return ((vmax - values) / (vmax - vmin)) * 100

scores_df['read_retention_score'] = normalize_score(scores_df['usable_reads'], higher_is_better=True)
scores_df['specificity_norm'] = normalize_score(scores_df['false_positive_rate'], higher_is_better=False)
scores_df['uniformity_score'] = normalize_score(scores_df['uniformity_cv'], higher_is_better=False)
scores_df['availability_score'] = normalize_score(scores_df['avg_marker_count'], higher_is_better=True)

# Calculate weighted overall score
# Since ALL k-mer sizes have excellent specificity (<0.2%), and ALL have plenty of markers,
# we should heavily prioritize read retention (the ONLY factor with big practical differences)
weights = {
    'read_retention': 0.70,  # CRITICAL - 15% difference between k=21 and k=41!
    'specificity': 0.20,     # All excellent (<0.2%), small practical impact
    'uniformity': 0.05,      # Minor factor
    'availability': 0.05     # Minor factor (all have sufficient markers)
}

scores_df['overall_score'] = (
    scores_df['read_retention_score'] * weights['read_retention'] +
    scores_df['specificity_norm'] * weights['specificity'] +
    scores_df['uniformity_score'] * weights['uniformity'] +
    scores_df['availability_score'] * weights['availability']
)

print("="*80)
print("COMPREHENSIVE K-MER EVALUATION SCORES")
print("="*80)
print(scores_df.to_string(index=False))
print("="*80)

# Create visualization
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

fig.suptitle('Comprehensive K-mer Size Evaluation & Recommendation',
             fontsize=17, fontweight='bold', y=0.98)

# Panel A: Overall Score (BIG - most important!)
ax1 = fig.add_subplot(gs[0, :2])

colors = ['#2ecc71' if score == scores_df['overall_score'].max() else '#3498db'
          for score in scores_df['overall_score']]
bars = ax1.bar(range(len(k_sizes)), scores_df['overall_score'],
               color=colors, edgecolor='black', linewidth=2)

ax1.set_xlabel('K-mer Size', fontweight='bold', fontsize=12)
ax1.set_ylabel('Overall Score (0-100)', fontweight='bold', fontsize=12)
ax1.set_title('A. Overall Performance Score ★ FINAL RECOMMENDATION',
              fontweight='bold', loc='left', fontsize=14, color='darkred')
ax1.set_xticks(range(len(k_sizes)))
ax1.set_xticklabels([f'k={k}' for k in k_sizes])
ax1.set_ylim(0, 110)
ax1.grid(axis='y', alpha=0.3)

# Add value labels
for bar, score in zip(bars, scores_df['overall_score']):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
            f'{score:.1f}',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

# Highlight winner
best_idx = scores_df['overall_score'].idxmax()
best_k = scores_df.loc[best_idx, 'k_size']
bars[best_idx].set_edgecolor('gold')
bars[best_idx].set_linewidth(4)

ax1.annotate(f'★ RECOMMENDED: k={int(best_k)}',
            xy=(best_idx, scores_df.loc[best_idx, 'overall_score']),
            xytext=(best_idx + 0.5, scores_df.loc[best_idx, 'overall_score'] - 10),
            arrowprops=dict(arrowstyle='->', color='gold', lw=3),
            fontsize=13, color='darkgreen', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='gold', alpha=0.4, edgecolor='darkgreen', linewidth=2))

# Panel B: Scoring weights explanation
ax2 = fig.add_subplot(gs[0, 2])
ax2.axis('off')

weight_text = """
SCORING CRITERIA
━━━━━━━━━━━━━━━━━━━━━
Optimized for ONT sequencing

📊 Read Retention: 70%
   ★★★ MOST CRITICAL
   81% vs 66% usable reads
   = 15% practical difference!

🎯 Specificity: 20%
   All <0.2% (excellent!)
   Minimal practical impact

📏 Uniformity: 5%
   All reasonably uniform

📦 Availability: 5%
   All have enough markers

━━━━━━━━━━━━━━━━━━━━━
★ Overall Score: Optimized
   for maximum read retention
"""

ax2.text(0.1, 0.5, weight_text, fontsize=10, family='monospace',
        verticalalignment='center',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=2))

# Panel C: Radar chart comparing k=21 vs k=41
ax3 = fig.add_subplot(gs[1, 0], projection='polar')

categories = ['Read\nRetention', 'Specificity', 'Uniformity', 'Availability']
k21_scores = scores_df[scores_df['k_size'] == 21][
    ['read_retention_score', 'specificity_norm', 'uniformity_score', 'availability_score']
].values[0]
k41_scores = scores_df[scores_df['k_size'] == 41][
    ['read_retention_score', 'specificity_norm', 'uniformity_score', 'availability_score']
].values[0]

angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
k21_scores_plot = np.concatenate((k21_scores, [k21_scores[0]]))
k41_scores_plot = np.concatenate((k41_scores, [k41_scores[0]]))
angles += angles[:1]

ax3.plot(angles, k21_scores_plot, 'o-', linewidth=2, label='k=21', color='#2ecc71')
ax3.fill(angles, k21_scores_plot, alpha=0.25, color='#2ecc71')
ax3.plot(angles, k41_scores_plot, 'o-', linewidth=2, label='k=41', color='#e74c3c')
ax3.fill(angles, k41_scores_plot, alpha=0.25, color='#e74c3c')

ax3.set_xticks(angles[:-1])
ax3.set_xticklabels(categories, fontsize=9)
ax3.set_ylim(0, 100)
ax3.set_title('C. k=21 vs k=41 Comparison', fontweight='bold', pad=20)
ax3.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
ax3.grid(True)

# Panel D: Individual component scores
ax4 = fig.add_subplot(gs[1, 1:])

x = np.arange(len(k_sizes))
width = 0.2

components = [
    ('read_retention_score', 'Read Retention', '#2ecc71'),
    ('specificity_norm', 'Specificity', '#3498db'),
    ('uniformity_score', 'Uniformity', '#f39c12'),
    ('availability_score', 'Availability', '#9b59b6')
]

for i, (col, label, color) in enumerate(components):
    offset = (i - 1.5) * width
    bars = ax4.bar(x + offset, scores_df[col], width, label=label, color=color, alpha=0.8)

ax4.set_xlabel('K-mer Size', fontweight='bold')
ax4.set_ylabel('Component Score (0-100)', fontweight='bold')
ax4.set_title('D. Individual Component Scores', fontweight='bold', loc='left')
ax4.set_xticks(x)
ax4.set_xticklabels([f'k={k}' for k in k_sizes])
ax4.legend(loc='upper right', ncol=2)
ax4.grid(axis='y', alpha=0.3)
ax4.set_ylim(0, 105)

# Panel E: Decision matrix table
ax5 = fig.add_subplot(gs[2, :])
ax5.axis('tight')
ax5.axis('off')

# Create table data
table_data = []
table_data.append(['K-mer\nSize', 'Usable\nReads (%)', 'False Pos.\nRate (%)',
                   'Marker\nCount', 'Uniformity\nCV (%)', 'Overall\nScore', 'Rank'])

for idx, row in scores_df.iterrows():
    rank = len(scores_df) - (scores_df['overall_score'] > row['overall_score']).sum()
    rank_str = f"#{rank}"
    if rank == 1:
        rank_str = "★ #1"

    table_data.append([
        f"k={int(row['k_size'])}",
        f"{row['usable_reads']:.2f}",
        f"{row['false_positive_rate']:.4f}",
        f"{int(row['avg_marker_count']):,}",
        f"{row['uniformity_cv']:.1f}",
        f"{row['overall_score']:.1f}",
        rank_str
    ])

table = ax5.table(cellText=table_data, cellLoc='center', loc='center',
                 bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.5)

# Style header row
for i in range(7):
    cell = table[(0, i)]
    cell.set_facecolor('#34495e')
    cell.set_text_props(weight='bold', color='white')

# Style data rows
# Calculate ranks for all scores
scores_df['rank'] = scores_df['overall_score'].rank(ascending=False, method='min')

for i in range(1, len(table_data)):
    rank = int(scores_df.iloc[i-1]['rank'])
    if rank == 1:
        row_color = '#d4edda'  # Light green for best
        for j in range(7):
            cell = table[(i, j)]
            cell.set_facecolor(row_color)
            cell.set_text_props(weight='bold')
    else:
        row_color = '#f8f9fa' if i % 2 == 0 else 'white'
        for j in range(7):
            table[(i, j)].set_facecolor(row_color)

ax5.set_title('E. Complete Decision Matrix', fontweight='bold', loc='left', fontsize=13, pad=10)

plt.savefig('final_results/04_final_recommendation.png', dpi=300, bbox_inches='tight')
plt.savefig('final_results/04_final_recommendation.pdf', bbox_inches='tight')
print(f"\n✓ Saved: final_results/04_final_recommendation.png")
print(f"✓ Saved: final_results/04_final_recommendation.pdf")

# Print final recommendation
print("\n" + "="*80)
print("🎯 FINAL RECOMMENDATION")
print("="*80)
best_row = scores_df.loc[best_idx]
print(f"\n★ RECOMMENDED K-MER SIZE: k={int(best_k)}")
print(f"\nOverall Score: {best_row['overall_score']:.1f}/100")
print(f"\nKey Metrics:")
print(f"  • Read Retention:       {best_row['usable_reads']:.2f}%")
print(f"  • False Positive Rate:  {best_row['false_positive_rate']:.4f}%")
print(f"  • Average Marker Count: {int(best_row['avg_marker_count']):,}")
print(f"  • Uniformity (CV):      {best_row['uniformity_cv']:.1f}%")

print(f"\n💡 Why k={int(best_k)} is best:")
if best_k == 21:
    print(f"  ✓ Highest read retention under ONT sequencing errors")
    print(f"  ✓ Excellent specificity (very low false positive rate)")
    print(f"  ✓ Best balance of all factors")
    print(f"  ✓ Shorter k-mers = fewer bases = less likely to contain errors")

print("\n" + "="*80)

# Save scores to CSV
scores_df.to_csv('final_results/comprehensive_scores.csv', index=False)
print(f"✓ Saved scores: final_results/comprehensive_scores.csv")
