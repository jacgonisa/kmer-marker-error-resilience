# Analysis Results - Figures and Visualizations

This document provides an overview of all generated figures from the comprehensive k-mer marker analysis.

---

## Figure 1: Marker Availability Analysis

![Marker Availability](final_results/01_marker_availability.png)

### Description
Four-panel analysis showing marker quantity and distribution across k-mer sizes.

**Panel A: Average Marker Count**
- Compares total k-mer counts between ARMS and CEN regions
- Shows clear increase in markers with larger k-mer sizes
- ARMS regions have significantly more markers than CEN (~100x more)

**Panel B: Marker Density (k-mers/Mb)**
- Density normalized per megabase
- ARMS: 100-220 million k-mers/Mb
- CEN: 90-970 thousand k-mers/Mb

**Panel C: Marker Count Distribution**
- Violin plots showing distribution across all databases
- ARMS markers show consistent high counts
- CEN markers more variable

**Panel D: Uniformity Across Databases**
- Coefficient of Variation (CV) by k-mer size
- Lower CV = more uniform distribution
- Longer k-mers (k=41) show better uniformity (CV ~30%)

### Key Finding
✅ **Larger k-mers provide more total markers but with diminishing returns**

---

## Figure 2: False Discovery Rate (FDR) Analysis

![FDR Analysis](final_results/02_cross_contamination.png)

### Description
Four-panel analysis assessing the risk of false positives from sequencing errors using **proper FDR calculation** (FDR = FP/(FP+TP)).

**Panel A: False Discovery Rate (FDR = FP/(FP+TP))**
- Shows absolute FDR for all k-mers
- All values < 0.2% - excellent performance!
- ARMS and CEN show similar low FDR
- **Critical metric**: Measures actual false positive rate among all k-mers

**Panel B: Conditional FDR (Among K-mers WITH Errors)**
- Shows FDR only for k-mers that contain errors
- ARMS: 93-94% FDR (high)
- CEN: 45-63% FDR (lower)
- **Interpretation**: When errors occur, they rarely maintain correct classification
- Most errors become "novel" k-mers (read loss, not false positive)

**Panel C: Per-Database FDR Heatmap**
- Detailed FDR values for each database
- Color scale: green (low FDR) to red (high FDR)
- Shows consistency across chromosomes

**Panel D: Database FDR Ranking**
- Databases ranked by mean FDR
- CEN markers generally have lower FDR (better error tolerance)
- ARMS markers more susceptible to errors but still excellent

### Key Finding
✅ **High conditional FDR (45-100%) but ~99% of errors become novel k-mers**
✅ **Absolute FDR < 0.2% for all k-mer sizes - excellent specificity!**

### Important Note on FDR Interpretation
- **Conditional FDR**: Among k-mers WITH errors, what % are false positives?
  - High values (45-100%) indicate most errors can't maintain classification
- **Absolute FDR**: Among ALL k-mers, what % become false positives?
  - Very low (<0.2%) because most errors become "novel" k-mers
- **Net result**: Errors cause read loss, not false discoveries

---

## Figure 3: Error Resilience and Read Retention

![Error Resilience](final_results/03_error_resilience.png)

### Description
Six-panel comprehensive analysis of read retention under realistic 1% per-base error rate (ONT R9-like).

**Panel A: Read Retention by K-mer Size** ⭐ MOST IMPORTANT
- k=21: **81.1% retention** (best!)
- k=25: 77.9% retention
- k=31: 73.4% retention
- k=35: 70.5% retention
- k=41: 66.3% retention (worst)
- **15% difference between k=21 and k=41!**

**Panel B: K-mers Affected by Errors**
- Shows percentage of k-mers containing ≥1 error
- k=21: 19% affected
- k=41: 34% affected
- Follows formula: P(error) = 1 - (0.99)^k

**Panel C: Error Tolerance by Region**
- ARMS vs CEN comparison
- Both regions show similar trends
- Shorter k-mers consistently better

**Panel D: Distribution of Error Tolerance**
- Box plots showing variation across databases
- Tight distributions = consistent behavior
- No major outliers

**Panel E: Error Outcomes**
- What happens to k-mers with errors?
- ~99% become "novel" (not in any database) → read loss
- ~0.06% remain correctly classified (error-tolerant)
- ~0.7% match wrong database (false positive)
- ~0% ambiguous (match multiple databases)

**Panel F: Practical Impact Visualization**
- With 1 million ONT reads:
- k=21: 811,000 usable reads
- k=41: 663,000 usable reads
- **148,000 more usable reads with k=21!**

### Key Finding
✅ **k=21 retains 15% more reads than k=41 under realistic ONT error rates**
✅ **Most errors cause read loss (novel k-mers), not false positives**

---

## Figure 4: Final Recommendation and Integrated Scoring

![Final Recommendation](final_results/04_final_recommendation.png)

### Description
Five-panel integrated analysis combining all metrics with weighted scoring system.

**Panel A: Overall Weighted Scores** ⭐ FINAL RECOMMENDATION
- k=21: **70/100** (Winner!)
- k=25: 61/100
- k=31: 49/100
- k=35: 42/100
- k=41: 30/100

**Panel B: Scoring Criteria Explanation**
- Read Retention: 40% weight (most critical)
- Specificity: 30% weight (FDR)
- Uniformity: 15% weight (CV)
- Availability: 15% weight (total markers)

**Panel C: Radar Chart (k=21 vs k=41)**
- Pentagon showing all 5 dimensions
- k=21 dominates in retention and resilience
- k=41 better in specificity and availability
- Clear visual trade-off

**Panel D: Individual Component Scores**
- Breakdown of each metric by k-mer size
- Shows how scores were normalized (0-100 scale)

**Panel E: Complete Decision Table**
- Raw values for all metrics
- Allows custom decision-making based on priorities

### Key Finding
✅ **k=21 recommended for ONT sequencing due to superior read retention**
✅ **Weighted scoring accounts for importance of each factor**

---

## Figure 5: Pentagon Radar Comparison

![Radar Comparison](final_results/05_radar_comparison.png)

### Description
Comprehensive pentagon radar charts comparing all five k-mer sizes across five key metrics.

**Metrics Displayed:**
1. **Usable Reads (%)** - Read retention after errors
2. **Specificity (%)** - 100 - FDR
3. **Marker Availability (millions)** - Total k-mer count
4. **Uniformity (%)** - 100 - CV (distribution consistency)
5. **Error Resilience (%)** - Same as usable reads

**Visual Interpretation:**
- Larger pentagon area = better overall performance
- k=21 shows large area in retention/resilience dimensions
- k=41 shows large area in specificity/markers dimensions
- k=31 provides middle ground

### Key Finding
✅ **Visual confirmation: k=21 dominates critical metrics (retention, resilience)**
✅ **k=31 offers balanced compromise if higher marker count needed**

---

## Figure 6: Coverage Loss Analysis

![Coverage Loss](final_results/06_coverage_loss_analysis.png)

### Description
Detailed analysis of how errors affect k-mer coverage per megabase.

Shows:
- Initial marker density (k-mers/Mb)
- Effective coverage after accounting for errors
- Coverage loss by k-mer size
- Regional differences (ARMS vs CEN)

### Key Finding
✅ **Shorter k-mers maintain higher effective coverage despite fewer total markers**

---

## Figure 7: Objective Comparison (Trade-off Analysis)

![Objective Comparison](final_results/07_objective_comparison.png)

### Description
Objective comparison showing all metrics side-by-side without arbitrary weighting.

Allows users to make decisions based on their own priorities:
- **Priority: Error resilience** → Choose k=21
- **Priority: Maximum markers** → Choose k=41
- **Priority: Balance** → Choose k=31

### Key Finding
✅ **No universal "best" k-mer - depends on use case and data quality**

---

## Summary Table: All Key Metrics

| K-mer | Read Retention | FDR (Conditional) | Error Impact | Total Markers | Recommendation |
|-------|---------------|-------------------|--------------|---------------|----------------|
| **k=21** | **81.10%** | 93.1% (ARMS), 62.9% (CEN) | 19.0% | 3.5M | ⭐ **Best for ONT R9** |
| k=25 | 77.93% | 93.6% (ARMS), 59.2% (CEN) | 22.2% | 4.0M | Conservative |
| **k=31** | 73.36% | 92.9% (ARMS), 55.7% (CEN) | 26.7% | 4.7M | ⭐ **Balanced** |
| k=35 | 70.47% | 93.9% (ARMS), 49.4% (CEN) | 29.6% | 5.1M | High coverage |
| k=41 | 66.34% | 93.8% (ARMS), 45.1% (CEN) | 33.7% | 5.6M | Maximum resolution |

---

## Scientific Interpretation

### Understanding FDR Results

The FDR analysis reveals an important pattern:

1. **High Conditional FDR (45-100%)**:
   - Among k-mers that contain errors AND still match a database, most are false positives
   - This means errors rarely preserve the correct database classification

2. **Low Absolute FDR (<0.2%)**:
   - Among ALL k-mers, very few become false positives
   - Why? Because ~99% of errors create "novel" k-mers that don't match any database

3. **Net Effect**:
   - Errors primarily cause **read loss** (novel k-mers)
   - False positives are rare (<0.2% of all k-mers)
   - This is actually good for specificity but bad for sensitivity

### Error Probability Mathematics

**Formula**: P(k-mer contains ≥1 error) = 1 - (0.99)^k

At 1% per-base error rate (ONT R9):
- k=21: P(error) = 1 - (0.99)^21 = **19.0%** → 81% retention
- k=31: P(error) = 1 - (0.99)^31 = **26.7%** → 73% retention
- k=41: P(error) = 1 - (0.99)^41 = **33.7%** → 66% retention

Each additional base adds ~1% chance of containing an error.

---

## Use Case Recommendations

### 1. Standard Arabidopsis Crossover Analysis
**Recommendation**: **k=31** (balanced)
- Good retention (73%)
- High marker count (4.7M)
- Excellent FDR control
- Best all-around choice

### 2. Low-Quality or Early ONT Data (R9)
**Recommendation**: **k=21** (error-resilient)
- Maximum retention (81%)
- Minimizes read loss
- Still excellent specificity
- Trade-off: Fewer total markers

### 3. High-Quality ONT Data (R10.4+)
**Recommendation**: **k=41** (maximum resolution)
- Maximum markers (5.6M)
- Best uniformity
- Best specificity
- Error rate lower in R10.4, so retention penalty less severe

### 4. Publication-Quality Results
**Recommendation**: Run **both k=31 and k=41**, compare results
- Demonstrates robustness
- Provides sensitivity analysis
- Shows results are not k-mer dependent

---

## Files Generated

All figures available in both formats:
- **PNG** (300 dpi, 500KB - 1.5MB) - For presentations and documents
- **PDF** (Vector, 30-80KB) - For publications and scalable graphics

```
final_results/
├── 01_marker_availability.png/pdf         (503KB / 33KB)
├── 02_cross_contamination.png/pdf         (898KB / 47KB) ⭐ Updated with FDR
├── 03_error_resilience.png/pdf            (594KB / 63KB)
├── 04_final_recommendation.png/pdf        (628KB / 60KB)
├── 05_radar_comparison.png/pdf            (1.4MB / 51KB)
├── 06_coverage_loss_analysis.png/pdf      (930KB / 65KB)
└── 07_objective_comparison.png/pdf        (1.1MB / 76KB)
```

---

## Citation

If you use these figures in publications, please cite:

```
Gonçalves, J. (2025). K-mer Marker Error Resilience Analysis for ONT Sequencing.
GitHub: https://github.com/jacgonisa/kmer-marker-error-resilience
```

And the KMC tool:
```
Deorowicz S, Kokot M, Grabowski S, Debudaj-Grabysz A.
KMC 2: Fast and resource-frugal k-mer counting.
Bioinformatics. 2015;31(10):1569-1576.
```

---

**Last Updated**: December 4, 2025
**Analysis Version**: 2.0 (FDR calculations implemented)
