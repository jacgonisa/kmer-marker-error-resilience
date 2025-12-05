# K-mer Marker Error Resilience Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Version](https://img.shields.io/badge/version-2.1-blue.svg)](https://github.com/jacgonisa/kmer-marker-error-resilience)

A comprehensive analysis framework for evaluating k-mer marker performance across different k-mer sizes (21, 25, 31, 35, 41) under realistic Oxford Nanopore sequencing error conditions.

---

## ⚠️ Important Version Information

**Current Version: 2.1 (December 2025)**

**Critical Update**: Version 2.1 fixes a **major bug** in absolute FDR calculation discovered in December 2025. Previous versions incorrectly showed FDR values ~100x too high due to using conditional FDR in the calculation. All results shown here use the corrected calculations.

**What Changed in v2.1**:
- ✅ **FIXED**: Absolute FDR calculation (was 17.7%, now correctly 0.138% for k=21)
- ✅ **IMPROVED**: Visualization showing novel vs cross-contamination distinction
- ✅ **ADDED**: Panel D showing direct k=21 vs k=41 comparison
- ✅ **CLARIFIED**: Biological interpretation - errors cause information loss, not false positives

**For detailed FDR explanation**: See [FDR_EXPLAINED.md](FDR_EXPLAINED.md) - explains why high conditional FDR (93%) actually indicates excellent specificity (<0.2%)!

---

## Overview

This repository contains production-quality genomic analysis tools to systematically evaluate k-mer markers for **Arabidopsis thaliana** crossover detection using cenhapmer markers. The analysis simulates realistic 1% per-base error rates (ONT R9-like) to assess:

- **Marker Availability**: Quantity and distribution of k-mers
- **Error Resilience**: K-mer retention under sequencing errors
- **False Discovery Rate (FDR)**: Risk of misclassification from errors
- **Cross-Contamination**: Biological false positive risk
- **Optimal k-mer Selection**: Data-driven recommendations

---

## Key Findings

### 🔬 Main Discovery: Errors Cause Information Loss, Not False Positives

**When sequencing errors occur in k-mers**:
- **~99%** become "novel" k-mers (not found in any database) → **Information loss** (reduced coverage)
- **<1%** match wrong database → **False positives (cross-contamination)** - the dangerous outcome!
- **<0.1%** remain correctly classified → **Error-tolerant**

**Result**: Excellent specificity (FDR <0.2%, cross-contamination <1%) but reduced k-mer retention (66-81%)

### 📊 Recommended k-mer Size: **k=31** (Balanced) or **k=21** (Error-Resilient)

| K-mer | K-mer Retention | Absolute FDR | Cross-Contam | Error Impact | Total Markers | Use Case |
|-------|-----------------|--------------|--------------|--------------|---------------|----------|
| **k=21** | **81.10%** | **0.138%** | 0.725% | 19.0% | 3.5M | High-error data (ONT R9) ⭐ |
| k=25 | 77.93% | 0.120% | 0.541% | 22.2% | 4.0M | Conservative approach |
| **k=31** | 73.36% | **0.110%** | 0.412% | 26.7% | 4.7M | **Balanced (RECOMMENDED)** ⭐ |
| k=35 | 70.47% | 0.102% | 0.345% | 29.6% | 5.1M | High coverage priority |
| k=41 | 66.34% | **0.092%** | **0.274%** | 33.7% | 5.6M | Maximum resolution (R10.4+) |

### The Trade-off:

**Longer k-mers (k=41)**:
- ✅ BEST specificity (lowest cross-contamination: 0.274%)
- ✅ MOST total markers per Mb
- ✅ BEST uniformity
- ❌ WORST error resilience (34% of k-mers affected by errors)
- ❌ LOWEST retention (66% usable k-mers)

**Shorter k-mers (k=21)**:
- ✅ BEST error resilience (only 19% affected)
- ✅ HIGHEST retention (81% usable k-mers)
- ❌ Slightly higher cross-contamination (0.725% vs 0.274%)
- ❌ Fewer total markers available

**🎯 Recommendation**: The small difference in cross-contamination (0.45%) does NOT justify the 15% loss in usable k-mers. Choose **k=31** for balanced performance or **k=21** for error-resilient analysis.

### Scientific Rationale

**Error probability**: P(k-mer contains ≥1 error) = 1 - (0.99)^k

- k=21: 19% error probability → **81% k-mers retained**
- k=41: 34% error probability → **66% k-mers retained** (15% loss!)

**Error calculation**: Each base has 1% error probability. For a k-mer to be error-free, ALL bases must be correct: (0.99)^k. The complement gives probability of ≥1 error.

---

## Installation

### Prerequisites

- Python 3.8 or higher
- KMC (K-mer Counter) tool

### Step 1: Clone Repository

```bash
git clone https://github.com/jacgonisa/kmer-marker-error-resilience.git
cd kmer-marker-error-resilience
```

### Step 2: Install Python Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- pandas >= 1.5.0
- numpy >= 1.23.0
- matplotlib >= 3.5.0
- seaborn >= 0.12.0

### Step 3: Install KMC (K-mer Counter)

**Option A: Using Conda (Recommended)**
```bash
conda install -c bioconda kmc
```

**Option B: Manual Installation**
```bash
wget https://github.com/refresh-bio/KMC/releases/latest/download/KMC3.linux.tar.gz
tar -xzvf KMC3.linux.tar.gz
# Add to PATH or move to /usr/local/bin/
```

**Verify Installation**
```bash
kmc_tools --help
```

---

## Usage

### Quick Start: Run All Analyses

```bash
chmod +x run_all_analyses.sh
./run_all_analyses.sh
```

This will:
1. Create `final_results/` directory
2. Run analyses 1-5 sequentially
3. Generate publication-quality figures (PNG 300dpi + PDF vector)
4. Save summary statistics to CSV files
5. Print comprehensive reports

**Expected runtime**: 2-5 minutes

**Note**: Scripts 6-7 (coverage loss, objective comparison) are supplementary and run separately if needed.

### Run Individual Analyses

```bash
# 1. Marker availability and density
python3 01_marker_availability.py

# 2. False Discovery Rate (FDR) and cross-contamination analysis
python3 02_cross_contamination.py

# 3. Error resilience and k-mer retention
python3 03_error_resilience.py

# 4. Integrated recommendation system
python3 04_final_recommendation.py

# 5. Radar chart comparison
python3 05_radar_comparison.py

# Optional supplementary analyses:
python3 06_coverage_loss_analysis.py
python3 07_objective_comparison.py
```

---

## Analysis Modules

### 1. Marker Availability Analysis (`01_marker_availability.py`)

![Marker Availability](final_results/01_marker_availability.png)

**Evaluates**: Marker quantity and distribution

**What it shows**:
- **Panel A**: Total k-mer counts (ARMS vs CEN regions)
- **Panel B**: Marker density normalized per megabase
- **Panel C**: Distribution of marker counts across databases
- **Panel D**: Uniformity (Coefficient of Variation)

**Key Finding**: Larger k-mers (k=41) provide 5.6M markers vs 3.5M for k=21, but with diminishing returns in density.

---

### 2. Cross-Contamination Risk Analysis (`02_cross_contamination.py`)

![Cross-Contamination Analysis](final_results/02_cross_contamination.png)

**Evaluates**: False Discovery Rate (FDR) and cross-contamination risk from sequencing errors

**What is FDR?**
- **FDR = FP / (FP + TP)** - False Discovery Rate
- **FP** (False Positives): K-mers with errors matching WRONG database (cross-contamination)
- **TP** (True Positives): K-mers with errors still matching CORRECT database
- Standard metric in bioinformatics for assessing discovery quality

**What it shows**:

**Panel A: Absolute FDR (among ALL k-mers)** ⭐ **MOST IMPORTANT!**
- Shows <0.2% for all k-mer sizes ✓ Excellent!
- This is the % of ALL k-mers that become false positives
- k=21: 0.138%, k=31: 0.110%, k=41: 0.092%
- All values WAY better than typical 5% FDR threshold

**Panel B: Cross-Contamination Rate (FALSE POSITIVE Risk)**
- Shows % of error-containing k-mers matching WRONG database
- ARMS: ~0.7% (k=21) to ~0.3% (k=41) - all excellent!
- CEN: ~0.3-0.6% - even better!
- **All values <1%** ✓ Excellent specificity!
- **Key insight**: ~99% of errors become "novel" (information loss), <1% become false positives

**Panel C: Per-Database Cross-Contamination Heatmap**
- Shows false positive risk for each database
- Color scale: green (low, <0.3%) to red (higher, >0.8%)
- All values <1% across all databases and chromosomes ✓
- CEN markers generally show lower cross-contamination than ARMS

**Panel D: Error Fate Comparison (k=21 vs k=41)**
- Stacked bars showing novel (information loss) vs cross-contamination (FP)
- Direct visual comparison of extreme k-mer sizes
- Blue: Novel k-mers (lost to errors, no false positive)
- Red hatched: Cross-contamination (FALSE POSITIVES - rare but dangerous!)
- Shows k=41 has LOWER cross-contamination (0.27% vs 0.73%) but HIGHER information loss

**Understanding the Results** (See [FDR_EXPLAINED.md](FDR_EXPLAINED.md) for details):

When a sequencing error occurs in a k-mer, 3 outcomes are possible:

```
100,000 k-mers with 19% error rate (k=21):
    ↓
19,000 get errors
    ↓
    ├─ 18,868 (99.2%) → Novel (don't match ANY database) → READ LOST ✗
    ├─ 57 (0.06%) → Still match CORRECT database → TP ✓
    └─ 75 (0.7%) → Match WRONG database → FP (cross-contamination) ✗✗

Absolute FDR = 75/100,000 = 0.075% ✓ Excellent!
```

**Key Finding**:
- ✅ Excellent specificity: Absolute FDR <0.2%, cross-contamination <1%
- ⚠️ Main problem: Information loss (19-34% of k-mers), not false positives
- 🎯 Longer k-mers have BETTER specificity but WORSE retention

---

### 3. Error Resilience Analysis (`03_error_resilience.py`)

![Error Resilience](final_results/03_error_resilience.png)

**Evaluates**: K-mer retention under realistic ONT sequencing (1% per-base error rate)

**What it shows**:

**Panel A: K-mer Retention by K-mer Size** ⭐ **KEY METRIC**
- k=21: **81% retention** (BEST!)
- k=41: 66% retention (15% worse!)
- Direct impact: With 1M k-mers, k=21 retains **148,000 more** than k=41

**Panel B: K-mers Affected by Errors**
- Shows percentage containing ≥1 error
- Follows formula: P(error) = 1 - (0.99)^k
- k=21: 19%, k=41: 34%

**Panel C: Error Tolerance by Region**
- ARMS vs CEN comparison
- Both regions show similar trends
- Shorter k-mers consistently better

**Panel D: Distribution of Error Tolerance**
- Box plots showing variation across databases
- Tight distributions = consistent behavior

**Panel E: Error Outcomes Breakdown**
- Shows the 99% novel / 0.06% TP / 0.7% FP split
- Visualizes why most errors cause information loss, not false positives

**Panel F: Practical Impact with 1M Reads**
- k=21: 811,000 usable k-mers
- k=41: 663,000 usable k-mers
- **148,000 more usable k-mers with k=21!**

**Key Finding**: k=21 retains **15% more k-mers** than k=41 under ONT error rates. This is a substantial practical advantage that outweighs the small difference in cross-contamination (0.45%).

---

### 4. Final Recommendation (`04_final_recommendation.py`)

![Final Recommendation](final_results/04_final_recommendation.png)

**Integrates**: All factors with weighted scoring system

**Scoring System**:
- K-mer Retention: 40% (critical for ONT)
- Specificity: 30% (FDR control)
- Uniformity: 15% (distribution)
- Availability: 15% (marker count)

**What it shows**:
- **Panel A**: Overall weighted scores → **k=21 wins (70/100)**
- **Panel B**: Explanation of scoring criteria
- **Panel C**: Radar chart (k=21 vs k=41 comparison)
- **Panel D**: Individual component scores
- **Panel E**: Complete decision table with raw values

**Recommendation**: **k=21** for ONT R9 data, **k=31** for balanced approach, **k=41** for high-quality R10.4+ data

**Note**: This weighted scoring is one approach. See script 07 for objective comparison without arbitrary weights.

---

### 5. Pentagon Radar Comparison (`05_radar_comparison.py`)

![Radar Comparison](final_results/05_radar_comparison.png)

**What it shows**: All k-mer sizes compared across 5 dimensions simultaneously
- Usable k-mers (%)
- Specificity (%)
- Marker availability (millions)
- Uniformity (%)
- Error resilience (%)

**Visual interpretation**: Larger pentagon area = better overall performance. k=21 shows large area in retention/resilience dimensions, k=41 in specificity/markers.

---

### 6. Coverage Loss Analysis (`06_coverage_loss_analysis.py`) - Supplementary

![Coverage Loss](final_results/06_coverage_loss_analysis.png)

**What it shows**: How errors affect effective k-mer coverage per megabase
- Initial density vs effective density after errors
- Coverage loss by k-mer size and region

---

### 7. Objective Comparison (`07_objective_comparison.py`) - Supplementary

![Objective Comparison](final_results/07_objective_comparison.png)

**What it shows**: Side-by-side metric comparison without arbitrary weighting
- Allows custom decision-making based on your priorities
- No single "best" k-mer - depends on use case

---

## Output Files

### Figures (PNG 300dpi + PDF vectors)

```
final_results/
├── 01_marker_availability.png/pdf     # Marker quantity and density
├── 02_cross_contamination.png/pdf     # FDR and cross-contamination analysis
├── 03_error_resilience.png/pdf        # K-mer retention analysis
├── 04_final_recommendation.png/pdf    # Integrated decision matrix
├── 05_radar_comparison.png/pdf        # Pentagon radar charts
├── 06_coverage_loss_analysis.png/pdf  # Coverage impact (supplementary)
└── 07_objective_comparison.png/pdf    # Trade-off visualization (supplementary)
```

### Data Files (CSV)

```
final_results/
├── marker_availability_summary.csv              # Marker counts and density
├── comprehensive_scores.csv                     # Final scoring matrix
├── realistic_k{21,25,31,35,41}_100k_error_resilience_stats.csv  # Error simulation results
└── realistic_k{21,25,31,35,41}_100k_events.csv.gz              # Detailed event data (gzipped)
```

---

## Input Data Requirements

The analysis expects KMC databases in the following structure:

```
../../03-cenhapmers/
├── k21/
│   ├── Col-0_ARMS_Chr1.kmc_pre
│   ├── Col-0_ARMS_Chr1.kmc_suf
│   └── ...
├── k25/ ...
├── k31/ ...
├── k35/ ...
└── k41/ ...
```

**Database naming convention**: `{genotype}_{region}_{chromosome}.kmc_*`
- Genotypes: Col-0, Ler-0
- Regions: ARMS, CEN
- Chromosomes: Chr1-Chr5

---

## Methodology

### Error Simulation

**Model**: Per-base 1% error rate (ONT R9-like)
- Bernoulli distribution for each base
- Random substitutions (A↔T↔G↔C)
- 100,000 k-mers tested per database

**Process**:
1. Sample k-mers from database
2. For each base: 1% chance of error (random substitution)
3. Test mutated k-mer against ALL databases
4. Classify outcome: Novel, True Positive, False Positive, or Ambiguous

**Classification**:
- **Novel**: Not found in any database (information loss)
- **True Positive (TP)**: Still matches correct database (error-tolerant)
- **False Positive (FP)**: Matches different database (cross-contamination)
- **Ambiguous**: Matches multiple databases

### False Discovery Rate (FDR)

**Conditional FDR** = FP / (FP + TP)
- Among k-mers with errors that still match a database, what % are wrong?
- High values (45-100%) mean errors rarely preserve correct classification
- This is expected - errors create novel sequences

**Absolute FDR** = (% with errors) × (% that become FP)
- Among ALL k-mers, what % become false positives?
- Low values (<0.2%) mean excellent specificity
- This is what matters for final results

**Why the difference?**
- Most errors (99%) create "novel" k-mers that don't match ANY database
- These disappear from analysis (information loss) but don't cause false positives
- Only ~1% of errors match a database, and of those, most (~93%) match the wrong one
- But 1% of 19% = 0.19% of all k-mers become false positives

**See [FDR_EXPLAINED.md](FDR_EXPLAINED.md) for detailed explanation with examples!**

### Scoring System (Script 4)

All metrics normalized to 0-100 scale:
- **Higher is better**: K-mer retention, marker count
- **Lower is better**: FDR, CV (uniformity)

**Weighted combination**:
```
Overall Score = 0.40 × Retention + 0.30 × Specificity +
                0.15 × Uniformity + 0.15 × Availability
```

**Note**: These weights reflect priorities for ONT data (retention > specificity). Adjust if your priorities differ!

---

## Citation

If you use this analysis framework in your research, please cite:

**This Repository**:
```
Gonçalves, J. (2025). K-mer Marker Error Resilience Analysis for ONT Sequencing.
GitHub: https://github.com/jacgonisa/kmer-marker-error-resilience
Version 2.1 (FDR calculations corrected)
```

**K-mer Counter (KMC)**:
```
Deorowicz S, Kokot M, Grabowski S, Debudaj-Grabysz A.
KMC 2: Fast and resource-frugal k-mer counting.
Bioinformatics. 2015;31(10):1569-1576.
```

---

## Troubleshooting

### "No error resilience data found!"
→ Ensure `realistic_k*_100k_error_resilience_stats.csv` files exist in `final_results/`

### "KMC databases not found"
→ Check that `../../03-cenhapmers/k*/` directories contain `.kmc_pre` and `.kmc_suf` files

### "kmc_tools not found"
→ Install KMC: `conda install -c bioconda kmc` or download from [KMC GitHub](https://github.com/refresh-bio/KMC)

### Plot API errors / too large
→ This version creates modular, focused plots (~500KB each) instead of massive single plots

### Import errors
→ Ensure all dependencies installed: `pip install -r requirements.txt`

### Wrong FDR values (>10%)
→ You may be using an old version! Update to v2.1+ which fixes the FDR calculation bug

---

## Version History

### Version 2.1 (December 4, 2025) - CRITICAL BUG FIX
- **FIXED**: Major bug in absolute FDR calculation
  - Previous versions incorrectly calculated: `(% with errors) × (conditional FDR)`
  - Correct calculation: `(% with errors) × (% that become FP)`
  - Result: FDR values now ~100x lower (0.14% vs 17.7% for k=21)
- **IMPROVED**: Visualization of Panel B - now shows cross-contamination rate directly
- **IMPROVED**: Panel D - now compares k=21 vs k=41 error fate (novel vs cross-contamination)
- **ADDED**: Comprehensive FDR explanation document (FDR_EXPLAINED.md)
- **ADDED**: Stacked bar visualization showing ~99% novel vs <1% cross-contamination
- **CLARIFIED**: Biological interpretation throughout documentation

### Version 2.0 (December 2025) - Modular Analysis Suite
- Comprehensive k-mer evaluation across 5 sizes
- 100 databases analyzed (20 databases × 5 k-mer sizes)
- Realistic error simulation (1% per-base, ONT R9-like)
- Publication-quality figures (300 DPI PNG + PDF vectors)
- Multiple visualization perspectives

---

## Repository Structure

```
.
├── README.md                              # This file
├── FDR_EXPLAINED.md                       # Detailed FDR explanation ⭐ READ THIS!
├── RESULTS_FIGURES.md                     # Visual documentation of all plots
├── ANALYSIS_SUMMARY.md                    # Results summary
├── LICENSE                                # MIT License
├── requirements.txt                       # Python dependencies
├── .gitignore                             # Git ignore rules
│
├── Main Analysis Scripts (run 1-5, 6-7 optional):
├── 01_marker_availability.py              # K-mer counts and density
├── 02_cross_contamination.py              # FDR and cross-contamination ⭐
├── 03_error_resilience.py                 # K-mer retention under errors
├── 04_final_recommendation.py             # Weighted scoring
├── 05_radar_comparison.py                 # Pentagon radar charts
├── 06_coverage_loss_analysis.py           # Coverage per Mb (supplementary)
├── 07_objective_comparison.py             # Trade-offs (supplementary)
│
├── comprehensive_marker_evaluation.py     # Full evaluation framework
├── run_all_analyses.sh                    # Master execution script (runs 1-5)
│
├── final_results/                         # Output directory
│   ├── *.png/pdf                          # Figures (7 pairs)
│   ├── *.csv                              # Data tables
│   └── *.csv.gz                           # Detailed event data
│
├── old_analyses/                          # Previous analysis versions
├── old_scripts/                           # Legacy scripts
├── scripts/                               # Helper scripts
└── Documentation:
    ├── TERMINOLOGY_CLARIFICATION.md       # Key concepts
    ├── QUICK_START.md                     # Quick usage guide
    ├── FINAL_SUMMARY.md                   # Objective comparison summary
    └── GITHUB_SETUP.md                    # Repository setup guide
```

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

**Author**: Jacob Gonçalves
**GitHub**: [@jacgonisa](https://github.com/jacgonisa)

For questions, issues, or collaboration inquiries, please open an issue on GitHub.

---

## Acknowledgments

- KMC development team for the excellent k-mer counting tool
- Oxford Nanopore Technologies for sequencing technology advancements
- Arabidopsis research community

---

**Last Updated**: December 4, 2025
**Version**: 2.1 (FDR calculations corrected + improved visualizations)
