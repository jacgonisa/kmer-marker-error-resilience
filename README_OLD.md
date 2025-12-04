# K-mer Marker Comprehensive Analysis Suite

This directory contains a complete analysis framework for evaluating k-mer markers across different k-mer sizes (21, 25, 31, 35, 41).

## 📁 Directory Structure

```
ERROR_RESILIENCE_ANALYSIS/
├── README.md                          # This file
├── run_all_analyses.sh               # Master script to run all analyses
│
├── Analysis Scripts (run in order):
├── 01_marker_availability.py         # Marker counts and density
├── 02_cross_contamination.py         # False positive risk analysis
├── 03_error_resilience.py           # Read retention under errors
├── 04_final_recommendation.py       # Integrated decision matrix
│
├── final_results/                   # Output directory
│   ├── realistic_k*_100k_error_resilience_stats.csv  # Raw data
│   ├── 01_marker_availability.png/pdf
│   ├── 02_cross_contamination.png/pdf
│   ├── 03_error_resilience.png/pdf
│   └── 04_final_recommendation.png/pdf
│
├── old_analyses/                    # Previous analysis versions
└── scripts/                         # Helper scripts
```

## 🎯 What Each Analysis Does

### 1. Marker Availability Analysis (`01_marker_availability.py`)
**Purpose**: Evaluate marker quantity and distribution

**Metrics**:
- Total k-mer count per database
- Marker density (k-mers per Mb)
- Distribution uniformity (Coefficient of Variation)
- Coverage consistency across chromosomes/regions

**Output**: 4-panel figure showing:
- A: Average marker count (ARMS vs CEN)
- B: Marker density per megabase
- C: Distribution of k-mer counts
- D: Uniformity across databases

### 2. Cross-Contamination Analysis (`02_cross_contamination.py`)
**Purpose**: Assess false positive risk from sequencing errors

**Metrics**:
- Absolute false positive rate (% of ALL k-mers)
- Conditional false positive rate (% of k-mers WITH errors)
- Per-database vulnerability
- ARMS vs CEN comparison

**Key Question**: When sequencing errors occur, do they cause misclassification?

**Output**: 4-panel figure showing:
- A: Absolute false positive rates ⭐ MOST IMPORTANT
- B: Conditional rates (of k-mers with errors)
- C: Heatmap by database
- D: Database vulnerability ranking

### 3. Error Resilience Analysis (`03_error_resilience.py`)
**Purpose**: Quantify read retention under realistic ONT sequencing

**Simulation**: 1% per-base error rate (ONT R9-like)

**Metrics**:
- Percentage of usable reads after errors
- Error impact by k-mer size
- Error tolerance (ARMS vs CEN)
- Practical impact with 1M reads

**Output**: 6-panel figure showing:
- A: Read retention by k-mer size ⭐ KEY METRIC
- B: K-mers affected by errors
- C: Error tolerance by region
- D: Distribution of error tolerance
- E: Error outcomes (novel vs false positive)
- F: Practical impact visualization

### 4. Final Recommendation (`04_final_recommendation.py`)
**Purpose**: Integrate all factors to recommend optimal k-mer size

**Scoring System** (weighted):
- Read Retention: 40% (highest priority)
- Specificity: 30% (false positive rate)
- Uniformity: 15% (marker distribution)
- Availability: 15% (total marker count)

**Output**: Multi-panel decision matrix showing:
- A: Overall weighted scores ⭐ FINAL RECOMMENDATION
- B: Scoring criteria explanation
- C: Radar chart (k=21 vs k=41)
- D: Individual component scores
- E: Complete decision table

## 🚀 How to Run

### Run All Analyses (Recommended)
```bash
./run_all_analyses.sh
```

This will:
1. Create `final_results/` directory
2. Run all 4 analyses in sequence
3. Generate PNG (300 dpi) and PDF (vector) figures
4. Save summary statistics to CSV files
5. Print comprehensive report

**Expected runtime**: 2-5 minutes (depending on k-mer database size)

### Run Individual Analyses
```bash
# Individual scripts can be run separately:
python3 01_marker_availability.py
python3 02_cross_contamination.py
python3 03_error_resilience.py
python3 04_final_recommendation.py
```

## 📊 Key Findings Summary

Based on comprehensive analysis across 5 k-mer sizes:

### ✅ Recommended: **k=21**

**Reasons**:
1. **Best Read Retention**: ~81% usable reads vs ~66% for k=41 (15% improvement!)
2. **Excellent Specificity**: <0.2% false positive rate
3. **Lowest Error Impact**: Only 19% of k-mers affected by errors (vs 34% for k=41)
4. **Practical Advantage**: ~148,000 more usable reads per million ONT reads

### 📈 Detailed Comparison

| K-mer | Usable Reads | False Pos. Rate | K-mers Affected | Recommendation |
|-------|--------------|-----------------|-----------------|----------------|
| k=21  | 81.10%      | 0.14%          | 19.0%          | ⭐ **BEST**   |
| k=25  | 77.93%      | 0.13%          | 22.2%          | Good          |
| k=31  | 73.36%      | 0.11%          | 27.0%          | Acceptable    |
| k=35  | 70.47%      | 0.10%          | 30.1%          | Not optimal   |
| k=41  | 66.34%      | 0.09%          | 33.7%          | Poor retention|

### 🔬 Scientific Rationale

**Why shorter k-mers are better for ONT**:
1. Each base has ~1% error probability
2. Longer k-mers = more bases = higher probability of ≥1 error
3. P(error) = 1 - (0.99)^k
   - k=21: P(error) = 19%
   - k=41: P(error) = 34%
4. Errors usually create "novel" k-mers → read loss
5. Shorter k-mers minimize read loss while maintaining specificity

## 📁 Output Files

### Figures (PNG + PDF)
- `01_marker_availability.[png|pdf]` - Marker quantity and density
- `02_cross_contamination.[png|pdf]` - False positive risk
- `03_error_resilience.[png|pdf]` - Read retention analysis
- `04_final_recommendation.[png|pdf]` - Integrated decision matrix

### Data Files (CSV)
- `marker_availability_summary.csv` - Marker counts and density stats
- `comprehensive_scores.csv` - Final scoring matrix
- `realistic_k{21,25,31,35,41}_100k_error_resilience_stats.csv` - Raw error simulation data

## 🔧 Requirements

### Python Packages
```bash
pip install pandas matplotlib seaborn numpy
```

### External Tools
- `kmc_tools` (for counting k-mers in databases)

### Input Data
- KMC databases in `../../03-cenhapmers/k{21,25,31,35,41}/`
- Error resilience statistics in `final_results/realistic_k*_100k_error_resilience_stats.csv`

## 📝 Methodology

### Error Simulation
- **Error model**: Per-base 1% error rate (Bernoulli distribution)
- **Sample size**: 100,000 k-mers per database
- **Error types**: Random base substitutions (A↔T↔G↔C)
- **Classification**: Each mutated k-mer tested against all databases

### Scoring System
All metrics normalized to 0-100 scale:
- **Higher is better**: Read retention, marker count
- **Lower is better**: False positive rate, CV (uniformity)

Weighted combination:
```
Overall Score = 0.40 × ReadRetention + 0.30 × Specificity +
                0.15 × Uniformity + 0.15 × Availability
```

## 🐛 Troubleshooting

### "No error resilience data found!"
→ Ensure `realistic_k*_100k_error_resilience_stats.csv` files exist in `final_results/`

### "KMC databases not found"
→ Check that `../../03-cenhapmers/k*/` directories contain `.kmc_pre` and `.kmc_suf` files

### "kmc_tools not found"
→ Install KMC: `conda install -c bioconda kmc` or download from https://github.com/refresh-bio/KMC

### Plot too large / API error
→ This version creates smaller, focused plots (~500KB each) instead of one massive plot

## 📚 Citation

If using this analysis framework, please cite:
- KMC: Deorowicz et al. (2015) Bioinformatics
- Error simulation methodology: This study

## 📧 Contact

For questions or issues with this analysis suite, please consult your bioinformatics team or the original researcher.

---

**Last updated**: 2025-12-04
**Version**: 2.0 (modular analysis suite)
