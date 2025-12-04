# Quick Start Guide

## Run All Analyses (Recommended)

```bash
cd ERROR_RESILIENCE_ANALYSIS
./run_all_analyses.sh
```

**Runtime**: ~10-15 seconds

**Output**: 8 publication-quality figures (PNG + PDF) in `final_results/`

---

## Run Individual Analyses

```bash
# 1. Marker availability and density
python3 01_marker_availability.py

# 2. Cross-contamination risk
python3 02_cross_contamination.py

# 3. Error resilience and read retention
python3 03_error_resilience.py

# 4. Final recommendation
python3 04_final_recommendation.py
```

---

## Key Results (TL;DR)

### ⭐ Recommendation: **Use k=21**

**Why?**
- ✅ Best read retention: 81% vs 66% for k=41 (+148k reads per 1M)
- ✅ Excellent specificity: 0.16% false positive rate
- ✅ Optimal for ONT sequencing (1% error rate)

**Practical Impact with 1M ONT reads:**
- k=21: ~811,000 usable reads
- k=41: ~663,000 usable reads
- **You lose 15% of your data with k=41!**

---

## Generated Files

### Figures (in `final_results/`)
- `01_marker_availability.png` - Marker counts and density
- `02_cross_contamination.png` - False positive risk
- `03_error_resilience.png` - Read retention ⭐ KEY FIGURE
- `04_final_recommendation.png` - Decision matrix

### Data (in `final_results/`)
- `marker_availability_summary.csv` - K-mer counts
- `comprehensive_scores.csv` - Final scores
- `realistic_k*_error_resilience_stats.csv` - Error simulation results

---

## Troubleshooting

**"No error resilience data found!"**
→ Run the error simulation first (see main README.md)

**"KMC databases not found"**
→ Check that `../../../03-cenhapmers/k*/` exists

**"kmc_tools not found"**
→ Install KMC: `conda install -c bioconda kmc`

---

## Documentation

- **Full details**: See `README.md`
- **Results summary**: See `ANALYSIS_SUMMARY.md`
- **This guide**: `QUICK_START.md`

---

**Questions?** Check the full README.md or contact your bioinformatics team.
