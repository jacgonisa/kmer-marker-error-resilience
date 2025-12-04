# Terminology Clarification: Coverage Loss, Not Read Loss

## ⚠️ IMPORTANT CORRECTION

The earlier analyses used the term **"read retention"** which could be misleading. Here's the correct understanding:

---

## What Actually Happens with Sequencing Errors

### ❌ INCORRECT: "You lose reads"

### ✅ CORRECT: "You lose coverage/classification ability per megabase"

---

## The Correct Model

### 1. Reads Contain Many Overlapping K-mers

A single 150bp ONT read contains:
- **~129 overlapping 21-mers** (positions 1-21, 2-22, 3-23, ..., 130-150)
- **~109 overlapping 41-mers** (positions 1-41, 2-42, 3-43, ..., 110-150)

### 2. Errors Affect Local K-mer Windows

If there's a sequencing error at position 50:
- **k=21**: Affects k-mers at positions 30-50 (21 consecutive k-mers)
- **k=41**: Affects k-mers at positions 10-50 (41 consecutive k-mers)

The rest of the read's k-mers are still fine!

### 3. You Lose Classification Ability Locally

```
Read: ═══════════[ERROR]═══════════════════════
      ↓          ↓        ↓
      OK    21 lost k-mers    OK
                or
      OK    41 lost k-mers    OK
```

You **don't discard the entire read**. You lose the ability to **classify that genomic region** (that megabase).

---

## The Correct Metric: Error-Affected K-mers per Mb

### Definition

**Coverage Loss** = (% k-mers with errors) × (k-mer density per Mb)

This tells you: **How many classification attempts fail per megabase of genome?**

### Real Numbers (ARMS markers)

| K-mer | Total K-mers/Mb | Error Rate | Lost K-mers/Mb | Coverage Retained |
|-------|-----------------|------------|----------------|-------------------|
| k=21  | 115.4 million   | 19.0%      | **21.9 million**  | 93.5 million (81%) |
| k=41  | 187.6 million   | 33.7%      | **63.2 million**  | 124.4 million (66%) |

**Key Finding**: k=41 loses **2.9× MORE coverage per Mb** than k=21!

---

## Why Longer K-mers Lose More Coverage

### Mathematical Explanation

With 1% per-base error rate:
```
P(k-mer contains ≥1 error) = 1 - (0.99)^k

k=21: P(error) = 19.0%
k=41: P(error) = 33.7%
```

### Physical Explanation

1. **More bases** = more chances for an error
2. Each error affects **k consecutive k-mers** (overlapping window)
3. Longer k-mers = **larger affected window** per error
4. Result: **More classification attempts fail per megabase**

---

## Practical Impact

### Example: 10 Mb Genomic Region

**With k=21 ARMS markers:**
- Total classification attempts: 115.4M × 10 = 1,154 million
- Lost to errors: 21.9M × 10 = 219 million (19%)
- Usable: 935 million attempts ✓
- **Effective coverage: 8.1 Mb**

**With k=41 ARMS markers:**
- Total classification attempts: 187.6M × 10 = 1,876 million
- Lost to errors: 63.2M × 10 = 632 million (33.7%)
- Usable: 1,244 million attempts
- **Effective coverage: 6.6 Mb**

**You lose 1.5 Mb more effective coverage with k=41!**

---

## Updated Terminology

| Old Term (Misleading) | New Term (Correct) |
|-----------------------|-------------------|
| "Read retention" | "Coverage retention" or "K-mer retention" |
| "Reads lost" | "Coverage lost per Mb" |
| "Usable reads %" | "Usable k-mers %" or "Coverage retention %" |
| "X% of reads remain" | "X% of classification attempts succeed" |

---

## Why This Matters

### Biological Interpretation

You're not losing sequencing **reads** (the raw data), you're losing:
1. **Classification ability** in local genomic windows
2. **Genome coverage** for marker-based analysis
3. **Resolution** for detecting crossover events

### Practical Consequences

With k=41 vs k=21:
- ❌ You don't sequence fewer reads
- ❌ You don't throw away reads
- ✅ You lose 2.9× more classification events per Mb
- ✅ You need deeper sequencing to compensate
- ✅ You have lower effective genome coverage

---

## Conclusion

The **recommendation remains the same** (k=21 is better), but the reasoning is more accurate:

**k=21 is superior because it maintains higher genome coverage per megabase under ONT error rates, requiring less sequencing depth to achieve the same effective coverage as k=41.**

---

## New Figure: 06_coverage_loss_analysis

This figure shows the **corrected metric**:
- Panel A: Lost k-mers per Mb (THE KEY METRIC)
- Panel B: Total k-mer density per Mb
- Panel C: Error rate (% k-mers affected)
- Panel D: Remaining classification ability
- Panel E: Coverage composition (usable vs lost)
- Panel F: Key insights and interpretation

**File**: `final_results/06_coverage_loss_analysis.png`

---

**Last Updated**: 2025-12-04
**Status**: ✅ Terminology corrected, metrics updated
