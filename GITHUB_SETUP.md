# GitHub Repository Setup Guide

This guide will walk you through setting up this project as a GitHub repository.

## Prerequisites

- GitHub account (username: `jacgonisa`)
- Git installed on your system
- GitHub CLI (optional, but recommended)

## Step-by-Step Setup

### Step 1: Create GitHub Repository

**Option A: Using GitHub Website (Recommended)**

1. Go to https://github.com/new
2. Fill in the repository details:
   - **Repository name**: `kmer-marker-error-resilience` (or your preferred name)
   - **Description**: "Comprehensive k-mer marker analysis framework for evaluating error resilience in ONT sequencing data"
   - **Visibility**: Public (recommended for research) or Private
   - **DO NOT** initialize with README, .gitignore, or license (we already have these)
3. Click "Create repository"

**Option B: Using GitHub CLI**

```bash
# Install gh if not already installed
# Ubuntu/Debian: sudo apt install gh
# macOS: brew install gh

# Login to GitHub
gh auth login

# Create repository
gh repo create kmer-marker-error-resilience --public --description "Comprehensive k-mer marker analysis framework for evaluating error resilience in ONT sequencing data"
```

### Step 2: Initialize Local Git Repository

From the `ERROR_RESILIENCE_ANALYSIS` directory:

```bash
# Initialize git repository
git init

# Add all files to staging
git add .

# Check what will be committed
git status

# Create initial commit
git commit -m "Initial commit: K-mer marker error resilience analysis framework

- Complete analysis suite for k-mer sizes 21, 25, 31, 35, 41
- FDR (False Discovery Rate) calculations
- Error resilience simulation (1% per-base error rate)
- Publication-quality visualizations
- Comprehensive documentation
- Requirements and setup instructions"
```

### Step 3: Connect to GitHub Remote

Replace `YOUR_REPO_NAME` with the actual repository name you created:

```bash
# Add remote repository
git remote add origin https://github.com/jacgonisa/YOUR_REPO_NAME.git

# Verify remote
git remote -v

# Push to GitHub
git branch -M main
git push -u origin main
```

**If you chose a different default branch name (e.g., 'master'):**

```bash
git push -u origin master
```

### Step 4: Verify Upload

1. Go to `https://github.com/jacgonisa/YOUR_REPO_NAME`
2. Verify that all files are present
3. Check that README.md displays correctly

---

## What's Included in This Repository

### Essential Files
- ✅ `README.md` - Comprehensive documentation with installation, usage, and methodology
- ✅ `LICENSE` - MIT License
- ✅ `requirements.txt` - Python dependencies
- ✅ `.gitignore` - Excludes unnecessary files (Python cache, large binary files)

### Analysis Scripts (7 main scripts)
- ✅ `01_marker_availability.py` - K-mer counts and density
- ✅ `02_cross_contamination.py` - **FDR analysis (UPDATED!)**
- ✅ `03_error_resilience.py` - Read retention under errors
- ✅ `04_final_recommendation.py` - Integrated scoring
- ✅ `05_radar_comparison.py` - Pentagon radar charts
- ✅ `06_coverage_loss_analysis.py` - Coverage impact analysis
- ✅ `07_objective_comparison.py` - Trade-off visualization

### Supporting Files
- ✅ `run_all_analyses.sh` - Master execution script
- ✅ `comprehensive_marker_evaluation.py` - Full evaluation framework
- ✅ Documentation files (ANALYSIS_SUMMARY.md, FINAL_SUMMARY.md, etc.)

### Output Directory
- ✅ `final_results/` - Contains generated figures and data
  - Note: Large gzipped event files (`.csv.gz`) are included
  - Total size: ~48 MB

---

## Important Notes

### Large Files

The repository includes some larger files (~48MB total):
- `final_results/realistic_k*_events.csv.gz` (3-11 MB each)

**If GitHub rejects large files**, you have two options:

**Option 1: Keep Large Files (Use Git LFS)**

```bash
# Install Git LFS
git lfs install

# Track large files
git lfs track "final_results/*.csv.gz"

# Add .gitattributes
git add .gitattributes

# Commit and push
git commit -m "Add Git LFS tracking for large data files"
git push
```

**Option 2: Exclude Large Files**

Edit `.gitignore` to add:
```
final_results/*.csv.gz
```

Then:
```bash
git rm --cached final_results/*.csv.gz
git commit -m "Remove large data files from repository"
git push
```

### Updating FDR Calculations

**IMPORTANT**: The script `02_cross_contamination.py` has been updated to use **False Discovery Rate (FDR)** instead of False Positive Rate (FPR):

- **Old metric**: FPR = False Positives
- **New metric**: FDR = FP / (FP + TP)

This is scientifically more accurate for this type of analysis.

---

## After Repository is Created

### Add Repository Topics (Tags)

On GitHub repository page:
1. Click the settings gear icon next to "About"
2. Add topics:
   - `bioinformatics`
   - `genomics`
   - `kmer-analysis`
   - `oxford-nanopore`
   - `error-correction`
   - `arabidopsis`
   - `python`
   - `data-analysis`

### Create Release (Optional)

```bash
# Tag the initial version
git tag -a v1.0.0 -m "Initial release: K-mer marker error resilience analysis"
git push origin v1.0.0
```

Or use GitHub's Release feature:
1. Go to repository → Releases → Create a new release
2. Tag: `v1.0.0`
3. Title: `Initial Release - K-mer Marker Analysis v1.0.0`
4. Description: Copy from README key findings section

### Add Zenodo DOI (For Academic Citation)

1. Link repository to Zenodo: https://zenodo.org/account/settings/github/
2. Create a release on GitHub
3. Zenodo will automatically create a DOI
4. Add DOI badge to README

---

## Maintenance

### Updating the Repository

```bash
# Make changes to files
git add .
git commit -m "Description of changes"
git push
```

### Adding New Features

```bash
# Create a new branch
git checkout -b feature/new-analysis

# Make changes and commit
git add .
git commit -m "Add new analysis feature"

# Push branch
git push -u origin feature/new-analysis

# Create pull request on GitHub to merge into main
```

---

## Recommended Repository Settings

### Branch Protection (Optional but Recommended)

1. Go to Settings → Branches
2. Add rule for `main` branch:
   - ✅ Require pull request reviews before merging
   - ✅ Require status checks to pass before merging

### GitHub Actions (Optional - Future Enhancement)

Consider adding automated testing:
```yaml
# .github/workflows/test.yml
name: Run Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.8
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run analysis scripts
        run: ./run_all_analyses.sh
```

---

## Troubleshooting

### Permission Denied

If you get permission errors:
```bash
# Set up SSH key or use Personal Access Token
# See: https://docs.github.com/en/authentication
```

### Large File Errors

```
error: file too large
```

Solution: Use Git LFS (see "Large Files" section above)

### Merge Conflicts

If you edit files on GitHub and locally:
```bash
git pull --rebase origin main
# Resolve conflicts
git add .
git rebase --continue
git push
```

---

## Quick Reference

```bash
# Check status
git status

# View commit history
git log --oneline

# View remote repository
git remote -v

# Pull latest changes
git pull

# Create and switch to new branch
git checkout -b new-feature

# Switch back to main
git checkout main

# Delete branch
git branch -d branch-name
```

---

## Next Steps After Upload

1. ✅ Verify all files uploaded correctly
2. ✅ Test the README displays properly
3. ✅ Add repository topics/tags
4. ✅ Share repository URL with collaborators
5. ✅ Consider creating a release with DOI
6. ✅ Add to your CV/publications if applicable

---

## Support

If you encounter issues:
1. Check GitHub documentation: https://docs.github.com
2. GitHub community forum: https://github.community
3. Contact: Open an issue in the repository

---

**Repository URL** (after creation):
`https://github.com/jacgonisa/YOUR_REPO_NAME`

**Suggested Repository Names**:
- `kmer-marker-error-resilience`
- `kmer-marker-analysis`
- `ont-kmer-error-analysis`
- `arabidopsis-kmer-markers`

---

**Good luck with your repository! 🚀**
