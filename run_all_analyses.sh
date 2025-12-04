#!/bin/bash
#
# Master script to run all k-mer analysis plots in sequence
#

echo "════════════════════════════════════════════════════════════════════════"
echo "  COMPREHENSIVE K-MER MARKER ANALYSIS SUITE"
echo "════════════════════════════════════════════════════════════════════════"
echo ""

# Create output directory
mkdir -p final_results

# Track timing
START_TIME=$(date +%s)

echo "📊 [1/5] Analyzing marker availability and density..."
echo "────────────────────────────────────────────────────────────────────────"
python3 01_marker_availability.py
if [ $? -eq 0 ]; then
    echo "✓ Marker availability analysis complete!"
else
    echo "✗ ERROR in marker availability analysis!"
    exit 1
fi
echo ""

echo "🎯 [2/5] Analyzing cross-contamination risk..."
echo "────────────────────────────────────────────────────────────────────────"
python3 02_cross_contamination.py
if [ $? -eq 0 ]; then
    echo "✓ Cross-contamination analysis complete!"
else
    echo "✗ ERROR in cross-contamination analysis!"
    exit 1
fi
echo ""

echo "🧬 [3/5] Analyzing error resilience and read retention..."
echo "────────────────────────────────────────────────────────────────────────"
python3 03_error_resilience.py
if [ $? -eq 0 ]; then
    echo "✓ Error resilience analysis complete!"
else
    echo "✗ ERROR in error resilience analysis!"
    exit 1
fi
echo ""

echo "🏆 [4/5] Generating final recommendation..."
echo "────────────────────────────────────────────────────────────────────────"
python3 04_final_recommendation.py
if [ $? -eq 0 ]; then
    echo "✓ Final recommendation complete!"
else
    echo "✗ ERROR in final recommendation!"
    exit 1
fi
echo ""

echo "🌟 [5/5] Creating pentagon radar comparison..."
echo "────────────────────────────────────────────────────────────────────────"
python3 05_radar_comparison.py
if [ $? -eq 0 ]; then
    echo "✓ Pentagon radar comparison complete!"
else
    echo "✗ ERROR in radar comparison!"
    exit 1
fi
echo ""

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "════════════════════════════════════════════════════════════════════════"
echo "  ✓ ALL ANALYSES COMPLETE!"
echo "════════════════════════════════════════════════════════════════════════"
echo ""
echo "Runtime: ${RUNTIME} seconds"
echo ""
echo "Generated files in final_results/:"
echo "  • 01_marker_availability.png/pdf"
echo "  • 02_cross_contamination.png/pdf"
echo "  • 03_error_resilience.png/pdf"
echo "  • 04_final_recommendation.png/pdf"
echo "  • 05_radar_comparison.png/pdf ⭐ NEW!"
echo "  • marker_availability_summary.csv"
echo "  • comprehensive_scores.csv"
echo ""
echo "════════════════════════════════════════════════════════════════════════"
