#!/bin/bash
##############################################################################
# collect_review_figures.sh
# Finds the 6 poster figures from output/figures/ (including subdirs with
# spaces like "update before lydia" and "Final volcano plots") and copies
# them to poster/figures/ with clean names for the LaTeX poster.
#
# Prefers PDF (vector) over PNG (raster). Falls back to PNG if no PDF.
#
# Usage: bash poster/collect_review_figures.sh
##############################################################################
set -euo pipefail

FIGURE_SRC="output/figures"
POSTER_DIR="poster/figures"

mkdir -p "$POSTER_DIR"

# ---- Helper: find and copy a figure ----
# Args: $1 = search pattern, $2 = clean output name (without extension)
copy_figure() {
    local pattern="$1"
    local clean_name="$2"

    # Try PDF first (vector format — best for LaTeX)
    local found_pdf
    found_pdf=$(find "$FIGURE_SRC" -name "${pattern}*.pdf" -print -quit 2>/dev/null || true)

    if [ -n "$found_pdf" ] && [ -f "$found_pdf" ]; then
        cp "$found_pdf" "${POSTER_DIR}/${clean_name}.pdf"
        echo "  Copied: $(basename "$found_pdf") -> ${clean_name}.pdf"
        return 0
    fi

    # Fall back to PNG
    local found_png
    found_png=$(find "$FIGURE_SRC" -name "${pattern}*.png" -print -quit 2>/dev/null || true)

    if [ -n "$found_png" ] && [ -f "$found_png" ]; then
        cp "$found_png" "${POSTER_DIR}/${clean_name}.png"
        echo "  Copied: $(basename "$found_png") -> ${clean_name}.png"
        echo "  NOTE: Only PNG found (no PDF). Consider regenerating with save_figure()."
        return 0
    fi

    echo "  WARNING: Could not find '${pattern}*' in ${FIGURE_SRC}/"
    echo "    Searched recursively in all subdirectories."
    return 1
}

echo "============================================"
echo " Collecting review poster figures..."
echo "============================================"

# ---- Dot plots (3) ----
echo ""
echo "Dot plots:"
copy_figure "flagip_GO_validated_any_BP_dotplot" "dotplot_flagip_validated_BP"
copy_figure "network_go_comparison_BP"            "dotplot_network_go_BP"
copy_figure "targeted_GO_RA_shared_core_MF_dotplot" "dotplot_RA_shared_MF"

# ---- Volcano plots (3) ----
echo ""
echo "Volcano plots:"
copy_figure "flagip_overlap_volcano_BK467_TRIP4_vs_WT" "volcano_flagip_overlap"
copy_figure "lydia_network_volcano"                     "volcano_lydia_network"
copy_figure "targeted_volcano_BK504_RA_effect"          "volcano_RA_BK504"

echo ""
echo "============================================"
echo " Collection complete. Figures in poster/figures/"
echo " Now compile: cd poster && pdflatex poster_review.tex"
echo "============================================"
