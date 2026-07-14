#!/usr/bin/env bash
##############################################################################
# run_poster_steps.sh
#
# Runs each poster analysis step with caching:
#   - If output file already exists in output/figures/, skip the step
#   - Otherwise, run the R step
#
# Usage: bash poster/run_poster_steps.sh [--force]
#   --force  Ignore cache, regenerate all steps
#
##############################################################################
set -euo pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
FIGURE_DIR="output/figures"
FORCE=0

if [[ "${1:-}" == "--force" ]]; then
    FORCE=1
fi

# Step name → key output pattern (if this file exists, step is cached)
# Array of "step_name|output_pattern" strings
STEPS=(
    "targeted_volcanos|targeted_volcano_BK504_RA_effect"
    "flagip_volcano|flagip_overlap_volcano_BK467_TRIP4_vs_WT"
    "lydia_network_volcano|lydia_network_volcano"
    "flagip_validated_go|flagip_GO_validated_any_BP_dotplot"
    "network_go|network_go_comparison_BP"
    "ra_common|RA_common_enriched"
    "targeted_go|targeted_GO_RA_shared_core_MF_dotplot"
)

SKIPPED=0
FAILED=0

for entry in "${STEPS[@]}"; do
    step="${entry%%|*}"
    pattern="${entry##*|}"

    # Check cache: does output exist?
    if [ "$FORCE" -eq 0 ]; then
        found=$(find "$FIGURE_DIR" -name "${pattern}*.pdf" -print -quit 2>/dev/null || true)
        if [ -z "$found" ]; then
            found=$(find "$FIGURE_DIR" -name "${pattern}*.png" -print -quit 2>/dev/null || true)
        fi
        if [ -n "$found" ] && [ -f "$found" ]; then
            echo "  [CACHED] $step — output exists: $(basename "$found")"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    fi

    # Run the step
    echo ""
    echo "========================================="
    echo " Running: $step"
    echo "========================================="
    if "$RSCRIPT" R/run_step.R "$step"; then
        echo "  [OK] $step"
    else
        echo "  [FAILED] $step (exit $?)"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "========================================="
if [ "$SKIPPED" -gt 0 ]; then
    echo " $SKIPPED step(s) cached, $FAILED failed"
else
    echo " All steps completed, $FAILED failed"
fi
if [ "$FORCE" -eq 1 ]; then
    echo " (--force: all steps regenerated)"
fi
echo "========================================="

exit $FAILED
