###############################################################################
# 18_venn_overflow_examples.R
# Three approaches to the "small section number overflow" problem in Venn
# diagrams where one set is much smaller than the other.
#
# PROBLEM: When circle areas are proportional to set size (e.g., -RA = 80
# proteins, +RA = 450 proteins), the small "-RA only" section is so tiny
# that the count number overflows it.
#
# THREE SOLUTIONS:
#   Example 1: External labels with leader lines (RECOMMENDED)
#     - Preserves proportional circles
#     - Small region's count placed OUTSIDE with a thin connecting line
#     - Uses ext.text=TRUE with raised ext.percent threshold
#
#   Example 2: Equal (non-proportional) circles
#     - Both circles same size, counts shown by numbers not area
#     - All labels fit comfortably inside
#     - Trade-off: less visually honest about set sizes
#
#   Example 3: ggVennDiagram with boxed labels
#     - Modern ggplot2-based approach
#     - White-background boxes behind count labels
#     - Readable even when overlapping circle boundaries
#
# Uses ACTUAL pipeline data (real protein counts from RA effect comparison).
# Colors from GLOBAL_COLORS for consistency across all plots.
#
# Usage:
#   make venn-overflow-examples
###############################################################################
cat("\n=========================================\n")
cat(" Venn Overflow Examples (3 variants)\n")
cat("=========================================\n\n")

library(VennDiagram)
library(grid)

experiments <- load_all_experiments()

# Get actual RA comparison data
exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"

if (!exp_base %in% names(experiments) || !exp_ra %in% names(experiments)) {
  cat("ERROR: Required experiments not found.\n")
  cat("Need:", exp_base, "and", exp_ra, "\n")
  quit(status = 1)
}

set_a <- get_significant_genes(experiments[[exp_base]])
set_b <- get_significant_genes(experiments[[exp_ra]])

area_a <- length(set_a)
area_b <- length(set_b)
overlap <- length(intersect(set_a, set_b))

cat(sprintf("Set A (-RA): %d proteins\n", area_a))
cat(sprintf("Set B (+RA): %d proteins\n", area_b))
cat(sprintf("Overlap:     %d proteins\n\n", overlap))

commit_hash <- get_git_hash()
fill_a <- GLOBAL_COLORS[["venn_a_only"]]
fill_b <- GLOBAL_COLORS[["venn_b_only"]]

title_str <- "TRIP4 without vs with Retinoic Acid"

# Helper: save viewport to PNG + PDF
save_venn <- function(vp, prefix, subtitle) {
  safe_name <- sanitize_filename(prefix)
  versioned <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

  grDevices::png(png_path, width = 7, height = 6, units = "in", res = 300)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text(subtitle, 0.5, 0.04, gp = gpar(fontsize = 9, fontface = "italic"))
  popViewport()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 7, height = 6)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title_str, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  grid.text(subtitle, 0.5, 0.04, gp = gpar(fontsize = 9, fontface = "italic"))
  popViewport()
  grDevices::dev.off()

  cat(sprintf("  Saved: %s\n", basename(png_path)))
}

# =====================================================================
# EXAMPLE 1: External labels with leader lines (RECOMMENDED)
# =====================================================================
# Preserves proportional circles. Small region's count placed OUTSIDE
# the circle with a thin connecting line (leader line).
# ext.percent raised to 0.25 forces small regions to use external text.
# =====================================================================
cat("[Example 1] External labels with leader lines...\n")

vp1 <- draw.pairwise.venn(
  area1     = area_a,
  area2     = area_b,
  cross.area = overlap,
  category  = c("- RA", "+ RA"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  # KEY: External text for small regions
  ext.text  = TRUE,
  ext.percent = c(0.05, 0.25, 0.25),  # Raise threshold for small regions
  ext.pos   = c(0, 0),
  ext.dist  = c(0.1, 0.15),
  ext.line.lty = "solid",
  ext.length = c(0.95, 0.95),
  # Per-region font sizing: [area1-only, area2-only, intersection]
  cex       = c(1.4, 1.4, 1.0),
  fontface  = "bold",
  fontfamily = "sans",
  cat.cex   = 1.6,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  col       = "transparent",
  cat.pos   = c(220, 140),
  cat.dist  = c(0.06, 0.06),
  margin    = 0.08,
  ind       = FALSE
)
save_venn(vp1, "example1_venn_ext_text",
          "Example 1: Proportional circles + external labels with leader lines")

# =====================================================================
# EXAMPLE 2: Equal (non-proportional) circles
# =====================================================================
# Both circles same size. Counts shown by numbers, not circle area.
# All labels fit comfortably. Trade-off: less visually honest about sizes.
# =====================================================================
cat("\n[Example 2] Equal circles (non-proportional)...\n")

vp2 <- draw.pairwise.venn(
  area1     = area_a,
  area2     = area_b,
  cross.area = overlap,
  category  = c("- RA", "+ RA"),
  fill      = c(fill_a, fill_b),
  alpha     = rep(0.5, 2),
  # KEY: Disable proportional scaling
  euler.d   = FALSE,
  scaled    = FALSE,
  ext.text  = TRUE,
  ext.percent = c(0.05, 0.05, 0.05),
  cex       = 1.5,
  fontface  = "bold",
  fontfamily = "sans",
  cat.cex   = 1.6,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  col       = "transparent",
  cat.pos   = c(220, 140),
  cat.dist  = c(0.06, 0.06),
  margin    = 0.08,
  ind       = FALSE
)
save_venn(vp2, "example2_venn_equal_circles",
          "Example 2: Equal circles (non-proportional) — counts by number, not area")

# =====================================================================
# EXAMPLE 3: ggVennDiagram with boxed labels (if installed)
# =====================================================================
# Modern ggplot2-based approach. White-background boxes behind count labels.
# Readable even when overlapping circle boundaries.
# =====================================================================
cat("\n[Example 3] ggVennDiagram with boxed labels...\n")

has_ggvenn <- requireNamespace("ggVennDiagram", quietly = TRUE)

if (has_ggvenn) {
  library(ggVennDiagram)
  library(ggplot2)

  sets_list <- list(
    "- RA" = set_a,
    "+ RA" = set_b
  )

  p3 <- ggVennDiagram(
    sets_list,
    label       = "count",
    label_geom  = "label",    # Boxed background labels
    label_alpha = 0.8,
    label_size  = 5,
    edge_size   = 1,
    set_color   = "black",
    set_size    = 7
  ) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4292C6") +
    labs(title = title_str) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

  safe_name <- sanitize_filename("example3_venn_boxed_labels")
  versioned <- paste0(safe_name, "_", commit_hash)
  save_figure(p3, "example3_venn_boxed_labels", width = 7, height = 6)
} else {
  cat("  ggVennDiagram not installed — skipping Example 3.\n")
  cat("  Install with: install.packages('ggVennDiagram')\n")
  cat("  Example 3 would show boxed labels with white background.\n")
}

cat("\n=========================================\n")
cat(" Venn overflow examples complete!\n")
cat("=========================================\n")
cat("\nThree variants saved to output/figures/:\n")
cat("  example1 — Proportional circles + external labels (RECOMMENDED)\n")
cat("  example2 — Equal circles (non-proportional)\n")
cat("  example3 — Boxed labels (if ggVennDiagram installed)\n")
cat("\nRun: make open-venn-overflow-examples\n")
