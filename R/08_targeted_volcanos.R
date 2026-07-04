###############################################################################
# 08_targeted_volcanos.R
# Custom volcano plots per researcher's specific requirements.
#
# THREE PLOTS:
#   1. BK467: TRIP4 TurboID vs Wild Type
#   2. BK467: TRIP4 without RA vs TRIP4 with RA (hormone effect)
#   3. BK504: TRIP4 without RA vs TRIP4 with RA (hormone effect, replicate)
#
# COLORING SCHEME (for each plot):
#   - Light gray:  |log2FC| < 1 (did not reach 2-fold change threshold)
#   - Light pink:  log2FC >= 1 AND -log10(padj) > 1 (highly enriched)
#   - Light purple: Known interactors (the 27 genes in known_interactors.txt,
#                   excluding ASCC1-3/TRIP4 which get their own color)
#   - Dark purple: ASCC1, ASCC2, ASCC3, TRIP4 (the ASCC complex + bait protein)
#
# AXES:
#   X-axis: log2 fold change (logFC from the mass spec pipeline)
#   Y-axis: -log10(adjusted p-value)
#
# LABELS:
#   All highlighted proteins (pink, light purple, dark purple) are labeled
#   with their gene name using ggrepel (non-overlapping text labels).
#
# Usage:
#   Rscript R/run_step.R targeted_volcanos
#   or:
#   make targeted-volcano
###############################################################################

cat("\n=========================================\n")
cat(" Targeted Volcano Plots (Researcher Spec)\n")
cat("=========================================\n\n")

# ---- Load data ----
# Load all experiments, then pick the specific ones we need.
experiments <- load_all_experiments()

# Load known interactors list (27 genes from literature/previous experiments)
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

# ---- Define the highlight groups ----

# DARK PURPLE: The ASCC complex + TRIP4 itself.
# These are the bait protein (TRIP4) and its known complex partners (ASCC1-3).
# They get the highest visual priority because they are the core of the study.
ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")

# LIGHT PURPLE: Known interactors from the literature list,
# EXCLUDING the ASCC core (those get dark purple).
# setdiff() removes ASCC_CORE genes from the known interactors list
# so each gene only gets ONE color (dark purple wins for ASCC1-3/TRIP4).
known_ia_excl_core <- setdiff(known_interactors, ASCC_CORE)

# ---- Custom colors ----
# The researcher specified these colors:
#   Light gray:  not significant
#   Light pink:  highly enriched
#   Light purple: known interactors
#   Dark purple: ASCC core
TARGETED_COLORS <- c(
  "ascc_core"    = "#4A148C",   # Dark purple — ASCC1-3 + TRIP4
  "known_ia"     = "#CE93D8",   # Light purple — other known interactors
  "enriched"     = "#FFCDD2",   # Light pink — highly enriched proteins
  "nonsig"       = "#E0E0E0"    # Light gray — below threshold
)

# ---- Helper function: build one targeted volcano plot ----
# This function takes ONE experiment's data frame and creates a volcano
# plot with the specific coloring scheme requested.
#
# Parameters:
#   df:          Data frame with columns gene, log2FC, padj
#   title:       Plot title
#   fc_cutoff:   Fold change threshold (default 1.0 = 2-fold change)
#   pval_line:   -log10(padj) threshold (default 1.0 = padj < 0.1)
#
# Returns: a ggplot2 plot object
make_targeted_volcano <- function(df, title, fc_cutoff = 1.0, pval_line = 1.0) {

  # Create a working copy of the data
  toPlot <- df

  # ---- Step 1: Calculate -log10(adjusted p-value) for the Y-axis ----
  # This transforms very small p-values into readable numbers.
  # Example: padj = 0.00001 → -log10(0.00001) = 5.0 (high on Y-axis)
  #          padj = 0.5      → -log10(0.5)    = 0.3 (low on Y-axis)
  toPlot$neglog10p <- -log10(toPlot$padj)

  # ---- Step 2: Assign each protein to a category ----
  # We start with everyone as "nonsig" (light gray), then upgrade
  # specific proteins to higher-priority categories.
  #
  # Priority (highest wins):
  #   1. ascc_core (dark purple) — ASCC1-3 + TRIP4
  #   2. known_ia  (light purple) — known interactors from literature
  #   3. enriched  (light pink) — log2FC >= 1 AND -log10(padj) > 1
  #   4. nonsig    (light gray) — everything else

  # Start: everyone is light gray
  toPlot$category <- "nonsig"

  # Layer 1: Highly enriched proteins (light pink)
  # Condition: absolute fold change >= 1 (2-fold) AND -log10(padj) > 1
  # abs() ensures we catch BOTH up-regulated (positive) and down-regulated (negative)
  toPlot$category[abs(toPlot$log2FC) >= fc_cutoff & toPlot$neglog10p > pval_line] <- "enriched"

  # Layer 2: Known interactors (light purple) — overrides "enriched"
  # %in% checks if each gene name is in the known_ia_excl_core list
  toPlot$category[toPlot$gene %in% known_ia_excl_core] <- "known_ia"

  # Layer 3: ASCC core (dark purple) — highest priority, overrides everything
  toPlot$category[toPlot$gene %in% ASCC_CORE] <- "ascc_core"

  # Convert to factor with ordered levels (controls legend order)
  # factor() creates a categorical variable. levels defines the allowed
  # values and their order. This ensures the legend shows categories in
  # priority order (ascc_core first, nonsig last).
  toPlot$category <- factor(toPlot$category,
                            levels = c("ascc_core", "known_ia", "enriched", "nonsig"))

  # ---- Step 3: Decide which points to label ----
  # Label everything EXCEPT nonsig (gray points).
  # These are the proteins we want to identify by gene name.
  label_data <- toPlot[toPlot$category != "nonsig", ]

  # ---- Step 4: Build the plot ----
  # ggplot() initializes the plot. aes() maps data columns to visual properties:
  #   x = log2FC → horizontal position
  #   y = neglog10p → vertical position
  #   color = category → dot color determined by the category we assigned
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +

    # geom_point() draws one circle per protein.
    # alpha controls transparency (0=invisible, 1=fully opaque).
    # We use higher alpha for highlighted points, lower for gray background.
    ggplot2::geom_point(data = toPlot[toPlot$category == "nonsig", ],
                        alpha = 0.4, size = 1.0) +
    ggplot2::geom_point(data = toPlot[toPlot$category != "nonsig", ],
                        alpha = 0.8, size = 2.0) +

    # scale_color_manual() assigns our custom colors to each category.
    # drop = FALSE ensures all categories appear in the legend even if
    # some have zero proteins (e.g., if no known interactors are present).
    ggplot2::scale_color_manual(
      values = TARGETED_COLORS,
      labels = c(
        "ascc_core" = "ASCC1-3 + TRIP4",
        "known_ia"  = "Known interactors",
        "enriched"  = "Highly enriched",
        "nonsig"    = "Not significant"
      ),
      name = NULL,
      drop = FALSE
    ) +

    # Add text labels for all highlighted proteins using ggrepel.
    # geom_text_repel() automatically pushes labels apart so they don't overlap.
    #   max.overlaps: if more than 30 labels would overlap a single point,
    #   skip the extra ones (prevents unreadable label clusters).
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.8, fontface = "bold",
      max.overlaps = 30,
      show.legend = FALSE
    ) +

    # Dashed threshold lines:
    #   Horizontal line at y = pval_line (the -log10(padj) threshold)
    #   Vertical lines at x = ±fc_cutoff (the fold change thresholds)
    ggplot2::geom_hline(yintercept = pval_line, linetype = "dashed",
                        color = "grey50", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed",
                        color = "grey50", linewidth = 0.3) +

    # Axis labels using mathematical notation
    # expression() lets us use subscripts and italics in axis labels
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +

    # theme_bw(): clean black-and-white theme with grid lines
    # theme(): customize appearance
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text = ggplot2::element_text(colour = "black", size = 8),
      axis.title = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

# =====================================================================
# PLOT 1: BK467 — TRIP4 TurboID vs Wild Type
# =====================================================================
# This is the MAIN experiment: TRIP4-TurboID (with biotin) vs Wild Type (biotin only).
# Proteins enriched here are candidate TRIP4 interaction partners.
#
# The CSV file "BK467_TRIP4_vs_BK467_WT_diffEx_minProb.csv" contains
# the differential expression results: log2FC and padj for ~3000 proteins.
cat("[1/3] BK467: TRIP4 TurboID vs Wild Type...\n")

exp1 <- "BK467_TRIP4_vs_BK467_WT"
if (exp1 %in% names(experiments)) {
  p1 <- make_targeted_volcano(
    experiments[[exp1]],
    title = "BK467: TRIP4 TurboID vs Wild Type"
  )
  save_figure(p1, "targeted_volcano_BK467_TRIP4_vs_WT",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found: ", exp1, "\n")
  cat("  Available experiments:\n")
  for (n in names(experiments)) cat("    - ", n, "\n")
}

# =====================================================================
# PLOT 2: BK467 — TRIP4 without RA vs TRIP4 with RA
# =====================================================================
# This compares TRIP4-TurboID cells WITHOUT retinoic acid vs WITH RA (0.2µM).
# Proteins enriched here are those whose association with TRIP4 CHANGES
# upon hormone treatment — these may be RA-dependent interactors.
cat("\n[2/3] BK467: TRIP4 without RA vs TRIP4 with RA...\n")

exp2 <- "BK467_TRIP4_RA02_vs_BK467_TRIP4"
if (exp2 %in% names(experiments)) {
  p2 <- make_targeted_volcano(
    experiments[[exp2]],
    title = "BK467: TRIP4 without RA vs TRIP4 with RA (0.2 µM)"
  )
  save_figure(p2, "targeted_volcano_BK467_RA_effect",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found: ", exp2, "\n")
}

# =====================================================================
# PLOT 3: BK504 — TRIP4 without RA vs TRIP4 with RA
# =====================================================================
# Same comparison as Plot 2, but in the BK504 batch (biological replicate).
# Comparing BK467 and BK504 results shows which RA effects are reproducible.
cat("\n[3/3] BK504: TRIP4 without RA vs TRIP4 with RA...\n")

exp3 <- "BK504_TRIP4_RA04_vs_BK504_TRIP4"
if (exp3 %in% names(experiments)) {
  p3 <- make_targeted_volcano(
    experiments[[exp3]],
    title = "BK504: TRIP4 without RA vs TRIP4 with RA (0.4 µM)"
  )
  save_figure(p3, "targeted_volcano_BK504_RA_effect",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found: ", exp3, "\n")
}

cat("\n=========================================\n")
cat(" Targeted volcano plots complete!\n")
cat(sprintf(" Output: %s/\n", FIGURE_DIR))
cat("=========================================\n")
