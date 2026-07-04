###############################################################################
# 08_targeted_volcanos.R
# Custom volcano plots per researcher's specific requirements.
#
# THREE PLOTS with DIFFERENT rules per plot type:
#
# PLOT 1: TRIP4 TurboID vs Wild Type (MAIN experiment)
#   - All circles, same size
#   - WT-enriched proteins (negative log2FC): gray
#   - TRIP4-enriched proteins (positive log2FC, sig): dark orange
#   - Known interactors: colored but NOT labeled
#   - ASCC core: colored but NOT labeled
#   - Label: none for this plot
#
# PLOT 2: BK467 TRIP4 without RA vs with RA
#   - All circles, same size
#   - Enriched proteins: orange
#   - Do NOT label known interactors or ASCC core
#   - Label: top 20 proteins by combined significance (10 up + 10 down)
#
# PLOT 3: BK504 TRIP4 without RA vs with RA
#   - Same rules as Plot 2
#
# Usage:
#   make targeted-volcano
###############################################################################

cat("\n=========================================\n")
cat(" Targeted Volcano Plots\n")
cat("=========================================\n\n")

experiments <- load_all_experiments()
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
known_ia_excl_core <- setdiff(known_interactors, ASCC_CORE)

# Colors (all categories use circles now — researcher said no triangles)
TARGETED_COLORS <- c(
  "ascc_core"    = "#0072B2",   # Blue
  "known_ia"     = "#009E73",   # Green
  "enriched_up"  = "#D55E00",   # Vermillion (orange) — enriched in TRIP4
  "enriched_dn"  = "#B0B0B0",   # Gray — enriched in WT control
  "nonsig"       = "#D0D0D0"    # Light gray — not significant
)

# ---- Helper: RA comparison volcano (Plots 2 and 3) ----
# For RA comparison plots: label top 20 proteins by significance,
# do NOT label known interactors or ASCC core specifically.
make_ra_volcano <- function(df, title, n_top = 20) {
  toPlot <- df
  toPlot$neglog10p <- -log10(toPlot$padj)

  # Categories: enriched or nonsig (no separate known_ia/ascc categories)
  toPlot$category <- "nonsig"
  toPlot$category[abs(toPlot$log2FC) >= 1 & toPlot$neglog10p > 1] <- "enriched_up"
  toPlot$category <- factor(toPlot$category, levels = c("enriched_up", "nonsig"))

  # Find top N proteins by combined significance on each side
  # "Combined significance" = product of fold change magnitude and -log10(p)
  sig_only <- toPlot[abs(toPlot$log2FC) >= 1 & toPlot$neglog10p > 1, ]
  sig_only$combined_score <- abs(sig_only$log2FC) * sig_only$neglog10p

  # Top 10 up-regulated (positive log2FC)
  top_up <- sig_only[sig_only$log2FC > 0, ]
  top_up <- top_up[order(-top_up$combined_score), ]
  top_up <- head(top_up, n_top / 2)

  # Top 10 down-regulated (negative log2FC)
  top_dn <- sig_only[sig_only$log2FC < 0, ]
  top_dn <- top_dn[order(-top_dn$combined_score), ]
  top_dn <- head(top_dn, n_top / 2)

  label_data <- rbind(top_up, top_dn)

  RA_COLORS <- c(
    "enriched_up" = "#D55E00",
    "nonsig"      = "#D0D0D0"
  )

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = RA_COLORS,
      labels = c("enriched_up" = "Significant", "nonsig" = "Not significant"),
      name = NULL
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 25, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text = ggplot2::element_text(colour = "black", size = 8),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

# ---- Helper: Main TRIP4 vs WT volcano ----
# WT-enriched proteins = gray, only TRIP4-enriched = orange
# No labels on this plot
make_main_volcano <- function(df, title) {
  toPlot <- df
  toPlot$neglog10p <- -log10(toPlot$padj)

  # Three categories: enriched_up (TRIP4, orange), enriched_dn (WT, gray), nonsig
  toPlot$category <- "nonsig"
  toPlot$category[toPlot$log2FC >= 1 & toPlot$neglog10p > 1] <- "enriched_up"
  toPlot$category[toPlot$log2FC <= -1 & toPlot$neglog10p > 1] <- "enriched_dn"
  toPlot$category <- factor(toPlot$category,
                            levels = c("enriched_up", "enriched_dn", "nonsig"))

  MAIN_COLORS <- c(
    "enriched_up" = "#D55E00",   # Orange — enriched in TRIP4
    "enriched_dn" = "#B0B0B0",   # Gray — enriched in WT control
    "nonsig"      = "#D0D0D0"    # Light gray
  )

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = MAIN_COLORS,
      labels = c(
        "enriched_up" = "Enriched in TRIP4",
        "enriched_dn" = "Enriched in WT",
        "nonsig"      = "Not significant"
      ),
      name = NULL
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text = ggplot2::element_text(colour = "black", size = 8),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

# =====================================================================
# PLOT 1: TRIP4 TurboID vs Wild Type
# =====================================================================
cat("[1/3] TRIP4 TurboID vs Wild Type...\n")

exp1 <- "BK467_TRIP4_vs_BK467_WT"
if (exp1 %in% names(experiments)) {
  p1 <- make_main_volcano(
    experiments[[exp1]],
    title = "TRIP4 TurboID vs Wild Type"
  )
  save_figure(p1, "targeted_volcano_BK467_TRIP4_vs_WT",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found:", exp1, "\n")
}

# =====================================================================
# PLOT 2: BK467 TRIP4 without RA vs with RA
# =====================================================================
cat("\n[2/3] TRIP4 without RA vs with RA...\n")

exp2 <- "BK467_TRIP4_RA02_vs_BK467_TRIP4"
if (exp2 %in% names(experiments)) {
  p2 <- make_ra_volcano(
    experiments[[exp2]],
    title = "TRIP4 without RA vs with RA (0.2 uM)",
    n_top = 20
  )
  save_figure(p2, "targeted_volcano_BK467_RA_effect",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found:", exp2, "\n")
}

# =====================================================================
# PLOT 3: BK504 TRIP4 without RA vs with RA
# =====================================================================
cat("\n[3/3] BK504 TRIP4 without RA vs with RA...\n")

exp3 <- "BK504_TRIP4_RA04_vs_BK504_TRIP4"
if (exp3 %in% names(experiments)) {
  p3 <- make_ra_volcano(
    experiments[[exp3]],
    title = "TRIP4 without RA vs with RA (0.4 uM)",
    n_top = 20
  )
  save_figure(p3, "targeted_volcano_BK504_RA_effect",
              width = 8, height = 6)
} else {
  cat("  WARNING: Experiment not found:", exp3, "\n")
}

cat("\n=========================================\n")
cat(" Targeted volcano plots complete!\n")
cat(sprintf(" Output: %s/\n", FIGURE_DIR))
cat("=========================================\n")
