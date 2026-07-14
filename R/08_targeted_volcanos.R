###############################################################################
# 08_targeted_volcanos.R
# Custom volcano plots per researcher's specific requirements.
#
# THREE PLOTS with DIFFERENT rules per plot type:
#
# PLOT 1: TRIP4 TurboID vs Wild Type (MAIN experiment)
#   - All circles, same size
#   - ASCC complex (ASCC1, ASCC2, ASCC3, TRIP4): blue, LABELED
#   - Known interactors (excluding ASCC): green, LABELED
#   - WT-enriched proteins (negative log2FC): gray, no labels
#   - Other TRIP4-enriched proteins (positive, sig): vermillion, no labels
#   - Category names in legend: ASCC complex, Known interactors, etc.
#
# PLOT 2: BK467 TRIP4 without RA vs with RA
#   - All circles, same size
#   - Enriched proteins: orange, LABEL ALL of them
#   - Title: proper µM symbol
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

source("R/00_theme.R")

experiments <- load_all_experiments()
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
known_ia_excl_core <- setdiff(known_interactors, ASCC_CORE)

# ---- Helper: RA comparison volcano (Plots 2 and 3) ----
# Label top 20 significant genes (10 up + 10 down by combined score).
# Uses a white background behind labels for readability on any background color.
make_ra_volcano <- function(df, title, n_top = 60) {
  toPlot <- df
  toPlot$neglog10p <- -log10(toPlot$padj)

  # Use config thresholds: padj < P_VALUE_CUTOFF, |log2FC| > LOG2FC_CUTOFF
  # THREE categories per researcher feedback (Jul 12 voice message):
  #   enriched_up (positive log2FC) = orange "Enriched with RA"
  #   enriched_dn (negative log2FC) = blue "Enriched without RA"
  #   nonsig = grey "Not significant"
  toPlot$category <- "nonsig"
  toPlot$category[toPlot$log2FC > LOG2FC_CUTOFF &
                  toPlot$padj < P_VALUE_CUTOFF] <- "enriched_up"
  toPlot$category[toPlot$log2FC < -LOG2FC_CUTOFF &
                  toPlot$padj < P_VALUE_CUTOFF] <- "enriched_dn"
  toPlot$category <- factor(toPlot$category,
    levels = c("enriched_up", "enriched_dn", "nonsig"))

  # Find top N proteins by combined significance on each side
  sig_up <- toPlot[toPlot$category == "enriched_up", ]
  sig_up <- sig_up[!is.na(sig_up$log2FC) & !is.na(sig_up$neglog10p), ]
  sig_up$combined_score <- abs(sig_up$log2FC) * sig_up$neglog10p

  sig_dn <- toPlot[toPlot$category == "enriched_dn", ]
  sig_dn <- sig_dn[!is.na(sig_dn$log2FC) & !is.na(sig_dn$neglog10p), ]
  sig_dn$combined_score <- abs(sig_dn$log2FC) * sig_dn$neglog10p

  # Top N/2 up-regulated (positive log2FC) — enriched with RA
  top_up <- sig_up[order(-sig_up$combined_score), ]
  top_up <- head(top_up, n_top / 2)

  # Top N/2 down-regulated (negative log2FC) — enriched without RA
  top_dn <- sig_dn[order(-sig_dn$combined_score), ]
  top_dn <- head(top_dn, n_top / 2)

  label_data <- rbind(top_up, top_dn)

  # Aruna's requested colors: orange=with RA, blue=without RA, gray=nonsig
  RA_COLORS <- c(
    "enriched_up" = "#D55E00",   # Vermillion Orange — enriched with RA
    "enriched_dn" = "#0072B2",   # Deep Navy — enriched without RA
    "nonsig"      = "#D0D0D0"    # Light grey — not significant
  )

  RA_LABELS <- c(
    "enriched_up" = "Enriched with RA",
    "enriched_dn" = "Enriched without RA",
    "nonsig"      = "Not significant"
  )

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                             color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = RA_COLORS,
      labels = RA_LABELS,
      name = NULL, drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)  # Legend dot size
      )
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 50, show.legend = FALSE,
      bg.color = "white", bg.r = 0.15  # White background for readability
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(P_VALUE_CUTOFF),
      linetype = "dashed", color = "grey50", linewidth = 0.3
    ) +
    ggplot2::geom_vline(
      xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
      linetype = "dashed", color = "grey50", linewidth = 0.3
    ) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +
    theme_poster(font_size = 16) +
    ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3)
    )

  return(p)
}

# ---- Helper: Main TRIP4 vs WT volcano ----
# ASCC complex and known interactors highlighted and LABELED.
# WT-enriched = gray, other TRIP4-enriched = vermillion, not sig = light gray.
make_main_volcano <- function(df, title) {
  toPlot <- df
  toPlot$neglog10p <- -log10(toPlot$padj)

  # Start: not significant
  toPlot$category <- "nonsig"

  # WT-enriched (significant, negative log2FC)
  wt_idx <- toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC < -LOG2FC_CUTOFF
  toPlot$category[wt_idx] <- "enriched_dn"

  # TRIP4-enriched (significant, positive log2FC)
  trip4_idx <- toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC > LOG2FC_CUTOFF
  toPlot$category[trip4_idx] <- "enriched_up"

  # Known interactors (excluding ASCC core)
  toPlot$category[toPlot$gene %in% known_ia_excl_core] <- "known_ia"

  # ASCC complex (overwrites everything else)
  toPlot$category[toPlot$gene %in% ASCC_CORE] <- "ascc_core"

  # Factor with ordering (last = drawn on top)
  toPlot$category <- factor(toPlot$category,
    levels = c("ascc_core", "known_ia", "enriched_up", "enriched_dn", "nonsig"))

  MAIN_COLORS <- c(
    "ascc_core"   = GLOBAL_COLORS[["ascc_core"]],     # Blue — ASCC complex
    "known_ia"    = GLOBAL_COLORS[["known_ia"]],       # Green — known interactors
    "enriched_up" = GLOBAL_COLORS[["enriched_up"]],    # Vermillion — enriched in TRIP4
    "enriched_dn" = GLOBAL_COLORS[["enriched_dn"]],    # Gray — enriched in WT
    "nonsig"      = GLOBAL_COLORS[["nonsig"]]          # Light gray
  )

  MAIN_LABELS <- c(
    "ascc_core"   = "ASCC complex",
    "known_ia"    = "Known interactors",
    "enriched_up" = "Enriched in TRIP4",
    "enriched_dn" = "Enriched in WT",
    "nonsig"      = "Not significant"
  )

  # Label only ASCC complex and known interactors
  # BUT: don't label proteins enriched in WT (negative log2FC) — per researcher request
  label_data <- toPlot[toPlot$category %in% c("ascc_core", "known_ia") &
                       !is.na(toPlot$log2FC) & toPlot$log2FC > 0, ]

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                             color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = MAIN_COLORS,
      labels = MAIN_LABELS,
      name = NULL, drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)  # Legend dot size
      )
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 30, show.legend = FALSE,
      bg.color = "white", bg.r = 0.15
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(P_VALUE_CUTOFF),
      linetype = "dashed", color = "grey50", linewidth = 0.3
    ) +
    ggplot2::geom_vline(
      xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
      linetype = "dashed", color = "grey50", linewidth = 0.3
    ) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +
    theme_poster(font_size = 16) +
    ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3)
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
  # Use proper mu symbol \u00b5M for micromolar
  title2 <- sprintf("TRIP4 without RA vs with RA (0.2 %sM)", "\u00b5")
  p2 <- make_ra_volcano(
    experiments[[exp2]],
    title = title2
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
  title3 <- sprintf("TRIP4 without RA vs with RA (0.4 %sM)", "\u00b5")
  p3 <- make_ra_volcano(
    experiments[[exp3]],
    title = title3
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
