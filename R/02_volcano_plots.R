###############################################################################
# 02_volcano_plots.R
# Lydia-style volcano plots with multi-category highlighting.
#
# Categories overlaid on each volcano:
#   - ia: Known interactors (from known_interactors.txt)
#   - flagMulti / flagOnce: Flag IP hits (from Flag IP data)
#   - CRAC: RNA interactome hits (if CRAC data available)
#   - gp / dhx / ddx / LARPs: Gene families
#   - high: High-confidence hits (top combined significance)
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/02_volcano_plots.R")
###############################################################################

cat("\n=========================================\n")
cat(" Volcano Plot Generation (Lydia-style)\n")
cat("=========================================\n\n")

# ---- Load known interactors ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

# ---- Load all experiment data ----
cat("Loading experiment data...\n")
experiments <- load_all_experiments()

# ---- Build Flag IP hit lists for overlay ----
# Proteins significant in 2+ Flag conditions = flagMulti
# Proteins significant in exactly 1 Flag condition = flagOnce
flag_exp_names <- grep("^flag_", names(experiments), value = TRUE)
if (length(flag_exp_names) >= 2) {
  cat("\nBuilding Flag IP hit categories...\n")
  flag_sig_lists <- lapply(experiments[flag_exp_names], get_significant_genes)

  flag_all_genes <- unlist(flag_sig_lists)
  flag_counts <- table(flag_all_genes)
  flag_multi <- names(flag_counts)[flag_counts >= 2]
  flag_once <- names(flag_counts)[flag_counts == 1]
  cat(sprintf("  Flag IP hits: %d multi-condition, %d single-condition\n",
              length(flag_multi), length(flag_once)))
} else {
  flag_multi <- character(0)
  flag_once <- character(0)
}

# ---- Assign gene family categories ----
assign_gene_family <- function(genes) {
  fam <- rep(NA_character_, length(genes))
  fam[grepl("^DHX", genes)] <- "dhx"
  fam[grepl("^DDX", genes)] <- "ddx"
  fam[grepl("^LARP", genes)] <- "LARPs"
  fam[genes %in% GENE_FAMILIES$GPATCH] <- "gp"
  return(fam)
}

# =====================================================================
# Helper: Lydia-style volcano with category highlighting
# =====================================================================
plot_lydia_volcano <- function(df, title, known_ia = NULL,
                               flag_m = NULL, flag_o = NULL,
                               show_gene_families = TRUE,
                               highlight_high = TRUE) {

  toPlot <- df
  # Start with significance as category
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )

  # Layer 1: Gene families (only for significant proteins)
  if (show_gene_families) {
    fam <- assign_gene_family(toPlot$gene)
    fam_mask <- !is.na(fam) & toPlot$category == "TRUE"
    toPlot$category[fam_mask] <- fam[fam_mask]
  }

  # Layer 2: Flag IP hits
  if (!is.null(flag_m) && length(flag_m) > 0) {
    toPlot$category[toPlot$gene %in% flag_m] <- "flagMulti"
  }
  if (!is.null(flag_o) && length(flag_o) > 0) {
    toPlot$category[toPlot$gene %in% flag_o] <- "flagOnce"
  }

  # Layer 3: Known interactors (highest priority)
  if (!is.null(known_ia) && length(known_ia) > 0) {
    toPlot$category[toPlot$gene %in% known_ia] <- "ia"
  }

  # Layer 4: High-confidence hits
  if (highlight_high) {
    toPlot$category[!is.na(toPlot$padj) &
                    ((toPlot$log2FC > 2 & -log10(toPlot$padj) > 6) |
                     (toPlot$log2FC > 7 & -log10(toPlot$padj) > 2))] <- "high"
  }

  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  # Determine which points to label
  label_cats <- c("ia", "gp", "dhx", "ddx", "LARPs", "flagMulti", "flagOnce", "high")
  label_data <- toPlot[toPlot$category %in% label_cats, ]

  # Build plot
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    ggplot2::scale_color_manual(
      values = CATEGORY_COLORS,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = title,
      color = "Category"
    )

  # Add labels for highlighted points
  if (nrow(label_data) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 30,
      show.legend = FALSE
    )
  }

  # Threshold lines
  p <- p +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50")

  # Theme (Lydia's compact style)
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(colour = "black", size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(colour = "black", size = 8),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10)
    )

  return(p)
}

# =====================================================================
# Helper: Overlay two experiments on one volcano
# =====================================================================
plot_volcano_overlay <- function(df1, name1, df2, name2, label_genes, title) {
  df1$experiment <- name1
  df2$experiment <- name2
  df1$sig <- df1$padj < P_VALUE_CUTOFF & abs(df1$log2FC) > LOG2FC_CUTOFF
  df2$sig <- df2$padj < P_VALUE_CUTOFF & abs(df1$log2FC) > LOG2FC_CUTOFF

  combined <- rbind(
    df1[, c("gene", "log2FC", "padj", "experiment", "sig")],
    df2[, c("gene", "log2FC", "padj", "experiment", "sig")]
  )
  combined$neglog10p <- -log10(combined$padj)

  # Label only from one dataset (deduplicate)
  label_data <- combined[combined$gene %in% label_genes & combined$experiment == name1, ]

  p <- ggplot2::ggplot(combined, ggplot2::aes(x = log2FC, y = neglog10p, color = experiment)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 25, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = EXPERIMENT_COLORS) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = title, color = "Experiment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  return(p)
}

# =====================================================================
# MAIN: Generate all volcano plots
# =====================================================================

# ---- 1. Lydia-style volcano for each TurboID experiment ----
cat("\n[1/4] TurboID volcano plots (Lydia-style)...\n")
turbo_names <- grep("^turbo_", names(experiments), value = TRUE)
for (name in turbo_names) {
  p <- plot_lydia_volcano(
    experiments[[name]], name,
    known_ia = known_interactors,
    flag_m = flag_multi, flag_o = flag_once,
    show_gene_families = TRUE
  )
  save_figure(p, paste0("volcano_", name), width = 7, height = 5)
}

# ---- 2. Lydia-style volcano for each Flag IP experiment ----
cat("\n[2/4] Flag IP volcano plots (Lydia-style)...\n")
flag_names <- grep("^flag_", names(experiments), value = TRUE)
for (name in flag_names) {
  p <- plot_lydia_volcano(
    experiments[[name]], name,
    known_ia = known_interactors,
    show_gene_families = TRUE
  )
  save_figure(p, paste0("volcano_", name), width = 7, height = 5)
}

# ---- 3. Overlay: TurboID TRIP4 vs Flag IP C-Flag ----
cat("\n[3/4] TurboID vs Flag IP overlay...\n")
turbo_main <- "turbo_trip4_vs_wt"
flag_main <- "flag_cflag_vs_ctrl"
if (turbo_main %in% names(experiments) && flag_main %in% names(experiments)) {
  p <- plot_volcano_overlay(
    experiments[[turbo_main]], "TurboID_TRIP4",
    experiments[[flag_main]], "FlagIP_C",
    label_genes = known_interactors,
    title = "TurboID (TRIP4 vs WT) vs Flag IP (C-Flag vs Ctrl)"
  )
  save_figure(p, "volcano_overlay_turboid_flag", width = 10, height = 7)
}

# ---- 4. Overlay: TRIP4 vs TRIP4+RA (hormone effect) ----
cat("\n[4/4] RA treatment comparison overlay...\n")
turbo_base <- "turbo_trip4_vs_wt"
turbo_ra <- "turbo_RA_vs_wt"
if (turbo_base %in% names(experiments) && turbo_ra %in% names(experiments)) {
  p <- plot_volcano_overlay(
    experiments[[turbo_base]], "TurboID_TRIP4",
    experiments[[turbo_ra]], "TurboID_TRIP4_RA",
    label_genes = known_interactors,
    title = "TurboID: TRIP4 vs TRIP4 + Retinoic Acid"
  )
  save_figure(p, "volcano_overlay_RA_effect", width = 10, height = 7)
}

cat("\n=========================================\n")
cat(" Volcano plots complete!\n")
cat(sprintf(" Output: %s/\n", FIGURE_DIR))
cat("=========================================\n")
