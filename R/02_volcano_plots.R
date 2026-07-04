###############################################################################
# 02_volcano_plots.R
# Generates publication-quality volcano plots for each experiment.
# Also creates overlay plots comparing TurboID vs Flag co-IP.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/02_volcano_plots.R")
###############################################################################

cat("\n=========================================\n")
cat(" Volcano Plot Generation\n")
cat("=========================================\n\n")

# ---- Load known interactors ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

# ---- Helper: EnhancedVolcano for a single experiment ----
plot_enhanced_volcano <- function(df, experiment_name, label_genes = NULL) {
  cat(sprintf("\n  [%s] Generating EnhancedVolcano...\n", experiment_name))

  # Handle NAs
  df$padj[is.na(df$padj)] <- 1
  df$log2FC[is.na(df$log2FC)] <- 0

  select_genes <- intersect(label_genes, df$gene)

  p <- EnhancedVolcano::EnhancedVolcano(
    toptable    = df,
    lab         = df$gene,
    x           = "log2FC",
    y           = "padj",
    pCutoff     = P_VALUE_CUTOFF,
    FCcutoff    = LOG2FC_CUTOFF,
    selectLab   = if (length(select_genes) > 0) select_genes else NULL,
    labSize     = 3.0,
    labCol      = "black",
    labFace     = "bold",
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    colConnectors = "grey50",
    col = c("grey70", "grey50", "#4DBBD5", "#E64B35"),
    colAlpha    = 0.6,
    pointSize   = 1.5,
    cutoffLineType = "dashed",
    cutoffLineCol  = "grey60",
    cutoffLineWidth = 0.5,
    legendPosition = "right",
    title       = experiment_name,
    subtitle    = NULL,
    caption     = NULL
  )

  return(p)
}

# ---- Helper: Overlay two experiments on the same volcano ----
plot_volcano_overlay <- function(df1, name1, df2, name2, label_genes, title) {
  cat(sprintf("\n  [OVERLAY] %s vs %s...\n", name1, name2))

  # Prepare data
  df1$neglog10p <- -log10(df1$padj)
  df2$neglog10p <- -log10(df2$padj)
  df1$experiment <- name1
  df2$experiment <- name2

  combined <- rbind(
    df1[, c("gene", "log2FC", "neglog10p", "experiment")],
    df2[, c("gene", "log2FC", "neglog10p", "experiment")]
  )

  # Label data — deduplicate: show one label per gene
  label_data <- combined[combined$gene %in% label_genes, ]
  if (nrow(label_data) > 0) {
    # Prefer the row with higher significance for labeling
    label_data <- label_data[order(label_data$neglog10p, decreasing = TRUE), ]
    label_data <- label_data[!duplicated(label_data$gene), ]
  }

  p <- ggplot2::ggplot(combined,
    ggplot2::aes(x = log2FC, y = neglog10p, color = experiment)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.2) +
    # Emphasize labeled points
    ggplot2::geom_point(
      data = combined[combined$gene %in% label_genes, ],
      ggplot2::aes(fill = experiment),
      shape = 21, size = 3, alpha = 0.8, stroke = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 3, fontface = "bold",
      max.overlaps = 25,
      show.legend = FALSE,
      color = "black"
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = EXPERIMENT_COLORS[c(name1, name2)]) +
    ggplot2::scale_fill_manual(values = EXPERIMENT_COLORS[c(name1, name2)]) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title,
      color = "Experiment",
      fill = "Experiment"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  return(p)
}

# =====================================================================
# MAIN: Run all volcano plots
# =====================================================================

# ---- Load all experiment data ----
cat("Loading experiment data...\n")

experiments <- list()
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  name <- tools::file_path_sans_ext(basename(f))
  # Skip non-data files
  if (name == "known_interactors") next
  experiments[[name]] <- load_proteomics_csv(f)
}

if (length(experiments) == 0) {
  stop("No CSV files found in ", DATA_DIR)
}

# ---- 1. Individual EnhancedVolcano for each experiment ----
cat("\n[1/3] Generating individual volcano plots...\n")

for (name in names(experiments)) {
  df <- experiments[[name]]
  p <- plot_enhanced_volcano(df, name, label_genes = known_interactors)
  save_figure(p, paste0("volcano_", name), width = 8, height = 6)
}

# ---- 2. Overlay: TurboID vs Flag IP ----
cat("\n[2/3] Generating TurboID vs Flag IP overlay...\n")

if ("turboid" %in% names(experiments) && "flag_ip" %in% names(experiments)) {
  p_overlay <- plot_volcano_overlay(
    experiments$turboid, "TurboID",
    experiments$flag_ip, "Flag_IP",
    label_genes = known_interactors,
    title = "TurboID vs Flag Co-IP Comparison"
  )
  save_figure(p_overlay, "volcano_overlay_turboid_flag", width = 10, height = 8)
} else {
  cat("  (skipped: turboid or flag_ip data not found)\n")
}

# ---- 3. Overlay: WT vs POI vs POI+Hormone ----
cat("\n[3/3] Generating condition comparison overlays...\n")

if ("wt_vs_poi" %in% names(experiments) && "poi_vs_poi_hormone" %in% names(experiments)) {
  p_cond <- plot_volcano_overlay(
    experiments$wt_vs_poi, "WT_vs_POI",
    experiments$poi_vs_poi_hormone, "POI_vs_POIHormone",
    label_genes = known_interactors,
    title = "Wild Type vs POI vs POI + Hormone"
  )
  save_figure(p_cond, "volcano_overlay_conditions", width = 10, height = 8)
} else {
  cat("  (skipped: condition data not found)\n")
}

cat("\n=========================================\n")
cat(" Volcano plots complete!\n")
cat(sprintf(" Output: %s/\n", FIGURE_DIR))
cat("=========================================\n")
