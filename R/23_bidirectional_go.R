###############################################################################
# 23_bidirectional_go.R
# Bidirectional GO enrichment dot plot.
#
# WHAT THIS DOES (requested by Dr. Aruna, July 12 voice message):
#   Take an experiment, split significant proteins into up-regulated and
#   down-regulated sets. Run GO enrichment on each separately. Combine into
#   a SINGLE dot plot where:
#     - UP-regulated GO terms on the RIGHT (positive GeneRatio)
#     - DOWN-regulated GO terms on the LEFT (negative GeneRatio)
#     - Color = p-adjusted value
#     - Size = gene count
#
#   This gives a one-figure overview of what pathways go up vs down.
#
# Usage:
#   make bidirectional-go
###############################################################################

source('R/01_config.R', chdir = TRUE)
source('R/utils.R')

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

cat("\n=========================================\n")
cat(" Bidirectional GO Dot Plot\n")
cat("=========================================\n\n")

# ---- Configuration ----
# Default experiment: TRIP4 vs WT (main experiment)
BIDIR_EXP <- "BK467_TRIP4_vs_BK467_WT"
ONT_TO_RUN <- c("BP")  # Biological Process — most informative
TOP_N_PER_DIR <- 15    # Show top 15 terms per direction

# ---- Load experiment data ----
experiments <- load_all_experiments()
df <- find_experiment(experiments, BIDIR_EXP)

if (is.null(df)) {
  cat("ERROR:", BIDIR_EXP, "not found.\n")
  quit(status = 1)
}

cat(sprintf("  Experiment: %s (%d proteins)\n", BIDIR_EXP, nrow(df)))

# ---- Split into up and down gene sets ----
up_genes <- df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF &
                    !is.na(df$gene)]
dn_genes <- df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC < -LOG2FC_CUTOFF &
                    !is.na(df$gene)]

cat(sprintf("  Up-regulated:   %d genes (log2FC > %.1f)\n", length(up_genes), LOG2FC_CUTOFF))
cat(sprintf("  Down-regulated: %d genes (log2FC < -%.1f)\n", length(dn_genes), LOG2FC_CUTOFF))

if (length(up_genes) < 10) cat("  WARNING: Few up-genes for GO enrichment.\n")
if (length(dn_genes) < 10) cat("  WARNING: Few down-genes for GO enrichment.\n")

# ---- Helper: run enrichGO and return processed results ----
run_bidir_go <- function(gene_set, direction, ont) {
  if (length(gene_set) < 5) {
    cat(sprintf("  Skipping %s %s: too few genes (%d)\n", direction, ont, length(gene_set)))
    return(NULL)
  }

  cat(sprintf("  Running enrichGO: %s, %s (%d genes)...\n", direction, ont, length(gene_set)))

  ego <- enrichGO(
    gene         = gene_set,
    OrgDb        = org.Hs.eg.db,
    keyType      = KEYTYPE,
    ont          = ont,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat(sprintf("    No enriched terms found for %s %s\n", direction, ont))
    return(NULL)
  }

  # Simplify to remove redundant terms (if rrvgo available)
  ego_simplified <- tryCatch({
    simplify(ego, cutoff = 0.7)
  }, error = function(e) ego)

  res <- as.data.frame(ego_simplified)
  if (nrow(res) == 0) return(NULL)

  # Take top N by p.adjust
  res <- res[order(res$p.adjust), ]
  res <- head(res, TOP_N_PER_DIR)

  # Parse GeneRatio "a/b" → numeric
  res$GeneRatioNum <- sapply(res$GeneRatio, function(x) {
    parts <- strsplit(as.character(x), "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA
  })

  # Add direction sign: up = positive, down = negative
  res$signed_GeneRatio <- if (direction == "up") res$GeneRatioNum else -res$GeneRatioNum
  res$direction <- direction
  res$neg_log10_padj <- -log10(res$p.adjust)

  # Truncate long descriptions for readability
  res$short_Description <- sapply(res$Description, function(x) {
    if (nchar(x) > 55) paste0(substr(x, 1, 52), "...") else x
  })

  return(res)
}

# ---- Run GO for each ontology ----
for (ont in ONT_TO_RUN) {
  cat(sprintf("\n--- Ontology: %s ---\n", ont))

  up_res <- run_bidir_go(up_genes, "up", ont)
  dn_res <- run_bidir_go(dn_genes, "down", ont)

  combined <- rbind(up_res, dn_res)

  if (is.null(combined) || nrow(combined) == 0) {
    cat("  No enriched terms in either direction. Skipping plot.\n")
    next
  }

  # Sort: up at top, down at bottom (reversed for ggplot y-axis)
  combined <- combined[order(combined$signed_GeneRatio, decreasing = TRUE), ]

  # Factor Description to preserve order
  combined$short_Description <- factor(combined$short_Description,
                                        levels = rev(combined$short_Description))

  cat(sprintf("  Combined: %d terms (up=%d, down=%d)\n", nrow(combined),
              sum(combined$direction == "up"),
              sum(combined$direction == "down")))

  # ---- Build the bidirectional dot plot ----
  p <- ggplot2::ggplot(combined,
                        ggplot2::aes(x = signed_GeneRatio,
                                     y = short_Description,
                                     color = neg_log10_padj,
                                     size = Count)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient(
      low = "#F0E442",     # Warm Yellow (low significance)
      high = "#0072B2",    # Deep Navy (high significance)
      name = expression(-Log[10]~(adjusted~italic(p)~value))
    ) +
    ggplot2::scale_size_continuous(
      name = "Gene Count",
      range = c(2, 8)      # Bigger dots for poster
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Signed~Gene~Ratio~(left~"="~down-regulated~"|"~right~"="~up-regulated)),
      y = NULL,
      title = sprintf("GO Enrichment — %s (%s)\nUp-regulated (right) vs Down-regulated (left)",
                      gsub("_", " ", BIDIR_EXP), ont),
      subtitle = sprintf("Up: %d sig genes → %d terms | Down: %d sig genes → %d terms",
                         length(up_genes), sum(combined$direction == "up", na.rm = TRUE),
                         length(dn_genes), sum(combined$direction == "down", na.rm = TRUE))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, color = "grey30"),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(size = 8),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")
    )

  # Add direction annotations
  p <- p +
    ggplot2::annotate("text", x = max(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "UP", color = "#D55E00",
                      fontface = "bold", size = 4, vjust = 0) +
    ggplot2::annotate("text", x = min(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "DOWN", color = "#0072B2",
                      fontface = "bold", size = 4, vjust = 0)

  filename <- sprintf("bidirectional_go_%s_%s", gsub("_", "-", BIDIR_EXP), ont)
  save_figure(p, filename, width = 10, height = max(6, nrow(combined) * 0.35))

  # Also export the combined results table
  export_cols <- c("Description", "direction", "GeneRatio", "Count",
                   "pvalue", "p.adjust", "signed_GeneRatio")
  export_cols <- intersect(export_cols, names(combined))
  save_table(combined[, export_cols], filename)

  cat(sprintf("  Saved: %s\n", filename))
}

cat("\nDone. Output saved to output/figures/.\n")
