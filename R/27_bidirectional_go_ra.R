###############################################################################
# 27_bidirectional_go_ra.R
# Bidirectional GO plot comparing -RA vs +RA conditions.
#
# Instead of up/down within one experiment, this compares:
#   LEFT  = GO enrichment of TRIP4 TurboID WITHOUT RA (BK467_TRIP4_vs_BK467_WT)
#   RIGHT = GO enrichment of TRIP4 TurboID WITH RA (BK467_TRIP4_RA02_vs_BK467_WT)
#
# Universe = all proteins detected in the respective experiment.
# Foreground = log2FC >= 0.5 AND padj < 0.05 (enriched only, positive direction).
#
# Usage:
#   make bidirectional-go-ra
###############################################################################

source('R/01_config.R', chdir = TRUE)
source('R/utils.R')

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)

cat("\n=========================================\n")
cat(" Bidirectional GO: -RA vs +RA\n")
cat("=========================================\n\n")

# ---- Configuration ----
EXP_NO_RA <- "BK467_TRIP4_vs_BK467_WT"        # Without RA
EXP_WITH_RA <- "BK467_TRIP4_RA02_vs_BK467_WT"   # With RA
ONT_TO_RUN <- c("BP", "CC", "MF")
TOP_N <- 15

# ---- Load data ----
experiments <- load_all_experiments()

df_nora <- find_experiment(experiments, EXP_NO_RA)
df_ra   <- find_experiment(experiments, EXP_WITH_RA)

if (is.null(df_nora)) {
  cat("ERROR:", EXP_NO_RA, "not found.\n")
  quit(status = 1)
}
if (is.null(df_ra)) {
  cat("ERROR:", EXP_WITH_RA, "not found.\n")
  quit(status = 1)
}

cat(sprintf("  Without RA (%s): %d proteins\n", EXP_NO_RA, nrow(df_nora)))
cat(sprintf("  With RA (%s): %d proteins\n", EXP_WITH_RA, nrow(df_ra)))

# ---- Extract enriched genes (positive log2FC only) ----
# Using >= for floating point tolerance
nora_genes <- df_nora$gene[df_nora$padj < P_VALUE_CUTOFF &
                           df_nora$log2FC >= LOG2FC_CUTOFF &
                           !is.na(df_nora$gene)]
nora_genes <- unique(nora_genes)

ra_genes <- df_ra$gene[df_ra$padj < P_VALUE_CUTOFF &
                       df_ra$log2FC >= LOG2FC_CUTOFF &
                       !is.na(df_ra$gene)]
ra_genes <- unique(ra_genes)

# Universe for each
nora_universe <- unique(df_nora$gene[!is.na(df_nora$gene)])
ra_universe   <- unique(df_ra$gene[!is.na(df_ra$gene)])

cat(sprintf("  -RA enriched: %d genes (universe: %d)\n", length(nora_genes), length(nora_universe)))
cat(sprintf("  +RA enriched: %d genes (universe: %d)\n", length(ra_genes), length(ra_universe)))

# ---- Helper: run GO and return processed data frame ----
run_go_for_condition <- function(genes, universe, condition_label, ont) {
  if (length(genes) < 5) {
    cat(sprintf("  Skipping %s %s: too few genes (%d)\n", condition_label, ont, length(genes)))
    return(NULL)
  }

  cat(sprintf("  Running enrichGO: %s, %s (%d genes)...\n", condition_label, ont, length(genes)))

  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = KEYTYPE,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 2,
    maxGSSize     = 5000
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat(sprintf("    No enriched terms for %s %s\n", condition_label, ont))
    return(NULL)
  }

  # Simplify
  ego_simplified <- tryCatch({
    simplify(ego, cutoff = 0.7)
  }, error = function(e) ego)

  res <- as.data.frame(ego_simplified)
  if (nrow(res) == 0) return(NULL)

  res <- res[order(res$p.adjust), ]
  res <- head(res, TOP_N)

  # Parse GeneRatio
  res$GeneRatioNum <- sapply(res$GeneRatio, function(x) {
    parts <- strsplit(as.character(x), "/")[[1]]
    if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA
  })

  # Direction: -RA = negative (left), +RA = positive (right)
  res$signed_GeneRatio <- if (condition_label == "-RA")
    -res$GeneRatioNum else res$GeneRatioNum
  res$condition <- condition_label
  res$neg_log10_padj <- -log10(res$p.adjust)

  res$short_Description <- sapply(res$Description, function(x) {
    if (nchar(x) > 55) paste0(substr(x, 1, 52), "...") else x
  })

  return(res)
}

ONT_NAMES <- c("BP" = "Biological Process",
               "CC" = "Cellular Component",
               "MF" = "Molecular Function")

# ---- Run for each ontology ----
for (ont in ONT_TO_RUN) {
  cat(sprintf("\n--- Ontology: %s ---\n", ont))

  nora_res <- run_go_for_condition(nora_genes, nora_universe, "-RA", ont)
  ra_res   <- run_go_for_condition(ra_genes, ra_universe, "+RA", ont)

  combined <- rbind(nora_res, ra_res)

  if (is.null(combined) || nrow(combined) == 0) {
    cat("  No enriched terms in either condition. Skipping.\n")
    next
  }

  # Sort: +RA at top, -RA at bottom
  combined <- combined[order(combined$signed_GeneRatio, decreasing = TRUE), ]
  combined$short_Description <- make.unique(combined$short_Description)
  combined$short_Description <- factor(combined$short_Description,
                                        levels = rev(combined$short_Description))

  cat(sprintf("  Combined: %d terms (-RA=%d, +RA=%d)\n", nrow(combined),
              sum(combined$condition == "-RA"),
              sum(combined$condition == "+RA")))

  # ---- Build bidirectional dot plot ----
  p <- ggplot2::ggplot(combined,
                        ggplot2::aes(x = signed_GeneRatio,
                                     y = short_Description,
                                     color = neg_log10_padj,
                                     size = Count)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradientn(
      colors = c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"),
      name = expression(-Log[10]~(adjusted~italic(p)~value))
    ) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid",
                        color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Signed~Gene~Ratio~(left~"="-~RA~"|"~right~"="+~RA)),
      y = NULL,
      title = sprintf("GO Enrichment: -RA vs +RA (%s)", ONT_NAMES[ont]),
      subtitle = sprintf("-RA: %d sig genes → %d terms | +RA: %d sig genes → %d terms",
                         length(nora_genes), sum(combined$condition == "-RA", na.rm = TRUE),
                         length(ra_genes), sum(combined$condition == "+RA", na.rm = TRUE))
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

  # Direction annotations
  p <- p +
    ggplot2::annotate("text",
                      x = max(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "+ RA", color = "#D55E00",
                      fontface = "bold", size = 4, vjust = 0) +
    ggplot2::annotate("text",
                      x = min(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "- RA", color = "#0072B2",
                      fontface = "bold", size = 4, vjust = 0)

  filename <- sprintf("bidirectional_go_RA_%s", ont)
  save_figure(p, filename, width = 10,
              height = max(6, nrow(combined) * 0.35))

  # Export table
  export_cols <- c("Description", "condition", "GeneRatio", "Count",
                   "pvalue", "p.adjust", "signed_GeneRatio")
  export_cols <- intersect(export_cols, names(combined))
  save_table(combined[, export_cols], filename)

  cat(sprintf("  Saved: %s\n", filename))
}

cat("\n=========================================\n")
cat(" Bidirectional GO (-RA vs +RA) complete!\n")
cat("=========================================\n")
