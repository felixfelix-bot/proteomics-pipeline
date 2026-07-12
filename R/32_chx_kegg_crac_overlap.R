###############################################################################
# 32_chx_kegg_crac_overlap.R
# KEGG pathway enrichment on CHX-enriched/depleted sets + CRAC↔TurboID overlap.
#
#   1. KEGG enrichment on CHX-enriched and CHX-depleted protein sets
#   2. Export enriched/depleted protein lists to CSV (gene, log2FC, padj)
#   3. Find common genes between CRAC significant and TurboID TRIP4 vs WT
#   4. Volcano plot of TurboID TRIP4 vs WT, labeling CRAC↔TurboID common genes
#
# Usage:
#   make chx-kegg-crac
###############################################################################

cat("\n=========================================\n")
cat(" CHX KEGG + CRAC/TurboID Overlap\n")
cat("=========================================\n\n")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggrepel)

# ---- Load experiments ----
experiments <- load_all_experiments()

chx_exp <- "TRIP4_CHX_vs_TRIP4_DMSO"
turbo_exp <- "BK467_TRIP4_vs_BK467_WT"

df_chx   <- find_experiment(experiments, chx_exp)
df_turbo <- find_experiment(experiments, turbo_exp)

if (is.null(df_chx)) {
  cat("ERROR: CHX experiment '", chx_exp, "' not found.\n", sep = "")
  quit(status = 1)
}
if (is.null(df_turbo)) {
  cat("ERROR: TurboID experiment '", turbo_exp, "' not found.\n", sep = "")
  quit(status = 1)
}

cat(sprintf("  CHX vs DMSO:        %d proteins\n", nrow(df_chx)))
cat(sprintf("  TurboID TRIP4 vs WT: %d proteins\n", nrow(df_turbo)))

# =====================================================================
# Section 1: KEGG pathway enrichment on CHX-enriched and CHX-depleted
# =====================================================================
cat("\n--- Section 1: KEGG on CHX data ---\n\n")

# CHX-enriched: padj < 0.05 & log2FC >= 0.5 (positive = more abundant with CHX)
# CHX-depleted: padj < 0.05 & log2FC <= -0.5 (negative = less abundant with CHX)
chx_enriched_genes <- unique(df_chx$gene[df_chx$padj < P_VALUE_CUTOFF &
                                           df_chx$log2FC >= LOG2FC_CUTOFF &
                                           !is.na(df_chx$gene)])
chx_depleted_genes <- unique(df_chx$gene[df_chx$padj < P_VALUE_CUTOFF &
                                           df_chx$log2FC <= -LOG2FC_CUTOFF &
                                           !is.na(df_chx$gene)])

cat(sprintf("  CHX-enriched: %d genes\n", length(chx_enriched_genes)))
cat(sprintf("  CHX-depleted: %d genes\n", length(chx_depleted_genes)))

# Universe = all proteins in the CHX experiment
universe_genes <- unique(df_chx$gene[!is.na(df_chx$gene)])
cat(sprintf("  Universe: %d proteins\n", length(universe_genes)))

run_kegg_chx <- function(genes, set_label) {
  if (length(genes) < 5) {
    cat(sprintf("\n  [%s] Too few genes (%d). Skipping KEGG.\n", set_label, length(genes)))
    return(NULL)
  }

  cat(sprintf("\n  [%s] %d genes → KEGG enrichment\n", set_label, length(genes)))

  ekegg <- tryCatch({
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    universe_entrez <- bitr(universe_genes, fromType = "SYMBOL",
                            toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    ekegg <- enrichKEGG(
      gene = unique(entrez_map$ENTREZID),
      universe = unique(universe_entrez$ENTREZID),
      organism = "hsa", pAdjustMethod = "BH",
      pvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 5000
    )
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    ekegg
  }, error = function(e) {
    cat(sprintf("    ERROR (network?): %s\n", conditionMessage(e)))
    NULL
  })

  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    cat("    No enriched KEGG pathways\n")
    return(NULL)
  }

  res_df <- as.data.frame(ekegg)
  cat(sprintf("    %d enriched pathways\n", nrow(res_df)))

  save_table(res_df, sprintf("KEGG_CHX_%s", set_label))

  n_show <- min(15, nrow(res_df))
  p <- dotplot(ekegg, showCategory = n_show) +
    ggplot2::labs(
      title = sprintf("CHX %s — KEGG Pathways", gsub("_", " ", set_label)),
      x = "GeneRatio"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  save_figure(p, sprintf("KEGG_CHX_%s_dotplot", set_label),
              width = 10, height = max(6, n_show * 0.4))
}

run_kegg_chx(chx_enriched_genes, "enriched")
run_kegg_chx(chx_depleted_genes, "depleted")

# =====================================================================
# Section 1b: GO enrichment (BP, MF, CC) on CHX-enriched and CHX-depleted
# =====================================================================
cat("\n--- Section 1b: GO on CHX data ---\n\n")

ONT_LABELS <- c("BP" = "Biological Process",
                "MF" = "Molecular Function",
                "CC" = "Cellular Component")

run_go_chx <- function(genes, set_label) {
  if (length(genes) < 5) {
    cat(sprintf("\n  [%s] Too few genes (%d). Skipping GO.\n", set_label, length(genes)))
    return(NULL)
  }

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  [%s] [%s] enrichGO... ", set_label, ont))
    ego <- tryCatch({
      enrichGO(gene = genes, universe = universe_genes,
               OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = GO_PVALUE_CUTOFF, qvalueCutoff = GO_QVALUE_CUTOFF,
               minGSSize = 2, maxGSSize = 5000)
    }, error = function(e) NULL)

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      cat("No terms\n"); next
    }

    ego_s <- tryCatch(simplify(ego, cutoff = 0.7), error = function(e) ego)
    res <- as.data.frame(ego_s)
    cat(sprintf("%d terms\n", nrow(res)))

    prefix <- sanitize_filename(paste0("GO_CHX_", set_label, "_", ont))
    save_table(res, prefix)

    n_show <- min(20, nrow(res))
    go_title <- sprintf("CHX %s — %s", gsub("_", " ", set_label), ONT_LABELS[ont])
    p_dot <- dotplot(ego_s, showCategory = n_show) +
      ggplot2::labs(title = go_title, x = "GeneRatio") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = max(6, n_show * 0.4))

    n_bar <- min(15, nrow(res))
    p_bar <- barplot(ego_s, showCategory = n_bar, orderBy = "Count") +
      ggplot2::labs(title = go_title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = max(6, n_bar * 0.5))
  }
}

run_go_chx(chx_enriched_genes, "enriched")
run_go_chx(chx_depleted_genes, "depleted")

# =====================================================================
# Section 1c: CHX vs DMSO volcano (top 30 per side labeled)
# =====================================================================
cat("\n--- Section 1c: CHX vs DMSO volcano (top 30 per side) ---\n\n")

chx_vol <- df_chx
chx_vol$neglog10p <- -log10(chx_vol$padj)
chx_vol$category <- "nonsig"
chx_vol$category[chx_vol$padj < P_VALUE_CUTOFF & chx_vol$log2FC >= LOG2FC_CUTOFF]  <- "enriched"
chx_vol$category[chx_vol$padj < P_VALUE_CUTOFF & chx_vol$log2FC <= -LOG2FC_CUTOFF] <- "depleted"

# Top 30 by combined significance score on each side
sig_only <- chx_vol[chx_vol$category %in% c("enriched", "depleted") &
                     !is.na(chx_vol$log2FC) & !is.na(chx_vol$neglog10p), ]
sig_only$combined_score <- abs(sig_only$log2FC) * sig_only$neglog10p

top_up <- head(sig_only[sig_only$log2FC > 0, ][order(-sig_only[sig_only$log2FC > 0, ]$combined_score), ], 30)
top_dn <- head(sig_only[sig_only$log2FC < 0, ][order(-sig_only[sig_only$log2FC < 0, ]$combined_score), ], 30)
label_data_chx <- rbind(top_up, top_dn)
cat(sprintf("  Labeling %d enriched + %d depleted (top 30 each)\n",
            nrow(top_up), nrow(top_dn)))

CHX_VOL_COLORS <- c(
  "enriched" = GLOBAL_COLORS[["chx_enriched"]],
  "depleted" = GLOBAL_COLORS[["dmso_enriched"]],
  "nonsig"   = GLOBAL_COLORS[["nonsig"]]
)
CHX_VOL_LABELS <- c(
  "enriched" = "Enriched in CHX",
  "depleted" = "Enriched in DMSO",
  "nonsig"   = "Not significant"
)

p_chx_volcano <- ggplot2::ggplot(chx_vol,
  ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.8) +
  ggplot2::scale_color_manual(
    values = CHX_VOL_COLORS, labels = CHX_VOL_LABELS,
    name = NULL, drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_data_chx, ggplot2::aes(label = gene),
    size = 2.5, fontface = "bold", max.overlaps = 60,
    show.legend = FALSE, bg.color = "white", bg.r = 0.15,
    segment.color = "grey40", segment.size = 0.3,
    min.segment.length = 0.2
  ) +
  ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID: CHX vs DMSO (top 30 per side)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = "right",
    panel.grid.minor = ggplot2::element_blank()
  )

save_figure(p_chx_volcano, "volcano_chx_vs_dmso_top30", width = 8, height = 6)
cat("  Saved volcano_chx_vs_dmso_top30\n")

# =====================================================================
# Section 2: Export enriched/depleted protein lists to CSV
# =====================================================================
cat("\n--- Section 2: Export enriched/depleted CSVs ---\n\n")

chx_enriched_df <- df_chx[df_chx$padj < P_VALUE_CUTOFF &
                            df_chx$log2FC >= LOG2FC_CUTOFF &
                            !is.na(df_chx$gene) &
                            !is.na(df_chx$log2FC) &
                            !is.na(df_chx$padj),
                          c("gene", "log2FC", "padj")]
chx_enriched_df <- chx_enriched_df[order(-chx_enriched_df$log2FC), ]

chx_depleted_df <- df_chx[df_chx$padj < P_VALUE_CUTOFF &
                            df_chx$log2FC <= -LOG2FC_CUTOFF &
                            !is.na(df_chx$gene) &
                            !is.na(df_chx$log2FC) &
                            !is.na(df_chx$padj),
                          c("gene", "log2FC", "padj")]
chx_depleted_df <- chx_depleted_df[order(chx_depleted_df$log2FC), ]

save_table(chx_enriched_df, "CHX_enriched_proteins")
save_table(chx_depleted_df, "CHX_depleted_proteins")

cat(sprintf("  Exported CHX-enriched: %d proteins\n", nrow(chx_enriched_df)))
cat(sprintf("  Exported CHX-depleted: %d proteins\n", nrow(chx_depleted_df)))

# =====================================================================
# Section 3: CRAC ↔ TurboID overlap
# =====================================================================
cat("\n--- Section 3: CRAC/TurboID overlap ---\n\n")

crac_path <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")

if (!file.exists(crac_path)) {
  cat("WARNING: CRAC data file not found:", crac_path, "\n")
  cat("Skipping CRAC overlap and labeled volcano.\n")
} else {
  cat("Loading CRAC data...\n")
  crac_raw <- readr::read_csv(crac_path, show_col_types = FALSE)
  cat(sprintf("  Raw rows: %d\n", nrow(crac_raw)))

  crac_df <- data.frame(
    gene   = crac_raw[[CRAC_GENE_COL]],
    log2FC = crac_raw[[CRAC_LOG2FC_COL]],
    padj   = crac_raw[[CRAC_PADJ_COL]],
    stringsAsFactors = FALSE
  )

  crac_df <- crac_df[!is.na(crac_df$gene) & crac_df$gene != "" &
                       !is.na(crac_df$log2FC) & !is.na(crac_df$padj), ]
  cat(sprintf("  Clean rows: %d proteins\n", nrow(crac_df)))

  crac_sig <- crac_df$gene[crac_df$padj < P_VALUE_CUTOFF &
                             abs(crac_df$log2FC) > LOG2FC_CUTOFF]
  crac_sig <- unique(crac_sig[!is.na(crac_sig)])
  cat(sprintf("  CRAC significant: %d genes\n", length(crac_sig)))

  turbo_sig <- get_significant_genes(df_turbo, direction = "enriched")
  cat(sprintf("  TurboID TRIP4 enriched: %d genes\n", length(turbo_sig)))

  common_genes <- intersect(crac_sig, turbo_sig)
  cat(sprintf("  Common (CRAC ∩ TurboID): %d genes\n", length(common_genes)))

  if (length(common_genes) > 0) {
    cat(sprintf("    %s\n", paste(common_genes, collapse = ", ")))
  }

  # Export overlap CSV with gene, log2FC, padj from TurboID
  overlap_df <- df_turbo[df_turbo$gene %in% common_genes,
                         c("gene", "log2FC", "padj")]
  overlap_df$category <- "CRAC_TurboID_common"
  overlap_df <- overlap_df[order(-abs(overlap_df$log2FC)), ]
  save_table(overlap_df, "CRAC_TurboID_common_genes")

  cat(sprintf("  Exported overlap table: %d rows\n", nrow(overlap_df)))

  # =====================================================================
  # Section 4: Volcano plot of TurboID vs WT, labeling common genes
  # =====================================================================
  cat("\n--- Section 4: TurboID vs WT volcano (overlap labels) ---\n\n")

  toPlot <- df_turbo
  toPlot$neglog10p <- -log10(toPlot$padj)

  toPlot$category <- "nonsig"
  toPlot$category[toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC > LOG2FC_CUTOFF]  <- "enriched_up"
  toPlot$category[toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC < -LOG2FC_CUTOFF] <- "enriched_dn"
  toPlot$category <- factor(toPlot$category,
    levels = c("enriched_up", "enriched_dn", "nonsig"))

  # Label data = only the CRAC-TurboID common genes
  label_data <- toPlot[toPlot$gene %in% common_genes &
                         !is.na(toPlot$gene), ]

  PLOT_COLORS <- c(
    "enriched_up" = GLOBAL_COLORS[["enriched_up"]],
    "enriched_dn" = GLOBAL_COLORS[["enriched_dn"]],
    "nonsig"      = GLOBAL_COLORS[["nonsig"]]
  )

  PLOT_LABELS <- c(
    "enriched_up" = "Enriched in TRIP4",
    "enriched_dn" = "Enriched in WT",
    "nonsig"      = "Not significant"
  )

  p_volcano <- ggplot2::ggplot(toPlot,
    ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = PLOT_COLORS,
      labels = PLOT_LABELS,
      name = NULL, drop = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_data, ggplot2::aes(label = gene),
      size = 2.8, fontface = "bold", max.overlaps = 50,
      show.legend = FALSE, bg.color = "white", bg.r = 0.15,
      segment.color = "grey40", segment.size = 0.3,
      min.segment.length = 0.2,
      arrow = grid::arrow(length = grid::unit(0.01, "npc"))
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF),
                        linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
                        linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = "TurboID TRIP4 vs WT — CRAC↔TurboID common genes labeled"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      axis.title.x = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(colour = "black", size = 8),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  save_figure(p_volcano, "volcano_turboid_vs_wt_crac_overlap",
              width = 8, height = 6)
  cat("  Saved volcano_turboid_vs_wt_crac_overlap\n")
}

cat("\n=========================================\n")
cat(" CHX KEGG + CRAC/TurboID overlap complete!\n")
cat("=========================================\n")
