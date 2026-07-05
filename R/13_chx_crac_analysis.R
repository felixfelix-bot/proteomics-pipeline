###############################################################################
# 13_chx_crac_analysis.R
# Volcano plots + GO analysis for new datasets:
#
#   1. TRIP4 CHX vs TRIP4 DMSO (translation inhibitor effect)
#   2. TRIP4 CHX vs WT
#   3. TRIP4 DMSO vs WT
#   4. FLAG-TRIP4 CRAC data (RNA interactome — different column structure)
#
# The CHX/DMSO datasets use the same column format as the existing pipeline
# (Gene, logFC, adj.P.Val). They'll be loaded automatically by
# load_all_experiments() since they end in _diffEx_minProb.csv.
#
# The CRAC data has different columns (external_gene_name, logFC, FDR).
# It's loaded separately with a custom loader function.
#
# Usage:
#   make chx-crack-analysis
###############################################################################
cat("\n=========================================\n")
cat(" CHX/DMSO + CRAC Analysis\n")
cat("=========================================\n\n")

# ---- Load mass spec experiments (includes CHX/DMSO automatically) ----
experiments <- load_all_experiments()

# ---- CRAC title map (no underscores, descriptive) ----
CHX_TITLES <- list(
  "TRIP4_CHX_vs_TRIP4_DMSO" = "TRIP4 TurboID: CHX vs DMSO",
  "TRIP4_CHX_vs_WT"         = "TRIP4 TurboID: CHX vs Wild Type",
  "TRIP4_DMSO_vs_WT"        = "TRIP4 TurboID: DMSO vs Wild Type"
)

# ---- CHX/DMSO volcano plots ----
# Same style as the main targeted volcano: label top 20 by significance,
# WT-enriched = gray, TRIP4-enriched = orange
make_chx_volcano <- function(df, title, n_top = 20) {
  toPlot <- df
  toPlot$neglog10p <- -log10(toPlot$padj)

  toPlot$category <- "nonsig"
  toPlot$category[toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC > LOG2FC_CUTOFF] <- "enriched_up"
  toPlot$category[toPlot$padj < P_VALUE_CUTOFF & toPlot$log2FC < -LOG2FC_CUTOFF] <- "enriched_dn"
  toPlot$category <- factor(toPlot$category,
    levels = c("enriched_up", "enriched_dn", "nonsig"))

  # Top N labels by combined significance score
  sig_only <- toPlot[toPlot$category %in% c("enriched_up", "enriched_dn"), ]
  sig_only <- sig_only[!is.na(sig_only$log2FC) & !is.na(sig_only$neglog10p), ]
  sig_only$combined_score <- abs(sig_only$log2FC) * sig_only$neglog10p

  top_up <- head(sig_only[sig_only$log2FC > 0, ][order(-sig_only[sig_only$log2FC > 0, ]$combined_score), ], n_top / 2)
  top_dn <- head(sig_only[sig_only$log2FC < 0, ][order(-sig_only[sig_only$log2FC < 0, ]$combined_score), ], n_top / 2)
  label_data <- rbind(top_up, top_dn)

  CHX_COLORS <- c(
    "enriched_up" = GLOBAL_COLORS[["enriched_up"]],
    "enriched_dn" = GLOBAL_COLORS[["enriched_dn"]],
    "nonsig"      = GLOBAL_COLORS[["nonsig"]]
  )

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(
      values = CHX_COLORS,
      labels = c("enriched_up" = "Enriched in TRIP4", "enriched_dn" = "Enriched in control", "nonsig" = "Not significant"),
      name = NULL
    ) +
    ggrepel::geom_text_repel(
      data = label_data, ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold", max.overlaps = 50,
      show.legend = FALSE, bg.color = "white", bg.r = 0.15
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggplot2::labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~(adjusted~italic(p)~value)), title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
                   axis.title.x = ggplot2::element_text(hjust = 0.5),
                   axis.text = ggplot2::element_text(colour = "black", size = 8),
                   legend.position = "right", panel.grid.minor = ggplot2::element_blank())
  return(p)
}

cat("\n--- CHX/DMSO Volcano Plots ---\n\n")
for (exp_name in names(CHX_TITLES)) {
  if (exp_name %in% names(experiments)) {
    cat(sprintf("[%s] Creating volcano...\n", exp_name))
    p <- make_chx_volcano(experiments[[exp_name]], CHX_TITLES[[exp_name]])
    save_figure(p, paste0("volcano_", exp_name), width = 8, height = 6)
  } else {
    cat(sprintf("  WARNING: %s not found in data\n", exp_name))
  }
}

# ---- CHX/DMSO GO analysis ----
cat("\n--- CHX/DMSO GO Enrichment ---\n\n")
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

universe <- unique(unlist(lapply(experiments, function(df) df$gene)))

ONT_LABELS <- c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")

for (exp_name in names(CHX_TITLES)) {
  if (!exp_name %in% names(experiments)) next
  sig_genes <- get_significant_genes(experiments[[exp_name]])
  cat(sprintf("\n[%s] %d significant genes → GO analysis\n", exp_name, length(sig_genes)))
  if (length(sig_genes) < 5) { cat("  Skipped: < 5 genes\n"); next }

  for (ont in c("BP", "MF", "CC")) {
    result <- tryCatch({
      enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.05, universe = universe)
    }, error = function(e) NULL)

    if (is.null(result) || nrow(as.data.frame(result)) == 0) { cat(sprintf("  [%s] No enriched terms\n", ont)); next }

    result <- tryCatch(simplify(result, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) result)
    res_df <- as.data.frame(result)
    cat(sprintf("  [%s] %d terms after simplify\n", ont, nrow(res_df)))

    prefix <- sanitize_filename(paste0("GO_", exp_name, "_", ont))
    save_table(res_df, prefix)

    n_show <- min(20, nrow(res_df))
    fig_h <- max(7, n_show * 0.4)
    title_str <- paste0("GO analysis of ", CHX_TITLES[[exp_name]], " — ", ONT_LABELS[ont])

    p_dot <- dotplot(result, showCategory = n_show, title = title_str) +
      ggplot2::labs(x = "Gene Ratio", color = "p-adjusted value") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(size = 7))
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_h)

    n_bar <- min(15, nrow(res_df))
    fig_hb <- max(7, n_bar * 0.5)
    p_bar <- barplot(result, showCategory = n_bar, orderBy = "Count", title = title_str) +
      ggplot2::labs(color = "p-adjusted value", fill = "p-adjusted value") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(size = 7))
    save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = fig_hb)
  }
}

# =====================================================================
# CRAC DATA ANALYSIS (different column structure)
# =====================================================================
cat("\n\n=========================================\n")
cat(" CRAC Data Analysis\n")
cat("=========================================\n\n")

crac_path <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")

if (!file.exists(crac_path)) {
  cat("WARNING: CRAC data file not found:", crac_path, "\n")
  cat("Expected at: data/FLAG-TRIP4_list_CRACdata.csv\n")
} else {
  cat("Loading CRAC data...\n")
  crac_df <- readr::read_csv(crac_path, show_col_types = FALSE)
  cat(sprintf("  Raw rows: %d\n", nrow(crac_df)))

  # Hardcoded CRAC column names (confirmed from actual file headers)
  # CRAC data uses different names than the mass spec pipeline
  crac_gene_col <- "external_gene_name"
  crac_fc_col   <- "logFC"
  crac_padj_col <- "FDR"

  # Verify columns exist
  missing <- c()
  if (!(crac_gene_col %in% colnames(crac_df))) missing <- c(missing, crac_gene_col)
  if (!(crac_fc_col %in% colnames(crac_df)))   missing <- c(missing, crac_fc_col)
  if (!(crac_padj_col %in% colnames(crac_df))) missing <- c(missing, crac_padj_col)

  if (length(missing) > 0) {
    cat("\n  ERROR: Missing columns in CRAC data:", paste(missing, collapse=", "), "\n")
    cat("  Available columns:", paste(colnames(crac_df), collapse=", "), "\n")
  } else {
    # Build clean data frame with standardized names
    crac_clean <- data.frame(
      gene = crac_df[[crac_gene_col]],
      log2FC = crac_df[[crac_fc_col]],
      padj = crac_df[[crac_padj_col]],
      stringsAsFactors = FALSE
    )

    # Remove rows with missing gene names or NA values
    crac_clean <- crac_clean[!is.na(crac_clean$gene) & crac_clean$gene != "" &
                               !is.na(crac_clean$log2FC) & !is.na(crac_clean$padj), ]
    cat(sprintf("  Clean rows: %d proteins\n", nrow(crac_clean)))

    n_sig <- sum(crac_clean$padj < P_VALUE_CUTOFF & abs(crac_clean$log2FC) > LOG2FC_CUTOFF)
    n_up  <- sum(crac_clean$padj < P_VALUE_CUTOFF & crac_clean$log2FC > LOG2FC_CUTOFF)
    cat(sprintf("  Significant: %d (%d up, %d down)\n", n_sig, n_up, n_sig - n_up))

    # ---- CRAC volcano plot ----
    cat("\n  Creating CRAC volcano plot...\n")
    p_crac <- make_chx_volcano(crac_clean, "FLAG-TRIP4 CRAC: RNA Interactome")
    save_figure(p_crac, "volcano_CRAC_FLAG_TRIP4", width = 8, height = 6)

    # ---- CRAC GO analysis ----
    cat("\n  CRAC GO enrichment...\n")
    crac_sig <- crac_clean$gene[crac_clean$padj < P_VALUE_CUTOFF & abs(crac_clean$log2FC) > LOG2FC_CUTOFF]
    crac_sig <- unique(crac_sig[!is.na(crac_sig)])
    cat(sprintf("  %d significant genes for GO\n", length(crac_sig)))

    if (length(crac_sig) >= 5) {
      for (ont in c("BP", "MF", "CC")) {
        result <- tryCatch({
          enrichGO(gene = crac_sig, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05, universe = universe)
        }, error = function(e) NULL)

        if (is.null(result) || nrow(as.data.frame(result)) == 0) { cat(sprintf("    [%s] No enriched terms\n", ont)); next }

        result <- tryCatch(simplify(result, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) result)
        res_df <- as.data.frame(result)
        cat(sprintf("    [%s] %d terms after simplify\n", ont, nrow(res_df)))

        prefix <- sanitize_filename(paste0("GO_CRAC_FLAG_TRIP4_", ont))
        save_table(res_df, prefix)

        n_show <- min(20, nrow(res_df))
        fig_h <- max(7, n_show * 0.4)
        title_str <- paste0("GO analysis of FLAG-TRIP4 CRAC RNA Interactome — ", ONT_LABELS[ont])

        p_dot <- dotplot(result, showCategory = n_show, title = title_str) +
          ggplot2::labs(x = "Gene Ratio", color = "p-adjusted value") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(size = 7))
        save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_h)

        n_bar <- min(15, nrow(res_df))
        fig_hb <- max(7, n_bar * 0.5)
        p_bar <- barplot(result, showCategory = n_bar, orderBy = "Count", title = title_str) +
          ggplot2::labs(color = "p-adjusted value", fill = "p-adjusted value") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), axis.text.y = ggplot2::element_text(size = 7))
        save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = fig_hb)
      }
    }

    # Save CRAC significant gene list
    save_table(data.frame(gene = crac_sig, category = "CRAC_significant"),
               "CRAC_significant_genes")
  }
}

cat("\n=========================================\n")
cat(" CHX/DMSO + CRAC analysis complete!\n")
cat("=========================================\n")
