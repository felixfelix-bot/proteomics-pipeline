###############################################################################
# 19_gsea.R
# GSEA (Gene Set Enrichment Analysis) using the full ranked gene list.
#
# WHAT THIS DOES (vs ORA):
#   ORA (Over-Representation Analysis) uses ONLY significant genes and asks
#   "are they enriched in GO terms?" — requires an arbitrary p-value cutoff.
#
#   GSEA uses ALL genes ranked by a continuous metric (signed statistic) and
#   asks "are genes in this GO term clustered at the top or bottom?" — no
#   cutoff needed, more sensitive for subtle but coordinated changes.
#
# GENE RANKING:
#   We rank by: sign(logFC) * -log10(adj.P.Val)
#   This gives a signed score where:
#     + = up-regulated, high significance (top of ranked list)
#     - = down-regulated, high significance (bottom of ranked list)
#     0 = not significant (middle of ranked list)
#
# EXPERIMENTS ANALYZED:
#   1) TRIP4 TurboID vs WT — core TRIP4 interactome
#   2) RA vs base (+RA vs -RA, within TRIP4) — RA-dependent changes
#   3) CHX vs DMSO — RNA-dependent changes
#   4) CRAC — RNA interactome
#
# OUTPUT:
#   - Dotplot (top 20 terms per ontology)
#   - Ridgeplot (running enrichment scores)
#   - GSEA table CSV with NES, p-value, leading edge genes
#
# Usage:
#   make gsea
###############################################################################
cat("\n=========================================\n")
cat(" GSEA: Gene Set Enrichment Analysis\n")
cat(" (ranked gene list — Lydia style)\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

experiments <- load_all_experiments()

ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

# ---- Helper: run GSEA for one experiment ----
# GSEA takes a named, ranked numeric vector:
#   names = gene symbols
#   values = ranking metric (positive = up, negative = down)
run_gsea <- function(df, experiment_name) {
  cat(sprintf("\n--- %s (%d proteins) ---\n", experiment_name, nrow(df)))

  # ---- Build the ranked gene list ----
  # Rank by: sign(logFC) * -log10(adj.P.Val)
  # This gives a continuous score that captures BOTH direction AND significance.
  # Genes with large +logFC and small padj get high positive ranks.
  # Genes with large -logFC and small padj get low negative ranks.
  # Non-significant genes cluster near zero.
  df$rank_metric <- sign(df$logFC) * -log10(df$adj.P.Val)

  # Remove any rows with NA or Inf in the ranking metric
  df <- df[is.finite(df$rank_metric), ]

  # Sort by rank_metric descending (most up-regulated first)
  df <- df[order(df$rank_metric, decreasing = TRUE), ]

  gene_list <- df$rank_metric
  names(gene_list) <- df$gene

  cat(sprintf("  Ranked %d genes (min=%.2f, max=%.2f)\n",
              length(gene_list), min(gene_list), max(gene_list)))

  # ---- Run GSEA for each ontology (BP, MF, CC) ----
  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  [%s] %s — Running gseGO...\n", ont, ONT_LABELS[ont]))

    result <- tryCatch({
      gseGO(
        geneList      = gene_list,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        minGSSize     = 10,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        verbose       = FALSE
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      cat("    No enriched gene sets found\n")
      next
    }

    res_df <- as.data.frame(result)
    cat(sprintf("    Found %d enriched gene sets\n", nrow(res_df)))

    # Build filename prefix
    safe_exp <- sanitize_filename(gsub("[ _]", "_", experiment_name))
    prefix <- sanitize_filename(paste0("GSEA_", safe_exp, "_", ont))

    # Save full GSEA table
    save_table(res_df, prefix)

    # ---- Dotplot (top 20 terms by NES) ----
    n_show <- min(20, nrow(res_df))
    fig_height <- max(7, n_show * 0.4)

    dot_title <- paste0("GSEA: ", experiment_name, " — ", ONT_LABELS[ont])

    p_dot <- dotplot(result, showCategory = n_show,
                     title = dot_title,
                     label_format = 50) +
      ggplot2::labs(x = "Gene Ratio",
                    color = "p-adjusted value") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 12, height = fig_height)

    # ---- Ridgeplot (running enrichment score distribution) ----
    # Shows how the enrichment score runs across the ranked list
    fig_height_ridge <- max(7, n_show * 0.45)

    p_ridge <- ridgeplot(result, showCategory = n_show,
                         label_format = 50) +
      ggplot2::labs(
        title = dot_title,
        x = "Running Enrichment Score",
        fill = "p-adjusted value"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_ridge, paste0(prefix, "_ridgeplot"),
                width = 12, height = fig_height_ridge)

    # ---- GSEA running score plot (top 5 terms) ----
    # Shows the actual enrichment line going up/down across the ranked list
    if (nrow(res_df) >= 2) {
      n_gsea <- min(5, nrow(res_df))
      tryCatch({
        p_gsea <- enrichplot::gseaplot2(result,
                                         geneSetID = 1:n_gsea,
                                         title = dot_title)
        save_figure(p_gsea, paste0(prefix, "_enrichment"),
                    width = 12, height = 8)
      }, error = function(e) {
        cat(sprintf("    Enrichment plot skipped: %s\n", conditionMessage(e)))
      })
    }
  }
}

# =====================================================================
# EXPERIMENT 1: TRIP4 TurboID vs WT — core TRIP4 interactome
# =====================================================================
cat("\n[1/4] TRIP4 TurboID vs Wild Type...\n")
exp_name <- "turbo_trip4_vs_wt"
if (exp_name %in% names(EXPERIMENTS)) {
  actual_name <- EXPERIMENTS[[exp_name]]
  if (actual_name %in% names(experiments)) {
    run_gsea(experiments[[actual_name]], "TRIP4 vs WT")
  } else {
    cat(sprintf("  Experiment '%s' not found in data\n", actual_name))
  }
} else {
  cat("  Experiment not defined in config\n")
}

# =====================================================================
# EXPERIMENT 2: RA effect — TRIP4+RA vs TRIP4 (RA-dependent changes)
# =====================================================================
cat("\n[2/4] RA effect (TRIP4+RA vs TRIP4)...\n")
exp_name <- "turbo_trip4_ra_vs_trip4"
if (exp_name %in% names(EXPERIMENTS)) {
  actual_name <- EXPERIMENTS[[exp_name]]
  if (actual_name %in% names(experiments)) {
    run_gsea(experiments[[actual_name]], "TRIP4+RA vs TRIP4")
  } else {
    cat(sprintf("  Experiment '%s' not found in data\n", actual_name))
  }
} else {
  cat("  RA effect experiment not defined — checking alternative names...\n")
  # Fallback: try to find RA vs base experiment
  candidates <- names(experiments)[grep("RA02_vs_TRIP4$", names(experiments))]
  if (length(candidates) > 0) {
    run_gsea(experiments[[candidates[1]]], "TRIP4+RA vs TRIP4")
  } else {
    cat("  No suitable RA experiment found\n")
  }
}

# =====================================================================
# EXPERIMENT 3: CHX vs DMSO — RNA-dependent changes
# =====================================================================
cat("\n[3/4] CHX vs DMSO (RNA-dependency)...\n")
exp_name <- "chx_vs_dmso"
if (exp_name %in% names(EXPERIMENTS)) {
  actual_name <- EXPERIMENTS[[exp_name]]
  if (actual_name %in% names(experiments)) {
    run_gsea(experiments[[actual_name]], "CHX vs DMSO")
  } else {
    cat(sprintf("  Experiment '%s' not found in data\n", actual_name))
  }
} else {
  cat("  CHX vs DMSO experiment not defined — checking data...\n")
  candidates <- names(experiments)[grep("CHX_vs.*DMSO", names(experiments))]
  if (length(candidates) > 0) {
    run_gsea(experiments[[candidates[1]]], "CHX vs DMSO")
  } else {
    cat("  No suitable CHX vs DMSO experiment found\n")
  }
}

# =====================================================================
# EXPERIMENT 4: CRAC RNA interactome
# =====================================================================
cat("\n[4/4] CRAC RNA interactome...\n")
crac_file <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")
if (file.exists(crac_file)) {
  crac_df <- load_proteomics_csv(
    crac_file,
    gene_col = CRAC_GENE_COL,
    log2fc_col = CRAC_LOG2FC_COL,
    padj_col = CRAC_PADJ_COL,
    pval_col = CRAC_PVAL_COL
  )
  crac_df$adj.P.Val <- crac_df$adj.P.Val  # already named correctly
  run_gsea(crac_df, "CRAC RNA Interactome")
} else {
  cat("  CRAC data file not found at:\n")
  cat(sprintf("    %s\n", crac_file))
}

cat("\n=========================================\n")
cat(" GSEA complete!\n")
cat("=========================================\n")
