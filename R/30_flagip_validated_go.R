###############################################################################
# 30_flagip_validated_go.R
# GO + KEGG enrichment analysis on proteins validated by BOTH C-Flag and N-Flag IP.
#
# These are the highest-confidence TRIP4 interactors — found significant in
# TurboID TRIP4 AND confirmed by both C-Flag IP AND N-Flag IP.
#
# Usage:
#   make flagip-validated-go
###############################################################################

cat("\n=========================================\n")
cat(" Flag IP Validated Proteins — GO + KEGG\n")
cat("=========================================\n\n")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
source("R/00_theme.R")

experiments <- load_all_experiments()

# ---- Load experiments ----
turbo_exp <- "BK467_TRIP4_vs_BK467_WT"
cflag_exp <- "BK516_Cflag_vs_BK516_Ctrl"
nflag_exp <- "BK516_Nflag_vs_BK516_Ctrl"

df_turbo <- find_experiment(experiments, turbo_exp)
df_cflag <- find_experiment(experiments, cflag_exp)
df_nflag <- find_experiment(experiments, nflag_exp)

if (is.null(df_turbo) || is.null(df_cflag) || is.null(df_nflag)) {
  cat("ERROR: Missing experiments.\n")
  quit(status = 1)
}

# ---- Extract significant genes (enriched only, positive log2FC >= threshold) ----
get_enriched <- function(df) {
  unique(df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC >= LOG2FC_CUTOFF & !is.na(df$gene)])
}

turbo_sig <- get_enriched(df_turbo)
cflag_sig <- get_enriched(df_cflag)
nflag_sig <- get_enriched(df_nflag)

# ---- Validated by BOTH C-Flag and N-Flag ----
validated_both <- Reduce(intersect, list(turbo_sig, cflag_sig, nflag_sig))

# Also: validated by C-Flag OR N-Flag (union)
validated_any <- intersect(turbo_sig, union(cflag_sig, nflag_sig))

cat(sprintf("  TurboID TRIP4-enriched:     %d\n", length(turbo_sig)))
cat(sprintf("  C-Flag IP enriched:         %d\n", length(cflag_sig)))
cat(sprintf("  N-Flag IP enriched:         %d\n", length(nflag_sig)))
cat(sprintf("  Validated by BOTH (C+N):    %d\n", length(validated_both)))
cat(sprintf("  Validated by EITHER (C|N):  %d\n", length(validated_any)))

if (length(validated_both) > 0) {
  cat(sprintf("    Genes: %s\n", paste(validated_both, collapse = ", ")))
}

# ---- Export gene lists ----
write.table(data.frame(gene = validated_both),
            file = file.path(TABLE_DIR, "flagip_validated_both_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(data.frame(gene = validated_any),
            file = file.path(TABLE_DIR, "flagip_validated_any_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
save_table(data.frame(gene = validated_both, category = "validated_both"), "flagip_validated_both")
save_table(data.frame(gene = validated_any, category = "validated_any"), "flagip_validated_any")

# ---- Universe = all TurboID TRIP4 proteins ----
universe <- unique(df_turbo$gene[!is.na(df_turbo$gene)])
cat(sprintf("  Universe: %d proteins\n", length(universe)))

# =====================================================================
# GO + KEGG enrichment
# =====================================================================
ONT_NAMES <- c("BP" = "Biological Process",
               "CC" = "Cellular Component",
               "MF" = "Molecular Function")

run_validated_go <- function(genes, set_label) {
  if (length(genes) < 5) {
    cat(sprintf("\n  [%s] Too few genes (%d). Skipping.\n", set_label, length(genes)))
    return(NULL)
  }

  cat(sprintf("\n  [%s] %d genes\n", set_label, length(genes)))

  # ---- GO enrichment ----
  for (ont in c("BP", "CC", "MF")) {
    cat(sprintf("    [GO %s] enrichGO... ", ont))
    ego <- tryCatch({
      enrichGO(gene = genes, universe = universe,
               OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               minGSSize = 2, maxGSSize = 5000)
    }, error = function(e) { cat("ERROR\n"); NULL })

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      cat("No terms\n"); next
    }

    ego_s <- tryCatch(simplify(ego, cutoff = 0.7), error = function(e) ego)
    res <- as.data.frame(ego_s)
    cat(sprintf("%d terms\n", nrow(res)))

    save_table(res, sprintf("flagip_GO_%s_%s", set_label, ont))

    n_show <- min(15, nrow(res))
    p <- dotplot(ego_s, showCategory = n_show)
    p <- p +
      ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                     name = "p-adjusted value") +
      ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 12),
                                     breaks = make_size_breaks(p$data$Count)) +
      ggplot2::scale_y_discrete(labels = capitalize_first) +
      ggplot2::labs(
        title = sprintf("Flag IP %s — %s", gsub("_", " ", set_label), ONT_NAMES[ont]),
        x = "Gene Ratio"
      ) +
      theme_poster()
    save_figure(p, sprintf("flagip_GO_%s_%s_dotplot", set_label, ont),
                width = 18, height = max(14, n_show * 0.8))
  }

  # ---- KEGG enrichment ----
  cat("    [KEGG] enrichKEGG... ")
  ekegg <- tryCatch({
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    universe_entrez <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    ekegg <- enrichKEGG(
      gene = unique(entrez_map$ENTREZID),
      universe = unique(universe_entrez$ENTREZID),
      organism = "hsa", pAdjustMethod = "BH",
      pvalueCutoff = 0.05, minGSSize = 2, maxGSSize = 5000
    )
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    ekegg
  }, error = function(e) { cat("ERROR (network?)\n"); NULL })

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    cat(sprintf("%d pathways\n", nrow(as.data.frame(ekegg))))
    save_table(as.data.frame(ekegg), sprintf("flagip_KEGG_%s", set_label))
    n_show <- min(15, nrow(as.data.frame(ekegg)))
    p <- dotplot(ekegg, showCategory = n_show)
    p <- p +
      ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                     name = "p-adjusted value") +
      ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 12),
                                     breaks = make_size_breaks(p$data$Count)) +
      ggplot2::scale_y_discrete(labels = capitalize_first) +
      ggplot2::labs(title = sprintf("Flag IP %s — KEGG Pathways", gsub("_", " ", set_label)),
                    x = "Gene Ratio") +
      theme_poster()
    save_figure(p, sprintf("flagip_KEGG_%s_dotplot", set_label),
                width = 18, height = max(14, n_show * 0.8))
  } else {
    cat("No pathways\n")
  }
}

# Run for both gene sets
run_validated_go(validated_both, "validated_both")
run_validated_go(validated_any, "validated_any")

cat("\n=========================================\n")
cat(" Flag IP Validated GO complete!\n")
cat("=========================================\n")
