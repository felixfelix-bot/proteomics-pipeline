###############################################################################
# 12_go_network_volcano.R
# Highlights proteins from GO enrichment results on the TurboID volcano.
# Like Lydia's approach: color the GO-associated proteins, no text labels.
#
# WHAT THIS PLOTS:
#   1. Runs GO enrichment on TurboID TRIP4 vs WT significant proteins
#   2. Collects all genes that appear in enriched GO terms
#   3. Plots the TurboID volcano with those genes highlighted in orange
#   4. Everything else stays gray
#   5. NO text labels — just color highlighting
#
# Usage:
#   make go-network-volcano
###############################################################################
cat("\n=========================================\n")
cat(" GO Network Volcano (Lydia-style)\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)

experiments <- load_all_experiments()
turbo_main <- "BK467_TRIP4_vs_BK467_WT"

if (!turbo_main %in% names(experiments)) {
  cat("ERROR:", turbo_main, "not found.\n")
  quit(status = 1)
}

df <- experiments[[turbo_main]]
df$neglog10p <- -log10(df$padj)

# ---- Run GO enrichment to get the gene network ----
cat("Running GO enrichment to find network genes...\n")
turbo_sig <- get_significant_genes(df)
cat(sprintf("  %d significant proteins\n", length(turbo_sig)))

universe <- unique(unlist(lapply(experiments, function(x) x$gene)))

go_genes <- character(0)

for (ont in c("BP", "MF", "CC")) {
  cat(sprintf("  [%s]...\n", ont))
  result <- tryCatch({
    enrichGO(
      gene          = turbo_sig,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      universe      = universe
    )
  }, error = function(e) NULL)

  if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
    # Extract all genes from enriched GO terms
    # The "geneID" column contains semicolon-separated gene symbols
    res_df <- as.data.frame(result)
    all_genes_in_terms <- unlist(strsplit(res_df$geneID, "/"))
    go_genes <- union(go_genes, unique(all_genes_in_terms))
    cat(sprintf("    %d enriched terms, %d unique genes\n",
                nrow(res_df), length(unique(all_genes_in_terms))))
  }
}

cat(sprintf("\n  Total unique genes in GO network: %d\n", length(go_genes)))

# ---- Build the volcano ----
cat("\nBuilding volcano plot...\n")

# Categories: GO network genes (orange) vs everything else (gray)
df$category <- "nonsig"
df$category[df$gene %in% go_genes] <- "go_network"
df$category <- factor(df$category, levels = c("go_network", "nonsig"))

GO_COLORS <- c(
  "go_network" = "#D55E00",   # Orange — in GO enrichment network
  "nonsig"     = "#D0D0D0"    # Gray
)

# NO labels — researcher said just highlight, don't label
p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.8) +
  ggplot2::scale_color_manual(
    values = GO_COLORS,
    labels = c(
      "go_network" = "In GO enrichment network",
      "nonsig"     = "Other proteins"
    ),
    name = NULL
  ) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID: GO Network Proteins Highlighted"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = "right",
    panel.grid.minor = ggplot2::element_blank()
  )

save_figure(p, "go_network_volcano_TRIP4_vs_WT",
            width = 8, height = 6)

# Also save the GO gene list
save_table(data.frame(gene = go_genes, category = "GO_network"),
           "go_network_genes")

cat("\n=========================================\n")
cat(" GO network volcano complete!\n")
cat("=========================================\n")
