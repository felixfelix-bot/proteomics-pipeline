###############################################################################
# 04_go_enrichment.R
# Gene Ontology enrichment analysis using clusterProfiler.
# Performs ORA (over-representation analysis) and GSEA.
# Includes visualization: dotplot, barplot, cnetplot, and GO term reduction.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/04_go_enrichment.R")
###############################################################################

cat("\n=========================================\n")
cat(" GO Enrichment Analysis\n")
cat("=========================================\n\n")

# ---- Load annotation database ----
cat("Loading annotation database (org.Hs.eg.db)...\n")
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# ---- Load all experiment data ----
experiments <- list()
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  name <- tools::file_path_sans_ext(basename(f))
  if (name == "known_interactors") next
  experiments[[name]] <- load_proteomics_csv(f)
}

# ---- Define universe (all detected proteins across experiments) ----
universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n", length(universe)))

# ---- Extract significant gene sets ----
gene_sets <- lapply(experiments, get_significant_genes)

# =====================================================================
# Helper: Run ORA for one gene set and one ontology
# =====================================================================
run_ora <- function(genes, ontology, universe, experiment_name) {
  if (length(genes) < 5) {
    cat(sprintf("  [%s/%s] Skipped: only %d significant genes (need >=5)\n",
                experiment_name, ontology, length(genes)))
    return(NULL)
  }

  cat(sprintf("  [%s/%s] Running enrichGO with %d genes...\n",
              experiment_name, ontology, length(genes)))

  result <- tryCatch({
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      pAdjustMethod = GO_PADJUST_METHOD,
      pvalueCutoff  = GO_PVALUE_CUTOFF,
      qvalueCutoff  = GO_QVALUE_CUTOFF,
      universe      = universe
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat(sprintf("    No enriched terms found\n"))
    return(NULL)
  }

  cat(sprintf("    Found %d enriched GO terms\n", nrow(as.data.frame(result))))
  return(result)
}

# =====================================================================
# Helper: Visualize ORA results
# =====================================================================
visualize_ora <- function(ego, experiment_name, ontology) {
  if (is.null(ego)) return()

  prefix <- paste0("GO_", experiment_name, "_", ontology)
  res_df <- as.data.frame(ego)

  # Save results table
  save_table(res_df, prefix)

  # Dotplot (top 20)
  p_dot <- dotplot(ego, showCategory = 20, title = paste(experiment_name, ontology)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_dot, paste0(prefix, "_dotplot"), width = 9, height = 7)

  # Barplot (top 20)
  p_bar <- barplot(ego, showCategory = 20, title = paste(experiment_name, ontology)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_bar, paste0(prefix, "_barplot"), width = 9, height = 7)

  # Cnetplot (gene-concept network) — top 10 categories
  if (nrow(res_df) >= 3) {
    p_cnet <- tryCatch({
      cnetplot(ego, showCategory = 10,
               categorySize = "pvalue",
               foldChange = NULL) +
        ggplot2::ggtitle(paste(experiment_name, ontology))
    }, error = function(e) NULL)

    if (!is.null(p_cnet)) {
      save_figure(p_cnet, paste0(prefix, "_cnetplot"), width = 12, height = 10)
    }
  }
}

# =====================================================================
# MAIN: Run ORA for each experiment and each ontology
# =====================================================================
cat("\n[1/3] Running ORA (over-representation analysis)...\n")

all_ora_results <- list()

for (exp_name in names(gene_sets)) {
  cat(sprintf("\n  --- %s ---\n", exp_name))
  sig_genes <- gene_sets[[exp_name]]
  cat(sprintf("  Significant genes: %d\n", length(sig_genes)))

  for (ont in GO_ONTOLOGIES) {
    ego <- run_ora(sig_genes, ont, universe, exp_name)
    visualize_ora(ego, exp_name, ont)
    if (!is.null(ego)) {
      result_key <- paste(exp_name, ont, sep = "_")
      all_ora_results[[result_key]] <- ego
    }
  }
}

# =====================================================================
# GSEA (Gene Set Enrichment Analysis) — uses ranked gene list
# =====================================================================
cat("\n[2/3] Running GSEA (gene set enrichment analysis)...\n")

run_gsea_for_experiment <- function(df, experiment_name) {
  cat(sprintf("\n  --- %s GSEA ---\n", experiment_name))

  # Rank genes by signed -log10(padj): up-regulated at top, down at bottom
  df$metric <- -log10(df$padj + 1e-10) * sign(df$log2FC)
  df <- df[!is.na(df$metric), ]
  df <- df[order(df$metric, decreasing = TRUE), ]

  gene_list <- df$metric
  names(gene_list) <- df$gene

  # Remove duplicates and zeros
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- gene_list[gene_list != 0]

  gsea_result <- tryCatch({
    gseGO(
      geneList     = gene_list,
      OrgDb        = org.Hs.eg.db,
      keyType      = "SYMBOL",
      ont          = "BP",       # GSEA typically on Biological Process
      minGSSize    = 10,
      maxGSSize    = 500,
      pvalueCutoff = GO_PVALUE_CUTOFF,
      verbose      = FALSE
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    cat("    No enriched gene sets found\n")
    return(NULL)
  }

  cat(sprintf("    Found %d enriched gene sets\n", nrow(as.data.frame(gsea_result))))

  # Save table
  save_table(as.data.frame(gsea_result), paste0("GSEA_", experiment_name, "_BP"))

  # Dotplot
  p_gsea <- dotplot(gsea_result, showCategory = 20,
                    title = paste("GSEA:", experiment_name, "BP")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  save_figure(p_gsea, paste0("GSEA_", experiment_name, "_BP_dotplot"), width = 10, height = 8)

  # Ridgeplot (distribution of fold changes within enriched sets)
  if (nrow(as.data.frame(gsea_result)) >= 3) {
    p_ridge <- tryCatch({
      ridgeplot(gsea_result) +
        ggplot2::ggtitle(paste("GSEA Ridge:", experiment_name))
    }, error = function(e) NULL)

    if (!is.null(p_ridge)) {
      save_figure(p_ridge, paste0("GSEA_", experiment_name, "_BP_ridge"), width = 10, height = 8)
    }
  }

  return(gsea_result)
}

for (exp_name in names(experiments)) {
  run_gsea_for_experiment(experiments[[exp_name]], exp_name)
}

# =====================================================================
# GO Term Reduction (rrvgo) — collapse redundant enriched terms
# =====================================================================
cat("\n[3/3] Reducing GO term redundancy (rrvgo)...\n")

library(rrvgo)

for (key in names(all_ora_results)) {
  ego <- all_ora_results[[key]]
  res_df <- as.data.frame(ego)

  if (nrow(res_df) < 5) {
    cat(sprintf("  [%s] Skipped: only %d terms (need >=5)\n", key, nrow(res_df)))
    next
  }

  ont <- strsplit(key, "_")[[1]]
  ontology <- ont[length(ont)]

  cat(sprintf("  [%s] Reducing %d GO terms...\n", key, nrow(res_df)))

  sim_matrix <- tryCatch({
    calculateSimMatrix(
      res_df$ID,
      orgdb = org.Hs.eg.db,
      ont = ontology,
      method = "Rel"
    )
  }, error = function(e) {
    cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(sim_matrix)) next

  scores <- setNames(-log10(res_df$p.adjust), res_df$ID)
  reduced <- reduceSimMatrix(sim_matrix, scores, threshold = 0.7,
                              orgdb = org.Hs.eg.db)

  save_table(as.data.frame(reduced), paste0("reduced_", key))

  # Treemap
  p_tree <- tryCatch({
    treemapPlot(reduced) +
      ggplot2::ggtitle(paste("Reduced GO:", key))
  }, error = function(e) NULL)

  if (!is.null(p_tree)) {
    save_figure(p_tree, paste0("reduced_", key, "_treemap"), width = 10, height = 8)
  }

  cat(sprintf("    Reduced to %d representative terms\n", nrow(reduced)))
}

cat("\n=========================================\n")
cat(" GO enrichment analysis complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
