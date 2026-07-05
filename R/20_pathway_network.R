###############################################################################
# 20_pathway_network.R
# STRING-style protein interaction NETWORK MAPS for enriched GO pathways.
#
# WHAT THIS DOES:
#   The user wants to see actual interaction networks (nodes + edges) for the
#   proteins involved in each enriched GO pathway — not just dot/bar plots.
#
# FOR EACH EXPERIMENT:
#   1. Run ORA (enrichGO) to find enriched pathways (BP, MF, CC)
#   2. For the top 3 enriched GO terms, extract the gene lists
#   3. Map those genes to STRINGdb and retrieve interactions
#   4. Build an igraph network visualization
#   5. Color nodes by logFC (red=up, blue=down), size by -log10(padj)
#   6. Save each pathway network as a PDF
#
# This is the "STRING database style network map" Dr Aruna requested.
#
# Usage:
#   make pathway-network
#   make open-pathway-network
###############################################################################
cat("\n=========================================\n")
cat(" PATHWAY NETWORK MAPS\n")
cat(" STRING-style interaction networks for\n")
cat(" proteins in enriched GO pathways\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(STRINGdb)
library(igraph)
library(grid)

# ---- Initialize STRINGdb ----
cat("Initializing STRINGdb...\n")
string_db <- STRINGdb$new(
  version         = STRING_VERSION,
  species         = STRING_TAXON,
  score_threshold = 400,
  input_directory = ""
)

ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

experiments <- load_all_experiments()
all_names <- names(experiments)
universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n\n", length(universe)))

# ---- Helper: build a STRING network for genes in a GO term ----
build_go_network <- function(pathway_id, pathway_desc, go_type, gene_list,
                              experiment_df, exp_label, prefix) {
  cat(sprintf("  Building network for: %s (%d genes)\n",
              pathway_desc, length(gene_list)))

  # Need at least 3 genes for a meaningful network
  if (length(gene_list) < 3) {
    cat("    Too few genes for network\n")
    return(NULL)
  }

  # Create a temporary data frame with the pathway genes
  pathway_df <- data.frame(
    gene = gene_list,
    stringsAsFactors = FALSE
  )

  # Map to STRING
  mapped <- tryCatch({
    string_db$map(pathway_df, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    STRING mapping error: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(mapped) || nrow(mapped) < 3) {
    cat("    Too few mapped to STRING for network\n")
    return(NULL)
  }

  # Get interactions among mapped proteins
  int <- tryCatch({
    string_db$get_interactions(mapped$STRING_id)
  }, error = function(e) {
    cat(sprintf("    STRING interaction error: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(int) || nrow(int) < 2) {
    cat("    Too few interactions for network\n")
    return(NULL)
  }

  # Build igraph
  g <- graph_from_data_frame(
    int[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = unique(c(int$from, int$to)))
  )

  # Annotate vertices with logFC from the experiment data
  exp_map <- experiment_df[, c("gene", "logFC", "padj")]
  exp_map <- exp_map[!duplicated(exp_map$gene), ]

  v_genes <- mapped$gene[match(V(g)$name, mapped$STRING_id)]
  V(g)$gene_name <- ifelse(is.na(v_genes), V(g)$name, v_genes)
  V(g)$logFC <- exp_map$logFC[match(V(g)$gene_name, exp_map$gene)]
  V(g)$padj <- exp_map$padj[match(V(g)$gene_name, exp_map$gene)]
  V(g)$logFC[is.na(V(g)$logFC)] <- 0
  V(g)$padj[is.na(V(g)$padj)] <- 1
  V(g)$size <- pmin(25, pmax(5, -log10(V(g)$padj + 1e-10) * 3))

  # Color: red = up, blue = down, grey = neutral
  node_colors <- rep("#D0D0D0", vcount(g))
  node_colors[V(g)$logFC > 0.5] <- "#D55E00"  # up-regulated
  node_colors[V(g)$logFC < -0.5] <- "#0072B2"  # down-regulated
  V(g)$color <- node_colors

  # Label only well-connected or high-significance nodes
  deg <- degree(g)
  label_cutoff <- quantile(deg, probs = 0.7, na.rm = TRUE)
  label_nodes <- deg >= label_cutoff | (-log10(V(g)$padj + 1e-10)) > 2
  V(g)$label <- ifelse(label_nodes, V(g)$gene_name, "")

  # ---- Layout and plot ----
  l <- tryCatch({
    layout_with_fr(g, niter = 500)
  }, error = function(e) {
    layout_nicely(g)
  })

  # Save as PDF
  network_prefix <- paste0(prefix, "_network")
  title_line <- sprintf("%s\n(%s — %s)", pathway_desc, exp_label, go_type)

  pdf_width <- 10
  pdf_height <- 8

  pdf_path <- safe_filepath(FIGURE_DIR, paste0(network_prefix, "_", get_git_hash()), ".pdf")
  grDevices::pdf(pdf_path, width = pdf_width, height = pdf_height)

  # Set up plotting area
  par(mar = c(1, 1, 3, 1))

  # Plot the network
  plot(g,
       layout = l,
       vertex.size = V(g)$size,
       vertex.color = V(g)$color,
       vertex.frame.color = "#333333",
       vertex.frame.width = 0.5,
       vertex.label = V(g)$label,
       vertex.label.color = "#111111",
       vertex.label.font = 1,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       edge.color = "#999999",
       edge.width = 0.5,
       edge.curved = 0.1,
       main = title_line,
       cex.main = 0.8
  )

  # Add legend
  legend("bottomright",
         legend = c("Up-regulated (logFC > 0.5)", "Down-regulated (logFC < -0.5)", "No change"),
         fill = c("#D55E00", "#0072B2", "#D0D0D0"),
         cex = 0.6,
         box.lty = 0,
         bg = "white"
  )

  dev.off()
  cat(sprintf("  Saved: %s\n", basename(pdf_path)))

  # Also save a PNG version
  png_path <- safe_filepath(FIGURE_DIR, paste0(network_prefix, "_", get_git_hash()), ".png")
  grDevices::png(png_path, width = pdf_width, height = pdf_height, units = "in", res = 200)
  par(mar = c(1, 1, 3, 1))
  plot(g,
       layout = l,
       vertex.size = V(g)$size,
       vertex.color = V(g)$color,
       vertex.frame.color = "#333333",
       vertex.frame.width = 0.5,
       vertex.label = V(g)$label,
       vertex.label.color = "#111111",
       vertex.label.font = 1,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       edge.color = "#999999",
       edge.width = 0.5,
       edge.curved = 0.1,
       main = title_line,
       cex.main = 0.8
  )
  legend("bottomright",
         legend = c("Up-regulated (logFC > 0.5)", "Down-regulated (logFC < -0.5)", "No change"),
         fill = c("#D55E00", "#0072B2", "#D0D0D0"),
         cex = 0.6,
         box.lty = 0,
         bg = "white"
  )
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(png_path)))

  return(g)
}

# ---- Helper: run ORA + build networks for one experiment ----
analyze_experiment_pathways <- function(df, exp_name) {
  cat(sprintf("\n========== %s (%d proteins) ==========\n", exp_name, nrow(df)))

  sig_genes <- get_significant_genes(df)
  cat(sprintf("  Significant genes: %d\n", length(sig_genes)))

  if (length(sig_genes) < 10) {
    cat("  Too few significant genes for pathway analysis\n")
    return()
  }

  safe_exp <- sanitize_filename(gsub("[ _]", "_", exp_name))

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("\n  [%s] Running enrichGO...\n", ont))

    result <- tryCatch({
      enrichGO(
        gene          = sig_genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        universe      = universe
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      cat("    No enriched terms\n")
      next
    }

    res_df <- as.data.frame(result)
    cat(sprintf("    Found %d enriched terms\n", nrow(res_df)))

    # Take top 3 enriched pathways by p.adjust
    n_pathways <- min(3, nrow(res_df))

    for (j in 1:n_pathways) {
      pathway_id <- res_df$ID[j]
      pathway_desc <- res_df$Description[j]
      pathway_genes <- strsplit(as.character(res_df$geneID[j]), "/")[[1]]

      prefix <- paste0("pathway_net_", safe_exp, "_", ont, "_g", j)
      cat(sprintf("\n  [%s/%s] Top %d/%d: %s\n",
                  ONT_LABELS[ont], exp_name, j, n_pathways, pathway_desc))

      build_go_network(
        pathway_id    = pathway_id,
        pathway_desc  = pathway_desc,
        go_type       = ONT_LABELS[ont],
        gene_list     = pathway_genes,
        experiment_df = df,
        exp_label     = exp_name,
        prefix        = prefix
      )
    }
  }
}

# =====================================================================
# MAIN: Run pathway networks for ALL discovered experiments
# =====================================================================
for (i in seq_along(all_names)) {
  exp_name <- all_names[i]
  cat(sprintf("\n[%d/%d] %s\n", i, length(all_names), exp_name))
  analyze_experiment_pathways(experiments[[exp_name]], exp_name)
}

# =====================================================================
# ALSO: CRAC data
# =====================================================================
crac_file <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")
if (file.exists(crac_file)) {
  cat("\n[+] CRAC RNA Interactome\n")
  crac_df <- load_proteomics_csv(
    crac_file,
    gene_col = CRAC_GENE_COL,
    log2fc_col = CRAC_LOG2FC_COL,
    padj_col = CRAC_PADJ_COL
  )
  analyze_experiment_pathways(crac_df, "CRAC_RNA_Interactome")
} else {
  cat("\n  CRAC data file not found (skipping)\n")
}

cat("\n=========================================\n")
cat(" Pathway network maps complete!\n")
cat("=========================================\n")
