###############################################################################
# 16_string_network_targeted.R
# STRING protein-protein interaction network analysis for TRIP4 TurboID.
#
# WHAT THIS DOES:
#   1. Maps significant TRIP4-enriched proteins to the STRING database
#   2. Identifies "seed" proteins (high-confidence core hits)
#   3. Expands the network to find neighbors of seeds
#   4. Splits significant proteins into categories:
#      - Seed (highly enriched)
#      - In network (connected to seeds via STRING)
#      - Not in network (significant but isolated)
#   5. Saves gene lists for downstream GO enrichment
#   6. Generates network visualization
#
# Based on Lydia's 05_string_network.R, adapted for the targeted pipeline.
#
# Usage:
#   make string-network
###############################################################################
cat("\n=========================================\n")
cat(" STRING Network Analysis (Targeted)\n")
cat("=========================================\n\n")

library(STRINGdb)
library(igraph)

experiments <- load_all_experiments()

# ---- Initialize STRINGdb ----
cat("Initializing STRINGdb (version 12.0, human, score >= 400)...\n")
cat("  (first run downloads ~100MB, subsequent runs use cache)\n")

string_db <- STRINGdb$new(
  version         = STRING_VERSION,
  species         = STRING_TAXON,
  score_threshold = STRING_SCORE_THRESHOLD,
  input_directory = ""
)

# =====================================================================
# Helper: STRING network analysis for one experiment
# =====================================================================
# Adapted from Lydia's 05_string_network.R.
# Returns a list with: gene lists for each category, candidate table, graph.
analyze_string_network <- function(df, experiment_name,
                                   fc_threshold = 2, p_threshold = 6) {

  cat(sprintf("\n  [%s] Mapping proteins to STRING...\n", experiment_name))

  # ---- STEP 1: Map gene symbols to STRING IDs ----
  # STRINGdb$map() sometimes fails with "incorrect number of dimensions"
  # when the input data frame has extra columns it doesn't expect.
  # Fix: pass only the columns STRING needs (gene + log2FC + padj).
  map_input <- data.frame(
    gene = df$gene,
    log2FC = df$log2FC,
    padj = df$padj,
    stringsAsFactors = FALSE
  )
  map_input <- map_input[!is.na(map_input$gene) & map_input$gene != "", ]

  mapped <- tryCatch({
    string_db$map(map_input, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
    cat("    Retrying with simplified input...\n")
    # Fallback: just gene names, no extra columns
    simple_input <- data.frame(gene = unique(df$gene), stringsAsFactors = FALSE)
    simple_input <- simple_input[!is.na(simple_input$gene) & simple_input$gene != "", ]
    tryCatch({
      string_db$map(simple_input, "gene", removeUnmappedRows = TRUE)
    }, error = function(e2) {
      cat(sprintf("    Still failing: %s\n", conditionMessage(e2)))
      return(NULL)
    })
  })

  if (is.null(mapped) || nrow(mapped) == 0) {
    cat("    No proteins mapped to STRING\n")
    return(NULL)
  }
  cat(sprintf("    Mapped %d/%d proteins to STRING\n", nrow(mapped), nrow(df)))

  # ---- STEP 2: Define seed proteins (Lydia's strict criteria) ----
  # Seed = extremely confident hit:
  #   (a) log2FC > 2 AND -log10(padj) > 6  (4-fold change, p < 1e-6)
  #   OR
  #   (b) log2FC > 7 AND -log10(padj) > 2  (128-fold change, p < 0.01)
  high_mask <- !is.na(mapped$padj) & (
    (mapped$log2FC > fc_threshold & -log10(mapped$padj) > p_threshold) |
    (mapped$log2FC > 7 & -log10(mapped$padj) > 2)
  )

  seed_ids <- unique(mapped$gene[high_mask])
  cat(sprintf("    High-confidence seed hits: %d\n", length(seed_ids)))

  if (length(seed_ids) < 2) {
    cat("    Too few seed hits for network analysis\n")
    return(NULL)
  }

  seed_string_ids <- mapped$STRING_id[mapped$gene %in% seed_ids]

  # ---- STEP 3: Find neighbors of seed proteins ----
  cat("    Finding STRING neighbors of seeds...\n")
  nbrs_list <- lapply(seed_string_ids, function(x) string_db$get_neighbors(x))
  nbrs <- unique(unlist(nbrs_list))
  expanded_ids <- unique(c(seed_string_ids, nbrs))

  # ---- STEP 4: Get interactions among expanded set ----
  int_expanded <- string_db$get_interactions(expanded_ids)
  cat(sprintf("    Found %d interactions\n", nrow(int_expanded)))

  # Filter to edges touching at least one seed
  int_core <- int_expanded[
    int_expanded$from %in% seed_string_ids |
    int_expanded$to %in% seed_string_ids, ,
  ]

  if (nrow(int_core) < 2) {
    cat("    Too few interactions for network\n")
    return(NULL)
  }

  # ---- STEP 5: Build igraph network ----
  g <- graph_from_data_frame(
    int_core[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = unique(c(int_core$from, int_core$to)))
  )

  # Annotate vertices
  vertex_info <- mapped[, c("STRING_id", "gene", "log2FC", "padj")]
  vertex_info <- vertex_info[!duplicated(vertex_info$STRING_id), ]

  V(g)$gene_name <- vertex_info$gene[match(V(g)$name, vertex_info$STRING_id)]
  V(g)$log2FC <- vertex_info$log2FC[match(V(g)$name, vertex_info$STRING_id)]
  V(g)$padj <- vertex_info$padj[match(V(g)$name, vertex_info$STRING_id)]
  V(g)$is_seed <- V(g)$name %in% seed_string_ids

  # Count seed connections per node
  core_links <- sapply(V(g)$name, function(x) {
    nb <- neighbors(g, x)
    sum(nb$name %in% seed_string_ids)
  })
  V(g)$core_links <- core_links

  # ---- STEP 6: Flag proteins by network membership ----
  network_genes <- unique(V(g)$gene_name[!is.na(V(g)$gene_name)])

  # Significant genes NOT in network
  sig_genes <- get_significant_genes(df)
  sig_in_network <- intersect(sig_genes, network_genes)
  sig_not_in_network <- setdiff(sig_genes, network_genes)

  cat(sprintf("    Network members: %d proteins\n", length(network_genes)))
  cat(sprintf("    Significant in network: %d\n", length(sig_in_network)))
  cat(sprintf("    Significant NOT in network: %d\n", length(sig_not_in_network)))

  # ---- STEP 7: Save gene lists for GO enrichment ----
  save_table(data.frame(gene = seed_ids, category = "seed_high_confidence"),
             paste0("string_seeds_", experiment_name))
  save_table(data.frame(gene = sig_in_network, category = "in_network"),
             paste0("string_in_network_", experiment_name))
  save_table(data.frame(gene = sig_not_in_network, category = "not_in_network"),
             paste0("string_not_in_network_", experiment_name))

  # ---- STEP 8: Candidate interactors (non-seed network members) ----
  candidates <- data.frame(
    gene = V(g)$gene_name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    STRING_id = V(g)$name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    core_links = V(g)$core_links[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    log2FC = V(g)$log2FC[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    padj = V(g)$padj[!V(g)$is_seed & !is.na(V(g)$gene_name)]
  )
  candidates <- candidates[order(-candidates$core_links, -abs(candidates$log2FC)), ]
  save_table(candidates, paste0("string_candidates_", experiment_name))

  # ---- STEP 9: Network visualization ----
  cat("    Generating network plot...\n")

  V(g)$color <- ifelse(V(g)$is_seed, GLOBAL_COLORS[["enriched_up"]],
                       ifelse(!is.na(V(g)$gene_name), GLOBAL_COLORS[["known_ia"]], "grey70"))
  V(g)$size <- ifelse(V(g)$is_seed, 6, 3)
  V(g)$label <- ifelse(V(g)$is_seed | V(g)$core_links >= 2,
                       V(g)$gene_name, NA)

  set.seed(42)

  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(paste0("string_network_", experiment_name))
  versioned_name <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned_name, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned_name, ".pdf")

  grDevices::png(png_path, width = 12, height = 10, units = "in", res = FIG_DPI)
  plot(g,
       layout = layout_with_fr(g),
       vertex.frame.color = NA,
       edge.color = "grey80",
       edge.width = 0.5,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       main = paste("STRING Network:", experiment_name))
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 12, height = 10)
  plot(g,
       layout = layout_with_fr(g),
       vertex.frame.color = NA,
       edge.color = "grey80",
       edge.width = 0.5,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       main = paste("STRING Network:", experiment_name))
  grDevices::dev.off()

  cat(sprintf("    Saved: %s\n", basename(png_path)))

  return(list(
    candidates = candidates,
    seeds = seed_ids,
    in_network = sig_in_network,
    not_in_network = sig_not_in_network,
    graph = g
  ))
}

# =====================================================================
# Run STRING analysis on ALL experiments
# =====================================================================
all_string_results <- list()

for (exp_name in names(experiments)) {
  # Skip duplicate experiments (only process main copy)
  if (grepl("__[0-9]+$", exp_name)) next

  # Build display label (underscores → spaces for titles)
  exp_label <- gsub("_", " ", exp_name)

  exp_data <- experiments[[exp_name]]
  result <- analyze_string_network(exp_data, exp_label)

  if (!is.null(result)) {
    all_string_results[[exp_label]] <- result
  }
}

# Also process CRAC data
crac_file <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")
if (file.exists(crac_file)) {
  cat("\n--- CRAC RNA Interactome ---\n")
  crac_df <- load_proteomics_csv(
    crac_file,
    gene_col = CRAC_GENE_COL,
    log2fc_col = CRAC_LOG2FC_COL,
    padj_col = CRAC_PADJ_COL
  )
  result <- analyze_string_network(crac_df, "CRAC RNA Interactome")
  if (!is.null(result)) {
    all_string_results[["CRAC RNA Interactome"]] <- result
  }
}

# =====================================================================
# GO enrichment on STRING-categorized gene sets
# =====================================================================
cat("\n=========================================\n")
cat(" GO Enrichment by STRING Network Membership\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n\n", length(universe)))

ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

# For each experiment's STRING results, run GO on 3 gene sets:
#   1) In-network significant genes
#   2) Not-in-network significant genes
#   3) All significant genes (for comparison)
for (exp_label in names(all_string_results)) {
  result <- all_string_results[[exp_label]]

  gene_sets <- list(
    "in_network"     = result$in_network,
    "not_in_network" = result$not_in_network,
    "all_significant" = c(result$in_network, result$not_in_network)
  )

  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    cat(sprintf("\n--- %s / %s (%d genes) ---\n", exp_label, set_name, length(genes)))

    if (length(genes) < 5) {
      cat("  Skipped: fewer than 5 genes\n")
      next
    }

    for (ont in c("BP", "MF", "CC")) {
      cat(sprintf("  [%s] Running enrichGO...\n", ont))

      ego <- tryCatch({
        enrichGO(
          gene          = genes,
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

      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        cat("    No enriched terms found\n")
        next
      }

      ego <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
                      error = function(e) ego)
      res_df <- as.data.frame(ego)
      cat(sprintf("    %d terms after simplify\n", nrow(res_df)))

      prefix <- sanitize_filename(paste0("GO_STRING_", exp_label, "_", set_name, "_", ont))
      save_table(res_df, prefix)

      n_show <- min(20, nrow(res_df))
      fig_h <- max(7, n_show * 0.4)
      title_str <- paste0("GO: ", exp_label, " — ", set_name, " — ", ONT_LABELS[ont])

      p_dot <- dotplot(ego, showCategory = n_show, title = title_str) +
        ggplot2::labs(x = "Gene Ratio", color = "p-adjusted value") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                       axis.text.y = ggplot2::element_text(size = 7))
      save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_h)

      n_bar <- min(15, nrow(res_df))
      fig_hb <- max(7, n_bar * 0.5)
      p_bar <- barplot(ego, showCategory = n_bar, orderBy = "Count", title = title_str) +
        ggplot2::labs(color = "p-adjusted value", fill = "p-adjusted value") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                       axis.text.y = ggplot2::element_text(size = 7))
      save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = fig_hb)
    }
  }
}

cat("\n=========================================\n")
cat(" STRING network + GO analysis complete!\n")
cat("=========================================\n")
cat("\nOutputs:\n")
cat("  - Network visualization (PNG + PDF)\n")
cat("  - Gene lists: seeds, in-network, not-in-network (CSV)\n")
cat("  - Candidate interactors ranked by connectivity (CSV)\n")
cat("  - GO enrichment per network category (dotplot + barplot + CSV)\n")
cat("\nRun: make open-string-network\n")
