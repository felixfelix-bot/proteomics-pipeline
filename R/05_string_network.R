###############################################################################
# 05_string_network.R
# Maps significant proteins onto the STRING physical interaction network.
# Identifies candidate interactors connected to high-confidence hits.
#
# Based on Lydia's 20260521_stats_cutoffs_TurboID_mapToStringDBnetwork.R
#
# Requires: STRINGdb (already in install script), igraph
# Internet connection needed on first run (downloads STRING data).
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/05_string_network.R")
###############################################################################

cat("\n=========================================\n")
cat(" STRING Network Analysis\n")
cat("=========================================\n\n")

library(STRINGdb)
library(igraph)

# ---- Load data ----
experiments <- load_all_experiments()

# ---- Initialize STRINGdb ----
cat("Initializing STRINGdb (version 12.0, human, score >= 400)...\n")
cat("  (first run downloads ~100MB, subsequent runs use cache)\n")

string_db <- STRINGdb$new(
  version = STRING_VERSION,
  species = STRING_TAXON,
  score_threshold = STRING_SCORE_THRESHOLD,
  input_directory = ""
)

# ---- Helper: Run STRING network analysis for one experiment ----
analyze_string_network <- function(df, experiment_name,
                                   fc_threshold = 2, p_threshold = 6) {
  cat(sprintf("\n  [%s] Mapping proteins to STRING...\n", experiment_name))

  # Map gene symbols to STRING IDs
  mapped <- tryCatch({
    string_db$map(df, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(mapped) || nrow(mapped) == 0) {
    cat("    No proteins mapped to STRING\n")
    return(NULL)
  }
  cat(sprintf("    Mapped %d/%d proteins to STRING\n", nrow(mapped), nrow(df)))

  # Define high-confidence core hits
  # (Lydia's criteria: log2FC > 2 AND -log10(padj) > 6, OR log2FC > 7 AND -log10 > 2)
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

  # Get neighbors of seed proteins (1-hop expansion)
  cat("    Finding network neighbors...\n")
  nbrs_list <- lapply(seed_string_ids, function(x) string_db$get_neighbors(x))
  nbrs <- unique(unlist(nbrs_list))
  expanded_ids <- unique(c(seed_string_ids, nbrs))

  # Get interactions among expanded set
  int_expanded <- string_db$get_interactions(expanded_ids)
  cat(sprintf("    Found %d interactions\n", nrow(int_expanded)))

  # Keep only edges touching at least one seed protein
  int_core <- int_expanded[
    int_expanded$from %in% seed_string_ids |
    int_expanded$to %in% seed_string_ids, ,
  ]

  if (nrow(int_core) < 2) {
    cat("    Too few interactions for network\n")
    return(NULL)
  }

  # Build igraph network
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

  # Count connections to seed set
  core_links <- sapply(V(g)$name, function(x) {
    nb <- neighbors(g, x)
    sum(nb$name %in% seed_string_ids)
  })
  V(g)$core_links <- core_links

  # Flag proteins in the network (for volcano plot overlay)
  df$inNetwork <- FALSE
  network_genes <- unique(V(g)$gene_name[!is.na(V(g)$gene_name)])
  df$inNetwork[df$gene %in% network_genes] <- TRUE

  n_in_net <- sum(df$inNetwork[df$gene %in% mapped$gene])
  cat(sprintf("    Proteins in STRING network: %d\n", n_in_net))

  # Rank candidate interactors by connectivity
  candidates <- data.frame(
    gene = V(g)$gene_name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    STRING_id = V(g)$name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    core_links = V(g)$core_links[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    log2FC = V(g)$log2FC[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    padj = V(g)$padj[!V(g)$is_seed & !is.na(V(g)$gene_name)]
  )
  candidates <- candidates[order(-candidates$core_links, -abs(candidates$log2FC)), ]

  # ---- Plot network ----
  cat("    Generating network plot...\n")

  # Color: seeds = red, candidates = teal, other = grey
  V(g)$color <- ifelse(V(g)$is_seed, "#E64B35",
                       ifelse(!is.na(V(g)$gene_name), "#1b9e77", "grey70"))
  V(g)$size <- ifelse(V(g)$is_seed, 6, 3)
  V(g)$label <- ifelse(V(g)$is_seed | V(g)$core_links >= 2,
                       V(g)$gene_name, NA)

  set.seed(42)  # Reproducible layout
  png(file.path(FIGURE_DIR, paste0("string_network_", experiment_name, ".png")),
      width = 12, height = 10, units = "in", res = FIG_DPI)
  plot(g,
       layout = layout_with_fr(g),
       vertex.frame.color = NA,
       edge.color = "grey80",
       edge.width = 0.5,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       main = paste("STRING Network:", experiment_name))
  dev.off()
  cat(sprintf("    Saved: string_network_%s.png\n", experiment_name))

  # ---- Volcano with network overlay ----
  toPlot <- df
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )
  toPlot$category[!toPlot$inNetwork & toPlot$category == "TRUE"] <- "TRUE"
  toPlot$category[toPlot$inNetwork] <- "inNetwork"
  toPlot$category[high_mask] <- "high"
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  label_data <- toPlot[toPlot$category %in% c("high", "inNetwork"), ]

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.2) +
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, max.overlaps = 20, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = paste("STRING Network Overlay:", experiment_name)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_figure(p, paste0("volcano_string_", experiment_name), width = 8, height = 6)

  return(list(
    candidates = candidates,
    in_network = network_genes,
    graph = g
  ))
}

# =====================================================================
# MAIN: Run STRING analysis for key experiments
# =====================================================================

key_experiments <- c("turbo_trip4_vs_wt", "flag_cflag_vs_ctrl",
                     "turbo_RA_vs_wt", "flag_RA_cflag_vs_cflag")

for (exp_name in key_experiments) {
  if (exp_name %in% names(experiments)) {
    result <- analyze_string_network(experiments[[exp_name]], exp_name)
    if (!is.null(result)) {
      save_table(result$candidates, paste0("string_candidates_", exp_name))
      save_table(data.frame(gene = result$in_network, category = "in_network"),
                 paste0("string_innetwork_", exp_name))
    }
  } else {
    cat(sprintf("\n  [%s] Skipping (not found in data)\n", exp_name))
  }
}

cat("\n=========================================\n")
cat(" STRING network analysis complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
