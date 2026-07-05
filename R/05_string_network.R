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

# === WHAT IS STRING? =========================================================
# STRING (https://string-db.org) is a database of KNOWN and PREDICTED
# protein-protein interactions (PPIs). It collects evidence from:
#   - Experiments (lab-validated interactions)
#   - Databases (curated pathway databases like KEGG, Reactome)
#   - Text mining (published papers)
#   - Co-expression (genes turned on/off together)
#   Each interaction gets a "confidence score" (0-1000). Higher = more evidence.
# By mapping our significant proteins onto STRING, we can see if they form
# connected biological networks — and discover "neighbor" proteins that
# interact with our hits but weren't detected, which may be candidate
# interactors we missed.
# =============================================================================

# cat() prints text to the console (like print() but without line numbers/quotes).
# These lines create a visual header so you know which script is running.
cat("\n=========================================\n")
cat(" STRING Network Analysis\n")
cat("=========================================\n\n")

# library() loads an R package (must be installed first).
# STRINGdb: R package to query the STRING protein interaction database.
# igraph: R package for network/graph analysis — lets us build, visualize,
#         and analyze networks of nodes (proteins) connected by edges (interactions).
library(STRINGdb)
library(igraph)

# ---- Load data ----
# load_all_experiments() is a custom helper (defined in utils.R) that reads
# all the proteomics result files and returns them as a named list.
# Each list element is a data frame with columns: gene, log2FC, padj, etc.
# Think of a "list" as a container holding multiple data frames, accessed by name:
#   experiments[["turbo_trip4_vs_wt"]]  →  one data frame
experiments <- load_all_experiments()

# ---- Initialize STRINGdb ----
# This creates the connection to the STRING database.
cat("Initializing STRINGdb (version 12.0, human, score >= 400)...\n")
cat("  (first run downloads ~100MB, subsequent runs use cache)\n")

# STRINGdb$new() creates a new STRINGdb "object" (an R6 class — think of it as
# a database connection object with built-in methods you call using $).
# Parameters:
#   version         = which version of STRING to use (12.0 is recent)
#   species         = taxonomy ID for the organism (e.g., 9606 for human).
#                     STRING_TAXON is defined in 01_config.R.
#   score_threshold = minimum confidence score for interactions (400 = "medium").
#                     Only interactions scoring ≥ 400 are considered. This filters
#                     out low-confidence predictions. STRING_SCORE_THRESHOLD is in config.
#   input_directory = where to cache downloaded data ("" = temp location).
string_db <- STRINGdb$new(
  version = STRING_VERSION,
  species = STRING_TAXON,
  score_threshold = STRING_SCORE_THRESHOLD,
  input_directory = ""
)

# ---- Helper: Run STRING network analysis for one experiment ----
# This defines a FUNCTION (a reusable block of code).
# It takes:
#   df              = a data frame with one experiment's results (gene, log2FC, padj)
#   experiment_name = a label string (used in plot titles & file names)
#   fc_threshold    = minimum log2 fold change for "high confidence" hits (default 2)
#   p_threshold     = minimum -log10(padj) for "high confidence" hits (default 6)
analyze_string_network <- function(df, experiment_name,
                                   fc_threshold = 2, p_threshold = 6) {

  # sprintf() formats a string (like printf in other languages).
  # %s is replaced by the experiment name.
  cat(sprintf("\n  [%s] Mapping proteins to STRING...\n", experiment_name))

  # ---- STEP 1: Map gene symbols to STRING IDs ----
  # STRING uses its own internal identifiers (e.g., "9606.ENSP00000000233").
  # string_db$map() takes our data frame and adds a "STRING_id" column by
  # looking up each gene symbol (e.g., "TRIP4") in STRING's database.
  #   "gene"               = the column name in df containing gene symbols
  #   removeUnmappedRows   = TRUE drops any genes STRING doesn't recognize
  # tryCatch() is R's error handler: it tries the first block, and if an error
  # occurs, it runs the error handler instead of crashing the whole script.
  # This makes the script robust — one bad experiment won't stop everything.
  mapped <- tryCatch({
    string_db$map(df, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
    return(NULL)  # Return NULL (nothing) on error
  })

  # If mapping failed or returned zero rows, skip this experiment.
  # || means "OR", is.na() checks for missing values, nrow() counts rows.
  if (is.null(mapped) || nrow(mapped) == 0) {
    cat("    No proteins mapped to STRING\n")
    return(NULL)
  }
  # Report how many of our proteins STRING recognized.
  # %d is replaced by an integer value.
  cat(sprintf("    Mapped %d/%d proteins to STRING\n", nrow(mapped), nrow(df)))

  # ---- STEP 2: Define "seed" proteins (high-confidence core hits) ----
  # "Seed" proteins are our most confident hits — the ones we're most sure about.
  # We'll build a network around them. Lydia's criteria define a seed as either:
  #   (a) Very significant AND good fold change: log2FC > 2 AND -log10(padj) > 6
  #       (i.e., >4-fold change AND p-adj < 0.000001)
  #   OR
  #   (b) Huge fold change with modest significance: log2FC > 7 AND -log10(padj) > 2
  #       (i.e., >128-fold change AND p-adj < 0.01)
  # is.na(mapped$padj) checks if the adjusted p-value is missing.
  # The ! means "NOT", so !is.na() means "the p-value exists".
  # & means "AND" (both must be true). The parentheses group the logic.
  high_mask <- !is.na(mapped$padj) & (
    (mapped$log2FC > fc_threshold & -log10(mapped$padj) > p_threshold) |
    (mapped$log2FC > 7 & -log10(mapped$padj) > 2)
  )

  # Extract the unique gene names of seed proteins.
  # mapped$gene[high_mask] selects gene names where high_mask is TRUE.
  # unique() removes duplicates.
  seed_ids <- unique(mapped$gene[high_mask])
  cat(sprintf("    High-confidence seed hits: %d\n", length(seed_ids)))

  # We need at least 2 seeds to form a meaningful network.
  if (length(seed_ids) < 2) {
    cat("    Too few seed hits for network analysis\n")
    return(NULL)
  }

  # Get the STRING IDs corresponding to our seed gene names.
  # %in% checks membership: "is this value in that set?"
  seed_string_ids <- mapped$STRING_id[mapped$gene %in% seed_ids]

  # ---- STEP 3: Find "neighbors" of seed proteins (network expansion) ----
  # "Neighbors" are proteins that directly interact with a seed protein
  # in the STRING database (one "hop" away). These are CANDIDATE interactors
  # — proteins connected to our hits that we may not have detected directly.
  # string_db$get_neighbors() returns all STRING IDs that interact with a given ID.
  # lapply() applies a function to each element of a list/vector:
  #   for each seed_string_id, get its neighbors.
  nbrs_list <- lapply(seed_string_ids, function(x) string_db$get_neighbors(x))

  # unlist() flattens the list of neighbor vectors into one big vector.
  # unique() removes duplicates (a neighbor of multiple seeds appears once).
  nbrs <- unique(unlist(nbrs_list))

  # Combine seed IDs and neighbor IDs into one expanded set.
  # c() concatenates (combines) vectors.
  expanded_ids <- unique(c(seed_string_ids, nbrs))

  # ---- STEP 4: Get all interactions among the expanded set ----
  # string_db$get_interactions() returns a data frame of pairwise interactions
  # (columns "from" and "to" = STRING IDs of the two interacting proteins).
  int_expanded <- string_db$get_interactions(expanded_ids)
  cat(sprintf("    Found %d interactions\n", nrow(int_expanded)))

  # ---- STEP 5: Filter to keep only edges touching at least one seed ----
  # This focuses the network on connections relevant to our high-confidence hits.
  # The [, ,] syntax selects rows where the condition is TRUE, keeping all columns.
  int_core <- int_expanded[
    int_expanded$from %in% seed_string_ids |
    int_expanded$to %in% seed_string_ids, ,
  ]

  # Need at least 2 interactions to make a meaningful network.
  if (nrow(int_core) < 2) {
    cat("    Too few interactions for network\n")
    return(NULL)
  }

  # ---- STEP 6: Build an igraph network object ----
  # igraph is the main R package for network analysis.
  # graph_from_data_frame() converts a data frame into a graph (network) object.
  #   d        = a data frame where each row is an edge (from → to)
  #   directed = FALSE means interactions are bidirectional (undirected graph)
  #   vertices = a data frame listing all nodes (proteins) in the network
  # The result 'g' is a graph object we can analyze and plot.
  g <- graph_from_data_frame(
    int_core[, c("from", "to")],
    directed = FALSE,
    vertices = data.frame(name = unique(c(int_core$from, int_core$to)))
  )

  # ---- STEP 7: Annotate vertices (nodes) with our data ----
  # V(g) accesses all VERTICES (nodes/proteins) in graph g.
  # We add attributes (labels) to each node: gene name, fold change, p-value.
  # First, prepare a lookup table from our mapped data:
  vertex_info <- mapped[, c("STRING_id", "gene", "log2FC", "padj")]
  # Remove duplicate STRING IDs (keep one entry per protein):
  vertex_info <- vertex_info[!duplicated(vertex_info$STRING_id), ]

  # match(x, table) finds the position of each x in table.
  # So V(g)$name are STRING IDs, and we look up each one in vertex_info.
  # V(g)$gene_name = ... assigns a "gene_name" attribute to every vertex.
  V(g)$gene_name <- vertex_info$gene[match(V(g)$name, vertex_info$STRING_id)]
  V(g)$log2FC <- vertex_info$log2FC[match(V(g)$name, vertex_info$STRING_id)]
  V(g)$padj <- vertex_info$padj[match(V(g)$name, vertex_info$STRING_id)]
  # Mark which nodes are our high-confidence seeds:
  V(g)$is_seed <- V(g)$name %in% seed_string_ids

  # ---- STEP 8: Count how many seed connections each node has ----
  # For each node, count how many of its neighbors are seed proteins.
  # sapply() applies a function to each element and returns a vector of results.
  # neighbors(g, x) returns the neighbors (connected nodes) of node x in graph g.
  # nb$name gets their names, and %in% checks which are seeds.
  core_links <- sapply(V(g)$name, function(x) {
    nb <- neighbors(g, x)
    sum(nb$name %in% seed_string_ids)
  })
  # Store this count as a vertex attribute:
  V(g)$core_links <- core_links

  # ---- STEP 9: Flag proteins in the network (for volcano plot overlay) ----
  # Add a new column to the original data frame marking network membership.
  df$inNetwork <- FALSE
  # Get gene names of all nodes that have a gene name assigned:
  network_genes <- unique(V(g)$gene_name[!is.na(V(g)$gene_name)])
  # Mark proteins in df that appear in the network:
  df$inNetwork[df$gene %in% network_genes] <- TRUE

  # Count how many of our mapped proteins are in the network:
  n_in_net <- sum(df$inNetwork[df$gene %in% mapped$gene])
  cat(sprintf("    Proteins in STRING network: %d\n", n_in_net))

  # ---- STEP 10: Rank candidate interactors by connectivity ----
  # "Candidates" are NON-seed proteins in the network — potential interactors
  # we discovered through STRING but didn't detect directly.
  # We rank them by: (1) number of connections to seeds, (2) fold change magnitude.
  # data.frame() creates a new data frame from the selected columns.
  # The condition !V(g)$is_seed excludes seeds; !is.na(...) excludes unnamed nodes.
  candidates <- data.frame(
    gene = V(g)$gene_name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    STRING_id = V(g)$name[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    core_links = V(g)$core_links[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    log2FC = V(g)$log2FC[!V(g)$is_seed & !is.na(V(g)$gene_name)],
    padj = V(g)$padj[!V(g)$is_seed & !is.na(V(g)$gene_name)]
  )

  # order() returns the sort order; the leading - means DESCENDING.
  # So we sort by core_links (most connections first), then by |log2FC|.
  # candidates[order(...), ] reorders the rows.
  candidates <- candidates[order(-candidates$core_links, -abs(candidates$log2FC)), ]

  # ---- STEP 11: Plot the network ----
  cat("    Generating network plot...\n")

  # Assign visual attributes to vertices for plotting:
  # Color: seeds = red (#E64B35), named candidates = teal (#1b9e77), unnamed = grey
  # ifelse() is a vectorized if-else: ifelse(condition, value_if_true, value_if_false)
  V(g)$color <- ifelse(V(g)$is_seed, "#E64B35",
                       ifelse(!is.na(V(g)$gene_name), "#1b9e77", "grey70"))
  # Size: seeds are bigger (6), others smaller (3):
  V(g)$size <- ifelse(V(g)$is_seed, 6, 3)
  # Labels: show gene name for seeds OR nodes with ≥2 seed connections; hide others (NA):
  V(g)$label <- ifelse(V(g)$is_seed | V(g)$core_links >= 2,
                       V(g)$gene_name, NA)

  # set.seed(42) makes the random layout reproducible — every run gives the same plot.
  # (Force-directed layouts use randomness; the seed fixes it.)
  set.seed(42)  # Reproducible layout

  # png() opens a PNG graphics device — everything plotted goes to this file.
  # file.path() joins path components (cross-platform). FIGURE_DIR is from config.
  png(safe_filepath(FIGURE_DIR, paste0("string_network_", experiment_name), ".png"),
      width = 12, height = 10, units = "in", res = FIG_DPI)

  # plot() draws the network graph.
  # layout_with_fr() is the Fruchterman-Reingold force-directed layout algorithm:
  #   It simulates physical forces — nodes repel each other (like magnets),
  #   edges act as springs pulling connected nodes together.
  #   The result spreads nodes out naturally, clustering connected proteins.
  plot(g,
       layout = layout_with_fr(g),
       vertex.frame.color = NA,     # No border around nodes
       edge.color = "grey80",       # Light grey interaction lines
       edge.width = 0.5,
       vertex.label.cex = 0.6,      # Label text size (cex = character expansion)
       vertex.label.color = "black",
       main = paste("STRING Network:", experiment_name))

  # dev.off() closes the PNG device and saves the file.
  dev.off()
  cat(sprintf("    Saved: string_network_%s.png\n", experiment_name))

  # ---- STEP 12: Volcano plot with STRING network overlay ----
  # A volcano plot shows log2FC (x-axis) vs -log10(padj) (y-axis) for all proteins.
  # We color-code proteins that are in the STRING network.
  toPlot <- df

  # First, classify every protein as significant ("TRUE") or not ("FALSE"):
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )
  # Significant but NOT in network stays "TRUE":
  toPlot$category[!toPlot$inNetwork & toPlot$category == "TRUE"] <- "TRUE"
  # In network → "inNetwork" category:
  toPlot$category[toPlot$inNetwork] <- "inNetwork"
  # High-confidence seed → "high" category:
  toPlot$category[high_mask] <- "high"

  # factor() converts category to an ordered categorical variable.
  # levels = names(CATEGORY_COLORS) ensures consistent ordering & coloring.
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  # Select which proteins to label (high-confidence + network members):
  label_data <- toPlot[toPlot$category %in% c("high", "inNetwork"), ]

  # Build the ggplot (layered grammar of graphics):
  # ggplot2::ggplot() starts the plot; aes() maps data columns to visual properties.
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.2) +                    # Scatter points (alpha = transparency)
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +  # Custom colors
    ggrepel::geom_text_repel(                                          # Smart labels that avoid overlapping
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, max.overlaps = 20, show.legend = FALSE
    ) +
    # Dashed reference lines at the significance thresholds:
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    # Axis labels and title:
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = paste("STRING Network Overlay:", experiment_name)
    ) +
    ggplot2::theme_bw() +                                              # Clean black-and-white theme
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # save_figure() is a custom helper (from utils.R) that saves the plot.
  save_figure(p, paste0("volcano_string_", experiment_name), width = 8, height = 6)

  # Return a list of results from this function (available to the caller).
  # This includes the candidate table, network gene list, and the graph object.
  return(list(
    candidates = candidates,
    in_network = network_genes,
    graph = g
  ))
}

# =====================================================================
# MAIN: Run STRING analysis for key experiments
# =====================================================================
# Define which experiments to analyze (these names match keys in the experiments list):
key_experiments <- c("turbo_trip4_vs_wt", "flag_cflag_vs_ctrl",
                     "turbo_RA_vs_wt", "flag_RA_cflag_vs_cflag")

# Loop through each key experiment and run the STRING network analysis.
# for (exp_name in key_experiments) iterates over the vector.
for (exp_name in key_experiments) {
  # Check if this experiment exists in our data:
  if (exp_name %in% names(experiments)) {
    # Run the analysis function defined above.
    # experiments[[exp_name]] gets the data frame for this experiment.
    result <- analyze_string_network(experiments[[exp_name]], exp_name)

    # If results were produced (not NULL), save the candidate and network tables:
    if (!is.null(result)) {
      # save_table() is a custom helper (from utils.R) that writes a CSV/TSV file.
      save_table(result$candidates, paste0("string_candidates_", exp_name))
      save_table(data.frame(gene = result$in_network, category = "in_network"),
                 paste0("string_innetwork_", exp_name))
    }
  } else {
    # Print a skip message if the experiment doesn't exist:
    cat(sprintf("\n  [%s] Skipping (not found in data)\n", exp_name))
  }
}

# Print final summary with output locations:
cat("\n=========================================\n")
cat(" STRING network analysis complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
