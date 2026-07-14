###############################################################################
# 26_network_go_comparison.R
# GO enrichment comparison: "in network" vs "not in network" enriched proteins.
#
# This replicates Lydia's Panel B — two GO enrichment subplots showing
# how network-connected vs unconnected enriched proteins differ in
# their biological functions.
#
# Uses the SAME classification as R/21_lydia_network_volcano.R:
#   - Seeds = highly enriched proteins (Lydia's thresholds)
#   - Network = STRING physical neighbors of seeds (score > 250)
#   - Foreground = log2FC >= 1 AND padj <= 0.1 (TRIP4-enriched only)
#   - Universe = all proteins detected in this experiment
#
# Output:
#   - Side-by-side GO dotplots for BP, CC, MF ontologies
#   - CSV tables for each gene set
#
# Usage:
#   make network-go
###############################################################################

cat("\n=========================================\n")
cat(" Network GO Comparison (In-Net vs Not-In-Net)\n")
cat("=========================================\n\n")

library(STRINGdb)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)
source("R/00_theme.R")

# ---- Load data ----
experiments <- load_all_experiments()

turbo_exp_name <- "BK467_TRIP4_vs_BK467_WT"
if (!turbo_exp_name %in% names(experiments)) {
  cat("ERROR: Main experiment not found:", turbo_exp_name, "\n")
  quit(status = 1)
}

df <- experiments[[turbo_exp_name]]
df$neglog10p <- -log10(df$padj)
cat(sprintf("  Loaded TurboID TRIP4 vs WT: %d proteins\n", nrow(df)))

# ---- Load known interactors + ASCC ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)
ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")

# ---- Load STRING physical interactions ----
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))
} else {
  cat("ERROR: Physical links file not found at:", phys_file, "\n")
  quit(status = 1)
}

# ---- Lydia's helper functions ----
get_phys_interactions <- function(string_ids, phys_edges = phys) {
  phys_edges %>% filter(from %in% string_ids | to %in% string_ids)
}

get_phys_neighbors <- function(string_ids, phys_edges = phys) {
  partners <- phys_edges %>%
    filter(from %in% string_ids | to %in% string_ids) %>%
    transmute(partner = if_else(from %in% string_ids, to, from))
  unique(partners$partner)
}

# ---- STRINGdb for ID mapping ----
string_cache <- file.path(OUTPUT_DIR, "string_cache")
dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)

string_db <- STRINGdb$new(
  version = "12.0",
  species = STRING_TAXON,
  score_threshold = 400,
  input_directory = string_cache
)

df_for_map <- as.data.frame(df[, c("gene", "log2FC", "padj", "neglog10p")])
mapped <- string_db$map(df_for_map, "gene", removeUnmappedRows = TRUE)
cat(sprintf("  Mapped %d of %d proteins to STRING\n", nrow(mapped), nrow(df)))

# =====================================================================
# Build network classification (SAME as R/21)
# =====================================================================
FC_THRESHOLD <- 1
NEGLOG10P_THRESHOLD <- 1

# Seeds (Lydia's stringent thresholds)
seed_mask <- (!is.na(df$log2FC) & !is.na(df$neglog10p)) &
  ((df$log2FC > 7.5 & df$neglog10p > 3) |
   (df$log2FC > 2 & df$neglog10p > 6))
sig1_genes <- df$gene[seed_mask]
cat(sprintf("\n  Highly enriched seeds: %d proteins\n", length(sig1_genes)))

# Network expansion
core_mapped <- mapped %>% dplyr::filter(gene %in% sig1_genes)
seed_ids <- unique(core_mapped$STRING_id)
cat(sprintf("  Seeds mapped to STRING: %d\n", length(seed_ids)))

if (length(seed_ids) > 0) {
  phys_neighbors_seed <- get_phys_neighbors(seed_ids)
  expanded_ids <- unique(c(seed_ids, phys_neighbors_seed))
  int_expanded <- get_phys_interactions(expanded_ids)

  int_expanded2 <- int_expanded[
    (int_expanded$from %in% seed_ids | int_expanded$to %in% seed_ids) &
    int_expanded$combined_score > 250, ]

  nearby_strong <- unique(c(
    int_expanded2$from[int_expanded2$to %in% seed_ids & int_expanded2$combined_score > 700],
    int_expanded2$to[int_expanded2$from %in% seed_ids & int_expanded2$combined_score > 700]
  ))

  if (length(nearby_strong) > 0) {
    int_expanded1 <- int_expanded[
      int_expanded$combined_score > 700 &
      (int_expanded$from %in% nearby_strong | int_expanded$to %in% nearby_strong), ]
  } else {
    int_expanded1 <- int_expanded[0, ]
  }

  int_expanded <- rbind(int_expanded1, int_expanded2)

  a <- mapped[
    mapped$STRING_id %in% unique(c(int_expanded$from, int_expanded$to)) &
    !(mapped$STRING_id %in% seed_ids), ]
  cat(sprintf("  Network proteins (excl. seeds): %d\n", nrow(a)))

} else {
  cat("  WARNING: No seeds mapped. Skipping network.\n")
  a <- mapped[0, ]
}

# =====================================================================
# Split enriched proteins into two sets
# =====================================================================
# "In network" = enriched AND (in STRING network OR is a seed OR is known interactor)
# "Not in network" = enriched AND NOT in any of the above

network_genes <- unique(c(a$gene, sig1_genes))

# All enriched proteins (foreground)
enriched_mask <- df$padj <= 0.1 & df$log2FC >= FC_THRESHOLD & !is.na(df$gene)
enriched_genes <- df$gene[enriched_mask]

# Add known interactors and ASCC to the enriched set if they pass threshold
enriched_genes <- unique(c(enriched_genes,
                           df$gene[df$gene %in% known_interactors & df$log2FC >= FC_THRESHOLD],
                           df$gene[df$gene %in% ASCC_CORE & df$log2FC >= FC_THRESHOLD]))

# Split
in_network_genes <- enriched_genes[enriched_genes %in% network_genes]
not_in_network_genes <- enriched_genes[!enriched_genes %in% network_genes]

cat(sprintf("\n  Enriched proteins total: %d\n", length(enriched_genes)))
cat(sprintf("  In network: %d\n", length(in_network_genes)))
cat(sprintf("  Not in network: %d\n", length(not_in_network_genes)))

# Universe = all detected proteins in this experiment
universe_genes <- unique(df$gene[!is.na(df$gene)])
cat(sprintf("  Universe (all detected): %d\n", length(universe_genes)))

# =====================================================================
# Export gene lists
# =====================================================================
write.table(data.frame(gene = in_network_genes),
            file = file.path(TABLE_DIR, "network_go_in_network_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(data.frame(gene = not_in_network_genes),
            file = file.path(TABLE_DIR, "network_go_not_in_network_genes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# =====================================================================
# Run GO enrichment on each set
# =====================================================================
run_go_for_split <- function(genes, set_name, universe) {
  cat(sprintf("\n  Running GO for '%s' (%d genes)...\n", set_name, length(genes)))

  if (length(genes) < 5) {
    cat("    Too few genes (<5). Skipping.\n")
    return(list())
  }

  results <- list()
  for (ont in c("BP", "CC", "MF")) {
    cat(sprintf("    [%s] enrichGO... ", ont))

    ego <- tryCatch({
      enrichGO(
        gene          = genes,
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        minGSSize     = 2,
        maxGSSize     = 5000
      )
    }, error = function(e) {
      cat("ERROR: %s\n", conditionMessage(e))
      return(NULL)
    })

    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      cat("No enriched terms\n")
      results[[ont]] <- NULL
      next
    }

    # Simplify
    ego_simplified <- tryCatch({
      simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }, error = function(e) ego)

    res <- as.data.frame(ego_simplified)
    cat(sprintf("%d terms (after simplify)\n", nrow(res)))

    # Save table
    save_table(res, sprintf("network_go_%s_%s", set_name, ont))

    results[[ont]] <- ego_simplified
  }

  return(results)
}

cat("\n--- GO Enrichment: In Network ---")
go_in_net <- run_go_for_split(in_network_genes, "in_network", universe_genes)

cat("\n--- GO Enrichment: Not In Network ---")
go_not_net <- run_go_for_split(not_in_network_genes, "not_in_network", universe_genes)

# =====================================================================
# Plot In-Network GO enrichment (single plot — Not-In-Network excluded)
# =====================================================================
cat("\n--- Generating In-Network GO plots ---\n")

ONT_NAMES <- c("BP" = "Biological Process",
               "CC" = "Cellular Component",
               "MF" = "Molecular Function")

for (ont in c("BP", "CC", "MF")) {
  has_in <- !is.null(go_in_net[[ont]]) && nrow(as.data.frame(go_in_net[[ont]])) > 0

  if (!has_in) {
    cat(sprintf("  [%s] No enriched terms. Skipping.\n", ont))
    next
  }

  cat(sprintf("  [%s] Building plot... ", ont))

  n_show <- min(15, nrow(as.data.frame(go_in_net[[ont]])))
  p_in <- dotplot(go_in_net[[ont]], showCategory = n_show)
  cnt <- p_in$data$Count
  p_in <- p_in +
    ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                   name = "p-adjusted value") +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 10),
                                   breaks = make_size_breaks(cnt, n_breaks = 8),
                                   limits = c(min(cnt), max(cnt))) +
    ggplot2::guides(size = size_legend_guide()) +
    ggplot2::scale_y_discrete(labels = wrap_go_labels) +
    ggplot2::labs(x = "Gene Ratio") +
    theme_poster() +
    ggplot2::theme(axis.text =
      ggplot2::element_text(size = 22, color = "black"),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.text  = ggplot2::element_text(size = 18))

  filename <- sprintf("network_go_comparison_%s", ont)
  save_figure(p_in, filename, width = 18, height = max(14, n_show * 0.8))

  cat("done\n")
}

# =====================================================================
# KEGG pathway enrichment for in-network and not-in-network sets
# =====================================================================
cat("\n--- KEGG Pathway Enrichment ---\n")

run_kegg_for_split <- function(genes, set_name, universe) {
  cat(sprintf("\n  Running KEGG for '%s' (%d genes)...\n",
              set_name, length(genes)))

  if (length(genes) < 5) {
    cat("    Too few genes (<5). Skipping KEGG.\n")
    return(NULL)
  }

  tryCatch({
    # Convert SYMBOL to ENTREZID for KEGG
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    entrez_genes <- unique(entrez_map$ENTREZID)
    cat(sprintf("    Mapped %d/%d genes to ENTREZID\n",
                length(entrez_genes), length(genes)))

    # Convert universe SYMBOL to ENTREZID
    universe_map <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
    universe_entrez <- unique(universe_map$ENTREZID)

    ekegg <- enrichKEGG(
      gene          = entrez_genes,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      universe      = universe_entrez,
      minGSSize     = 2,
      maxGSSize     = 5000
    )

    if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
      cat("    No enriched KEGG pathways\n")
      return(NULL)
    }

    # Convert ENTREZID back to readable gene symbols
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

    kegg_df <- as.data.frame(ekegg)
    cat(sprintf("    Found %d enriched KEGG pathways\n", nrow(kegg_df)))

    save_table(kegg_df, sprintf("network_KEGG_%s", set_name))

    # Dotplot
    n_show <- min(20, nrow(kegg_df))
    fig_height <- max(7, n_show * 0.4)

    p <- enrichplot::dotplot(ekegg, showCategory = n_show)
    cnt_k <- p$data$Count
    p <- p +
      ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                     name = "p-adjusted value") +
      ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 10),
                                     breaks = make_size_breaks(cnt_k, n_breaks = 8),
                                     limits = c(min(cnt_k), max(cnt_k))) +
      ggplot2::scale_y_discrete(labels = wrap_go_labels) +
      ggplot2::guides(size = size_legend_guide()) +
      ggplot2::labs(x = "Gene Ratio") +
      theme_poster() +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 18, face = "bold"),
                     legend.text  = ggplot2::element_text(size = 18))
    save_figure(p, sprintf("network_KEGG_%s_dotplot", set_name),
                width = 18, height = max(14, n_show * 0.8))

  }, error = function(e) {
    cat(sprintf("    KEGG error (may need internet): %s\n",
                conditionMessage(e)))
  })
}

# KEGG for in-network set
run_kegg_for_split(in_network_genes, "in_network", universe_genes)

# KEGG for not-in-network set
run_kegg_for_split(not_in_network_genes, "not_in_network", universe_genes)

cat("\n=========================================\n")
cat(" Network GO Comparison complete!\n")
cat("=========================================\n")
cat(sprintf("  In network: %d genes\n", length(in_network_genes)))
cat(sprintf("  Not in network: %d genes\n", length(not_in_network_genes)))
cat("  Output: output/figures/network_go_comparison_*.png\n")
cat("=========================================\n")
