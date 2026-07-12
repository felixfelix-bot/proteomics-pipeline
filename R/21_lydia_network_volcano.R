###############################################################################
# 21_lydia_network_volcano.R
# Lydia-style volcano: STRING physical interaction network overlay.
#
# WHAT THIS DOES (matching Lydia's actual code exactly):
#   1. Load TurboID TRIP4 vs WT data + known interactors
#   2. Map proteins to STRINGdb
#   3. Define SEEDS: high-confidence hits (Lydia's stringent thresholds)
#   4. Find physical neighbors of seeds
#   5. Mark each protein: in the seed-neighbor network or not
#   6. Plot volcano colored by network membership:
#        "high"  = seed protein (stringent threshold)  → red
#        TRUE    = in STRING physical network           → teal #1b9e77
#        FALSE   = not in network                       → grey60
#   7. Label known interactors + ASCC complex only
#
# This is SIMPLE on purpose. Lydia's plot has 3 colors, not 9.
# Gene families (DHX/DDX/etc) have their own script (R/06_gene_families.R).
#
# Usage:
#   make lydia-volcano
###############################################################################

library(STRINGdb)
library(igraph)
library(dplyr)

cat("\n=========================================\n")
cat(" Lydia-Style Network Volcano\n")
cat("=========================================\n\n")

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

# ---- Load known interactors ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
known_ia_excl_core <- setdiff(known_interactors, ASCC_CORE)

# ---- Define high-confidence SEED proteins (Lydia's stringent thresholds) ----
# Seed = (log2FC > 2 AND -log10(padj) > 6) OR (log2FC > 7 AND -log10(padj) > 2)
seed_mask <- (!is.na(df$log2FC) & !is.na(df$padj)) &
  ((df$log2FC > 2 & df$neglog10p > 6) |
   (df$log2FC > 7 & df$neglog10p > 2))

seed_genes <- df$gene[seed_mask]
cat(sprintf("  High-confidence seeds: %d proteins\n", length(seed_genes)))
if (length(seed_genes) > 0) {
  cat(sprintf("    Sample seeds: %s\n", paste(head(seed_genes, 15), collapse = ", ")))
}

# ---- STRING network integration ----
cat("\n--- STRING Network Analysis ---\n")

# Use local physical links file if available (avoids ~100MB download on Windows)
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

string_cache <- file.path(OUTPUT_DIR, "string_cache")
dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))
} else {
  cat("No local physical links file. Interactions will come from STRINGdb online.\n")
}

# Still need STRINGdb for gene name → STRING_id mapping
string_db <- STRINGdb$new(
  version = "12.0",
  species = STRING_TAXON,
  score_threshold = 400,
  input_directory = string_cache
)

# Map our proteins to STRING
df_for_map <- as.data.frame(df[, c("gene", "log2FC", "padj", "neglog10p")])
mapped <- string_db$map(df_for_map, "gene", removeUnmappedRows = TRUE)
cat(sprintf("  Mapped %d of %d proteins to STRING\n", nrow(mapped), nrow(df)))

# Get seed STRING IDs
core_mapped <- mapped %>% dplyr::filter(gene %in% seed_genes)
seed_ids <- unique(core_mapped$STRING_id)
cat(sprintf("  Seeds mapped to STRING: %d\n", length(seed_ids)))

# ---- Lydia's neighbor expansion ----
if (length(seed_ids) > 0) {
  cat("  Expanding network via physical interactions...\n")

  # Use local phys file if loaded
  if (exists("phys")) {
    all_edges <- phys
  } else {
    all_edges <- string_db$get_interactions()
    if (!"combined_score" %in% names(all_edges)) {
      sc <- grep("score", names(all_edges), ignore.case = TRUE, value = TRUE)
      if (length(sc) > 0) names(all_edges)[names(all_edges) == sc[1]] <- "combined_score"
    }
  }

  # Direct interactions with seeds (score > 250)
  direct_edges <- all_edges[
    (all_edges$from %in% seed_ids | all_edges$to %in% seed_ids) &
    all_edges$combined_score > 250, ]

  # Partners directly connected to seeds
  direct_partners <- unique(c(
    direct_edges$from[!(direct_edges$from %in% seed_ids)],
    direct_edges$to[!(direct_edges$to %in% seed_ids)]
  ))

  # Strong neighbors (score > 700 with seed)
  strong_edges <- direct_edges[direct_edges$combined_score > 700, ]
  nearby_strong <- unique(c(
    strong_edges$from[strong_edges$to %in% seed_ids],
    strong_edges$to[strong_edges$from %in% seed_ids]
  ))

  # Secondary interactions among strong neighbors (score > 700)
  if (length(nearby_strong) > 0) {
    secondary_edges <- all_edges[
      all_edges$combined_score > 700 &
      (all_edges$from %in% nearby_strong | all_edges$to %in% nearby_strong), ]
  } else {
    secondary_edges <- direct_edges[0, ]
  }

  all_net_edges <- unique(rbind(direct_edges, secondary_edges))
  all_net_ids <- unique(c(all_net_edges$from, all_net_edges$to))
  network_partner_ids <- setdiff(all_net_ids, seed_ids)

  # Map STRING IDs back to gene names
  network_genes <- mapped %>%
    dplyr::filter(STRING_id %in% network_partner_ids) %>%
    dplyr::pull(gene) %>%
    unique()

  cat(sprintf("  Network neighbors: %d proteins\n", length(network_genes)))

  # ---- Mark inNetwork: Lydia's exact approach ----
  # Values: "high" for seeds, TRUE for network neighbors, FALSE for everything else
  df$inNetwork <- "FALSE"
  df$inNetwork[df$gene %in% network_genes] <- TRUE
  df$inNetwork[seed_mask] <- "high"

  cat(sprintf("  In network (total): %d\n", sum(df$inNetwork != "FALSE")))

} else {
  cat("  WARNING: No seeds mapped to STRING. Skipping network analysis.\n")
  df$inNetwork <- "FALSE"
}

# ---- THE VOLCANO (Lydia's exact style) ----
# Three colors only: high=red, TRUE=teal, FALSE=grey
cat("\n--- Generating plot ---\n")

toPlot <- df

# Color scale matching Lydia's code exactly
LYDIA_NET_COLORS <- c(
  "high"  = "red",         # Lydia's exact color for seeds
  "TRUE"  = "#1b9e77",     # Lydia's exact teal for network members
  "FALSE" = "grey60"       # Lydia's exact grey for everything else
)

LYDIA_NET_LABELS <- c(
  "high"  = "High-confidence seed",
  "TRUE"  = "In STRING network",
  "FALSE" = "Not in network"
)

# Label known interactors + ASCC that are TRIP4-enriched (positive log2FC)
label_data <- toPlot[toPlot$gene %in% c(known_interactors, ASCC_CORE) &
                     !is.na(toPlot$log2FC) & toPlot$log2FC > 0, ]

p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                           color = inNetwork)) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::scale_color_manual(
    values = LYDIA_NET_COLORS,
    labels = LYDIA_NET_LABELS,
    name = NULL, drop = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 2.5, fontface = "bold",
    hjust = 0, nudge_x = 0.1,
    max.overlaps = 30, show.legend = FALSE
  ) +
  ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID vs WT — STRING Physical Network"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    axis.title.x = ggplot2::element_text(hjust = 0.5),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")
  )

save_figure(p, "lydia_network_volcano")

# ---- Export gene lists ----
cat("\n--- Exporting gene lists ---\n")

sig_genes <- df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF &
                     !is.na(df$gene)]
write.table(data.frame(gene = sig_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_all.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

net_genes <- df$gene[df$inNetwork != "FALSE" & !is.na(df$gene)]
write.table(data.frame(gene = net_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_network.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat(sprintf("  Exported: enriched_all (%d), network (%d)\n",
            length(sig_genes), length(net_genes)))

cat("\nDone. Output saved to output/figures/.\n")
