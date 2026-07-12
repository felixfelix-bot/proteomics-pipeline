###############################################################################
# 21_lydia_network_volcano.R
# Lydia's STRING physical interaction network overlay on TRIP4 vs WT volcano.
#
# This is a faithful adaptation of Lydia's code:
#   2026_05_22_stats_cutoffs_turboID_map_to_STRING
#
# LYDIA'S EXACT APPROACH (traced from her code line by line):
#   1. Define seeds: (log2FC > 2 AND -log10(q) > 6) OR (log2FC > 7.5 AND -log10(q) > 3)
#   2. Map ALL proteins to STRINGdb
#   3. Get physical neighbors of seeds using local phys file
#   4. Build interaction network:
#      a. expanded_ids = seeds + ALL their physical neighbors
#      b. int_expanded = ALL interactions among expanded_ids
#      c. int_expanded2 = direct interactions with seeds (score > 250)
#      d. nearby_strong = partners with seed at score > 700
#      e. int_expanded1 = interactions among nearby_strong (score > 700)
#      f. Final network = int_expanded1 + int_expanded2
#   5. Mark each protein: inNetwork = TRUE (neighbor) / FALSE (not in network)
#   6. Overwrite seeds: inNetwork = "high"
#   7. Plot: color by inNetwork (3 colors only)
#      "high" = red, TRUE = #1b9e77 (teal), FALSE = grey60
#   8. Labels: known interactors + high-confidence seeds
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

# ---- Load STRING physical interactions ----
# Lydia loads this once at the top and uses it for all queries
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))
} else {
  cat("ERROR: Physical links file not found at:", phys_file, "\n")
  cat("This file is required. Place 9606.protein.physical.links.v12.0.txt in data/\n")
  quit(status = 1)
}

# ---- Lydia's helper functions (exact copy) ----
get_phys_interactions <- function(string_ids, phys_edges = phys) {
  phys_edges %>%
    filter(from %in% string_ids | to %in% string_ids)
}

get_phys_neighbors <- function(string_ids, phys_edges = phys) {
  partners <- phys_edges %>%
    filter(from %in% string_ids | to %in% string_ids) %>%
    transmute(partner = if_else(from %in% string_ids, to, from))
  unique(partner$partner)
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

# Map our proteins to STRING (as.data.frame prevents dimension error)
df_for_map <- as.data.frame(df[, c("gene", "log2FC", "padj", "neglog10p")])
mapped <- string_db$map(df_for_map, "gene", removeUnmappedRows = TRUE)
cat(sprintf("  Mapped %d of %d proteins to STRING\n", nrow(mapped), nrow(df)))

# =====================================================================
# STEP 1: Define seeds (Lydia's EXACT thresholds for entry 1 = TRIP4 vs WT)
# =====================================================================
# Lydia's code:
#   sig1 <- Gene.name[!is.na(q) & ((log2FC>2 & -log10(q)>6) | (log2FC>7.5 & -log10(q)>3))]
#
# Note: 7.5 and 3, NOT 7 and 2. This is the WT comparison (most stringent).

seed_mask <- (!is.na(df$log2FC) & !is.na(df$padj)) &
  ((df$log2FC > 2 & df$neglog10p > 6) |
   (df$log2FC > 7.5 & df$neglog10p > 3))

sig1 <- df$gene[seed_mask]
cat(sprintf("\n  Seeds (Lydia's threshold): %d proteins\n", length(sig1)))
if (length(sig1) > 0) {
  cat(sprintf("    Sample: %s\n", paste(head(sig1, 15), collapse = ", ")))
}

# =====================================================================
# STEP 2: Network expansion (Lydia's EXACT code, traced line by line)
# =====================================================================

# core_genes = sig1 (the seeds)
# core_mapped = seeds mapped to STRING IDs
core_mapped <- mapped %>% dplyr::filter(gene %in% sig1)
seed_ids <- unique(core_mapped$STRING_id)
cat(sprintf("  Seeds mapped to STRING: %d\n", length(seed_ids)))

if (length(seed_ids) > 0) {
  cat("  Expanding network...\n")

  # 1. Get ALL physical neighbors of seeds
  phys_neighbors_seed <- get_phys_neighbors(seed_ids)
  cat(sprintf("    Physical neighbors of seeds: %d\n", length(phys_neighbors_seed)))

  # 2. Expand set: seeds + ALL their neighbors
  expanded_ids <- unique(c(seed_ids, phys_neighbors_seed))

  # 3. Get ALL interactions among seeds + neighbors
  #    (NOT the entire phys database — only interactions WITHIN the expanded set)
  int_expanded <- get_phys_interactions(expanded_ids)
  cat(sprintf("    Interactions among expanded set: %d\n", nrow(int_expanded)))

  # 4. Direct interactions with seeds (score > 250)
  #    Lydia: int_expanded2 <- int_expanded[(from or to in seed_ids) & combined_score > 250, ]
  int_expanded2 <- int_expanded[
    (int_expanded$from %in% seed_ids | int_expanded$to %in% seed_ids) &
    int_expanded$combined_score > 250, ]
  cat(sprintf("    Direct interactions with seeds (score>250): %d\n", nrow(int_expanded2)))

  # 5. Strong neighbors (interacting with seed at score > 700)
  #    Lydia: nearby_strong = partners of seeds where edge score > 700
  nearby_strong <- unique(c(
    int_expanded2$from[int_expanded2$to %in% seed_ids & int_expanded2$combined_score > 700],
    int_expanded2$to[int_expanded2$from %in% seed_ids & int_expanded2$combined_score > 700]
  ))
  cat(sprintf("    Strong neighbors (score>700 with seed): %d\n", length(nearby_strong)))

  # 6. Secondary interactions among strong neighbors (score > 700)
  #    Lydia uses int_expanded (NOT the full phys database) — only interactions
  #    WITHIN the expanded set. This is the key difference from my previous version.
  if (length(nearby_strong) > 0) {
    int_expanded1 <- int_expanded[
      int_expanded$combined_score > 700 &
      (int_expanded$from %in% nearby_strong | int_expanded$to %in% nearby_strong), ]
    cat(sprintf("    Secondary interactions (score>700 among strong): %d\n", nrow(int_expanded1)))
  } else {
    int_expanded1 <- int_expanded[0, ]
  }

  # 7. Combine: final network edges
  int_expanded <- rbind(int_expanded1, int_expanded2)
  cat(sprintf("    Combined network edges: %d\n", nrow(int_expanded)))

  # 8. Network proteins (excluding seeds themselves)
  #    Lydia: a = mapped proteins in network but NOT seeds
  a <- mapped[
    mapped$STRING_id %in% unique(c(int_expanded$from, int_expanded$to)) &
    !(mapped$STRING_id %in% seed_ids), ]
  cat(sprintf("    Network proteins (excl. seeds): %d\n", nrow(a)))

} else {
  cat("  WARNING: No seeds mapped to STRING.\n")
  a <- mapped[0, ]
}

# =====================================================================
# STEP 3: Mark inNetwork (Lydia's EXACT approach)
# =====================================================================

# First: set TRUE/FALSE for all proteins
df$inNetwork <- FALSE
df$inNetwork[df$gene %in% a$gene & !is.na(df$gene)] <- TRUE

# Then: overwrite seeds with "high"
# Lydia: toPlot$inNetwork[(log2FC>2 & -log10(q)>6) | (log2FC>7.5 & -log10(q)>3)] <- "high"
df$inNetwork[seed_mask] <- "high"

# Report
cat(sprintf("\n  inNetwork breakdown:\n"))
cat(sprintf("    'high' (seeds): %d\n", sum(df$inNetwork == "high")))
cat(sprintf("    TRUE (network neighbors): %d\n", sum(df$inNetwork == TRUE)))
cat(sprintf("    FALSE (not in network): %d\n", sum(df$inNetwork == FALSE)))

# =====================================================================
# STEP 4: The Volcano (Lydia's EXACT style)
# =====================================================================
cat("\n--- Generating plot ---\n")

toPlot <- df

# Lydia's EXACT colors
LYDIA_COLORS <- c(
  "high"  = "red",
  "TRUE"  = "#1b9e77",
  "FALSE" = "grey60"
)

LYDIA_LABELS <- c(
  "high"  = "High-confidence seed",
  "TRUE"  = "In STRING network",
  "FALSE" = "Not in network"
)

# Labels: known interactors + ASCC + seeds (per Aruna: no gene families)
# Lydia labels: sig == "ia" | sig == "gp" | sig == "dhx" | sig == "ddx" | sig == "high"
# Aruna doesn't want gene families, so: known interactors + ASCC + high seeds
label_data <- toPlot[
  (toPlot$gene %in% known_interactors | toPlot$gene %in% ASCC_CORE |
   toPlot$inNetwork == "high") &
  !is.na(toPlot$log2FC) & toPlot$log2FC > 0, ]

# Lydia's plot style: alpha=0.4, geom_text hjust=0 nudge_x=0.1 size=2
p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                           color = inNetwork)) +
  ggplot2::geom_point(alpha = 0.4, size = 1) +
  ggplot2::scale_color_manual(
    values = LYDIA_COLORS,
    labels = LYDIA_LABELS,
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
    hjust = 0,
    nudge_x = 0.1,
    max.overlaps = 30, show.legend = FALSE,
    bg.color = "white", bg.r = 0.15
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

# =====================================================================
# STEP 5: Export gene lists (Lydia's write.table calls)
# =====================================================================
cat("\n--- Exporting gene lists ---\n")

# All significant (enriched) genes
sig_genes <- df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF &
                     !is.na(df$gene)]
write.table(data.frame(gene = sig_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_all.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Network genes (in STRING interaction network, including seeds)
net_genes <- df$gene[df$inNetwork != FALSE & !is.na(df$gene)]
write.table(data.frame(gene = net_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_network.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat(sprintf("  Exported: enriched_all (%d), network (%d)\n",
            length(sig_genes), length(net_genes)))

cat("\nDone. Output saved to output/figures/.\n")
