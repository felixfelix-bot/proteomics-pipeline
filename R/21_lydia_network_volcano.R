###############################################################################
# 21_lydia_network_volcano.R
# Lydia-style volcano: STRING physical interaction network overlay.
#
# WHAT THIS DOES (adapted from Lydia's reference code):
#   1. Load TurboID TRIP4 vs WT experiment data
#   2. Load known interactors from known_interactors.txt
#   3. Map all proteins to STRINGdb (physical interactions, score >= 400)
#   4. Define SEED proteins: high-confidence hits
#      (log2FC > 2 AND -log10(padj) > 6) OR (log2FC > 7 AND -log10(padj) > 2)
#   5. Get physical neighbors of seeds from STRING
#   6. Build interaction network: seeds + direct neighbors + secondary neighbors
#      Score filtering: > 250 for direct, > 700 for secondary
#   7. Mark each protein: inNetwork = TRUE/FALSE
#   8. Generate THREE volcano views:
#      A. Known interactors highlighted + STRING network overlay
#      B. Colored by network membership only (Lydia's "volcano_wPhysIa")
#      C. Colored by combined sig_network (Lydia's "volcano_wPhysIa_wIa")
#
# VOLCANO CATEGORIES (Lydia's priority order, adapted to our palette):
#   "ascc"     = ASCC complex (TRIP4, ASCC1-3)              → Navy #0072B2
#   "ia"       = Known interactors (from database)           → Green #009E73
#   "high"     = High-confidence seed (stringent threshold)  → Red #D55E00
#   "net_true" = In STRING network (neighbor of seed)        → Teal #1b9e77
#   "sig"      = Significant but not in network              → Vermillion #D55E00
#   "nonsig"   = Not significant                             → Grey #D0D0D0
#
# LEGEND: dots are LARGE for poster visibility. Dot colors match legend colors.
#
# Usage:
#   make lydia-volcano
###############################################################################

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

# ---- Define Lydia's categories ----
# Priority order (last assignment wins, matching Lydia's code):
# nonsig → sig → wt_enriched → dhx → ddx → gpatch → larp → ia → ascc → high

df$lydia_cat <- "nonsig"

# Significant proteins
sig_idx <- df$padj < P_VALUE_CUTOFF & abs(df$log2FC) > LOG2FC_CUTOFF
df$lydia_cat[sig_idx] <- "sig"

# WT-enriched (significant, negative log2FC)
wt_idx <- df$padj < P_VALUE_CUTOFF & df$log2FC < -LOG2FC_CUTOFF
df$lydia_cat[wt_idx] <- "wt_enriched"

# Gene families (Lydia highlights these — RNA-binding helicases)
df$lydia_cat[grepl("^DHX", df$gene)] <- "dhx"
df$lydia_cat[grepl("^DDX", df$gene)] <- "ddx"
df$lydia_cat[grepl("^GPATCH", df$gene)] <- "gpatch"
df$lydia_cat[grepl("^LARP", df$gene)] <- "larp"

# Known interactors (overwrites sig/wt/families regardless of significance)
df$lydia_cat[df$gene %in% known_ia_excl_core] <- "ia"

# ASCC complex (highest priority)
df$lydia_cat[df$gene %in% ASCC_CORE] <- "ascc"

# ---- Define high-confidence SEED proteins (Lydia's stringent thresholds) ----
# Seed = (log2FC > 2 AND -log10(padj) > 6) OR (log2FC > 7 AND -log10(padj) > 2)
# These are the proteins Lydia uses as starting points for STRING neighbor expansion
seed_mask <- (!is.na(df$log2FC) & !is.na(df$padj)) &
  ((df$log2FC > 2 & df$neglog10p > 6) |
   (df$log2FC > 7 & df$neglog10p > 2))

df$lydia_cat[seed_mask] <- "high"

seed_genes <- df$gene[seed_mask]
cat(sprintf("  High-confidence seeds: %d proteins\n", length(seed_genes)))
if (length(seed_genes) > 0) {
  cat(sprintf("    Sample seeds: %s\n", paste(head(seed_genes, 15), collapse = ", ")))
}

# ---- STRING network integration ----
cat("\n--- STRING Network Analysis ---\n")

# Use local physical links file if available (provided by researcher).
# This avoids the ~100MB STRINGdb download that fails on Windows/OneDrive.
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))

  # Still need STRINGdb for ID mapping (gene name → STRING_id)
  string_cache <- file.path(OUTPUT_DIR, "string_cache")
  dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)
  string_db <- STRINGdb$new(
    version = "12.0",
    species = STRING_TAXON,
    score_threshold = 400,
    input_directory = string_cache
  )
} else {
  cat("No local physical links file found. Using STRINGdb online...\n")
  cat("(Place 9606.protein.physical.links.v12.0.txt in data/ for offline use)\n")
  string_cache <- file.path(OUTPUT_DIR, "string_cache")
  dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)
  string_db <- STRINGdb$new(
    version = "12.0",
    species = STRING_TAXON,
    score_threshold = 400,
    input_directory = string_cache
  )
}

# Map our proteins to STRING
# as.data.frame() prevents the "incorrect dimensions" error with tibbles
df_for_map <- as.data.frame(df[, c("gene", "log2FC", "padj", "neglog10p", "lydia_cat")])
mapped <- string_db$map(df_for_map, "gene", removeUnmappedRows = TRUE)
cat(sprintf("  Mapped %d of %d proteins to STRING\n", nrow(mapped), nrow(df)))

# Get seed STRING IDs
core_mapped <- mapped %>%
  dplyr::filter(gene %in% seed_genes)
seed_ids <- unique(core_mapped$STRING_id)
cat(sprintf("  Seeds mapped to STRING: %d\n", length(seed_ids)))

# Lydia's neighbor expansion approach:
# 1. Get all physical interactions touching seeds
# 2. Filter: direct interactions with seeds (score > 250)
# 3. Secondary: interactions among neighbors of seeds (score > 700)

if (length(seed_ids) > 0) {
  cat("  Expanding network via physical interactions...\n")

  # Use local physical links file if loaded, otherwise query STRINGdb
  if (exists("phys")) {
    all_interactions <- phys
    cat(sprintf("  Using local physical links: %d edges\n", nrow(all_interactions)))
  } else {
    all_interactions <- string_db$get_interactions()
    if (!"combined_score" %in% names(all_interactions)) {
      score_col <- grep("score", names(all_interactions), ignore.case = TRUE, value = TRUE)
      if (length(score_col) > 0) {
        names(all_interactions)[names(all_interactions) == score_col[1]] <- "combined_score"
      }
    }
  }

  # Helper: get interactions touching any of the given IDs
  get_interactions_for <- function(ids, edges) {
    edges[edges$from %in% ids | edges$to %in% ids, ]
  }

  # Direct interactions with seeds (score > 250)
  direct_edges <- all_interactions[
    (all_interactions$from %in% seed_ids | all_interactions$to %in% seed_ids) &
    all_interactions$combined_score > 250, ]

  # IDs directly connected to seeds (excluding seeds themselves)
  direct_partners <- unique(c(
    direct_edges$from[!(direct_edges$from %in% seed_ids)],
    direct_edges$to[!(direct_edges$to %in% seed_ids)]
  ))

  # Strong neighbors (score > 700 with a seed)
  strong_edges <- direct_edges[direct_edges$combined_score > 700, ]
  nearby_strong <- unique(c(
    strong_edges$from[strong_edges$to %in% seed_ids],
    strong_edges$to[strong_edges$from %in% seed_ids]
  ))

  # Secondary interactions among strong neighbors (score > 700)
  if (length(nearby_strong) > 0) {
    secondary_edges <- all_interactions[
      all_interactions$combined_score > 700 &
      (all_interactions$from %in% nearby_strong | all_interactions$to %in% nearby_strong), ]
  } else {
    secondary_edges <- direct_edges[0, ]  # empty
  }

  # Combined interaction network
  all_net_edges <- unique(rbind(direct_edges, secondary_edges))
  all_net_ids <- unique(c(all_net_edges$from, all_net_edges$to))

  # Mark which proteins are in the network (excluding seeds themselves)
  network_partner_ids <- setdiff(all_net_ids, seed_ids)

  # Map STRING IDs back to gene names
  network_genes <- mapped %>%
    dplyr::filter(STRING_id %in% network_partner_ids) %>%
    dplyr::pull(gene) %>%
    unique()

  cat(sprintf("  Network proteins (neighbors): %d\n", length(network_genes)))

  # Mark inNetwork on our data
  df$inNetwork <- FALSE
  df$inNetwork[df$gene %in% network_genes] <- TRUE
  df$inNetwork[df$gene %in% seed_genes] <- TRUE  # Seeds are also "in network"

  cat(sprintf("  Proteins in network (total): %d\n", sum(df$inNetwork, na.rm = TRUE)))

  # Create combined category: sig_network (like Lydia's paste(sig, inNetwork))
  df$sig_network <- paste(df$lydia_cat,
                          ifelse(df$inNetwork, "net", "nonet"),
                          sep = "_")

} else {
  cat("  WARNING: No seeds mapped to STRING. Skipping network analysis.\n")
  df$inNetwork <- FALSE
  df$sig_network <- paste(df$lydia_cat, "nonet", sep = "_")
}

# ---- PLOT A: Network membership volcano (Lydia's primary view) ----
# Colored by inNetwork: high=red, in network=teal, not in network=grey
cat("\n--- Generating plots ---\n")

toPlot <- df
# Override category for network visualization
toPlot$net_cat <- "not_in_network"
toPlot$net_cat[toPlot$inNetwork & toPlot$lydia_cat != "high"] <- "in_network"
toPlot$net_cat[toPlot$lydia_cat == "high"] <- "seed_high"

# Ensure known interactors + gene families visible even if in network
toPlot$net_cat[toPlot$lydia_cat == "ia"] <- "known_ia"
toPlot$net_cat[toPlot$lydia_cat == "ascc"] <- "ascc"
toPlot$net_cat[toPlot$lydia_cat == "dhx"] <- "dhx"
toPlot$net_cat[toPlot$lydia_cat == "ddx"] <- "ddx"
toPlot$net_cat[toPlot$lydia_cat == "gpatch"] <- "gpatch"
toPlot$net_cat[toPlot$lydia_cat == "larp"] <- "larp"

toPlot$net_cat <- factor(toPlot$net_cat,
  levels = c("ascc", "known_ia", "dhx", "ddx", "gpatch", "larp",
             "seed_high", "in_network", "not_in_network"))

NET_COLORS <- c(
  "ascc"            = "#0072B2",   # Navy — ASCC complex
  "known_ia"        = "#009E73",   # Green — known interactors
  "dhx"             = "#66a61e",   # Green-brown — DHX helicases
  "ddx"             = "#1f78b4",   # Blue — DDX helicases
  "gpatch"          = "#e6ab02",   # Gold — GPATCH family
  "larp"            = "#a6761d",   # Bronze — LARP family
  "seed_high"       = "#D55E00",   # Vermillion — high-confidence seed
  "in_network"      = "#1b9e77",   # Teal — in STRING network
  "not_in_network"  = "#D0D0D0"    # Light grey
)

NET_LABELS <- c(
  "ascc"            = "ASCC complex",
  "known_ia"        = "Known interactors",
  "dhx"             = "DHX helicases",
  "ddx"             = "DDX helicases",
  "gpatch"          = "GPATCH family",
  "larp"            = "LARP family",
  "seed_high"       = "High-confidence (STRING seed)",
  "in_network"      = "In STRING network",
  "not_in_network"  = "Not in network"
)

# Label only known interactors and ASCC (not network neighbors, per Aruna's request)
label_data <- toPlot[toPlot$lydia_cat %in% c("ia", "ascc") &
                     !is.na(toPlot$log2FC) & toPlot$log2FC > 0, ]

p_net <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                               color = net_cat)) +
  ggplot2::geom_point(alpha = 0.5, size = 1.2) +
  ggplot2::scale_color_manual(
    values = NET_COLORS,
    labels = NET_LABELS,
    name = NULL, drop = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5, alpha = 1)  # BIGGER legend dots for poster
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 3, fontface = "bold",
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
    title = "TRIP4 TurboID vs WT — STRING Physical Network Overlay"
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

save_figure(p_net, "lydia_network_volcano")

# ---- PLOT B: Combined sig_network volcano (Lydia's second view) ----
# Shows the interaction between significance and network membership
# Each category shows whether it's also in the STRING network

# Define colors for combined categories
# Only show categories that actually appear
combined_cats <- unique(toPlot$sig_network)
combined_cats <- combined_cats[order(combined_cats)]

# Build a color mapping for combined categories
# Network=True variants get saturated colors, Network=False get lighter
SIG_NET_COLORS <- c()

# Assign colors based on base category + network status
for (cc in combined_cats) {
  parts <- strsplit(cc, "_")[[1]]
  net_status <- parts[length(parts)]
  base_cat <- paste(parts[-length(parts)], collapse = "_")

  if (net_status == "net") {
    # In network — use saturated color
    base_color <- switch(base_cat,
      "ascc" = "#0072B2",
      "ia" = "#009E73",
      "high" = "#D55E00",
      "sig" = "#56B4E9",
      "wt" = "#8C8C8C",
      "nonsig" = "#B0B0B0",
      "#999999")
  } else {
    # Not in network — use lighter/grey version
    base_color <- switch(base_cat,
      "ascc" = "#7BAFD4",
      "ia" = "#7EC8A0",
      "high" = "#E8A577",
      "sig" = "#B0D8E8",
      "wt" = "#C4C4C4",
      "nonsig" = "#D8D8D8",
      "#D0D0D0")
  }
  SIG_NET_COLORS[cc] <- base_color
}

# Human-readable labels for combined categories
SIG_NET_LABELS <- c()
for (cc in combined_cats) {
  parts <- strsplit(cc, "_")[[1]]
  net_status <- parts[length(parts)]
  base_cat <- paste(parts[-length(parts)], collapse = "_")

  base_label <- switch(base_cat,
    "ascc" = "ASCC complex",
    "ia" = "Known interactors",
    "high" = "High-confidence seed",
    "sig" = "Significant",
    "wt" = "WT-enriched",
    "nonsig" = "Not significant",
    base_cat)

  net_label <- ifelse(net_status == "net", "in network", "not in network")
  SIG_NET_LABELS[cc] <- paste0(base_label, " (", net_label, ")")
}

p_combined <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = neglog10p,
                                                     color = sig_network)) +
  ggplot2::geom_point(alpha = 0.5, size = 1.2) +
  ggplot2::scale_color_manual(
    values = SIG_NET_COLORS,
    labels = SIG_NET_LABELS,
    name = NULL, drop = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 4, alpha = 1)
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 3, fontface = "bold",
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
    title = "TRIP4 TurboID vs WT — Significance + Network Combined"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    axis.title.x = ggplot2::element_text(hjust = 0.5),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 8),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")
  )

save_figure(p_combined, "lydia_combined_volcano")

# ---- Export gene lists (like Lydia's write.table calls) ----
# These are useful for ShinyGO and other tools
cat("\n--- Exporting gene lists ---\n")

# All significant (enriched) genes
sig_genes <- df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF &
                     !is.na(df$gene)]
write.table(data.frame(gene = sig_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_all.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Network genes (in STRING interaction network)
net_genes <- df$gene[df$inNetwork & !is.na(df$gene)]
write.table(data.frame(gene = net_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_network.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Significant but NOT in network
nonnet_sig <- df$gene[!df$inNetwork & df$padj < P_VALUE_CUTOFF &
                      df$log2FC > LOG2FC_CUTOFF & !is.na(df$gene)]
write.table(data.frame(gene = nonnet_sig),
            file = file.path(TABLE_DIR, "lydia_enriched_notInNetwork.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat(sprintf("  Exported: enriched_all (%d), network (%d), notInNetwork (%d)\n",
            length(sig_genes), length(net_genes), length(nonnet_sig)))

cat("\nDone. Output saved to output/figures/.\n")
