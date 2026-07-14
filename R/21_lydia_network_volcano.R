###############################################################################
# 21_lydia_network_volcano.R
# Lydia's STRING physical interaction network overlay on TRIP4 vs WT volcano.
#
# Per Aruna's July 12 voice message:
#   - Everything with log2FC < 1 = GREEN (background, ignored for analysis)
#   - Only TRIP4-enriched proteins (log2FC >= 1) get categorized
#   - Each category gets its own DISTINCT color:
#       1. ASC complex (TRIP4, ASCC1-3)
#       2. Verified interaction & in network
#       3. Verified interaction & not in network
#       4. Highly enriched (stringent seeds)
#       5. In network of highly enriched proteins
#       6. Not assigned to interaction network
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

# ---- Load known interactors ----
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

# ---- Lydia's helper functions (verbatim) ----
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
# ARUNA'S THRESHOLDS (from voice message):
#   Enriched = log2FC >= 1 AND -log10(padj) >= 1 (i.e., padj <= 0.1)
#   Background = everything else = GREEN
# =====================================================================
FC_THRESHOLD <- 0.999  # log2FC >= 1 (fold change >= 2x)
                         # Using 0.999 for floating point tolerance
                         # (R may store "1.0" from CSV as 0.9999999)
                         # Threshold line on plot still drawn at 1.0
NEGLOG10P_THRESHOLD <- 1  # -log10(padj) >= 1 (padj <= 0.1)

# =====================================================================
# STEP 1: Build sig column
# =====================================================================
# Priority order (highest wins):
#   ascc > ia > high > TRUE > FALSE
#
# Only applies to TRIP4-enriched proteins (log2FC >= 1).
# Everything with log2FC < 1 stays FALSE = green background.

df$sig <- FALSE  # default: background (green)

# Highly enriched seeds (Lydia's stringent threshold)
df$sig[!is.na(df$log2FC) & !is.na(df$neglog10p) &
       df$log2FC >= FC_THRESHOLD &
       ((df$log2FC > 7.5 & df$neglog10p > 3) |
        (df$log2FC > 2 & df$neglog10p > 6))] <- "high"

# Significant AND TRIP4-enriched (log2FC >= 1, padj <= 0.1)
sig_idx <- df$padj <= 0.1 & df$log2FC >= FC_THRESHOLD & df$sig == FALSE
df$sig[sig_idx] <- TRUE

# Known interactors (override — only if TRIP4-enriched)
df$sig[df$gene %in% known_interactors &
       df$log2FC >= FC_THRESHOLD &
       df$sig %in% c("TRUE", "high")] <- "ia"

# ASCC complex (highest priority)
df$sig[df$gene %in% ASCC_CORE &
       df$log2FC >= FC_THRESHOLD] <- "ascc"

# Seeds for network expansion (same as sig == "high")
sig1_genes <- df$gene[df$sig == "high"]
cat(sprintf("\n  Highly enriched seeds: %d proteins\n", length(sig1_genes)))

# =====================================================================
# STEP 2: STRING network expansion (Lydia's exact code)
# =====================================================================
core_mapped <- mapped %>% dplyr::filter(gene %in% sig1_genes)
seed_ids <- unique(core_mapped$STRING_id)
cat(sprintf("  Seeds mapped to STRING: %d\n", length(seed_ids)))

if (length(seed_ids) > 0) {
  cat("  Expanding network...\n")

  # 1. Get ALL physical neighbors of seeds
  phys_neighbors_seed <- get_phys_neighbors(seed_ids)
  expanded_ids <- unique(c(seed_ids, phys_neighbors_seed))

  # 2. Get ALL interactions among expanded set
  int_expanded <- get_phys_interactions(expanded_ids)

  # 3. Direct interactions with seeds (score > 250)
  int_expanded2 <- int_expanded[
    (int_expanded$from %in% seed_ids | int_expanded$to %in% seed_ids) &
    int_expanded$combined_score > 250, ]

  # 4. Strong neighbors (score > 700 with seed)
  nearby_strong <- unique(c(
    int_expanded2$from[int_expanded2$to %in% seed_ids & int_expanded2$combined_score > 700],
    int_expanded2$to[int_expanded2$from %in% seed_ids & int_expanded2$combined_score > 700]
  ))

  # 5. Secondary interactions among strong neighbors (score > 700)
  if (length(nearby_strong) > 0) {
    int_expanded1 <- int_expanded[
      int_expanded$combined_score > 700 &
      (int_expanded$from %in% nearby_strong | int_expanded$to %in% nearby_strong), ]
  } else {
    int_expanded1 <- int_expanded[0, ]
  }

  # 6. Final network
  int_expanded <- rbind(int_expanded1, int_expanded2)

  # Network proteins (excluding seeds)
  a <- mapped[
    mapped$STRING_id %in% unique(c(int_expanded$from, int_expanded$to)) &
    !(mapped$STRING_id %in% seed_ids), ]
  cat(sprintf("  Network proteins (excl. seeds): %d\n", nrow(a)))

} else {
  cat("  WARNING: No seeds mapped. Skipping network.\n")
  a <- mapped[0, ]
}

# =====================================================================
# STEP 3: Build inNetwork column
# =====================================================================
# Only TRIP4-enriched proteins (log2FC >= 1) get network membership.
# Seeds = "high", network neighbors = TRUE, everything else = FALSE.
df$inNetwork <- FALSE

# Network neighbors (proteins in STRING network of seeds)
df$inNetwork[df$gene %in% a$gene &
             !is.na(df$gene) &
             df$log2FC >= FC_THRESHOLD] <- TRUE

# Seeds get "high" inNetwork status
df$inNetwork[df$sig == "high"] <- "high"

# =====================================================================
# STEP 4: Create sig_network = paste(sig, inNetwork) — THE COMBINED COLUMN
# =====================================================================
df$sig_network <- paste(df$sig, df$inNetwork, sep = "_")

cat("\n  sig_network categories:\n")
for (cat_name in sort(unique(df$sig_network))) {
  count <- sum(df$sig_network == cat_name, na.rm = TRUE)
  cat(sprintf("    %-20s: %d proteins\n", cat_name, count))
}

# =====================================================================
# STEP 5: Map sig_network values to DISTINCT colors + labels
# =====================================================================
# Per Aruna: each enriched category gets its own distinct color.
# Background (log2FC < 1) = green.

SIG_NET_COLORS <- c(
  # ASCC complex (only if TRIP4-enriched) — BLUE per Aruna
  "ascc_FALSE"   = "#0072B2",
  "ascc_TRUE"    = "#0072B2",
  "ascc_high"    = "#0072B2",

  # Known interactors — GOLD (not in net) / teal (in net)
  # Changed from light purple — too close to dark purple
  "ia_FALSE"     = "#FFD92F",
  "ia_TRUE"      = "#1B9E77",
  "ia_high"      = "#1B9E77",

  # Highly enriched seeds — red
  "high_high"    = "#E41A1C",
  "high_TRUE"    = "#E41A1C",
  "high_FALSE"   = "#E41A1C",

  # Significant in network — PURPLE (changed from blue, since ASCC is now blue)
  "TRUE_TRUE"    = "#984EA3",
  "TRUE_high"    = "#984EA3",

  # Significant not in network — orange
  "TRUE_FALSE"   = "#FF7F00",

  # Background (not enriched, log2FC < 1) — grey
  "FALSE_FALSE"  = "grey60",
  "FALSE_TRUE"   = "grey60",
  "FALSE_high"   = "grey60"
)

# Human-readable labels
SIG_NET_LABELS <- c(
  "ascc_FALSE"   = "ASCC complex",
  "ascc_TRUE"    = "ASCC complex",
  "ascc_high"    = "ASCC complex",

  "ia_FALSE"     = "Verified interaction & not in network",
  "ia_TRUE"      = "Verified interaction & in network",
  "ia_high"      = "Verified interaction & in network",

  "high_high"    = "Highly enriched",
  "high_TRUE"    = "Highly enriched",
  "high_FALSE"   = "Highly enriched",

  "TRUE_TRUE"    = "In network of highly enriched proteins",
  "TRUE_high"    = "In network of highly enriched proteins",

  "TRUE_FALSE"   = "Not assigned to interaction network",

  "FALSE_FALSE"  = "Not enriched",
  "FALSE_TRUE"   = "Not enriched",
  "FALSE_high"   = "Not enriched"
)

# Get all unique sig_network values present in data
all_cats <- sort(unique(df$sig_network))
present_colors <- SIG_NET_COLORS[all_cats]

# Build legend in Aruna's specified order:
# 1. ASCC complex (blue)
# 2. Verified interaction & in network (teal)
# 3. Verified interaction & not in network (gold)
# 4. Highly enriched (red)
# 5. In network of highly enriched proteins (dark purple)
# 6. Not assigned to interaction network (orange)
# 7. Not enriched (grey)
legend_order <- c(
  "ascc_FALSE", "ascc_TRUE", "ascc_high",          # 1. ASCC
  "ia_TRUE", "ia_high",                              # 2. Verified & in network
  "ia_FALSE",                                        # 3. Verified & not in network
  "high_high", "high_TRUE", "high_FALSE",            # 4. Highly enriched
  "TRUE_TRUE", "TRUE_high",                          # 5. In network of highly enriched
  "TRUE_FALSE",                                      # 6. Not assigned
  "FALSE_FALSE", "FALSE_TRUE", "FALSE_high"          # 7. Not enriched
)
# Keep only categories present in data, deduplicate by label
legend_breaks <- legend_order[legend_order %in% all_cats]
legend_labels <- SIG_NET_LABELS[legend_breaks]
dedup_mask <- !duplicated(legend_labels)
legend_breaks <- legend_breaks[dedup_mask]
legend_labels <- legend_labels[dedup_mask]

# Add annotation: % of enriched proteins in network
enriched_total <- sum(df$sig != FALSE & !is.na(df$gene) & df$log2FC >= FC_THRESHOLD, na.rm = TRUE)
enriched_in_net <- sum(df$inNetwork != FALSE & df$sig != FALSE &
                       !is.na(df$gene) & df$log2FC >= FC_THRESHOLD, na.rm = TRUE)
pct_in_net <- if (enriched_total > 0) round(100 * enriched_in_net / enriched_total) else 0

# =====================================================================
# STEP 6: The Volcano
# =====================================================================
cat("\n--- Generating plot ---\n")

# Labels: ASCC complex + verified interactors (including in-network).
# NOT highly enriched (Aruna: don't label those).
label_data <- df[
  ((df$gene %in% known_interactors & df$log2FC >= FC_THRESHOLD) |
   (df$gene %in% ASCC_CORE & df$log2FC >= FC_THRESHOLD)) &
  !is.na(df$log2FC), ]

p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p,
                                       color = sig_network)) +
  ggplot2::geom_point(alpha = 0.5, size = 2.5) +
  ggplot2::scale_color_manual(
    values = present_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    name = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 5, alpha = 1)
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 5.0, fontface = "bold",
    hjust = 0, nudge_x = 0.1,
    max.overlaps = 30, show.legend = FALSE,
    bg.color = "white", bg.r = 0.15
  ) +
  # Aruna's thresholds: log2FC = 1, -log10(padj) = 1
  ggplot2::geom_hline(yintercept = NEGLOG10P_THRESHOLD,
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-1, 1),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID vs WT — STRING Physical Network",
    caption = sprintf("%d of %d enriched proteins in network (%d%%)", enriched_in_net, enriched_total, pct_in_net)
  ) +
  theme_poster(font_size = 24)

save_figure(p, "lydia_network_volcano", width = 16, height = 12)

# ---- Export gene lists ----
cat("\n--- Exporting gene lists ---\n")

# Foreground: TRIP4-enriched proteins (log2FC >= 1 AND padj <= 0.1)
sig_genes <- df$gene[df$padj <= 0.1 & df$log2FC >= FC_THRESHOLD &
                     !is.na(df$gene)]
write.table(data.frame(gene = sig_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_all.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

net_genes <- df$gene[df$inNetwork != FALSE & !is.na(df$gene) &
                     df$log2FC >= FC_THRESHOLD]
write.table(data.frame(gene = net_genes),
            file = file.path(TABLE_DIR, "lydia_enriched_network.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

cat(sprintf("  Exported: enriched_all (%d), network (%d)\n",
            length(sig_genes), length(net_genes)))

cat("\nDone. Output saved to output/figures/.\n")
