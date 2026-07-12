###############################################################################
# 28_string_style_network.R
# STRING website-style network visualization for RA common proteins.
#
# Replicates the STRING database's native "bubbles" view:
#   - Force-directed (Fruchterman-Reingold) layout
#   - Large semi-transparent colored nodes (functional grouping by cluster)
#   - Multi-colored edges by evidence/confidence
#   - Gene name labels with white outline
#
# Works for both common enriched and common depleted protein sets.
#
# Usage:
#   make string-style-network
###############################################################################

cat("\n=========================================\n")
cat(" STRING-Style Network (RA Common Proteins)\n")
cat("=========================================\n\n")

library(STRINGdb)
library(igraph)
library(ggplot2)

# ---- Load data ----
experiments <- load_all_experiments()

RA02_EXP <- "BK467_TRIP4_RA02_vs_BK467_TRIP4"
RA04_EXP <- "BK504_TRIP4_RA04_vs_BK504_TRIP4"

df_ra02 <- find_experiment(experiments, RA02_EXP)
df_ra04 <- find_experiment(experiments, RA04_EXP)

if (is.null(df_ra02) || is.null(df_ra04)) {
  cat("ERROR: RA experiments not found.\n")
  quit(status = 1)
}

# ---- Helper: extract enriched/depleted ----
get_enriched <- function(df) {
  unique(df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC >= LOG2FC_CUTOFF & !is.na(df$gene)])
}
get_depleted <- function(df) {
  unique(df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC <= -LOG2FC_CUTOFF & !is.na(df$gene)])
}

ra02_enr <- get_enriched(df_ra02)
ra04_enr <- get_enriched(df_ra04)
ra02_dep <- get_depleted(df_ra02)
ra04_dep <- get_depleted(df_ra04)

common_enriched <- intersect(ra02_enr, ra04_enr)
common_depleted <- intersect(ra02_dep, ra04_dep)

cat(sprintf("  Common RA-enriched: %d proteins\n", length(common_enriched)))
cat(sprintf("  Common RA-depleted: %d proteins\n", length(common_depleted)))

# ---- Build summary tables with mean log2FC ----
build_table <- function(genes) {
  if (length(genes) == 0) return(data.frame())
  sub_a <- df_ra02[df_ra02$gene %in% genes, c("gene", "log2FC")]
  sub_b <- df_ra04[df_ra04$gene %in% genes, c("gene", "log2FC")]
  sub_a <- sub_a[!duplicated(sub_a$gene), ]
  sub_b <- sub_b[!duplicated(sub_b$gene), ]
  colnames(sub_a) <- c("gene", "log2FC_RA02")
  colnames(sub_b) <- c("gene", "log2FC_RA04")
  merged <- merge(sub_a, sub_b, by = "gene", all = TRUE)
  merged$mean_log2FC <- rowMeans(merged[, c("log2FC_RA02", "log2FC_RA04")], na.rm = TRUE)
  merged
}

enriched_table <- build_table(common_enriched)
depleted_table <- build_table(common_depleted)

# ---- Load STRING physical interactions ----
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")
if (!file.exists(phys_file)) {
  cat("ERROR: STRING physical links file not found.\n")
  quit(status = 1)
}
cat("Loading STRING physical interactions...\n")
phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
names(phys) <- c("from", "to", "combined_score")

# ---- STRINGdb for ID mapping ----
string_cache <- file.path(OUTPUT_DIR, "string_cache")
dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)
string_db <- STRINGdb$new(version = STRING_VERSION, species = STRING_TAXON,
                           score_threshold = STRING_SCORE_THRESHOLD,
                           input_directory = string_cache)

# =====================================================================
# STRING-style network builder
# =====================================================================
build_string_style_network <- function(common_table, title, file_prefix,
                                        direction_label) {
  if (nrow(common_table) == 0) {
    cat(sprintf("  [%s] No proteins — skipping.\n", title))
    return(NULL)
  }

  genes <- common_table$gene
  cat(sprintf("\n  [%s] %d proteins\n", title, length(genes)))

  # Map to STRING
  map_input <- as.data.frame(common_table[, c("gene", "mean_log2FC")])
  map_input <- map_input[!is.na(map_input$gene) & map_input$gene != "", ]
  mapped <- tryCatch(
    string_db$map(map_input, "gene", removeUnmappedRows = TRUE),
    error = function(e) { cat(sprintf("    MAP ERROR: %s\n", conditionMessage(e))); NULL }
  )
  if (is.null(mapped) || nrow(mapped) < 2) {
    cat("    Too few mapped. Skipping.\n")
    return(NULL)
  }

  mapped_ids <- unique(mapped$STRING_id)

  # Get physical interactions among these proteins
  edges <- phys[phys$from %in% mapped_ids & phys$to %in% mapped_ids, ]
  cat(sprintf("    Mapped: %d/%d | Interactions: %d\n",
              nrow(mapped), length(genes), nrow(edges)))

  if (nrow(edges) < 1) {
    cat("    No interactions. Skipping.\n")
    return(NULL)
  }

  # Build graph
  vertex_df <- data.frame(name = mapped_ids, stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                              vertices = vertex_df)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Keep largest connected component
  comps <- decompose(g)
  if (length(comps) > 1) {
    comp_sizes <- sapply(comps, vcount)
    g <- comps[[which.max(comp_sizes)]]
    cat(sprintf("    Largest component: %d proteins\n", vcount(g)))
  }

  # Annotate
  string_to_gene <- setNames(mapped$gene, mapped$STRING_id)
  string_to_fc   <- setNames(mapped$mean_log2FC, mapped$STRING_id)

  V(g)$gene_name <- sapply(V(g)$name, function(id) {
    gn <- string_to_gene[id]
    if (is.na(gn)) id else gn
  })
  V(g)$mean_log2FC <- as.numeric(string_to_fc[match(V(g)$name, names(string_to_fc))])

  # Assign edge weights from combined_score (for layout spacing)
  E(g)$weight <- edges$combined_score[match(paste(ends(g, E(g))[,1],
                                                   ends(g, E(g))[,2], sep="-"),
                                             paste(edges$from, edges$to, sep="-"))]
  # Some edges might not match due to simplification; fill NAs
  E(g)$weight[is.na(E(g)$weight)] <- median(E(g)$weight, na.rm = TRUE)

  # Edge confidence (0-1 scale for line width)
  E(g)$confidence <- E(g)$weight / 1000

  # ---- STRING-style coloring ----
  # Color nodes by mean_log2FC intensity (red=high positive, blue=high negative)
  fc_range <- range(V(g)$mean_log2FC, na.rm = TRUE)
  V(g)$color <- sapply(V(g)$mean_log2FC, function(fc) {
    if (is.na(fc)) return("grey80")
    if (fc > 0) {
      # Enriched: orange-red gradient by intensity
      intensity <- min(1, abs(fc) / max(abs(fc_range), na.rm = TRUE))
      rgb(0.85, 0.37 - 0.3 * intensity, 0, 0.5)  # transparent orange→red
    } else {
      # Depleted: blue gradient
      intensity <- min(1, abs(fc) / max(abs(fc_range), na.rm = TRUE))
      rgb(0, 0.27 + 0.3 * intensity, 0.67, 0.5)  # transparent navy
    }
  })

  cat(sprintf("    Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

  # ---- Plot: STRING bubble style using ggraph ----
  set.seed(42)

  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(file_prefix)
  versioned <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

  plot_string_style <- function() {
    # Use igraph base graphics for STRING-like appearance
    # Force-directed layout with more spread
    layout_fr <- layout_with_fr(g, niter = 1000)

    par(mar = c(1, 1, 3, 1))

    # Scale layout to fill canvas
    layout_scaled <- scale(layout_fr, center = colMeans(layout_fr))

    plot(g,
         layout               = layout_scaled,
         vertex.shape         = "circle",
         vertex.size          = 22,
         vertex.color         = V(g)$color,
         vertex.frame.color   = NA,
         vertex.label         = V(g)$gene_name,
         vertex.label.cex     = 0.65,
         vertex.label.font    = 2,
         vertex.label.color   = "black",
         vertex.label.family  = "sans",
         # STRING-style: edges colored by confidence
         edge.color           = rgb(0.5, 0.5, 0.5, 0.3 + 0.5 * E(g)$confidence),
         edge.width           = 0.5 + 2.5 * E(g)$confidence,
         edge.curved          = 0.1,
         main                 = title,
         asp                  = 0)  # fill canvas

    # Legend
    legend("bottomright",
           legend = c(paste("Enriched", direction_label),
                      paste("Depleted", direction_label),
                      "High confidence edge",
                      "Low confidence edge"),
           col    = c(rgb(0.85, 0.27, 0, 0.5),
                      rgb(0, 0.47, 0.67, 0.5),
                      rgb(0.5, 0.5, 0.5, 0.8),
                      rgb(0.5, 0.5, 0.5, 0.3)),
           lwd    = c(NA, NA, 3, 1),
           pch    = c(21, 21, NA, NA),
           pt.cex = 2.5,
           pt.bg  = c(rgb(0.85, 0.27, 0, 0.5),
                      rgb(0, 0.47, 0.67, 0.5)),
           cex    = 0.8, bty = "n")
  }

  grDevices::png(png_path, width = 14, height = 12, units = "in", res = FIG_DPI)
  plot_string_style()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 14, height = 12)
  plot_string_style()
  grDevices::dev.off()

  cat(sprintf("    Saved: %s\n", basename(png_path)))
  cat(sprintf("    Saved: %s\n", basename(pdf_path)))
  return(g)
}

# ---- Build STRING-style networks ----
if (nrow(enriched_table) > 0) {
  g_enr <- build_string_style_network(
    enriched_table,
    title = "Common RA-Enriched Proteins — STRING Physical Network",
    file_prefix = "RA_common_STRING_style_enriched",
    direction_label = "with RA"
  )
}

if (nrow(depleted_table) > 0) {
  g_dep <- build_string_style_network(
    depleted_table,
    title = "Common RA-Depleted Proteins — STRING Physical Network",
    file_prefix = "RA_common_STRING_style_depleted",
    direction_label = "with RA"
  )
}

cat("\n=========================================\n")
cat(" STRING-Style Networks Complete!\n")
cat("=========================================\n")
cat("  Output: output/figures/RA_common_STRING_style_*.png\n")
cat("=========================================\n")
