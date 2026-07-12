###############################################################################
# 22_ra_common_analysis.R
# Retinoic Acid (RA) common protein analysis across two concentrations.
#
# WHAT THIS DOES:
#   Retinoic acid was applied at two concentrations in two TRIP4 TurboID batches:
#     - BK467: 0.2 µM RA  (experiment: BK467_TRIP4_RA02_vs_BK467_TRIP4)
#     - BK504: 0.4 µM RA  (experiment: BK504_TRIP4_RA04_vs_BK504_TRIP4)
#   Each compares "+RA" vs "-RA" within the SAME TRIP4 cell line, so they
#   isolate the RA treatment effect from the TRIP4-vs-WT effect.
#
#   This script finds proteins that respond CONSISTENTLY to RA across BOTH
#   concentrations — i.e. the INTERSECTION of significant hits in both:
#     1. Common RA-ENRICHED proteins (positive log2FC, padj < 0.05 in BOTH)
#     2. Common RA-DEPLETED proteins (negative log2FC, padj < 0.05 in BOTH)
#     3. STRING physical-interaction network of each common set (Lydia's method)
#     4. Venn diagram: RA02-enriched vs RA04-enriched (overlap visualization)
#     5. Gene lists exported as plain-text files for ShinyGO / external tools
#
# STRING network approach follows Lydia's local physical-links method
# (see R/21_lydia_network_volcano.R): use the local
# 9606.protein.physical.links.v12.0.txt file to build the induced subnetwork
# among the common proteins, avoiding the ~100MB STRINGdb download.
#
# Usage:
#   make ra-common           (full analysis + networks + venn + gene lists)
#   make ra-common-network   (alias — runs the same comprehensive analysis)
###############################################################################

cat("\n=========================================\n")
cat(" RA Common Protein Analysis (RA02 + RA04)\n")
cat("=========================================\n\n")

library(STRINGdb)
library(igraph)
library(VennDiagram)
library(grid)

# ---- Load all experiments ----
experiments <- load_all_experiments()

# ---- Define the two RA-effect experiments ----
# These are "+RA vs -RA within TRIP4" comparisons — the pure RA treatment effect.
RA02_EXP <- "BK467_TRIP4_RA02_vs_BK467_TRIP4"   # 0.2 µM RA, batch BK467
RA04_EXP <- "BK504_TRIP4_RA04_vs_BK504_TRIP4"   # 0.4 µM RA, batch BK504

df_ra02 <- find_experiment(experiments, RA02_EXP)
df_ra04 <- find_experiment(experiments, RA04_EXP)

if (is.null(df_ra02)) {
  cat("ERROR: Experiment not found:", RA02_EXP, "\n")
  cat("Available experiments:\n")
  for (n in names(experiments)) cat("  -", n, "\n")
  quit(status = 1)
}
if (is.null(df_ra04)) {
  cat("ERROR: Experiment not found:", RA04_EXP, "\n")
  cat("Available experiments:\n")
  for (n in names(experiments)) cat("  -", n, "\n")
  quit(status = 1)
}

cat(sprintf("  RA02 (%s): %d proteins\n", RA02_EXP, nrow(df_ra02)))
cat(sprintf("  RA04 (%s): %d proteins\n", RA04_EXP, nrow(df_ra04)))

# =====================================================================
# Helper: extract significant genes by DIRECTION
# =====================================================================
# get_significant_genes() (from utils.R) returns BOTH up and down combined.
# Here we need direction-specific sets for the intersection analysis.

get_enriched_genes <- function(df,
                               padj_cutoff = P_VALUE_CUTOFF,
                               log2fc_cutoff = LOG2FC_CUTOFF) {
  # Enriched with RA = significant AND positive log2FC (higher abundance with RA)
  genes <- df$gene[df$padj < padj_cutoff & df$log2FC > log2fc_cutoff]
  genes <- unique(genes[!is.na(genes)])
  return(genes)
}

get_depleted_genes <- function(df,
                               padj_cutoff = P_VALUE_CUTOFF,
                               log2fc_cutoff = LOG2FC_CUTOFF) {
  # Depleted with RA = significant AND negative log2FC (lower abundance with RA)
  genes <- df$gene[df$padj < padj_cutoff & df$log2FC < -log2fc_cutoff]
  genes <- unique(genes[!is.na(genes)])
  return(genes)
}

# =====================================================================
# STEP 1: Common RA-enriched proteins (intersection of both concentrations)
# =====================================================================
cat("\n--- Step 1: Common RA-enriched proteins ---\n")

ra02_enriched <- get_enriched_genes(df_ra02)
ra04_enriched <- get_enriched_genes(df_ra04)

cat(sprintf("  RA02-enriched (0.2 µM): %d proteins\n", length(ra02_enriched)))
cat(sprintf("  RA04-enriched (0.4 µM): %d proteins\n", length(ra04_enriched)))

common_enriched <- intersect(ra02_enriched, ra04_enriched)
cat(sprintf("  COMMON enriched (both): %d proteins\n", length(common_enriched)))
if (length(common_enriched) > 0) {
  cat(sprintf("    %s\n", paste(common_enriched, collapse = ", ")))
}

# =====================================================================
# STEP 2: Common RA-depleted proteins (intersection of both concentrations)
# =====================================================================
cat("\n--- Step 2: Common RA-depleted proteins ---\n")

ra02_depleted <- get_depleted_genes(df_ra02)
ra04_depleted <- get_depleted_genes(df_ra04)

cat(sprintf("  RA02-depleted (0.2 µM): %d proteins\n", length(ra02_depleted)))
cat(sprintf("  RA04-depleted (0.4 µM): %d proteins\n", length(ra04_depleted)))

common_depleted <- intersect(ra02_depleted, ra04_depleted)
cat(sprintf("  COMMON depleted (both): %d proteins\n", length(common_depleted)))
if (length(common_depleted) > 0) {
  cat(sprintf("    %s\n", paste(common_depleted, collapse = ", ")))
}

# =====================================================================
# STEP 3: Build summary tables with log2FC from both experiments
# =====================================================================
cat("\n--- Step 3: Build summary tables ---\n")

build_common_table <- function(genes, df_a, df_b, label_a, label_b) {
  if (length(genes) == 0) return(data.frame())
  sub_a <- df_a[df_a$gene %in% genes, c("gene", "log2FC", "padj")]
  sub_b <- df_b[df_b$gene %in% genes, c("gene", "log2FC", "padj")]
  colnames(sub_a) <- c("gene", paste0("log2FC_", label_a), paste0("padj_", label_a))
  colnames(sub_b) <- c("gene", paste0("log2FC_", label_b), paste0("padj_", label_b))
  # Deduplicate (in case a gene appears multiple times in one experiment)
  sub_a <- sub_a[!duplicated(sub_a$gene), ]
  sub_b <- sub_b[!duplicated(sub_b$gene), ]
  merged <- merge(sub_a, sub_b, by = "gene", all = TRUE)
  merged$mean_log2FC <- rowMeans(merged[, c(
    paste0("log2FC_", label_a), paste0("log2FC_", label_b))], na.rm = TRUE)
  merged <- merged[order(-abs(merged$mean_log2FC)), ]
  return(merged)
}

common_enriched_table <- build_common_table(
  common_enriched, df_ra02, df_ra04, "RA02", "RA04")
common_depleted_table <- build_common_table(
  common_depleted, df_ra02, df_ra04, "RA02", "RA04")

if (nrow(common_enriched_table) > 0) {
  save_table(common_enriched_table, "RA_common_enriched_RA02_RA04")
}
if (nrow(common_depleted_table) > 0) {
  save_table(common_depleted_table, "RA_common_depleted_RA02_RA04")
}

# ---- Export gene lists as plain-text files (for ShinyGO, DAVID, etc.) ----
export_gene_list <- function(genes, filename) {
  if (length(genes) == 0) {
    cat(sprintf("  Skipped (empty set): %s\n", filename))
    return(invisible(NULL))
  }
  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(filename)
  versioned <- paste0(safe_name, "_", commit_hash)
  filepath <- file.path(TABLE_DIR, paste0(versioned, ".txt"))
  write.table(data.frame(gene = genes), file = filepath,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Saved: %s (%d genes)\n", basename(filepath), length(genes)))
  invisible(filepath)
}

export_gene_list(common_enriched, "RA_common_enriched_genes")
export_gene_list(common_depleted, "RA_common_depleted_genes")
export_gene_list(ra02_enriched, "RA02_enriched_genes")
export_gene_list(ra04_enriched, "RA04_enriched_genes")
export_gene_list(ra02_depleted, "RA02_depleted_genes")
export_gene_list(ra04_depleted, "RA04_depleted_genes")

# =====================================================================
# STEP 4: Venn diagram — RA02 vs RA04 (enriched and depleted)
# =====================================================================
cat("\n--- Step 4: Venn diagrams (RA02 vs RA04) ---\n")

make_ra_venn <- function(set_a, set_b, label_a, label_b, title, file_prefix) {
  overlap_count <- length(intersect(set_a, set_b))
  cat(sprintf("  %s: %d | %s: %d | overlap: %d\n",
              label_a, length(set_a), label_b, length(set_b), overlap_count))

  vp <- draw.pairwise.venn(
    area1          = length(set_a),
    area2          = length(set_b),
    cross.area     = overlap_count,
    category       = c(label_a, label_b),
    fill           = c(GLOBAL_COLORS[["venn_a_only"]], GLOBAL_COLORS[["venn_b_only"]]),
    alpha          = rep(0.5, 2),
    cat.cex        = 1.4,
    cex            = 1.8,
    fontfamily     = "sans",
    cat.fontfamily = "sans",
    col            = "transparent",
    cat.pos        = c(-30, 30),
    cat.dist       = c(0.06, 0.06),
    margin         = 0.08,
    ind            = FALSE
  )

  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(file_prefix)
  versioned <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

  # Draw to both PNG and PDF
  for (ext_path in list(png_path, pdf_path)) {
    if (grepl("\\.png$", ext_path)) {
      grDevices::png(ext_path, width = 7, height = 6, units = "in", res = FIG_DPI)
    } else {
      grDevices::pdf(ext_path, width = 7, height = 6)
    }
    grid.draw(vp)
    pushViewport(viewport())
    grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 13, fontface = "bold"))
    popViewport()
    grDevices::dev.off()
  }
  cat(sprintf("  Saved: %s\n  Saved: %s\n", basename(png_path), basename(pdf_path)))
}

# Enriched overlap (primary, per Aruna's request)
make_ra_venn(ra02_enriched, ra04_enriched,
             label_a = "RA02 (0.2uM)",
             label_b = "RA04 (0.4uM)",
             title = "RA-Enriched Proteins: 0.2uM vs 0.4uM",
             file_prefix = "RA_common_venn_enriched_RA02_vs_RA04")

# Depleted overlap (for completeness)
make_ra_venn(ra02_depleted, ra04_depleted,
             label_a = "RA02 (0.2uM)",
             label_b = "RA04 (0.4uM)",
             title = "RA-Depleted Proteins: 0.2uM vs 0.4uM",
             file_prefix = "RA_common_venn_depleted_RA02_vs_RA04")

# =====================================================================
# STEP 5: STRING physical-interaction networks (Lydia's method)
# =====================================================================
cat("\n--- Step 5: STRING physical-interaction networks ---\n")
cat("(Using Lydia's local physical-links approach)\n\n")

# ---- Load local STRING physical links file (avoids ~100MB download) ----
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")
phys <- NULL

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))
} else {
  cat("NOTE: Local physical links file not found at:\n")
  cat(sprintf("  %s\n", phys_file))
  cat("  Will query STRINGdb online for interactions instead.\n")
}

# ---- Initialize STRINGdb for gene-name -> STRING-ID mapping ----
cat("Initializing STRINGdb (v", STRING_VERSION, ", human, score >= ",
    STRING_SCORE_THRESHOLD, ")...\n", sep = "")
string_cache <- file.path(OUTPUT_DIR, "string_cache")
dir.create(string_cache, recursive = TRUE, showWarnings = FALSE)

string_db <- STRINGdb$new(
  version         = STRING_VERSION,
  species         = STRING_TAXON,
  score_threshold = STRING_SCORE_THRESHOLD,
  input_directory = string_cache
)

# ---- Helper: build & plot a STRING network for a set of common proteins ----
# Lydia's method (induced subnetwork among the common proteins):
#   1. Map gene symbols -> STRING IDs
#   2. Find physical interactions among those proteins (induced subnetwork)
#   3. Build igraph with force-directed (Fruchterman-Reingold) layout
#   4. Color nodes by log2FC direction (vermillion = up, navy = down)
build_ra_string_network <- function(common_table, title, file_prefix) {

  if (nrow(common_table) == 0) {
    cat(sprintf("  [%s] No proteins to map — skipping.\n", title))
    return(NULL)
  }

  genes <- common_table$gene
  cat(sprintf("\n  [%s] %d proteins to map to STRING\n", title, length(genes)))

  # Map genes -> STRING IDs.
  # as.data.frame() avoids the "incorrect number of dimensions" tibble error.
  map_input <- as.data.frame(common_table[, c("gene", "mean_log2FC")])
  map_input <- map_input[!is.na(map_input$gene) & map_input$gene != "", ]

  mapped <- tryCatch({
    string_db$map(map_input, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
    cat("    Retrying with gene names only...\n")
    simple_input <- as.data.frame(data.frame(gene = unique(map_input$gene),
                                             stringsAsFactors = FALSE))
    tryCatch({
      string_db$map(simple_input, "gene", removeUnmappedRows = TRUE)
    }, error = function(e2) {
      cat(sprintf("    Still failing: %s\n", conditionMessage(e2)))
      return(NULL)
    })
  })

  if (is.null(mapped) || nrow(mapped) == 0) {
    cat("    No proteins mapped to STRING — skipping network.\n")
    return(NULL)
  }
  cat(sprintf("    Mapped %d / %d proteins to STRING\n",
              nrow(mapped), length(genes)))

  mapped_ids <- unique(mapped$STRING_id[!is.na(mapped$STRING_id)])
  if (length(mapped_ids) < 2) {
    cat("    Fewer than 2 mapped proteins — no network possible.\n")
    return(NULL)
  }

  # Find physical interactions among the mapped proteins (induced subnetwork)
  if (!is.null(phys)) {
    edges <- phys[phys$from %in% mapped_ids & phys$to %in% mapped_ids, ]
    cat(sprintf("    Physical interactions among mapped proteins: %d edges\n",
                nrow(edges)))
  } else {
    # Fallback: query STRINGdb online
    interactions <- string_db$get_interactions(mapped_ids)
    edges <- interactions[interactions$from %in% mapped_ids &
                          interactions$to   %in% mapped_ids, ]
    cat(sprintf("    STRINGdb interactions among mapped proteins: %d edges\n",
                nrow(edges)))
  }

  if (nrow(edges) < 1) {
    cat("    No interactions among these proteins — skipping network plot.\n")
    cat("    (These common proteins do not physically interact in STRING.)\n")
    return(NULL)
  }

  # ---- Build igraph (include ALL mapped proteins as vertices) ----
  vertex_df <- data.frame(name = mapped_ids, stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                             vertices = vertex_df)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Annotate vertices with gene names and mean log2FC
  string_to_gene <- setNames(mapped$gene, mapped$STRING_id)
  string_to_fc   <- setNames(mapped$mean_log2FC, mapped$STRING_id)

  V(g)$gene_name <- sapply(V(g)$name, function(id) {
    gn <- string_to_gene[id]
    if (is.na(gn)) id else gn
  })
  V(g)$mean_log2FC <- as.numeric(string_to_fc[match(V(g)$name, names(string_to_fc))])

  # ---- Color nodes by log2FC direction ----
  # Vermillion (#D55E00) = up/enriched, Navy (#0072B2) = down/depleted
  V(g)$color <- ifelse(!is.na(V(g)$mean_log2FC) & V(g)$mean_log2FC > 0,
                       GLOBAL_COLORS[["enriched_up"]],   # vermillion
                       GLOBAL_COLORS[["ascc_core"]])     # navy
  V(g)$size        <- 5
  V(g)$label       <- V(g)$gene_name
  V(g)$frame.color <- NA

  cat(sprintf("    Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

  # ---- Plot with force-directed (Fruchterman-Reingold) layout ----
  set.seed(42)  # reproducible layout
  layout_fr <- layout_with_fr(g)

  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(file_prefix)
  versioned <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

  plot_network <- function() {
    plot(g,
         layout             = layout_fr,
         vertex.label.cex   = 0.6,
         vertex.label.font  = 2,
         vertex.label.color = "black",
         edge.color         = "grey80",
         edge.width         = 0.5,
         main               = title)
    legend("topright",
           legend = c("Up (enriched with RA)", "Down (depleted with RA)"),
           col    = c(GLOBAL_COLORS[["enriched_up"]], GLOBAL_COLORS[["ascc_core"]]),
           pch    = 19, pt.cex = 2, cex = 0.9, bty = "n")
  }

  grDevices::png(png_path, width = 12, height = 10, units = "in", res = FIG_DPI)
  plot_network()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 12, height = 10)
  plot_network()
  grDevices::dev.off()

  cat(sprintf("    Saved: %s\n", basename(png_path)))
  cat(sprintf("    Saved: %s\n", basename(pdf_path)))

  return(g)
}

# ---- Network of common RA-enriched proteins ----
g_enriched <- build_ra_string_network(
  common_enriched_table,
  title = "Common RA-Enriched Proteins - STRING Physical Network (RA02 & RA04)",
  file_prefix = "RA_common_string_network_enriched"
)

# ---- Network of common RA-depleted proteins ----
g_depleted <- build_ra_string_network(
  common_depleted_table,
  title = "Common RA-Depleted Proteins - STRING Physical Network (RA02 & RA04)",
  file_prefix = "RA_common_string_network_depleted"
)

# =====================================================================
# STEP 6: Circular network plots (Lydia's Panel C style)
# =====================================================================
cat("\n--- Step 6: Circular network plots (Lydia Panel C style) ---\n")

build_circular_network <- function(common_table, title, file_prefix) {
  if (is.null(common_table) || nrow(common_table) == 0) {
    cat(sprintf("  [%s] No proteins — skipping circular plot.\n", title))
    return(NULL)
  }

  # Re-use the graph if STRING mapped successfully
  genes <- common_table$gene
  cat(sprintf("\n  [%s] %d proteins for circular layout\n", title, length(genes)))

  map_input <- as.data.frame(common_table[, c("gene", "mean_log2FC")])
  map_input <- map_input[!is.na(map_input$gene) & map_input$gene != "", ]

  mapped <- tryCatch({
    string_db$map(map_input, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    cat(sprintf("    ERROR mapping: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(mapped) || nrow(mapped) < 2) {
    cat("    Fewer than 2 mapped — skipping.\n")
    return(NULL)
  }

  mapped_ids <- unique(mapped$STRING_id)

  # Get physical interactions among these proteins
  if (!is.null(phys)) {
    edges <- phys[phys$from %in% mapped_ids & phys$to %in% mapped_ids, ]
  } else {
    interactions <- string_db$get_interactions(mapped_ids)
    edges <- interactions[interactions$from %in% mapped_ids &
                          interactions$to %in% mapped_ids, ]
  }

  if (nrow(edges) < 1) {
    cat("    No interactions — proteins don't form a network.\n")
    cat("    Still generating circular protein list (no edges).\n")
  }

  # Build graph — include ALL mapped proteins as vertices
  vertex_df <- data.frame(name = mapped_ids, stringsAsFactors = FALSE)
  if (nrow(edges) > 0) {
    g <- graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                               vertices = vertex_df)
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  } else {
    g <- make_graph("empty", n = length(mapped_ids))
    V(g)$name <- mapped_ids
  }

  # Add isolated vertices that have no edges
  isolated <- mapped_ids[!mapped_ids %in% c(edges$from, edges$to)]
  for (id in isolated) {
    g <- g + vertices(id)
  }

  # Annotate
  string_to_gene <- setNames(mapped$gene, mapped$STRING_id)
  string_to_fc   <- setNames(mapped$mean_log2FC, mapped$STRING_id)

  V(g)$gene_name <- sapply(V(g)$name, function(id) {
    gn <- string_to_gene[id]
    if (is.na(gn)) id else gn
  })
  V(g)$mean_log2FC <- as.numeric(string_to_fc[match(V(g)$name, names(string_to_fc))])

  # Color: vermillion = enriched, navy = depleted
  V(g)$color <- ifelse(!is.na(V(g)$mean_log2FC) & V(g)$mean_log2FC > 0,
                       GLOBAL_COLORS[["enriched_up"]],
                       GLOBAL_COLORS[["ascc_core"]])
  V(g)$size        <- 6
  V(g)$label       <- V(g)$gene_name
  V(g)$frame.color <- NA

  cat(sprintf("    Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

  # Circular layout (Lydia Panel C style) — labels OUTSIDE circle
  set.seed(42)
  layout_circle <- layout_in_circle(g, order = order(V(g)$gene_name))

  # Compute radial angles so labels point outward from center
  angles <- atan2(layout_circle[,2], layout_circle[,1])
  # Place labels radially outside: degree = -angle (igraph convention)
  label_degrees <- -angles

  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(file_prefix)
  versioned <- paste0(safe_name, "_", commit_hash)
  png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
  pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

  plot_circular <- function() {
    plot(g,
         layout             = layout_circle,
         vertex.label.cex   = 0.55,
         vertex.label.font  = 2,
         vertex.label.color = "black",
         vertex.label.dist  = 0,      # distance from node (0 = on edge)
         vertex.label.degree = label_degrees,  # radially outward
         edge.color         = "grey80",
         edge.width         = 0.4,
         main               = title)
    legend("bottomright",
           legend = c("Enriched with RA", "Depleted with RA"),
           col    = c(GLOBAL_COLORS[["enriched_up"]], GLOBAL_COLORS[["ascc_core"]]),
           pch    = 19, pt.cex = 2, cex = 0.9, bty = "n")
  }

  grDevices::png(png_path, width = 12, height = 12, units = "in", res = FIG_DPI)
  plot_circular()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 12, height = 12)
  plot_circular()
  grDevices::dev.off()

  cat(sprintf("    Saved: %s\n", basename(png_path)))
  cat(sprintf("    Saved: %s\n", basename(pdf_path)))
  return(g)
}

# Circular network for common enriched
g_enriched_circle <- build_circular_network(
  common_enriched_table,
  title = "Common RA-Enriched Proteins — Circular Network (RA02 & RA04)",
  file_prefix = "RA_common_circular_enriched"
)

# Circular network for common depleted
g_depleted_circle <- build_circular_network(
  common_depleted_table,
  title = "Common RA-Depleted Proteins — Circular Network (RA02 & RA04)",
  file_prefix = "RA_common_circular_depleted"
)

# =====================================================================
# STEP 7: GO enrichment analysis
# =====================================================================
cat("\n--- Step 7: GO enrichment analysis ---\n")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

ONT_NAMES <- c("BP" = "Biological Process",
               "CC" = "Cellular Component",
               "MF" = "Molecular Function")

# Universe = all proteins detected in RA02 experiment
go_universe <- unique(df_ra02$gene[!is.na(df_ra02$gene)])
cat(sprintf("  GO universe: %d proteins\n", length(go_universe)))

run_common_go <- function(genes, set_label) {
  cat(sprintf("\n  [%s] %d genes\n", set_label, length(genes)))

  if (length(genes) < 5) {
    cat("    Too few genes (<5). Skipping GO.\n")
    return(NULL)
  }

  for (ont in c("BP", "CC", "MF")) {
    cat(sprintf("    [%s] enrichGO... ", ont))

    ego <- tryCatch({
      enrichGO(
        gene          = genes,
        universe      = go_universe,
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
      next
    }

    ego_simplified <- tryCatch({
      simplify(ego, cutoff = 0.7)
    }, error = function(e) ego)

    res <- as.data.frame(ego_simplified)
    cat(sprintf("%d terms\n", nrow(res)))

    # Save table
    save_table(res, sprintf("RA_common_GO_%s_%s", set_label, ont))

    # Dotplot
    n_show <- min(15, nrow(res))
    fig_height <- max(6, n_show * 0.4)

    p <- dotplot(ego_simplified, showCategory = n_show) +
      ggplot2::labs(
        title = sprintf("Common RA-%s — %s",
                        set_label, ONT_NAMES[ont]),
        x = "GeneRatio", color = "p.adjust"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    fig_name <- sprintf("RA_common_GO_%s_%s_dotplot", set_label, ont)
    save_figure(p, fig_name, width = 10, height = fig_height)
  }
}

# GO for common enriched
run_common_go(common_enriched, "enriched")

# GO for common depleted
run_common_go(common_depleted, "depleted")

# =====================================================================
# STEP 8: Bidirectional GO plot — common enriched (right) vs depleted (left)
# =====================================================================
cat("\n--- Step 8: Bidirectional GO (enriched vs depleted) ---\n")

run_bidir_common_go <- function(genes, direction, ont) {
  if (length(genes) < 5) {
    cat(sprintf("  Skipping %s %s: too few genes (%d)\n", direction, ont, length(genes)))
    return(NULL)
  }

  cat(sprintf("  Running enrichGO: %s, %s (%d genes)...\n", direction, ont, length(genes)))

  ego <- enrichGO(
    gene          = genes,
    universe      = go_universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 2,
    maxGSSize     = 5000
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)

  ego_s <- tryCatch(simplify(ego, cutoff = 0.7), error = function(e) ego)
  res <- as.data.frame(ego_s)
  if (nrow(res) == 0) return(NULL)

  res <- res[order(res$p.adjust), ]
  res <- head(res, 15)

  res$GeneRatioNum <- sapply(res$GeneRatio, function(x) {
    p <- strsplit(as.character(x), "/")[[1]]
    if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2]) else NA
  })

  res$signed_GeneRatio <- if (direction == "enriched") res$GeneRatioNum else -res$GeneRatioNum
  res$direction <- direction
  res$neg_log10_padj <- -log10(res$p.adjust)
  res$short_Description <- sapply(res$Description, function(x) {
    if (nchar(x) > 55) paste0(substr(x, 1, 52), "...") else x
  })
  return(res)
}

for (ont in c("BP", "CC", "MF")) {
  cat(sprintf("\n  [%s] bidirectional...\n", ont))

  en_res <- run_bidir_common_go(common_enriched, "enriched", ont)
  de_res <- run_bidir_common_go(common_depleted, "depleted", ont)

  combined <- rbind(en_res, de_res)
  if (is.null(combined) || nrow(combined) == 0) {
    cat("    No terms. Skipping.\n")
    next
  }

  combined <- combined[order(combined$signed_GeneRatio, decreasing = TRUE), ]
  combined$short_Description <- factor(combined$short_Description,
                                        levels = rev(combined$short_Description))

  p <- ggplot2::ggplot(combined,
                        ggplot2::aes(x = signed_GeneRatio,
                                     y = short_Description,
                                     color = neg_log10_padj,
                                     size = Count)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradientn(
      colors = c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"),
      name = expression(-Log[10]~(adjusted~italic(p)~value))
    ) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid",
                        color = "grey50", linewidth = 0.3) +
    ggplot2::labs(
      x = expression(Signed~Gene~Ratio~(left~"="-~depleted~"|"~right~"="+~enriched)),
      y = NULL,
      title = sprintf("Common RA Proteins — Bidirectional GO (%s)", ONT_NAMES[ont]),
      subtitle = sprintf("Enriched: %d genes → %d terms | Depleted: %d genes → %d terms",
                         length(common_enriched),
                         sum(combined$direction == "enriched"),
                         length(common_depleted),
                         sum(combined$direction == "depleted"))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, color = "grey30"),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")
    )

  p <- p +
    ggplot2::annotate("text",
                      x = max(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "ENRICHED", color = "#D55E00",
                      fontface = "bold", size = 4, vjust = 0) +
    ggplot2::annotate("text",
                      x = min(combined$signed_GeneRatio, na.rm = TRUE) * 0.8,
                      y = 0.5, label = "DEPLETED", color = "#0072B2",
                      fontface = "bold", size = 4, vjust = 0)

  fig_name <- sprintf("RA_common_bidirectional_GO_%s", ont)
  save_figure(p, fig_name, width = 10,
              height = max(6, nrow(combined) * 0.35))
  cat(sprintf("    Saved: %s\n", fig_name))
}
cat("\n=========================================\n")
cat(" RA Common Protein Analysis Complete!\n")
cat("=========================================\n")
cat(sprintf("  Common RA-enriched: %d proteins\n", length(common_enriched)))
cat(sprintf("  Common RA-depleted: %d proteins\n", length(common_depleted)))
cat("\nOutputs:\n")
cat("  - Venn diagrams: output/figures/\n")
cat("  - STRING force-directed networks: output/figures/\n")
cat("  - Circular networks (Lydia Panel C): output/figures/\n")
cat("  - GO dotplots (BP, CC, MF): output/figures/\n")
cat("  - Summary tables (CSV): output/tables/\n")
cat("  - Gene lists (TXT for ShinyGO): output/tables/\n")
cat("\n")
