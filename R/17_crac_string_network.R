###############################################################################
# 17_crac_string_network.R
# STRING protein-protein interaction network for CRAC RNA interactome data.
#
# WHAT THIS DOES:
#   1. Loads the CRAC data (FLAG-TRIP4 RNA interactome)
#   2. Identifies significant CRAC hits (TRIP4-bound RNAs/proteins)
#   3. Maps them to the STRING database for protein-protein interactions
#   4. Identifies "seed" proteins using Lydia's criteria:
#      - log2FC > 2 AND p < 0.000001 (very high confidence)
#      - OR log2FC > 7 (extreme fold change)
#   5. Creates a network visualization showing:
#      - Seeds (large orange nodes)
#      - Candidates (smaller green nodes)
#      - Connection edges from STRING
#
# OUTPUT:
#   - Network plot (PNG + PDF)
#   - Candidate table (CSV ranked by connections)
#
# Usage:
#   make crac-network
###############################################################################
cat("\n=========================================\n")
cat(" CRAC STRING Network Analysis\n")
cat("=========================================\n\n")

# ---- Load CRAC data ----
crac_path <- file.path(DATA_DIR, "FLAG-TRIP4_list_CRACdata.csv")

if (!file.exists(crac_path)) {
  cat("ERROR: CRAC data file not found:", crac_path, "\n")
  cat("Expected at: data/FLAG-TRIP4_list_CRACdata.csv\n")
  quit(status = 1)
}

cat("Loading CRAC data...\n")
crac_df <- readr::read_csv(crac_path, show_col_types = FALSE)
cat(sprintf("  Raw rows: %d\n", nrow(crac_df)))

# CRAC data uses different column names than mass spec pipeline
crac_gene_col <- CRAC_GENE_COL    # "external_gene_name"
crac_fc_col   <- CRAC_LOG2FC_COL  # "logFC"
crac_padj_col <- CRAC_PADJ_COL    # "FDR"

# Verify columns exist
missing <- c()
if (!(crac_gene_col %in% colnames(crac_df))) missing <- c(missing, crac_gene_col)
if (!(crac_fc_col %in% colnames(crac_df)))   missing <- c(missing, crac_fc_col)
if (!(crac_padj_col %in% colnames(crac_df))) missing <- c(missing, crac_padj_col)

if (length(missing) > 0) {
  cat("ERROR: Missing columns in CRAC data:", paste(missing, collapse = ", "), "\n")
  cat("Available columns:", paste(colnames(crac_df), collapse = ", "), "\n")
  quit(status = 1)
}

# Build clean data frame with standardized column names
crac_clean <- data.frame(
  gene = crac_df[[crac_gene_col]],
  log2FC = crac_df[[crac_fc_col]],
  padj = crac_df[[crac_padj_col]],
  stringsAsFactors = FALSE
)

# Remove rows with missing gene names or NA values
crac_clean <- crac_clean[!is.na(crac_clean$gene) & crac_clean$gene != "" &
                           !is.na(crac_clean$log2FC) & !is.na(crac_clean$padj), ]
cat(sprintf("  Clean rows: %d proteins\n", nrow(crac_clean)))

# ---- Identify significant CRAC hits ----
# Significant = padj < 0.05 AND |log2FC| > 0.5
crac_sig <- crac_clean[crac_clean$padj < P_VALUE_CUTOFF &
                         abs(crac_clean$log2FC) > LOG2FC_CUTOFF, ]
crac_sig_genes <- unique(crac_sig$gene)
cat(sprintf("\n  Significant CRAC hits: %d proteins\n", length(crac_sig_genes)))

if (length(crac_sig_genes) < 5) {
  cat("  WARNING: Fewer than 5 significant CRAC hits. Network may be sparse.\n")
}

# ---- Classify seeds (Lydia's criteria) ----
# Seeds = very high confidence hits:
#   log2FC > 2 AND p < 0.000001
#   OR log2FC > 7
SEED_PADJ <- 0.000001   # 1e-6 — very stringent
SEED_LOG2FC <- 2         # 4-fold change minimum
SEED_LOG2FC_EXTREME <- 7 # extreme fold change

is_seed <- (crac_sig$padj < SEED_PADJ & crac_sig$log2FC > SEED_LOG2FC) |
           (crac_sig$log2FC > SEED_LOG2FC_EXTREME)

seed_genes <- unique(crac_sig$gene[is_seed])
candidate_genes <- setdiff(crac_sig_genes, seed_genes)

cat(sprintf("  Seeds (high-confidence): %d\n", length(seed_genes)))
cat(sprintf("  Candidates:              %d\n", length(candidate_genes)))

if (length(seed_genes) > 0) {
  cat(sprintf("    Seeds: %s\n", paste(seed_genes, collapse = ", ")))
}

# ---- Map to STRING database ----
cat("\nMapping to STRING database (human, v", STRING_VERSION, ")...\n", sep = "")
cat("  (First run downloads ~100MB — this is normal)\n")

# Sanitize cache dir to prevent colon issues on Windows/OneDrive
string_cache_dir <- file.path(OUTPUT_DIR, "string_cache")
dir.create(string_cache_dir, showWarnings = FALSE, recursive = TRUE)

string_map <- tryCatch({
  STRINGdb::STRINGdb$new(
    version = STRING_VERSION,
    species = STRING_TAXON,
    score_threshold = STRING_SCORE_THRESHOLD,
    input_directory = string_cache_dir
  )
}, error = function(e) {
  cat("ERROR: Cannot connect to STRING database:\n", conditionMessage(e), "\n")
  quit(status = 1)
})

# Map CRAC genes to STRING identifiers
cat("Mapping CRAC gene symbols to STRING IDs...\n")
string_hits <- string_map$map(crac_sig, "gene", removeUnmappedRows = TRUE)

n_mapped <- sum(!is.na(string_hits$STRING_id))
cat(sprintf("  Mapped: %d / %d genes to STRING\n", n_mapped, nrow(crac_sig)))

if (n_mapped < 3) {
  cat("  WARNING: Too few genes mapped to STRING for meaningful network.\n")
  cat("  Saving gene list only.\n")
  save_table(data.frame(gene = crac_sig_genes, category = "CRAC_significant"),
             "CRAC_network_genes")
  quit(status = 0)
}

mapped_genes <- string_hits$STRING_id[!is.na(string_hits$STRING_id)]

# ---- Get STRING interactions ----
# Use local physical links file if available (avoids 100MB download on Windows)
phys_file <- file.path(DATA_DIR, "9606.protein.physical.links.v12.0.txt")

if (file.exists(phys_file)) {
  cat("Loading local STRING physical interactions file...\n")
  phys <- read.table(phys_file, header = TRUE, stringsAsFactors = FALSE)
  names(phys) <- c("from", "to", "combined_score")
  cat(sprintf("  Loaded %d physical interactions\n", nrow(phys)))
  # Filter to interactions involving our mapped genes
  interactions <- phys[phys$from %in% mapped_genes | phys$to %in% mapped_genes, ]
} else {
  cat("Querying STRING for protein interactions (online)...\n")
  interactions <- string_map$get_interactions(mapped_genes)
}
cat(sprintf("  Found %d interactions\n", nrow(interactions)))

if (nrow(interactions) == 0) {
  cat("  No interactions found. Saving gene list only.\n")
  save_table(data.frame(gene = crac_sig_genes, category = "CRAC_significant"),
             "CRAC_network_genes")
  quit(status = 0)
}

# ---- Build igraph network ----
cat("\nBuilding network visualization...\n")
library(igraph)

# Create edge list from STRING interactions
edges <- data.frame(
  from = interactions$from,
  to   = interactions$to,
  weight = interactions$combined_score / 1000  # Normalize to 0-1
)

# Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Mark seeds vs candidates
string_to_gene <- setNames(string_hits$gene, string_hits$STRING_id)
V(g)$name_display <- sapply(V(g)$name, function(id) {
  gene <- string_to_gene[id]
  if (is.na(gene)) id else gene
})

# Determine node type
seed_string_ids <- string_hits$STRING_id[string_hits$gene %in% seed_genes]
V(g)$is_seed <- V(g)$name %in% seed_string_ids

# Node colors: seeds = orange, candidates = green
V(g)$color <- ifelse(V(g)$is_seed,
                     GLOBAL_COLORS[["enriched_up"]],   # Vermillion orange
                     GLOBAL_COLORS[["known_ia"]])       # Bluish green
V(g)$size <- ifelse(V(g)$is_seed, 8, 4)

# ---- Network plot ----
# Use Fruchterman-Reingold layout for natural-looking network
set.seed(42)  # Reproducible layout
layout_fr <- layout_with_fr(g)

# Label only seeds (too many labels otherwise)
label_vec <- ifelse(V(g)$is_seed, V(g)$name_display, NA)

# Save as PNG
commit_hash <- get_git_hash()
png_path <- safe_filepath(FIGURE_DIR, paste0("crac_string_network_", commit_hash), ".png")
pdf_path <- safe_filepath(FIGURE_DIR, paste0("crac_string_network_", commit_hash), ".pdf")

grDevices::png(png_path, width = 12, height = 10, units = "in", res = 300)
plot(g,
     layout = layout_fr,
     vertex.label = label_vec,
     vertex.label.cex = 0.7,
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = "grey70",
     edge.width = 0.5,
     edge.curved = 0.1,
     main = "FLAG-TRIP4 CRAC RNA Interactome — STRING Network")
legend("topright",
       legend = c("Seeds (high-confidence)", "Candidates"),
       col = c(GLOBAL_COLORS[["enriched_up"]], GLOBAL_COLORS[["known_ia"]]),
       pch = 19, pt.cex = c(2, 1.5), cex = 0.9)
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 12, height = 10)
plot(g,
     layout = layout_fr,
     vertex.label = label_vec,
     vertex.label.cex = 0.7,
     vertex.label.font = 2,
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = "grey70",
     edge.width = 0.5,
     edge.curved = 0.1,
     main = "FLAG-TRIP4 CRAC RNA Interactome — STRING Network")
legend("topright",
       legend = c("Seeds (high-confidence)", "Candidates"),
       col = c(GLOBAL_COLORS[["enriched_up"]], GLOBAL_COLORS[["known_ia"]]),
       pch = 19, pt.cex = c(2, 1.5), cex = 0.9)
grDevices::dev.off()

cat(sprintf("  Saved: %s\n", basename(png_path)))
cat(sprintf("  Saved: %s\n", basename(pdf_path)))

# ---- Candidate table ranked by connections ----
cat("\nBuilding candidate table...\n")
degree_df <- data.frame(
  gene = V(g)$name_display,
  string_id = V(g)$name,
  is_seed = V(g)$is_seed,
  connections = degree(g),
  stringsAsFactors = FALSE
)
degree_df <- degree_df[order(-degree_df$connections), ]

# Merge back log2FC and padj from CRAC data
degree_df <- merge(degree_df, crac_clean[, c("gene", "log2FC", "padj")],
                   by = "gene", all.x = TRUE)
degree_df <- degree_df[order(-degree_df$connections), ]

save_table(degree_df, "CRAC_network_candidates")

cat(sprintf("\n  Top 10 candidates by connections:\n"))
for (i in seq_len(min(10, nrow(degree_df)))) {
  row <- degree_df[i, ]
  cat(sprintf("    %2d. %-12s (connections=%d, log2FC=%.2f, padj=%.2g, %s)\n",
              i, row$gene, row$connections, row$log2FC, row$padj,
              ifelse(row$is_seed, "SEED", "candidate")))
}

cat("\n=========================================\n")
cat(" CRAC STRING network analysis complete!\n")
cat("=========================================\n")
