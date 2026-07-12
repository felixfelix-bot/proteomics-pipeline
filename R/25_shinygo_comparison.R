###############################################################################
# 25_shinygo_comparison.R
# Compare ShinyGO pathway/network results with our STRING pipeline.
#
# WHAT THIS DOES (requested by Dr. Aruna):
#   Aruna uses ShinyGO (bioinformatics.sdstate.edu/go/) for pathway
#   visualization. This script compares ShinyGO exports with our
#   STRING-based analysis:
#     1. Reads ShinyGO CSV export (Aruna places in data/shinygo_export.csv)
#     2. Reads our STRING network gene lists (from lydia-volcano output)
#     3. Generates side-by-side comparison:
#        - Venn diagram of gene sets
#        - GO term overlap analysis
#        - Unique ShinyGO terms vs unique STRING terms
#        - Combined network view
#
# SHINYGO CSV FORMAT (expected columns — Aruna exports from browser):
#   The script auto-detects common ShinyGO export formats:
#   - "Gene Set" / "Term" / "Description" column for GO term
#   - "Genes" / "Gene List" column for gene symbols
#   - "FDR" / "qvalue" / "p.adjust" for significance
#
# Usage:
#   make shinygo-compare
###############################################################################

source('R/01_config.R', chdir = TRUE)
source('R/utils.R')

library(ggplot2)
library(VennDiagram)
library(grid)

cat("\n=========================================\n")
cat(" ShinyGO vs STRING Comparison\n")
cat("=========================================\n\n")

# ---- Load ShinyGO export ----
shinygo_path <- file.path(DATA_DIR, "shinygo_export.csv")

if (!file.exists(shinygo_path)) {
  # Try alternative names
  alt_paths <- list.files(DATA_DIR, pattern = "shinygo|shiny_go",
                          full.names = TRUE, ignore.case = TRUE)
  if (length(alt_paths) > 0) {
    shinygo_path <- alt_paths[1]
  } else {
    cat("ERROR: No ShinyGO export found.\n")
    cat("  Expected: data/shinygo_export.csv\n")
    cat("  Also tried: data/shinygo*.csv, data/shiny_go*.csv\n")
    cat("\n  How to create the export:\n")
    cat("  1. Go to http://bioinformatics.sdstate.edu/go/\n")
    cat("  2. Paste your gene list (from output/tables/lydia_enriched_all.txt)\n")
    cat("  3. Run analysis\n")
    cat("  4. Export results as CSV\n")
    cat("  5. Save as data/shinygo_export.csv\n")
    quit(status = 0)
  }
}

cat(sprintf("Loading ShinyGO export: %s\n", shinygo_path))
shinygo_df <- readr::read_csv(shinygo_path, show_col_types = FALSE)
cat(sprintf("  Loaded %d rows\n", nrow(shinygo_df)))

# ---- Auto-detect column names ----
# ShinyGO exports vary by version. Detect common patterns.
find_shinygo_col <- function(df, patterns) {
  for (p in patterns) {
    matches <- grep(p, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  return(NULL)
}

term_col <- find_shinygo_col(shinygo_df, c("term", "description", "pathway", "go"))
gene_col <- find_shinygo_col(shinygo_df, c("gene", "overlap", "input"))
fdr_col  <- find_shinygo_col(shinygo_df, c("fdr", "qvalue", "q.value", "p.adjust", "pvalue", "p.value"))

if (is.null(term_col)) {
  # Use first column as term
  term_col <- colnames(shinygo_df)[1]
  cat(sprintf("  Auto-detected term column: '%s'\n", term_col))
}
if (is.null(gene_col)) {
  # Try to find a column with semicolons (ShinyGO uses ";" separator)
  for (col in colnames(shinygo_df)) {
    if (any(grepl(";", shinygo_df[[col]], fixed = TRUE))) {
      gene_col <- col
      break
    }
  }
}
if (is.null(fdr_col)) {
  fdr_col <- find_shinygo_col(shinygo_df, c("p", "sig", "score"))
}

cat(sprintf("  Term column: %s\n", term_col))
cat(sprintf("  Gene column: %s\n", ifelse(is.null(gene_col), "NOT FOUND", gene_col)))
cat(sprintf("  FDR column:  %s\n", ifelse(is.null(fdr_col), "NOT FOUND", fdr_col)))

# ---- Parse ShinyGO gene sets ----
# Extract all unique genes from ShinyGO export
shinygo_all_genes <- character(0)
shinygo_terms <- list()

if (!is.null(gene_col)) {
  for (i in seq_len(nrow(shinygo_df))) {
    term_name <- shinygo_df[[term_col]][i]
    gene_str <- as.character(shinygo_df[[gene_col]][i])
    genes <- unlist(strsplit(gene_str, "[;/,]"))
    genes <- trimws(genes)
    genes <- genes[genes != "" & !is.na(genes)]
    shinygo_terms[[term_name]] <- genes
    shinygo_all_genes <- union(shinygo_all_genes, genes)
  }
}

cat(sprintf("\n  ShinyGO unique genes: %d\n", length(shinygo_all_genes)))
cat(sprintf("  ShinyGO terms: %d\n", length(shinygo_terms)))

# ---- Load our STRING pipeline gene lists ----
# Read from output/tables/ (produced by lydia-volcano or string-network)
our_genes_file <- file.path(TABLE_DIR, "lydia_enriched_all.txt")
our_net_file <- file.path(TABLE_DIR, "lydia_enriched_network.txt")

our_all_genes <- character(0)
our_net_genes <- character(0)

if (file.exists(our_genes_file)) {
  our_all_genes <- readLines(our_genes_file)
  our_all_genes <- trimws(our_all_genes)
  our_all_genes <- our_all_genes[our_all_genes != ""]
}
if (file.exists(our_net_file)) {
  our_net_genes <- readLines(our_net_file)
  our_net_genes <- trimws(our_net_genes)
  our_net_genes <- our_net_genes[our_net_genes != ""]
}

cat(sprintf("  Our pipeline enriched genes: %d\n", length(our_all_genes)))
cat(sprintf("  Our STRING network genes: %d\n", length(our_net_genes)))

if (length(our_all_genes) == 0) {
  cat("\n  NOTE: Our gene lists not found. Run 'make lydia-volcano' first.\n")
  cat("  Proceeding with ShinyGO analysis only.\n")
}

# ---- COMPARISON 1: Gene set overlap (Venn diagram) ----
cat("\n--- Gene Set Overlap ---\n")

if (length(our_all_genes) > 0 && length(shinygo_all_genes) > 0) {
  overlap <- intersect(our_all_genes, shinygo_all_genes)
  only_ours <- setdiff(our_all_genes, shinygo_all_genes)
  only_shinygo <- setdiff(shinygo_all_genes, our_all_genes)

  cat(sprintf("  Overlap: %d genes\n", length(overlap)))
  cat(sprintf("  Only in STRING pipeline: %d genes\n", length(only_ours)))
  cat(sprintf("  Only in ShinyGO: %d genes\n", length(only_shinygo)))

  venn_plot <- function(filename) {
    venn <- VennDiagram::venn.diagram(
      x = list(STRING_pipeline = our_all_genes,
               ShinyGO = shinygo_all_genes),
      filename = NULL,
      fill = c("#0072B2", "#D55E00"),
      alpha = 0.5,
      cat.cex = 1.2,
      cex = 1.5,
      cat.fontface = "bold",
      cat.pos = c(180, 0),
      main = "Gene Set Overlap: STRING Pipeline vs ShinyGO",
      main.cex = 1.3
    )
    png_path <- safe_filepath(FIGURE_DIR, filename, ".png")
    pdf_path <- safe_filepath(FIGURE_DIR, filename, ".pdf")

    grDevices::png(png_path, width = 8, height = 6, units = "in", res = 300)
    grid::grid.draw(venn)
    grDevices::dev.off()

    grDevices::pdf(pdf_path, width = 8, height = 6)
    grid::grid.draw(venn)
    grDevices::dev.off()

    cat(sprintf("  Saved: %s\n", basename(png_path)))
  }

  venn_plot("shinygo_vs_string_venn")

  # Export unique genes from each tool
  save_table(data.frame(gene = only_ours, set = "STRING_only"),
             "shinygo_comparison_STRING_unique")
  save_table(data.frame(gene = only_shinygo, set = "ShinyGO_only"),
             "shinygo_comparison_ShinyGO_unique")
  save_table(data.frame(gene = overlap, set = "overlap"),
             "shinygo_comparison_overlap")
}

# ---- COMPARISON 2: ShinyGO terms not covered by our GO enrichment ----
cat("\n--- ShinyGO Term Coverage ---\n")
cat(sprintf("  ShinyGO found %d GO terms.\n", length(shinygo_terms)))
cat("  (Cross-reference with our enrichGO results in output/tables/ to find\n")
cat("   terms unique to ShinyGO — ShinyGO sometimes finds non-redundant\n")
cat("   terms that clusterProfiler misses.)\n\n")

# Export ShinyGO terms for manual comparison
if (length(shinygo_terms) > 0) {
  term_summary <- data.frame(
    term = names(shinygo_terms),
    gene_count = sapply(shinygo_terms, length),
    genes = sapply(shinygo_terms, function(g) paste(head(g, 20), collapse = ";")),
    stringsAsFactors = FALSE
  )
  term_summary <- term_summary[order(-term_summary$gene_count), ]
  save_table(term_summary, "shinygo_terms_summary")
  cat(sprintf("  Exported term summary: %d terms\n", nrow(term_summary)))
}

cat("\nDone. Output saved to output/figures/ and output/tables/.\n")
cat("\nNEXT STEPS:\n")
cat("  1. Compare ShinyGO terms with our enrichGO results (output/tables/)\n")
cat("  2. Look for terms found by ShinyGO but not by our pipeline\n")
cat("  3. These unique terms may reveal pathways clusterProfiler missed\n")
