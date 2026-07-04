###############################################################################
# 01_config.R
# Central configuration — thresholds, paths, and parameters.
# Edit these values to customize the analysis for your data.
###############################################################################

# ---- Significance thresholds ----
P_VALUE_CUTOFF    <- 0.05   # adjusted p-value (FDR) threshold
LOG2FC_CUTOFF     <- 1.0    # |log2 fold change| threshold (1.0 = 2-fold)

# ---- File paths ----
# Place your CSV files in the data/ directory
DATA_DIR   <- here::here("data")
OUTPUT_DIR <- here::here("output")
FIGURE_DIR <- here::here("output", "figures")
TABLE_DIR  <- here::here("output", "tables")

# Create output directories if they don't exist
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Organism annotation ----
ORGDB <- "org.Hs.eg.db"   # Human annotation database
KEYTYPE <- "SYMBOL"        # Input gene IDs are gene symbols (e.g., TP53, BRCA1)

# ---- Column name mapping ----
# Edit these if your CSV columns have different names
COL_GENE   <- "gene"      # Column with gene symbols
COL_LOG2FC <- "log2FC"    # Column with log2 fold change
COL_PADJ   <- "padj"      # Column with adjusted p-value
COL_PVAL   <- "pvalue"    # Column with raw p-value (if available)

# ---- Plot aesthetics ----
# Colors for different experiments on overlay plots
EXPERIMENT_COLORS <- c(
  "WT_vs_POI"        = "#999999",  # grey
  "POI_vs_POIHormone" = "#E69F00",  # orange
  "TurboID"          = "#E64B35",  # red
  "Flag_IP"          = "#4DBBD5",  # blue
  "CRAC_RNA"         = "#00A087"   # teal
)

# Figure dimensions and DPI for publication quality
FIG_WIDTH  <- 8
FIG_HEIGHT <- 6
FIG_DPI    <- 300

# ---- GO enrichment parameters ----
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.05
GO_PADJUST_METHOD <- "BH"     # Benjamini-Hochberg FDR
GO_ONTOLOGIES <- c("BP", "MF", "CC")  # Biological Process, Molecular Function, Cellular Component

# ---- Helper: load and validate a CSV file ----
load_proteomics_csv <- function(filepath,
                                gene_col   = COL_GENE,
                                log2fc_col = COL_LOG2FC,
                                padj_col   = COL_PADJ) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: %s", filepath))
  }

  df <- readr::read_csv(filepath, show_col_types = FALSE)

  # Check required columns exist
  missing <- setdiff(c(gene_col, log2fc_col, padj_col), colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s\nAvailable columns: %s",
      basename(filepath),
      paste(missing, collapse = ", "),
      paste(colnames(df), collapse = ", ")
    ))
  }

  # Select and rename for consistency
  df <- df[, c(gene_col, log2fc_col, padj_col)]
  colnames(df) <- c("gene", "log2FC", "padj")

  # Clean: remove NAs in essential columns
  df <- df[!is.na(df$gene) & df$gene != "", ]

  cat(sprintf("  Loaded %s: %d proteins\n", basename(filepath), nrow(df)))
  return(df)
}

# ---- Helper: extract significant genes ----
get_significant_genes <- function(df,
                                  padj_cutoff = P_VALUE_CUTOFF,
                                  log2fc_cutoff = LOG2FC_CUTOFF) {
  sig <- df$gene[df$padj < padj_cutoff & abs(df$log2FC) > log2fc_cutoff]
  sig <- unique(sig[!is.na(sig)])
  return(sig)
}

cat("[Config] Configuration loaded.\n")
cat(sprintf("  Organism: Human (org.Hs.eg.db)\n"))
cat(sprintf("  Significance: padj < %.2f AND |log2FC| > %.1f\n", P_VALUE_CUTOFF, LOG2FC_CUTOFF))
cat("\n")
