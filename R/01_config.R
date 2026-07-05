###############################################################################
# 01_config.R
# Central configuration — thresholds, paths, parameters, experiment definitions.
# Customized for TRIP4/ASCC proteomics study.
#
# This file is loaded FIRST by every analysis script. It defines:
#   - Statistical thresholds (what counts as "significant")
#   - File paths (where to find data, where to save output)
#   - Experiment definitions (maps short names to real mass spec filenames)
#   - Color schemes for plots
#   - Two key helper functions: load_proteomics_csv() and get_significant_genes()
#
# All variables defined here become GLOBAL — they're available to every
# script that sources this file.
###############################################################################

# ---- Significance thresholds ----
# These define what we consider a "significant" protein.
#
# P_VALUE_CUTOFF: The adjusted p-value (False Discovery Rate) threshold.
#   A protein is "significant" if its adj.P.Val is below this number.
#   0.05 means "we accept a 5% chance that this protein is a false positive."
#   The "adjusted" part means it's already been corrected for multiple testing
#   (we tested ~3000 proteins, so without correction, many would look
#   significant just by chance).
#
# LOG2FC_CUTOFF: The minimum fold change (magnitude) for significance.
#   log2FC of 1.0 = 2-fold change, 2.0 = 4-fold change, 0.5 = ~1.4-fold change.
#   Lydia used 0.5 (more permissive) rather than 1.0 (more strict).
#   We use abs() so both up-regulated (positive) and down-regulated (negative)
#   proteins count.
P_VALUE_CUTOFF    <- 0.05
LOG2FC_CUTOFF     <- 0.5

# ---- File paths ----
# here::here() builds paths relative to the project root.
# The "here" package automatically finds the project root by looking for
# files like .git, .Rproj, or the directory structure.
# This makes the code work regardless of what folder R is running from.
#
# DATA_DIR:   Where the input CSV files live (mass spec output)
# OUTPUT_DIR: Where all generated results go
# FIGURE_DIR: PNG and PDF plots go here
# TABLE_DIR:  CSV result tables go here
DATA_DIR   <- here::here("data")
OUTPUT_DIR <- here::here("output")
FIGURE_DIR <- here::here("output", "figures")
TABLE_DIR  <- here::here("output", "tables")

# dir.create() makes a folder if it doesn't exist.
# recursive = TRUE: also creates parent folders (output/ → output/figures/)
# showWarnings = FALSE: don't warn if the folder already exists
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Organism annotation ----
# ORGDB: The annotation database package. org.Hs.eg.db contains mappings
#   between gene symbols, Entrez IDs, GO terms, etc. for HUMAN genes.
#   "Hs" = Homo sapiens. Other options: org.Mm.eg.db (mouse), org.Dm.eg.db (fly).
KEYTYPE <- "SYMBOL"     # Our gene identifiers are gene symbols (e.g., "TRIP4")
STRING_TAXON <- 9606    # Human NCBI taxon ID (used by STRING database)

# ---- CRAC data column names (different from mass spec CSVs) ----
# CRAC (Cross-linking and Analysis of cDNAs) data uses different column
# names than the DIA-NN/limma pipeline output. The Excel export has:
#   external_gene_name → gene symbol (maps to COL_GENE)
#   logFC              → log fold change (same name, maps to COL_LOG2FC)
#   FDR                → false discovery rate (maps to COL_PADJ)
CRAC_GENE_COL   <- "external_gene_name"
CRAC_LOG2FC_COL <- "logFC"
CRAC_PADJ_COL   <- "FDR"
CRAC_PVAL_COL   <- "PValue"

# ---- Column name mapping (from mass spec output CSV) ----
# These tell the pipeline which columns to use from the CSV files.
# The mass spectrometry pipeline outputs CSVs with specific column names.
# If your column names differ, change them here.
COL_GENE   <- "Gene"       # Column containing gene symbols (e.g., "TRIP4")
COL_LOG2FC <- "logFC"      # Column with log2 fold change values
COL_PADJ   <- "adj.P.Val"  # Column with adjusted p-values (FDR)
COL_PVAL   <- "P.Val"      # Column with raw p-values (less commonly used)

# ---- Experiment definitions ----
# This list maps SHORT names (left side) to the ACTUAL filenames from
# the mass spectrometry output (right side, without the _diffEx_minProb suffix).
#
# In R, list() creates a named list. You access elements by name:
#   EXPERIMENTS$turbo_trip4_vs_wt  →  "BK467_TRIP4_vs_BK467_WT"
#
# The pipeline uses the SHORT names (left) in its logic to identify
# which experiment is which. The ACTUAL names (right) are matched
# against the discovered CSV filenames.
#
# Naming convention:
#   turbo_ = TurboID experiments (proximity labeling in HeLa cells)
#   flag_  = Flag IP experiments (co-immunoprecipitation in HEK293 cells)
#   RA     = Retinoic Acid treated samples
EXPERIMENTS <- list(
  # TurboID (HeLa, BK467) — proximity labeling experiment
  turbo_trip4_vs_wt      = "BK467_TRIP4_vs_BK467_WT",
  turbo_RA_vs_trip4      = "BK467_TRIP4_RA02_vs_BK467_TRIP4",
  turbo_RA_vs_wt         = "BK467_TRIP4_RA02_vs_BK467_WT",
  turbo_cross_467_504    = "BK467_TRIP4_vs_BK504_TRIP4",
  turbo_RA_cross         = "BK467_TRIP4_RA02_vs_BK504_TRIP4_RA04",
  # TurboID (HeLa, BK504) — biological replicate batch
  turbo_RA04_vs_trip4    = "BK504_TRIP4_RA04_vs_BK504_TRIP4",
  turbo_RA04_vs_wt       = "BK504_TRIP4_RA04_vs_BK467_WT",
  # Flag IP (HEK293, BK516) — co-immunoprecipitation experiment
  flag_cflag_vs_ctrl     = "BK516_Cflag_vs_BK516_Ctrl",
  flag_nflag_vs_ctrl     = "BK516_Nflag_vs_BK516_Ctrl",
  flag_cflag_vs_nflag    = "BK516_Cflag_vs_BK516_Nflag",
  # Flag IP + RA (HEK293, BK523) — retinoic acid treatment
  flag_RA_cflag_vs_cflag = "BK523_Cflag_RA04_vs_BK516_Cflag",
  flag_RA_cflag_vs_ctrl  = "BK523_Cflag_RA04_vs_BK523_Ctrl_RA04",
  flag_RA_cflag_vs_nflag = "BK523_Cflag_RA04_vs_BK523_Nflag_RA04",
  flag_RA_nflag_vs_nflag = "BK523_Nflag_RA04_vs_BK516_Nflag",
  flag_RA_nflag_vs_ctrl  = "BK523_Nflag_RA04_vs_BK523_Ctrl_RA04",
  # CHX/DMSO TurboID (TRIP4 translation inhibitor experiments)
  chx_vs_dmso            = "TRIP4_CHX_vs_TRIP4_DMSO",
  chx_vs_wt              = "TRIP4_CHX_vs_WT",
  dmso_vs_wt             = "TRIP4_DMSO_vs_WT"
)

# ---- Key comparisons (used for Venn diagrams, overlay plots, GO) ----
# Not all experiments get the full treatment — these are the 4 "main"
# comparisons that drive the key figures.
MAIN_COMPARISONS <- list(
  turbotrip4  = "turbo_trip4_vs_wt",      # Main TurboID experiment
  flagC       = "flag_cflag_vs_ctrl",     # Main Flag IP experiment
  turbotrip4RA = "turbo_RA_vs_wt",        # TurboID with retinoic acid
  flagRA      = "flag_RA_cflag_vs_cflag"  # Flag IP with retinoic acid
)

# ---- Global color scheme (consistent across ALL plots) ----
# Each category gets the SAME color wherever it appears — this is a
# critical scientific best practice. Okabe-Ito / Wong Color Universal
# Design palette for colorblind safety (Wong B, Nature Methods 2011).
#
# Usage: GLOBAL_COLORS["enriched_up"]  →  "#D55E00"
GLOBAL_COLORS <- c(
  # Main volcano categories
  "ascc_core"     = "#0072B2",   # Blue — ASCC complex
  "known_ia"      = "#009E73",   # Bluish-green — known interactors
  "enriched_up"   = "#D55E00",   # Vermillion/orange — enriched in TRIP4
  "enriched_dn"   = "#B0B0B0",   # Gray — enriched in WT
  "nonsig"        = "#D0D0D0",   # Light gray — not significant

  # Flag IP validation categories (DISTINCT hues, not shades of orange)
  "flag_both"     = "#F0E442",   # Yellow — validated by C-Flag + N-Flag
  "flag_c_only"   = "#56B4E9",   # Sky blue — validated by C-Flag only
  "flag_n_only"   = "#CC79A7",   # Pink — validated by N-Flag only

  # CHX/DMSO specific
  "chx_enriched"  = "#D55E00",   # Orange — enriched in CHX
  "dmso_enriched" = "#7B3294",   # Purple — enriched in DMSO

  # Venn diagram fill colors (solid, not count-based)
  "venn_overlap"  = "#08519C",   # Dark blue — common to both sets
  "venn_a_only"   = "#6BAED6",   # Light blue — set A unique region
  "venn_b_only"   = "#FDAE6B"    # Light salmon — set B unique region
)

# Legacy color scheme (Lydia's style, retained for 02_volcano_plots.R)
CATEGORY_COLORS <- c(
  "ia"           = "#d95f02",   # Known interactors — orange
  "flagMulti"    = "#7570b3",   # Flag IP hit in 2+ conditions — purple
  "flagOnce"     = "#1f78b4",   # Flag IP hit in 1 condition — blue
  "CRAC"         = "#e7298a",   # CRAC RNA hit — pink
  "gp"           = "#7570b3",   # GPATCH family — purple
  "LARPs"        = "#a6761d",   # LARP family — brown
  "dhx"          = "#66a61e",   # DHX helicases — green
  "ddx"          = "#1f78b4",   # DDX helicases — blue
  "high"         = "red",       # High-confidence hits — red
  "inNetwork"    = "#1b9e77",   # STRING network hit — teal
  "TRUE"         = "#1b9e77",   # Significant but uncategorized — teal
  "FALSE"        = "grey60"     # Not significant — grey
)

# Colors for multi-experiment overlay plots
EXPERIMENT_COLORS <- c(
  "TurboID_TRIP4"     = "#E64B35",   # Red
  "TurboID_TRIP4_RA"  = "#E69F00",   # Orange
  "FlagIP_C"          = "#4DBBD5",   # Blue
  "FlagIP_C_RA"       = "#3C5488"    # Dark blue
)

# Default figure dimensions (in inches) and resolution (dots per inch)
FIG_WIDTH  <- 8
FIG_HEIGHT <- 6
FIG_DPI    <- 300   # 300 DPI = publication quality

# ---- Gene families (from Lydia's scripts) ----
# These gene families are highlighted specially in volcano plots.
# They are RNA-binding proteins, relevant to TRIP4's role.
#
# GPATCH: A specific list of genes (hardcoded because there's no simple
#   naming pattern to match them automatically).
# DHX, DDX, LARP: Matched by PREFIX. NULL means "match by regular expression
#   ^DHX, ^DDX, ^LARP later in the code" rather than listing every member.
GENE_FAMILIES <- list(
  GPATCH = c("AGGF1", "CMTR1", "GPATCH1", "GPATCH2", "GPATCH3", "GPATCH4",
             "GPATCH8", "GPATCH11", "PINX1", "SUGP2", "SUGP1", "NKRF",
             "GPKOW", "GPANK1", "RBM17", "RBM10", "RBM5", "RBM6",
             "SON", "ZGPAT", "CHERP", "TFIP11"),
  DHX = NULL,  # Matched by prefix "^DHX" (e.g., DHX9, DHX15, DHX29...)
  DDX = NULL,  # Matched by prefix "^DDX" (e.g., DDX3, DDX5, DDX17...)
  LARP = NULL  # Matched by prefix "^LARP" (e.g., LARP1, LARP4, LARP7...)
)

# ---- GO enrichment parameters ----
# BP = Biological Process (what does the protein DO)
# MF = Molecular Function (what does the protein BIND/CATALYZE)
# CC = Cellular Component (WHERE in the cell is the protein)
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.05
GO_PADJUST_METHOD <- "BH"             # Benjamini-Hochberg FDR correction
GO_ONTOLOGIES <- c("BP", "MF", "CC")  # Run all three domains

# ---- STRING network parameters ----
STRING_VERSION <- "12.0"
STRING_SCORE_THRESHOLD <- 400  # Medium confidence (0-1000 scale)

# ---- Helper function: fuzzy column name matching ----
# Tries a list of candidate column names against what's actually in the CSV.
# Returns the first match found (case-insensitive).
# This handles minor naming variations between different software versions.
#
# function() defines a new function in R.
# The syntax is: function(parameters) { body }
find_column <- function(available_cols, candidates) {
  # First pass: try exact match (case-sensitive)
  for (cand in candidates) {
    # %in% checks if cand exists in the available_cols vector
    if (cand %in% available_cols) return(cand)
  }
  # Second pass: try case-insensitive match
  # tolower() converts text to lowercase for comparison
  for (cand in candidates) {
    # which() returns the INDEX (position) where the condition is TRUE
    match_idx <- which(tolower(available_cols) == tolower(cand))
    if (length(match_idx) > 0) return(available_cols[match_idx[1]])
  }
  # If nothing matched, return NULL (R's version of "nothing")
  return(NULL)
}

# Custom operator: "null coalescing" — if the left side is NULL, use the right.
# This is like the ?? operator in C# or the || operator in JavaScript.
# We use it for readable default values: actual_gene %||% "NOT FOUND"
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Helper function: load and validate a CSV file ----
# This function:
#   1. Reads the CSV file
#   2. Finds the gene, logFC, and padj columns (with fuzzy matching)
#   3. Reports what it found
#   4. If columns are missing, prints a detailed error
#   5. Returns a clean data frame with exactly 3 columns: gene, log2FC, padj
#
# Parameters with = signs are DEFAULT values — if you don't specify them
# when calling the function, these defaults are used.
load_proteomics_csv <- function(filepath,
                                gene_col = COL_GENE,
                                log2fc_col = COL_LOG2FC,
                                padj_col = COL_PADJ) {
  # Check if the file exists before trying to read it
  if (!file.exists(filepath)) {
    # stop() is R's "throw an error" — it stops execution with a message.
    stop(sprintf("File not found: %s", filepath))
  }

  # readr::read_csv() reads a CSV file into a "data frame" (a table, like
  # an Excel sheet). Each column is a variable, each row is a protein.
  # The :: syntax means "use the read_csv function FROM the readr package."
  # show_col_types = FALSE suppresses a harmless message about column types.
  df <- readr::read_csv(filepath, show_col_types = FALSE)
  cols <- colnames(df)   # Get the column names as a character vector

  # Try to find the actual column names using fuzzy matching.
  # c() combines multiple values into a vector — these are all the
  # candidate names to try, in order of preference.
  actual_gene <- find_column(cols, c(gene_col, "Gene", "gene", "Gene.name",
                                      "Gene.symbol", "Gene.names", "SYMBOL",
                                      "Gene.Symbol"))
  actual_log2fc <- find_column(cols, c(log2fc_col, "logFC", "log2FC",
                                        "Log2FC", "log2.fold.change",
                                        "log2FoldChange"))
  actual_padj <- find_column(cols, c(padj_col, "adj.P.Val", "adj.P.value",
                                      "FDR", "padj", "q.value", "adj.PVal"))

  # Print what was matched (so the user can verify correctness)
  cat(sprintf("  Columns: gene='%s', logFC='%s', padj='%s'\n",
              actual_gene %||% "NOT FOUND",
              actual_log2fc %||% "NOT FOUND",
              actual_padj %||% "NOT FOUND"))

  # Build a list of what's missing (if anything)
  # c() with no arguments creates an empty vector that we append to
  missing <- c()
  if (is.null(actual_gene))   missing <- c(missing, sprintf("gene column (expected '%s')", gene_col))
  if (is.null(actual_log2fc)) missing <- c(missing, sprintf("logFC column (expected '%s')", log2fc_col))
  if (is.null(actual_padj))   missing <- c(missing, sprintf("padj column (expected '%s')", padj_col))

  # If any required columns are missing, print a detailed error and stop
  if (length(missing) > 0) {
    cat("\n  ========================================\n")
    cat("  COLUMN ERROR in file:\n")
    cat(sprintf("    %s\n", basename(filepath)))
    cat("  ========================================\n")
    cat("\n  Missing required columns:\n")
    # for loops in R: for (item in collection) { do something with item }
    for (m in missing) cat("    - ", m, "\n", sep = "")
    cat("\n  Available columns in this file:\n")
    for (c in cols) cat("    - ", c, "\n", sep = "")
    cat("\n  Tip: If your column names differ, edit R/01_config.R\n")
    cat("       and update COL_GENE, COL_LOG2FC, COL_PADJ.\n\n")
    stop("Required columns not found in ", basename(filepath))
  }

  # Extract only the 3 columns we need and rename them to standard names
  keep_cols <- c(actual_gene, actual_log2fc, actual_padj)
  rename_map <- c("gene", "log2FC", "padj")    # Standardized names
  names(rename_map) <- keep_cols                # Map old → new

  # df[, keep_cols] selects only the specified columns from the data frame
  # drop = FALSE keeps it as a data frame even if only 1 column is selected
  df <- df[, keep_cols, drop = FALSE]
  colnames(df) <- rename_map[keep_cols]

  # Remove rows where the gene name is missing (NA) or empty string
  # is.na() returns TRUE for missing values
  # & is the logical AND operator
  # != is "not equal to"
  df <- df[!is.na(df$gene) & df$gene != "", ]
  cat(sprintf("  Loaded %s: %d proteins\n", basename(filepath), nrow(df)))
  return(df)
}

# ---- Helper function: extract significant genes ----
# Given a data frame with gene/log2FC/padj columns, returns just the names
# of genes that pass BOTH the p-value AND fold-change thresholds.
#
# This is used everywhere — volcano plots, Venn diagrams, GO enrichment,
# STRING network all need to know "which proteins are significant?"
get_significant_genes <- function(df,
                                  padj_cutoff = P_VALUE_CUTOFF,
                                  log2fc_cutoff = LOG2FC_CUTOFF) {
  # df$gene[condition] selects gene names where the condition is TRUE.
  # The condition: padj < cutoff AND absolute value of log2FC > cutoff.
  # abs() = absolute value (handles both up- and down-regulation).
  sig <- df$gene[df$padj < padj_cutoff & abs(df$log2FC) > log2fc_cutoff]
  # Remove duplicates and missing values
  sig <- unique(sig[!is.na(sig)])
  return(sig)
}

# ---- Print configuration summary ----
# This runs automatically when the file is sourced.
cat("[Config] TRIP4/ASCC proteomics study configuration loaded.\n")
cat(sprintf("  Organism: Human (org.Hs.eg.db, STRING taxon %d)\n", STRING_TAXON))
cat(sprintf("  Significance: adj.P.Val < %.2f AND |logFC| > %.1f\n", P_VALUE_CUTOFF, LOG2FC_CUTOFF))
cat(sprintf("  Experiments defined: %d\n", length(EXPERIMENTS)))
cat("\n")
