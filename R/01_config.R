###############################################################################
# 01_config.R
# Central configuration — thresholds, paths, parameters, experiment definitions.
# Customized for TRIP4/ASCC proteomics study.
###############################################################################

# ---- Significance thresholds ----
P_VALUE_CUTOFF    <- 0.05   # adjusted p-value (FDR) threshold
LOG2FC_CUTOFF     <- 0.5    # |log2 fold change| threshold (Lydia uses 0.5, not 1.0)

# ---- File paths ----
DATA_DIR   <- here::here("data")
OUTPUT_DIR <- here::here("output")
FIGURE_DIR <- here::here("output", "figures")
TABLE_DIR  <- here::here("output", "tables")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Organism annotation ----
ORGDB <- "org.Hs.eg.db"
KEYTYPE <- "SYMBOL"
STRING_TAXON <- 9606   # Human NCBI taxon ID for STRINGdb

# ---- Column name mapping (from mass spec output CSV) ----
COL_GENE   <- "Gene"
COL_LOG2FC <- "logFC"
COL_PADJ   <- "adj.P.Val"
COL_PVAL   <- "P.Val"

# ---- Experiment definitions ----
# Maps a short name to the CSV filename prefix
# The pipeline loads any *.csv in data/ but these names are used for
# Venn diagrams, overlay plots, and output file naming.
EXPERIMENTS <- list(
  # TurboID (HeLa, BK467)
  turbo_trip4_vs_wt      = "BK467_TRIP4_vs_BK467_WT",
  turbo_RA_vs_trip4      = "BK467_TRIP4_RA02_vs_BK467_TRIP4",
  turbo_RA_vs_wt         = "BK467_TRIP4_RA02_vs_BK467_WT",
  # TurboID (HeLa, BK504)
  turbo_RA04_vs_trip4    = "BK504_TRIP4_RA04_vs_BK504_TRIP4",
  turbo_RA04_vs_wt       = "BK504_TRIP4_RA04_vs_BK467_WT",
  # Flag IP (HEK293, BK516)
  flag_cflag_vs_ctrl     = "BK516_Cflag_vs_BK516_Ctrl",
  flag_nflag_vs_ctrl     = "BK516_Nflag_vs_BK516_Ctrl",
  flag_cflag_vs_nflag    = "BK516_Cflag_vs_BK516_Nflag",
  # Flag IP + RA (HEK293, BK523)
  flag_RA_cflag_vs_cflag = "BK523_Cflag_RA04_vs_BK516_Cflag",
  flag_RA_cflag_vs_ctrl  = "BK523_Cflag_RA04_vs_BK523_Ctrl_RA04"
)

# ---- Key comparisons (used for Venn, overlay, GO) ----
MAIN_COMPARISONS <- list(
  turbotrip4  = "turbo_trip4_vs_wt",
  flagC       = "flag_cflag_vs_ctrl",
  turbotrip4RA = "turbo_RA_vs_wt",
  flagRA      = "flag_RA_cflag_vs_cflag"
)

# ---- Plot aesthetics ----
# Lydia-style color scheme for volcano plot categories
CATEGORY_COLORS <- c(
  "ia"           = "#d95f02",   # known interactors — orange
  "flagMulti"    = "#7570b3",   # Flag IP hit in 2+ conditions — purple
  "flagOnce"     = "#1f78b4",   # Flag IP hit in 1 condition — blue
  "CRAC"         = "#e7298a",   # CRAC RNA hit — pink
  "gp"           = "#7570b3",   # GPATCH family — purple
  "LARPs"        = "#a6761d",   # LARP family — brown
  "dhx"          = "#66a61e",   # DHX helicases — green
  "ddx"          = "#1f78b4",   # DDX helicases — blue
  "high"         = "red",       # High-confidence hits — red
  "inNetwork"    = "#1b9e77",   # STRING network hit — teal
  "TRUE"         = "#1b9e77",   # Significant (unlabeled) — teal
  "FALSE"        = "grey60"     # Not significant — grey
)

EXPERIMENT_COLORS <- c(
  "TurboID_TRIP4"     = "#E64B35",
  "TurboID_TRIP4_RA"  = "#E69F00",
  "FlagIP_C"          = "#4DBBD5",
  "FlagIP_C_RA"       = "#3C5488"
)

# Figure dimensions
FIG_WIDTH  <- 8
FIG_HEIGHT <- 6
FIG_DPI    <- 300

# ---- Gene families (from Lydia's scripts) ----
GENE_FAMILIES <- list(
  GPATCH = c("AGGF1", "CMTR1", "GPATCH1", "GPATCH2", "GPATCH3", "GPATCH4",
             "GPATCH8", "GPATCH11", "PINX1", "SUGP2", "SUGP1", "NKRF",
             "GPKOW", "GPANK1", "RBM17", "RBM10", "RBM5", "RBM6",
             "SON", "ZGPAT", "CHERP", "TFIP11"),
  DHX = NULL,  # Matched by prefix "^DHX"
  DDX = NULL,  # Matched by prefix "^DDX"
  LARP = NULL  # Matched by prefix "^LARP"
)

# ---- GO enrichment parameters ----
GO_PVALUE_CUTOFF <- 0.05
GO_QVALUE_CUTOFF <- 0.05
GO_PADJUST_METHOD <- "BH"
GO_ONTOLOGIES <- c("BP", "MF", "CC")

# ---- STRING network parameters ----
STRING_VERSION <- "12.0"
STRING_SCORE_THRESHOLD <- 400  # medium confidence

# ---- Helper: load and validate a CSV file ----
load_proteomics_csv <- function(filepath,
                                gene_col = COL_GENE,
                                log2fc_col = COL_LOG2FC,
                                padj_col = COL_PADJ) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: %s", filepath))
  }

  df <- readr::read_csv(filepath, show_col_types = FALSE)

  missing <- setdiff(c(gene_col, log2fc_col, padj_col), colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s\nAvailable columns: %s",
      basename(filepath),
      paste(missing, collapse = ", "),
      paste(colnames(df), collapse = ", ")
    ))
  }

  keep_cols <- c(gene_col, log2fc_col, padj_col)
  rename_map <- c("gene", "log2FC", "padj")
  names(rename_map) <- keep_cols

  df <- df[, keep_cols, drop = FALSE]
  colnames(df) <- rename_map[keep_cols]

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

cat("[Config] TRIP4/ASCC proteomics study configuration loaded.\n")
cat(sprintf("  Organism: Human (org.Hs.eg.db, STRING taxon %d)\n", STRING_TAXON))
cat(sprintf("  Significance: adj.P.Val < %.2f AND |logFC| > %.1f\n", P_VALUE_CUTOFF, LOG2FC_CUTOFF))
cat(sprintf("  Experiments defined: %d\n", length(EXPERIMENTS)))
cat("\n")
