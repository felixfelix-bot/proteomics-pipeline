library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(viridis)
library(dplyr)
library(readr)

## ---- 1. USER SETTINGS -- change these per dataset --------------------------

# Path to working directory
setwd("C:/Users/MS-Auswertung-Pro/Data/Benno_data/Astral Projects data analysis (raw data backuped)/2026-03-27_AZ_BK533_AA_TurboID/analysis_minProb_TRIP4_CHX_vs_TRIP4_HeLa/")

# Path to your differential expression csv
input_csv <- "TRIP4_CHX_vs_TRIP4_HeLa_diffEx_minProb.csv"

# Folder where results (tables + plots) will be written
output_dir <- "GO_enrichment_results"

direction <- "up" #"up", "down", or "both"

use_custom_background <- FALSE # FALSE for default whole proteome, TRUE for only proteins in this experiment

# Column names in your diffex table
col_uniprot   <- "UniProt.ID"
col_logfc     <- "logFC"
col_adjpval   <- "adj.P.Val"
col_is_sig    <- "is_significant"

logfc_cutoff  <- 1
adjp_cutoff   <- 0.05

# Plot settings (used in Section 6)
show <- "FoldEnrichment"   # "FoldEnrichment" or "GeneRatio"
n    <- 20                 # number of top terms to show per plot

# GO settings
ontologies    <- c("BP", "MF", "CC")   # Biological Process, Molecular Function, Cellular Component
pAdjustMethod <- "BH"                  # Benjamini-Hochberg adjusted p-value
pvalueCutoff  <- 0.05
minGSSize     <- 2
maxGSSize     <- 5000
kegg_organism <- "hsa"   # KEGG code for Homo sapiens

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)


## ---- 2. Read and prepare data ----------------------------------------------

diffex <- read_csv(paste0("../",input_csv), show_col_types = FALSE)
diffex[[col_uniprot]] <- sub(";.*", "", diffex[[col_uniprot]])

# Drop rows without a UniProt ID (can't be mapped) and de-duplicate
diffex <- diffex %>%
  filter(!is.na(.data[[col_uniprot]]), .data[[col_uniprot]] != "") %>%
  distinct(.data[[col_uniprot]], .keep_all = TRUE)

# Map UniProt -> Entrez ID and gene SYMBOL (needed for org.Hs.eg.db)
id_map <- bitr(diffex[[col_uniprot]],
               fromType = "UNIPROT",
               toType   = c("ENTREZID", "SYMBOL"),
               OrgDb    = org.Hs.eg.db)

cat(sprintf("Mapped %d / %d UniProt IDs to Entrez IDs (%.1f%%)\n",
            nrow(id_map), nrow(diffex), 100 * nrow(id_map) / nrow(diffex)))

diffex <- diffex %>%
  left_join(id_map, by = setNames("UNIPROT", col_uniprot))

# Custom background/universe = every protein quantified in this experiment
# that was successfully mapped to an Entrez ID. Only used if
# use_custom_background = TRUE (see Section 1); otherwise stays unused and
# enrichGO()/enrichKEGG() fall back to their default whole-proteome background.
universe_entrez <- unique(na.omit(diffex$ENTREZID))

if (use_custom_background) {
  cat(sprintf("Using custom background: %d detected proteins\n", length(universe_entrez)))
} else {
  cat("Using default whole-proteome background\n")
}

## ---- 3. Define the gene set(s) to test, based on `direction` ---------------

stopifnot(direction %in% c("up", "down", "both"))

sig_col <- diffex[[col_adjpval]] < adjp_cutoff & abs(diffex[[col_logfc]]) > logfc_cutoff

up_genes   <- character(0)
down_genes <- character(0)

if (direction %in% c("up", "both")) {
  up_genes <- diffex %>%
    filter(sig_col, .data[[col_logfc]] > 0, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>%
    unique()
  
  cat(sprintf("Up-regulated ('enriched') proteins (mapped): %d\n", length(up_genes)))
  if (length(up_genes) < 5) warning("Very few up-regulated genes - enrichment results may not be meaningful.")
}

if (direction %in% c("down", "both")) {
  down_genes <- diffex %>%
    filter(sig_col, .data[[col_logfc]] < 0, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>%
    unique()
  
  cat(sprintf("Down-regulated ('depleted') proteins (mapped): %d\n", length(down_genes)))
  if (length(down_genes) < 5) warning("Very few down-regulated genes - enrichment results may not be meaningful.")
}

## ---- 4. Run GO and KEGG enrichment ------------------------------------------

# GO enrichment against the default (whole-proteome) background.
# `universe` is intentionally left unspecified.
run_go_enrichment <- function(gene_list, ont, label) {
  if (length(gene_list) < 3) {
    message(sprintf("Skipping %s (%s): fewer than 3 genes.", label, ont))
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = gene_list,
    universe      = if (use_custom_background) universe_entrez else NULL,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff  = pvalueCutoff,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    message(sprintf("No enriched GO terms found for %s (%s).", label, ont))
    return(ego)
  }
  
  # Reduce redundancy among enriched GO terms
  tryCatch(
    simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
}

# KEGG pathway enrichment against the default (whole-genome/proteome)
# background. Also leaves `universe` unspecified for the same reason.
run_kegg_enrichment <- function(gene_list, organism, label) {
  if (length(gene_list) < 3) {
    message(sprintf("Skipping KEGG (%s): fewer than 3 genes.", label))
    return(NULL)
  }
  
  ekegg <- tryCatch(
    enrichKEGG(
      gene          = gene_list,
      universe      = if (use_custom_background) universe_entrez else NULL,
      organism      = organism,
      keyType       = "kegg",       # KEGG uses Entrez IDs for "hsa"
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize
    ),
    error = function(e) {
      message(sprintf("KEGG query failed (%s): %s", label, conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(ekegg) || nrow(ekegg@result) == 0) {
    message(sprintf("No enriched KEGG pathways found for %s.", label))
    return(ekegg)
  }
  
  # Convert Entrez IDs to gene symbols in the result table for readability
  tryCatch(setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
           error = function(e) ekegg)
}

results <- list()

if (length(up_genes) > 0) {
  for (ont in ontologies) {
    results[[paste0("up_", ont)]] <- run_go_enrichment(up_genes, ont, "UP")
  }
  results[["up_KEGG"]] <- run_kegg_enrichment(up_genes, kegg_organism, "UP")
}

if (length(down_genes) > 0) {
  for (ont in ontologies) {
    results[[paste0("down_", ont)]] <- run_go_enrichment(down_genes, ont, "DOWN")
  }
  results[["down_KEGG"]] <- run_kegg_enrichment(down_genes, kegg_organism, "DOWN")
}

## ---- 5. Save result tables --------------------------------------------------

for (name in names(results)) {
  ego <- results[[name]]
  if (!is.null(ego) && nrow(ego@result) > 0) {
    out_path <- file.path(paste0("GO_", name, ".csv"))
    write.csv(ego@result, out_path, row.names = FALSE)
    cat(sprintf("Wrote: %s (%d terms)\n", out_path, nrow(ego@result)))
  }
}

## ---- 6. Plots ---------------------------------------------------------------
## Reads whichever result CSVs exist for the chosen `direction` (missing ones,
## e.g. down_* files when direction = "up", are simply skipped) and produces
## one lollipop plot per ontology/KEGG, per direction present.

read_result_csv <- function(name) {
  path <- paste0("GO_", name, ".csv")
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path)
  if (nrow(df) == 0) return(NULL)
  df
}

bp_up   <- read_result_csv("up_BP")
cc_up   <- read_result_csv("up_CC")
mf_up   <- read_result_csv("up_MF")
kegg_up <- read_result_csv("up_KEGG")

bp_down   <- read_result_csv("down_BP")
cc_down   <- read_result_csv("down_CC")
mf_down   <- read_result_csv("down_MF")
kegg_down <- read_result_csv("down_KEGG")

if (!is.null(kegg_up))   kegg_up   <- kegg_up[kegg_up$p.adjust   < pvalueCutoff, ]
if (!is.null(kegg_down)) kegg_down <- kegg_down[kegg_down$p.adjust < pvalueCutoff, ]

make_lollipop <- function(df, title, filename) {
  
  if (is.null(df) || nrow(df) == 0) {
    message(sprintf("Nothing to plot for '%s' - skipping.", title))
    return(invisible(NULL))
  }
  
  # Transform p-value
  df$p.adjust <- -log10(df$p.adjust)
  
  # Select top terms
  df_top <- df %>%
    arrange(desc(p.adjust)) %>%
    slice_head(n = n)
  
  # Plot
  p <- ggplot(
    df_top,
    aes(
      x = .data[[show]],
      y = reorder(Description, .data[[show]])
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = .data[[show]],
        yend = reorder(Description, .data[[show]]),
        colour = p.adjust
      ),
      linewidth = 1
    ) +
    geom_point(
      aes(
        size = Count,
        colour = p.adjust
      )
    ) +
    scale_colour_viridis_c(
      option = "D",
      direction = -1,
      name = "-log10(adj. p-value)"
    ) +
    scale_size(
      name = "Gene count"
    ) +
    labs(
      x = show,
      y = NULL,
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )
  
  # Save PNG
  ggsave(
    filename = filename,
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

## Called one by one (not in a loop) so you keep full control over each plot's
## title/filename. Any of these is safely skipped (with a message) if that
## direction wasn't run - e.g. all *_down_* calls do nothing when
## direction = "up".

bp_up_p   <- make_lollipop(bp_up,   "Enrichment: Biological Process (UP)", "GO_BP_up_lollipop.png")
cc_up_p   <- make_lollipop(cc_up,   "Enrichment: Cellular Component (UP)", "GO_CC_up_lollipop.png")
mf_up_p   <- make_lollipop(mf_up,   "Enrichment: Molecular Function (UP)", "GO_MF_up_lollipop.png")
kegg_up_p <- make_lollipop(kegg_up, "Enrichment: KEGG pathways (UP)",      "KEGG_up_lollipop.png")

bp_down_p   <- make_lollipop(bp_down,   "Enrichment: Biological Process (DOWN)", "GO_BP_down_lollipop.png")
cc_down_p   <- make_lollipop(cc_down,   "Enrichment: Cellular Component (DOWN)", "GO_CC_down_lollipop.png")
mf_down_p   <- make_lollipop(mf_down,   "Enrichment: Molecular Function (DOWN)", "GO_MF_down_lollipop.png")
kegg_down_p <- make_lollipop(kegg_down, "Enrichment: KEGG pathways (DOWN)",      "KEGG_down_lollipop.png")

# Display (each prints only if it was actually generated)
bp_up_p; cc_up_p; mf_up_p; kegg_up_p
bp_down_p; cc_down_p; mf_down_p; kegg_down_p