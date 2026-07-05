###############################################################################
# 11_targeted_go.R
# GO enrichment analysis on specific gene sets.
#
# THREE GENE SETS:
#   1) TurboID TRIP4 vs WT significant
#   2) RA shared core (sig in both -RA and +RA)
#   3) RA-gained (sig in +RA but not -RA)
#
# For each gene set, runs enrichGO for all three ontologies:
#   BP = Biological Process — what biological pathway the proteins participate in
#   MF = Molecular Function — what molecular activity the proteins have
#   CC = Cellular Component — where in the cell the proteins are located
#
# OUTPUT:
#   - Dotplot: sorted by GeneRatio (x-axis), dot size = count, color = padj
#   - Barplot: sorted by Count (high at top), color = padj
#   - CSV table of all enriched terms
#
# Changes from previous version:
#   - Barplots sorted by Count (capital C, highest at top)
#   - Titles are descriptive (no underscores), include ontology name
#   - Redundant GO terms (e.g., multiple ribosome terms) are collapsed
#     via clusterProfiler::simplify() with cutoff=0.7
#   - "GeneRatio" label fixed to "Gene Ratio" with space
#
# Usage:
#   make targeted-go
###############################################################################
cat("\n=========================================\n")
cat(" Targeted GO Enrichment Analysis\n")
cat("=========================================\n\n")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

experiments <- load_all_experiments()

universe <- unique(unlist(lapply(experiments, function(df) df$gene)))
cat(sprintf("Background universe: %d unique proteins\n\n", length(universe)))

# ---- Human-readable names for each GO ontology ----
ONT_LABELS <- c(
  "BP" = "Biological Process",
  "MF" = "Molecular Function",
  "CC" = "Cellular Component"
)

# ---- Human-readable names for experiment sets (NO underscores) ----
SET_TITLES <- list(
  "TurboID_TRIP4_vs_WT" = "GO analysis of TRIP4 TurboID vs Wild Type",
  "RA_shared_core"      = "GO analysis of RA Shared Core Interactome",
  "RA_gained"           = "GO analysis of RA-Dependent Interactors"
)

# ---- Helper: run GO enrichment for one gene set ----
run_targeted_go <- function(genes, set_name, universe) {
  cat(sprintf("\n--- %s (%d genes) ---\n", set_name, length(genes)))

  if (length(genes) < 5) {
    cat("  Skipped: fewer than 5 genes\n")
    return(NULL)
  }

  # Get human-readable title for this set
  set_title <- SET_TITLES[[set_name]]
  if (is.null(set_title)) {
    # Fallback: replace underscores with spaces
    set_title <- gsub("_", " ", set_name)
  }

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  [%s] %s — Running enrichGO...\n", ont, ONT_LABELS[ont]))

    result <- tryCatch({
      enrichGO(
        gene          = genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        universe      = universe
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      return(NULL)
    })

    if (is.null(result) || nrow(as.data.frame(result)) == 0) {
      cat("    No enriched terms found\n")
      next
    }

    res_df <- as.data.frame(result)
    cat(sprintf("    Found %d enriched GO terms\n", nrow(res_df)))

    # ---- Simplify: group redundant GO terms (e.g., multiple ribosome terms) ----
    # clusterProfiler::simplify() uses GO graph structure to collapse
    # similar child terms into parent terms. cutoff=0.7 is the standard
    # threshold (terms with >0.7 semantic similarity are collapsed).
    result_simple <- tryCatch({
      simplify(result, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }, error = function(e) {
      cat(sprintf("    simplify() note: %s (using full result)\n",
                  conditionMessage(e)))
      result  # Fallback to unsimplified
    })

    simple_df <- as.data.frame(result_simple)
    cat(sprintf("    After simplification: %d terms (was %d)\n",
                nrow(simple_df), nrow(res_df)))

    # Build a filename-safe prefix (underscores are OK in filenames)
    prefix <- sanitize_filename(paste0("targeted_GO_", set_name, "_", ont))

    save_table(simple_df, prefix)

    # ---- Dotplot ----
    # Sorted by GeneRatio (default — user confirmed this looks good)
    n_show <- min(20, nrow(simple_df))
    fig_height <- max(7, n_show * 0.4)

    # Build descriptive title: "GO analysis of TRIP4 TurboID vs Wild Type — Biological Process"
    dot_title <- paste0(set_title, " — ", ONT_LABELS[ont])

    p_dot <- dotplot(result_simple, showCategory = n_show,
                     title = dot_title) +
      ggplot2::labs(x = "Gene Ratio") +  # Fix: space between Gene and Ratio
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_height)

    # ---- Barplot ----
    # Sorted by Count (highest at top, lowest at bottom).
    # NOTE: orderBy = "Count" (capital C) — matches the column name in
    # the enrichResult object. Lowercase "count" silently fails matching,
    # defaulting to p.adjust ordering. See enrichplot::barplot docs.
    n_show_bar <- min(15, nrow(simple_df))
    fig_height_bar <- max(7, n_show_bar * 0.5)

    p_bar <- barplot(result_simple, showCategory = n_show_bar,
                     orderBy = "Count",  # Capital C — sorts by gene count
                     title = dot_title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_bar, paste0(prefix, "_barplot"), width = 10, height = fig_height_bar)
  }
}

# =====================================================================
# GENE SET 1: TurboID TRIP4 vs WT significant
# =====================================================================
cat("[1/3] TurboID TRIP4 vs WT...\n")
turbo_exp <- "BK467_TRIP4_vs_BK467_WT"
if (turbo_exp %in% names(experiments)) {
  turbo_sig <- get_significant_genes(experiments[[turbo_exp]])
  run_targeted_go(turbo_sig, "TurboID_TRIP4_vs_WT", universe)
}

# =====================================================================
# GENE SET 2: RA effect — shared (core interactome)
# =====================================================================
cat("\n[2/3] RA shared (core interactome)...\n")
exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"
if (exp_base %in% names(experiments) && exp_ra %in% names(experiments)) {
  base_sig <- get_significant_genes(experiments[[exp_base]])
  ra_sig   <- get_significant_genes(experiments[[exp_ra]])
  shared   <- intersect(base_sig, ra_sig)
  run_targeted_go(shared, "RA_shared_core", universe)
}

# =====================================================================
# GENE SET 3: RA-gained (RA-dependent interactors)
# =====================================================================
cat("\n[3/3] RA-gained (RA-dependent)...\n")
if (exp_base %in% names(experiments) && exp_ra %in% names(experiments)) {
  ra_gained <- setdiff(ra_sig, base_sig)
  run_targeted_go(ra_gained, "RA_gained", universe)
}

cat("\n=========================================\n")
cat(" Targeted GO enrichment complete!\n")
cat("=========================================\n")
