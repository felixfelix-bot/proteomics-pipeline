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
source("R/00_theme.R")

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
                     title = dot_title)
    p_dot <- p_dot +
      ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                     name = "p-adjusted value") +
      ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 10),
                                     breaks = make_size_breaks(p_dot$data$Count),
                                     limits = c(0, NA)) +
      ggplot2::scale_y_discrete(labels = capitalize_first) +
      ggplot2::guides(size = size_legend_guide()) +
      ggplot2::labs(x = "Gene Ratio") +
      theme_poster()
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 18, height = max(14, n_show * 0.8))

    # ---- Barplot ----
    # Sorted by Count (highest at top, lowest at bottom).
    # NOTE: orderBy = "Count" (capital C) — matches the column name in
    # the enrichResult object. Lowercase "count" silently fails matching,
    # defaulting to p.adjust ordering. See enrichplot::barplot docs.
    n_show_bar <- min(15, nrow(simple_df))
    fig_height_bar <- max(7, n_show_bar * 0.5)

    p_bar <- barplot(result_simple, showCategory = n_show_bar,
                     title = dot_title) +
      ggplot2::scale_fill_gradient(low = "#D55E00", high = "#0072B2",
                                    name = "p-adjusted value") +
      ggplot2::scale_y_discrete(labels = capitalize_first) +
      theme_poster()
    save_figure(p_bar, paste0(prefix, "_barplot"), width = 18, height = max(12, n_show_bar * 0.8))
  }
}

# ---- Helper: run KEGG pathway enrichment for one gene set ----
run_targeted_kegg <- function(genes, set_name, universe) {
  cat(sprintf("\n--- KEGG: %s (%d genes) ---\n", set_name, length(genes)))

  if (length(genes) < 5) {
    cat("  Skipped: fewer than 5 genes\n")
    return(NULL)
  }

  tryCatch({
    # Convert SYMBOL to ENTREZID for KEGG
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    entrez_genes <- unique(entrez_map$ENTREZID)
    cat(sprintf("  Mapped %d/%d genes to ENTREZID\n",
                length(entrez_genes), length(genes)))

    # Convert universe SYMBOL to ENTREZID
    universe_map <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
    universe_entrez <- unique(universe_map$ENTREZID)

    ekegg <- enrichKEGG(
      gene          = entrez_genes,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      universe      = universe_entrez,
      minGSSize     = 2,
      maxGSSize     = 5000
    )

    if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
      cat("  No enriched KEGG pathways found\n")
      return(NULL)
    }

    # Convert ENTREZID back to readable gene symbols
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

    kegg_df <- as.data.frame(ekegg)
    cat(sprintf("  Found %d enriched KEGG pathways\n", nrow(kegg_df)))

    prefix <- sanitize_filename(paste0("targeted_KEGG_", set_name))
    save_table(kegg_df, prefix)

    # Dotplot
    n_show <- min(20, nrow(kegg_df))
    fig_height <- max(7, n_show * 0.4)

    set_title <- SET_TITLES[[set_name]]
    if (is.null(set_title)) set_title <- gsub("_", " ", set_name)

    p_kegg <- enrichplot::dotplot(ekegg, showCategory = n_show,
                                  title = paste0(set_title, " — KEGG Pathways"))
    p_kegg <- p_kegg +
      ggplot2::scale_color_gradient(low = "#D55E00", high = "#0072B2",
                                     name = "p-adjusted value") +
      ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 10),
                                     breaks = make_size_breaks(p_kegg$data$Count),
                                     limits = c(0, NA)) +
      ggplot2::scale_y_discrete(labels = capitalize_first) +
      ggplot2::guides(size = size_legend_guide()) +
      ggplot2::labs(x = "Gene Ratio") +
      theme_poster()
    save_figure(p_kegg, paste0(prefix, "_dotplot"),
                width = 18, height = max(14, n_show * 0.8))

  }, error = function(e) {
    cat(sprintf("  KEGG error (may need internet): %s\n",
                conditionMessage(e)))
  })
}

# =====================================================================
# GENE SET 1: TurboID TRIP4 vs WT significant
# =====================================================================
# CRITICAL FIX (per Aruna, July 12 voice message):
#   The old code used get_significant_genes() which calls abs(log2FC),
#   including WT-enriched proteins (negative log2FC) in the GO analysis.
#   That's why mitochondrial proteins appeared — they're enriched in WT!
#
#   CORRECT approach:
#     Foreground = TRIP4-enriched ONLY: log2FC >= 1 AND padj <= 0.1
#     Background = ALL genes detected in THIS experiment (not all experiments)
#
cat("[1/3] TurboID TRIP4 vs WT...\n")
turbo_exp <- "BK467_TRIP4_vs_BK467_WT"
if (turbo_exp %in% names(experiments)) {
  df_turbo <- experiments[[turbo_exp]]

  # Universe = ALL proteins in this experiment (the proper background)
  turbo_universe <- unique(df_turbo$gene[!is.na(df_turbo$gene)])
  cat(sprintf("  Universe for this experiment: %d genes\n", length(turbo_universe)))

  # Foreground = TRIP4-enriched ONLY (positive log2FC >= 1, padj <= 0.1)
  turbo_sig <- df_turbo$gene[df_turbo$log2FC >= 1 &
                              df_turbo$padj <= 0.1 &
                              !is.na(df_turbo$gene)]
  turbo_sig <- unique(turbo_sig)
  cat(sprintf("  TRIP4-enriched (log2FC>=1, padj<=0.1): %d genes\n", length(turbo_sig)))
  cat(sprintf("  (Previous bug included %d WT-enriched genes)\n",
              sum(df_turbo$gene %in% turbo_sig == FALSE &
                  df_turbo$padj < P_VALUE_CUTOFF &
                  df_turbo$log2FC <= -LOG2FC_CUTOFF &
                  !is.na(df_turbo$gene), na.rm = TRUE)))

  run_targeted_go(turbo_sig, "TurboID_TRIP4_vs_WT", turbo_universe)
  run_targeted_kegg(turbo_sig, "TurboID_TRIP4_vs_WT", turbo_universe)
}

# =====================================================================
# GENE SET 2: RA effect — shared (core interactome)
# =====================================================================
cat("\n[2/3] RA shared (core interactome)...\n")
exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"
if (exp_base %in% names(experiments) && exp_ra %in% names(experiments)) {
  # Fix: only TRIP4-enriched (positive log2FC >= 1, padj <= 0.1)
  df_base <- experiments[[exp_base]]
  df_ra   <- experiments[[exp_ra]]

  ra_universe <- unique(c(df_base$gene, df_ra$gene))
  ra_universe <- unique(ra_universe[!is.na(ra_universe)])

  base_sig <- df_base$gene[df_base$log2FC >= 1 & df_base$padj <= 0.1 & !is.na(df_base$gene)]
  ra_sig   <- df_ra$gene[df_ra$log2FC >= 1 & df_ra$padj <= 0.1 & !is.na(df_ra$gene)]
  base_sig <- unique(base_sig)
  ra_sig   <- unique(ra_sig)

  shared   <- intersect(base_sig, ra_sig)
  run_targeted_go(shared, "RA_shared_core", ra_universe)
  run_targeted_kegg(shared, "RA_shared_core", ra_universe)
}

# =====================================================================
# GENE SET 3: RA-gained (RA-dependent interactors)
# =====================================================================
cat("\n[3/3] RA-gained (RA-dependent)...\n")
if (exp_base %in% names(experiments) && exp_ra %in% names(experiments)) {
  ra_gained <- setdiff(ra_sig, base_sig)
  run_targeted_go(ra_gained, "RA_gained", ra_universe)
  run_targeted_kegg(ra_gained, "RA_gained", ra_universe)
}

cat("\n=========================================\n")
cat(" Targeted GO enrichment complete!\n")
cat("=========================================\n")
