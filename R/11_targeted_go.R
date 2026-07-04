###############################################################################
# 11_targeted_go.R
# GO enrichment analysis on specific gene sets.
#
# CHANGES:
#   - Barplots sorted by count (high to low) via orderBy
#   - Taller figures to prevent label overlap
#   - Shorter labels via label wrap
#   - Fewer categories in barplots (15 instead of 20) for readability
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

# ---- Helper: run GO enrichment for one gene set ----
run_targeted_go <- function(genes, set_name, universe) {
  cat(sprintf("\n--- %s (%d genes) ---\n", set_name, length(genes)))

  if (length(genes) < 5) {
    cat("  Skipped: fewer than 5 genes\n")
    return(NULL)
  }

  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  [%s] Running enrichGO...\n", ont))

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

    prefix <- sanitize_filename(paste0("targeted_GO_", set_name, "_", ont))

    save_table(res_df, prefix)

    # Dotplot (taller figure for readability)
    n_show <- min(20, nrow(res_df))
    fig_height <- max(7, n_show * 0.4)  # Scale height with number of terms

    p_dot <- dotplot(result, showCategory = n_show,
                     title = paste(set_name, "-", ont)) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 7)
      )
    save_figure(p_dot, paste0(prefix, "_dotplot"), width = 10, height = fig_height)

    # Barplot — sorted by count, taller figure
    n_show_bar <- min(15, nrow(res_df))
    fig_height_bar <- max(7, n_show_bar * 0.5)

    p_bar <- barplot(result, showCategory = n_show_bar,
                     orderBy = "count",
                     title = paste(set_name, "-", ont)) +
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
