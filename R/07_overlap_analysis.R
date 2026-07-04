###############################################################################
# 07_overlap_analysis.R
# Cross-experiment overlap analysis:
#   - Flag IP hits overlaid on TurboID volcano
#   - RA-specific changes identified
#   - CRAC RNA interactome overlap (if data available)
#
# Based on Lydia's 20260424_Overlap_TurboID_wOthers.R
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/07_overlap_analysis.R")
###############################################################################

cat("\n=========================================\n")
cat(" Cross-Experiment Overlap Analysis\n")
cat("=========================================\n\n")

# ---- Load data ----
experiments <- list()
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)
for (f in csv_files) {
  name <- tools::file_path_sans_ext(basename(f))
  if (name == "known_interactors") next
  experiments[[name]] <- load_proteomics_csv(f)
}

# ---- Load known interactors ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

# =====================================================================
# 1. Flag IP hits overlaid on TurboID volcano
# =====================================================================
cat("[1/3] Flag IP overlap on TurboID volcano...\n")

turbo_main <- "turbo_trip4_vs_wt"
flag_exp_names <- grep("^flag_", names(experiments), value = TRUE)

if (turbo_main %in% names(experiments) && length(flag_exp_names) > 0) {
  turbo_df <- experiments[[turbo_main]]

  # Build Flag IP hit categories
  flag_sig_lists <- lapply(experiments[flag_exp_names], get_significant_genes)
  flag_all_genes <- unlist(flag_sig_lists)
  flag_counts <- table(flag_all_genes)
  flag_multi <- names(flag_counts)[flag_counts >= 2]
  flag_once <- names(flag_counts)[flag_counts == 1]

  toPlot <- turbo_df
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )
  toPlot$category[toPlot$gene %in% flag_once] <- "flagOnce"
  toPlot$category[toPlot$gene %in% flag_multi] <- "flagMulti"
  toPlot$category[toPlot$gene %in% known_interactors] <- "ia"
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  label_data <- toPlot[toPlot$category %in% c("ia", "flagMulti", "flagOnce"), ]

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 30, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = "TurboID (TRIP4 vs WT) with Flag IP Hits Overlaid",
      color = "Category"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_figure(p, "overlap_turboid_with_flagip", width = 8, height = 6)

  # Save extracted gene lists
  save_table(data.frame(gene = flag_multi, category = "FlagIP_multi"),
             "flagIP_hits_multi")
  save_table(data.frame(gene = flag_once, category = "FlagIP_once"),
             "flagIP_hits_once")
}

# =====================================================================
# 2. RA-specific changes (retinoic acid effect)
# =====================================================================
cat("\n[2/3] RA-specific change analysis...\n")

# TurboID: compare TRIP4 vs WT (base) vs TRIP4+RA vs WT
turbo_base <- "turbo_trip4_vs_wt"
turbo_ra <- "turbo_RA_vs_wt"

if (turbo_base %in% names(experiments) && turbo_ra %in% names(experiments)) {
  base_sig <- get_significant_genes(experiments[[turbo_base]])
  ra_sig <- get_significant_genes(experiments[[turbo_ra]])

  # RA-specific = significant in RA but not in base
  ra_specific <- setdiff(ra_sig, base_sig)
  # Shared = significant in both
  shared <- intersect(base_sig, ra_sig)
  # Base-only = significant in base but lost with RA
  base_only <- setdiff(base_sig, ra_sig)

  cat(sprintf("  Shared (TRIP4 + TRIP4+RA): %d proteins\n", length(shared)))
  cat(sprintf("  RA-specific: %d proteins\n", length(ra_specific)))
  cat(sprintf("  Base-only (lost with RA): %d proteins\n", length(base_only)))

  # Venn diagram
  venn_list <- list(
    "TRIP4" = base_sig,
    "TRIP4 + RA" = ra_sig
  )
  p_venn <- ggVennDiagram::ggVennDiagram(venn_list, label = "count", set_size = 5) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("RA Effect: TRIP4 vs TRIP4+RA Interactome") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_figure(p_venn, "venn_RA_effect", width = 6, height = 6)

  # Save gene lists
  save_table(data.frame(gene = ra_specific, category = "RA_specific"), "RA_specific_genes")
  save_table(data.frame(gene = shared, category = "shared_TRIP4_RA"), "shared_TRIP4_RA")
  save_table(data.frame(gene = base_only, category = "base_only_lost_with_RA"), "base_only_genes")
}

# Flag IP RA analysis
flag_base <- "flag_cflag_vs_ctrl"
flag_ra <- "flag_RA_cflag_vs_cflag"

if (flag_base %in% names(experiments) && flag_ra %in% names(experiments)) {
  flagbase_sig <- get_significant_genes(experiments[[flag_base]])
  flagra_sig <- get_significant_genes(experiments[[flag_ra]])

  flag_ra_specific <- setdiff(flagra_sig, flagbase_sig)
  flag_shared <- intersect(flagbase_sig, flagra_sig)

  cat(sprintf("  Flag IP shared (Ctrl + RA): %d proteins\n", length(flag_shared)))
  cat(sprintf("  Flag IP RA-specific: %d proteins\n", length(flag_ra_specific)))

  venn_list <- list(
    "Flag IP" = flagbase_sig,
    "Flag IP + RA" = flagra_sig
  )
  p_venn <- ggVennDiagram::ggVennDiagram(venn_list, label = "count", set_size = 5) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("RA Effect: Flag IP vs Flag IP+RA Interactome") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_figure(p_venn, "venn_RA_effect_flagip", width = 6, height = 6)
}

# =====================================================================
# 3. CRAC overlap (if CRAC data available)
# =====================================================================
cat("\n[3/3] CRAC RNA overlap (if data available)...\n")

crac_file <- list.files(DATA_DIR, pattern = "*crac*|*CRAC*", full.names = TRUE, ignore.case = TRUE)
if (length(crac_file) > 0) {
  cat("  CRAC data found — running overlap analysis...\n")

  # Read CRAC data (expect gene names + significance info)
  for (cf in crac_file) {
    cat(sprintf("    Processing: %s\n", basename(cf)))
    # CRAC files may have different formats — this is a best-effort read
    crac_df <- tryCatch(
      readr::read_tsv(cf, show_col_types = FALSE),
      error = function(e) readr::read_csv(cf, show_col_types = FALSE)
    )
    cat(sprintf("    Columns: %s\n", paste(colnames(crac_df), collapse = ", ")))

    # Try to find gene name column
    gene_col <- grep("gene|Gene|external_gene|symbol", colnames(crac_df), value = TRUE)[1]
    if (is.na(gene_col)) next

    crac_genes <- unique(crac_df[[gene_col]][!is.na(crac_df[[gene_col]])])

    # Overlay CRAC hits on TurboID volcano
    if ("turbo_trip4_vs_wt" %in% names(experiments)) {
      turbo_df <- experiments$turbo_trip4_vs_wt
      turbo_df$category <- ifelse(
        turbo_df$padj < P_VALUE_CUTOFF & abs(turbo_df$log2FC) > LOG2FC_CUTOFF,
        "TRUE", "FALSE"
      )
      turbo_df$category[turbo_df$gene %in% crac_genes] <- "CRAC"
      turbo_df$category[turbo_df$gene %in% known_interactors] <- "ia"
      turbo_df$category <- factor(turbo_df$category, levels = names(CATEGORY_COLORS))

      label_data <- turbo_df[turbo_df$category %in% c("CRAC", "ia"), ]

      p <- ggplot2::ggplot(turbo_df, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
        ggplot2::geom_point(alpha = 0.3, size = 1.2) +
        ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = gene),
          size = 2.5, fontface = "bold",
          max.overlaps = 25, show.legend = FALSE
        ) +
        ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
        ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
        ggplot2::labs(
          x = expression(Log[2]~Fold~Change),
          y = expression(-Log[10]~(adj.~italic(p)~value)),
          title = "TurboID with CRAC RNA Hits Overlaid",
          color = "Category"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

      save_figure(p, "overlap_turboid_with_crac", width = 8, height = 6)
    }
  }
} else {
  cat("  No CRAC data files found. Skipping CRAC overlap.\n")
  cat("  To enable: place CRAC data in data/ directory.\n")
}

cat("\n=========================================\n")
cat(" Overlap analysis complete!\n")
cat("=========================================\n")
