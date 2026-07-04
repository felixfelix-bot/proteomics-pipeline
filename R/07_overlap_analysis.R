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
#   source("R/07_overlap_analysis.R"
###############################################################################

# === WHAT DOES THIS SCRIPT DO? ================================================
# This script compares results across DIFFERENT experiments to find proteins
# that show up consistently (or specifically) across methods. Cross-method
# validation is powerful: if a protein is a hit in both TurboID AND Flag IP,
# we're much more confident it's a real TRIP4 interactor.
#
# Three types of overlap analysis:
#   1. Flag IP vs TurboID: Are Flag IP hits also significant in TurboID?
#      (Two different methods measuring protein-protein interactions.)
#   2. RA effect: How does adding retinoic acid (RA) change the interactome?
#      RA is a signaling molecule that may alter TRIP4's binding partners.
#   3. CRAC overlap: CRAC (UV cross-linking and analysis of cDNA) identifies
#      proteins that directly contact RNA. Overlap with our hits suggests
#      RNA-binding interactors.
# =============================================================================

# Print a header to the console so you know which script is running:
cat("\n=========================================\n")
cat(" Cross-Experiment Overlap Analysis\n")
cat("=========================================\n\n")

# ---- Load data ----
# load_all_experiments() is a custom helper (from utils.R) that reads all
# proteomics result files. Returns a named list of data frames.
experiments <- load_all_experiments()

# ---- Load known interactors ----
# "Known interactors" are proteins previously reported to interact with TRIP4
# (from literature or databases). Loading them lets us check if our data
# confirms prior findings.
# file.path() joins path components: DATA_DIR/known_interactors.txt
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
# load_known_interactors() is a custom helper (from utils.R) that reads a
# text file of gene names and returns them as a character vector.
known_interactors <- load_known_interactors(interactors_file)

# =====================================================================
# 1. Flag IP hits overlaid on TurboID volcano
# =====================================================================
# This section takes proteins found significant in Flag IP experiments and
# shows where they land on the TurboID volcano plot. If a Flag IP hit also
# shows a strong TurboID signal, that's strong cross-method validation.
cat("[1/3] Flag IP overlap on TurboID volcano...\n")

# The main TurboID experiment (TRIP4 vs wild-type):
turbo_main <- "turbo_trip4_vs_wt"

# Find all Flag IP experiment names (those starting with "flag_"):
# grep() with value=TRUE returns matching strings.
flag_exp_names <- grep("^flag_", names(experiments), value = TRUE)

# Only proceed if we have both TurboID data AND Flag IP data:
if (turbo_main %in% names(experiments) && length(flag_exp_names) > 0) {
  # Get the TurboID data frame:
  turbo_df <- experiments[[turbo_main]]

  # For each Flag IP experiment, get its list of significant gene names.
  # lapply() applies get_significant_genes() to each Flag IP data frame.
  # get_significant_genes() is a custom helper (from utils.R).
  flag_sig_lists <- lapply(experiments[flag_exp_names], get_significant_genes)

  # unlist() combines all gene lists into one big vector (may have duplicates):
  flag_all_genes <- unlist(flag_sig_lists)

  # table() counts how many times each gene appears across Flag IP experiments.
  # Example: if "TRIP4" is significant in 3 Flag IP experiments, table shows TRIP4=3.
  flag_counts <- table(flag_all_genes)

  # Genes significant in ≥2 Flag IP experiments (more reproducible/robust hits):
  flag_multi <- names(flag_counts)[flag_counts >= 2]
  # Genes significant in only 1 Flag IP experiment (less reproducible):
  flag_once <- names(flag_counts)[flag_counts == 1]

  # ---- Build the volcano plot with Flag IP overlay ----
  toPlot <- turbo_df

  # Classify each protein as significant ("TRUE") or not ("FALSE"):
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )

  # Override categories for Flag IP hits (these take visual priority):
  # Flag IP hit in one experiment → "flagOnce"
  toPlot$category[toPlot$gene %in% flag_once] <- "flagOnce"
  # Flag IP hit in multiple experiments → "flagMulti" (higher confidence)
  toPlot$category[toPlot$gene %in% flag_multi] <- "flagMulti"
  # Known literature interactor → "ia" (interactor annotation)
  toPlot$category[toPlot$gene %in% known_interactors] <- "ia"

  # Convert to factor for consistent coloring (CATEGORY_COLORS from config):
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  # Select proteins to label (Flag IP hits + known interactors):
  label_data <- toPlot[toPlot$category %in% c("ia", "flagMulti", "flagOnce"), ]

  # Build the volcano plot:
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 30, show.legend = FALSE
    ) +
    # Dashed significance threshold lines:
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

  # Save the plot:
  save_figure(p, "overlap_turboid_with_flagip", width = 8, height = 6)

  # Save the gene lists as separate tables for downstream use:
  save_table(data.frame(gene = flag_multi, category = "FlagIP_multi"),
             "flagIP_hits_multi")
  save_table(data.frame(gene = flag_once, category = "FlagIP_once"),
             "flagIP_hits_once")
}

# =====================================================================
# 2. RA-specific changes (retinoic acid effect)
# =====================================================================
# RA (retinoic acid) is a signaling molecule (derivative of vitamin A) that
# regulates gene expression. We compare:
#   - TRIP4 vs WT (baseline): what interacts with TRIP4 normally?
#   - TRIP4+RA vs WT: what interacts with TRIP4 when RA is present?
# Proteins that change between these conditions may be RA-dependent
# interactors — potentially biologically interesting.
cat("\n[2/3] RA-specific change analysis...\n")

# Define the two TurboID experiments to compare:
turbo_base <- "turbo_trip4_vs_wt"    # Baseline (no RA)
turbo_ra <- "turbo_RA_vs_wt"         # With RA treatment

# Only proceed if both experiments exist:
if (turbo_base %in% names(experiments) && turbo_ra %in% names(experiments)) {
  # Get significant gene lists for each condition:
  base_sig <- get_significant_genes(experiments[[turbo_base]])
  ra_sig <- get_significant_genes(experiments[[turbo_ra]])

  # setdiff(A, B) returns elements in A but NOT in B (set difference).
  # RA-specific = significant in RA condition but NOT in baseline.
  # These proteins only interact with TRIP4 when RA is present.
  ra_specific <- setdiff(ra_sig, base_sig)

  # intersect(A, B) returns elements in BOTH A and B.
  # Shared = significant in both conditions (robust interactors regardless of RA).
  shared <- intersect(base_sig, ra_sig)

  # Base-only = significant in baseline but NOT with RA.
  # These interactors are LOST when RA is added — the interaction depends on
  # absence of RA.
  base_only <- setdiff(base_sig, ra_sig)

  # Report the counts:
  cat(sprintf("  Shared (TRIP4 + TRIP4+RA): %d proteins\n", length(shared)))
  cat(sprintf("  RA-specific: %d proteins\n", length(ra_specific)))
  cat(sprintf("  Base-only (lost with RA): %d proteins\n", length(base_only)))

  # Create a Venn diagram to visualize the overlap:
  # A Venn diagram shows overlapping sets as intersecting circles.
  venn_list <- list(
    "TRIP4" = base_sig,
    "TRIP4 + RA" = ra_sig
  )

  # ggVennDiagram::ggVennDiagram() creates a Venn diagram as a ggplot object.
  # label = "count" shows the number of genes in each region.
  p_venn <- ggVennDiagram::ggVennDiagram(venn_list, label = "count", set_size = 5) +
    # scale_fill_gradient() colors the regions from light blue to dark blue:
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("RA Effect: TRIP4 vs TRIP4+RA Interactome") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # Save the Venn diagram:
  save_figure(p_venn, "venn_RA_effect", width = 6, height = 6)

  # Save each gene category as a separate table:
  save_table(data.frame(gene = ra_specific, category = "RA_specific"), "RA_specific_genes")
  save_table(data.frame(gene = shared, category = "shared_TRIP4_RA"), "shared_TRIP4_RA")
  save_table(data.frame(gene = base_only, category = "base_only_lost_with_RA"), "base_only_genes")
}

# ---- Flag IP RA analysis ----
# Same comparison but for Flag IP experiments (a different method).
flag_base <- "flag_cflag_vs_ctrl"       # Flag IP baseline (no RA)
flag_ra <- "flag_RA_cflag_vs_cflag"     # Flag IP with RA

if (flag_base %in% names(experiments) && flag_ra %in% names(experiments)) {
  # Get significant genes for each Flag IP condition:
  flagbase_sig <- get_significant_genes(experiments[[flag_base]])
  flagra_sig <- get_significant_genes(experiments[[flag_ra]])

  # Same set operations as above:
  # RA-specific = in RA condition but not baseline:
  flag_ra_specific <- setdiff(flagra_sig, flagbase_sig)
  # Shared = in both conditions:
  flag_shared <- intersect(flagbase_sig, flagra_sig)

  cat(sprintf("  Flag IP shared (Ctrl + RA): %d proteins\n", length(flag_shared)))
  cat(sprintf("  Flag IP RA-specific: %d proteins\n", length(flag_ra_specific)))

  # Create Venn diagram for Flag IP RA comparison:
  venn_list <- list(
    "Flag IP" = flagbase_sig,
    "Flag IP + RA" = flagra_sig
  )
  p_venn <- ggVennDiagram::ggVennDiagram(venn_list, label = "count", set_size = 5) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("RA Effect: Flag IP vs Flag IP+RA Interactome") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # Save the Venn diagram:
  save_figure(p_venn, "venn_RA_effect_flagip", width = 6, height = 6)
}

# =====================================================================
# 3. CRAC overlap (if CRAC data available)
# =====================================================================
# CRAC = Cross-linking and Analysis of CDNAs. It's a method that uses UV light
# to "freeze" protein-RNA interactions in living cells, then identifies which
# proteins are directly bound to RNA. If our TRIP4 interactors overlap with
# CRAC-identified RNA-binding proteins, it suggests those interactors are
# involved in RNA-related processes.
#
# This section only runs if CRAC data files exist in the data/ directory.
cat("\n[3/3] CRAC RNA overlap (if data available)...\n")

# list.files() searches for files matching a pattern.
# pattern = "*crac*|*CRAC*" matches any filename containing "crac" (case-insensitive).
# full.names = TRUE returns full paths (not just file names).
# ignore.case = TRUE makes the pattern case-insensitive.
crac_file <- list.files(DATA_DIR, pattern = "*crac*|*CRAC*", full.names = TRUE, ignore.case = TRUE)

# Only proceed if at least one CRAC file was found:
if (length(crac_file) > 0) {
  cat("  CRAC data found — running overlap analysis...\n")

  # Process each CRAC file found:
  for (cf in crac_file) {
    cat(sprintf("    Processing: %s\n", basename(cf)))  # basename() strips the directory path

    # ---- Read the CRAC data file ----
    # CRAC files may come in different formats (TSV or CSV), so we try both.
    # tryCatch() is R's error handler:
    #   - First, try readr::read_tsv() (tab-separated values).
    #   - If that fails (error), fall back to readr::read_csv() (comma-separated).
    # This makes the code ROBUST — it won't crash on unexpected file formats.
    # show_col_types = FALSE suppresses the column type messages.
    crac_df <- tryCatch(
      readr::read_tsv(cf, show_col_types = FALSE),
      error = function(e) readr::read_csv(cf, show_col_types = FALSE)
    )
    # Print the column names so you can see the file's structure:
    cat(sprintf("    Columns: %s\n", paste(colnames(crac_df), collapse = ", ")))

    # ---- Find the gene name column ----
    # Different CRAC files may name their gene column differently.
    # grep() searches column names for any of these patterns:
    #   "gene", "Gene", "external_gene", "symbol"
    # [1] takes the first match (in case multiple match).
    gene_col <- grep("gene|Gene|external_gene|symbol", colnames(crac_df), value = TRUE)[1]

    # If no gene column was found, skip this file:
    if (is.na(gene_col)) next

    # Extract unique gene names from the CRAC data:
    # crac_df[[gene_col]] gets the column as a vector.
    # [!is.na(...)] removes missing values.
    # unique() removes duplicates.
    crac_genes <- unique(crac_df[[gene_col]][!is.na(crac_df[[gene_col]])])

    # ---- Overlay CRAC hits on TurboID volcano ----
    # Show which TurboID significant proteins are also CRAC RNA-binding hits.
    if ("turbo_trip4_vs_wt" %in% names(experiments)) {
      turbo_df <- experiments$turbo_trip4_vs_wt

      # Classify proteins as significant or not:
      turbo_df$category <- ifelse(
        turbo_df$padj < P_VALUE_CUTOFF & abs(turbo_df$log2FC) > LOG2FC_CUTOFF,
        "TRUE", "FALSE"
      )
      # Mark CRAC RNA-binding proteins:
      turbo_df$category[turbo_df$gene %in% crac_genes] <- "CRAC"
      # Mark known literature interactors:
      turbo_df$category[turbo_df$gene %in% known_interactors] <- "ia"
      # Convert to factor for coloring:
      turbo_df$category <- factor(turbo_df$category, levels = names(CATEGORY_COLORS))

      # Select CRAC and known-interactor proteins for labeling:
      label_data <- turbo_df[turbo_df$category %in% c("CRAC", "ia"), ]

      # Build the volcano plot:
      p <- ggplot2::ggplot(turbo_df, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
        ggplot2::geom_point(alpha = 0.3, size = 1.2) +
        ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = gene),
          size = 2.5, fontface = "bold",
          max.overlaps = 25, show.legend = FALSE
        ) +
        # Dashed significance threshold lines:
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

      # Save the plot:
      save_figure(p, "overlap_turboid_with_crac", width = 8, height = 6)
    }
  }
} else {
  # If no CRAC files were found, inform the user:
  cat("  No CRAC data files found. Skipping CRAC overlap.\n")
  cat("  To enable: place CRAC data in data/ directory.\n")
}

# Print completion message:
cat("\n=========================================\n")
cat(" Overlap analysis complete!\n")
cat("=========================================\n")
