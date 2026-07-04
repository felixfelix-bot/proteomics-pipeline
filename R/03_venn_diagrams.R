###############################################################################
# 03_venn_diagrams.R
# Generates Venn diagrams comparing protein sets across experiments.
# Extracts overlap and unique sets for downstream GO analysis.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/03_venn_diagrams.R")
###############################################################################

cat("\n=========================================\n")
cat(" Venn Diagram Analysis\n")
cat("=========================================\n\n")

# ---- Load all experiment data ----
experiments <- list()
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  name <- tools::file_path_sans_ext(basename(f))
  if (name == "known_interactors") next
  experiments[[name]] <- load_proteomics_csv(f)
}

# ---- Extract significant gene sets ----
cat("Extracting significant gene sets...\n")
gene_sets <- lapply(experiments, function(df) {
  genes <- get_significant_genes(df)
  cat(sprintf("  %d significant genes\n", length(genes)))
  return(genes)
})

# Give friendly names for display
set_names <- names(gene_sets)
display_names <- gsub("_", " ", set_names)

# =====================================================================
# 1. Venn diagram: TurboID vs Flag IP
# =====================================================================
cat("\n[1/2] TurboID vs Flag IP Venn diagram...\n")

if ("turboid" %in% names(gene_sets) && "flag_ip" %in% names(gene_sets)) {
  venn_list_1 <- list(
    "TurboID"  = gene_sets$turboid,
    "Flag IP"  = gene_sets$flag_ip
  )

  p_venn1 <- ggVennDiagram::ggVennDiagram(
    venn_list_1,
    label = "count",
    set_size = 5,
    label_size = 4
  ) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("TurboID vs Flag Co-IP Overlap") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

  save_figure(p_venn1, "venn_turboid_vs_flagip", width = 6, height = 6)

  # Extract overlaps
  turbo_flag_shared <- intersect(gene_sets$turboid, gene_sets$flag_ip)
  turbo_only <- setdiff(gene_sets$turboid, gene_sets$flag_ip)
  flag_only <- setdiff(gene_sets$flag_ip, gene_sets$turboid)

  cat(sprintf("  Shared (TurboID ∩ Flag IP): %d proteins\n", length(turbo_flag_shared)))
  cat(sprintf("  Unique to TurboID: %d proteins\n", length(turbo_only)))
  cat(sprintf("  Unique to Flag IP: %d proteins\n", length(flag_only)))

  # Save extracted sets
  save_table(data.frame(gene = turbo_flag_shared, category = "TurboID_AND_FlagIP"),
             "overlap_turboid_flagip")
  save_table(data.frame(gene = turbo_only, category = "TurboID_only"),
             "unique_turboid")
  save_table(data.frame(gene = flag_only, category = "FlagIP_only"),
             "unique_flagip")
}

# =====================================================================
# 2. Venn diagram: Protein Interactome vs CRAC RNA Interactome
# =====================================================================
cat("\n[2/2] Protein Interactome vs CRAC RNA Venn diagram...\n")

# For 3+ sets, we can do a combined Venn
if ("crac_rna" %in% names(gene_sets)) {
  # If we have TurboID + Flag IP + CRAC, do a 3-set Venn
  sets_available <- intersect(c("turboid", "flag_ip", "crac_rna"), names(gene_sets))

  if (length(sets_available) >= 3) {
    venn_list_2 <- list(
      "TurboID"   = gene_sets$turboid,
      "Flag IP"   = gene_sets$flag_ip,
      "CRAC RNA"  = gene_sets$crac_rna
    )

    p_venn2 <- ggVennDiagram::ggVennDiagram(
      venn_list_2,
      label = "count",
      set_size = 4,
      label_size = 3.5
    ) +
      ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggplot2::ggtitle("Protein Interactome vs CRAC RNA Interactome") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

    save_figure(p_venn2, "venn_turboid_flagip_crac", width = 7, height = 7)

    # Extract overlaps
    all_three <- Reduce(intersect, gene_sets[sets_available])
    turbo_crac <- intersect(gene_sets$turboid, gene_sets$crac_rna)
    flag_crac <- intersect(gene_sets$flag_ip, gene_sets$crac_rna)

    cat(sprintf("  All three overlap: %d proteins\n", length(all_three)))
    cat(sprintf("  TurboID ∩ CRAC: %d proteins\n", length(turbo_crac)))
    cat(sprintf("  Flag IP ∩ CRAC: %d proteins\n", length(flag_crac)))

    save_table(data.frame(gene = all_three, category = "ALL_THREE_overlap"),
               "overlap_all_three")
    save_table(data.frame(gene = turbo_crac, category = "TurboID_AND_CRAC"),
               "overlap_turboid_crac")
    save_table(data.frame(gene = flag_crac, category = "FlagIP_AND_CRAC"),
               "overlap_flagip_crac")

  } else if (length(sets_available) == 2) {
    # Just do TurboID vs CRAC or whatever two we have
    names_friendly <- gsub("_", " ", sets_available)
    venn_list_2 <- setNames(gene_sets[sets_available], names_friendly)

    p_venn2 <- ggVennDiagram::ggVennDiagram(
      venn_list_2,
      label = "count",
      set_size = 5,
      label_size = 4
    ) +
      ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      ggplot2::ggtitle("Protein Interactome vs CRAC RNA Interactome") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

    save_figure(p_venn2, "venn_protein_vs_rna", width = 6, height = 6)
  }
}

# =====================================================================
# 3. UpSet plot for all experiments (handles >4 sets gracefully)
# =====================================================================
if (length(gene_sets) >= 2) {
  cat("\n[BONUS] Generating UpSet plot for all experiments...\n")

  all_genes <- unique(unlist(gene_sets))
  binary_matrix <- data.frame(gene = all_genes)

  for (name in names(gene_sets)) {
    binary_matrix[[name]] <- as.integer(all_genes %in% gene_sets[[name]])
  }

  # Convert to the format UpSetR expects
  upset_data <- binary_matrix[, -1]  # remove gene column
  rownames(upset_data) <- binary_matrix$gene

  png(file.path(FIGURE_DIR, "upset_all_experiments.png"),
      width = 10, height = 6, units = "in", res = FIG_DPI)
  print(UpSetR::upset(
    upset_data,
    sets = names(gene_sets),
    order.by = "freq",
    main.bar.color = "grey30",
    sets.bar.color = unname(EXPERIMENT_COLORS[names(gene_sets)]),
    point.size = 3.5,
    line.size = 1.2,
    text.scale = 1.1
  ))
  dev.off()
  cat("  Saved: upset_all_experiments.png\n")
}

cat("\n=========================================\n")
cat(" Venn diagrams complete!\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))
cat("=========================================\n")
