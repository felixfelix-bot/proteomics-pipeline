###############################################################################
# 10_targeted_venns.R
# Two specific Venn diagrams per researcher request.
#
# VENN 1: TRIP4 without RA vs TRIP4 with RA (RA effect in TurboID)
#   - Set A: Significant proteins in BK467_TRIP4_vs_BK467_WT (TRIP4, no RA)
#   - Set B: Significant proteins in BK467_TRIP4_RA02_vs_BK467_WT (TRIP4 + RA)
#   - Shows which interactions are GAINED or LOST upon retinoic acid treatment
#
# VENN 2: TurboID vs C-Flag vs N-Flag (cross-method comparison)
#   - Set A: TurboID significant (BK467_TRIP4_vs_BK467_WT)
#   - Set B: C-Flag significant (BK516_Cflag_vs_BK516_Ctrl)
#   - Set C: N-Flag significant (BK516_Nflag_vs_BK516_Ctrl)
#   - Shows which proteins are found by multiple independent methods
#
# Usage:
#   make targeted-venn
###############################################################################
cat("\n=========================================\n")
cat(" Targeted Venn Diagrams\n")
cat("=========================================\n\n")

experiments <- load_all_experiments()

# ---- Extract significant gene sets ----
cat("Extracting significant gene sets...\n")
gene_sets <- lapply(experiments, function(df) {
  get_significant_genes(df)
})

# =====================================================================
# VENN 1: TRIP4 without RA vs TRIP4 with RA
# =====================================================================
cat("\n[1/2] Venn: TRIP4 -RA vs TRIP4 +RA (BK467)...\n")

exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"

if (exp_base %in% names(gene_sets) && exp_ra %in% names(gene_sets)) {
  set_a <- gene_sets[[exp_base]]
  set_b <- gene_sets[[exp_ra]]

  cat(sprintf("  TRIP4 -RA: %d significant proteins\n", length(set_a)))
  cat(sprintf("  TRIP4 +RA: %d significant proteins\n", length(set_b)))

  venn1 <- list(
    "TRIP4 (-RA)" = set_a,
    "TRIP4 (+RA)" = set_b
  )

  p1 <- ggVennDiagram::ggVennDiagram(venn1, label = "count",
                                     set_size = 5, label_size = 5) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("RA Effect: TRIP4 without vs with Retinoic Acid (BK467)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12))

  save_figure(p1, "targeted_venn_RA_effect_BK467", width = 6, height = 6)

  # Extract and save the sets
  shared <- intersect(set_a, set_b)
  ra_gained <- setdiff(set_b, set_a)
  ra_lost <- setdiff(set_a, set_b)

  cat(sprintf("  Shared (both): %d proteins\n", length(shared)))
  cat(sprintf("  RA-gained: %d proteins\n", length(ra_gained)))
  cat(sprintf("  RA-lost: %d proteins\n", length(ra_lost)))

  save_table(data.frame(gene = shared, category = "RA_shared"),
             "RA_effect_shared")
  save_table(data.frame(gene = ra_gained, category = "RA_gained"),
             "RA_effect_gained")
  save_table(data.frame(gene = ra_lost, category = "RA_lost"),
             "RA_effect_lost")
} else {
  cat("  WARNING: Missing experiments for RA effect Venn.\n")
}

# =====================================================================
# VENN 2: TurboID vs C-Flag vs N-Flag
# =====================================================================
cat("\n[2/2] Venn: TurboID vs C-Flag vs N-Flag...\n")

exp_turbo <- "BK467_TRIP4_vs_BK467_WT"
exp_cflag <- "BK516_Cflag_vs_BK516_Ctrl"
exp_nflag <- "BK516_Nflag_vs_BK516_Ctrl"

sets_found <- c()
if (exp_turbo %in% names(gene_sets)) sets_found <- c(sets_found, exp_turbo)
if (exp_cflag %in% names(gene_sets)) sets_found <- c(sets_found, exp_cflag)
if (exp_nflag %in% names(gene_sets)) sets_found <- c(sets_found, exp_nflag)

if (length(sets_found) >= 2) {
  venn2 <- list()
  if (exp_turbo %in% sets_found) venn2[["TurboID"]] <- gene_sets[[exp_turbo]]
  if (exp_cflag %in% sets_found) venn2[["C-Flag IP"]] <- gene_sets[[exp_cflag]]
  if (exp_nflag %in% sets_found) venn2[["N-Flag IP"]] <- gene_sets[[exp_nflag]]

  for (n in names(venn2)) {
    cat(sprintf("  %s: %d significant proteins\n", n, length(venn2[[n]])))
  }

  p2 <- ggVennDiagram::ggVennDiagram(venn2, label = "count",
                                     set_size = 4, label_size = 4) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("TurboID vs C-Flag vs N-Flag IP Overlap") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12))

  save_figure(p2, "targeted_venn_turboid_cflag_nflag", width = 7, height = 7)

  # If all 3 sets available, extract the triple overlap
  if (length(sets_found) == 3) {
    all_three <- Reduce(intersect, venn2)
    cat(sprintf("  In all three methods: %d proteins\n", length(all_three)))
    save_table(data.frame(gene = all_three, category = "all_three_methods"),
               "overlap_turboid_cflag_nflag")
  }
} else {
  cat("  WARNING: Need at least 2 of TurboID/C-Flag/N-Flag experiments.\n")
}

cat("\n=========================================\n")
cat(" Targeted Venn diagrams complete!\n")
cat("=========================================\n")
