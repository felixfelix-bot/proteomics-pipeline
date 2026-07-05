###############################################################################
# 10_targeted_venns.R
# Two specific Venn diagrams per researcher request.
#
# VENN 1: TRIP4 without RA vs TRIP4 with RA
#   - Labels: "- RA" and "+ RA" (since title already says TRIP4)
#   - Labels positioned further left to avoid overlapping circles
#   - Extra margin to prevent PDF clipping
#
# VENN 2: TurboID vs Flag IP (2-way, not 3-way)
#   - Combines C-Flag + N-Flag into one "Flag IP" set (union)
#   - No count legend — manual fill colors
#   - Overlap = dark blue, non-overlap = same lighter color
#   - Title: "TurboID vs Flag IP Overlap"
#
# Usage:
#   make targeted-venn
###############################################################################
cat("\n=========================================\n")
cat(" Targeted Venn Diagrams\n")
cat("=========================================\n\n")

experiments <- load_all_experiments()

cat("Extracting significant gene sets...\n")
gene_sets <- lapply(experiments, function(df) {
  get_significant_genes(df)
})

# =====================================================================
# VENN 1: TRIP4 without RA vs TRIP4 with RA
# =====================================================================
cat("\n[1/2] Venn: TRIP4 -RA vs TRIP4 +RA...\n")

exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"

if (exp_base %in% names(gene_sets) && exp_ra %in% names(gene_sets)) {
  set_a <- gene_sets[[exp_base]]
  set_b <- gene_sets[[exp_ra]]

  cat(sprintf("  -RA: %d significant proteins\n", length(set_a)))
  cat(sprintf("  +RA: %d significant proteins\n", length(set_b)))

  # Create Venn with shortened labels
  # "TRIP4" is in the title, so just "- RA" and "+ RA"
  venn1 <- list(
    "- RA" = set_a,
    "+ RA" = set_b
  )

  # Build Venn diagram with manual positioning of set labels
  # set_label_x/y: position labels away from the circles
  p1 <- ggVennDiagram::ggVennDiagram(
    venn1,
    label = "count",
    set_size = 5,
    label_size = 5,
    label_alpha = 0,
    # Move set labels to the left, away from circles
    set_label_x = c(-0.35, -0.35),  # Both labels shifted left
    set_label_y = c(0.25, -0.25)    # One above, one below center
  ) +
    ggplot2::scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggplot2::ggtitle("TRIP4 without vs with Retinoic Acid") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      # Ensure the plot has enough space around it to prevent label clipping
      plot.margin = ggplot2::margin(10, 10, 10, 20, "pt")
    )

  save_figure(p1, "targeted_venn_RA_effect_BK467",
              width = 7, height = 6)  # Wider to prevent clipping

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
# VENN 2: TurboID vs Flag IP (2-way — C-Flag and N-Flag combined)
# =====================================================================
cat("\n[2/2] Venn: TurboID vs Flag IP...\n")

exp_turbo <- "BK467_TRIP4_vs_BK467_WT"
exp_cflag <- "BK516_Cflag_vs_BK516_Ctrl"
exp_nflag <- "BK516_Nflag_vs_BK516_Ctrl"

if (exp_turbo %in% names(gene_sets)) {
  turbo_sig <- gene_sets[[exp_turbo]]
  cat(sprintf("  TurboID: %d significant proteins\n", length(turbo_sig)))

  # Combine C-Flag and N-Flag into one Flag IP set (union)
  flag_sig <- character(0)
  if (exp_cflag %in% names(gene_sets)) {
    flag_sig <- union(flag_sig, gene_sets[[exp_cflag]])
  }
  if (exp_nflag %in% names(gene_sets)) {
    flag_sig <- union(flag_sig, gene_sets[[exp_nflag]])
  }

  if (length(flag_sig) > 0) {
    cat(sprintf("  Flag IP (C-Flag + N-Flag combined): %d significant proteins\n",
                length(flag_sig)))

    venn2 <- list(
      "TurboID" = turbo_sig,
      "Flag IP" = flag_sig
    )

    # 2-way Venn with manual colors:
    #   - Overlap (intersection) = dark blue #0072B2
    #   - Non-overlapping parts = same lighter blue #B3CDE3
    #   - No count legend (no gradient fill)
    p2 <- ggVennDiagram::ggVennDiagram(
      venn2,
      label = "count",
      set_size = 5,
      label_size = 5,
      label_alpha = 0,
      # Position labels left of circles to avoid clipping
      set_label_x = c(-0.30, -0.30),
      set_label_y = c(0.20, -0.20)
    ) +
      # Manual 2-color scheme: non-overlap = light blue, overlap = dark blue
      ggplot2::scale_fill_gradientn(
        colors = c("#B3CDE3", "#0072B2"),
        guide = "none"  # Remove the gradient legend entirely
      ) +
      ggplot2::ggtitle("TurboID vs Flag IP Overlap") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        plot.margin = ggplot2::margin(10, 10, 10, 20, "pt")
      )

    save_figure(p2, "targeted_venn_turboid_flagip",
                width = 6, height = 6)

    # Also save the overlap gene list
    turbo_flag_common <- intersect(turbo_sig, flag_sig)
    cat(sprintf("  Common to both: %d proteins\n", length(turbo_flag_common)))
    save_table(data.frame(gene = turbo_flag_common,
                          category = "TurboID_and_FlagIP"),
               "overlap_turboid_flagip")

  } else {
    cat("  WARNING: No Flag IP data found.\n")
  }
} else {
  cat("  WARNING: TurboID experiment not found.\n")
}

cat("\n=========================================\n")
cat(" Targeted Venn diagrams complete!\n")
cat("=========================================\n")
