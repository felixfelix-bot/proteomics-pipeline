###############################################################################
# 10_targeted_venns.R
# Two specific Venn diagrams per researcher request.
#
# VENN 1: TRIP4 without RA vs TRIP4 with RA
# VENN 2: TurboID vs Flag IP (2-way — C-Flag and N-Flag combined)
#
# Uses VennDiagram::draw.pairwise.venn for clean, publication-quality
# 2-set Venns with manual fill colors and alpha-blended overlap.
# Colors from GLOBAL_COLORS for consistency across all plots.
#
# Usage:
#   make targeted-venn
###############################################################################
cat("\n=========================================\n")
cat(" Targeted Venn Diagrams\n")
cat("=========================================\n\n")

library(VennDiagram)
library(grid)

experiments <- load_all_experiments()

cat("Extracting significant gene sets...\n")
gene_sets <- lapply(experiments, function(df) {
  get_significant_genes(df)
})

# ---- Helper: draw a 2-set Venn with solid colors ----
# Uses VennDiagram::draw.pairwise.venn for clean, standard Venns.
# Two circles get distinct fill colors; the overlap is alpha-blended.
# Colors come from GLOBAL_COLORS for cross-plot consistency.
make_two_set_venn <- function(set_a, set_b, label_a, label_b, title, file_prefix) {
  cat(sprintf("  %s: %d significant proteins\n", label_a, length(set_a)))
  cat(sprintf("  %s: %d significant proteins\n", label_b, length(set_b)))
  overlap_count <- length(intersect(set_a, set_b))
  cat(sprintf("  Overlap: %d proteins\n", overlap_count))

  # Build the Venn using VennDiagram (base R grid graphics)
  # fill: two colors for the circles; alpha = 0.5 blends them in overlap
  # col = "transparent": no circle borders (clean modern look)
  vp <- draw.pairwise.venn(
    area1     = length(set_a),
    area2     = length(set_b),
    cross.area = overlap_count,
    category  = c(label_a, label_b),
    fill      = c(GLOBAL_COLORS[["venn_a_only"]], GLOBAL_COLORS[["venn_b_only"]]),
    alpha     = rep(0.5, 2),
    cat.cex   = 1.6,
    cex       = 2.0,
    fontfamily = "sans",
    cat.fontfamily = "sans",
    col       = "transparent",
    cat.pos   = c(-30, 30),      # Labels at upper edges, away from circle centers
    cat.dist  = c(0.06, 0.06),   # Pushed further out to avoid overlapping circles
    margin    = 0.08,
    ind       = FALSE
  )

  # Add title using grid.text
  pushViewport(viewport())
  grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()

  # Save as PNG and PDF using grid utilities
  commit_hash <- get_git_hash()
  safe_name <- sanitize_filename(file_prefix)
  versioned_name <- paste0(safe_name, "_", commit_hash)

  png_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".png"))
  pdf_path <- file.path(FIGURE_DIR, paste0(versioned_name, ".pdf"))

  # Re-draw to file
  grDevices::png(png_path, width = 7, height = 6, units = "in", res = 300)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  grDevices::dev.off()

  grDevices::pdf(pdf_path, width = 7, height = 6)
  grid.draw(vp)
  pushViewport(viewport())
  grid.text(title, 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  grDevices::dev.off()

  cat(sprintf("  Saved: %s\n  Saved: %s\n", basename(png_path), basename(pdf_path)))

  return(overlap_count)
}

# =====================================================================
# VENN 1: TRIP4 without RA vs TRIP4 with RA
# =====================================================================
cat("\n[1/2] Venn: TRIP4 -RA vs TRIP4 +RA...\n")

exp_base <- "BK467_TRIP4_vs_BK467_WT"
exp_ra   <- "BK467_TRIP4_RA02_vs_BK467_WT"

if (exp_base %in% names(gene_sets) && exp_ra %in% names(gene_sets)) {
  set_a <- gene_sets[[exp_base]]
  set_b <- gene_sets[[exp_ra]]

  overlap <- make_two_set_venn(
    set_a, set_b,
    label_a = "- RA",
    label_b = "+ RA",
    title = "TRIP4 without vs with Retinoic Acid",
    file_prefix = "targeted_venn_RA_effect_BK467"
  )

  ra_gained <- setdiff(set_b, set_a)
  ra_lost <- setdiff(set_a, set_b)

  cat(sprintf("  Shared (both): %d proteins\n", overlap))
  cat(sprintf("  RA-gained: %d proteins\n", length(ra_gained)))
  cat(sprintf("  RA-lost: %d proteins\n", length(ra_lost)))

  save_table(data.frame(gene = intersect(set_a, set_b), category = "RA_shared"),
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

  # Combine C-Flag and N-Flag into one Flag IP set (union)
  flag_sig <- character(0)
  if (exp_cflag %in% names(gene_sets)) {
    flag_sig <- union(flag_sig, gene_sets[[exp_cflag]])
  }
  if (exp_nflag %in% names(gene_sets)) {
    flag_sig <- union(flag_sig, gene_sets[[exp_nflag]])
  }

  if (length(flag_sig) > 0) {
    overlap <- make_two_set_venn(
      turbo_sig, flag_sig,
      label_a = "TurboID",
      label_b = "Flag IP",
      title = "TurboID vs Flag IP Overlap",
      file_prefix = "targeted_venn_turboid_flagip"
    )

    cat(sprintf("  Common to both: %d proteins\n", overlap))
    save_table(data.frame(gene = intersect(turbo_sig, flag_sig),
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
