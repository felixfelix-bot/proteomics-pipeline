###############################################################################
# 09_flagip_volcano.R
# TurboID TRIP4 vs WT volcano with Flag IP common hits highlighted.
#
# CHANGES:
#   - All circles, same size (no triangles/squares/diamonds)
#   - Do NOT label known interactors
#   - ONLY label proteins that are significant in TurboID AND C-Flag AND N-Flag
#   - These triple-validated hits are the highest-confidence interactors
#
# COLORING:
#   - Light gray:     Not significant
#   - Orange circles: Significant in TurboID (or Flag IP but not all three)
#   - Blue circles:   Significant in all three: TurboID + C-Flag + N-Flag
#
# Usage:
#   make flagip-volcano
###############################################################################
cat("\n=========================================\n")
cat(" Flag IP Overlap Volcano\n")
cat("=========================================\n\n")

experiments <- load_all_experiments()

# ---- Find proteins significant in BOTH Flag IP experiments ----
cat("Finding common Flag IP hits...\n")

flag_c <- "BK516_Cflag_vs_BK516_Ctrl"
flag_n <- "BK516_Nflag_vs_BK516_Ctrl"

if (flag_c %in% names(experiments) && flag_n %in% names(experiments)) {
  cflag_sig <- get_significant_genes(experiments[[flag_c]])
  nflag_sig <- get_significant_genes(experiments[[flag_n]])
  flag_common <- intersect(cflag_sig, nflag_sig)

  cat(sprintf("  C-flag significant: %d proteins\n", length(cflag_sig)))
  cat(sprintf("  N-flag significant: %d proteins\n", length(nflag_sig)))
  cat(sprintf("  Common (in BOTH): %d proteins\n", length(flag_common)))

  save_table(data.frame(gene = flag_common, category = "FlagIP_common"),
             "flagip_common_significant")
} else {
  cat("  WARNING: Flag IP experiments not found. Using empty list.\n")
  flag_common <- character(0)
}

# ---- Build the volcano ----
cat("\nBuilding volcano plot...\n")
turbo_main <- "BK467_TRIP4_vs_BK467_WT"

if (!turbo_main %in% names(experiments)) {
  cat("ERROR: Main experiment not found:", turbo_main, "\n")
  quit(status = 1)
}

df <- experiments[[turbo_main]]
df$neglog10p <- -log10(df$padj)

# Categories: only flag_common (blue) and everything else (gray/orange)
# Do NOT separately highlight known interactors
df$category <- "nonsig"
df$category[abs(df$log2FC) >= 1 & df$neglog10p > 1] <- "enriched"
df$category[df$gene %in% flag_common] <- "flag_common"
df$category <- factor(df$category,
                      levels = c("flag_common", "enriched", "nonsig"))

# ONLY label the triple-validated hits (flag_common proteins)
label_data <- df[df$category == "flag_common", ]

FLAGIP_COLORS <- c(
  "flag_common" = "#0072B2",   # Blue — common to TurboID + C-Flag + N-Flag
  "enriched"    = "#D55E00",   # Orange — significant in TurboID only
  "nonsig"      = "#D0D0D0"    # Gray
)

# All circles now — researcher said no triangles
p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.8) +
  ggplot2::scale_color_manual(
    values = FLAGIP_COLORS,
    labels = c(
      "flag_common" = "Common: TurboID + C-Flag + N-Flag",
      "enriched"    = "Significant in TurboID",
      "nonsig"      = "Not significant"
    ),
    name = NULL, drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 2.5, fontface = "bold",
    max.overlaps = 30, show.legend = FALSE
  ) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID vs WT with Common Flag IP Hits"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = "right",
    panel.grid.minor = ggplot2::element_blank()
  )

save_figure(p, "flagip_overlap_volcano_BK467_TRIP4_vs_WT",
            width = 8, height = 6)

cat("\n=========================================\n")
cat(" Flag IP overlap volcano complete!\n")
cat("=========================================\n")
