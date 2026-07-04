###############################################################################
# 09_flagip_volcano.R
# TurboID TRIP4 vs WT volcano with Flag IP common hits labeled.
#
# WHAT THIS PLOTS:
#   A volcano plot of BK467_TRIP4_vs_BK467_WT (the main TurboID experiment).
#   On this volcano, we highlight proteins that are significantly enriched
#   in BOTH Flag IP experiments (C-flag vs Ctrl AND N-flag vs Ctrl).
#   These overlapping proteins are high-confidence TRIP4 interactors
#   validated by two independent methods.
#
# COLORING (color-blind safe Okabe-Ito):
#   - Gray circles:     Not significant
#   - Orange triangles: Significant in TurboID but NOT in both Flag IP
#   - Green squares:    Known interactors (from literature list)
#   - Blue diamonds:    Common Flag IP hits (sig in BOTH C-flag AND N-flag)
#
# LABELS: Only green (known interactors) and blue (Flag IP common) are labeled.
#
# Usage:
#   make flagip-volcano
###############################################################################
cat("\n=========================================\n")
cat(" Flag IP Overlap Volcano\n")
cat("=========================================\n\n")

experiments <- load_all_experiments()
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)

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

# Assign categories (priority: flag_common > known_ia > enriched > nonsig)
df$category <- "nonsig"
df$category[abs(df$log2FC) >= 1 & df$neglog10p > 1] <- "enriched"
df$category[df$gene %in% known_interactors] <- "known_ia"
df$category[df$gene %in% flag_common] <- "flag_common"
df$category <- factor(df$category,
                      levels = c("flag_common", "known_ia", "enriched", "nonsig"))

# Colors (Okabe-Ito color-blind safe)
FLAGIP_COLORS <- c(
  "flag_common" = "#0072B2",   # Blue + Diamond
  "known_ia"    = "#009E73",   # Green + Square
  "enriched"    = "#D55E00",   # Vermillion + Triangle
  "nonsig"      = "#B0B0B0"    # Gray + Circle
)
FLAGIP_SHAPES <- c(
  "flag_common" = 18,          # Diamond
  "known_ia"    = 15,          # Square
  "enriched"    = 17,          # Triangle
  "nonsig"      = 16           # Circle
)

# Only label known interactors and flag IP common hits
label_data <- df[df$category %in% c("flag_common", "known_ia"), ]

p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p,
                                      color = category, shape = category)) +
  ggplot2::geom_point(data = df[df$category == "nonsig", ],
                      alpha = 0.3, size = 0.6) +
  ggplot2::geom_point(data = df[df$category != "nonsig", ],
                      alpha = 0.8, size = 1.5) +
  ggplot2::scale_color_manual(
    values = FLAGIP_COLORS,
    labels = c(
      "flag_common" = "Common Flag IP (C + N)",
      "known_ia"    = "Known interactors",
      "enriched"    = "Significant in TurboID",
      "nonsig"      = "Not significant"
    ),
    name = NULL, drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = FLAGIP_SHAPES,
    labels = c(
      "flag_common" = "Common Flag IP (C + N)",
      "known_ia"    = "Known interactors",
      "enriched"    = "Significant in TurboID",
      "nonsig"      = "Not significant"
    ),
    name = NULL, drop = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 2.8, fontface = "bold",
    max.overlaps = 30, show.legend = FALSE
  ) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TurboID TRIP4 vs WT with Common Flag IP Hits"
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
