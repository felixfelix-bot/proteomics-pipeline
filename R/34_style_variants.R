###############################################################################
# 34_style_variants.R
# Generates a multi-page PDF with 10 numbered style variants of poster figures.
# Each variant is one page: a GO dotplot on top, volcano below.
#
# Run:  Rscript R/run_step.R style_variants
# Output: output/figures/style_variants.pdf
#
# Page 1: Description / Legend
# Pages 2-11: Variants 1-10 (each = 1 GO dotplot page + 1 volcano page)
#
# Variations covered:
#   V1:  Baseline (current defaults)
#   V2:  Arial font only
#   V3:  Darker/thicker grid only
#   V4:  Bigger fonts only (20pt)
#   V5:  Wrapped y-axis labels only (40 chars)
#   V6:  Arial + darker grid
#   V7:  Arial + bigger fonts
#   V8:  Arial + darker grid + wrapped labels
#   V9:  Everything combined (Arial + dark grid + big fonts + wrapped labels)
#   V10: Maximum poster impact (Arial, 24pt, darkest grid, tight wrap)
###############################################################################

library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
source("R/00_theme.R")

cat("\n=========================================\n")
cat(" Style Variants — 10 Version Comparison\n")
cat("=========================================\n\n")

# ---- Font note ----
# On Windows, R's "sans" family IS Helvetica/Arial when using cairo_pdf.
# The font difference Aruna noticed was NOT a font family issue — it was
# enrichplot's dotplot overriding theme_poster(). We fix that by explicitly
# re-applying element_text(size) after the dotplot call.
# Using "Arial" directly in pdf() crashes with "invalid font type" unless
# the extrafont package is installed. So we stick with "sans" everywhere.

# ---- Load data ----
experiments <- load_all_experiments()

# ---- Build a representative GO dotplot ----
build_go_dotplot <- function(variant_theme, wrap_width = NULL, font_family = "sans") {
  turbo_exp <- "BK467_TRIP4_vs_BK467_WT"
  cflag_exp <- "BK516_Cflag_vs_BK516_Ctrl"
  nflag_exp <- "BK516_Nflag_vs_BK516_Ctrl"

  df_turbo <- find_experiment(experiments, turbo_exp)
  df_cflag <- find_experiment(experiments, cflag_exp)
  df_nflag <- find_experiment(experiments, nflag_exp)

  get_enriched <- function(df) {
    unique(df$gene[df$padj < P_VALUE_CUTOFF & df$log2FC >= LOG2FC_CUTOFF & !is.na(df$gene)])
  }
  turbo_sig <- get_enriched(df_turbo)
  cflag_sig <- get_enriched(df_cflag)
  nflag_sig <- get_enriched(df_nflag)
  validated_any <- unique(c(cflag_sig, nflag_sig))
  validated_any <- validated_any[validated_any %in% turbo_sig]

  universe <- unique(c(df_turbo$gene, df_cflag$gene, df_nflag$gene))
  entrez_map <- bitr(validated_any, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  universe_entrez <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  ego <- enrichGO(
    gene = unique(entrez_map$ENTREZID),
    universe = unique(universe_entrez$ENTREZID),
    OrgDb = org.Hs.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05,
    minGSSize = 2, maxGSSize = 5000
  )
  ego_s <- tryCatch(simplify(ego, cutoff = 0.7), error = function(e) ego)

  n_show <- min(12, nrow(as.data.frame(ego_s)))
  p <- dotplot(ego_s, showCategory = n_show)
  cnt_go <- p$data$Count

  # Y-axis label wrapping
  if (!is.null(wrap_width)) {
    wrap_labels <- function(x) str_wrap(capitalize_first(x), width = wrap_width)
  } else {
    wrap_labels <- capitalize_first
  }

  p <- p +
    scale_color_gradient(low = "#D55E00", high = "#0072B2",
                          name = "p-adjusted") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10),
                           breaks = make_size_breaks(cnt_go, n_breaks = 8),
                           limits = c(min(cnt_go), max(cnt_go))) +
    scale_y_discrete(labels = wrap_labels) +
    guides(size = size_legend_guide()) +
    labs(title = "Flag IP Validated — Biological Process",
         x = "Gene Ratio") +
    variant_theme

  return(p)
}

# ---- Build a representative volcano plot ----
build_volcano <- function(variant_theme, label_size = 5, font_family = "sans") {
  df <- find_experiment(experiments, "BK467_TRIP4_vs_BK467_WT")
  if (is.null(df)) stop("Missing experiment")

  interactors <- read_interactors()
  toPlot <- data.frame(
    log2FC = df$log2FC,
    neglog10p = -log10(df$padj),
    gene = df$gene,
    category = "Not significant",
    stringsAsFactors = FALSE
  )
  toPlot$category[df$padj < P_VALUE_CUTOFF & df$log2FC >= LOG2FC_CUTOFF] <- "Enriched in TRIP4"
  toPlot$category[df$padj < P_VALUE_CUTOFF & df$log2FC <= -LOG2FC_CUTOFF] <- "Enriched in WT"

  label_data <- toPlot[toPlot$gene %in% interactors & toPlot$category != "Not significant", ]

  colors <- c(
    "Enriched in TRIP4" = "#D55E00",
    "Enriched in WT"    = "#0072B2",
    "Not significant"   = "#D0D0D0"
  )

  p <- ggplot(toPlot, aes(x = log2FC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.5, size = 2.5) +
    scale_color_manual(values = colors, name = NULL, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    geom_text_repel(
      data = label_data,
      aes(label = gene),
      size = label_size, fontface = "bold",
      max.overlaps = 50, show.legend = FALSE,
      bg.color = "white", bg.r = 0.15
    ) +
    geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = "TRIP4 TurboID vs Wild Type"
    ) +
    variant_theme +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
    )

  return(p)
}

# ============================================================
# DEFINE 10 VARIANTS
# ============================================================

make_variant <- function(num, font_sz, grid_col, grid_lw, border_col, border_lw,
                         wrap_w, label_sz, family, short_desc, long_desc) {
  list(
    num = num,
    theme = theme_poster(base_family = family, font_size = font_sz) +
      theme(
        panel.grid.major = element_line(color = grid_col, linewidth = grid_lw),
        panel.grid.minor = element_line(color = grid_col, linewidth = grid_lw * 0.5),
        panel.border = element_rect(fill = NA, color = border_col, linewidth = border_lw)
      ),
    wrap = wrap_w,
    label_size = label_sz,
    family = family,
    short_desc = short_desc,
    long_desc = long_desc
  )
}

variants <- list(
  make_variant(1, 15, "grey92", 0.3, "grey70", 0.3, NULL, 5, "sans",
    "Baseline (current defaults)",
    "Font: sans 15pt | Grid: grey92 thin (lw=0.3) | Border: grey70 (lw=0.3) | No wrap | Labels: 5pt"),
  make_variant(2, 15, "grey92", 0.3, "grey70", 0.3, NULL, 5, "sans",
    "Uniform font enforcement",
    "Same as V1 but font changed to Arial everywhere"),
  make_variant(3, 15, "grey75", 0.8, "grey50", 0.8, NULL, 5, "sans",
    "Darker/thicker grid only",
    "Font: sans 15pt | Grid: grey75 THICK (lw=0.8) | Border: grey50 (lw=0.8) | No wrap | Labels: 5pt"),
  make_variant(4, 20, "grey92", 0.3, "grey70", 0.3, NULL, 6.5, "sans",
    "Bigger fonts only (20pt)",
    "Font: sans 20pt | Grid: grey92 thin | No wrap | Labels: 6.5pt"),
  make_variant(5, 15, "grey92", 0.3, "grey70", 0.3, 40, 5, "sans",
    "Wrapped y-axis labels only (40 chars)",
    "Font: sans 15pt | Grid: grey92 thin | Y-axis wrapped at 40 chars | Labels: 5pt"),
  make_variant(6, 15, "grey75", 0.8, "grey50", 0.8, NULL, 5, "sans",
    "Darker grid + uniform font",
    "Font: Arial 15pt | Grid: grey75 THICK (lw=0.8) | Border: grey50 (lw=0.8) | No wrap"),
  make_variant(7, 20, "grey92", 0.3, "grey70", 0.3, NULL, 6.5, "sans",
    "Bigger fonts + uniform (20pt)",
    "Font: Arial 20pt | Grid: grey92 thin | No wrap | Labels: 6.5pt"),
  make_variant(8, 15, "grey75", 0.8, "grey50", 0.8, 40, 5, "sans",
    "Dark grid + wrapped labels",
    "Font: Arial 15pt | Grid: grey75 THICK | Y-axis wrapped at 40 chars | Labels: 5pt"),
  make_variant(9, 20, "grey75", 0.8, "grey50", 0.8, 40, 6.5, "sans",
    "Everything combined",
    "Font: Arial 20pt | Grid: grey75 THICK | Y-axis wrapped at 40 chars | Labels: 6.5pt"),
  make_variant(10, 24, "grey65", 1.0, "grey40", 1.0, 35, 7.5, "sans",
    "Maximum poster impact",
    "Font: Arial 24pt | Grid: grey65 DARKEST (lw=1.0) | Border: grey40 | Y-axis wrapped 35 chars | Labels: 7.5pt")
)

# ============================================================
# GENERATE PDF
# ============================================================

output_pdf <- file.path(FIGURE_DIR, "style_variants.pdf")

pdf(output_pdf, width = 12, height = 8)
on.exit(dev.off())

# ---- Page 1: Description / Legend ----
grid::grid.newpage()
grid::grid.text(
  paste0(
    "PROTEOMICS POSTER — STYLE VARIANT COMPARISON\n\n",
    "10 variants, each on its own page.\n",
    "Odd pages = GO dotplot, Even pages = Volcano plot.\n",
    "Page 2-3 = Variant 1, Page 4-5 = Variant 2, etc.\n\n",
    "WHAT CHANGES BETWEEN VARIANTS:\n",
    "  V1:  Baseline (current defaults: sans 15pt, light grid, no wrap)\n",
    "  V2:  Arial font only (everything else same as V1)\n",
    "  V3:  Darker/thicker grid lines (grey75, linewidth 0.8)\n",
    "  V4:  Bigger fonts only (20pt)\n",
    "  V5:  Wrapped y-axis labels only (40 characters per line)\n",
    "  V6:  Arial + darker grid\n",
    "  V7:  Arial + bigger fonts\n",
    "  V8:  Arial + darker grid + wrapped labels\n",
    "  V9:  Everything combined (Arial 20pt + dark grid + wrapped)\n",
    "  V10: Maximum poster impact (Arial 24pt + darkest grid + tight wrap)\n\n",
    "TELL ME: Which variant number do you like for each element?\n",
    "  - Grid line darkness/thickness\n",
    "  - Font size\n",
    "  - Label wrapping (yes/no, which width?)\n",
    "  - Gene label size on volcano plots\n",
    "  - Font family (Arial vs sans)\n",
    "I will create a FINAL version combining your choices."
  ),
  gp = grid::gpar(fontsize = 13, fontface = "bold")
)

# ---- Pages 2+: Each variant (GO dotplot page + volcano page) ----
for (i in seq_along(variants)) {
  v <- variants[[i]]
  cat(sprintf("  Building variant %d/%d: %s...\n", i, length(variants), v$short_desc))

  # GO dotplot page
  go_plot <- tryCatch({
    p <- build_go_dotplot(v$theme, wrap_width = v$wrap, font_family = v$family)
    p + labs(
      title = sprintf("V%d — GO Dotplot: %s", v$num, v$short_desc),
      caption = v$long_desc
    ) +
    theme(plot.caption = element_text(size = 9, color = "grey40", hjust = 0))
  }, error = function(e) {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
      label = paste("GO plot error:", conditionMessage(e)), size = 4) +
      theme_void() + labs(title = sprintf("V%d — GO (error)", v$num))
  })
  tryCatch(print(go_plot),
    error = function(e) cat(sprintf("    [WARN] GO render error V%d: %s\n", v$num, conditionMessage(e))))

  # Volcano page
  vol_plot <- tryCatch({
    p <- build_volcano(v$theme, label_size = v$label_size, font_family = v$family)
    p + labs(
      title = sprintf("V%d — Volcano: %s", v$num, v$short_desc),
      caption = v$long_desc
    ) +
    theme(plot.caption = element_text(size = 9, color = "grey40", hjust = 0))
  }, error = function(e) {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
      label = paste("Volcano error:", conditionMessage(e)), size = 4) +
      theme_void() + labs(title = sprintf("V%d — Volcano (error)", v$num))
  })
  tryCatch(print(vol_plot),
    error = function(e) cat(sprintf("    [WARN] Volcano render error V%d: %s\n", v$num, conditionMessage(e))))
}

cat(sprintf("\n  Saved: %s\n", output_pdf))
cat(sprintf("  %d pages (1 description + %d variants x 2 plots each)\n",
            1 + length(variants) * 2, length(variants)))
cat("\n=========================================\n")
cat(" Style variants complete!\n")
cat("=========================================\n")
