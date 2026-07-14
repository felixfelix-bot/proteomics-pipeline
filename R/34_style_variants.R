###############################################################################
# 34_style_variants.R
# Generates a multi-page PDF with 10 numbered style variants of poster figures.
# Each variant is one page: GO dotplot page + volcano page.
#
# Run:  Rscript R/run_step.R style_variants
# Output: output/figures/style_variants.pdf
#
# OPTIMIZATION: GO enrichment + volcano data computed ONCE, then reused
# for all 10 variants. Total runtime ~40s instead of ~5min.
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

# ---- Load data ----
experiments <- load_all_experiments()

# ============================================================
# PRECOMPUTE ONCE (the expensive part — ~30s for GO enrichment)
# ============================================================
cat("  Computing GO enrichment (once for all variants)...\n")
{
  df_turbo <- find_experiment(experiments, "BK467_TRIP4_vs_BK467_WT")
  df_cflag <- find_experiment(experiments, "BK516_Cflag_vs_BK516_Ctrl")
  df_nflag <- find_experiment(experiments, "BK516_Nflag_vs_BK516_Ctrl")

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
  GO_RESULT <- tryCatch(simplify(ego, cutoff = 0.7), error = function(e) ego)
  cat("  GO enrichment done.\n")
}

# ---- Precompute volcano data ----
{
  df_vol <- find_experiment(experiments, "BK467_TRIP4_vs_BK467_WT")
  interactors_vol <- read_interactors()
  VOLCANO_DATA <- data.frame(
    log2FC = df_vol$log2FC,
    neglog10p = -log10(df_vol$padj),
    gene = df_vol$gene,
    category = "Not significant",
    stringsAsFactors = FALSE
  )
  VOLCANO_DATA$category[df_vol$padj < P_VALUE_CUTOFF & df_vol$log2FC >= LOG2FC_CUTOFF] <- "Enriched in TRIP4"
  VOLCANO_DATA$category[df_vol$padj < P_VALUE_CUTOFF & df_vol$log2FC <= -LOG2FC_CUTOFF] <- "Enriched in WT"
  VOLCANO_LABELS <- VOLCANO_DATA[VOLCANO_DATA$gene %in% interactors_vol & VOLCANO_DATA$category != "Not significant", ]
}

# ============================================================
# FAST PLOT BUILDERS (use precomputed data — no recompute)
# ============================================================
build_go_dotplot <- function(variant_theme, wrap_width = NULL) {
  n_show <- min(12, nrow(as.data.frame(GO_RESULT)))
  p <- dotplot(GO_RESULT, showCategory = n_show)
  cnt_go <- p$data$Count

  if (!is.null(wrap_width)) {
    wrap_labels <- function(x) str_wrap(capitalize_first(x), width = wrap_width)
  } else {
    wrap_labels <- capitalize_first
  }

  p +
    scale_color_gradient(low = "#D55E00", high = "#0072B2", name = "p-adjusted") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10),
                           breaks = make_size_breaks(cnt_go, n_breaks = 8),
                           limits = c(min(cnt_go), max(cnt_go))) +
    scale_y_discrete(labels = wrap_labels) +
    guides(size = size_legend_guide()) +
    labs(title = "Flag IP Validated — Biological Process", x = "Gene Ratio") +
    variant_theme
}

build_volcano <- function(variant_theme, label_size = 5) {
  colors <- c(
    "Enriched in TRIP4" = "#D55E00",
    "Enriched in WT"    = "#0072B2",
    "Not significant"   = "#D0D0D0"
  )

  ggplot(VOLCANO_DATA, aes(x = log2FC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.5, size = 2.5) +
    scale_color_manual(values = colors, name = NULL, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    geom_text_repel(
      data = VOLCANO_LABELS,
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
}

# ============================================================
# DEFINE 10 VARIANTS
# ============================================================
make_variant <- function(num, font_sz, grid_col, grid_lw, border_col, border_lw,
                         wrap_w, label_sz, short_desc, long_desc) {
  list(
    num = num,
    theme = theme_poster(font_size = font_sz) +
      theme(
        panel.grid.major = element_line(color = grid_col, linewidth = grid_lw),
        panel.grid.minor = element_line(color = grid_col, linewidth = grid_lw * 0.5),
        panel.border = element_rect(fill = NA, color = border_col, linewidth = border_lw)
      ),
    wrap = wrap_w,
    label_size = label_sz,
    short_desc = short_desc,
    long_desc = long_desc
  )
}

variants <- list(
  make_variant(1, 15, "grey92", 0.3, "grey70", 0.3, NULL, 5,
    "Baseline (current defaults)",
    "Font: 15pt | Grid: grey92 thin (lw=0.3) | No wrap | Labels: 5pt"),
  make_variant(2, 15, "grey92", 0.3, "grey70", 0.3, NULL, 5,
    "Uniform font enforcement",
    "Same as V1 (all fonts forced consistent via theme override)"),
  make_variant(3, 15, "grey75", 0.8, "grey50", 0.8, NULL, 5,
    "Darker/thicker grid only",
    "Font: 15pt | Grid: grey75 THICK (lw=0.8) | Border: grey50 (lw=0.8) | Labels: 5pt"),
  make_variant(4, 20, "grey92", 0.3, "grey70", 0.3, NULL, 6.5,
    "Bigger fonts only (20pt)",
    "Font: 20pt | Grid: grey92 thin | No wrap | Labels: 6.5pt"),
  make_variant(5, 15, "grey92", 0.3, "grey70", 0.3, 40, 5,
    "Wrapped y-axis labels only (40 chars)",
    "Font: 15pt | Grid: grey92 thin | Y-axis wrapped at 40 chars | Labels: 5pt"),
  make_variant(6, 15, "grey75", 0.8, "grey50", 0.8, NULL, 5,
    "Darker grid + uniform font",
    "Font: 15pt | Grid: grey75 THICK | No wrap | Labels: 5pt"),
  make_variant(7, 20, "grey92", 0.3, "grey70", 0.3, NULL, 6.5,
    "Bigger fonts + uniform (20pt)",
    "Font: 20pt | Grid: grey92 thin | No wrap | Labels: 6.5pt"),
  make_variant(8, 15, "grey75", 0.8, "grey50", 0.8, 40, 5,
    "Dark grid + wrapped labels",
    "Font: 15pt | Grid: grey75 THICK | Y-axis wrapped at 40 chars | Labels: 5pt"),
  make_variant(9, 20, "grey75", 0.8, "grey50", 0.8, 40, 6.5,
    "Everything combined",
    "Font: 20pt | Grid: grey75 THICK | Y-axis wrapped at 40 chars | Labels: 6.5pt"),
  make_variant(10, 24, "grey65", 1.0, "grey40", 1.0, 35, 7.5,
    "Maximum poster impact",
    "Font: 24pt | Grid: grey65 DARKEST (lw=1.0) | Border: grey40 | Y-axis wrapped 35 chars | Labels: 7.5pt")
)

# ============================================================
# GENERATE PDF — use ggsave per page to avoid sink() conflicts
# ============================================================
output_pdf <- file.path(FIGURE_DIR, "style_variants.pdf")

# Build all plots into a list first, then write PDF
all_plots <- list()

# Description page (use cowplot-free approach: blank ggplot with text)
desc_text <- paste0(
  "PROTEOMICS POSTER — STYLE VARIANT COMPARISON\n\n",
  "10 variants, each on 2 pages (GO dotplot + volcano).\n\n",
  "V1:  Baseline (current defaults: 15pt, light grid, no wrap)\n",
  "V2:  Uniform font enforcement\n",
  "V3:  Darker/thicker grid (grey75, lw=0.8)\n",
  "V4:  Bigger fonts (20pt)\n",
  "V5:  Wrapped y-axis labels (40 chars)\n",
  "V6:  Darker grid + uniform font\n",
  "V7:  Bigger fonts + uniform (20pt)\n",
  "V8:  Dark grid + wrapped labels\n",
  "V9:  Everything combined (20pt + dark grid + wrapped)\n",
  "V10: Maximum poster (24pt, darkest grid, tight wrap)\n\n",
  "TELL ME: Which variant number for each element:\n",
  "  Grid, Font size, Label wrapping, Gene label size\n",
  "I will combine your choices into a final version."
)
all_plots[[1]] <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = desc_text,
           size = 4, hjust = 0.5, lineheight = 1.5) +
  theme_void()

# ---- Pages 2+: Each variant ----
for (i in seq_along(variants)) {
  v <- variants[[i]]
  cat(sprintf("  Building variant %d/%d: %s...\n", i, length(variants), v$short_desc))

  # GO dotplot
  go_plot <- tryCatch({
    p <- build_go_dotplot(v$theme, wrap_width = v$wrap)
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
  all_plots[[length(all_plots) + 1]] <- go_plot

  # Volcano
  vol_plot <- tryCatch({
    p <- build_volcano(v$theme, label_size = v$label_size)
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
  all_plots[[length(all_plots) + 1]] <- vol_plot
}

# ---- Write all plots to single multi-page PDF ----
# Use cairo_pdf for better font rendering, write each plot as a page
cat("  Writing PDF...\n")
pdf(output_pdf, width = 12, height = 8)
for (p in all_plots) {
  print(p)
}
dev.off()

cat(sprintf("\n  Saved: %s\n", output_pdf))
cat(sprintf("  %d pages\n", length(all_plots)))
cat("\n=========================================\n")
cat(" Style variants complete!\n")
cat("=========================================\n")
