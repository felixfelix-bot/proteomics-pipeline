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
  interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
  interactors_vol <- load_known_interactors(interactors_file)
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
build_go_dotplot <- function(variant_theme, wrap_width = 30) {
  n_show <- min(12, nrow(as.data.frame(GO_RESULT)))
  p <- dotplot(GO_RESULT, showCategory = n_show)
  cnt_go <- p$data$Count

  wrap_labels <- function(x) str_wrap(capitalize_first(x), width = wrap_width)

  p +
    scale_color_gradient(low = "#D55E00", high = "#0072B2", name = "p-adjusted") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10),
                           breaks = make_size_breaks(cnt_go, n_breaks = 8),
                           limits = c(min(cnt_go), max(cnt_go))) +
    scale_y_discrete(labels = wrap_labels) +
    guides(size = size_legend_guide()) +
    labs(x = "Gene Ratio") +
    variant_theme
}

build_volcano <- function(variant_theme, label_size = 5) {
  colors <- c(
    "Enriched in TRIP4" = "#D55E00",
    "Enriched in WT"    = "#0072B2",
    "Not significant"   = "#D0D0D0"
  )

  # Clip x-axis: outlier at ~-10 wastes space, clip at -7
  x_min <- -7

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
    coord_cartesian(xlim = c(x_min, NA)) +
    labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value))
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
  make_variant(1, 15, "grey75", 0.6, "grey50", 0.5, 30, 5,
    "Baseline (dark grid + wrapped labels)",
    "Font: 15pt | Grid: grey75 (lw=0.6) | Wrap: 30ch | Labels: 5pt"),
  make_variant(2, 18, "grey75", 0.6, "grey50", 0.5, 30, 6,
    "Bigger fonts (18pt)",
    "Font: 18pt | Grid: grey75 (lw=0.6) | Wrap: 30ch | Labels: 6pt"),
  make_variant(3, 15, "grey60", 0.8, "grey40", 0.8, 30, 5,
    "Darker grid (grey60)",
    "Font: 15pt | Grid: grey60 THICK (lw=0.8) | Wrap: 30ch | Labels: 5pt"),
  make_variant(4, 20, "grey60", 0.8, "grey40", 0.8, 30, 6.5,
    "Bigger fonts + dark grid (20pt)",
    "Font: 20pt | Grid: grey60 (lw=0.8) | Wrap: 30ch | Labels: 6.5pt"),
  make_variant(5, 15, "grey75", 0.6, "grey50", 0.5, 22, 5,
    "Tighter wrap (22 chars)",
    "Font: 15pt | Grid: grey75 | Wrap: 22ch (tighter) | Labels: 5pt"),
  make_variant(6, 18, "grey60", 0.8, "grey40", 0.8, 25, 6,
    "Medium combo (18pt + dark grid + 25ch wrap)",
    "Font: 18pt | Grid: grey60 | Wrap: 25ch | Labels: 6pt"),
  make_variant(7, 22, "grey60", 0.8, "grey40", 0.8, 30, 7,
    "Poster-ready (22pt)",
    "Font: 22pt | Grid: grey60 | Wrap: 30ch | Labels: 7pt"),
  make_variant(8, 18, "grey50", 1.0, "grey30", 1.0, 25, 6,
    "Darkest grid (grey50)",
    "Font: 18pt | Grid: grey50 DARKEST (lw=1.0) | Wrap: 25ch | Labels: 6pt")
)

# ============================================================
# GENERATE PDF — use ggsave per page to avoid sink() conflicts
# ============================================================
commit_hash <- get_git_hash()
output_pdf <- file.path(FIGURE_DIR, paste0("style_variants_", commit_hash, ".pdf"))

# Build all plots into a list first, then write PDF
all_plots <- list()

# Description page (use cowplot-free approach: blank ggplot with text)
desc_text <- paste0(
  "PROTEOMICS POSTER — STYLE VARIANT COMPARISON\n\n",
  "8 variants, each on 2 pages (GO dotplot + volcano).\n",
  "ALL variants: dark grid + wrapped labels (30ch default).\n",
  "Volcano x-axis clipped at -7 (outlier removed).\n\n",
  "V1:  Baseline — 15pt, grey75 grid, 30ch wrap\n",
  "V2:  Bigger fonts — 18pt\n",
  "V3:  Darker grid — grey60, lw=0.8\n",
  "V4:  Bigger + dark — 20pt, grey60\n",
  "V5:  Tighter wrap — 22 chars\n",
  "V6:  Medium combo — 18pt, grey60, 25ch wrap\n",
  "V7:  Poster-ready — 22pt\n",
  "V8:  Darkest grid — grey50, lw=1.0\n\n",
  "TELL ME: Which variant number for each element:\n",
  "  Grid darkness, Font size, Label wrapping\n",
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
cat("  Writing PDF...\n")
pdf(output_pdf, width = 12, height = 8)
tryCatch({
  for (p in all_plots) {
    tryCatch(print(p),
      error = function(e) {
        cat(sprintf("    [WARN] render error: %s\n", conditionMessage(e)))
        # Print a blank placeholder so the page isn't missing
        print(ggplot() + theme_void())
      })
  }
}, finally = {
  dev.off()
})

cat(sprintf("\n  Saved: %s\n", output_pdf))
cat(sprintf("  %d pages\n", length(all_plots)))
cat("\n=========================================\n")
cat(" Style variants complete!\n")
cat("=========================================\n")
