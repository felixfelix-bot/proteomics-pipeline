###############################################################################
# 34_style_variants.R
# Generates a multi-page PDF comparing different styling options for
# poster figures. Each page shows one variant with a description header.
#
# Run:  Rscript R/run_step.R style_variants
# Output: output/figures/style_variants.pdf
#
# The script picks representative figures (one GO dotplot + one volcano)
# and renders them in 6 different style variants so Aruna can pick
# which elements she likes best.
###############################################################################

library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
source("R/00_theme.R")

cat("\n=========================================\n")
cat(" Style Variants — Multi-Version Comparison\n")
cat("=========================================\n\n")

# ---- Load data (same as other scripts) ----
experiments <- load_all_experiments()

# ---- Build a representative GO dotplot with long y-axis labels ----
build_go_dotplot <- function(variant_theme, wrap_width = NULL, font_scale = 1) {
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

# ---- Build a representative volcano plot with gene labels ----
build_volcano <- function(variant_theme, label_size = 5) {
  df <- find_experiment(experiments, "BK467_TRIP4_vs_BK467_WT")
  if (is.null(df)) stop("Missing experiment")

  interactors <- read_interactors()
  sig <- df[df$padj < P_VALUE_CUTOFF & abs(df$log2FC) >= LOG2FC_CUTOFF & !is.na(df$gene), ]
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
  labels <- c(
    "Enriched in TRIP4" = "Enriched in TRIP4",
    "Enriched in WT"    = "Enriched in WT",
    "Not significant"   = "Not significant"
  )

  p <- ggplot(toPlot, aes(x = log2FC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.5, size = 2.5) +
    scale_color_manual(values = colors, labels = labels, name = NULL, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    geom_text_repel(
      data = label_data,
      aes(label = gene),
      size = label_size, fontface = "bold",
      max.overlaps = 50, show.legend = FALSE,
      bg.color = "white", bg.r = 0.15
    ) +
    geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = "TRIP4 TurboID vs Wild Type"
    ) +
    variant_theme +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      coord_cartesian(clip = "off")  # not valid here, but keeping for compatibility
    )

  return(p)
}

# ---- Helper: add title annotation to each variant ----
annotate_variant <- function(plot, variant_num, description) {
  title_text <- sprintf("Variant %d: %s", variant_num, description)
  plot + labs(title = title_text, caption = description)
}

# ---- Define 6 theme variants ----

# Variant 1: BASELINE (current theme_poster defaults)
v1_theme <- theme_poster(font_size = 15)
v1_desc <- "Baseline (current: font=15pt, grid=grey92 thin, no wrap)"
v1_wrap <- NULL
v1_label_size <- 5

# Variant 2: DARKER THICKER GRID only
v2_theme <- theme_poster(font_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey75", linewidth = 0.8),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.4),
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.8)
  )
v2_desc <- "Darker/thicker grid (grid=grey75 lw=0.8, border=grey50 lw=0.8)"
v2_wrap <- NULL
v2_label_size <- 5

# Variant 3: BIGGER FONTS only
v3_theme <- theme_poster(font_size = 20)
v3_desc <- "Bigger fonts (font=20pt, grid=grey92 thin, no wrap)"
v3_wrap <- NULL
v3_label_size <- 6.5

# Variant 4: WRAPPED Y-AXIS only
v4_theme <- theme_poster(font_size = 15)
v4_desc <- "Wrapped y-axis labels (40 chars, font=15pt)"
v4_wrap <- 40
v4_label_size <- 5

# Variant 5: ALL COMBINED (best guess)
v5_theme <- theme_poster(font_size = 20) +
  theme(
    panel.grid.major = element_line(color = "grey75", linewidth = 0.8),
    panel.grid.minor = element_line(color = "grey85", linewidth = 0.4),
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.8)
  )
v5_desc <- "All combined: bigger fonts (20pt) + darker grid + wrapped labels (40ch)"
v5_wrap <- 40
v5_label_size <- 6.5

# Variant 6: POSTER-READY (even bigger, for reading from distance)
v6_theme <- theme_poster(font_size = 24) +
  theme(
    panel.grid.major = element_line(color = "grey65", linewidth = 1.0),
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.5),
    panel.border = element_rect(fill = NA, color = "grey40", linewidth = 1.0)
  )
v6_desc <- "Poster-ready: font=24pt, grid=grey65 lw=1.0, labels wrapped (35ch)"
v6_wrap <- 35
v6_label_size <- 7.5

# ---- Generate all variants ----
variants <- list(
  list(theme = v1_theme, desc = v1_desc, wrap = v1_wrap, label = v1_label_size),
  list(theme = v2_theme, desc = v2_desc, wrap = v2_wrap, label = v2_label_size),
  list(theme = v3_theme, desc = v3_desc, wrap = v3_wrap, label = v3_label_size),
  list(theme = v4_theme, desc = v4_desc, wrap = v4_wrap, label = v4_label_size),
  list(theme = v5_theme, desc = v5_desc, wrap = v5_wrap, label = v5_label_size),
  list(theme = v6_theme, desc = v6_desc, wrap = v6_wrap, label = v6_label_size)
)

output_pdf <- file.path(FIGURE_DIR, "style_variants.pdf")

# Open PDF device
pdf(output_pdf, width = 14, height = 10)
on.exit(dev.off())

# ---- Page 0: Description / Legend ----
grid::grid.newpage()
grid::grid.text(
  paste0(
    "PROTEOMICS POSTER — STYLE VARIANT COMPARISON\n\n",
    "6 variants shown. Each page = one variant (GO dotplot + volcano side by side).\n\n",
    "Variant 1: BASELINE — current defaults (font=15, grid=grey92 thin, no wrap)\n",
    "Variant 2: DARKER GRID — grid=grey75 lw=0.8, border=grey50 lw=0.8\n",
    "Variant 3: BIGGER FONTS — font=20pt\n",
    "Variant 4: WRAPPED LABELS — y-axis wrapped at 40 characters\n",
    "Variant 5: ALL COMBINED — font=20 + darker grid + wrapped (40ch)\n",
    "Variant 6: POSTER-READY — font=24, darkest grid, wrapped (35ch)\n\n",
    "Tell me which variant number you like for EACH element:\n",
    "  Grid lines, Font size, Label wrapping, Gene label size\n",
    "I will create a final version combining your choices."
  ),
  gp = grid::gpar(fontsize = 14, fontface = "bold")
)

# ---- Generate each variant page ----
for (i in seq_along(variants)) {
  v <- variants[[i]]
  cat(sprintf("  Building variant %d/%d...\n", i, length(variants)))

  go_plot <- tryCatch({
    p <- build_go_dotplot(v$theme, wrap_width = v$wrap)
    p + labs(title = sprintf("Variant %d — GO Dotplot\n%s", i, v$desc))
  }, error = function(e) {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
      label = paste("GO plot error:", conditionMessage(e)), size = 4) +
      theme_void() + labs(title = sprintf("Variant %d — GO (error)", i))
  })
  print(go_plot)

  vol_plot <- tryCatch({
    p <- build_volcano(v$theme, label_size = v$label)
    p + labs(title = sprintf("Variant %d — Volcano\n%s", i, v$desc))
  }, error = function(e) {
    ggplot() + annotate("text", x = 0.5, y = 0.5,
      label = paste("Volcano error:", conditionMessage(e)), size = 4) +
      theme_void() + labs(title = sprintf("Variant %d — Volcano (error)", i))
  })
  print(vol_plot)
}

cat(sprintf("\n  Saved: %s\n", output_pdf))
cat(sprintf("  %d pages (1 description + %d variants)\n", length(variants) + 1, length(variants)))
cat("\n=========================================\n")
cat(" Style variants complete!\n")
cat("=========================================\n")
