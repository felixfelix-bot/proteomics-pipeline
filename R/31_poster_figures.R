###############################################################################
# 31_poster_figures.R
# Generate ALL figures for the LaTeX poster with STANDARDIZED FONT SIZES.
#
# This script:
#   1. Sources R/00_theme.R (defines theme_poster with uniform fonts)
#   2. Loads experiment data
#   3. For GO/KEGG: reads result tables from output/tables/ (already computed
#      by targeted-go, etc.) and re-creates dotplots/barplots with theme_poster()
#   4. For volcano/Venn: regenerates from raw data with theme_poster()
#   5. Saves ALL figures as PDF to poster/figures/ (vector format for LaTeX)
#
# Prerequisites: Run `make aruna-fast` FIRST to generate GO/KEGG tables.
#
# Usage:
#   make poster      (runs aruna-fast, then this script, then LaTeX)
#   make poster-figures  (runs this script only)
###############################################################################
cat("\n=========================================\n")
cat(" Poster Figures — Arial 20, Black, GO-dot/bar, synced legends\n")
cat("=========================================\n\n")

# ---- Source theme (must come BEFORE any plotting) ----
# This defines theme_poster(), POSTER_FONTS, save_poster_figure()
source("R/00_theme.R")

library(ggplot2)
library(ggrepel)

# ---- Output directory ----
POSTER_FIG_DIR <- "poster/figures"
dir.create(POSTER_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: save directly to poster/figures/ (no git hash, clean names)
save_poster <- function(plot, name, width = 8, height = 6) {
  pdf_path <- file.path(POSTER_FIG_DIR, paste0(name, ".pdf"))
  png_path <- file.path(POSTER_FIG_DIR, paste0(name, ".png"))
  svg_path <- file.path(POSTER_FIG_DIR, paste0(name, ".svg"))
  ggsave(pdf_path, plot = plot, width = width, height = height,
         device = grDevices::cairo_pdf)
  ggsave(png_path, plot = plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(svg_path, plot = plot, width = width, height = height, bg = "white")
  cat(sprintf("  [POSTER] %s.pdf + .png + .svg (%.0f x %.0f in)\n", name, width, height))
}

# ---- capitalize_first() is defined in R/00_theme.R (sourced above) ----

# ---- Load experiments ----
experiments <- load_all_experiments()
turbo_name <- "BK467_TRIP4_vs_BK467_WT"
df_turbo <- find_experiment(experiments, turbo_name)
if (is.null(df_turbo)) {
  cat("ERROR: Cannot find", turbo_name, "\n")
  quit(status = 1)
}
df_turbo$neglog10p <- -log10(df_turbo$padj)
cat(sprintf("  Loaded %d proteins from %s\n", nrow(df_turbo), turbo_name))

# Known interactors + ASCC core
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- if (file.exists(interactors_file)) {
  readLines(interactors_file)
} else {
  character(0)
}
ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
FC_THRESH <- 0.999

cat("\n--- Generating poster figures ---\n\n")

#=============================================================================
# SECTION A — VOLCANO PLOTS (5 total)
#=============================================================================

#---------------------------------------------------------------------------
# A1. VOLCANO PLOT — TurboID TRIP4 vs WT (existing, updated fonts)
#---------------------------------------------------------------------------
cat("[A1] TurboID TRIP4 vs WT volcano...\n")

df_volcano <- df_turbo
df_volcano$category <- "Not enriched"
df_volcano$category[!is.na(df_volcano$log2FC) &
  df_volcano$log2FC >= FC_THRESH &
  df_volcano$padj < P_VALUE_CUTOFF] <- "Enriched"
df_volcano$category[df_volcano$gene %in% ASCC_CORE &
  df_volcano$log2FC >= FC_THRESH] <- "ASCC complex"
df_volcano$category[df_volcano$gene %in% known_interactors &
  df_volcano$log2FC >= FC_THRESH] <- "Known interactor"

df_volcano$category <- factor(df_volcano$category,
  levels = c("ASCC complex", "Known interactor", "Enriched", "Not enriched"))

volcano_colors <- c(
  "ASCC complex"    = "#0072B2",
  "Known interactor" = "#009E73",
  "Enriched"        = "#D55E00",
  "Not enriched"    = "grey80"
)

label_genes <- unique(c(ASCC_CORE,
  intersect(known_interactors,
    df_volcano$gene[df_volcano$category != "Not enriched"])))

p_volcano <- ggplot(df_volcano,
    aes(x = log2FC, y = neglog10p, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colors, name = NULL) +
  geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = LOG2FC_CUTOFF, linetype = "dashed", color = "grey50") +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TurboID TRIP4 vs Wild Type"
  ) +
  geom_text_repel(
    data = df_volcano[df_volcano$gene %in% label_genes &
      !is.na(df_volcano$log2FC), ],
    aes(label = gene), size = 4, fontface = "bold",
    max.overlaps = 20, show.legend = FALSE, color = "black"
  ) +
  theme_poster() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "grey80"))

save_poster(p_volcano, "volcano_turbo", width = 16, height = 12)

#---------------------------------------------------------------------------
# A2. FLAG IP OVERLAY VOLCANO (existing, updated fonts)
#---------------------------------------------------------------------------
cat("[A2] Flag IP overlap volcano...\n")

cflag_exp <- find_experiment(experiments, "BK516_Cflag_vs_BK516_Ctrl")
nflag_exp <- find_experiment(experiments, "BK516_Nflag_vs_BK516_Ctrl")

cflag_sig <- if (!is.null(cflag_exp)) get_significant_genes(cflag_exp) else character(0)
nflag_sig <- if (!is.null(nflag_exp)) get_significant_genes(nflag_exp) else character(0)

df_flag <- df_turbo
df_flag$flag_cat <- "TRIP4 only"
turbo_sig <- df_flag$gene[df_flag$padj < P_VALUE_CUTOFF & df_flag$log2FC >= FC_THRESH]
df_flag$flag_cat[df_flag$gene %in% intersect(cflag_sig, nflag_sig) &
  df_flag$gene %in% turbo_sig] <- "Both C+N Flag"
df_flag$flag_cat[df_flag$gene %in% setdiff(cflag_sig, nflag_sig) &
  df_flag$gene %in% turbo_sig] <- "C-Flag only"
df_flag$flag_cat[df_flag$gene %in% setdiff(nflag_sig, cflag_sig) &
  df_flag$gene %in% turbo_sig] <- "N-Flag only"
df_flag$flag_cat[!(df_flag$gene %in% turbo_sig)] <- "Not enriched"

df_flag$flag_cat <- factor(df_flag$flag_cat,
  levels = c("Both C+N Flag", "C-Flag only", "N-Flag only",
             "TRIP4 only", "Not enriched"))

flag_colors <- c(
  "Both C+N Flag" = "#009E73",
  "C-Flag only"  = "#0072B2",
  "N-Flag only"  = "#CC79A7",
  "TRIP4 only"   = "#D55E00",
  "Not enriched" = "grey80"
)

label_flag <- intersect(intersect(cflag_sig, nflag_sig), turbo_sig)

p_flag <- ggplot(df_flag,
    aes(x = log2FC, y = neglog10p, color = flag_cat)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = flag_colors, name = NULL) +
  geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = LOG2FC_CUTOFF, linetype = "dashed", color = "grey50") +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "Flag IP Validation Overlay"
  ) +
  geom_text_repel(
    data = df_flag[df_flag$gene %in% label_flag & !is.na(df_flag$log2FC), ],
    aes(label = gene), size = 4, fontface = "bold",
    max.overlaps = 15, show.legend = FALSE, color = "black"
  ) +
  theme_poster() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "grey80"))

save_poster(p_flag, "volcano_flagip", width = 16, height = 12)

#---------------------------------------------------------------------------
# A3-A5. CHX/DMSO VOLCANO PLOTS (new)
# Generic helper: produce a colored volcano for a CHX/DMSO comparison.
#   - enriched (CHX-side, log2FC > 0) in chx_enriched color
#   - depleted (DMSO/WT-side, log2FC < 0) in dmso_enriched color
#   - label top 15 per side by combined significance score
#---------------------------------------------------------------------------

# Internal helper: build one CHX/DMSO volcano from a single experiment data frame.
make_chx_volcano <- function(df, experiment_name, plot_title, out_name) {
  if (is.null(df)) {
    cat(sprintf("  SKIP (experiment not loaded): %s\n", experiment_name))
    return(invisible(NULL))
  }
  cat(sprintf("[CHX] Building volcano for %s ...\n", experiment_name))

  df <- df[!is.na(df$log2FC) & !is.na(df$padj), ]
  df$neglog10p <- -log10(df$padj)

  df$category <- "Not significant"
  df$category[df$padj < P_VALUE_CUTOFF & df$log2FC >= LOG2FC_CUTOFF]  <- "Enriched (CHX side)"
  df$category[df$padj < P_VALUE_CUTOFF & df$log2FC <= -LOG2FC_CUTOFF] <- "Depleted (DMSO side)"

  df$category <- factor(df$category,
    levels = c("Enriched (CHX side)", "Depleted (DMSO side)", "Not significant"))

  chx_colors <- c(
    "Enriched (CHX side)"   = GLOBAL_COLORS["chx_enriched"],
    "Depleted (DMSO side)"  = GLOBAL_COLORS["dmso_enriched"],
    "Not significant"       = "grey80"
  )

  # Combined significance score for top-N labeling
  df$score <- abs(df$log2FC) * df$neglog10p
  top_up   <- head(df$gene[df$category == "Enriched (CHX side)"][
    order(df$score[df$category == "Enriched (CHX side)"], decreasing = TRUE)], 15)
  top_dn   <- head(df$gene[df$category == "Depleted (DMSO side)"][
    order(df$score[df$category == "Depleted (DMSO side)"], decreasing = TRUE)], 15)
  top_genes <- unique(c(top_up, top_dn))

  p <- ggplot(df, aes(x = log2FC, y = neglog10p, color = category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = chx_colors, name = NULL) +
    geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
               linetype = "dashed", color = "grey50") +
    labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = plot_title
    ) +
    geom_text_repel(
      data = df[df$gene %in% top_genes, ],
      aes(label = gene), size = 4, fontface = "bold",
      max.overlaps = 20, show.legend = FALSE, color = "black"
    ) +
    theme_poster() +
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", color = "grey80"))

  save_poster(p, out_name, width = 16, height = 12)
}

# A3. CHX vs DMSO
make_chx_volcano(
  df             = find_experiment(experiments, "TRIP4_CHX_vs_TRIP4_DMSO"),
  experiment_name = "TRIP4_CHX_vs_TRIP4_DMSO",
  plot_title      = "CHX vs DMSO (TRIP4)",
  out_name        = "volcano_chx_vs_dmso"
)

# A4. CHX vs WT
make_chx_volcano(
  df             = find_experiment(experiments, "TRIP4_CHX_vs_WT"),
  experiment_name = "TRIP4_CHX_vs_WT",
  plot_title      = "CHX vs Wild Type",
  out_name        = "volcano_chx_vs_wt"
)

# A5. DMSO vs WT
make_chx_volcano(
  df             = find_experiment(experiments, "TRIP4_DMSO_vs_WT"),
  experiment_name = "TRIP4_DMSO_vs_WT",
  plot_title      = "DMSO vs Wild Type",
  out_name        = "volcano_dmso_vs_wt"
)

#=============================================================================
# VENN DIAGRAM — TurboID vs C-Flag vs N-Flag (existing, fonts updated)
#=============================================================================
cat("[VENN] TurboID vs C-Flag vs N-Flag...\n")

venn_data <- list(
  "TurboID TRIP4" = turbo_sig,
  "C-Flag IP" = cflag_sig,
  "N-Flag IP" = nflag_sig
)

p_venn <- tryCatch({
  library(ggVennDiagram)
  p <- ggVennDiagram(venn_data,
    label_alpha = 0.7,
    edge_size = 0.5,
    label = "count") +
    scale_fill_gradient(low = "white", high = "#0072B2") +
    labs(title = "TurboID vs Flag IP Overlap") +
    theme_poster() +
    theme(legend.position = "right")
  p
}, error = function(e) {
  cat("  ggVennDiagram not available, trying VennDiagram...\n")
  library(VennDiagram)
  venn_file <- file.path(POSTER_FIG_DIR, "venn_overlap.pdf")
  venn.plot <- venn.diagram(
    x = venn_data,
    filename = venn_file,
    col = c("#D55E00", "#0072B2", "#CC79A7"),
    fill = c("#D55E00", "#0072B2", "#CC79A7"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.5,
    fontfamily = "Arial",
    cat.fontfamily = "Arial",
    main = "TurboID vs Flag IP Overlap",
    main.cex = 1.5,
    main.fontfamily = "Arial"
  )
  NULL
})

if (!is.null(p_venn)) {
  save_poster(p_venn, "venn_overlap", width = 16, height = 12)
}

#=============================================================================
# SECTION B — GO DOT PLOTS AND BAR PLOTS
# Reads GO/KEGG result CSVs from output/tables/, creates BOTH a dotplot
# and a barplot for each table found.
#=============================================================================

# Find the most recent table file matching a pattern
find_latest_table <- function(pattern) {
  files <- list.files(TABLE_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

# ---- Helper: standardized dotplot from a GO/KEGG result table ----
#  - Capitalizes first letter of each GO term
#  - BIG dots: size range c(4, 10)
#  - Synced legends: size override matches actual dot range
make_poster_dotplot <- function(table_path, title_text, max_terms = 15) {
  if (!file.exists(table_path)) {
    cat(sprintf("  SKIP (table not found): %s\n", basename(table_path)))
    return(NULL)
  }

  df <- read.csv(table_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    cat(sprintf("  SKIP (empty table): %s\n", basename(table_path)))
    return(NULL)
  }

  # Parse GeneRatio (stored as "3/10" string -> numeric)
  if ("GeneRatio" %in% names(df)) {
    df$ratio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(strsplit(as.character(x), "/")[[1]])
      if (length(parts) == 2 && parts[2] > 0) parts[1] / parts[2] else NA
    })
  } else {
    cat(sprintf("  SKIP (no GeneRatio column): %s\n", basename(table_path)))
    return(NULL)
  }

  # Capitalize GO terms
  if ("Description" %in% names(df)) {
    df$Description <- capitalize_first(df$Description)
  }

  # Sort by p.adjust (most significant first), keep top N
  df <- df[order(df$p.adjust), ]
  df <- head(df, max_terms)
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = ratio, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#D55E00", high = "#0072B2",
                         name = "p-adjusted") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10),
                          breaks = make_size_breaks(df$Count)) +
    labs(
      x = "Gene Ratio",
      y = NULL,
      title = title_text
    ) +
    guides(
      size  = guide_legend(title.position = "top",
                           override.aes = list(color = "#D55E00")),
      color = guide_colorbar(title.position = "top",
                             barwidth = 1.5, barheight = 8)
    ) +
    theme_poster(font_size = 24) +
    theme(legend.position = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.box = "vertical")
}

# ---- Helper: standardized barplot from a GO/KEGG result table ----
#  - Horizontal bar plot (y = Description, x = Count)
#  - Bar fill colored by p.adjust gradient (same colors as dotplot)
#  - Capitalizes GO terms
#  - theme_poster()
make_poster_barplot <- function(table_path, title_text, max_terms = 15) {
  if (!file.exists(table_path)) {
    cat(sprintf("  SKIP (table not found): %s\n", basename(table_path)))
    return(NULL)
  }

  df <- read.csv(table_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    cat(sprintf("  SKIP (empty table): %s\n", basename(table_path)))
    return(NULL)
  }

  if (!("Count" %in% names(df))) {
    cat(sprintf("  SKIP (no Count column): %s\n", basename(table_path)))
    return(NULL)
  }

  # Capitalize GO terms
  if ("Description" %in% names(df)) {
    df$Description <- capitalize_first(df$Description)
  }

  # Sort by p.adjust (most significant first), keep top N
  df <- df[order(df$p.adjust), ]
  df <- head(df, max_terms)
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "#D55E00", high = "#0072B2",
                        name = "p-adjusted") +
    labs(
      x = "Gene Count",
      y = NULL,
      title = title_text
    ) +
    guides(
      fill = guide_colorbar(title.position = "top",
                            barwidth = 1.5, barheight = 8)
    ) +
    theme_poster() +
    theme(legend.position = "right")
}

# ---- Generate dotplots AND barplots for every available GO/KEGG table ----
# For each (pattern, output_prefix, title) entry below, we look up the most
# recent matching CSV in TABLE_DIR. If found, we create BOTH a dotplot and a
# barplot with theme_poster().
go_table_specs <- list(
  # TurboID targeted GO/KEGG
  list(pattern = "targeted_GO_turbo_trip4.*BP",   prefix = "go_bp_turbo",
       title = "GO Biological Process — TurboID TRIP4"),
  list(pattern = "targeted_GO_turbo_trip4.*MF",   prefix = "go_mf_turbo",
       title = "GO Molecular Function — TurboID TRIP4"),
  list(pattern = "targeted_GO_turbo_trip4.*CC",   prefix = "go_cc_turbo",
       title = "GO Cellular Component — TurboID TRIP4"),
  list(pattern = "targeted_KEGG_turbo_trip4",     prefix = "kegg_turbo",
       title = "KEGG Pathways — TurboID TRIP4"),

  # CHX-enriched GO
  list(pattern = "GO_CHX_enriched.*BP", prefix = "go_bp_chx_enriched",
       title = "GO BP — CHX Enriched"),
  list(pattern = "GO_CHX_enriched.*MF", prefix = "go_mf_chx_enriched",
       title = "GO MF — CHX Enriched"),
  list(pattern = "GO_CHX_enriched.*CC", prefix = "go_cc_chx_enriched",
       title = "GO CC — CHX Enriched"),

  # CHX-depleted GO (= DMSO-enriched side)
  list(pattern = "GO_CHX_depleted.*BP", prefix = "go_bp_chx_depleted",
       title = "GO BP — CHX Depleted (DMSO side)"),
  list(pattern = "GO_CHX_depleted.*MF", prefix = "go_mf_chx_depleted",
       title = "GO MF — CHX Depleted (DMSO side)"),
  list(pattern = "GO_CHX_depleted.*CC", prefix = "go_cc_chx_depleted",
       title = "GO CC — CHX Depleted (DMSO side)"),

  # Flag-IP validated GO (both C+N Flag)
  list(pattern = "flagip_GO_validated_both.*BP", prefix = "go_bp_validated",
       title = "GO BP — Validated by Both C-Flag + N-Flag"),
  list(pattern = "flagip_GO_validated_both.*MF", prefix = "go_mf_validated",
       title = "GO MF — Validated by Both C-Flag + N-Flag"),
  list(pattern = "flagip_GO_validated_both.*CC", prefix = "go_cc_validated",
       title = "GO CC — Validated by Both C-Flag + N-Flag"),

  # Any additional KEGG tables for CHX/DMSO comparisons
  list(pattern = "KEGG_CHX_enriched",   prefix = "kegg_chx_enriched",
       title = "KEGG Pathways — CHX Enriched"),
  list(pattern = "KEGG_CHX_depleted",   prefix = "kegg_chx_depleted",
       title = "KEGG Pathways — CHX Depleted (DMSO side)")
)

cat("\n[GO/KEGG] Generating dotplots and barplots...\n")
go_count <- 0
for (spec in go_table_specs) {
  table_file <- find_latest_table(spec$pattern)
  if (is.null(table_file)) {
    cat(sprintf("  SKIP (no table matching %s)\n", spec$pattern))
    next
  }

  # Dotplot
  p_dot <- make_poster_dotplot(table_file, spec$title)
  if (!is.null(p_dot)) {
    save_poster(p_dot, paste0(spec$prefix, "_dot"), width = 18, height = 14)
    go_count <- go_count + 1
  }

  # Barplot
  p_bar <- make_poster_barplot(table_file, spec$title)
  if (!is.null(p_bar)) {
    save_poster(p_bar, paste0(spec$prefix, "_bar"), width = 18, height = 14)
    go_count <- go_count + 1
  }
}
cat(sprintf("  Generated %d GO/KEGG dot/bar plots.\n", go_count))

#=============================================================================
# DONE
#=============================================================================
cat("\n=========================================\n")
cat(" Poster figures complete!\n")
cat(sprintf(" Output: %s/*.pdf\n", POSTER_FIG_DIR))
cat("=========================================\n")
cat("\nNext steps:\n")
cat("  1. Compile poster: cd poster && pdflatex poster_simple.tex\n")
cat("  2. Or use the full poster: cd poster && pdflatex poster.tex\n")
cat("  3. For STRING network figures: make lydia-volcano string-style-network\n")
cat("     then copy PDFs from output/figures/ to poster/figures/\n")
