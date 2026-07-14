###############################################################################
# 00_theme.R
# Global plot theme — STANDARDIZED FONT SIZES for poster-quality figures.
#
# PROBLEM: When figures are generated at different sizes (8x6 vs 10x7 vs 12x8)
# and then resized for a poster, font sizes become inconsistent.
#
# SOLUTION: This file defines a SINGLE theme applied to ALL plots. Every figure
# uses the same base font size regardless of output dimensions. When exported
# as PDF (vector format), text stays crisp at any poster size.
#
# USAGE: Source this file in every R script BEFORE creating plots:
#   source("R/00_theme.R")
#
# Then use:
#   my_plot + theme_poster()        # for standard poster figures
#   my_plot + theme_publication()   # for journal figures (slightly smaller)
#
# KEY: Export figures as PDF for posters. PDF = vector = no quality loss.
###############################################################################

# ---- Standardized font sizes (in points) ----
# These are the ONLY font size values used across the entire pipeline.
# Changing a number here changes it EVERYWHERE.
POSTER_FONTS <- list(
  title       = 20,
  subtitle    = 20,
  axis_title  = 20,
  axis_text   = 20,
  legend_title = 20,
  legend_text  = 20,
  facet_title  = 20,
  annotate    = 20
)

# Helper: create size breaks that cover the FULL data range (min to max).
# Produces n_breaks circles with CONSTANT step size from min to max.
# The smallest legend dot = actual data minimum (e.g. 1, never 0).
#
# CRITICAL: The caller MUST also pass limits = c(min(counts), max(counts))
# to scale_size_continuous() so the size domain matches the breaks exactly.
# Otherwise ggplot extends the domain below min and the dot sizes don't match.
#
# Usage:
#   cnt <- p$data$Count
#   scale_size_continuous(name = "Gene Count", range = c(3, 10),
#                         breaks = make_size_breaks(cnt, n_breaks = 8),
#                         limits = c(min(cnt), max(cnt)))
make_size_breaks <- function(counts, n_breaks = 8) {
  cmin <- min(counts, na.rm = TRUE)
  cmax <- max(counts, na.rm = TRUE)
  # Constant-step integer sequence from min to max
  unique(round(seq(cmin, cmax, length.out = n_breaks)))
}

# Standard guides() call for size legend — ensures legend circles
# render consistently (grey fill, not black) regardless of the
# color scale override applied to the dotplot.
# Usage: p + guides(size = size_legend_guide())
size_legend_guide <- function() {
  ggplot2::guide_legend(
    override.aes = list(color = "grey50", fill = NA),
    title.position = "top",
    order = 2
  )
}
# clusterProfiler returns terms like "negative regulation of..."
# This capitalizes the very first character.
# Shared across all GO scripts (R/11, R/23, R/30, R/31).
capitalize_first <- function(x) {
  sapply(x, function(s) {
    s <- as.character(s)
    if (nchar(s) > 0) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s))) else s
  })
}

# Wrap + capitalize labels for dotplot y-axis (GO terms are very long)
# Usage: scale_y_discrete(labels = wrap_go_labels)
wrap_go_labels <- function(x, width = 35) {
  stringr::str_wrap(capitalize_first(x), width = width)
}

# Publication fonts (slightly smaller for journal figures)
PUBLICATION_FONTS <- list(
  title       = 12,
  subtitle    = 10,
  axis_title  = 11,
  axis_text   = 9,
  legend_title = 10,
  legend_text  = 8,
  facet_title  = 10,
  annotate    = 7
)

# ---- Poster theme (for all poster figures) ----
# Applies consistent font sizes to every ggplot element.
# Usage: my_plot + theme_poster()              # default 20pt
#        my_plot + theme_poster(font_size = 18) # smaller (volcanos)
# ---- Poster theme (FINAL — based on Aruna's voice feedback Jul 15) ----
# Settings derived from style_variants review:
#   Grid: grey65, lw=0.8 (V3/V6 darkness level, between grey75 and grey50)
#   Font: 18pt axes/ticks/legend (V9 "close to perfect" at 20, V7 "good" at 20,
#         but split difference since 20pt protein labels were too big)
#   Titles: blank (LaTeX provides panel labels)
theme_poster <- function(base_family = "sans", font_size = 18) {
  ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title.x  = ggplot2::element_text(size = font_size,
                                             margin = ggplot2::margin(t = 6),
                                             color = "black"),
      axis.title.y  = ggplot2::element_text(size = font_size,
                                             margin = ggplot2::margin(r = 6),
                                             color = "black"),
      axis.text.x   = ggplot2::element_text(size = font_size, color = "black"),
      axis.text.y   = ggplot2::element_text(size = font_size, color = "black"),
      legend.title  = ggplot2::element_text(size = font_size, face = "bold",
                                            color = "black"),
      legend.text   = ggplot2::element_text(size = font_size,
                                            color = "black"),
      legend.position  = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.key       = ggplot2::element_rect(fill = "white", color = NA),
      strip.text    = ggplot2::element_text(size = font_size, face = "bold",
                                            color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
      # Grid: grey65, lw=0.8 (Aruna praised V3/V6 darkness)
      panel.grid.major = ggplot2::element_line(color = "grey65", linewidth = 0.8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, color = "grey40", linewidth = 0.6),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

# ---- Publication theme (for journal figures) ----
theme_publication <- function(base_family = "sans") {
  ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = PUBLICATION_FONTS$title,
                                             hjust = 0.5, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(size = PUBLICATION_FONTS$subtitle, hjust = 0.5),
      axis.title.x  = ggplot2::element_text(size = PUBLICATION_FONTS$axis_title),
      axis.title.y  = ggplot2::element_text(size = PUBLICATION_FONTS$axis_title),
      axis.text.x   = ggplot2::element_text(size = PUBLICATION_FONTS$axis_text, color = "black"),
      axis.text.y   = ggplot2::element_text(size = PUBLICATION_FONTS$axis_text, color = "black"),
      legend.title  = ggplot2::element_text(size = PUBLICATION_FONTS$legend_title, face = "bold"),
      legend.text   = ggplot2::element_text(size = PUBLICATION_FONTS$legend_text),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      strip.text    = ggplot2::element_text(size = PUBLICATION_FONTS$facet_title, face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, color = "grey70", linewidth = 0.3),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

# ---- Helper: Save figure in BOTH PNG and PDF (for posters) ----
# PDF = vector format, text stays sharp at ANY size in LaTeX.
# PNG = backup for software that can't handle PDF.
#
# Usage: save_poster_figure(my_plot, "volcano_main", width=8, height=6)
# Creates: output/figures/volcano_main.pdf + .png + .svg
save_poster_figure <- function(plot, name, width = 8, height = 6, dpi = 300) {
  # PDF (vector — for LaTeX poster)
  pdf_path <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  device = grDevices::cairo_pdf)

  # PNG (raster — backup)
  png_path <- file.path(FIGURE_DIR, paste0(name, ".png"))
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")

  # SVG (vector — editable in Illustrator/Inkscape)
  svg_path <- file.path(FIGURE_DIR, paste0(name, ".svg"))
  ggplot2::ggsave(svg_path, plot = plot, width = width, height = height,
                  bg = "white")

  cat(sprintf("  Saved: %s (.pdf + .png + .svg)\n", name))
}

# ---- Standardized figure dimensions ----
# Use these constants instead of ad-hoc width/height values.
# This ensures all figures have consistent aspect ratios.
FIG_DIMS <- list(
  single_col = c(width = 6, height = 4.5),     # Single column in poster
  double_col = c(width = 12, height = 4.5),     # Full width
  square     = c(width = 7, height = 7),        # Square (network plots)
  tall       = c(width = 7, height = 10),       # Tall (many GO terms)
  wide       = c(width = 12, height = 6),       # Wide (volcano + Venn side-by-side)
  poster_full = c(width = 11, height = 8.5)     # Full poster panel
)

cat("[Theme] Standardized font theme loaded. Use theme_poster() or theme_publication().\n")
cat("         Export with save_poster_figure() for PDF+PNG vector output.\n")
