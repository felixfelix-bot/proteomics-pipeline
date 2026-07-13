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
  title       = 24,
  subtitle    = 24,
  axis_title  = 24,
  axis_text   = 24,
  legend_title = 24,
  legend_text  = 24,
  facet_title  = 24,
  annotate    = 24
)

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
# Usage: my_plot + theme_poster()
theme_poster <- function(base_family = "Arial") {
  ggplot2::theme_minimal(base_family = base_family) +
    ggplot2::theme(
      # Titles
      plot.title    = ggplot2::element_text(face = "bold", size = POSTER_FONTS$title,
                                             hjust = 0.5, margin = ggplot2::margin(b = 8),
                                             color = "black"),
      plot.subtitle = ggplot2::element_text(size = POSTER_FONTS$subtitle,
                                             hjust = 0.5, margin = ggplot2::margin(b = 6),
                                             color = "black"),
      # Axes
      axis.title.x  = ggplot2::element_text(size = POSTER_FONTS$axis_title,
                                             margin = ggplot2::margin(t = 6),
                                             color = "black"),
      axis.title.y  = ggplot2::element_text(size = POSTER_FONTS$axis_title,
                                             margin = ggplot2::margin(r = 6),
                                             color = "black"),
      axis.text.x   = ggplot2::element_text(size = POSTER_FONTS$axis_text, color = "black"),
      axis.text.y   = ggplot2::element_text(size = POSTER_FONTS$axis_text, color = "black"),
      # Legend
      legend.title  = ggplot2::element_text(size = POSTER_FONTS$legend_title, face = "bold",
                                            color = "black"),
      legend.text   = ggplot2::element_text(size = POSTER_FONTS$legend_text,
                                            color = "black"),
      legend.position  = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3),
      legend.key       = ggplot2::element_rect(fill = "white", color = NA),
      # Facets
      strip.text    = ggplot2::element_text(size = POSTER_FONTS$facet_title, face = "bold",
                                            color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
      # Panel
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(fill = NA, color = "grey70", linewidth = 0.3),
      # Margins (compact for poster density)
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
# Creates: output/figures/volcano_main.pdf AND .png
save_poster_figure <- function(plot, name, width = 8, height = 6, dpi = 300) {
  # PDF (vector — for LaTeX poster)
  pdf_path <- file.path(FIGURE_DIR, paste0(name, ".pdf"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height,
                  device = grDevices::cairo_pdf)

  # PNG (raster — backup)
  png_path <- file.path(FIGURE_DIR, paste0(name, ".png"))
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")

  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
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
