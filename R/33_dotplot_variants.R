###############################################################################
# 33_dotplot_variants.R
# Regenerates the 3 key GO BP dotplots with different font sizes and label
# treatments — directly from the saved GO enrichment CSV tables.
#
# NO need to re-run enrichGO. Reads the tables that were already saved by:
#   R/30_flagip_validated_go.R  → flagip_GO_validated_both_BP
#   R/26_network_go_comparison.R → network_go_in_network_BP
#   R/11_targeted_go.R           → targeted_GO_TurboID_TRIP4_vs_WT_BP
#
# For EACH of the 3 plots, generates TWO sets of 3 font sizes:
#
#   SET 1 — Full GO terms (with y-axis text wrapping):
#     <name>_full_terms_tickSIZE20_xlabelSIZE20.png
#     <name>_full_terms_tickSIZE24_xlabelSIZE24.png
#     <name>_full_terms_tickSIZE28_xlabelSIZE28.png
#
#   SET 2 — Short GO terms only (max 3 words):
#     <name>_short_terms_tickSIZE20_xlabelSIZE20.png
#     <name>_short_terms_tickSIZE24_xlabelSIZE24.png
#     <name>_short_terms_tickSIZE28_xlabelSIZE28.png
#
# Total: 3 plots × 2 variants × 3 sizes = 18 PNG files.
#
# Usage:
#   make dotplot-variants
#
# Prerequisite: Run make targeted-go, make network-go, make flagip-validated-go
# at least once so the CSV tables exist in output/tables/.
###############################################################################

cat("\n=========================================\n")
cat(" GO Dotplot Font-Size Variants\n")
cat("=========================================\n\n")

library(ggplot2)
library(stringr)

source("R/00_theme.R")

# ---- Configuration ----
# Font sizes: 18, 20, 24.
FONT_SIZES <- c(18, 20, 24)
# Y-axis label wrapping width (characters per line)
WRAP_WIDTH <- 45
# Max words for "short" GO terms
MAX_SHORT_WORDS <- 3

# ---- Plot definitions ----
# Each entry maps to a CSV table in output/tables/
PLOT_DEFS <- list(
  list(
    table_pattern = "flagip_GO_validated_both_BP",
    plot_prefix   = "flagip_GO_validated_both_BP_dotplot",
    title         = "Flag IP Validated (Both) — Biological Process"
  ),
  list(
    table_pattern = "network_go_in_network_BP",
    plot_prefix   = "network_go_comparison_BP",
    title         = "GO Enrichment: In-Network — Biological Process"
  ),
  list(
    table_pattern = "targeted_GO_TurboID_TRIP4_vs_WT_BP",
    plot_prefix   = "targeted_GO_TurboID_TRIP4_vs_WT_BP_dotplot",
    title         = "GO analysis of TRIP4 TurboID vs Wild Type — Biological Process"
  )
)

# ---- Helper: find the most recent CSV table matching a pattern ----
find_table <- function(pattern) {
  # CSV files are named: <pattern>_<commithash>.csv
  files <- list.files(TABLE_DIR, pattern = paste0("^", pattern, "_.*\\.csv$"),
                      full.names = TRUE)
  if (length(files) == 0) {
    cat(sprintf("  WARNING: No CSV table matching '%s' in %s\n", pattern, TABLE_DIR))
    return(NULL)
  }
  # Use the most recently modified file
  most_recent <- files[order(file.info(files)$mtime, decreasing = TRUE)][1]
  cat(sprintf("  Using: %s\n", basename(most_recent)))
  most_recent
}

# ---- Helper: parse GeneRatio string ("5/120") to numeric ----
parse_gene_ratio <- function(ratio_str) {
  parts <- strsplit(as.character(ratio_str), "/")
  sapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]))
}

# ---- Helper: count words in a string ----
word_count <- function(x) {
  sapply(strsplit(as.character(x), "\\s+"), length)
}

# ---- Helper: wrap long labels for y-axis ----
wrap_label <- function(labels, width = WRAP_WIDTH) {
  sapply(labels, function(lab) {
    lab <- as.character(lab)
    # Capitalize first letter
    lab <- paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab)))
    # Wrap if longer than width
    str_wrap(lab, width = width)
  })
}

# ---- Helper: build a dotplot from a data frame ----
# Recreates what enrichplot::dotplot does, but from raw CSV data.
make_dotplot <- function(df, font_size, use_short_terms, title_text) {
  # Filter to short terms if requested
  if (use_short_terms) {
    df <- df[word_count(df$Description) <= MAX_SHORT_WORDS, ]
    if (nrow(df) == 0) {
      cat("    WARNING: No short GO terms found. Skipping.\n")
      return(NULL)
    }
  }

  # Parse GeneRatio from "k/N" to numeric
  df$GeneRatioNum <- parse_gene_ratio(df$GeneRatio)

  # Sort by GeneRatio descending (same as enrichplot default)
  df <- df[order(df$GeneRatioNum, decreasing = TRUE), ]

  # Use Description as factor, ordered by GeneRatio (bottom = highest ratio)
  df$Description <- factor(df$Description, levels = rev(df$Description))

  # Label wrapping function for y-axis
  label_func <- if (use_short_terms) {
    function(labels) capitalize_first(labels)
  } else {
    function(labels) wrap_label(labels, width = WRAP_WIDTH)
  }

  # Build the plot
  cnt <- df$Count
  size_breaks <- make_size_breaks(cnt, n_breaks = 8)

  p <- ggplot(df, aes(x = GeneRatioNum, y = Description,
                       color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#D55E00", high = "#0072B2",
                         name = "p-adjusted value") +
    scale_size_continuous(name = "Gene Count", range = c(3, 10),
                          breaks = size_breaks,
                          limits = c(min(cnt), max(cnt))) +
    scale_y_discrete(labels = label_func) +
    guides(size = size_legend_guide()) +
    labs(title = title_text, x = "Gene Ratio") +
    theme_poster(font_size = font_size) +
    theme(
      axis.text = element_text(size = font_size, color = "black"),
      axis.title.x = element_text(size = font_size, color = "black"),
      plot.title = element_text(size = font_size, face = "bold",
                                 hjust = 0.5, color = "black"),
      panel.grid.major = element_line(color = "grey75", linewidth = 0.4)
    )

  return(p)
}

# ---- Helper: save a plot with descriptive filename ----
save_variant <- function(plot, plot_prefix, variant_label, font_size, n_terms) {
  # Filename format: <prefix>_<variant>_tickSIZE<size>_xlabelSIZE<size>
  filename <- sprintf("%s_%s_tickSIZE%d_xlabelSIZE%d",
                      plot_prefix, variant_label, font_size, font_size)
  fig_height <- max(14, n_terms * 0.8)
  save_figure(plot, filename, width = 18, height = fig_height)
}

# =====================================================================
# Main: generate all 18 variants
# =====================================================================

for (def in PLOT_DEFS) {
  cat(sprintf("\n--- Processing: %s ---\n", def$plot_prefix))

  csv_path <- find_table(def$table_pattern)
  if (is.null(csv_path)) {
    cat("  Skipping (no data).\n")
    next
  }

  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d GO terms\n", nrow(df)))

  # ---- SET 1: Full GO terms (with wrapping) ----
  cat("  Generating full-term variants (3 sizes)...\n")
  for (sz in FONT_SIZES) {
    p <- make_dotplot(df, font_size = sz, use_short_terms = FALSE,
                      title_text = def$title)
    if (!is.null(p)) {
      save_variant(p, def$plot_prefix, "full_terms", sz, nrow(df))
    }
  }

  # ---- SET 2: Short GO terms only (max 3 words) ----
  short_df <- df[word_count(df$Description) <= MAX_SHORT_WORDS, ]
  cat(sprintf("  Short GO terms (<= %d words): %d of %d\n",
              MAX_SHORT_WORDS, nrow(short_df), nrow(df)))

  if (nrow(short_df) > 0) {
    cat("  Generating short-term variants (3 sizes)...\n")
    for (sz in FONT_SIZES) {
      p <- make_dotplot(df, font_size = sz, use_short_terms = TRUE,
                        title_text = def$title)
      if (!is.null(p)) {
        save_variant(p, def$plot_prefix, "short_terms_max3words", sz,
                     nrow(short_df))
      }
    }
  } else {
    cat("  No short GO terms available. Skipping short-term variants.\n")
  }
}

cat("\n=========================================\n")
cat(" Dotplot variants complete!\n")
cat(sprintf(" Output: %s\n", FIGURE_DIR))
cat(" Total: 18 PNG + PDF files (3 plots x 2 variants x 3 sizes)\n")
cat("=========================================\n")
