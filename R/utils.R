###############################################################################
# utils.R
# Shared utility functions used across all analysis scripts.
###############################################################################

# ---- Save a ggplot object to PNG and PDF ----
save_figure <- function(plot, filename, width = FIG_WIDTH, height = FIG_HEIGHT,
                        dpi = FIG_DPI) {
  filepath_png <- file.path(FIGURE_DIR, paste0(filename, ".png"))
  filepath_pdf <- file.path(FIGURE_DIR, paste0(filename, ".pdf"))

  ggplot2::ggsave(filepath_png, plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filepath_pdf, plot, width = width, height = height)

  cat(sprintf("  Saved: %s\n  Saved: %s\n", basename(filepath_png), basename(filepath_pdf)))
}

# ---- Save a data frame to CSV ----
save_table <- function(df, filename) {
  filepath <- file.path(TABLE_DIR, paste0(filename, ".csv"))
  write.csv(df, filepath, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", basename(filepath)))
  invisible(filepath)
}

# ---- Load known interactors from file ----
load_known_interactors <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("Known interactors file not found: ", filepath)
    return(character(0))
  }
  genes <- readLines(filepath)
  genes <- trimws(genes)
  genes <- genes[genes != "" & !grepl("^#", genes)]
  cat(sprintf("  Loaded %d known interactors\n", length(genes)))
  return(genes)
}

# ---- Create a labeled volcano plot using ggplot2 (for overlay or custom) ----
# This is the custom version for multi-experiment overlays.
make_volcano_ggplot <- function(df, title = "", label_genes = NULL,
                                experiment_name = NULL,
                                padj_cutoff = P_VALUE_CUTOFF,
                                log2fc_cutoff = LOG2FC_CUTOFF) {

  df$neglog10p <- -log10(df$padj)
  df$category <- with(df, ifelse(
    padj < padj_cutoff & abs(log2FC) > log2fc_cutoff,
    ifelse(log2FC > 0, "Up", "Down"),
    "NS"
  ))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggplot2::scale_color_manual(
      values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70"),
      labels = c("Down", "Not significant", "Up"),
      name = NULL
    ) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adjusted~italic(p)~value)),
      title = title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )

  # Label specific genes if requested
  if (!is.null(label_genes) && length(label_genes) > 0) {
    label_data <- df[df$gene %in% label_genes, ]
    if (nrow(label_data) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = gene),
          size = 3, fontface = "bold",
          max.overlaps = 25,
          show.legend = FALSE,
          color = "black"
        )
    }
  }

  return(p)
}

cat("[utils] Utility functions loaded.\n")
