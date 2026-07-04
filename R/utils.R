###############################################################################
# utils.R
# Shared utility functions used across all analysis scripts.
#
# This file defines helper functions that EVERY analysis script needs:
#   - save_figure():      Save a plot as PNG + PDF
#   - save_table():       Save a data table as CSV
#   - load_known_interactors(): Read the known interactors text file
#   - discover_diffex_csvs():   Find all mass spec data files recursively
#   - load_all_experiments():   Load all data files into memory
#   - make_volcano_ggplot():    Create a basic volcano plot (used by some scripts)
###############################################################################

# ---- Save a ggplot object to PNG and PDF ----
# ggplot is R's main plotting system. You build a plot by adding layers
# with the + operator, then save it with ggsave().
#
# We save BOTH PNG (for viewing) and PDF (for publication-quality figures).
# The function takes:
#   plot:     The ggplot object to save
#   filename: Name WITHOUT extension (we add .png and .pdf)
#   width/height: Figure size in inches (default from config)
#   dpi:      Resolution (300 = publication quality)
save_figure <- function(plot, filename, width = FIG_WIDTH, height = FIG_HEIGHT,
                        dpi = FIG_DPI) {
  # file.path() joins path segments correctly for the current OS
  # paste0() concatenates strings without spaces
  filepath_png <- file.path(FIGURE_DIR, paste0(filename, ".png"))
  filepath_pdf <- file.path(FIGURE_DIR, paste0(filename, ".pdf"))

  # ggsave() writes the plot to a file. It auto-detects format from extension.
  # The :: syntax means "use ggsave from the ggplot2 package."
  ggplot2::ggsave(filepath_png, plot, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filepath_pdf, plot, width = width, height = height)

  # basename() strips the directory path, showing just the filename
  cat(sprintf("  Saved: %s\n  Saved: %s\n", basename(filepath_png), basename(filepath_pdf)))
}

# ---- Save a data frame to CSV ----
# Used to save extracted gene lists, GO enrichment results, etc.
save_table <- function(df, filename) {
  filepath <- file.path(TABLE_DIR, paste0(filename, ".csv"))
  # write.csv() exports a data frame to CSV format.
  # row.names = FALSE: don't add an extra column with row numbers.
  write.csv(df, filepath, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", basename(filepath)))
  # invisible() returns a value without printing it.
  # This is a convention in R — functions that produce "side effects"
  # (like writing a file) return their result invisibly so they don't
  # clutter the console when called.
  invisible(filepath)
}

# ---- Load known interactors from text file ----
# Reads data/known_interactors.txt — a plain text file with one gene name
# per line. These are proteins already known to interact with TRIP4
# (from literature or previous experiments).
#
# They get highlighted specially on volcano plots (orange color).
load_known_interactors <- function(filepath) {
  # Check if file exists; if not, warn and return empty list
  if (!file.exists(filepath)) {
    # warning() prints a warning message but does NOT stop execution
    # (unlike stop() which does)
    warning("Known interactors file not found: ", filepath)
    return(character(0))   # character(0) = empty character vector
  }
  # readLines() reads a text file into a character vector (one string per line)
  genes <- readLines(filepath)
  # trimws() removes leading/trailing whitespace from each line
  genes <- trimws(genes)
  # Filter: keep only non-empty lines that don't start with # (comment lines)
  # grepl("^#", x) returns TRUE if x starts with #
  # ! negates it: keep lines that do NOT start with #
  genes <- genes[genes != "" & !grepl("^#", genes)]
  cat(sprintf("  Loaded %d known interactors\n", length(genes)))
  return(genes)
}

# ---- Create a basic labeled volcano plot using ggplot2 ----
# This is the SIMPLE version (Up/Down/NS coloring only).
# The Lydia-style version with multi-category highlighting is in
# 02_volcano_plots.R — that one adds known interactors, Flag IP hits,
# gene families, etc.
#
# This simpler version is used for overlay plots where we compare
# two experiments on the same axes.
#
# A volcano plot has:
#   X-axis: log2 fold change (how much protein abundance changed)
#   Y-axis: -log10(adjusted p-value) (statistical confidence)
#   Each dot = one protein
make_volcano_ggplot <- function(df, title = "", label_genes = NULL,
                                experiment_name = NULL,
                                padj_cutoff = P_VALUE_CUTOFF,
                                log2fc_cutoff = LOG2FC_CUTOFF) {

  # Calculate -log10(padj) for the Y-axis.
  # The $ operator accesses a named column: df$padj means "the padj column
  # from the data frame df."
  # -log10() transforms very small p-values (0.00001) into large numbers (5),
  # making them easier to see on the plot.
  df$neglog10p <- -log10(df$padj)

  # Classify each protein into "Up" (significantly increased), "Down"
  # (significantly decreased), or "NS" (not significant).
  # with() evaluates an expression in the context of the data frame,
  # so we can write padj instead of df$padj.
  # ifelse() is vectorized: it applies to every row at once.
  df$category <- with(df, ifelse(
    padj < padj_cutoff & abs(log2FC) > log2fc_cutoff,
    ifelse(log2FC > 0, "Up", "Down"),  # Positive log2FC = Up, negative = Down
    "NS"                                # Failed thresholds = Not Significant
  ))

  # Build the plot layer by layer using ggplot2's + syntax:
  #
  # ggplot(data, aes(...)): Initialize with data and "aesthetic mappings"
  #   aes() maps data columns to visual properties:
  #   - x = log2FC → X-axis position
  #   - y = neglog10p → Y-axis position
  #   - color = category → dot color determined by Up/Down/NS
  #
  # geom_point(): Draw dots (one per protein)
  #   alpha = 0.5: semi-transparent (so overlapping dots are visible)
  #   size = 1.5: dot size
  #
  # scale_color_manual(): Define custom colors for each category
  #
  # geom_hline/vline: Draw threshold lines (dashed, grey)
  #   hline = horizontal line at y = -log10(0.05) ≈ 1.3
  #   vline = vertical lines at x = ±0.5 (fold change threshold)
  #
  # labs(): Axis labels and title
  #   expression() creates mathematical notation (subscripts, italics)
  #
  # theme_classic(): Clean white background with only axis lines
  # theme(): Fine-tune appearance (centered bold title, legend position)
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

  # Add text labels for specific genes (e.g., known interactors)
  # ggrepel prevents labels from overlapping each other.
  # Only label genes that are in the label_genes list AND present in the data.
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

# ---- Discover all diffEx_minProb.csv files recursively ----
# This function scans the data/ directory (and ALL its subdirectories)
# for files ending in "_diffEx_minProb.csv". These are the standard
# output from the mass spectrometry differential expression analysis.
#
# The function:
#   1. Finds all matching files recursively
#   2. Strips the "_diffEx_minProb" suffix to get clean experiment names
#   3. Handles duplicates (loads BOTH with __1, __2 suffixes)
#   4. Returns a NAMED CHARACTER VECTOR: names → file paths
discover_diffex_csvs <- function(data_dir = DATA_DIR) {
  # list.files() lists files in a directory.
  #   pattern: regular expression — "_diffEx_minProb\\.csv$" means
  #     "ends with _diffEx_minProb.csv" (\\ escapes the dot, $ = end)
  #   full.names = TRUE: return full paths, not just filenames
  #   recursive = TRUE: search inside ALL subdirectories
  #   ignore.case = TRUE: match regardless of case
  all_csvs <- list.files(data_dir, pattern = "_diffEx_minProb\\.csv$",
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE)

  # If no files found, print a helpful error and stop
  if (length(all_csvs) == 0) {
    cat("\n  ========================================\n")
    cat("  DATA ERROR\n")
    cat("  ========================================\n")
    cat("\n  No *_diffEx_minProb.csv files found in:\n")
    cat(sprintf("    %s\n", data_dir))
    cat("\n  The pipeline expects mass spectrometry output files\n")
    cat("  ending in '_diffEx_minProb.csv'.\n")
    cat("\n  Expected directory structure:\n")
    cat("    data/\n")
    cat("      <experiment_folder>/\n")
    cat("        <subfolder>/\n")
    cat("          *_diffEx_minProb.csv   <-- these are loaded\n\n")
    stop("No diffEx_minProb.csv files found in data/")
  }

  # Extract experiment names from filenames.
  # basename() gets just the filename (removes directory path).
  # gsub() replaces text using a pattern: strip "_diffEx_minProb.csv"
  # from the end, leaving just the experiment name.
  # Example: "BK467_TRIP4_vs_BK467_WT_diffEx_minProb.csv"
  #       → "BK467_TRIP4_vs_BK467_WT"
  exp_names <- basename(all_csvs)
  exp_names <- gsub("_diffEx_minProb\\.csv$", "", exp_names, ignore.case = TRUE)

  # Handle duplicates (same experiment name found in different folders).
  # Example: BK516_Cflag_vs_BK516_Ctrl might exist in two folders
  # (one normal, one with "(n=5)" suffix). We load BOTH and disambiguate.
  dup <- duplicated(exp_names)
  if (any(dup)) {
    dup_names <- unique(exp_names[dup])
    cat("  [INFO] Duplicate experiment names found - loading ALL copies:\n")
    for (d in dup_names) {
      paths <- all_csvs[exp_names == d]
      cat(sprintf("    '%s' appears %d times:\n", d, length(paths)))
      for (i in seq_along(paths)) {
        cat(sprintf("      [%d] %s\n", i, paths[i]))
      }
    }
    cat("  (Each will get a unique name: _1, _2, etc.)\n\n")

    # Rename ALL copies: append __1, __2, __3 to disambiguate
    for (d in dup_names) {
      idx <- which(exp_names == d)
      for (i in seq_along(idx)) {
        # paste0() concatenates strings without spaces
        exp_names[idx[i]] <- paste0(exp_names[idx[i]], "__", i)
      }
    }
  }

  # Assign the experiment names as the names of the character vector.
  # After this: all_csvs["BK467_TRIP4_vs_BK467_WT"] returns the full path.
  names(all_csvs) <- exp_names

  # Print a summary of what was found
  cat(sprintf("  Found %d experiment files:\n", length(all_csvs)))
  for (n in exp_names) {
    cat(sprintf("    - %s\n", n))
  }
  cat("\n")

  return(all_csvs)
}

# ---- Load all experiments from data directory ----
# Convenience function: discovers all CSV files, then loads each one.
# Returns a NAMED LIST where each element is a data frame:
#   experiments$BK467_TRIP4_vs_BK467_WT  → data frame with 3000 proteins
#   experiments$BK516_Cflag_vs_BK516_Ctrl → data frame with 3000 proteins
#   ...
#
# list() creates an empty list. Lists in R can hold any type of object,
# including other lists and data frames. Think of it like a Python dict
# or a JavaScript object.
load_all_experiments <- function(data_dir = DATA_DIR) {
  # Step 1: Find all CSV files
  csv_paths <- discover_diffex_csvs(data_dir)

  # Step 2: Load each one into the list
  experiments <- list()
  for (name in names(csv_paths)) {
    # Load the CSV using the function from 01_config.R
    experiments[[name]] <- load_proteomics_csv(csv_paths[[name]])
  }

  cat(sprintf("  Loaded %d experiments successfully.\n\n", length(experiments)))
  return(experiments)
}

# This line runs when the file is sourced — confirms it loaded successfully.
cat("[utils] Utility functions loaded.\n")
