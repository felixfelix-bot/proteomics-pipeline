###############################################################################
# 06_gene_families.R
# Highlights gene families (GPATCH, DHX, DDX, LARP) on volcano plots.
# Based on Lydia's 20260424_Overlap_TurboID_wOthers.R gene family highlighting.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/06_gene_families.R"
###############################################################################

# === WHAT ARE THESE GENE FAMILIES? ============================================
# GPATCH, DHX, DDX, and LARP are families of RNA-BINDING PROTEINS — proteins
# that interact with RNA molecules in the cell. Since TRIP4 is involved in
# RNA-related processes, we want to see if members of these families show up
# in our proteomics data.
#
#   GPATCH = proteins with a "G-patch" domain (found in RNA-processing proteins).
#            Members are listed explicitly in GENE_FAMILIES$GPATCH (from config).
#   DHX    = DEAH-box helicases (unwind RNA/DNA; e.g., DHX9, DHX30).
#            Identified by name prefix "DHX".
#   DDX    = DEAD-box helicases (another family of RNA unwinding enzymes;
#            e.g., DDX3X, DDX5). Identified by name prefix "DDX".
#   LARP   = La-related proteins (bind RNA, involved in translation/RNA stability;
#             e.g., LARP1, LARP4B). Identified by name prefix "LARP".
#
# Highlighting these on volcano plots shows whether known RNA-binding proteins
# are among our significant hits — a biological signal worth following up.
# =============================================================================

# Print a header to the console so you know which script is running:
cat("\n=========================================\n")
cat(" Gene Family Highlighting\n")
cat("=========================================\n\n")

# ---- Load data ----
# load_all_experiments() is a custom helper (from utils.R) that reads all
# proteomics result files. Returns a named list of data frames.
experiments <- load_all_experiments()

# ---- Helper: classify a gene into its family ----
# This function takes a single gene name (e.g., "DDX5") and returns its
# family abbreviation ("ddx"), or NA if it doesn't belong to any tracked family.
# return() immediately exits the function with the given value.
get_gene_family <- function(gene) {
  # %in% checks if the gene is in the explicit GPATCH list (from config):
  if (gene %in% GENE_FAMILIES$GPATCH) return("gp")
  # grepl() tests if a string matches a pattern (regular expression).
  # "^DHX" means "starts with DHX" (^ = start of string).
  # This catches DHX9, DHX30, DHX36, etc.
  if (grepl("^DHX", gene)) return("dhx")
  # Similarly for DDX family:
  if (grepl("^DDX", gene)) return("ddx")
  # And LARP family:
  if (grepl("^LARP", gene)) return("LARPs")
  # NA_character_ = "missing text value" (NA for character/string data).
  return(NA_character_)
}

# ---- Helper: volcano with gene families highlighted ----
# This function creates a volcano plot (fold change vs significance) where
# members of our gene families are colored and labeled distinctly.
#   df    = a data frame with gene, log2FC, padj columns
#   title = a string for the plot title
plot_gene_family_volcano <- function(df, title) {
  # Make a copy of the data to work with (avoid modifying the original):
  toPlot <- df

  # Classify each protein as significant ("TRUE") or not ("FALSE"):
  # P_VALUE_CUTOFF and LOG2FC_CUTOFF are defined in 01_config.R.
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )

  # sapply() applies a function to EACH ELEMENT of a vector, one at a time,
  # and returns a vector of results. Here, it calls get_gene_family() on
  # every gene name in toPlot$gene.
  # (Contrast with lapply(), which returns a list; sapply() simplifies to a vector.)
  fams <- sapply(toPlot$gene, get_gene_family)

  # !is.na(fams) identifies genes that belong to a tracked family:
  fam_mask <- !is.na(fams)

  # Override the category for family members with their family abbreviation.
  # This means family members get special colors regardless of significance.
  toPlot$category[fam_mask] <- fams[fam_mask]

  # factor() converts the category column to a categorical variable with
  # fixed levels. CATEGORY_COLORS (from config) defines the color for each category.
  # drop = FALSE ensures all levels appear in the legend even if some are absent.
  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  # Select only family-member rows for labeling (so they're easy to spot):
  label_data <- toPlot[!is.na(fams), ]

  # Build the volcano plot using ggplot2's layered grammar of graphics:
  # ggplot2::ggplot() initializes the plot.
  # aes() maps data columns to visual properties (x, y, color).
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    # geom_point() draws the scatter points (alpha = transparency, 0.3 = 30% opaque):
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    # scale_color_manual() assigns specific colors from CATEGORY_COLORS:
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    # ggrepel::geom_text_repel() adds text labels that automatically avoid
    # overlapping each other — much better than plain text labels:
    ggrepel::geom_text_repel(
      data = label_data,                           # Only label family members
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",               # Bold text for visibility
      max.overlaps = 25, show.legend = FALSE
    ) +
    # Dashed horizontal line at the p-value significance threshold:
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    # Dashed vertical lines at the fold change thresholds (both directions):
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    # Axis labels and plot title:
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = paste("Gene Families:", title)
    ) +
    # theme_bw() = clean black-and-white background:
    ggplot2::theme_bw() +
    # Customize text sizes for readability:
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8)
    )

  # Return the plot object (doesn't display it yet — caller decides to save):
  return(p)
}

# =====================================================================
# MAIN: Generate gene family volcanos for TurboID experiments
# =====================================================================

# grep() searches for pattern matches. Here we find experiment names that
# start with "turbo_" (the ^ means "start of string").
# value = TRUE returns the matched strings themselves (not their positions).
turbo_names <- grep("^turbo_", names(experiments), value = TRUE)

# Loop through each TurboID experiment:
for (name in turbo_names) {
  cat(sprintf("  [%s]\n", name))

  # Classify all genes in this experiment and count family members:
  # sapply() applies get_gene_family to each gene name.
  fams <- sapply(experiments[[name]]$gene, get_gene_family)

  # table() creates a FREQUENCY COUNT — it counts how many times each
  # unique value appears. Here, it counts genes per family.
  # Example: table(c("dhx","ddx","dhx")) → dhx=2, ddx=1
  # fams[!is.na(fams)] excludes genes that aren't in any family.
  fam_counts <- table(fams[!is.na(fams)])

  # If any family members were found, print a summary:
  # sprintf("%s=%d", ...) formats as "familyName=count".
  # paste(..., collapse=", ") joins them with commas.
  if (length(fam_counts) > 0) {
    cat(sprintf("    Gene families found: %s\n",
                paste(sprintf("%s=%d", names(fam_counts), fam_counts), collapse = ", ")))
  }

  # Create the volcano plot and save it:
  p <- plot_gene_family_volcano(experiments[[name]], name)
  save_figure(p, paste0("genefamily_", name), width = 7, height = 5)
}

# Also do Flag IP experiments (experiments starting with "flag_"):
flag_names <- grep("^flag_", names(experiments), value = TRUE)
for (name in flag_names) {
  p <- plot_gene_family_volcano(experiments[[name]], name)
  save_figure(p, paste0("genefamily_", name), width = 7, height = 5)
}

# ---- Summary table: which family members are significant where ----
# This builds a CROSS-EXPERIMENT summary: for each gene-family member, it shows
# whether that gene was significant (TRUE/FALSE) in each experiment.
# Think of it as a matrix: rows = genes, columns = experiments.
cat("\nGenerating gene family summary table...\n")

# Step 1: Collect all unique gene names across ALL experiments.
# lapply() applies a function to each experiment's data frame, extracting $gene.
# unlist() flattens the resulting list into one vector.
# unique() removes duplicates.
all_genes <- unique(unlist(lapply(experiments, function(df) df$gene)))

# Step 2: Keep only genes that belong to a tracked family:
# sapply() applies get_gene_family to each gene; !is.na() keeps family members.
fam_genes <- all_genes[!is.na(sapply(all_genes, get_gene_family))]

# Step 3: Build the summary table (only if we found family members):
if (length(fam_genes) > 0) {
  summary_list <- list()  # Initialize an empty list to collect rows

  # For each family-member gene:
  for (gene in fam_genes) {
    # Start building a row: gene name + its family:
    row <- list(gene = gene, family = get_gene_family(gene))

    # For each experiment, check if this gene is significant:
    for (exp_name in names(experiments)) {
      row[[exp_name]] <- FALSE  # Default: not significant
      df <- experiments[[exp_name]]

      # Find this gene in the experiment's data frame:
      # df[df$gene == gene, ] selects rows where gene matches.
      hit <- df[df$gene == gene, ]

      # If the gene exists in this experiment, check significance:
      if (nrow(hit) > 0) {
        # TRUE if both p-value and fold change pass thresholds:
        row[[exp_name]] <- hit$padj < P_VALUE_CUTOFF & abs(hit$log2FC) > LOG2FC_CUTOFF
      }
    }

    # Convert the row (a list) to a one-row data frame and store it:
    summary_list[[gene]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }

  # do.call(rbind, list) stacks all one-row data frames into one table.
  # rbind = "row bind" = stack rows vertically.
  fam_summary <- do.call(rbind, summary_list)

  # Sort the table: first by family, then by gene name (alphabetical):
  fam_summary <- fam_summary[order(fam_summary$family, fam_summary$gene), ]

  # Save the summary table (CSV/TSV file via the custom helper):
  save_table(fam_summary, "gene_family_significance_summary")
}

# Print completion message:
cat("\n=========================================\n")
cat(" Gene family highlighting complete!\n")
cat("=========================================\n")
