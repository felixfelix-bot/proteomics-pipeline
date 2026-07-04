###############################################################################
# 06_gene_families.R
# Highlights gene families (GPATCH, DHX, DDX, LARP) on volcano plots.
# Based on Lydia's 20260424_Overlap_TurboID_wOthers.R gene family highlighting.
#
# Usage:
#   source("R/01_config.R")
#   source("R/utils.R")
#   source("R/06_gene_families.R")
###############################################################################

cat("\n=========================================\n")
cat(" Gene Family Highlighting\n")
cat("=========================================\n\n")

# ---- Load data ----
experiments <- list()
csv_files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)
for (f in csv_files) {
  name <- tools::file_path_sans_ext(basename(f))
  if (name == "known_interactors") next
  experiments[[name]] <- load_proteomics_csv(f)
}

# ---- Helper: classify a gene into its family ----
get_gene_family <- function(gene) {
  if (gene %in% GENE_FAMILIES$GPATCH) return("gp")
  if (grepl("^DHX", gene)) return("dhx")
  if (grepl("^DDX", gene)) return("ddx")
  if (grepl("^LARP", gene)) return("LARPs")
  return(NA_character_)
}

# ---- Helper: volcano with gene families highlighted ----
plot_gene_family_volcano <- function(df, title) {
  toPlot <- df
  toPlot$category <- ifelse(
    toPlot$padj < P_VALUE_CUTOFF & abs(toPlot$log2FC) > LOG2FC_CUTOFF,
    "TRUE", "FALSE"
  )

  # Assign gene families (overrides sig category)
  fams <- sapply(toPlot$gene, get_gene_family)
  fam_mask <- !is.na(fams)
  toPlot$category[fam_mask] <- fams[fam_mask]

  toPlot$category <- factor(toPlot$category, levels = names(CATEGORY_COLORS))

  label_data <- toPlot[!is.na(fams), ]

  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x = log2FC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.3, size = 1.2) +
    ggplot2::scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 2.5, fontface = "bold",
      max.overlaps = 25, show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~(adj.~italic(p)~value)),
      title = paste("Gene Families:", title)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8)
    )

  return(p)
}

# =====================================================================
# MAIN: Generate gene family volcanos for TurboID experiments
# =====================================================================
turbo_names <- grep("^turbo_", names(experiments), value = TRUE)

for (name in turbo_names) {
  cat(sprintf("  [%s]\n", name))

  # Count family members present
  fams <- sapply(experiments[[name]]$gene, get_gene_family)
  fam_counts <- table(fams[!is.na(fams)])
  if (length(fam_counts) > 0) {
    cat(sprintf("    Gene families found: %s\n",
                paste(sprintf("%s=%d", names(fam_counts), fam_counts), collapse = ", ")))
  }

  p <- plot_gene_family_volcano(experiments[[name]], name)
  save_figure(p, paste0("genefamily_", name), width = 7, height = 5)
}

# Also do Flag IP experiments
flag_names <- grep("^flag_", names(experiments), value = TRUE)
for (name in flag_names) {
  p <- plot_gene_family_volcano(experiments[[name]], name)
  save_figure(p, paste0("genefamily_", name), width = 7, height = 5)
}

# ---- Summary table: which family members are significant where ----
cat("\nGenerating gene family summary table...\n")

all_genes <- unique(unlist(lapply(experiments, function(df) df$gene)))
fam_genes <- all_genes[!is.na(sapply(all_genes, get_gene_family))]

if (length(fam_genes) > 0) {
  summary_list <- list()
  for (gene in fam_genes) {
    row <- list(gene = gene, family = get_gene_family(gene))
    for (exp_name in names(experiments)) {
      row[[exp_name]] <- FALSE
      df <- experiments[[exp_name]]
      hit <- df[df$gene == gene, ]
      if (nrow(hit) > 0) {
        row[[exp_name]] <- hit$padj < P_VALUE_CUTOFF & abs(hit$log2FC) > LOG2FC_CUTOFF
      }
    }
    summary_list[[gene]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  fam_summary <- do.call(rbind, summary_list)
  fam_summary <- fam_summary[order(fam_summary$family, fam_summary$gene), ]
  save_table(fam_summary, "gene_family_significance_summary")
}

cat("\n=========================================\n")
cat(" Gene family highlighting complete!\n")
cat("=========================================\n")
