###############################################################################
# 29_chx_volcano_venn.R
# Volcano plot + Venn diagram for CHX vs DMSO comparison.
#
# Volcano: CHX vs DMSO, enriched (orange) vs depleted (blue) vs non-sig (grey)
# Venn: CHX vs DMSO enriched vs depleted overlap
#
# Usage:
#   make chx-volcano-venn
###############################################################################

cat("\n=========================================\n")
cat(" CHX vs DMSO — Volcano + Venn\n")
cat("=========================================\n\n")

library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(grid)

experiments <- load_all_experiments()

CHX_EXP <- "TRIP4_CHX_vs_TRIP4_DMSO"
df <- find_experiment(experiments, CHX_EXP)

if (is.null(df)) {
  cat("ERROR:", CHX_EXP, "not found.\n")
  quit(status = 1)
}

df$neglog10p <- -log10(df$padj)
cat(sprintf("  Loaded %s: %d proteins\n", CHX_EXP, nrow(df)))

# ---- Classification ----
df$category <- "Not significant"
df$category[df$padj < P_VALUE_CUTOFF & df$log2FC >= 1] <- "Enriched in CHX"
df$category[df$padj < P_VALUE_CUTOFF & df$log2FC <= -1] <- "Enriched in DMSO"

cat(sprintf("  Enriched in CHX:  %d\n", sum(df$category == "Enriched in CHX")))
cat(sprintf("  Enriched in DMSO: %d\n", sum(df$category == "Enriched in DMSO")))
cat(sprintf("  Not significant:  %d\n", sum(df$category == "Not significant")))

# ---- Known interactors for labeling ----
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known_interactors <- load_known_interactors(interactors_file)
ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
label_genes <- unique(c(known_interactors, ASCC_CORE))

label_data <- df[df$gene %in% label_genes &
                 df$category != "Not significant" &
                 !is.na(df$gene), ]

# =====================================================================
# VOLCANO PLOT
# =====================================================================
cat("\n--- Generating volcano plot ---\n")

VOLCANO_COLORS <- c(
  "Enriched in CHX"  = "#D55E00",   # Vermillion orange
  "Enriched in DMSO" = "#0072B2",   # Navy blue
  "Not significant"  = "grey70"
)

p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
  ggplot2::geom_point(alpha = 0.5, size = 1) +
  ggplot2::scale_color_manual(values = VOLCANO_COLORS, name = NULL) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 2.5, fontface = "bold",
    max.overlaps = 30, show.legend = FALSE,
    bg.color = "white", bg.r = 0.15
  ) +
  ggplot2::geom_hline(yintercept = -log10(P_VALUE_CUTOFF),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = c(-1, 1),
                      linetype = "dashed", color = "grey50", linewidth = 0.3) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID: CHX vs DMSO",
    caption = sprintf("CHX-enriched: %d | DMSO-enriched: %d | n=%d",
                      sum(df$category == "Enriched in CHX"),
                      sum(df$category == "Enriched in DMSO"),
                      nrow(df))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 11),
    plot.caption = ggplot2::element_text(hjust = 0.5, size = 8, color = "grey30"),
    axis.text = ggplot2::element_text(colour = "black", size = 8),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.text = ggplot2::element_text(size = 8),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")
  )

save_figure(p, "chx_vs_dmso_volcano", width = 8, height = 6)

# =====================================================================
# VENN DIAGRAM: CHX-enriched vs DMSO-enriched
# =====================================================================
cat("\n--- Generating Venn diagram ---\n")

chx_genes <- df$gene[df$category == "Enriched in CHX" & !is.na(df$gene)]
dmso_genes <- df$gene[df$category == "Enriched in DMSO" & !is.na(df$gene)]

overlap <- length(intersect(chx_genes, dmso_genes))
cat(sprintf("  CHX-enriched: %d | DMSO-enriched: %d | Overlap: %d\n",
            length(chx_genes), length(dmso_genes), overlap))

vp <- draw.pairwise.venn(
  area1     = length(chx_genes),
  area2     = length(dmso_genes),
  cross.area = overlap,
  category  = c("Enriched in CHX", "Enriched in DMSO"),
  fill      = c("#D55E00", "#0072B2"),
  alpha     = rep(0.5, 2),
  cat.cex   = 1.4,
  cex       = 1.8,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  col       = "transparent",
  cat.pos   = c(-30, 30),
  cat.dist  = c(0.06, 0.06),
  margin    = 0.08,
  ind       = FALSE
)

commit_hash <- get_git_hash()
safe_name <- sanitize_filename("chx_vs_dmso_venn")
versioned <- paste0(safe_name, "_", commit_hash)
png_path <- safe_filepath(FIGURE_DIR, versioned, ".png")
pdf_path <- safe_filepath(FIGURE_DIR, versioned, ".pdf")

# Draw to PNG
grDevices::png(png_path, width = 7, height = 6, units = "in", res = FIG_DPI)
grid.draw(vp)
pushViewport(viewport())
grid.text("CHX vs DMSO Overlap", 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

# Draw to PDF
grDevices::pdf(pdf_path, width = 7, height = 6)
grid.draw(vp)
pushViewport(viewport())
grid.text("CHX vs DMSO Overlap", 0.5, 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()
grDevices::dev.off()

cat(sprintf("  Saved: %s\n", basename(png_path)))
cat(sprintf("  Saved: %s\n", basename(pdf_path)))

# Export gene lists
save_table(data.frame(gene = chx_genes, category = "CHX_enriched"), "chx_volcano_enriched")
save_table(data.frame(gene = dmso_genes, category = "DMSO_enriched"), "chx_volcano_dmso_enriched")

overlap_genes <- intersect(chx_genes, dmso_genes)
if (length(overlap_genes) > 0) {
  save_table(data.frame(gene = overlap_genes, category = "CHX_DMSO_overlap"), "chx_dmso_overlap")
} else {
  cat("  No overlap between CHX and DMSO enriched sets.\n")
}

cat("\nDone. Output saved to output/figures/.\n")
