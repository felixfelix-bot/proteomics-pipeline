###############################################################################
# 31_poster_figures.R
# Generate ALL figures for the LaTeX poster with STANDARDIZED FONT SIZES.
#
# This script:
#   1. Sources R/00_theme.R (defines theme_poster with uniform fonts)
#   2. Loads experiment data
#   3. For GO/KEGG: reads result tables from output/tables/ (already computed
#      by targeted-go, etc.) and re-creates dotplots with theme_poster()
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
cat(" Poster Figures — Standardized Fonts\n")
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
  ggsave(pdf_path, plot = plot, width = width, height = height,
         device = grDevices::cairo_pdf)
  ggsave(png_path, plot = plot, width = width, height = height,
         dpi = 300, bg = "white")
  cat(sprintf("  [POSTER] %s.pdf (%.0f x %.0f in)\n", name, width, height))
}

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
# 1. VOLCANO PLOT — TurboID TRIP4 vs WT
#=============================================================================
cat("[1/6] Volcano plot...\n")

df_volcano <- df_turbo
df_volcano$category <- "Not enriched"
df_volcano$category[!is.na(df_volcano$log2FC) &
  df_volcano$log2FC >= FC_THRESH &
  df_volcano$padj < 0.05] <- "Enriched"
df_volcano$category[df_volcano$gene %in% ASCC_CORE &
  df_volcano$log2FC >= FC_THRESH] <- "ASCC complex"
df_volcano$category[df_volcano$gene %in% known_interactors &
  df_volcano$log2FC >= FC_THRESH] <- "Known interactor"

# Factor with consistent ordering
df_volcano$category <- factor(df_volcano$category,
  levels = c("ASCC complex", "Known interactor", "Enriched", "Not enriched"))

# Colors matching pipeline
volcano_colors <- c(
  "ASCC complex" = "#0072B2",
  "Known interactor" = "#009E73",
  "Enriched" = "#D55E00",
  "Not enriched" = "grey80"
)

# Label ASCC + known interactors that are enriched
label_genes <- unique(c(ASCC_CORE,
  intersect(known_interactors,
    df_volcano$gene[df_volcano$category != "Not enriched"])))

p_volcano <- ggplot(df_volcano,
    aes(x = log2FC, y = neglog10p, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colors, name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TurboID TRIP4 vs Wild Type"
  ) +
  geom_text_repel(
    data = df_volcano[df_volcano$gene %in% label_genes &
      !is.na(df_volcano$log2FC), ],
    aes(label = gene), size = 2.5, fontface = "bold",
    max.overlaps = 20, show.legend = FALSE, color = "black"
  ) +
  theme_poster() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "grey80"))

save_poster(p_volcano, "volcano_turbo", width = 7, height = 6)

#=============================================================================
# 2. FLAG IP OVERLAP VOLCANO
#=============================================================================
cat("[2/6] Flag IP overlap volcano...\n")

cflag_exp <- find_experiment(experiments, "BK516_Cflag_vs_BK516_Ctrl")
nflag_exp <- find_experiment(experiments, "BK516_Nflag_vs_BK516_Ctrl")

cflag_sig <- if (!is.null(cflag_exp)) get_significant_genes(cflag_exp) else character(0)
nflag_sig <- if (!is.null(nflag_exp)) get_significant_genes(nflag_exp) else character(0)

df_flag <- df_turbo
df_flag$flag_cat <- "TRIP4 only"
turbo_sig <- df_flag$gene[df_flag$padj < 0.05 & df_flag$log2FC >= FC_THRESH]
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
  "C-Flag only" = "#0072B2",
  "N-Flag only" = "#CC79A7",
  "TRIP4 only" = "#D55E00",
  "Not enriched" = "grey80"
)

# Label validated proteins (both C+N)
label_flag <- intersect(intersect(cflag_sig, nflag_sig), turbo_sig)

p_flag <- ggplot(df_flag,
    aes(x = log2FC, y = neglog10p, color = flag_cat)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = flag_colors, name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "Flag IP Validation Overlay"
  ) +
  geom_text_repel(
    data = df_flag[df_flag$gene %in% label_flag & !is.na(df_flag$log2FC), ],
    aes(label = gene), size = 2.5, fontface = "bold",
    max.overlaps = 15, show.legend = FALSE, color = "black"
  ) +
  theme_poster() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "grey80"))

save_poster(p_flag, "volcano_flagip", width = 7, height = 6)

#=============================================================================
# 3. VENN DIAGRAM — TurboID vs C-Flag vs N-Flag
#=============================================================================
cat("[3/6] Venn diagram...\n")

# Use VennDiagram package if available, otherwise ggVennDiagram
venn_data <- list(
  "TurboID TRIP4" = turbo_sig,
  "C-Flag IP" = cflag_sig,
  "N-Flag IP" = nflag_sig
)

# Try ggVennDiagram (cleaner, ggplot-based, theme-compatible)
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
  # Fallback: static VennDiagram
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
    fontfamily = "sans",
    cat.fontfamily = "sans",
    main = "TurboID vs Flag IP Overlap",
    main.cex = 1.5,
    main.fontfamily = "sans"
  )
  NULL  # Return NULL — figure already saved by venn.diagram
})

if (!is.null(p_venn)) {
  save_poster(p_venn, "venn_overlap", width = 7, height = 6)
}

#=============================================================================
# 4-6. GO + KEGG DOTPLOTS — read from saved tables, re-plot with theme_poster
#=============================================================================
# This reads the CSV tables generated by targeted-go (R/11) and
# flagip-validated-go (R/30) and creates clean dotplots.
#
# If tables don't exist yet, prints a warning.

# Helper: create a standardized dotplot from a GO/KEGG result table
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

  # Parse GeneRatio (stored as "3/10" string → numeric)
  if ("GeneRatio" %in% names(df)) {
    df$ratio <- sapply(df$GeneRatio, function(x) {
      parts <- as.numeric(strsplit(as.character(x), "/")[[1]])
      if (length(parts) == 2 && parts[2] > 0) parts[1] / parts[2] else NA
    })
  } else {
    return(NULL)
  }

  # Sort by p.adjust (most significant at top)
  df <- df[order(df$p.adjust), ]
  df <- head(df, max_terms)
  df$Description <- factor(df$Description, levels = rev(df$Description))

  ggplot(df, aes(x = ratio, y = Description,
                 size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#D55E00", high = "#0072B2",
                         name = "p-adjusted") +
    scale_size_continuous(name = "Gene Count", range = c(2, 6)) +
    labs(
      x = "Gene Ratio",
      y = NULL,
      title = title_text
    ) +
    theme_poster() +
    theme(axis.text.y = element_text(size = POSTER_FONTS$axis_text),
          legend.position = "right")
}

# Find the most recent table file matching a pattern
find_latest_table <- function(pattern) {
  files <- list.files(TABLE_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[order(file.info(files)$mtime, decreasing = TRUE)][1]
}

# ---- GO BP dotplot (TurboID) ----
cat("[4/6] GO BP dotplot...\n")
go_bp_file <- find_latest_table("targeted_GO_turbo_trip4.*BP")
if (!is.null(go_bp_file)) {
  p_go <- make_poster_dotplot(go_bp_file,
    "GO Biological Process — TurboID TRIP4")
  if (!is.null(p_go)) save_poster(p_go, "go_bp_turbo", width = 8, height = 6)
} else {
  cat("  SKIP: Run 'make targeted-go' first\n")
}

# ---- KEGG dotplot (TurboID) ----
cat("[5/6] KEGG dotplot...\n")
kegg_file <- find_latest_table("targeted_KEGG_turbo_trip4")
if (!is.null(kegg_file)) {
  p_kegg <- make_poster_dotplot(kegg_file,
    "KEGG Pathways — TurboID TRIP4")
  if (!is.null(p_kegg)) save_poster(p_kegg, "kegg_turbo", width = 8, height = 6)
} else {
  cat("  SKIP: Run 'make targeted-go' first (KEGG runs automatically)\n")
}

# ---- Validated GO BP dotplot (C-Flag + N-Flag) ----
cat("[6/6] Validated GO BP dotplot...\n")
val_go_file <- find_latest_table("flagip_GO_validated_both.*BP")
if (!is.null(val_go_file)) {
  p_val <- make_poster_dotplot(val_go_file,
    "GO BP — Validated by Both C-Flag + N-Flag")
  if (!is.null(p_val)) save_poster(p_val, "go_bp_validated", width = 8, height = 6)
} else {
  cat("  SKIP: Run 'make flagip-validated-go' first\n")
}

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
