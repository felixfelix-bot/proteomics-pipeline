###############################################################################
# 09_flagip_volcano.R
# TRIP4 TurboID vs WT volcano — Flag IP overlap highlighting.
#
# WHAT THIS PLOT SHOWS:
#   Start with the TRIP4 TurboID vs WT volcano.
#   Only genes significantly enriched IN TRIP4 (positive log2FC) are candidates.
#   For each TRIP4-enriched gene, check: was it also found as significant
#   in the C-Flag and/or N-Flag experiments?
#
# CATEGORIES (priority order — later overwrites earlier):
#   "in_both"      = Sig in TurboID TRIP4 AND C-Flag AND N-Flag    → Blue
#   "in_cflag"     = Sig in TurboID TRIP4 AND C-Flag (not N-Flag)  → Green
#   "in_nflag"     = Sig in TurboID TRIP4 AND N-Flag (not C-Flag)  → Vermillion
#   "trip4_only"   = Sig in TurboID TRIP4 only (no Flag IP support) → Pale orange
#   "wt_enriched"  = Sig but enriched in WT (negative log2FC)      → Gray
#   "nonsig"       = Not significant in TurboID                    → Light gray
#
# Only the Flag IP-validated categories (in_both, in_cflag, in_nflag)
# get text labels — trip4_only and wt_enriched do not.
#
# Usage:
#   make flagip-volcano
#
# FIXED: 2026-07-05 — was incorrectly doing C-flag∩N-flag only without
#   requiring TurboID TRIP4 significance. Now starts from TurboID TRIP4-
#   enriched genes and cross-references against Flag IP experiments.
###############################################################################
cat("\n=========================================\n")
cat(" Flag IP Overlap Volcano (FIXED)\n")
cat("=========================================\n\n")

source("R/00_theme.R")

experiments <- load_all_experiments()

# ---- DIAGNOSTIC: Print ALL experiment names found in data ----
# This helps debug when experiment names don't match expected values.
# The names come from the CSV filenames (minus _diffEx_minProb.csv suffix).
cat("\n--- EXPERIMENT DIAGNOSTICS ---\n")
cat(sprintf("Total experiments loaded: %d\n", length(experiments)))
cat("All experiment names:\n")
for (nm in names(experiments)) {
  cat(sprintf("  - %s (%d proteins)\n", nm, nrow(experiments[[nm]])))
}
cat("--- END DIAGNOSTICS ---\n\n")

# ---- Step 1: Find genes significant in each Flag IP experiment ----
cat("Loading Flag IP experiment data...\n")

flag_c <- "BK516_Cflag_vs_BK516_Ctrl"
flag_n <- "BK516_Nflag_vs_BK516_Ctrl"

cflag_sig <- character(0)
nflag_sig <- character(0)

# Use find_experiment() to handle dedup suffixes (BK516_Cflag_vs_BK516_Ctrl__1)
cflag_exp <- find_experiment(experiments, flag_c)
if (!is.null(cflag_exp)) {
  cflag_sig <- get_significant_genes(cflag_exp)
  cat(sprintf("  C-flag significant:   %d proteins\n", length(cflag_sig)))
  # Print sample gene names for debugging
  if (length(cflag_sig) > 0) {
    cat(sprintf("    Sample genes (first 10): %s\n",
                paste(head(cflag_sig, 10), collapse = ", ")))
  }
} else {
  cat("  WARNING: C-flag experiment not found:", flag_c, "\n")
  cat("  Tried also:", paste0(flag_c, "__1"), ",", paste0(flag_c, "__2"), "\n")
  cat("  Check the experiment names listed above.\n")
}

nflag_exp <- find_experiment(experiments, flag_n)
if (!is.null(nflag_exp)) {
  nflag_sig <- get_significant_genes(nflag_exp)
  cat(sprintf("  N-flag significant:   %d proteins\n", length(nflag_sig)))
  if (length(nflag_sig) > 0) {
    cat(sprintf("    Sample genes (first 10): %s\n",
                paste(head(nflag_sig, 10), collapse = ", ")))
  }
} else {
  cat("  WARNING: N-flag experiment not found:", flag_n, "\n")
  cat("  Tried also:", paste0(flag_n, "__1"), ",", paste0(flag_n, "__2"), "\n")
  cat("  Check the experiment names listed above.\n")
}

# ---- Step 2: Load TurboID TRIP4 vs WT data ----
cat("\nLoading TurboID TRIP4 vs WT data...\n")
turbo_main <- "BK467_TRIP4_vs_BK467_WT"

if (!turbo_main %in% names(experiments)) {
  cat("ERROR: Main experiment not found:", turbo_main, "\n")
  quit(status = 1)
}

df <- experiments[[turbo_main]]
df$neglog10p <- -log10(df$padj)
cat(sprintf("  TurboID dataset: %d proteins\n", nrow(df)))

# ---- Step 3: Classify each protein ----
# Significance thresholds from config (padj < 0.05, |log2FC| > 0.5)
# But for TRIP4 enrichment, we only care about POSITIVE log2FC
cat("\nClassifying proteins...\n")

# Start: all are "nonsig"
df$category <- "nonsig"

# WT-enriched (significant, negative log2FC) → gray, no labels
wt_idx <- df$padj < P_VALUE_CUTOFF & df$log2FC < -LOG2FC_CUTOFF
df$category[wt_idx] <- "wt_enriched"

# TRIP4-enriched (significant, positive log2FC) — these are the candidates
trip4_idx <- df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF
df$category[trip4_idx] <- "trip4_only"

# Now check TRIP4-enriched genes against Flag IP data
# in_nflag: in N-flag but NOT in C-flag
in_nflag <- trip4_idx & df$gene %in% nflag_sig & !(df$gene %in% cflag_sig)
df$category[in_nflag] <- "in_nflag"

# in_cflag: in C-flag but NOT in N-flag
in_cflag <- trip4_idx & df$gene %in% cflag_sig & !(df$gene %in% nflag_sig)
df$category[in_cflag] <- "in_cflag"

# in_both: in BOTH C-flag and N-flag
in_both <- trip4_idx & df$gene %in% cflag_sig & df$gene %in% nflag_sig
df$category[in_both] <- "in_both"

# Report counts
cat(sprintf("  TRIP4-enriched (total):     %d proteins\n", sum(trip4_idx)))
cat(sprintf("    -> Also in C-flag only:    %d\n", sum(in_cflag)))
cat(sprintf("    -> Also in N-flag only:    %d\n", sum(in_nflag)))
cat(sprintf("    -> Also in BOTH:           %d\n", sum(in_both)))
cat(sprintf("    -> TRIP4 only (no Flag IP): %d\n",
            sum(trip4_idx) - sum(in_cflag | in_nflag | in_both)))
cat(sprintf("  WT-enriched:                 %d proteins\n", sum(wt_idx)))
cat(sprintf("  Not significant:             %d proteins\n",
            nrow(df) - sum(trip4_idx) - sum(wt_idx)))

# ---- DIAGNOSTIC: If overlap is zero, explain why ----
if (sum(in_both) == 0 && sum(in_cflag) == 0 && sum(in_nflag) == 0) {
  cat("\n  ========================================\n")
  cat("  OVERLAP DIAGNOSTIC: No Flag IP overlap found!\n")
  cat("  ========================================\n")
  cat(sprintf("  TRIP4-enriched genes: %d\n", sum(trip4_idx)))
  cat(sprintf("  C-flag significant:   %d\n", length(cflag_sig)))
  cat(sprintf("  N-flag significant:   %d\n", length(nflag_sig)))

  # Check raw intersection (without TRIP4 filter)
  trip4_genes <- df$gene[trip4_idx]
  raw_cflag_overlap <- intersect(trip4_genes, cflag_sig)
  raw_nflag_overlap <- intersect(trip4_genes, nflag_sig)
  cat(sprintf("\n  Raw overlap TRIP4 x C-flag: %d genes\n", length(raw_cflag_overlap)))
  cat(sprintf("  Raw overlap TRIP4 x N-flag: %d genes\n", length(raw_nflag_overlap)))

  if (length(raw_cflag_overlap) > 0) {
    cat(sprintf("    Genes: %s\n", paste(head(raw_cflag_overlap, 20), collapse = ", ")))
  }
  if (length(raw_nflag_overlap) > 0) {
    cat(sprintf("    Genes: %s\n", paste(head(raw_nflag_overlap, 20), collapse = ", ")))
  }

  # Check if gene names use different cases
  cat("\n  Checking gene name format...\n")
  cat(sprintf("    TurboID sample: %s\n",
              paste(head(df$gene[trip4_idx], 5), collapse = ", ")))
  if (length(cflag_sig) > 0) {
    cat(sprintf("    C-flag sample:  %s\n",
                paste(head(cflag_sig, 5), collapse = ", ")))
  }
  if (length(nflag_sig) > 0) {
    cat(sprintf("    N-flag sample:  %s\n",
                paste(head(nflag_sig, 5), collapse = ", ")))
  }
  cat("  ========================================\n\n")
}

# Factor with explicit ordering (priority: later levels drawn on top)
df$category <- factor(df$category,
  levels = c("in_both", "in_cflag", "in_nflag",
             "trip4_only", "wt_enriched", "nonsig"))

# ---- Step 4: Only label Flag IP-validated hits ----
label_data <- df[df$category %in% c("in_both", "in_cflag", "in_nflag"), ]

# Also save the Flag IP validated gene lists
validate_dir <- file.path(TABLE_DIR, "flagip_validation")
dir.create(validate_dir, showWarnings = FALSE, recursive = TRUE)

if (sum(in_both) > 0) {
  vp <- safe_filepath(validate_dir, "trip4_validated_by_both_flagIP", ".csv")
  write.csv(data.frame(gene = df$gene[in_both]), vp, row.names = FALSE)
}
if (sum(in_cflag) > 0) {
  vp <- safe_filepath(validate_dir, "trip4_validated_by_cflag_only", ".csv")
  write.csv(data.frame(gene = df$gene[in_cflag]), vp, row.names = FALSE)
}
if (sum(in_nflag) > 0) {
  vp <- safe_filepath(validate_dir, "trip4_validated_by_nflag_only", ".csv")
  write.csv(data.frame(gene = df$gene[in_nflag]), vp, row.names = FALSE)
}

# ---- Step 5: Build the volcano plot ----
# Colors from GLOBAL_COLORS for consistency across all plots
FLAGIP_COLORS <- c(
  "in_both"    = GLOBAL_COLORS[["flag_both"]],     # Yellow — both Flag IP methods
  "in_cflag"   = GLOBAL_COLORS[["flag_c_only"]],   # Sky blue — C-Flag only
  "in_nflag"   = GLOBAL_COLORS[["flag_n_only"]],   # Pink — N-Flag only
  "trip4_only" = GLOBAL_COLORS[["enriched_up"]],   # Vermillion — TRIP4-enriched only
  "wt_enriched" = GLOBAL_COLORS[["enriched_dn"]],  # Gray
  "nonsig"      = GLOBAL_COLORS[["nonsig"]]        # Light gray
)

FLAGIP_LABELS <- c(
  "in_both"    = "Validated by C-Flag + N-Flag",
  "in_cflag"   = "Validated by C-Flag only",
  "in_nflag"   = "Validated by N-Flag only",
  "trip4_only" = "TRIP4-enriched only",
  "wt_enriched" = "Enriched in WT",
  "nonsig"      = "Not significant"
)

cat("\nBuilding volcano plot...\n")

p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = category)) +
  ggplot2::geom_point(alpha = 0.5, size = 0.8) +
  ggplot2::scale_color_manual(
    values = FLAGIP_COLORS,
    labels = FLAGIP_LABELS,
    name = NULL, drop = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1)  # Legend dot size
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    ggplot2::aes(label = gene),
    size = 2.5, fontface = "bold",
    max.overlaps = 30, show.legend = FALSE,
    bg.color = "white", bg.r = 0.15
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(P_VALUE_CUTOFF),
    linetype = "dashed", color = "grey50", linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
    linetype = "dashed", color = "grey50", linewidth = 0.3
  ) +
  ggplot2::labs(
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~(adjusted~italic(p)~value)),
    title = "TRIP4 TurboID vs WT — Flag IP Validation"
  ) +
  theme_poster(font_size = 14) +
  ggplot2::theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

save_figure(p, "flagip_overlap_volcano_BK467_TRIP4_vs_WT",
            width = 16, height = 12)

cat("\n=========================================\n")
cat(" Flag IP overlap volcano complete!\n")
cat("=========================================\n")
