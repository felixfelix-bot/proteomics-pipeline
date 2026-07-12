###############################################################################
# diagnostics.R — Prints structural info about the data (NO actual values).
# Run: Rscript R/diagnostics.R
# Safe to share — outputs only counts, ranges, category sizes.
###############################################################################

source('R/01_config.R')
source('R/utils.R')

experiments <- load_all_experiments()
df <- experiments[["BK467_TRIP4_vs_BK467_WT"]]

df$neglog10p <- -log10(df$padj)

cat("=========================================\n")
cat(" DATA DIAGNOSTICS (safe to share)\n")
cat("=========================================\n\n")

cat(sprintf("Total proteins: %d\n", nrow(df)))
cat(sprintf("With valid gene names: %d\n", sum(!is.na(df$gene) & df$gene != "")))
cat(sprintf("With valid log2FC: %d\n", sum(!is.na(df$log2FC))))
cat(sprintf("With valid padj: %d\n", sum(!is.na(df$padj))))

cat("\n--- log2FC distribution ---\n")
cat(sprintf("  Range: %.2f to %.2f\n", min(df$log2FC, na.rm=TRUE), max(df$log2FC, na.rm=TRUE)))
cat(sprintf("  log2FC >= 1: %d proteins\n", sum(df$log2FC >= 1, na.rm=TRUE)))
cat(sprintf("  log2FC >= 0.5: %d proteins\n", sum(df$log2FC >= 0.5, na.rm=TRUE)))
cat(sprintf("  log2FC <= -1: %d proteins\n", sum(df$log2FC <= -1, na.rm=TRUE)))
cat(sprintf("  log2FC <= -0.5: %d proteins\n", sum(df$log2FC <= -0.5, na.rm=TRUE)))

cat("\n--- padj distribution ---\n")
cat(sprintf("  Range: %.2e to %.2e\n", min(df$padj, na.rm=TRUE), max(df$padj, na.rm=TRUE)))
cat(sprintf("  padj <= 0.05: %d proteins\n", sum(df$padj <= 0.05, na.rm=TRUE)))
cat(sprintf("  padj <= 0.1: %d proteins\n", sum(df$padj <= 0.1, na.rm=TRUE)))

cat("\n--- Combined thresholds ---\n")
cat(sprintf("  log2FC>=1 & padj<=0.1: %d proteins\n",
    sum(df$log2FC >= 1 & df$padj <= 0.1, na.rm=TRUE)))
cat(sprintf("  log2FC>=0.5 & padj<=0.05: %d proteins\n",
    sum(df$log2FC >= 0.5 & df$padj <= 0.05, na.rm=TRUE)))
cat(sprintf("  log2FC>=1 & padj<=0.05: %d proteins\n",
    sum(df$log2FC >= 1 & df$padj <= 0.05, na.rm=TRUE)))

cat("\n--- Lydia's seed thresholds ---\n")
cat(sprintf("  log2FC>7.5 & neglog10p>3: %d proteins\n",
    sum(df$log2FC > 7.5 & df$neglog10p > 3, na.rm=TRUE)))
cat(sprintf("  log2FC>2 & neglog10p>6: %d proteins\n",
    sum(df$log2FC > 2 & df$neglog10p > 6, na.rm=TRUE)))
cat(sprintf("  Seeds (union): %d proteins\n",
    sum((df$log2FC > 7.5 & df$neglog10p > 3) |
        (df$log2FC > 2 & df$neglog10p > 6), na.rm=TRUE)))

cat("\n--- Known interactors present ---\n")
interactors_file <- file.path(DATA_DIR, "known_interactors.txt")
known <- load_known_interactors(interactors_file)
present <- known[known %in% df$gene]
cat(sprintf("  Known interactors file: %d genes\n", length(known)))
cat(sprintf("  Present in this experiment: %d genes\n", length(present)))
cat(sprintf("  Of those, log2FC >= 1: %d\n",
    sum(present %in% df$gene[df$log2FC >= 1])))

ASCC_CORE <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
ascc_present <- ASCC_CORE[ASCC_CORE %in% df$gene]
cat(sprintf("\n  ASCC complex members present: %s\n",
    paste(ascc_present, collapse=", ")))

cat("\n=========================================\n")
cat(" END DIAGNOSTICS\n")
cat("=========================================\n")
