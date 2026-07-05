###############################################################################
# print_context.R
# Prints debugging context to the terminal — column names, file structure,
# overlapping gene sets — without modifying any data.
#
# Use this to understand what's in the CSV files and why certain gene
# overlaps might be empty. The output can be reviewed for sensitive data
# before sharing with the AI.
#
# Usage:
#   make context
###############################################################################

# Source config and utils (same setup as run_all.R)
source("R/01_config.R")
source("R/utils.R")

cat("\n=========================================\n")
cat(" PROTEOMICS PIPELINE — DATA CONTEXT\n")
cat("=========================================\n\n")

# ---- Load config and data ----
cat("[1] Loading experiment data...\n")
experiments <- load_all_experiments()

# ---- Show file info ----
for (name in names(experiments)) {
  df <- experiments[[name]]
  cat(sprintf("\n  --- %s ---\n", name))
  cat(sprintf("  Rows: %d\n", nrow(df)))
  cat(sprintf("  Columns: %s\n", paste(colnames(df), collapse = ", ")))

  n_sig <- sum(df$padj < P_VALUE_CUTOFF & abs(df$log2FC) > LOG2FC_CUTOFF)
  n_up <- sum(df$padj < P_VALUE_CUTOFF & df$log2FC > LOG2FC_CUTOFF)
  n_dn <- sum(df$padj < P_VALUE_CUTOFF & df$log2FC < -LOG2FC_CUTOFF)
  cat(sprintf("  Significant: %d (%d up, %d down)\n", n_sig, n_up, n_dn))

  # Show first few and last few gene names (NOT the actual values — just
  # the structure) to illustrate naming conventions
  cat(sprintf("  Gene names (first 5): %s\n",
              paste(head(df$gene, 5), collapse = ", ")))
}

# ---- Overlap analysis for flag IP volcano ----
cat("\n\n[2] Flag IP overlap analysis...\n")
turbo_exp <- "BK467_TRIP4_vs_BK467_WT"
cflag_exp <- "BK516_Cflag_vs_BK516_Ctrl"
nflag_exp <- "BK516_Nflag_vs_BK516_Ctrl"

for (exp_name in c(turbo_exp, cflag_exp, nflag_exp)) {
  if (exp_name %in% names(experiments)) {
    cat(sprintf("  ✓ %s: FOUND\n", exp_name))
  } else {
    cat(sprintf("  ✗ %s: NOT FOUND\n", exp_name))
    cat(sprintf("    Available names with similar pattern:\n"))
    for (avail in names(experiments)) {
      if (grepl("BK467|BK516|TRIP4|Cflag|Nflag|Ctrl", avail)) {
        cat(sprintf("      - %s\n", avail))
      }
    }
  }
}

if (all(c(turbo_exp, cflag_exp, nflag_exp) %in% names(experiments))) {
  turbo_sig <- experiments[[turbo_exp]]$gene[
    experiments[[turbo_exp]]$padj < P_VALUE_CUTOFF &
    experiments[[turbo_exp]]$log2FC > LOG2FC_CUTOFF
  ]
  cflag_sig <- experiments[[cflag_exp]]$gene[
    experiments[[cflag_exp]]$padj < P_VALUE_CUTOFF
  ]
  nflag_sig <- experiments[[nflag_exp]]$gene[
    experiments[[nflag_exp]]$padj < P_VALUE_CUTOFF
  ]

  cat(sprintf("\n  TRIP4-enriched (positive log2FC, sig): %d\n", length(turbo_sig)))
  cat(sprintf("  C-flag significant (any direction):    %d\n", length(cflag_sig)))
  cat(sprintf("  N-flag significant (any direction):    %d\n", length(nflag_sig)))

  in_cflag_only <- intersect(turbo_sig, setdiff(cflag_sig, nflag_sig))
  in_nflag_only <- intersect(turbo_sig, setdiff(nflag_sig, cflag_sig))
  in_both <- intersect(turbo_sig, intersect(cflag_sig, nflag_sig))

  cat(sprintf("\n  TRIP4-enriched AND C-flag only: %d genes\n",
              length(in_cflag_only)))
  cat(sprintf("  TRIP4-enriched AND N-flag only: %d genes\n",
              length(in_nflag_only)))
  cat(sprintf("  TRIP4-enriched AND both C+N:    %d genes\n",
              length(in_both)))

  if (length(in_cflag_only) > 0 && length(in_cflag_only) <= 20) {
    cat(sprintf("  Gene names (C-flag only): %s\n",
                paste(in_cflag_only, collapse = ", ")))
  } else if (length(in_cflag_only) > 20) {
    cat(sprintf("  Gene names (C-flag only, first 20): %s\n",
                paste(head(in_cflag_only, 20), collapse = ", ")))
  }

  if (length(in_nflag_only) > 0 && length(in_nflag_only) <= 20) {
    cat(sprintf("  Gene names (N-flag only): %s\n",
                paste(in_nflag_only, collapse = ", ")))
  } else if (length(in_nflag_only) > 20) {
    cat(sprintf("  Gene names (N-flag only, first 20): %s\n",
                paste(head(in_nflag_only, 20), collapse = ", ")))
  }

  if (length(in_both) > 0) {
    cat(sprintf("  Gene names (both C+N): %s\n",
                paste(in_both, collapse = ", ")))
  }

  # Check: are the genes that we expect to find actually in the data?
  ascc_genes <- c("TRIP4", "ASCC1", "ASCC2", "ASCC3")
  cat(sprintf("\n  ASCC complex genes in TurboID data:\n"))
  for (g in ascc_genes) {
    present <- g %in% experiments[[turbo_exp]]$gene
    sig_in_turbo <- g %in% turbo_sig
    sig_in_cflag <- if (g %in% experiments[[cflag_exp]]$gene) {
      experiments[[cflag_exp]]$padj[experiments[[cflag_exp]]$gene == g] < P_VALUE_CUTOFF
    } else NA
    sig_in_nflag <- if (g %in% experiments[[nflag_exp]]$gene) {
      experiments[[nflag_exp]]$padj[experiments[[nflag_exp]]$gene == g] < P_VALUE_CUTOFF
    } else NA
    cat(sprintf("    %s: in TurboID=%s, sig TurboID=%s, sig Cflag=%s, sig Nflag=%s\n",
                g, present, sig_in_turbo,
                if (is.na(sig_in_cflag)) "NA" else sig_in_cflag,
                if (is.na(sig_in_nflag)) "NA" else sig_in_nflag))
  }
}

cat("\n=========================================\n")
cat(" Context dump complete.\n")
cat("=========================================\n")
cat("\nShare this output with the AI to help diagnose overlap issues.\n")
cat("If any gene names above are sensitive, redact them before sharing.\n")
