###############################################################################
# generate_synthetic_data.R
# Generates realistic synthetic human proteomics data for testing the pipeline.
# Produces CSV files matching the expected format: gene, log2FC, padj
#
# Also generates a known interactors list for volcano plot labeling.
#
# Usage:
#   Rscript generate_synthetic_data.R
############################################################################---

# ---- Parameters ----
N_PROTEINS <- 3000          # Total proteins detected per experiment
N_SIGNIFICANT <- 400        # Proteins passing significance threshold
N_KNOWN_INTERACTORS <- 25   # Known interaction partners to highlight

# Real human gene symbols for known interactors (well-studied proteins)
KNOWN_INTERACTORS <- c(
  "TP53", "MYC", "BRCA1", "BRCA2", "EGFR", "AKT1", "MTOR", "PIK3CA",
  "RB1", "CCND1", "CDK4", "CDK6", "MDM2", "CDKN2A", "ATM", "ATR",
  "CHEK1", "CHEK2", "XRCC5", "XRCC6", "HSP90AA1", "HSPA8", "UBC",
  "SUMO1", "RPS27A"
)

# Random human gene symbols (broad sampling across pathways)
RANDOM_GENES <- c(
  paste0("GENE", sprintf("%04d", 1:(N_PROTEINS - length(KNOWN_INTERACTORS)))),
  KNOWN_INTERACTORS
)

set.seed(42)  # Reproducible synthetic data

# ---- Helper: generate a single experiment's proteomics results ----
generate_experiment <- function(name, n_sig, up_regulated_frac = 0.5) {
  n <- length(RANDOM_GENES)

  # Most proteins: not significant (low FC, high p-value)
  log2fc <- rnorm(n, mean = 0, sd = 0.4)
  padj <- runif(n, min = 0.06, max = 1.0)

  # Make some proteins significant
  sig_indices <- sample(1:n, n_sig)
  n_up <- round(n_sig * up_regulated_frac)
  up_idx <- sample(sig_indices, n_up)
  down_idx <- setdiff(sig_indices, up_idx)

  # Up-regulated: positive log2FC, low p-value
  log2fc[up_idx] <- runif(length(up_idx), min = 1.0, max = 4.5)
  padj[up_idx]   <- runif(length(up_idx), min = 0.00001, max = 0.049)

  # Down-regulated: negative log2FC, low p-value
  log2fc[down_idx] <- runif(length(down_idx), min = -4.5, max = -1.0)
  padj[down_idx]   <- runif(length(down_idx), min = 0.00001, max = 0.049)

  # Ensure known interactors are significant in at least some experiments
  known_idx <- which(RANDOM_GENES %in% KNOWN_INTERACTORS)
  known_sig <- sample(known_idx, max(1, round(length(known_idx) * 0.6)))
  log2fc[known_sig] <- rnorm(length(known_sig), mean = 2.5, sd = 1.0)
  padj[known_sig]   <- runif(length(known_sig), min = 0.00001, max = 0.04)

  df <- data.frame(
    Gene = RANDOM_GENES,
    Entry.Name = paste0(RANDOM_GENES, "_HUMAN"),
    UniProt.ID = paste0("P", sprintf("%05d", sample(1:99999, n, replace = TRUE))),
    Description = paste0("Uncharacterized protein ", RANDOM_GENES),
    logFC = round(log2fc, 4),
    P.Val = round(pmin(padj * runif(n, 1, 5), 1), 6),
    adj.P.Val = round(padj, 6),
    log.P.Val = round(-log10(pmin(padj * runif(n, 1, 5), 1)), 4),
    log.adj.P.Val = round(-log10(padj), 4),
    n_data_points_BK516_Ctrl = sample(2:6, n, replace = TRUE),
    n_data_points_BK516_Cflag = sample(2:6, n, replace = TRUE),
    imputed_BK516_Cflag = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.1, 0.9)),
    imputed_BK516_Ctrl = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.1, 0.9)),
    is_significant = (padj < 0.05 & abs(log2fc) > 1)
  )

  # Shuffle rows
  df <- df[sample(nrow(df)), ]
  rownames(df) <- NULL

  return(df)
}

# ---- Generate all experiments ----
cat("Generating synthetic proteomics data...\n\n")

experiments <- list(
  wt_vs_poi          = generate_experiment("WT_vs_POI",        350, 0.55),
  poi_vs_poi_hormone  = generate_experiment("POI_vs_POIHormone", 300, 0.60),
  turboid            = generate_experiment("TurboID",           200, 0.65),
  flag_ip            = generate_experiment("Flag_IP",           180, 0.58),
  crac_rna           = generate_experiment("CRAC_RNA",          220, 0.50)
)

# ---- Write to CSV ----
output_dir <- here::here("data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(experiments)) {
  filepath <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(experiments[[name]], filepath, row.names = FALSE)
  n_sig <- sum(experiments[[name]]$adj.P.Val < 0.05 & abs(experiments[[name]]$logFC) > 1)
  cat(sprintf("  %s.csv: %d proteins (%d significant)\n", name, nrow(experiments[[name]]), n_sig))
}

# ---- Write known interactors list ----
interactors_path <- file.path(output_dir, "known_interactors.txt")
writeLines(KNOWN_INTERACTORS, interactors_path)
cat(sprintf("  known_interactors.txt: %d genes\n", length(KNOWN_INTERACTORS)))

cat("\nDone. All synthetic data written to data/\n")
cat(sprintf("To test the full pipeline: source('R/01_config.R') then run analysis scripts\n"))
