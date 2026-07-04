###############################################################################
# generate_synthetic_data.R
# Generates realistic synthetic human proteomics data matching the exact
# CSV format from the TRIP4/ASCC mass spectrometry study.
#
# Produces one CSV per experiment with real BK experiment names.
#
# Usage:
#   Rscript generate_synthetic_data.R
###############################################################################

# ---- Parameters ----
N_PROTEINS <- 3000
N_SIGNIFICANT_BASE <- 400

# Real human gene symbols for known interactors
KNOWN_INTERACTORS <- c(
  "ASCC1", "ASCC2", "ASCC3", "TRIP4", "MED1", "MED15", "MED23",
  "PABPC1", "ZNF598", "FMR1", "NEK1", "NEK9", "NEK6", "ANKRD28",
  "PKN1", "PKN3", "PPP6R1", "WNK1", "STK3", "UFSP2", "DDRGK1",
  "UFL1", "EP300", "NCOA1", "NCOA7", "NCOA2", "NCOA3"
)

# Gene families for highlighting
GPATCH_GENES <- c("AGGF1", "CMTR1", "GPATCH1", "GPATCH2", "GPATCH3",
                  "GPATCH4", "GPATCH8", "SUGP1", "SUGP2", "GPKOW")

set.seed(42)

# Pool of gene symbols: random + known + gene families + some DHX/DDX/LARP
RANDOM_GENES <- c(
  paste0("GENE", sprintf("%04d", 1:(N_PROTEINS - 100))),
  KNOWN_INTERACTORS,
  GPATCH_GENES,
  paste0("DHX", c(9, 15, 29, 30, 36, 38, 58)),
  paste0("DDX", c(3, 5, 17, 21, 41, 46, 50)),
  paste0("LARP", c(1, 4, 5, 7))
)
RANDOM_GENES <- unique(RANDOM_GENES)

# ---- Helper: generate a single experiment ----
generate_experiment <- function(n_sig, up_frac = 0.6, known_boost = TRUE) {
  n <- length(RANDOM_GENES)

  log2fc <- rnorm(n, mean = 0, sd = 0.4)
  padj <- runif(n, min = 0.06, max = 1.0)

  sig_indices <- sample(1:n, min(n_sig, n))
  n_up <- round(length(sig_indices) * up_frac)
  up_idx <- sample(sig_indices, n_up)
  down_idx <- setdiff(sig_indices, up_idx)

  log2fc[up_idx] <- runif(length(up_idx), min = 0.5, max = 4.5)
  padj[up_idx]   <- runif(length(up_idx), min = 0.00001, max = 0.049)
  log2fc[down_idx] <- runif(length(down_idx), min = -4.5, max = -0.5)
  padj[down_idx]   <- runif(length(down_idx), min = 0.00001, max = 0.049)

  if (known_boost) {
    known_idx <- which(RANDOM_GENES %in% KNOWN_INTERACTORS)
    known_sig <- sample(known_idx, max(1, round(length(known_idx) * 0.7)))
    log2fc[known_sig] <- rnorm(length(known_sig), mean = 2.5, sd = 1.0)
    padj[known_sig]   <- runif(length(known_sig), min = 0.00001, max = 0.04)
  }

  df <- data.frame(
    Gene = RANDOM_GENES,
    Entry.Name = paste0(RANDOM_GENES, "_HUMAN"),
    UniProt.ID = paste0("P", sprintf("%05d", sample(1:99999, n, replace = TRUE))),
    Description = paste0("Protein ", RANDOM_GENES),
    logFC = round(log2fc, 4),
    P.Val = round(pmin(padj * runif(n, 1, 5), 1), 6),
    adj.P.Val = round(padj, 6),
    log.P.Val = round(-log10(pmin(padj * runif(n, 1, 5), 1)), 4),
    log.adj.P.Val = round(-log10(padj), 4),
    n_data_points_BK516_Ctrl = sample(2:6, n, replace = TRUE),
    n_data_points_BK516_Cflag = sample(2:6, n, replace = TRUE),
    imputed_BK516_Cflag = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.1, 0.9)),
    imputed_BK516_Ctrl = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.1, 0.9)),
    is_significant = (padj < 0.05 & abs(log2fc) > 0.5)
  )

  df <- df[sample(nrow(df)), ]
  rownames(df) <- NULL
  return(df)
}

# ---- Generate all experiments with real names ----
cat("Generating synthetic TRIP4/ASCC proteomics data...\n\n")

experiments <- list(
  # TurboID (HeLa, BK467)
  BK467_TRIP4_vs_BK467_WT         = generate_experiment(672, 0.65),
  BK467_TRIP4_RA02_vs_BK467_TRIP4 = generate_experiment(102, 0.50),
  BK467_TRIP4_RA02_vs_BK467_WT    = generate_experiment(550, 0.60),
  # TurboID (HeLa, BK504)
  BK504_TRIP4_RA04_vs_BK504_TRIP4 = generate_experiment(150, 0.55),
  BK504_TRIP4_RA04_vs_BK467_WT    = generate_experiment(480, 0.62),
  # Flag IP (HEK293, BK516)
  BK516_Cflag_vs_BK516_Ctrl       = generate_experiment(320, 0.60),
  BK516_Nflag_vs_BK516_Ctrl       = generate_experiment(280, 0.58),
  BK516_Cflag_vs_BK516_Nflag      = generate_experiment(90, 0.52),
  # Flag IP + RA (HEK293, BK523)
  BK523_Cflag_RA04_vs_BK516_Cflag = generate_experiment(210, 0.55),
  BK523_Cflag_RA04_vs_BK523_Ctrl_RA04 = generate_experiment(340, 0.60)
)

# ---- Write to CSV ----
output_dir <- here::here("data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(experiments)) {
  filepath <- file.path(output_dir, paste0(name, "_diffEx_minProb.csv"))
  write.csv(experiments[[name]], filepath, row.names = FALSE)
  n_sig <- sum(experiments[[name]]$adj.P.Val < 0.05 & abs(experiments[[name]]$logFC) > 0.5)
  cat(sprintf("  %s.csv: %d proteins (%d significant)\n", name,
              nrow(experiments[[name]]), n_sig))
}

# ---- Write known interactors ----
interactors_path <- file.path(output_dir, "known_interactors.txt")
writeLines(KNOWN_INTERACTORS, interactors_path)
cat(sprintf("\n  known_interactors.txt: %d genes\n", length(KNOWN_INTERACTORS)))

cat("\nDone. All synthetic data written to data/\n")
