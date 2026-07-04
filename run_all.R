###############################################################################
# run_all.R — Master pipeline: runs all analysis modules in sequence.
#
# Usage:
#   Rscript run_all.R          # run on real data in data/
#   Rscript run_all.R --test   # generate synthetic data first, then run
###############################################################################

# Parse args
args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

# Set project root (where this script lives)
project_root <- getwd()
Sys.setenv(PROJECT_ROOT = project_root)

# Source config + utils
source(file.path(project_root, "R", "01_config.R"))
source(file.path(project_root, "R", "utils.R"))

cat("\n")
cat("=========================================\n")
cat(" Proteomics Analysis Pipeline\n")
cat("=========================================\n")
cat(sprintf(" Mode: %s\n", ifelse(test_mode, "TEST (synthetic data)", "REAL DATA")))
cat(sprintf(" Data directory: %s\n", DATA_DIR))
cat("\n")

if (test_mode) {
  cat("Generating synthetic test data...\n")
  source(file.path(project_root, "R", "generate_synthetic_data.R"))
  cat("\n")
}

# Run each module
cat("Step 1/3: Volcano plots...\n")
source(file.path(project_root, "R", "02_volcano_plots.R"))

cat("\nStep 2/3: Venn diagrams...\n")
source(file.path(project_root, "R", "03_venn_diagrams.R"))

cat("\nStep 3/3: GO enrichment...\n")
source(file.path(project_root, "R", "04_go_enrichment.R"))

cat("\n=========================================\n")
cat(" PIPELINE COMPLETE!\n")
cat("=========================================\n")
cat(sprintf(" All figures in: %s/\n", FIGURE_DIR))
cat(sprintf(" All tables in:  %s/\n", TABLE_DIR))
cat("=========================================\n")
