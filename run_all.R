###############################################################################
# run_all.R — Master pipeline: runs all analysis modules in sequence.
#
# Usage:
#   Rscript run_all.R          # run on real data in data/
#   Rscript run_all.R --test   # generate synthetic data first, then run
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
test_mode <- "--test" %in% args

project_root <- getwd()
Sys.setenv(PROJECT_ROOT = project_root)

# Set up logging (saves all output to output/logs/ with commit hash)
source(file.path(project_root, "R", "setup_logging.R"))
script_label <- ifelse(test_mode, "run_all_test", "run_all")
setup_logging(script_name = script_label)

source(file.path(project_root, "R", "01_config.R"))
source(file.path(project_root, "R", "utils.R"))

# Ensure logging is stopped even on error
.on_pipeline_error <- function() {
  try(stop_logging(success = FALSE), silent = TRUE)
}

cat("\n")
cat("=========================================\n")
cat(" TRIP4/ASCC Proteomics Analysis Pipeline\n")
cat("=========================================\n")
cat(sprintf(" Mode: %s\n", ifelse(test_mode, "TEST (synthetic data)", "REAL DATA")))
cat(sprintf(" Data directory: %s\n", DATA_DIR))
cat("\n")

if (test_mode) {
  cat("Generating synthetic test data...\n")
  source(file.path(project_root, "R", "generate_synthetic_data.R"))
  cat("\n")
}

# Core analysis
cat("Step 1/7: Volcano plots (Lydia-style)...\n")
source(file.path(project_root, "R", "02_volcano_plots.R"))

cat("\nStep 2/7: Venn diagrams + set extraction...\n")
source(file.path(project_root, "R", "03_venn_diagrams.R"))

cat("\nStep 3/7: GO enrichment...\n")
source(file.path(project_root, "R", "04_go_enrichment.R"))

# Advanced analysis (Lydia-style)
cat("\nStep 4/7: STRING network analysis...\n")
tryCatch(source(file.path(project_root, "R", "05_string_network.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

cat("\nStep 5/7: Gene family highlighting...\n")
tryCatch(source(file.path(project_root, "R", "06_gene_families.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

cat("\nStep 6/7: Cross-experiment overlap...\n")
tryCatch(source(file.path(project_root, "R", "07_overlap_analysis.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

cat("\nStep 7/7: Summary...\n")
cat("=========================================\n")
cat(" PIPELINE COMPLETE!\n")
cat("=========================================\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))

# List output files
fig_files <- list.files(FIGURE_DIR, pattern = "\\.(png|pdf)$")
tab_files <- list.files(TABLE_DIR, pattern = "\\.csv$")
cat(sprintf("\n  %d figures generated\n", length(fig_files)))
cat(sprintf("  %d tables generated\n", length(tab_files)))
cat("=========================================\n")

# Stop logging and print summary
stop_logging(success = TRUE)
