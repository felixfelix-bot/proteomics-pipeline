###############################################################################
# run_all.R — Master pipeline: runs all analysis modules in sequence.
#
# This is the MAIN entry point. When you type "make all" or
# "Rscript run_all.R", this file runs. It loads the configuration,
# sets up logging, and then runs each analysis step one by one.
#
# Usage:
#   Rscript run_all.R          # run on real data in data/
#   Rscript run_all.R --test   # generate synthetic (fake) data first, then run
###############################################################################

# ---- Get command-line arguments ----
# commandArgs() reads whatever you typed after "Rscript run_all.R" on the
# command line. trailingOnly = TRUE means "only get the stuff after the
# script name", ignoring R's own internal arguments.
# For example: "Rscript run_all.R --test" gives args = c("--test")
args <- commandArgs(trailingOnly = TRUE)

# The %in% operator checks if a value exists inside a vector.
# Here we check if "--test" was passed. Result: TRUE or FALSE.
# If TRUE, we'll generate fake data for testing instead of using real data.
test_mode <- "--test" %in% args

# ---- Determine the project root directory ----
# getwd() returns the "current working directory" — the folder R is
# running from. We store this so all other scripts know where to find
# the data/, R/, and output/ folders.
project_root <- getwd()

# Sys.setenv() sets an "environment variable" that other programs can read.
# We store the project root here so any subprocess can find it.
Sys.setenv(PROJECT_ROOT = project_root)

# ---- Set up logging ----
# source() is how R loads and runs another R file. It's like #include in C
# or import in Python — it reads the file and executes all the code in it.
# Any functions defined in that file become available to use.
#
# Here we load setup_logging.R, which defines the functions setup_logging()
# and stop_logging(). After sourcing, we can call them below.
source(file.path(project_root, "R", "setup_logging.R"))

# file.path() joins folder names into a path that works on any OS.
# On Windows: "C:\...\R\setup_logging.R"
# On Linux:   "/home/.../R/setup_logging.R"

# ifelse() is R's ternary operator: ifelse(condition, value_if_true, value_if_false)
# We pick a different label for the log filename depending on mode.
script_label <- ifelse(test_mode, "run_all_test", "run_all")

# Start logging. From this point on, everything printed to the console
# is ALSO saved to a log file in output/logs/. The log filename includes
# the timestamp and git commit hash for traceability.
setup_logging(script_name = script_label)

# ---- Load configuration and shared utilities ----
# 01_config.R defines all thresholds (p-value cutoff, fold change cutoff),
# file paths, experiment names, color schemes, and two key functions:
#   - load_proteomics_csv(): reads a CSV and validates its columns
#   - get_significant_genes(): filters to only significant proteins
source(file.path(project_root, "R", "01_config.R"))

# utils.R defines shared helper functions used by ALL analysis scripts:
#   - save_figure(): saves a plot as both PNG and PDF
#   - save_table(): saves a data frame as CSV
#   - discover_diffex_csvs(): finds all data files recursively
#   - load_all_experiments(): loads all CSVs into a named list
source(file.path(project_root, "R", "utils.R"))

# ---- Error safety net ----
# Define a cleanup function that will run if the pipeline crashes.
# stop_logging() closes the log file properly and marks the status as FAILED.
# try(..., silent = TRUE) means "try this, and if it errors, don't crash —
# just silently ignore the error." We don't want the cleanup itself to crash.
.on_pipeline_error <- function() {
  try(stop_logging(success = FALSE), silent = TRUE)
}

# ---- Print a header banner ----
# cat() is R's version of print/console.log. It prints text to the screen.
# \n means "new line". We use cat() for plain text output.
# sprintf() is like printf in C or f-strings in Python — it inserts values
# into a string. %s = string, %d = integer, %.2f = 2 decimal places.
cat("\n")
cat("=========================================\n")
cat(" TRIP4/ASCC Proteomics Analysis Pipeline\n")
cat("=========================================\n")
cat(sprintf(" Mode: %s\n", ifelse(test_mode, "TEST (synthetic data)", "REAL DATA")))
cat(sprintf(" Data directory: %s\n", DATA_DIR))
cat("\n")

# ---- Generate synthetic test data if in test mode ----
# In test mode, we generate fake proteomics data so you can see what
# the pipeline does without needing real data. The fake data has the
# exact same format as real mass spec output.
if (test_mode) {
  cat("Generating synthetic test data...\n")
  source(file.path(project_root, "R", "generate_synthetic_data.R"))
  cat("\n")
}

# =====================================================================
# STEP 1: Volcano Plots
# =====================================================================
# A volcano plot shows all proteins as dots, with:
#   X-axis: log2 fold change (how much the protein changed)
#   Y-axis: -log10 p-value (how confident we are it really changed)
# Proteins in the top-right and top-left corners are the most interesting:
# they changed a lot AND we're very confident about it.
#
# Lydia's style adds color-coded categories on top: known interactors
# (orange), Flag IP hits (blue/purple), gene families (green/brown), etc.
cat("Step 1/7: Volcano plots (Lydia-style)...\n")
source(file.path(project_root, "R", "02_volcano_plots.R"))

# =====================================================================
# STEP 2: Venn Diagrams
# =====================================================================
# Venn diagrams show the OVERLAP between different experiments.
# For example: "Which proteins are significant in BOTH TurboID AND Flag IP?"
# The overlapping region contains proteins found by both methods — these
# are high-confidence interactors.
cat("\nStep 2/7: Venn diagrams + set extraction...\n")
source(file.path(project_root, "R", "03_venn_diagrams.R"))

# =====================================================================
# STEP 3: GO Enrichment Analysis
# =====================================================================
# GO (Gene Ontology) enrichment asks: "Among the significant proteins,
# are certain biological pathways over-represented?"
# For example: if 30 of our significant proteins are involved in
# transcription, but only 5% of all human proteins are, that's enrichment.
# This tells us WHAT the protein of interest (TRIP4) might be doing.
cat("\nStep 3/7: GO enrichment...\n")
source(file.path(project_root, "R", "04_go_enrichment.R"))

# =====================================================================
# STEP 4: STRING Network Analysis
# =====================================================================
# STRING is a database of known protein-protein interactions.
# We map our significant proteins onto the STRING network to see:
#   - Do our hits cluster together (suggesting a real biological complex)?
#   - What other proteins are connected to our hits (candidate interactors)?
#
# tryCatch() is R's try/catch error handler. It runs the first expression,
# and if it errors, runs the second expression instead of crashing.
# Steps 4-6 are wrapped in tryCatch because they depend on internet
# access (STRING downloads) and external databases that might fail.
# If one step fails, the pipeline continues with the next step.
cat("\nStep 4/7: STRING network analysis...\n")
tryCatch(source(file.path(project_root, "R", "05_string_network.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

# =====================================================================
# STEP 5: Gene Family Highlighting
# =====================================================================
# Highlights members of specific gene families (GPATCH, DHX, DDX, LARP)
# on volcano plots. These families are RNA-binding proteins involved in
# RNA processing — relevant to TRIP4's suspected role in transcription.
cat("\nStep 5/7: Gene family highlighting...\n")
tryCatch(source(file.path(project_root, "R", "06_gene_families.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

# =====================================================================
# STEP 6: Cross-Experiment Overlap
# =====================================================================
# Compares results across experiments:
#   - Overlays Flag IP hits on the TurboID volcano plot
#   - Identifies RA (retinoic acid) specific changes
#   - Checks for CRAC (RNA interactome) data overlap
# This step answers: "Which proteins show up in MULTIPLE methods?"
cat("\nStep 6/7: Cross-experiment overlap...\n")
tryCatch(source(file.path(project_root, "R", "07_overlap_analysis.R")),
         error = function(e) cat("  SKIPPED:", conditionMessage(e), "\n"))

# =====================================================================
# STEP 7: Summary
# =====================================================================
# Count how many files were generated and print a summary.
cat("\nStep 7/7: Summary...\n")
cat("=========================================\n")
cat(" PIPELINE COMPLETE!\n")
cat("=========================================\n")
cat(sprintf(" Figures: %s/\n", FIGURE_DIR))
cat(sprintf(" Tables:  %s/\n", TABLE_DIR))

# list.files() lists all files in a directory.
# pattern = "\\.(png|pdf)$" means "files ending in .png or .pdf"
# The \\ is R's way of escaping the dot (a dot means "any character" in regex).
# $ means "end of string" — so we match ".png" at the very end of the filename.
fig_files <- list.files(FIGURE_DIR, pattern = "\\.(png|pdf)$")
tab_files <- list.files(TABLE_DIR, pattern = "\\.csv$")

# length() counts how many items are in a vector.
cat(sprintf("\n  %d figures generated\n", length(fig_files)))
cat(sprintf("  %d tables generated\n", length(tab_files)))
cat("=========================================\n")

# ---- Stop logging ----
# Close the log file and print the final summary (including "COMPLETED" status).
# The log file is now complete and can be shared for review.
stop_logging(success = TRUE)
