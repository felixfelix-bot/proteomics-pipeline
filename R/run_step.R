###############################################################################
# run_step.R
# Runs a single analysis step with full logging.
#
# Usage:
#   Rscript R/run_step.R <step_name>
#
# Where step_name is one of:
#   volcano, venn, go, string, families, overlap
###############################################################################

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Usage: Rscript R/run_step.R <step_name>\n")
    cat("  Steps: volcano, venn, go, string, families, overlap\n")
    quit(status = 1)
  }

  step <- args[1]

  valid_steps <- c("volcano", "venn", "go", "string", "families", "overlap",
                   "targeted_volcanos", "flagip_volcano", "targeted_venns",
                   "targeted_go", "go_network_volcano", "chx_crac_analysis",
                   "venn_examples", "venn_label_examples", "string_network",
                   "crac_string_network", "venn_overflow_examples", "gsea",
                   "pathway_network", "lydia_network_volcano",
                   "chx_common_analysis")
  if (!(step %in% valid_steps)) {
    cat(sprintf("Error: Unknown step '%s'\n", step))
    cat(sprintf("  Valid steps: %s\n", paste(valid_steps, collapse = ", ")))
    quit(status = 1)
  }

  step_scripts <- list(
    volcano  = "R/02_volcano_plots.R",
    venn     = "R/03_venn_diagrams.R",
    go       = "R/04_go_enrichment.R",
    string   = "R/05_string_network.R",
    families = "R/06_gene_families.R",
    overlap  = "R/07_overlap_analysis.R",
    targeted_volcanos = "R/08_targeted_volcanos.R",
    flagip_volcano    = "R/09_flagip_volcano.R",
    targeted_venns    = "R/10_targeted_venns.R",
    targeted_go       = "R/11_targeted_go.R",
    go_network_volcano = "R/12_go_network_volcano.R",
    chx_crac_analysis  = "R/13_chx_crac_analysis.R",
    venn_examples      = "R/14_venn_examples.R",
    venn_label_examples = "R/15_venn_label_examples.R",
    string_network     = "R/16_string_network_targeted.R",
    crac_string_network = "R/17_crac_string_network.R",
    venn_overflow_examples = "R/18_venn_overflow_examples.R",
    gsea                   = "R/19_gsea.R",
    pathway_network        = "R/20_pathway_network.R",
    lydia_network_volcano  = "R/21_lydia_network_volcano.R",
    chx_common_analysis    = "R/24_chx_common_analysis.R"
  )

  # Set up logging (captures all output to file + console)
  sys.source("R/setup_logging.R", envir = .GlobalEnv)
  setup_logging(script_name = step)

  # Load config and utils into the GLOBAL environment
  # This ensures all functions are visible to subsequent source() calls.
  # source() with local=TRUE evaluates in parent.frame() (main's env),
  # but subsequent source() calls create sibling environments that can't
  # see each other. Using sys.source() into .GlobalEnv avoids this.
  project_root <- getwd()
  Sys.setenv(PROJECT_ROOT = project_root)
  sys.source(file.path(project_root, "R", "01_config.R"), envir = .GlobalEnv)
  sys.source(file.path(project_root, "R", "utils.R"), envir = .GlobalEnv)

  cat("\n=========================================\n")
  cat(sprintf(" Running step: %s\n", step))
  cat("=========================================\n\n")

  # Run the step — also source into .GlobalEnv so it can see all functions
  success <- TRUE
  tryCatch({
    sys.source(file.path(project_root, step_scripts[[step]]), envir = .GlobalEnv)
  }, error = function(e) {
    success <<- FALSE
    cat("\n\nFATAL ERROR in step '", step, "':\n", sep = "")
    cat("  ", conditionMessage(e), "\n", sep = "")
  })

  stop_logging(success = success)

  if (!success) quit(status = 1)
}

main()
