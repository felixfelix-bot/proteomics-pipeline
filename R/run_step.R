###############################################################################
# run_step.R
# Runs a single analysis step with full logging.
#
# Usage:
#   Rscript R/run_step.R <step_name> [--force]
#
# Where step_name is one of:
#   volcano, venn, go, string, families, overlap, ...
#
# Skip logic: if output/.stamps/<step>.stamp matches the current MD5 of
# the step script + config + utils, the step is skipped ("up to date").
# Use --force to override and always re-run.
###############################################################################

# ---- Compute a content hash for a step (script + shared deps) ----
# Returns NA if any file is missing.
compute_step_hash <- function(script_path) {
  deps <- c(script_path, "R/01_config.R", "R/utils.R", "R/setup_logging.R")
  hashes <- tools::md5sum(deps)
  paste(hashes, collapse = "")
}

# ---- Check if a step should be skipped (already completed, unchanged) ----
# Returns TRUE if the step can be skipped.
should_skip <- function(step, script_path) {
  stamp_dir <- file.path("output", ".stamps")
  stamp_file <- file.path(stamp_dir, step)

  if (!file.exists(stamp_file)) return(FALSE)

  current_hash <- compute_step_hash(script_path)
  if (is.na(current_hash)) return(FALSE)

  saved_hash <- tryCatch(
    trimws(readLines(stamp_file, n = 1, warn = FALSE)),
    error = function(e) ""
  )

  if (nchar(saved_hash) > 0 && saved_hash == current_hash) {
    return(TRUE)
  }
  return(FALSE)
}

# ---- Write a stamp after successful completion ----
write_stamp <- function(step, script_path) {
  stamp_dir <- file.path("output", ".stamps")
  dir.create(stamp_dir, recursive = TRUE, showWarnings = FALSE)
  stamp_file <- file.path(stamp_dir, step)
  current_hash <- compute_step_hash(script_path)
  writeLines(current_hash, stamp_file)
  cat(sprintf("[Cache] Stamp written: %s\n", stamp_file))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse --force flag
  force <- "--force" %in% args
  args <- args[args != "--force"]

  if (length(args) < 1) {
    cat("Usage: Rscript R/run_step.R <step_name> [--force]\n")
    cat("  Steps: volcano, venn, go, string, families, overlap, ...\n")
    cat("  --force: re-run even if output already exists\n")
    quit(status = 1)
  }

  step <- args[1]

  valid_steps <- c("volcano", "venn", "go", "string", "families", "overlap",
                   "targeted_volcanos", "flagip_volcano", "targeted_venns",
                   "targeted_go", "go_network_volcano", "chx_crac_analysis",
                   "venn_examples", "venn_label_examples", "string_network",
                   "crac_string_network", "venn_overflow_examples", "gsea",
                   "pathway_network", "lydia_network_volcano",
                   "ra_common", "bidirectional_go",
                   "chx_common_analysis", "shinygo_comparison",
                   "diagnostics", "network_go", "bidirectional_go_ra",
                   "string_style_network", "chx_volcano_venn",
                   "flagip_validated_go", "poster_figures",
                   "chx_kegg_crac_overlap", "dotplot_variants",
                   "style_variants")
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
    ra_common              = "R/22_ra_common_analysis.R",
    bidirectional_go       = "R/23_bidirectional_go.R",
    chx_common_analysis    = "R/24_chx_common_analysis.R",
    shinygo_comparison     = "R/25_shinygo_comparison.R",
    diagnostics            = "R/diagnostics.R",
    network_go             = "R/26_network_go_comparison.R",
    bidirectional_go_ra    = "R/27_bidirectional_go_ra.R",
    string_style_network   = "R/28_string_style_network.R",
    chx_volcano_venn       = "R/29_chx_volcano_venn.R",
    flagip_validated_go    = "R/30_flagip_validated_go.R",
    poster_figures         = "R/31_poster_figures.R",
    chx_kegg_crac_overlap  = "R/32_chx_kegg_crac_overlap.R",
    dotplot_variants       = "R/33_dotplot_variants.R",
    style_variants         = "R/34_style_variants.R"
  )

  # ---- Skip check: if step already ran with same code, skip it ----
  script_path <- step_scripts[[step]]
  if (!force && should_skip(step, script_path)) {
    cat(sprintf("[Cache] Skipping '%s' — up to date. Use --force to re-run.\n", step))
    quit(status = 0)
  }

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

  # ---- Write stamp on success so next run can skip ----
  if (success) {
    write_stamp(step, script_path)
  } else {
    quit(status = 1)
  }
}

main()
