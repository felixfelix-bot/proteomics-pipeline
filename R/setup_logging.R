###############################################################################
# setup_logging.R
# Sets up logging: all stdout/stderr goes to BOTH the console AND a log file.
# The log filename includes timestamp + git commit hash for traceability.
#
# Call this at the TOP of run_all.R or run_step.R, BEFORE sourcing any
# analysis scripts. The sink() must be stopped at the end (handled by
# stop_logging() or on.exit()).
#
# After sourcing, these global variables are available:
#   LOG_FILE    - path to the log file
#   LOG_COMMIT  - git commit hash
###############################################################################

setup_logging <- function(script_name = "run") {
  # Get git info
  git_commit <- tryCatch(
    trimws(system2("git", c("rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE)),
    error = function(e) "unknown"
  )
  git_branch <- tryCatch(
    trimws(system2("git", c("rev-parse", "--abbrev-ref", "HEAD"),
                   stdout = TRUE, stderr = FALSE)),
    error = function(e) "unknown"
  )

  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_dir <- file.path("output", "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(log_dir,
    sprintf("%s_%s_%s.log", timestamp_str, git_commit, script_name))

  # Write header to log file
  header <- c(
    "=========================================",
    "PROTEOMICS PIPELINE RUN LOG",
    "=========================================",
    sprintf("Script:      %s", script_name),
    sprintf("Timestamp:   %s", timestamp_str),
    sprintf("Git commit:  %s", git_commit),
    sprintf("Git branch:  %s", git_branch),
    sprintf("R version:   %s", R.version.string),
    sprintf("Working dir: %s", getwd()),
    sprintf("Platform:    %s", R.version$platform),
    "=========================================",
    ""
  )
  writeLines(header, log_file)

  # Open connections for sink
  log_con <- file(log_file, open = "at")  # append text
  msg_con <- file(log_file, open = "at")  # separate con for messages

  # sink with split=TRUE: output goes to BOTH console and file
  sink(log_con, split = TRUE)
  sink(msg_con, type = "message")  # capture messages/warnings too

  # Store handles for cleanup
  .GlobalEnv$.log_con <- log_con
  .GlobalEnv$.msg_con <- msg_con
  .GlobalEnv$LOG_FILE <- log_file
  .GlobalEnv$LOG_COMMIT <- git_commit
  .GlobalEnv$LOG_TIMESTAMP <- timestamp_str

  cat("[Logging] All output is being saved to:\n")
  cat(sprintf("          %s\n", log_file))
  cat(sprintf("          (commit: %s)\n\n", git_commit))

  return(invisible(log_file))
}

stop_logging <- function(success = TRUE) {
  # Flush and close sinks
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)

  # Close connections
  if (exists(".log_con", envir = .GlobalEnv)) {
    try(close(.GlobalEnv$.log_con), silent = TRUE)
    rm(".log_con", envir = .GlobalEnv)
  }
  if (exists(".msg_con", envir = .GlobalEnv)) {
    try(close(.GlobalEnv$.msg_con), silent = TRUE)
    rm(".msg_con", envir = .GlobalEnv)
  }

  # Print summary to console (after sinks are closed, so it only goes to console)
  if (exists("LOG_FILE", envir = .GlobalEnv)) {
    cat("\n=========================================\n")
    cat(sprintf("Log saved to: %s\n", .GlobalEnv$LOG_FILE))
    cat(sprintf("Git commit:   %s\n", .GlobalEnv$LOG_COMMIT))
    if (success) {
      cat("Status:       COMPLETED\n")
    } else {
      cat("Status:       FAILED\n")
    }
    cat("=========================================\n")

    # Also append status to the log file
    try({
      cat("\n=========================================\n",
          file = .GlobalEnv$LOG_FILE, append = TRUE)
      cat(sprintf("Status: %s\n", ifelse(success, "COMPLETED", "FAILED")),
          file = .GlobalEnv$LOG_FILE, append = TRUE)
      cat("=========================================\n",
          file = .GlobalEnv$LOG_FILE, append = TRUE)
    }, silent = TRUE)
  }
}
