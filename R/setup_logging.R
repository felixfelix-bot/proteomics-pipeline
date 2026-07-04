###############################################################################
# setup_logging.R
# Sets up logging: all stdout/stderr goes to BOTH the console AND a log file.
# The log filename includes timestamp + git commit hash for traceability.
#
# HOW IT WORKS:
# R has a function called sink() that "redirects" output. Normally, cat()
# and print() send text to your screen (the console). sink() can redirect
# that text to a file instead. With split = TRUE, it sends text to BOTH
# the file AND the screen simultaneously — so you can watch the output
# live AND have it saved for later.
#
# This is important because the pipeline produces a LOT of output, and
# having a saved log lets you:
#   1. Review what happened after the fact
#   2. Share the log with collaborators for debugging
#   3. Know exactly which code version (git commit) produced these results
#
# Two functions are defined here:
#   setup_logging() — call at the START to begin capturing output
#   stop_logging()  — call at the END to close the log file properly
###############################################################################

# ---- setup_logging: Begin capturing all output to a log file ----
# Parameters:
#   script_name: A label for the log filename (e.g., "run_all", "volcano")
#
# Returns: The path to the log file (invisibly)
setup_logging <- function(script_name = "run") {
  # ---- Get the git commit hash ----
  # system2() runs an external command (like subprocess in Python).
  # Here we run "git rev-parse --short HEAD" to get the short commit hash.
  # stdout = TRUE captures the output, stderr = FALSE ignores errors.
  # trimws() removes leading/trailing whitespace from the result.
  #
  # tryCatch() is R's try/catch — if git isn't available (or this isn't
  # a git repo), it returns "unknown" instead of crashing.
  git_commit <- tryCatch(
    trimws(system2("git", c("rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE)),
    error = function(e) "unknown"
  )

  # Also get the branch name (e.g., "main", "develop")
  git_branch <- tryCatch(
    trimws(system2("git", c("rev-parse", "--abbrev-ref", "HEAD"),
                   stdout = TRUE, stderr = FALSE)),
    error = function(e) "unknown"
  )

  # ---- Create a timestamped log filename ----
  # Sys.time() returns the current date and time.
  # format() converts it to a string in the format: 20260705_144530
  # (YearMonthDay_HourMinuteSecond)
  timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Create the logs directory if it doesn't exist
  log_dir <- file.path("output", "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # Build the log filename: timestamp_commit_script.log
  # Example: 20260705_144530_a1b2c3d_run_all.log
  log_file <- file.path(log_dir,
    sprintf("%s_%s_%s.log", timestamp_str, git_commit, script_name))

  # ---- Write a header to the log file ----
  # This header goes at the very top of the log file with metadata.
  # c() combines strings into a vector. writeLines() writes them to the file.
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

  # ---- Set up output redirection (sink) ----
  # file() opens a "connection" (like opening a file handle in Python).
  #   open = "at" means "append text mode" — add to end of file.
  # We open TWO connections: one for normal output, one for messages/warnings.
  log_con <- file(log_file, open = "at")
  msg_con <- file(log_file, open = "at")

  # sink() redirects R's output stream.
  #   sink(log_con, split = TRUE):
  #     All output from cat(), print(), etc. goes to log_con (the file)
  #     AND also to the console (split = TRUE means "send to both").
  #
  #   sink(msg_con, type = "message"):
  #     Redirects the MESSAGE stream (warnings, package messages) to the
  #     same file. type = "message" targets stderr instead of stdout.
  sink(log_con, split = TRUE)
  sink(msg_con, type = "message")

  # ---- Store connection handles globally for cleanup later ----
  # .GlobalEnv is R's "global environment" — the main workspace.
  # We store the connections and metadata there so stop_logging() can
  # access them later (functions in R don't automatically see variables
  # from the calling scope).
  .GlobalEnv$.log_con <- log_con
  .GlobalEnv$.msg_con <- msg_con
  .GlobalEnv$LOG_FILE <- log_file
  .GlobalEnv$LOG_COMMIT <- git_commit
  .GlobalEnv$LOG_TIMESTAMP <- timestamp_str

  # Print a message so the user knows logging started
  cat("[Logging] All output is being saved to:\n")
  cat(sprintf("          %s\n", log_file))
  cat(sprintf("          (commit: %s)\n\n", git_commit))

  # invisible() returns a value without printing it to the console.
  # We return the log file path in case the caller wants it.
  return(invisible(log_file))
}

# ---- stop_logging: Close the log file and print final status ----
# Parameters:
#   success: TRUE if the pipeline completed, FALSE if it crashed
#
# This function MUST be called at the end to properly close the log file.
# If it's not called, the file might be incomplete.
stop_logging <- function(success = TRUE) {
  # ---- Close the sink connections ----
  # sink() with no arguments STOPS redirection (returns output to console).
  # type = "message" stops the message/warning redirection first.
  # We wrap in try() in case the sinks were already closed.
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)

  # ---- Close the file connections ----
  # exists() checks if a variable exists. We check before closing to
  # avoid errors if setup_logging() was never called.
  if (exists(".log_con", envir = .GlobalEnv)) {
    try(close(.GlobalEnv$.log_con), silent = TRUE)
    # rm() removes a variable from the environment
    rm(".log_con", envir = .GlobalEnv)
  }
  if (exists(".msg_con", envir = .GlobalEnv)) {
    try(close(.GlobalEnv$.msg_con), silent = TRUE)
    rm(".msg_con", envir = .GlobalEnv)
  }

  # ---- Print the final summary ----
  # This prints AFTER sinks are closed, so it goes ONLY to the console
  # (not to the log file, which is already being written to separately below).
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

    # ---- Also append the status to the log file itself ----
    # cat() with file = ... and append = TRUE writes directly to the file
    # (bypassing the now-closed sink).
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
