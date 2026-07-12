#!/usr/bin/env Rscript
# clean_old_outputs.R
# Deletes output files that don't contain the current git commit hash.
# Called by Makefile before each analysis step.
# Cross-platform (works on Windows, macOS, Linux).

# Get current git hash
current_hash <- tryCatch({
  hash <- system("git rev-parse --short HEAD", intern = TRUE)
  trimws(hash)
}, error = function(e) "nogit")

cat(sprintf("  Cleaning old outputs (keeping hash: %s)...\n", current_hash))

# Directories to clean
dirs_to_clean <- c(
  file.path("output", "figures"),
  file.path("output", "tables")
)

# File extensions to consider
extensions <- c("png", "pdf", "csv", "txt")

removed <- 0
kept <- 0

for (dir_path in dirs_to_clean) {
  if (!dir.exists(dir_path)) next

  files <- list.files(dir_path, full.names = TRUE, recursive = FALSE)
  files <- files[tools::file_ext(files) %in% extensions]

  for (f in files) {
    if (grepl(current_hash, f, fixed = TRUE)) {
      kept <- kept + 1
    } else {
      file.remove(f)
      removed <- removed + 1
    }
  }
}

cat(sprintf("  Removed %d old file(s), kept %d current file(s).\n", removed, kept))
