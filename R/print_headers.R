###############################################################################
# print_headers.R
# Prints just the column names of every CSV file in data/.
# No file sizes, no row counts — just: filename → column names.
#
# Usage:
#   make headers
###############################################################################
source("R/01_config.R")

cat("\n=========================================\n")
cat(" CSV COLUMN HEADERS\n")
cat("=========================================\n\n")

all_csvs <- list.files(DATA_DIR, pattern = "\\.csv$",
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)

if (length(all_csvs) == 0) {
  cat("No CSV files found in:", DATA_DIR, "\n")
  quit(status = 0)
}

for (filepath in all_csvs) {
  rel_path <- sub(paste0(getwd(), "/"), "", filepath)
  cols <- tryCatch({
    header_line <- readLines(filepath, n = 1)
    # Parse CSV header to extract clean column names
    # Remove surrounding quotes, split by comma
    cols <- strsplit(header_line, ",")[[1]]
    cols <- gsub('^"|"$', "", cols)  # Strip quotes
    cols <- trimws(cols)             # Strip whitespace
    cols
  }, error = function(e) "(unreadable)")

  cat(sprintf("%s\n", rel_path))
  cat(sprintf("  → %s\n\n", paste(cols, collapse = " | ")))
}

cat(sprintf("Total: %d files\n", length(all_csvs)))
