###############################################################################
# list_data.R
# Lists ALL CSV files in the data/ directory (recursively) with full paths,
# file sizes, and row counts. Also shows the header (first line) of each
# file so the AI can understand the column structure without seeing data.
#
# Output is designed to be pasted directly into the AI chat.
# NO data values are printed — only file paths, sizes, and headers.
#
# Usage:
#   make list-data
###############################################################################
source("R/01_config.R")

cat("\n=========================================\n")
cat(" DATA DIRECTORY — CSV FILE INVENTORY\n")
cat("=========================================\n\n")

# Find ALL CSV files recursively in data/
all_csvs <- list.files(DATA_DIR, pattern = "\\.csv$",
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)

if (length(all_csvs) == 0) {
  cat("No CSV files found in:", DATA_DIR, "\n")
  cat("Make sure your mass spec data is in subdirectories under data/\n")
  quit(status = 0)
}

cat(sprintf("Found %d CSV files:\n\n", length(all_csvs)))

for (i in seq_along(all_csvs)) {
  filepath <- all_csvs[i]
  rel_path <- filepath
  # Show relative to project root for readability
  rel_path <- sub(paste0(getwd(), "/"), "", filepath)

  # Get file size
  size_bytes <- file.size(filepath)
  size_mb <- round(size_bytes / 1024 / 1024, 2)

  # Read header line only (first line)
  header_line <- tryCatch({
    readLines(filepath, n = 1)
  }, error = function(e) {
    "(unreadable)"
  })

  # Count rows (minus header)
  n_rows <- tryCatch({
    length(count.fields(filepath, sep = ",")) - 1
  }, error = function(e) {
    NA
  })

  cat(sprintf("[%d] %s\n", i, rel_path))
  cat(sprintf("    Size: %.2f MB | Rows: %s\n", size_mb,
              format(n_rows, big.mark = ",")))
  cat(sprintf("    Header: %s\n\n", header_line))
}

cat("=========================================\n")
cat(sprintf("Total: %d CSV files\n", length(all_csvs)))
cat("=========================================\n")
cat("\nPaste the above into the AI chat to populate its context\n")
cat("with your data file structure. No data values are shown —\n")
cat("only paths, sizes, row counts, and column headers.\n")
