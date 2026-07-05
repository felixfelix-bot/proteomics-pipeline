###############################################################################
# list_data.R
# Lists ALL CSV AND Excel files in the data/ directory (recursively) with
# full paths, file sizes, and row counts. Also shows the header (first
# line) of each file so the AI can understand the column structure.
#
# If any Excel files (.xlsx, .xls) are found, they are automatically
# exported as CSV alongside the original, so the pipeline can consume them.
#
# Output is designed to be pasted directly into the AI chat.
# NO data values are printed — only file paths, sizes, and headers.
#
# Usage:
#   make list-data
###############################################################################
source("R/01_config.R")

cat("\n=========================================\n")
cat(" DATA DIRECTORY — FILE INVENTORY\n")
cat("=========================================\n\n")

# ---- Step 1: Find and convert any Excel files ----
excel_files <- list.files(DATA_DIR, pattern = "\\.(xlsx|xls)$",
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = TRUE)

if (length(excel_files) > 0) {
  cat(sprintf("Found %d Excel file(s) — converting to CSV...\n\n", length(excel_files)))

  has_readxl <- requireNamespace("readxl", quietly = TRUE)
  if (!has_readxl) {
    cat("  WARNING: 'readxl' package not installed. Cannot auto-convert Excel.\n")
    cat("  Please export Excel files as CSV manually, or install readxl:\n")
    cat("    Rscript -e \"install.packages('readxl')\"\n\n")
  } else {
    library(readxl)
    for (xls_path in excel_files) {
      cat(sprintf("  Converting: %s\n", basename(xls_path)))

      # Get all sheet names
      sheets <- tryCatch(excel_sheets(xls_path), error = function(e) {
        cat(sprintf("    ERROR reading sheets: %s\n", conditionMessage(e)))
        return(NULL)
      })

      if (is.null(sheets)) next

      cat(sprintf("    Sheets: %s\n", paste(sheets, collapse = ", ")))

      # Convert each sheet to CSV (first sheet gets no suffix, others get _sheetname)
      for (i in seq_along(sheets)) {
        sheet_name <- sheets[i]
        df <- tryCatch(read_excel(xls_path, sheet = i), error = function(e) {
          cat(sprintf("    ERROR reading sheet '%s': %s\n", sheet_name, conditionMessage(e)))
          return(NULL)
        })
        if (is.null(df)) next

        # Build CSV filename next to the Excel file
        xls_dir <- dirname(xls_path)
        xls_base <- sub("\\.(xlsx|xls)$", "", basename(xls_path), ignore.case = TRUE)

        if (length(sheets) == 1) {
          csv_name <- paste0(xls_base, ".csv")
        } else {
          safe_sheet <- gsub("[^A-Za-z0-9_]", "_", sheet_name)
          csv_name <- paste0(xls_base, "_", safe_sheet, ".csv")
        }
        csv_path <- file.path(xls_dir, csv_name)

        write.csv(df, csv_path, row.names = FALSE)
        cat(sprintf("    → %s (%d rows, %d cols)\n", csv_name, nrow(df), ncol(df)))
      }
    }
    cat("\n")
  }
} else {
  cat("No Excel files found.\n\n")
}

# ---- Step 2: List ALL CSV files ----
all_csvs <- list.files(DATA_DIR, pattern = "\\.csv$",
                        full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)

if (length(all_csvs) == 0) {
  cat("No CSV files found in:", DATA_DIR, "\n")
  quit(status = 0)
}

cat(sprintf("Found %d CSV files:\n\n", length(all_csvs)))

for (i in seq_along(all_csvs)) {
  filepath <- all_csvs[i]
  rel_path <- sub(paste0(getwd(), "/"), "", filepath)

  size_bytes <- file.size(filepath)
  size_mb <- round(size_bytes / 1024 / 1024, 2)

  header_line <- tryCatch({
    readLines(filepath, n = 1)
  }, error = function(e) "(unreadable)")

  n_rows <- tryCatch({
    length(count.fields(filepath, sep = ",")) - 1
  }, error = function(e) NA)

  cat(sprintf("[%d] %s\n", i, rel_path))
  cat(sprintf("    Size: %.2f MB | Rows: %s\n", size_mb,
              format(n_rows, big.mark = ",")))
  cat(sprintf("    Header: %s\n\n", header_line))
}

# ---- Step 3: List remaining Excel files (for reference) ----
if (length(excel_files) > 0) {
  cat("Also found Excel files (converted to CSV above):\n\n")
  for (f in excel_files) {
    rel <- sub(paste0(getwd(), "/"), "", f)
    cat(sprintf("  %s (%.2f MB)\n", rel, round(file.size(f) / 1024 / 1024, 2)))
  }
  cat("\n")
}

cat("=========================================\n")
cat(sprintf("Total: %d CSV files, %d Excel files\n", length(all_csvs), length(excel_files)))
cat("=========================================\n")
cat("\nPaste the above into the AI chat to populate its context\n")
cat("with your data file structure. No data values are shown —\n")
cat("only paths, sizes, row counts, and column headers.\n")
