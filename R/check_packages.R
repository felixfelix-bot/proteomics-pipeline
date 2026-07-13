#############################################################################
# check_packages.R — Verify all required packages are installed and loadable.
# Used by: make check
#############################################################################

pkgs <- c(
  "ggplot2", "ggrepel", "EnhancedVolcano", "ggVennDiagram",
  "VennDiagram", "UpSetR", "dplyr", "readr", "tibble", "tidyr",
  "gprofiler2", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "STRINGdb", "ComplexHeatmap", "igraph"
)

# Optional packages — warn but don't block the pipeline
optional_pkgs <- c("rrvgo")

cat("\n")
cat("=========================================\n")
cat(" Package Check\n")
cat("=========================================\n\n")

loaded <- sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))

for (p in names(loaded)) {
  if (loaded[[p]]) {
    cat("  [OK]      ", p, "\n", sep = "")
  } else {
    cat("  [MISSING] ", p, "\n", sep = "")
  }
}

cat(sprintf("\n  %d/%d required packages installed\n", sum(loaded), length(loaded)))

# Check optional packages (warn but don't fail)
opt_loaded <- sapply(optional_pkgs, function(p) requireNamespace(p, quietly = TRUE))
for (p in names(opt_loaded)) {
  if (!opt_loaded[[p]]) {
    cat(sprintf("  [OPTIONAL] %s not installed (not needed for most analyses)\n", p))
  }
}

if (!all(loaded)) {
  cat("\n  Some packages are missing. Run:\n")
  cat("    make install\n")
  cat("  or:\n")
  cat("    Rscript R/00_install_packages.R\n\n")
  quit(status = 1)
} else {
  cat("  All packages ready!\n\n")
}
