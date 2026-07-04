###############################################################################
# 00_install_packages.R
# Installs all R packages required for the proteomics analysis pipeline.
# Safe to run multiple times — only installs missing packages.
#
# Run from R console or command line:
#   Rscript 00_install_packages.R
###############################################################################

cat("\n")
cat("=========================================\n")
cat(" Proteomics Pipeline — Package Installer\n")
cat("=========================================\n\n")

# ---- Helper function ----
is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_if_missing <- function(pkg, type = "cran") {
  if (is_installed(pkg)) {
    cat(sprintf("  [OK]      %s (already installed)\n", pkg))
    return(invisible(FALSE))
  }

  cat(sprintf("  [INSTALL] %s (%s)...\n", pkg, type))
  tryCatch({
    if (type == "bioc") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    cat(sprintf("           -> %s installed successfully\n", pkg))
    invisible(TRUE)
  }, error = function(e) {
    cat(sprintf("           -> FAILED: %s\n", conditionMessage(e)))
    invisible(FALSE)
  })
}

# ---- Install BiocManager first ----
cat("[1] Checking BiocManager...\n")
if (!is_installed("BiocManager")) {
  install.packages("BiocManager")
}
cat(sprintf("  BiocManager version: %s\n\n", as.character(packageVersion("BiocManager"))))

# ---- CRAN packages ----
cat("[2] Installing CRAN packages...\n\n")

cran_packages <- c(
  # Plotting
  "ggplot2",         # Core plotting engine
  "ggrepel",         # Non-overlapping text labels for plots
  "scales",          # Axis formatting
  "RColorBrewer",    # Color palettes

  # Volcano plots
  "EnhancedVolcano", # Publication-ready volcano plots (also on Bioconductor)

  # Venn diagrams
  "ggVennDiagram",   # 2-7 set Venn diagrams
  "VennDiagram",     # Classic Venn (backup)
  "UpSetR",          # UpSet plots for >4 sets

  # Data wrangling
  "dplyr",           # Data manipulation
  "readr",           # Fast CSV reading
  "tibble",          # Modern data frames
  "tidyr",           # Data reshaping

  # GO enrichment (online fallback)
  "gprofiler2",      # g:Profiler enrichment (queries web server)

  # Utility
  "config",          # Configuration management
  "here"             # Project-relative paths
)

for (pkg in cran_packages) {
  install_if_missing(pkg, "cran")
}

# ---- Bioconductor packages ----
cat("\n[3] Installing Bioconductor packages...\n\n")

bioc_packages <- c(
  "clusterProfiler",  # GO/KEGG enrichment — the gold standard
  "org.Hs.eg.db",     # Human annotation database
  "enrichplot",       # Visualization for enrichment results
  "DOSE",             # Disease ontology (dependency, explicit for clarity)
  "STRINGdb",         # STRING protein interaction network
  "rrvgo",            # Reduce GO term redundancy
  "ComplexHeatmap"    # Heatmap visualization
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, "bioc")
}

# ---- Verify all packages load ----
cat("\n[4] Verifying all packages can be loaded...\n\n")

all_packages <- c(cran_packages, bioc_packages)
loaded <- c()
failed <- c()

for (pkg in all_packages) {
  if (is_installed(pkg)) {
    loaded <- c(loaded, pkg)
  } else {
    failed <- c(failed, pkg)
  }
}

cat(sprintf("  Successfully loaded: %d/%d packages\n", length(loaded), length(all_packages)))

if (length(failed) > 0) {
  cat("\n  [WARNING] The following packages FAILED to install:\n")
  for (f in failed) {
    cat(sprintf("    - %s\n", f))
  }
  cat("\n  Try installing manually:\n")
  cat("    BiocManager::install(c(\"", paste(failed, collapse = '", "'), "\"))\n", sep = "")
} else {
  cat("\n  [SUCCESS] All packages installed and verified!\n")
}

cat("\n=========================================\n")
cat(" Installation complete.\n")
cat("=========================================\n\n")
