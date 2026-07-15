#!/bin/bash
# Bootstrap script: Install R + proteomics pipeline deps via apt prebuilt
# NO conda, NO source compilation of large packages, NO ansible.
#
# ~90% of packages are apt prebuilt (zero compilation).
# 4 small Bioconductor packages (clusterProfiler, enrichplot, rrvgo,
# EnhancedVolcano) compile from source — their heavy deps are already
# via apt, so it's fast.
#
# Usage:
#   chmod +x bootstrap.sh && ./bootstrap.sh
#
# Tested on: Ubuntu 22.04 / 24.04 / 26.04
# Time: ~5-10 min (downloads prebuilt .deb files)

# NOTE: No 'set -e' — we want to continue even if one package fails,
# and report all failures at the end.

echo "=========================================="
echo " Proteomics Pipeline — Fast Install"
echo "=========================================="
echo ""

# ---- Step 1: apt prebuilt packages ----
echo "[1/3] Installing R + prebuilt packages via apt..."
echo "      (sudo password needed for apt only)"

sudo apt-get update -y

# R base + build deps + LaTeX (for PDF output)
sudo apt-get install -y \
  r-base r-base-dev \
  libxml2-dev libcurl4-openssl-dev libssl-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev \
  libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev libxt-dev \
  libbz2-dev libzstd-dev liblzma-dev libgit2-dev libssh2-1-dev \
  libgdal-dev libgeos-dev libproj-dev libudunits2-dev \
  libglpk-dev libgmp-dev libmpfr-dev \
  texlive-latex-base texlive-latex-extra texlive-fonts-recommended \
  texinfo qpdf pandoc make git curl wget

# CRAN packages (apt prebuilt — no compilation)
sudo apt-get install -y \
  r-cran-ggplot2 r-cran-ggrepel r-cran-scales r-cran-rcolorbrewer \
  r-cran-venndiagram r-cran-upsetr r-cran-dplyr r-cran-readr \
  r-cran-tibble r-cran-tidyr r-cran-gprofiler2 r-cran-here \
  r-cran-stringr r-cran-patchwork r-cran-igraph

# Bioconductor packages (apt prebuilt — no compilation)
sudo apt-get install -y \
  r-bioc-org.hs.eg.db r-bioc-stringdb r-bioc-complexheatmap \
  r-bioc-gosemsim r-bioc-delayedarray r-bioc-sparsearray

# Remove apt BiocVersion pin (it pins to 3.20, but R 4.5 needs Bioc 3.22)
# This lets BiocManager auto-detect the correct version.
sudo apt-get remove -y r-bioc-biocversion 2>/dev/null || true

echo ""
echo "[1/3] apt prebuilt packages done."
echo ""

# ---- Step 2: Create user R library ----
echo "[2/3] Setting up R library..."

R_LIBDIR="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
mkdir -p "$R_LIBDIR"
export R_LIBS_USER="$R_LIBDIR"

R_VERSION=$(Rscript -e 'cat(paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))' 2>/dev/null)
echo "      Detected R version: $R_VERSION"

echo ""

# ---- Step 3: Install remaining packages via Rscript ----
echo "[3/3] Installing remaining packages via Rscript..."
echo "      (ggVennDiagram, config from CRAN)"
echo "      (EnhancedVolcano, clusterProfiler, enrichplot, rrvgo from Bioconductor)"
echo ""

Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

# ---- CRAN packages ----
cran_pkgs <- c("ggVennDiagram", "config")
for (p in cran_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (CRAN)\n", p))
    install.packages(p)
  }
}

# ---- Bioconductor packages ----
# ALL of these are Bioconductor, NOT CRAN.
# EnhancedVolcano was incorrectly tried from CRAN before — that fails.
bioc_pkgs <- c("EnhancedVolcano", "clusterProfiler", "enrichplot", "rrvgo")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Determine correct Bioconductor version for this R
# R 4.5 -> Bioc 3.22, R 4.4 -> Bioc 3.20
r_short <- paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
bioc_ver <- switch(r_short,
  "4.5" = "3.22",
  "4.4" = "3.20",
  "4.3" = "3.18",
  as.character(BiocManager::version())  # fallback: auto-detect
)
cat(sprintf("  R %s -> Bioconductor %s\n", r_short, bioc_ver))

# Clean stale 00LOCK dirs before installing
lib_dir <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", r_short)
if (dir.exists(lib_dir)) {
  locks <- list.files(lib_dir, pattern = "^00LOCK", full.names = TRUE)
  if (length(locks) > 0) {
    cat("  Cleaning stale 00LOCK directories...\n")
    unlink(locks, recursive = TRUE)
  }
}

for (p in bioc_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (Bioc %s)\n", p, bioc_ver))
    BiocManager::install(p, version = bioc_ver, update = FALSE, ask = FALSE, force = TRUE)
  }
}

# ---- Verify all packages ----
all_pkgs <- c(cran_pkgs, bioc_pkgs,
  "ggplot2", "ggrepel", "scales", "RColorBrewer", "VennDiagram", "UpSetR",
  "dplyr", "readr", "tibble", "tidyr", "gprofiler2", "here", "stringr",
  "patchwork", "igraph", "org.Hs.eg.db", "STRINGdb", "ComplexHeatmap",
  "GOSemSim", "DOSE")

failed <- c()
for (p in all_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [FAIL] %s\n", p))
    failed <- c(failed, p)
  } else {
    cat(sprintf("  [OK]   %s %s\n", p, as.character(packageVersion(p))))
  }
}

if (length(failed) > 0) {
  cat("\n  [FAILED]:", paste(failed, collapse = ", "), "\n")
  cat("  Try running this script again — it skips already-installed packages.\n")
} else {
  cat("\n  [SUCCESS] All", length(all_pkgs), "packages installed\n")
}
'

echo ""
echo "=========================================="
echo " Install complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Copy CSV data files into data/"
echo "  2. make check        (verify all packages)"
echo "  3. make test         (test on synthetic data)"
echo "  4. make all          (run full pipeline on real data)"
echo ""