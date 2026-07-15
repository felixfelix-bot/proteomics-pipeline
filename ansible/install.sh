#!/bin/bash
# ============================================================
# Proteomics Pipeline — Linux Quick Install (no conda, no ansible)
#
# Uses Ubuntu's prebuilt R packages for ~90% of dependencies.
# Only 3 small Bioconductor packages compile from source
# (clusterProfiler, enrichplot, rrvgo) — their deps are already
# prebuilt via apt, so it's fast.
#
# Usage:
#   chmod +x ansible/install.sh
#   ./ansible/install.sh
#
# Total time: ~5-10 min
# Requires: Ubuntu + sudo
# ============================================================

set -e

echo "=========================================="
echo " Proteomics Pipeline — Quick Install"
echo "=========================================="
echo ""

# ---- Step 1: apt prebuilt packages (no compilation) ----
echo "[1/3] Installing prebuilt R packages via apt..."
echo "  (sudo password may be required)"

sudo apt-get update -qq
sudo apt-get install -y -qq \
  r-base r-base-dev \
  libxml2-dev libcurl4-openssl-dev libssl-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libcairo2-dev libxt-dev libbz2-dev libzstd-dev liblzma-dev \
  libgit2-dev libssh2-1-dev libgdal-dev libgeos-dev libproj-dev \
  libudunits2-dev libglpk-dev libgmp-dev libmpfr-dev \
  texlive-latex-base texlive-latex-extra texlive-fonts-recommended \
  texinfo qpdf pandoc make git curl wget \
  \
  r-cran-ggplot2 r-cran-ggrepel r-cran-scales r-cran-rcolorbrewer \
  r-cran-venndiagram r-cran-upsetr r-cran-dplyr r-cran-readr \
  r-cran-tibble r-cran-tidyr r-cran-gprofiler2 r-cran-here \
  r-cran-stringr r-cran-patchwork r-cran-igraph \
  \
  r-bioc-org.hs.eg.db r-bioc-stringdb r-bioc-complexheatmap \
  r-bioc-gosemsim r-bioc-delayedarray r-bioc-sparsearray \
  2>&1 | tail -5

echo "  Done."
echo ""

# ---- Step 2: Install remaining packages via Rscript ----
echo "[2/3] Installing remaining Bioconductor packages via Rscript..."
echo "  (clusterProfiler, enrichplot, rrvgo, EnhancedVolcano, ggVennDiagram, config)"
echo "  These compile from source but their deps are already installed — fast."
echo ""

# Clean any stale lock files
find ~/R/x86_64-pc-linux-gnu-library/ -name '00LOCK*' -exec rm -rf {} + 2>/dev/null || true

Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- list(
  cran = c("ggVennDiagram", "config"),
  bioc = c("EnhancedVolcano", "clusterProfiler", "enrichplot", "rrvgo")
)

for (p in pkgs$cran) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (CRAN)\n", p))
    install.packages(p)
  }
}

for (p in pkgs$bioc) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (Bioconductor)\n", p))
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

failed <- c()
for (p in c(pkgs$cran, pkgs$bioc)) {
  if (!requireNamespace(p, quietly = TRUE)) failed <- c(failed, p)
}
if (length(failed) > 0) {
  cat("\n  [FAILED]:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("\n  [OK] All remaining packages installed\n")
}
'

echo ""

# ---- Step 3: Verify ----
echo "[3/3] Verifying all packages..."
echo ""

Rscript -e '
pkgs <- c("ggplot2","ggrepel","scales","RColorBrewer","ggVennDiagram",
          "VennDiagram","UpSetR","dplyr","readr","tibble","tidyr",
          "gprofiler2","config","here","stringr","patchwork","igraph",
          "clusterProfiler","enrichplot","DOSE","org.Hs.eg.db",
          "STRINGdb","rrvgo","ComplexHeatmap","EnhancedVolcano")
failed <- c()
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [FAIL] %s\n", p))
    failed <- c(failed, p)
  } else {
    cat(sprintf("  [OK]   %s\n", p))
  }
}
cat(sprintf("\n  %d/%d packages installed\n", length(pkgs) - length(failed), length(pkgs)))
if (length(failed) > 0) {
  cat("  Missing:", paste(failed, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("  All packages ready!\n")
}
'

echo ""
echo "=========================================="
echo " Install complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Copy CSV data files into data/"
echo "  2. make test   (test on synthetic data)"
echo "  3. make all    (run full pipeline)"
echo ""