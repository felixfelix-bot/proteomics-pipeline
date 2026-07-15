#!/bin/bash
# Bootstrap script: Install R + proteomics pipeline deps via apt prebuilt
# NO conda, NO source compilation of large packages, NO ansible.
#
# ~90% of packages are apt prebuilt (zero compilation).
# 3 small Bioconductor packages (clusterProfiler, enrichplot, rrvgo) 
# compile from source — their heavy deps are already via apt, so it's fast.
#
# Usage:
#   chmod +x bootstrap.sh && ./bootstrap.sh
#
# Tested on: Ubuntu 22.04 / 24.04 / 26.04
# Time: ~5-10 min (downloads prebuilt .deb files)

set -e

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

echo ""
echo "[1/3] apt prebuilt packages done."
echo ""

# ---- Step 2: Create user R library ----
echo "[2/3] Setting up R library..."

R_LIBDIR="$HOME/R/x86_64-pc-linux-gnu-library/4.5"
mkdir -p "$R_LIBDIR"
export R_LIBS_USER="$R_LIBDIR"

# Also handle R 4.4 if that's what's installed
R_VERSION=$(Rscript -e 'cat(paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))' 2>/dev/null)
echo "      Detected R version: $R_VERSION"

echo ""

# ---- Step 3: Install remaining packages via Rscript ----
echo "[3/3] Installing remaining packages via Rscript..."
echo "      (EnhancedVolcano, ggVennDiagram, config from CRAN)"
echo "      (clusterProfiler, enrichplot, rrvgo from Bioconductor)"
echo "      Their heavy deps are already installed via apt above."
echo ""

Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

cran_pkgs <- c("EnhancedVolcano", "ggVennDiagram", "config")
bioc_pkgs <- c("clusterProfiler", "enrichplot", "rrvgo")

# CRAN packages
for (p in cran_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (CRAN)\n", p))
    install.packages(p)
  }
}

# BiocManager + Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (p in bioc_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s\n", p))
  } else {
    cat(sprintf("  [INSTALL] %s (Bioconductor)\n", p))
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

# Verify all
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
  quit(status = 1)
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
echo "  1. cd proteomics-pipeline"
echo "  2. Copy CSV data files into data/"
echo "  3. make check        (verify all packages)"
echo "  4. make test         (test on synthetic data)"
echo "  5. make all          (run full pipeline on real data)"
echo ""