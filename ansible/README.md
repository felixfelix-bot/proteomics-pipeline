# Proteomics Pipeline — Ansible Setup

Installs R, all system libraries, LaTeX, and R/Bioconductor packages
needed to run the TRIP4/ASCC proteomics analysis pipeline on Ubuntu.

## Quick Start (on the target Ubuntu machine)

```bash
# Install Ansible if not already present
sudo apt install ansible

# Clone the repo
git clone https://github.com/c03rad0r/proteomics-pipeline.git
cd proteomics-pipeline

# Run the playbook (installs everything — takes ~20-40 min)
ansible-playbook ansible/setup.yml
```

## What It Installs

**System packages (apt):**
- `r-base`, `r-base-dev` — R runtime + dev headers
- `libxml2-dev`, `libcurl4-openssl-dev`, `libssl-dev` — R package build deps
- `libcairo2-dev`, `libxt-dev` — plotting backends
- `libgdal-dev`, `libgeos-dev`, `libproj-dev`, `libudunits2-dev` — spatial/GDAL
- `texlive-latex-base`, `texlive-latex-extra`, `texlive-fonts-recommended` — LaTeX for PDF output
- `pandoc`, `qpdf` — document conversion
- `make`, `git`, `curl`, `wget`

**R packages (CRAN):**
- `ggplot2`, `ggrepel`, `scales`, `RColorBrewer`
- `ggVennDiagram`, `VennDiagram`, `UpSetR`
- `dplyr`, `readr`, `tibble`, `tidyr`, `stringr`
- `patchwork`, `igraph`
- `gprofiler2`, `config`, `here`

**R packages (Bioconductor):**
- `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `DOSE`
- `EnhancedVolcano`, `STRINGdb`, `rrvgo`, `ComplexHeatmap`

## After Installation

```bash
# Copy Aruna's data files into data/
cp /path/to/aruna-data/*.csv data/

# Generate plots
make targeted-go     # GO dotplots (the ones with wrap variants)
make volcano          # volcano plots
make all              # full pipeline
```

## Running on Remote Machines

Create `ansible/inventory.ini`:

```ini
[ubuntu]
192.168.1.100 ansible_user=youruser

[ubuntu:vars]
ansible_python_interpreter=/usr/bin/python3
```

Then:

```bash
ansible-playbook ansible/setup.yml -i ansible/inventory.ini
```

## Supported Ubuntu Versions

Tested on Ubuntu 26.04 LTS (Resolute Raccoon).
Should work on Ubuntu 22.04+ (R 4.1+ from apt).