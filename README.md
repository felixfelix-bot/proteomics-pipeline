# Proteomics Analysis Pipeline — TRIP4/ASCC Study

R-based proteomics data analysis pipeline for TRIP4/ASCC mass spectrometry data.
Implements Lydia's analytical approach: multi-category volcano plots, STRING network 
mapping, gene family highlighting, cross-experiment overlap, and GO enrichment.

## Study Overview

**Protein of interest:** TRIP4 (with ASCC complex interactors)  
**Organism:** Human (*Homo sapiens*)  
**Experimental systems:**
- **TurboID** (HeLa): proximity labeling — BK467 (WT/TRIP4/±RA), BK504 (TRIP4±RA)
- **Flag IP** (HEK293): co-immunoprecipitation — BK516 (C-Flag/N-Flag/Ctrl), BK523 (+RA)
- **CRAC** (HEK293): UV cross-linking RNA interactome (if data available)
- **RA** = Retinoic acid treatment

## Quick Start

### Windows — one command

After cloning this repo, open PowerShell and run:

```powershell
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
.\ansible\bootstrap-windows-all.ps1
```

This installs R, Rtools, Git, and all R packages automatically (~30 min). When it says DONE, copy your CSV files to `data\` and run:

```powershell
Rscript.exe run_all.R
```

For full instructions and troubleshooting see **[Windows Install Guide](docs/WINDOWS_INSTALL.md)**.

### Linux — Fast Install

**Prerequisites:** Ubuntu 22.04/24.04/26.04, Ansible 2.10+

#### Option A: Conda — no sudo needed, all precompiled binaries (~5-10 min)

```bash
# 1. Install Ansible (if not already installed)
sudo apt-get update && sudo apt-get install -y ansible

# 2. Clone the repo
git clone https://github.com/c03rad0r/proteomics-pipeline.git
cd proteomics-pipeline

# 3. Run the conda playbook — auto-installs miniconda if not found
#    ALL R packages installed as precompiled binaries, zero compilation
ansible-playbook ansible/setup.yml

# 4. Activate the conda env (run this in each new terminal)
conda activate r-proteomics

# 5. Copy your CSV data files into data/
cp /path/to/your/data/*.csv data/

# 6. Run the pipeline
make test          # test on synthetic data (generates fake data, runs all steps)
make all           # run full pipeline on real data
make targeted-go   # GO dotplots only
```

#### Option B: apt prebuilt + minimal Bioconductor compile (~5-10 min, needs sudo)

```bash
# Same prerequisites as above (Ansible + git clone)
# Uses apt prebuilt R packages for ~90% of deps.
# Only 4-5 small Bioconductor packages compile from source
# (clusterProfiler, enrichplot, rrvgo, EnhancedVolcano, config).
ansible-playbook ansible/setup-prebuilt.yml --connection=local --ask-become-pass

# Verify
make check
# Expected: "17/17 required packages installed — All packages ready!"

# Copy data and run
cp /path/to/your/data/*.csv data/
make test          # test on synthetic data
make all           # run full pipeline on real data
make targeted-go   # GO dotplots only
```

#### Option C: Full source compile (fallback only, ~30+ min)

```bash
# Slowest — compiles ALL R packages from source. Only use if A and B fail.
ansible-playbook ansible/setup.yml -e install_method=sudo --ask-become-pass
```

**If `make check` fails** (packages missing), install them manually:

```bash
# Check which packages are missing
Rscript R/check_packages.R

# Install missing CRAN packages
Rscript -e 'install.packages("PACKAGE_NAME", repos="https://cloud.r-project.org")'

# Install missing Bioconductor packages
Rscript -e 'BiocManager::install("PACKAGE_NAME", update=FALSE, ask=FALSE)'

# If you get 00LOCK errors, clean them up:
find ~/R/x86_64-pc-linux-gnu-library/ -name '00LOCK*' -exec rm -rf {} + 2>/dev/null
# Then retry the install
```

## Analysis Modules

| Script | Description |
|--------|-------------|
| `02_volcano_plots.R` | Lydia-style multi-category volcano plots (interactors, gene families, Flag IP overlay) |
| `03_venn_diagrams.R` | Venn diagrams: TurboID vs Flag IP, RA-specific changes |
| `04_go_enrichment.R` | GO enrichment (clusterProfiler ORA + GSEA + rrvgo) |
| `05_string_network.R` | STRING network mapping, candidate interactor identification |
| `06_gene_families.R` | GPATCH/DHX/DDX/LARP family highlighting on volcanos |
| `07_overlap_analysis.R` | Cross-experiment overlap (Flag IP on TurboID, CRAC, RA effects) |
| `08_targeted_volcanos.R` | Custom volcanos: TRIP4 vs WT (ASCC+interactors labeled), RA effect (orange/blue) |
| `09_flagip_volcano.R` | Flag IP validation volcano (C-Flag/N-Flag/both categories) |
| `13_chx_crac_analysis.R` | CHX/DMSO + CRAC analysis |
| `16_string_network_targeted.R` | Lydia-style direct PPI STRING network (seeds + neighbors) |
| `17_crac_string_network.R` | CRAC RNA interactome STRING network |
| `19_gsea.R` | GSEA enrichment (ranked gene list, all experiments) |
| `20_pathway_network.R` | STRING network maps for enriched GO pathways |
| `21_lydia_network_volcano.R` | Lydia volcano w/ STRING physical network overlay + gene families |
| `22_ra_common_analysis.R` | Common RA-enriched/depleted across RA02+RA04 concentrations |
| `23_bidirectional_go.R` | Bidirectional GO dot plot (up=right, down=left, single figure) |
| `24_chx_common_analysis.R` | CHX-enriched/depleted sets + STRING networks + GO |
| `25_shinygo_comparison.R` | Compare ShinyGO exports with our STRING pipeline results |

## Key Make Targets

| Command | What it produces |
|---------|-----------------|
| `make all` | Full pipeline |
| `make targeted-volcano` | TRIP4 vs WT + RA effect volcanos |
| `make flagip-volcano` | Flag IP validation volcano |
| `make lydia-volcano` | STRING network overlay volcano (Lydia's method) |
| `make ra-common` | Common RA proteins across both concentrations + networks + Venn |
| `make bidirectional-go` | Up-right/down-left GO dot plot |
| `make chx-common` | CHX-enriched/depleted + STRING + GO |
| `make crac-network` | CRAC RNA interactome STRING network |
| `make gsea` | GSEA enrichment on all experiments |
| `make string-network` | Lydia-style direct PPI STRING network |
| `make shinygo-compare` | Compare ShinyGO export with our STRING results |
| `make all-volcano` | All volcano plots in one go |

**Full architecture with flowcharts, timing diagrams, and verification checklist:**
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)

## Data Format

CSV files with these columns (from DIA-NN / LFQ pipeline output):

| Column | Description |
|--------|-------------|
| `Gene` | Gene symbol (e.g., TRIP4, ASCC1) |
| `Entry.Name` | UniProt entry name |
| `UniProt.ID` | UniProt accession |
| `logFC` | Log fold change |
| `P.Val` | Raw p-value |
| `adj.P.Val` | Adjusted p-value (FDR) |
| `is_significant` | Pre-computed significance flag |

Only `Gene`, `logFC`, and `adj.P.Val` are required by the pipeline.

## Key Comparisons

- **TurboID main:** BK467_TRIP4_vs_BK467_WT
- **Flag IP main:** BK516_Cflag_vs_BK516_Ctrl
- **RA effect TurboID:** BK467_TRIP4_RA02_vs_BK467_TRIP4
- **RA effect Flag IP:** BK523_Cflag_RA04_vs_BK516_Cflag

## Significance Thresholds

- adj.P.Val < 0.05
- |logFC| > 0.5 (matching Lydia's scripts)

## Code/Data Separation

This repo contains **only code**. Actual proteomics data is never committed.
The `.gitignore` excludes all CSV/data files except README.md.
