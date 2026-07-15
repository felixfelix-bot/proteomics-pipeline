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

### Linux — Fast Install (RECOMMENDED)

Three options, fastest first:

#### Option A: Conda — precompiled binaries, no sudo, no compilation (~5 min)

```bash
# Prerequisite: miniconda installed (https://docs.conda.io/en/latest/miniconda.html)
# One command — installs R + ALL packages as precompiled binaries
ansible-playbook ansible/setup.yml -e install_method=conda

# Activate the conda env
conda activate r-proteomics

# Copy CSV data files into data/, then run:
make test          # test pipeline on synthetic data
make all           # run full pipeline on real data
make targeted-go   # generate GO dotplots only
```

#### Option B: apt prebuilt + minimal Bioconductor (~5-10 min, needs sudo)

```bash
# Uses Ubuntu's prebuilt R packages for ~90% of deps.
# Only 4-5 small Bioconductor packages compile from source.
ansible-playbook ansible/setup-prebuilt.yml --connection=local --ask-become-pass

# Copy CSV data files into data/, then run:
make test          # test pipeline on synthetic data
make all           # run full pipeline on real data
make targeted-go   # generate GO dotplots only
```

#### Option C: apt + compile everything from source (~30+ min, needs sudo)

```bash
# Slowest — compiles ALL R packages from source. Only use if A and B fail.
ansible-playbook ansible/setup.yml -e install_method=sudo --ask-become-pass

# Copy CSV data files into data/, then run:
make test          # test pipeline on synthetic data
make all           # run full pipeline on real data
make targeted-go   # generate GO dotplots only
```

### Verifying the install

After any install method, verify all packages load:

```bash
Rscript -e 'pkgs <- c("ggplot2","clusterProfiler","enrichplot","STRINGdb","ComplexHeatmap","EnhancedVolcano","rrvgo","DOSE","org.Hs.eg.db","igraph","ggVennDiagram"); for (p in pkgs) cat(sprintf("  %s  %s %s\n", ifelse(requireNamespace(p,quietly=TRUE),"[OK]","[FAIL]"), p, ifelse(requireNamespace(p,quietly=TRUE),as.character(packageVersion(p)),"MISSING")))'
```

All should show [OK]. Then test on synthetic data:

```bash
make test    # generates synthetic data, runs full pipeline, outputs to output/figures/
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
