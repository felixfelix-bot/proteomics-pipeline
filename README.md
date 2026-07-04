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

### On Windows (researcher's laptop)

```powershell
# 1. Install R from cran.r-project.org
# 2. Clone the repo
git clone https://github.com/c03rad0r/proteomics-pipeline.git
cd proteomics-pipeline

# 3. Install R packages (15-30 min first time)
Rscript.exe R\00_install_packages.R

# 4. Place CSV files in data/ folder
#    (files named like BK467_TRIP4_vs_BK467_WT.csv etc.)

# 5. Run the pipeline
Rscript.exe run_all.R

# Or test with synthetic data first:
Rscript.exe run_all.R --test
```

### On Linux (with micromamba)

```bash
micromamba create -y -n r-env -c conda-forge r-base r-essentials
micromamba run -n r-env Rscript R/00_install_packages.R
micromamba run -n r-env Rscript run_all.R
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
