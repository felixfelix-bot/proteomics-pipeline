# Proteomics Analysis Pipeline

R-based proteomics data analysis pipeline for mass spectrometry data.
Generates volcano plots, Venn diagrams, and Gene Ontology (GO) enrichment analysis.

## Target

- **Organism:** Human (*Homo sapiens*)
- **Input:** CSV files with gene symbols, log2 fold change, p-adjusted values
- **Output:** Publication-quality figures + enrichment result tables

## Quick Start

### Option 1: Full Ansible setup (recommended)

**On Linux control machine (e.g., the AI laptop):**

```bash
# Install Ansible
chmod +x ansible/bootstrap-linux.sh
./ansible/bootstrap-linux.sh

# Edit inventory to point to your Windows machine
vim ansible/inventory.yml

# Run the playbook
cd ansible
ansible-playbook -i inventory.yml setup.yml
```

**On Windows target machine (one-time WinRM setup):**

Run PowerShell as Administrator:
```powershell
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
.\ansible\bootstrap-windows.ps1
```

### Option 2: Local R install only

If you just want the R environment without Ansible:

```bash
# Linux
Rscript R/00_install_packages.R

# Windows (from R or Rscript)
Rscript.exe R\00_install_packages.R
```

## Project Structure

```
proteomics-pipeline/
├── ansible/                  # Infrastructure as code
│   ├── bootstrap-linux.sh    # Installs Ansible on Linux control machine
│   ├── bootstrap-windows.ps1 # Configures WinRM on Windows target
│   ├── inventory.yml         # Target machine configuration
│   ├── setup.yml             # Main Ansible playbook
│   └── roles/
│       ├── r-base/           # R language installation
│       ├── r-packages/       # Bioconductor + CRAN R packages
│       └── vscode/           # VS Code + R extensions
├── R/                        # Analysis scripts
│   ├── 00_install_packages.R # Package installer
│   ├── 01_config.R           # Configuration (thresholds, paths)
│   ├── 02_volcano_plots.R    # Volcano plot generation
│   ├── 03_venn_diagrams.R    # Venn diagram + set extraction
│   ├── 04_go_enrichment.R    # GO enrichment analysis
│   ├── utils.R               # Shared utility functions
│   └── generate_synthetic_data.R  # Test data for CI
├── data/                     # CSV data files (gitignored)
├── tests/                    # Test suite
└── README.md
```

## Data Requirements

Place CSV files in `data/`. Expected columns:

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol (e.g., TP53, BRCA1) |
| `log2FC` | Log2 fold change |
| `padj` | Adjusted p-value (FDR) |

## Analysis Modules

### 1. Volcano Plots
- EnhancedVolcano for single-experiment publication plots
- ggplot2 overlay for TurboID vs Flag co-IP comparison
- Labels known interactors

### 2. Venn Diagrams
- ggVennDiagram for 2-4 set comparisons
- UpSetR for >4 sets
- Programmatic overlap extraction for downstream analysis

### 3. GO Enrichment
- clusterProfiler + enrichGO (ORA)
- clusterProfiler + gseGO (GSEA)
- rrvgo for GO term redundancy reduction
- gprofiler2 as cross-check
- org.Hs.eg.db for human annotation

## Code/Data Separation

This repo contains **only code**. Actual proteomics data is never committed.
Researchers run the analysis scripts locally with their own data files.
