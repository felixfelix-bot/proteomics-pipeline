# Windows Installation Guide

Complete setup instructions for running the TRIP4/ASCC proteomics analysis pipeline
on a Windows laptop. Takes 30-45 minutes total.

---

## Option A: One-Click Script (Recommended)

Does everything automatically: installs R, Rtools, Git, clones the repo, and installs all R packages.

### Step 1: Open PowerShell as Administrator

- Click Start, type "PowerShell"
- Right-click "Windows PowerShell" → "Run as Administrator"

### Step 2: Allow scripts

```powershell
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
```

Type `Y` when prompted.

### Step 3: Download and run the bootstrap script

```powershell
curl -L -o bootstrap.ps1 https://raw.githubusercontent.com/c03rad0r/proteomics-pipeline/main/ansible/bootstrap-windows-all.ps1
.\bootstrap.ps1
```

The script will:
1. Install Chocolatey (Windows package manager)
2. Install Git
3. Install R (statistical computing language)
4. Install Rtools (C++ compiler needed for Bioconductor packages)
5. Clone the proteomics-pipeline repository to `C:\proteomics-pipeline`
6. Install all required R packages (this takes 20-40 minutes)

**Go get coffee.** Step 6 is a long compilation step and there's no way to speed it up.

When it finishes, you'll see "Setup Complete!" in green.

### Step 4: Test with synthetic data

```powershell
cd C:\proteomics-pipeline
Rscript.exe run_all.R --test
```

This generates fake test data and runs the entire pipeline. If it finishes without
errors, everything is installed correctly. You'll see output figures in
`output\figures\` and tables in `output\tables\`.

### Step 5: Run on real data

1. Copy your CSV files to `C:\proteomics-pipeline\data\`
   - Files should be named like: `BK467_TRIP4_vs_BK467_WT.csv`
   - The `known_interactors.txt` file is already in `data\`
2. Run the pipeline:
```powershell
cd C:\proteomics-pipeline
Rscript.exe run_all.R
```

---

## Option B: Manual Install (if the script doesn't work)

### Step 1: Install R

1. Go to https://cran.r-project.org/bin/windows/base/
2. Download the latest R installer (R-4.x.x for Windows)
3. Run the installer, click Next through all defaults
4. Note the install path (usually `C:\Program Files\R\R-4.x.x\`)

### Step 2: Install Rtools

**This is required.** Without Rtools, Bioconductor packages cannot be compiled.

1. Go to https://cran.r-project.org/bin/windows/Rtools/
2. Download the Rtools version matching your R version:
   - R 4.4 → Rtools43
   - R 4.5 → Rtools44
3. Run the installer, click Next through all defaults

### Step 3: Install Git

1. Go to https://git-scm.com/
2. Download and install with default settings

### Step 4: Clone the repository

Open PowerShell or Command Prompt:
```powershell
cd C:\
git clone https://github.com/c03rad0r/proteomics-pipeline.git
cd proteomics-pipeline
```

### Step 5: Install R packages

```powershell
Rscript.exe R\00_install_packages.R
```

This takes 20-40 minutes. It installs:
- **CRAN packages:** ggplot2, ggrepel, ggVennDiagram, UpSetR, gprofiler2, etc.
- **Bioconductor packages:** clusterProfiler, org.Hs.eg.db, EnhancedVolcano,
  enrichplot, STRINGdb, rrvgo, ComplexHeatmap

You'll see `[OK]` next to packages that are already installed and `[INSTALL]`
for ones being compiled. If some fail, re-run the script — it only installs
what's missing.

If you get errors about "Bioconductor version requires R version":
```powershell
Rscript.exe -e "BiocManager::install(version = '3.22', ask=FALSE); BiocManager::install(c('clusterProfiler','org.Hs.eg.db','enrichplot','DOSE','rrvgo'), ask=FALSE)"
```

### Step 6: Verify installation

```powershell
Rscript.exe -e "library(clusterProfiler); library(org.Hs.eg.db); library(EnhancedVolcano); cat('All packages loaded OK\n')"
```

If you see "All packages loaded OK" — you're ready.

### Step 7: Test and run

Same as Steps 4 and 5 from Option A above.

---

## Troubleshooting

### "Rscript is not recognized as a command"

R is not in your PATH. Use the full path instead:
```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" R\00_install_packages.R
```
(Adjust the version number to match what you installed.)

### "Bioconductor version requires R version"

Your R and Bioconductor versions are mismatched. Run:
```powershell
Rscript.exe -e "BiocManager::install(version = '3.22', ask=FALSE, update=TRUE, checkBuilt=TRUE)"
```
Then re-run `R\00_install_packages.R`.

### Package compilation fails with "could not find tools"

Rtools is not installed or not in PATH. Reinstall Rtools (Step 2 above) and
restart your PowerShell window.

### "ERROR: dependency 'xxx' is not available"

Some packages need other packages as prerequisites. Re-run the installer:
```powershell
Rscript.exe R\00_install_packages.R
```
The script is idempotent — it installs what's missing and skips what's present.

### The pipeline runs but produces no figures

Check that your CSV files are in the `data\` folder and have the correct column
names: `Gene`, `logFC`, `adj.P.Val`. See `data\README.md` for details.

---

## Data File Naming

Place CSV files in `C:\proteomics-pipeline\data\` with these names:

| File | Experiment |
|------|-----------|
| `BK467_TRIP4_vs_BK467_WT.csv` | TurboID: TRIP4 vs Wild Type |
| `BK467_TRIP4_RA02_vs_BK467_TRIP4.csv` | TurboID: Effect of retinoic acid on TRIP4 |
| `BK467_TRIP4_RA02_vs_BK467_WT.csv` | TurboID: TRIP4+RA vs Wild Type |
| `BK516_Cflag_vs_BK516_Ctrl.csv` | Flag IP: C-Flag vs Control |
| `BK516_Nflag_vs_BK516_Ctrl.csv` | Flag IP: N-Flag vs Control |
| `BK523_Cflag_RA04_vs_BK516_Cflag.csv` | Flag IP: RA effect on C-Flag interactome |

Files use only the `_diffEx_minProb.csv` output from the LFQ pipeline.
Do NOT use `imputatedMatrix` or `originalMatrix` files.

---

## Getting Help

If something goes wrong:
1. Copy the full error message
2. Note which step failed
3. Send the error text to the chat

We can diagnose from the error output.
