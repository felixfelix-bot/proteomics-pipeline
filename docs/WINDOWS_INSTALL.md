# Windows Install

## One command

Clone the repo, then run the bootstrap script:

```powershell
git clone https://github.com/c03rad0r/proteomics-pipeline.git
cd proteomics-pipeline
.\ansible\bootstrap-windows-all.ps1
```

The script auto-elevates to Admin, installs R + Rtools + Git, and installs all R packages. Takes ~30 min. When it says DONE, you're ready.

## Then run the pipeline

```powershell
# Copy your CSV files into the data\ folder first, then:
Rscript.exe run_all.R
```

To test with fake data first: `Rscript.exe run_all.R --test`

## Data files

Copy your `_diffEx_minProb.csv` files into `data\`:

- `BK467_TRIP4_vs_BK467_WT.csv`
- `BK516_Cflag_vs_BK516_Ctrl.csv`
- etc.

`known_interactors.txt` is already included.

## Troubleshooting

- **"Rscript not found"** → Restart PowerShell, or use full path: `& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"`
- **"Bioconductor version requires R version"** → Run: `Rscript.exe -e "BiocManager::install(version='3.22',ask=FALSE,update=TRUE,checkBuilt=TRUE)"` then re-run the package installer
- **Package compile fails** → Ensure Rtools is installed (the bootstrap does this, but if it failed, download from cran.r-project.org/bin/windows/Rtools/)
