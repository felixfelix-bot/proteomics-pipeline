# Windows Install

## One command

Open PowerShell and paste this:

```
irm https://raw.githubusercontent.com/c03rad0r/proteomics-pipeline/main/ansible/bootstrap-windows-all.ps1 | iex
```

It auto-elevates to Admin, installs R + Rtools + Git, clones the repo to `C:\proteomics-pipeline`, and installs all packages. Takes ~30 min. When it says DONE, you're ready.

## Then run the pipeline

```powershell
cd C:\proteomics-pipeline

# Copy your CSV files into the data\ folder first, then:
Rscript.exe run_all.R
```

To test with fake data first: `Rscript.exe run_all.R --test`

## Data files

Copy your `_diffEx_minProb.csv` files into `C:\proteomics-pipeline\data\`:

- `BK467_TRIP4_vs_BK467_WT.csv`
- `BK516_Cflag_vs_BK516_Ctrl.csv`
- etc.

`known_interactors.txt` is already included.

## Troubleshooting

- **"Rscript not found"** → Restart PowerShell, or use full path: `& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"`
- **"Bioconductor version requires R version"** → Run: `Rscript.exe -e "BiocManager::install(version='3.22',ask=FALSE,update=TRUE,checkBuilt=TRUE)"` then re-run the package installer
- **Package compile fails** → Ensure Rtools is installed (the bootstrap does this, but if it failed, download from cran.r-project.org/bin/windows/Rtools/)
