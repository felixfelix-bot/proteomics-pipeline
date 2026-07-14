# OneDrive Recovery — Handover Document

**Purpose:** Resolve git repository corruption caused by OneDrive uninstall on
Dr. Aruna's Windows laptop. This document contains all context needed by a
third-party helper who has NOT seen the conversation history.

---

## The Problem

Dr. Aruna's proteomics analysis pipeline (`proteomics-pipeline`) was cloned
inside an OneDrive-synced folder:

```
C:\Users\aruna\OneDrive\Dokumente\GitHub\proteomics-pipeline
```

OneDrive's "Files On-Demand" feature marks files as "online-only" — placeholder
stubs that aren't actually on disk. When OneDrive was uninstalled (due to disk
space issues), these stubs became orphaned. Files may appear to exist in
Explorer but contain zero bytes or unreadable data.

This affects TWO things:
1. **Git repository integrity** — `.git/objects/` files may be corrupt stubs
2. **Data files** — the experimental CSV files (NOT on GitHub, gitignored)
   may be corrupt stubs

The code itself is safe — it's on GitHub. The risk is local data loss and
broken git operations.

---

## What's at Risk

### Data Files (ONLY on Aruna's laptop — NOT on GitHub)

These files are gitignored and exist ONLY locally. If they're corrupt, they
must be recovered from onedrive.live.com before that expires.

| File/Pattern | Location | What It Is |
|---|---|---|
| `*_diffEx_minProb.csv` | `data/subfolders/` | Mass spectrometry output (~15 files, ~1.4MB each). TurboID, Flag IP, CHX/DMSO experiments. |
| `FLAG-TRIP4_list_CRACdata.csv` | `data/` | CRAC RNA interactome data (different column names: `external_gene_name`, `logFC`, `FDR`) |
| `9606.protein.physical.links.v12.0.txt` | `data/` | STRING database physical interactions (~65MB, 1.5M edges). Can be re-downloaded from STRING if lost. |
| `known_interactors.txt` | `data/` | List of ~76 known TRIP4 interactors. **THIS IS ON GITHUB** (gitignored exception) — safe, can be recovered via git. |
| `shinygo_export.csv` | `data/` (optional) | ShinyGO comparison export. Only exists if Aruna ran `make shinygo-compare`. Not critical. |

**CRITICAL:** The `*_diffEx_minProb.csv` files and `FLAG-TRIP4_list_CRACdata.csv`
are the experimental data. They CANNOT be regenerated. They must be recovered
from OneDrive cloud (onedrive.live.com) or from Aruna's backup if she has one.

### Git Repository

The code is fully backed up on GitHub:
```
https://github.com/c03rad0r/proteomics-pipeline.git
```
Repository is PRIVATE. Access requires c08r4d0r's GitHub credentials.

If git operations fail (log, pull, push), the simplest fix is to delete the
local repo and clone fresh. No code will be lost.

---

## Diagnostic Steps (Do These First)

### Step 1: Test Git Integrity

Open PowerShell in the repo folder and run:

```powershell
cd C:\Users\aruna\OneDrive\Dokumente\GitHub\proteomics-pipeline

# Does git work at all?
git log --oneline -5

# Full integrity check
git fsck --full
```

**If both succeed with no errors:** Git is fine. The repo survived the
OneDrive uninstall. Skip to "Relocate the Repo" below.

**If git log fails or fsck reports errors:** Git objects are corrupted.
Proceed to "Recovery Procedure" below.

### Step 2: Check Data File Integrity

Check if the data CSV files are real or zero-byte stubs:

```powershell
# List all CSV files with their sizes
Get-ChildItem -Path data\ -Filter *.csv -Recurse | Select-Object FullName, Length

# Check the CRAC file specifically
Get-Item data\FLAG-TRIP4_list_CRACdata.csv | Select-Object Length

# Check the STRING file
Get-Item data\9606.protein.physical.links.v12.0.txt -ErrorAction SilentlyContinue | Select-Object Length
```

**Red flags:**
- Any CSV file showing `Length = 0` → corrupted stub
- Any CSV file showing `Length < 1000` (less than 1KB) → likely a stub (real files are ~1.4MB)
- Missing files entirely → were online-only, never downloaded

**If files are intact:** Great. Back them up before doing anything else.

### Step 3: Back Up Data NOW

Before any further operations, copy all data files to a safe location:

```powershell
# Create backup on Desktop
$backup = "$env:USERPROFILE\Desktop\proteomics_data_backup_$(Get-Date -Format 'yyyyMMdd')"

# Copy ALL data files (preserves folder structure)
Copy-Item -Path data\* -Destination $backup -Recurse

Write-Host "Backup created at: $backup"
Write-Host "Verify the file sizes match before proceeding."
```

---

## Recovery Procedure

### Scenario A: Git Is Corrupted, Data Is Intact

1. **Back up data** (Step 3 above)
2. **Delete the old repo:**
   ```powershell
   cd C:\Users\aruna
   Remove-Item -Recurse -Force "OneDrive\Dokumente\GitHub\proteomics-pipeline"
   ```
3. **Clone fresh to a clean path** (NOT inside OneDrive):
   ```powershell
   mkdir C:\Users\aruna\GitHub
   cd C:\Users\aruna\GitHub
   git clone https://github.com/c03rad0r/proteomics-pipeline.git
   ```
4. **Restore data files:**
   ```powershell
   Copy-Item -Path "$env:USERPROFILE\Desktop\proteomics_data_backup_*\*" -Destination C:\Users\aruna\GitHub\proteomics-pipeline\data\ -Recurse
   ```
5. **Verify:**
   ```powershell
   cd C:\Users\aruna\GitHub\proteomics-pipeline
   git log --oneline -5
   .\run.ps1 diagnostics
   ```

### Scenario B: Data Files Are Corrupted / Missing

1. **Try OneDrive cloud recovery:** Go to onedrive.live.com, sign in with
   Aruna's Microsoft account, navigate to Documents → GitHub →
   proteomics-pipeline → data, and download the CSV files.

2. **If OneDrive cloud is unavailable** (account closed, files purged):
   - The `*_diffEx_minProb.csv` files must be re-exported from the original
     mass spectrometry analysis (DIA-NN + limma pipeline). Aruna or her
     collaborator Benno would have the raw data.
   - `FLAG-TRIP4_list_CRACdata.csv` must be re-exported from the CRAC
     analysis pipeline.
   - `9606.protein.physical.links.v12.0.txt` can be re-downloaded from
     https://string-db.org/cgi/download (select Homo sapiens, physical links,
     v12.0)
   - `known_interactors.txt` is on GitHub — recover with:
     ```powershell
     git show HEAD:data/known_interactors.txt > data\known_interactors.txt
     ```

3. **Once data is recovered**, follow Scenario A steps 3-5 to clone fresh and
   restore.

### Scenario C: Everything Is Fine (Git + Data Intact)

Still relocate to avoid future issues:

1. Back up data (Step 3 above)
2. Move the entire folder:
   ```powershell
   Move-Item "C:\Users\aruna\OneDrive\Dokumente\GitHub\proteomics-pipeline" "C:\Users\aruna\GitHub\proteomics-pipeline"
   ```
3. Or clone fresh (cleaner — avoids any residual OneDrive metadata):
   Follow Scenario A steps 2-5.

---

## New Directory Structure (Target)

```
C:\Users\aruna\GitHub\proteomics-pipeline\
├── .git\                          ← healthy, no OneDrive stubs
├── R\                             ← analysis scripts (from GitHub)
├── data\                          ← experimental CSVs (restored from backup)
│   ├── <experiment_folders>\
│   │   └── *_diffEx_minProb.csv   ← mass spec data (~15 files)
│   ├── FLAG-TRIP4_list_CRACdata.csv  ← CRAC RNA data
│   ├── 9606.protein.physical.links.v12.0.txt  ← STRING database (65MB)
│   └── known_interactors.txt      ← on GitHub, always recoverable
├── output\                        ← generated by pipeline (gitignored)
├── Makefile
├── run.ps1                        ← PowerShell runner (replaces make.exe)
├── docs\
│   ├── HANDOVER.md                ← full maintenance guide
│   └── HANDOVER_MAKE_TARGETS.md   ← Makefile patterns for AI/devs
└── ...
```

---

## Additional Context: make.exe Is Also Blocked

Aruna's Windows machine has Application Control policy that blocks
`make.exe`. The repo includes `run.ps1` — a PowerShell script that replaces
make entirely. After recovery, Aruna runs:

```powershell
cd C:\Users\aruna\GitHub\proteomics-pipeline
.\run.ps1 aruna-fast       # instead of: make aruna-fast
.\run.ps1 aruna-all        # instead of: make aruna-all
.\run.ps1 diagnostics      # verify data is loaded correctly
.\run.ps1 help             # see all targets
```

Do NOT install WSL. The PowerShell runner handles everything.

---

## R Installation

R must be installed and on PATH. Check:

```powershell
Rscript --version
```

If not found, R is likely at:
```
C:\Program Files\R\R-4.6.1\bin\Rscript.exe
```

The `run.ps1` script auto-detects the R path. If it can't find R, hardcode
the path at the top of `run.ps1` (look for the `$RSCRIPT` variable).

---

## Verification Checklist

After recovery, verify each item:

```
[ ] 1. git log --oneline -5 works (shows recent commits)
[ ] 2. git fsck --full reports no errors
[ ] 3. data/ folder contains *_diffEx_minProb.csv files (check sizes > 100KB)
[ ] 4. data/FLAG-TRIP4_list_CRACdata.csv exists and is > 0 bytes
[ ] 5. data/known_interactors.txt exists (28 lines, ~76 genes)
[ ] 6. .\run.ps1 diagnostics runs and prints experiment counts
[ ] 7. .\run.ps1 aruna-fast completes and generates figures in output/figures/
[ ] 8. Repo lives at C:\Users\aruna\GitHub\ (NOT OneDrive path)
```

---

## Key People

| Role | Who | What They Do |
|------|-----|-------------|
| Biochemist | Dr. Aruna | Runs the pipeline on her Windows laptop. Has the experimental data. Not a programmer. |
| Coordinator | c08r4d0r | Coordinates between groups. GitHub repo owner. Technical but not an R developer. |
| Third Group | (this handover's recipient) | Resolving the OneDrive/git corruption issue. |

---

## What NOT to Do

1. **Do NOT install WSL** — the run.ps1 PowerShell script replaces make.exe
2. **Do NOT re-clone without backing up data first** — the CSV files are NOT on GitHub
3. **Do NOT put the new repo inside any cloud-sync folder** (OneDrive, Google Drive, Dropbox, iCloud)
4. **Do NOT push data files to GitHub** — they're gitignored for a reason (sensitive research data)
5. **Do NOT delete data\known_interactors.txt** — it's the only data file tracked by git
6. **Do NOT run make.exe** — it's blocked by Application Control policy. Use run.ps1
7. **Do NOT push to ngit (nostr git)** — this is a private repo. GitHub only.

---

## Technical Details for Helpers

### Why This Pipeline Exists

This is an R-based proteomics analysis pipeline for TRIP4/ASCC protein
interaction research. It takes mass spectrometry CSV data and produces
publication-quality figures: volcano plots, Venn diagrams, GO enrichment dot
plots, STRING protein interaction networks.

### Key Architecture Points

- **Data/code separation:** Code is on GitHub. Data stays local (gitignored).
  Aruna runs `git pull` then `.\run.ps1 aruna-fast` to generate figures.
- **run_step.R dispatcher:** All make/run.ps1 targets call
  `Rscript R/run_step.R <step_name>`, which loads config + utils, then runs
  the target script. See `docs/HANDOVER_MAKE_TARGETS.md` for full details.
- **33 R scripts** in `R/` directory, each producing specific figure types.
- **Full maintenance guide** in `docs/HANDOVER.md`.

### If R Package Reinstallation Is Needed

If the fresh clone needs packages installed:
```powershell
.\run.ps1 install
```
This runs `R/00_install_packages.R` which installs ~20 packages including
Bioconductor packages (EnhancedVolcano, clusterProfiler, STRINGdb). Takes
30+ minutes on first install.
