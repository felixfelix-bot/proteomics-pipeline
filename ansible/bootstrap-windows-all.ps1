###############################################################################
# ONE-COMMAND WINDOWS SETUP FOR PROTEOMICS PIPELINE
# Resilient — safe to run multiple times.
#
# Usage:
#   Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
#   .\ansible\bootstrap-windows-all.ps1
###############################################################################

$ErrorActionPreference = "Continue"  # Don't die on non-critical errors

function Step($n, $msg) { Write-Host "`n[$n] $msg" -ForegroundColor Cyan }
function OK($msg)       { Write-Host "  [OK] $msg" -ForegroundColor Green }
function Info($msg)     { Write-Host "  $msg" -ForegroundColor Gray }
function Warn($msg)     { Write-Host "  [WARN] $msg" -ForegroundColor Yellow }
function Skip($msg)     { Write-Host "  [SKIP] $msg already present" -ForegroundColor DarkGray }

# --- Self-elevate to Administrator ---
if (-not ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
    Write-Host "Requesting Administrator privileges..." -ForegroundColor Yellow
    $scriptPath = $MyInvocation.MyCommand.Path
    if (-not $scriptPath) { $scriptPath = Join-Path $PSScriptRoot "bootstrap-windows-all.ps1" }
    Start-Process powershell -Verb RunAs -ArgumentList "-NoProfile -ExecutionPolicy Bypass -File `"$scriptPath`""
    exit
}

# --- Helper: refresh PATH from registry (fixes "choco not found" after install) ---
function Refresh-Path {
    $env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
}

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host " Proteomics Pipeline - Auto Setup" -ForegroundColor Cyan
Write-Host " (Safe to re-run)" -ForegroundColor DarkCyan
Write-Host "========================================" -ForegroundColor Cyan

# =====================================================================
# 1. Chocolatey
# =====================================================================
Step "1/5" "Checking Chocolatey..."
Refresh-Path
if (Get-Command choco -ErrorAction SilentlyContinue) {
    Skip "Chocolatey"
} else {
    Info "Installing Chocolatey..."
    Set-ExecutionPolicy Bypass -Scope Process -Force
    [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
    try {
        Invoke-Expression ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))
        Refresh-Path
        OK "Chocolatey installed"
    } catch {
        # Chocolatey files might exist from partial install — check anyway
        Refresh-Path
        if (Get-Command choco -ErrorAction SilentlyContinue) {
            OK "Chocolatey available"
        } else {
            Warn "Chocolatey install failed, but continuing (may already have R/Rtools)"
        }
    }
}

# =====================================================================
# 2. Git + R + Rtools (skip if already installed)
# =====================================================================
Step "2/5" "Installing Git, R, and Rtools..."
Refresh-Path

# Git
if (Get-Command git -ErrorAction SilentlyContinue) {
    Skip "Git"
} elseif (Get-Command choco -ErrorAction SilentlyContinue) {
    Info "Installing Git..."
    choco install git -y --no-progress 2>$null
    Refresh-Path
    OK "Git installed"
} else {
    Warn "Cannot install Git (no Chocolatey). Install manually if needed."
}

# R — check common locations before trying to install
$Rscript = $null
$rPaths = @("C:\Program Files\R", "${env:LOCALAPPDATA}\Programs\R")
foreach ($basePath in $rPaths) {
    if (Test-Path $basePath) {
        $rDir = Get-ChildItem $basePath -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
        if ($rDir) {
            $candidate = Join-Path $rDir.FullName "bin\Rscript.exe"
            if (Test-Path $candidate) {
                $Rscript = $candidate
                break
            }
        }
    }
}
# Also check PATH
if (-not $Rscript) {
    $rscriptCmd = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rscriptCmd) { $Rscript = $rscriptCmd.Source }
}

if ($Rscript) {
    $rVer = (Split-Path (Split-Path $Rscript) -Leaf)
    Skip "R ($rVer)"
    OK "Using: $Rscript"
} elseif (Get-Command choco -ErrorAction SilentlyContinue) {
    Info "Installing R..."
    choco install r.project -y --no-progress 2>$null
    Refresh-Path
    $rDir = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
    if ($rDir) {
        $Rscript = Join-Path $rDir.FullName "bin\Rscript.exe"
        OK "R installed: $($rDir.Name)"
    } else {
        Warn "R install may have failed. Check C:\Program Files\R\"
    }
} else {
    Warn "Cannot install R (no Chocolatey). Install manually from cran.r-project.org"
}

# Rtools — check before installing
$rtoolsFound = $false
foreach ($rt in @("C:\rtools45", "C:\rtools44", "C:\rtools43")) {
    if (Test-Path $rt) { $rtoolsFound = $true; break }
}
if ($rtoolsFound) {
    Skip "Rtools"
} elseif (Get-Command choco -ErrorAction SilentlyContinue) {
    Info "Installing Rtools..."
    choco install rtools -y --no-progress 2>$null
    Refresh-Path
    if (Test-Path "C:\rtools45" -or (Test-Path "C:\rtools44") -or (Test-Path "C:\rtools43")) {
        OK "Rtools installed"
    } else {
        Warn "Rtools install may have failed. Download from cran.r-project.org/bin/windows/Rtools/"
    }
} else {
    Warn "Cannot install Rtools (no Chocolatey)."
}

# =====================================================================
# 3. Locate the repo (already cloned OR clone fresh)
# =====================================================================
Step "3/5" "Locating repository..."

# If we're already inside the repo (running from cloned copy), use here
$RepoDir = $null
if (Test-Path "R\00_install_packages.R") {
    $RepoDir = (Get-Item .).FullName
    OK "Already in repo: $RepoDir"
} elseif (Test-Path "C:\proteomics-pipeline\R\00_install_packages.R") {
    $RepoDir = "C:\proteomics-pipeline"
    OK "Repo found: $RepoDir"
} else {
    # Try to clone (public repos work without auth; private needs credentials)
    Info "Cloning repository to C:\proteomics-pipeline..."
    try {
        git clone https://github.com/c03rad0r/proteomics-pipeline.git "C:\proteomics-pipeline" 2>&1 | Out-Null
        if (Test-Path "C:\proteomics-pipeline\R\00_install_packages.R") {
            $RepoDir = "C:\proteomics-pipeline"
            OK "Repository cloned"
        } else {
            Warn "Clone failed (private repo — run from inside your existing clone instead)"
            Warn "Example: cd to your GitHub folder, then run this script"
        }
    } catch {
        Warn "Clone failed. If you already have the repo, cd into it and re-run this script."
    }
}

if (-not $RepoDir) {
    # Last resort: search common locations
    $searchPaths = @(
        "$env:USERPROFILE\OneDrive\Dokumente\GitHub\proteomics-pipeline",
        "$env:USERPROFILE\Documents\GitHub\proteomics-pipeline",
        "$env:USERPROFILE\GitHub\proteomics-pipeline",
        "C:\proteomics-pipeline"
    )
    foreach ($sp in $searchPaths) {
        if (Test-Path "$sp\R\00_install_packages.R") {
            $RepoDir = $sp
            OK "Found repo: $RepoDir"
            break
        }
    }
}

if (-not $RepoDir) {
    Write-Host "`n  [ERROR] Could not find the repository." -ForegroundColor Red
    Write-Host "  Options:" -ForegroundColor Yellow
    Write-Host "    1. cd into your cloned repo folder and re-run this script" -ForegroundColor Yellow
    Write-Host "    2. Clone manually: git clone https://github.com/c03rad0r/proteomics-pipeline.git" -ForegroundColor Yellow
    Write-Host "`n  Skipping to package install only (if Rscript is available)..." -ForegroundColor Yellow
} else {
    # Pull latest
    Push-Location $RepoDir
    try { git pull 2>&1 | Out-Null; OK "Pulled latest" } catch { Warn "git pull skipped" }
    Pop-Location
}

# =====================================================================
# 4. Install R packages (the long step)
# =====================================================================
Step "4/5" "Installing R packages..."

if (-not $Rscript) {
    # Try to find it one more time
    Refresh-Path
    $Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
}
if (-not $Rscript -and $RepoDir) {
    # Search filesystem
    $found = Get-ChildItem "C:\Program Files\R" -Recurse -Filter "Rscript.exe" -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($found) { $Rscript = $found.FullName }
}

if ($Rscript -and (Test-Path $Rscript)) {
    Info "Using Rscript: $Rscript"
    if ($RepoDir) {
        Info "Running package installer (20-40 minutes)..."
        Push-Location $RepoDir
        & $Rscript "R\00_install_packages.R"
        Pop-Location
        OK "Package installer finished"
    } else {
        # No repo found but R exists — just install packages inline
        Info "Installing core packages inline..."
        & $Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('clusterProfiler','org.Hs.eg.db','EnhancedVolcano','enrichplot','DOSE','STRINGdb','rrvgo','ComplexHeatmap','ggplot2','ggrepel','ggVennDiagram','UpSetR','gprofiler2','igraph'), ask=FALSE, update=FALSE)"
        OK "Packages installed"
    }
} else {
    Write-Host "  [ERROR] Rscript.exe not found." -ForegroundColor Red
    Write-Host "  Install R manually from https://cran.r-project.org/bin/windows/base/" -ForegroundColor Yellow
}

# =====================================================================
# 5. MiKTeX (LaTeX) + SumatraPDF (lightweight PDF viewer)
# =====================================================================
Step "5/7" "Installing MiKTeX (LaTeX)..."
Refresh-Path

if (Get-Command pdflatex -ErrorAction SilentlyContinue) {
    Skip "MiKTeX/pdflatex"
} elseif (Get-Command choco -ErrorAction SilentlyContinue) {
    Info "Installing MiKTeX via Chocolatey..."
    choco install miktex -y --no-progress 2>$null
    Refresh-Path
    if (Get-Command pdflatex -ErrorAction SilentlyContinue) {
        OK "MiKTeX installed"
    } else {
        Warn "MiKTeX install may have failed. Download from https://miktex.org/download"
    }
} else {
    Warn "Cannot install MiKTeX (no Chocolatey). Download from https://miktex.org/download"
}

# ---------------------------------------------------------------------
# 5b. SumatraPDF — lightweight PDF viewer (Windows equivalent of evince)
#     Auto-reloads PDF on recompile, supports command-line open.
# ---------------------------------------------------------------------
Step "5b/7" "Installing SumatraPDF (PDF viewer)..."
$sumatra = Get-Command SumatraPDF -ErrorAction SilentlyContinue
if (-not $sumatra) {
    $sumatraExe = "${env:ProgramFiles}\SumatraPDF\SumatraPDF.exe"
    if (Test-Path $sumatraExe) {
        $sumatra = $true
    }
}
if ($sumatra) {
    Skip "SumatraPDF"
} elseif (Get-Command choco -ErrorAction SilentlyContinue) {
    Info "Installing SumatraPDF via Chocolatey..."
    choco install sumatrapdf -y --no-progress 2>$null
    Refresh-Path
    $sumatraExe = "${env:ProgramFiles}\SumatraPDF\SumatraPDF.exe"
    if ((Get-Command SumatraPDF -ErrorAction SilentlyContinue) -or (Test-Path $sumatraExe)) {
        OK "SumatraPDF installed"
    } else {
        Warn "SumatraPDF install may have failed. Download from https://sumatrapdfreader.org/"
    }
} else {
    Warn "Cannot install SumatraPDF (no Chocolatey). Download from https://sumatrapdfreader.org/"
}

# =====================================================================
# 6. Verification
# =====================================================================
Step "6/7" "Verification..."

if ($Rscript -and (Test-Path $Rscript)) {
    & $Rscript -e "cat('R:', R.version.string, '\n'); pkgs <- c('ggplot2','clusterProfiler','EnhancedVolcano','org.Hs.eg.db','igraph','ggVennDiagram'); for (p in pkgs) cat(sprintf('  %-25s %s', p, ifelse(requireNamespace(p, quietly=TRUE), 'OK', 'MISSING')), '\n')"
}

# =====================================================================
# Done
# =====================================================================
Write-Host "`n========================================" -ForegroundColor Green
Write-Host " Setup Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
if ($RepoDir) {
    Write-Host "  Repo:    $RepoDir" -ForegroundColor White
}
Write-Host "  Next:    Copy CSV files to data\" -ForegroundColor White
Write-Host "  Test:    Rscript.exe run_all.R --test" -ForegroundColor White
Write-Host "  Real:    Rscript.exe run_all.R" -ForegroundColor White
Write-Host ""
