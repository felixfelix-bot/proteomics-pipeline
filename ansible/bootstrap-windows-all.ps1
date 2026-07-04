###############################################################################
# bootstrap-windows-all.ps1
# ALL-IN-ONE setup for the researcher's Windows laptop.
# Installs everything needed to run the proteomics pipeline:
#   1. R + Rtools (via winget/chocolatey)
#   2. Git (if not present)
#   3. Clones the repo
#   4. Installs all R packages
#
# Usage:
#   Open PowerShell as Administrator
#   Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
#   .\bootstrap-windows-all.ps1
#
# Total time: ~30-45 minutes (mostly R package compilation)
###############################################################################

$ErrorActionPreference = "Stop"

function Write-Step($msg) {
    Write-Host "`n[STEP] $msg" -ForegroundColor Cyan
}
function Write-OK($msg) {
    Write-Host "  [OK] $msg" -ForegroundColor Green
}
function Write-Info($msg) {
    Write-Host "  $msg" -ForegroundColor Gray
}
function Write-Warn($msg) {
    Write-Host "  [WARNING] $msg" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host " Proteomics Pipeline - Full Windows Setup" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# Check Administrator
$isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $isAdmin) {
    Write-Warn "Not running as Administrator. Some installs may fail."
    Write-Host "  Right-click PowerShell -> 'Run as Administrator'" -ForegroundColor Yellow
    $continue = Read-Host "Continue anyway? (y/N)"
    if ($continue -ne "y") { exit 1 }
}

# ---- STEP 1: Install Chocolatey (if not present) ----
Write-Step "1/6: Checking Chocolatey..."
if (-not (Get-Command choco -ErrorAction SilentlyContinue)) {
    Write-Info "Installing Chocolatey..."
    Set-ExecutionPolicy Bypass -Scope Process -Force
    [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
    Invoke-Expression ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))
    Write-OK "Chocolatey installed"
} else {
    Write-OK "Chocolatey already installed"
}

# ---- STEP 2: Install Git ----
Write-Step "2/6: Installing Git..."
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    choco install git -y --no-progress
    Write-OK "Git installed"
    # Refresh PATH
    $env:Path = [System.Environment]::GetEnvironmentVariable("Path", "Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path", "User")
} else {
    Write-OK "Git already installed"
}

# ---- STEP 3: Install R ----
Write-Step "3/6: Installing R..."
$rInstalled = $false
$rscriptPath = $null

# Check common locations
$rPaths = @(
    "C:\Program Files\R",
    "${env:LOCALAPPDATA}\Programs\R"
)

foreach ($basePath in $rPaths) {
    if (Test-Path $basePath) {
        $rVersion = Get-ChildItem $basePath -Directory | Sort-Object Name -Descending | Select-Object -First 1
        if ($rVersion) {
            $rscriptPath = Join-Path $rVersion.FullName "bin\Rscript.exe"
            if (Test-Path $rscriptPath) {
                $rInstalled = $true
                Write-OK "R found: $($rVersion.Name)"
                break
            }
        }
    }
}

if (-not $rInstalled) {
    Write-Info "Installing R via Chocolatey..."
    choco install r.project -y --no-progress
    # Find the installed R
    $rDir = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
    if ($rDir) {
        $rscriptPath = Join-Path $rDir.FullName "bin\Rscript.exe"
        Write-OK "R installed: $($rDir.Name)"
    } else {
        Write-Host "  [ERROR] R installation failed. Install manually from cran.r-project.org" -ForegroundColor Red
        exit 1
    }
}

Write-Info "Rscript path: $rscriptPath"

# ---- STEP 4: Install Rtools ----
Write-Step "4/6: Installing Rtools..."
$rtoolsPath = "C:\rtools43"
if (-not (Test-Path $rtoolsPath)) {
    Write-Info "Installing Rtools (needed for compiling R packages)..."
    choco install rtools -y --no-progress
    if (Test-Path "C:\rtools43") {
        Write-OK "Rtools installed"
    } else {
        # Try rtools44
        if (Test-Path "C:\rtools44") {
            $rtoolsPath = "C:\rtools44"
            Write-OK "Rtools44 installed"
        } else {
            Write-Warn "Rtools may not have installed correctly."
            Write-Host "  Download manually from https://cran.r-project.org/bin/windows/Rtools/" -ForegroundColor Yellow
        }
    }
} else {
    Write-OK "Rtools already installed"
}

# ---- STEP 5: Clone the repo ----
Write-Step "5/6: Cloning proteomics-pipeline..."
$projectDir = "C:\proteomics-pipeline"
if (Test-Path $projectDir) {
    Write-Info "Directory exists. Pulling latest..."
    Push-Location $projectDir
    git pull
    Pop-Location
} else {
    git clone https://github.com/c03rad0r/proteomics-pipeline.git $projectDir
    Write-OK "Repository cloned to $projectDir"
}

# ---- STEP 6: Install R packages ----
Write-Step "6/6: Installing R packages (this takes 20-40 minutes)..."
Write-Info "Running: Rscript R\00_install_packages.R"
Write-Info "This downloads and compiles Bioconductor packages."
Write-Info "Go get coffee. Seriously."

Push-Location $projectDir
& $rscriptPath "R\00_install_packages.R"
$installExit = $LASTEXITCODE
Pop-Location

if ($installExit -eq 0) {
    Write-OK "R packages installed successfully!"
} else {
    Write-Warn "Package installer exited with code $installExit"
    Write-Info "Some packages may have failed. Check the output above."
    Write-Info "You can re-run: Rscript R\00_install_packages.R"
}

# ---- DONE ----
Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host " Setup Complete!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "Project directory: $projectDir" -ForegroundColor White
Write-Host ""
Write-Host "Next steps:" -ForegroundColor White
Write-Host "  1. Copy your CSV files to: $projectDir\data\" -ForegroundColor Gray
Write-Host "  2. Test with synthetic data:" -ForegroundColor Gray
Write-Host "     cd $projectDir" -ForegroundColor Gray
Write-Host "     Rscript.exe run_all.R --test" -ForegroundColor Gray
Write-Host "  3. Run on real data:" -ForegroundColor Gray
Write-Host "     Rscript.exe run_all.R" -ForegroundColor Gray
Write-Host ""

# Verify R works
Write-Info "Verifying R installation..."
& $rscriptPath -e "cat('R version:', R.version.string, '\n')"
