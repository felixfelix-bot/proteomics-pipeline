###############################################################################
# ONE-COMMAND WINDOWS SETUP FOR PROTEOMICS PIPELINE
#
# Run from any PowerShell window (no admin needed — it self-elevates):
#
#   irm https://raw.githubusercontent.com/c03rad0r/proteomics-pipeline/main/ansible/bootstrap-windows-all.ps1 | iex
#
# That's it. It installs R, Rtools, Git, clones the repo, installs all packages.
# Takes ~30-45 min. Go get coffee.
###############################################################################

# --- Self-elevate to Administrator if not already ---
if (-not ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
    Write-Host "Requesting Administrator privileges..." -ForegroundColor Yellow
    Start-Process powershell -Verb RunAs -ArgumentList "-NoProfile -ExecutionPolicy Bypass -Command `"irm https://raw.githubusercontent.com/c03rad0r/proteomics-pipeline/main/ansible/bootstrap-windows-all.ps1 | iex`""
    exit
}

$ErrorActionPreference = "Stop"
$ProjectDir = "C:\proteomics-pipeline"

function Step($n, $msg) { Write-Host "`n[$n] $msg" -ForegroundColor Cyan }
function OK($msg)       { Write-Host "  [OK] $msg" -ForegroundColor Green }
function Info($msg)     { Write-Host "  $msg" -ForegroundColor Gray }
function Fail($msg)     { Write-Host "  [FAIL] $msg" -ForegroundColor Red; exit 1 }

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host " Proteomics Pipeline - Auto Setup" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# 1. Chocolatey
Step "1/5" "Installing Chocolatey..."
if (-not (Get-Command choco -ErrorAction SilentlyContinue)) {
    Set-ExecutionPolicy Bypass -Scope Process -Force
    [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072
    try { Invoke-Expression ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1')) }
    catch { Fail "Chocolatey install failed: $_" }
}
OK "Chocolatey ready"

# 2. Git + R + Rtools
Step "2/5" "Installing Git, R, and Rtools..."
if (-not (Get-Command git -ErrorAction SilentlyContinue)) { choco install git -y --no-progress }
choco install r.project -y --no-progress
choco install rtools -y --no-progress
$env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
OK "Git, R, Rtools installed"

# 3. Clone repo
Step "3/5" "Cloning repository..."
if (Test-Path $ProjectDir) { Push-Location $ProjectDir; git pull; Pop-Location }
else { git clone https://github.com/c03rad0r/proteomics-pipeline.git $ProjectDir }
OK "Repository at $ProjectDir"

# 4. Find Rscript
$Rscript = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1 | ForEach-Object { Join-Path $_.FullName "bin\Rscript.exe" }
if (-not (Test-Path $Rscript)) { Fail "Rscript.exe not found. Check R installation." }
Info "Using: $Rscript"

# 5. Install R packages (the long step)
Step "4/5" "Installing R packages (20-40 minutes)..."
Info "Compiling Bioconductor packages from source. This is normal."
Push-Location $ProjectDir
& $Rscript "R\00_install_packages.R"
Pop-Location
OK "R packages installed"

# 6. Quick test
Step "5/5" "Quick verification..."
& $Rscript -e "cat('R:', R.version.string, '\n'); for (p in c('ggplot2','clusterProfiler','EnhancedVolcano','org.Hs.eg.db','igraph')) cat(p, ':', ifelse(requireNamespace(p,quietly=TRUE),'OK','MISSING'), '\n')"

Write-Host "`n========================================" -ForegroundColor Green
Write-Host " DONE! Next steps:" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host "  1. Copy CSV files to: $ProjectDir\data\" -ForegroundColor White
Write-Host "  2. Test:  cd $ProjectDir" -ForegroundColor White
Write-Host "            Rscript.exe run_all.R --test" -ForegroundColor White
Write-Host "  3. Real:  Rscript.exe run_all.R" -ForegroundColor White
Write-Host ""
