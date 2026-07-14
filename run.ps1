##############################################################################
## run.ps1 — PowerShell replacement for `make` on Windows
##
## USAGE:
##   .\run.ps1 aruna-fast       (same as: make aruna-fast)
##   .\run.ps1 aruna-slow       (same as: make aruna-slow)
##   .\run.ps1 aruna-all        (same as: make aruna-all)
##   .\run.ps1 targeted-volcano (run a single analysis step)
##   .\run.ps1 poster-figures   (generate poster-styled figures)
##
## WHY THIS EXISTS:
##   Windows Application Control policy blocks make.exe.
##   This script calls Rscript directly, bypassing make entirely.
##   No WSL needed, no make.exe needed, no policy exceptions needed.
##
## PREREQUISITE:
##   R must be installed and Rscript.exe must be on your PATH.
##   If Rscript is not on PATH, set the full path below:
##   $RSCRIPT = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
##############################################################################

param(
    [Parameter(Mandatory=$true, Position=0)]
    [string]$Target,

    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$RemainingArgs
)

# ---- Rscript path ----
# Auto-detect, or hardcode if auto-detect fails:
$RSCRIPT = $null

# Try PATH first
$RSCRIPT = (Get-Command Rscript.exe -ErrorAction SilentlyContinue).Source

# If not on PATH, try known install location
if (-not $RSCRIPT) {
    $rVersions = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue |
                 Sort-Object Name -Descending |
                 Select-Object -First 1
    if ($rVersions) {
        $RSCRIPT = Join-Path $rVersions.FullName "bin\Rscript.exe"
    }
}

if (-not $RSCRIPT -or -not (Test-Path $RSCRIPT)) {
    Write-Host "ERROR: Rscript.exe not found." -ForegroundColor Red
    Write-Host "  Install R from https://cran.r-project.org/bin/windows/base/" -ForegroundColor Yellow
    Write-Host "  Or hardcode the path at the top of this script." -ForegroundColor Yellow
    exit 1
}

Write-Host "Using Rscript: $RSCRIPT" -ForegroundColor Cyan

# ---- Target definitions ----
# Maps make target names to their Rscript calls.
# Group targets run multiple sub-targets in sequence.

$groupTargets = @{
    "aruna-fast" = @(
        "targeted_volcanos",
        "flagip_volcano",
        "targeted_venns",
        "targeted_go",
        "venn",
        "bidirectional_go",
        "flagip_validated_go",
        "poster_figures",
        "lydia_network_volcano"
    )
    "aruna-slow" = @(
        "lydia_network_volcano",
        "string_network",
        "crac_string_network",
        "ra_common",
        "chx_common_analysis",
        "gsea"
    )
    "aruna-all" = @(
        "targeted_volcanos",
        "flagip_volcano",
        "targeted_venns",
        "targeted_go",
        "venn",
        "bidirectional_go",
        "flagip_validated_go",
        "poster_figures",
        "lydia_network_volcano",
        "string_network",
        "crac_string_network",
        "ra_common",
        "chx_common_analysis",
        "gsea",
        "chx_kegg_crac_overlap",
        "chx_volcano_venn"
    )
    "all-volcano" = @(
        "targeted_volcanos",
        "flagip_volcano",
        "lydia_network_volcano"
    )
}

# Single-step targets that map directly to run_step.R step names
# (Make target name with hyphens → step name with underscores)
$singleTargets = @{
    "volcano"              = "volcano"
    "venn"                 = "venn"
    "go"                   = "go"
    "string-network"       = "string_network"
    "crac-network"         = "crac_string_network"
    "targeted-volcano"     = "targeted_volcanos"
    "flagip-volcano"       = "flagip_volcano"
    "targeted-venn"        = "targeted_venns"
    "targeted-go"          = "targeted_go"
    "bidirectional-go"     = "bidirectional_go"
    "bidirectional-go-ra"  = "bidirectional_go_ra"
    "flagip-validated-go"  = "flagip_validated_go"
    "string-style-network" = "string_style_network"
    "network-go"           = "network_go"
    "ra-common"            = "ra_common"
    "chx-common"           = "chx_common_analysis"
    "chx-volcano-venn"     = "chx_volcano_venn"
    "chx-kegg-crac"        = "chx_kegg_crac_overlap"
    "gsea"                 = "gsea"
    "lydia-volcano"        = "lydia_network_volcano"
    "pathway-network"      = "pathway_network"
    "poster-figures"       = "poster_figures"
    "diagnostics"          = "diagnostics"
}

# ---- Helper: run a single step ----
function Invoke-Step {
    param([string]$StepName)

    Write-Host ""
    Write-Host "=========================================" -ForegroundColor Cyan
    Write-Host " Running: $StepName" -ForegroundColor Cyan
    Write-Host "=========================================" -ForegroundColor Cyan
    Write-Host ""

    & $RSCRIPT "R/run_step.R" $StepName

    if ($LASTEXITCODE -ne 0) {
        Write-Host "FAILED: $StepName (exit code $LASTEXITCODE)" -ForegroundColor Red
        return $false
    }
    return $true
}

# ---- Check packages first ----
Write-Host "Checking R packages..." -ForegroundColor Cyan
& $RSCRIPT "R/check_packages.R" 2>&1 | Write-Host
if ($LASTEXITCODE -ne 0) {
    Write-Host "Package check failed. Run: .\run.ps1 install" -ForegroundColor Red
    exit 1
}

# ---- Dispatch ----
switch ($Target) {
    "install" {
        Write-Host "Installing R packages..." -ForegroundColor Cyan
        & $RSCRIPT "R/00_install_packages.R"
        break
    }

    { $groupTargets.ContainsKey($_) } {
        $steps = $groupTargets[$_]
        $failed = @()
        foreach ($step in $steps) {
            $ok = Invoke-Step -StepName $step
            if (-not $ok) { $failed += $step }
        }
        Write-Host ""
        if ($failed.Count -gt 0) {
            Write-Host "COMPLETED WITH ERRORS. Failed steps: $($failed -join ', ')" -ForegroundColor Yellow
            exit 1
        } else {
            Write-Host "ALL STEPS COMPLETED SUCCESSFULLY." -ForegroundColor Green
            Write-Host "Check output/figures/ for results." -ForegroundColor Green
        }
        break
    }

    { $singleTargets.ContainsKey($_) } {
        $stepName = $singleTargets[$_]
        $ok = Invoke-Step -StepName $stepName
        if ($ok) {
            Write-Host "Step completed. Check output/figures/" -ForegroundColor Green
        }
        break
    }

    "help" {
        Write-Host "TRIP4/ASCC Proteomics Pipeline — PowerShell Runner"
        Write-Host ""
        Write-Host "Group targets:"
        Write-Host "  .\run.ps1 aruna-fast       Quick analyses + poster figures"
        Write-Host "  .\run.ps1 aruna-slow        STRING networks + GSEA (slow)"
        Write-Host "  .\run.ps1 aruna-all         Everything"
        Write-Host ""
        Write-Host "Single steps:"
        Write-Host "  .\run.ps1 targeted-volcano   Volcano plots"
        Write-Host "  .\run.ps1 targeted-go        GO enrichment"
        Write-Host "  .\run.ps1 targeted-venn      Venn diagrams"
        Write-Host "  .\run.ps1 bidirectional-go   Bidirectional GO dot plot"
        Write-Host "  .\run.ps1 lydia-volcano      Lydia STRING overlay volcano"
        Write-Host "  .\run.ps1 string-network     STRING PPI network"
        Write-Host "  .\run.ps1 flagip-validated-go  Flag-validated GO+KEGG"
        Write-Host "  .\run.ps1 poster-figures     Poster-styled figures"
        Write-Host "  .\run.ps1 diagnostics        Structural data summary"
        Write-Host ""
        Write-Host "Setup:"
        Write-Host "  .\run.ps1 install            Install R packages"
        break
    }

    default {
        Write-Host "Unknown target: $Target" -ForegroundColor Red
        Write-Host "Run: .\run.ps1 help" -ForegroundColor Yellow
        exit 1
    }
}
