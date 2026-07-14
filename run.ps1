##############################################################################
## run.ps1 - PowerShell replacement for `make` on Windows
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
    $rVersions = Get-ChildItem 'C:\Program Files\R' -Directory -ErrorAction SilentlyContinue |
                 Sort-Object Name -Descending |
                 Select-Object -First 1
    if ($rVersions) {
        $RSCRIPT = Join-Path $rVersions.FullName 'bin\Rscript.exe'
    }
}

if (-not $RSCRIPT -or -not (Test-Path $RSCRIPT)) {
    Write-Host 'ERROR: Rscript.exe not found.' -ForegroundColor Red
    Write-Host '  Install R from https://cran.r-project.org/bin/windows/base/' -ForegroundColor Yellow
    Write-Host '  Or hardcode the path at the top of this script.' -ForegroundColor Yellow
    exit 1
}

Write-Host "Using Rscript: $RSCRIPT" -ForegroundColor Cyan

# ---- Target definitions ----
# Maps make target names to their Rscript calls.
# Group targets run multiple sub-targets in sequence.

$groupTargets = @{
    'aruna-fast' = @(
        'targeted_volcanos',
        'flagip_volcano',
        'targeted_venns',
        'targeted_go',
        'venn',
        'bidirectional_go',
        'flagip_validated_go',
        'poster_figures',
        'lydia_network_volcano'
    )
    'aruna-slow' = @(
        'lydia_network_volcano',
        'string_network',
        'crac_string_network',
        'ra_common',
        'chx_common_analysis',
        'gsea'
    )
    'aruna-all' = @(
        'targeted_volcanos',
        'flagip_volcano',
        'targeted_venns',
        'targeted_go',
        'venn',
        'bidirectional_go',
        'flagip_validated_go',
        'poster_figures',
        'lydia_network_volcano',
        'string_network',
        'crac_string_network',
        'ra_common',
        'chx_common_analysis',
        'gsea',
        'chx_kegg_crac_overlap',
        'chx_volcano_venn'
    )
    'all-volcano' = @(
        'targeted_volcanos',
        'flagip_volcano',
        'lydia_network_volcano'
    )
}

# Single-step targets that map directly to run_step.R step names
# (Make target name with hyphens -> step name with underscores)
$singleTargets = @{
    'volcano'              = 'volcano'
    'venn'                 = 'venn'
    'go'                   = 'go'
    'string-network'       = 'string_network'
    'crac-network'         = 'crac_string_network'
    'targeted-volcano'     = 'targeted_volcanos'
    'flagip-volcano'       = 'flagip_volcano'
    'targeted-venn'        = 'targeted_venns'
    'targeted-go'          = 'targeted_go'
    'bidirectional-go'     = 'bidirectional_go'
    'bidirectional-go-ra'  = 'bidirectional_go_ra'
    'flagip-validated-go'  = 'flagip_validated_go'
    'string-style-network' = 'string_style_network'
    'network-go'           = 'network_go'
    'ra-common'            = 'ra_common'
    'chx-common'           = 'chx_common_analysis'
    'chx-volcano-venn'     = 'chx_volcano_venn'
    'chx-kegg-crac'        = 'chx_kegg_crac_overlap'
    'gsea'                 = 'gsea'
    'lydia-volcano'        = 'lydia_network_volcano'
    'pathway-network'      = 'pathway_network'
    'poster-figures'       = 'poster_figures'
    'dotplot-variants'     = 'dotplot_variants'
    'style-variants'       = 'style_variants'
    'diagnostics'          = 'diagnostics'
}

# ---- Helper: run a single step ----
function Invoke-Step {
    param([string]$StepName, [switch]$Force)

    Write-Host ''
    Write-Host '=========================================' -ForegroundColor Cyan
    Write-Host " Running: $StepName" -ForegroundColor Cyan
    Write-Host '=========================================' -ForegroundColor Cyan
    Write-Host ''

    $cmdArgs = @('R/run_step.R', $StepName)
    if ($Force) { $cmdArgs += '--force' }

    # Capture both stdout and stderr so we can show the error
    $output = & $RSCRIPT @cmdArgs 2>&1
    $output | Write-Host

    if ($LASTEXITCODE -ne 0) {
        Write-Host "FAILED: $StepName (exit code $LASTEXITCODE)" -ForegroundColor Red
        return $false
    }
    return $true
}

# ---- Helper: collect figures + compile poster ----
function Invoke-CollectAndCompile {
    # ---- Collect figures from output/figures/ (RECURSIVE) ----
    $figDir = 'output\figures'
    $posterDir = 'poster\figures'
    if (-not (Test-Path $posterDir)) { New-Item -ItemType Directory -Path $posterDir | Out-Null }

    # Pattern -> clean name mapping (clean names match poster_review.tex \includegraphics)
    $figMap = @{
        'flagip_GO_validated_any_BP_dotplot'       = 'dotplot_flagip_validated_BP'
        'network_go_comparison_BP'                 = 'dotplot_network_go_BP'
        'targeted_GO_RA_shared_core_MF_dotplot'    = 'dotplot_RA_shared_MF'
        'flagip_overlap_volcano_BK467_TRIP4_vs_WT' = 'volcano_flagip_overlap'
        'lydia_network_volcano'                    = 'volcano_lydia_network'
        'targeted_volcano_BK504_RA_effect'         = 'volcano_RA_BK504'
    }

    Write-Host ''
    Write-Host 'Collecting figures (recursive search)...' -ForegroundColor Cyan
    foreach ($entry in $figMap.GetEnumerator()) {
        $pattern = $entry.Key
        $cleanName = $entry.Value

        # Search recursively for PDF (preferred) then PNG
        $found = Get-ChildItem -Path $figDir -Recurse -Filter "$($pattern)*.pdf" -ErrorAction SilentlyContinue |
                 Sort-Object LastWriteTime -Descending | Select-Object -First 1
        if (-not $found) {
            $found = Get-ChildItem -Path $figDir -Recurse -Filter "$($pattern)*.png" -ErrorAction SilentlyContinue |
                     Sort-Object LastWriteTime -Descending | Select-Object -First 1
        }

        if ($found) {
            $ext = $found.Extension  # .pdf or .png
            Copy-Item $found.FullName "$posterDir\$cleanName$ext" -Force
            Write-Host "  $($found.Name) -> $cleanName$ext"
        } else {
            Write-Host "  WARNING: Not found: $pattern" -ForegroundColor Yellow
        }
    }

    # ---- Compile LaTeX poster (if pdflatex is available) ----
    $pdflatex = Get-Command pdflatex -ErrorAction SilentlyContinue
    if (-not $pdflatex) {
        Write-Host ''
        Write-Host '=========================================' -ForegroundColor Yellow
        Write-Host ' FIGURES COLLECTED — pdflatex not installed' -ForegroundColor Yellow
        Write-Host '=========================================' -ForegroundColor Yellow
        Write-Host '  Figures are in poster\figures\' -ForegroundColor White
        Write-Host ''
        Write-Host '  To compile the poster PDF, install MiKTeX:' -ForegroundColor White
        Write-Host '    1. Download: https://miktex.org/download' -ForegroundColor White
        Write-Host '    2. Install (accept defaults)' -ForegroundColor White
        Write-Host '    3. Run: .\run.ps1 poster-review-only' -ForegroundColor White
        Write-Host ''
        Write-Host '  OR: send the figures to someone with LaTeX installed.' -ForegroundColor White
    } else {
        Write-Host ''
        Write-Host 'Compiling poster_review.tex ...' -ForegroundColor Cyan
        Push-Location 'poster'
        & pdflatex -interaction=nonstopmode poster_review.tex 2>&1 | Out-Null
        & pdflatex -interaction=nonstopmode poster_review.tex 2>&1 | Out-Null
        Pop-Location

        if (Test-Path 'poster\poster_review.pdf') {
            # Rename with commit hash so feedback references the right version
            $gitHash = (git rev-parse --short HEAD 2>$null).Trim()
            if ($gitHash) {
                $versionedPdf = "poster\poster_review_${gitHash}.pdf"
                Copy-Item 'poster\poster_review.pdf' $versionedPdf -Force
                Write-Host ''
                Write-Host '=========================================' -ForegroundColor Green
                Write-Host " POSTER REVIEW: $versionedPdf" -ForegroundColor Green
                Write-Host " (commit: $gitHash)" -ForegroundColor White
                Write-Host '=========================================' -ForegroundColor Green
                Start-Process $versionedPdf
            } else {
                Write-Host ''
                Write-Host '=========================================' -ForegroundColor Green
                Write-Host ' POSTER REVIEW: poster\poster_review.pdf' -ForegroundColor Green
                Write-Host '=========================================' -ForegroundColor Green
                Start-Process 'poster\poster_review.pdf'
            }
        } else {
            Write-Host 'LaTeX compilation failed.' -ForegroundColor Red
        }
    }
}

# ---- Helper: check data files exist ----
function Test-DataFiles {
    $dataFiles = Get-ChildItem -Path 'data' -Recurse -Filter '*diffEx_minProb*.csv' -ErrorAction SilentlyContinue
    if (-not $dataFiles -or $dataFiles.Count -eq 0) {
        Write-Host '=========================================' -ForegroundColor Red
        Write-Host ' ERROR: No data files found in data\' -ForegroundColor Red
        Write-Host '=========================================' -ForegroundColor Red
        Write-Host ''
        Write-Host 'The pipeline needs CSV files ending in _diffEx_minProb.csv' -ForegroundColor Yellow
        Write-Host 'in the data\ directory (or subdirectories).' -ForegroundColor Yellow
        Write-Host ''
        Write-Host 'Copy your mass spec output files to data\ and try again.' -ForegroundColor Yellow
        Write-Host ''
        $dataContents = Get-ChildItem -Path 'data' -Recurse -ErrorAction SilentlyContinue
        if ($dataContents) {
            Write-Host 'Current contents of data\:' -ForegroundColor Gray
            $dataContents | ForEach-Object { Write-Host "  $($_.FullName.Replace((Get-Location).Path, ''))" }
        } else {
            Write-Host 'data\ directory is empty or does not exist.' -ForegroundColor Gray
        }
        return $false
    }
    Write-Host "Found $($dataFiles.Count) data file(s) in data\" -ForegroundColor Green
    Write-Host ''
    return $true
}

# ---- Check packages first ----
Write-Host 'Checking R packages...' -ForegroundColor Cyan
& $RSCRIPT 'R/check_packages.R' 2>&1 | Write-Host
if ($LASTEXITCODE -ne 0) {
    Write-Host 'Package check failed. Run: .\run.ps1 install' -ForegroundColor Red
    exit 1
}

# ---- Dispatch ----
$force = $RemainingArgs -contains '--force'

switch ($Target) {
    'install' {
        Write-Host 'Installing R packages...' -ForegroundColor Cyan
        & $RSCRIPT 'R/00_install_packages.R'
        break
    }

    'setup' {
        Write-Host '=========================================' -ForegroundColor Cyan
        Write-Host ' FULL SETUP: Git + R + Rtools + R packages + MiKTeX' -ForegroundColor Cyan
        Write-Host '=========================================' -ForegroundColor Cyan
        Write-Host ''
        Write-Host 'This requires Administrator privileges. UAC prompt will appear.' -ForegroundColor Yellow
        Write-Host 'Takes 20-30 minutes on first run.' -ForegroundColor White
        Write-Host ''
        $scriptPath = Join-Path $PSScriptRoot 'ansible\bootstrap-windows-all.ps1'
        Start-Process powershell -Verb RunAs -Wait -ArgumentList "-NoProfile -ExecutionPolicy Bypass -File `"$scriptPath`""
        Write-Host ''
        Write-Host 'Setup script finished. Restart PowerShell to pick up new PATH entries.' -ForegroundColor Green
        break
    }

    { $groupTargets.ContainsKey($_) } {
        $steps = $groupTargets[$_]
        $failed = @()
        foreach ($step in $steps) {
            $ok = Invoke-Step -StepName $step -Force:$force
            if (-not $ok) { $failed += $step }
        }
        Write-Host ''
        if ($failed.Count -gt 0) {
            Write-Host "COMPLETED WITH ERRORS. Failed steps: $($failed -join ', ')" -ForegroundColor Yellow
            exit 1
        } else {
            Write-Host 'ALL STEPS COMPLETED SUCCESSFULLY.' -ForegroundColor Green
            Write-Host 'Check output/figures/ for results.' -ForegroundColor Green
        }
        break
    }

    { $singleTargets.ContainsKey($_) } {
        $stepName = $singleTargets[$_]
        $ok = Invoke-Step -StepName $stepName -Force:$force
        if ($ok) {
            Write-Host 'Step completed. Check output/figures/' -ForegroundColor Green
        }
        break
    }

    'poster-review' {
        # Step 1: Check data files
        if (-not (Test-DataFiles)) { break }

        # Step 2: Generate ONLY poster figures (run_step.R handles caching)
        # ra_common removed — its outputs aren't in the poster figMap.
        # Saves ~10 min on first run (no STRING networks, Venn, circular, etc.)
        $posterSteps = @(
            'targeted_volcanos',
            'flagip_volcano',
            'lydia_network_volcano',
            'flagip_validated_go',
            'network_go',
            'targeted_go'
        )
        $failed = @()
        foreach ($step in $posterSteps) {
            $ok = Invoke-Step -StepName $step -Force:$force
            if (-not $ok) { $failed += $step }
        }
        Write-Host ''
        if ($failed.Count -gt 0) {
            Write-Host "Some steps failed: $($failed -join ', ')" -ForegroundColor Yellow
            Write-Host 'Continuing with available figures...' -ForegroundColor Yellow
        }

        # Step 3: Collect + compile
        Invoke-CollectAndCompile
        break
    }

    'poster-review-force' {
        # Same as poster-review but ignores cache — regenerates everything
        if (-not (Test-DataFiles)) { break }
        Write-Host "Regenerating ALL figures (no cache)..." -ForegroundColor Green

        $posterSteps = @(
            'targeted_volcanos',
            'flagip_volcano',
            'lydia_network_volcano',
            'flagip_validated_go',
            'network_go',
            'targeted_go'
        )
        $failed = @()
        foreach ($step in $posterSteps) {
            $ok = Invoke-Step -StepName $step -Force
            if (-not $ok) { $failed += $step }
        }
        if ($failed.Count -gt 0) {
            Write-Host "Some steps failed: $($failed -join ', ')" -ForegroundColor Yellow
        }

        Invoke-CollectAndCompile
        break
    }

    'poster-review-only' {
        # Skip figure regeneration — just collect + compile (fast iteration)
        Invoke-CollectAndCompile
        break
    }

    'help' {
        Write-Host 'TRIP4/ASCC Proteomics Pipeline - PowerShell Runner'
        Write-Host ''
        Write-Host 'Group targets:'
        Write-Host '  .\run.ps1 aruna-fast       Quick analyses + poster figures'
        Write-Host '  .\run.ps1 aruna-slow        STRING networks + GSEA (slow)'
        Write-Host '  .\run.ps1 aruna-all         Everything'
        Write-Host ''
        Write-Host 'Single steps:'
        Write-Host '  .\run.ps1 targeted-volcano   Volcano plots'
        Write-Host '  .\run.ps1 targeted-go        GO enrichment'
        Write-Host '  .\run.ps1 targeted-venn      Venn diagrams'
        Write-Host '  .\run.ps1 bidirectional-go   Bidirectional GO dot plot'
        Write-Host '  .\run.ps1 lydia-volcano      Lydia STRING overlay volcano'
        Write-Host '  .\run.ps1 string-network     STRING PPI network'
        Write-Host '  .\run.ps1 flagip-validated-go  Flag-validated GO+KEGG'
        Write-Host '  .\run.ps1 poster-figures     Poster-styled figures'
        Write-Host '  .\run.ps1 poster-review       All 6 poster figures + compile (cached)'
        Write-Host '  .\run.ps1 poster-review-only   Just collect + compile (fast, skip R)'
        Write-Host '  .\run.ps1 poster-review-force  Regenerate ALL figures + compile (no cache)'
        Write-Host '  .\run.ps1 dotplot-variants   GO dotplot font-size variants'
        Write-Host '  .\run.ps1 style-variants     12 style variants in one PDF (pick your favorite)'
        Write-Host '  .\run.ps1 diagnostics        Structural data summary'
        Write-Host ''
        Write-Host 'Caching:'
        Write-Host '  Steps are skipped if code unchanged since last run.'
        Write-Host '  Add --force to any target to bypass cache:'
        Write-Host '    .\run.ps1 poster-review --force'
        Write-Host '    .\run.ps1 targeted-go --force'
        Write-Host ''
        Write-Host 'Setup:'
        Write-Host '  .\run.ps1 setup              Full setup: Git + R + Rtools + packages + MiKTeX'
        Write-Host '  .\run.ps1 install            Install R packages only'
        break
    }

    default {
        Write-Host "Unknown target: $Target" -ForegroundColor Red
        Write-Host 'Run: .\run.ps1 help' -ForegroundColor Yellow
        exit 1
    }
}
