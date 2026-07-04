###############################################################################
# bootstrap-windows.ps1
# Configures WinRM on Windows so Ansible can manage it remotely.
# Run this ONCE on the Windows target machine, as Administrator.
###############################################################################
#  Usage:
#    1. Open PowerShell as Administrator
#    2. Set-ExecutionPolicy RemoteSigned -Scope CurrentUser
#    3. .\bootstrap-windows.ps1
###############################################################################

$ErrorActionPreference = "Stop"

Write-Host ""
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host " Proteomics Pipeline - WinRM Bootstrap" -ForegroundColor Cyan
Write-Host " Windows Target Machine Setup" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

# Check for Administrator privileges
$isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
if (-not $isAdmin) {
    Write-Host "[ERROR] This script must be run as Administrator." -ForegroundColor Red
    Write-Host "   Right-click PowerShell -> 'Run as Administrator'" -ForegroundColor Yellow
    exit 1
}

Write-Host "[1/5] Configuring WinRM service..." -ForegroundColor Yellow

# Enable WinRM
$winrmQuickConfig = winrm quickconfig -quiet 2>&1
if ($LASTEXITCODE -ne 0) {
    # Sometimes needs the service started first
    Set-Service -Name WinRM -StartupType Automatic
    Start-Service -Name WinRM
    winrm quickconfig -quiet 2>&1 | Out-Null
}

Write-Host "[2/5] Setting WinRM service to auto-start..." -ForegroundColor Yellow
Set-Service -Name WinRM -StartupType Automatic

Write-Host "[3/5] Configuring WinRM for NTLM authentication..." -ForegroundColor Yellow

# Allow NTLM auth (needed for Ansible)
winrm set winrm/config/service/Auth '@{Basic="true";NTLM="true";CredSSP="true"}' 2>&1 | Out-Null
winrm set winrm/config/client/Auth '@{Basic="true";NTLM="true";CredSSP="true"}' 2>&1 | Out-Null

# Set AllowUnencrypted for NTLM (Ansible requirement)
winrm set winrm/config/service '@{AllowUnencrypted="true"}' 2>&1 | Out-Null
winrm set winrm/config/client '@{AllowUnencrypted="true"}' 2>&1 | Out-Null

Write-Host "[4/5] Opening firewall for WinRM (port 5985 HTTP)..." -ForegroundColor Yellow

# Open firewall for WinRM HTTP (port 5985) and HTTPS (port 5986)
netsh advfirewall firewall add rule name="WinRM HTTP" dir=in action=allow protocol=TCP localport=5985 2>&1 | Out-Null
netsh advfirewall firewall add rule name="WinRM HTTPS" dir=in action=allow protocol=TCP localport=5986 2>&1 | Out-Null

Write-Host "[5/5] Testing WinRM listener..." -ForegroundColor Yellow

# Test
$listener = winrm enumerate winrm/config/listener 2>&1
if ($listener -match "HTTP") {
    Write-Host ""
    Write-Host "[SUCCESS] WinRM is configured and listening!" -ForegroundColor Green
} else {
    Write-Host "[WARNING] WinRM listener not detected on HTTP. Check configuration." -ForegroundColor Yellow
}

# Get this machine's IP for reference
$ip = (Get-NetIPAddress -AddressFamily IPv4 | Where-Object { $_.InterfaceAlias -notlike "*Loopback*" -and $_.IPAddress -notlike "169.*" } | Select-Object -First 1).IPAddress

Write-Host ""
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host " WinRM Setup Complete!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "This machine IP: $ip" -ForegroundColor White
Write-Host ""
Write-Host "Next steps (from the Linux control machine):" -ForegroundColor White
Write-Host "  1. Edit ansible/inventory.yml:" -ForegroundColor Gray
Write-Host "     ansible_host: $ip" -ForegroundColor Gray
Write-Host "     ansible_user: $env:USERNAME" -ForegroundColor Gray
Write-Host "     ansible_password: <your-windows-password>" -ForegroundColor Gray
Write-Host ""
Write-Host "  2. Run the playbook:" -ForegroundColor Gray
Write-Host "     cd ansible && ansible-playbook -i inventory.yml setup.yml" -ForegroundColor Gray
Write-Host ""

# Also recommend creating the working directory
$workDir = "C:\proteomics-pipeline"
if (-not (Test-Path $workDir)) {
    New-Item -ItemType Directory -Path $workDir -Force | Out-Null
    Write-Host "Created working directory: $workDir" -ForegroundColor Green
}
