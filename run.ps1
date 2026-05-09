# Node Runner PowerShell launcher - pins the bundled venv Python.
# See run.bat for why this exists.

$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$py   = Join-Path $here 'venv\Scripts\python.exe'

if (-not (Test-Path $py)) {
    Write-Host "Could not find venv Python at $py." -ForegroundColor Red
    Write-Host "Create the venv first:" -ForegroundColor Yellow
    Write-Host "  python -m venv venv"
    Write-Host "  venv\Scripts\pip install -r requirements.txt"
    exit 1
}

& $py (Join-Path $here 'run.py') @args
