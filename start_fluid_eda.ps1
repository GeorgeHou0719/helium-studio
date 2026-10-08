param([int]$Port = 5004)
$ErrorActionPreference = 'Stop'
Push-Location -LiteralPath $PSScriptRoot
try {
    $edaPython = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
    if (-not (Test-Path -LiteralPath $edaPython)) {
        python -m venv .venv
        if ($LASTEXITCODE -ne 0) { throw 'Unable to create .venv. Install Python 3.11 or newer.' }
        & $edaPython -m pip install -e .
        if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
    }
    & $edaPython -m models.fluid_eda.app --port $Port
    if ($LASTEXITCODE -ne 0) { throw 'Helium Studio exited with an error.' }
} finally {
    Pop-Location
}
