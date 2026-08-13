$ErrorActionPreference = "Stop"
$Command = Join-Path $PSScriptRoot ".venv\Scripts\echelle.exe"
if (-not (Test-Path -LiteralPath $Command -PathType Leaf)) {
    throw "The NIFS kit is not installed. Run .\install.ps1 offline (or online) first."
}
& $Command @args
exit $LASTEXITCODE
