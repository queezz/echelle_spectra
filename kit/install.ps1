[CmdletBinding()]
param(
    [ValidateSet("offline", "online")]
    [string]$Mode = "offline"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
$KitRoot = $PSScriptRoot
$env:UV_CACHE_DIR = Join-Path $KitRoot ".cache\uv"

if ((Get-Content -Raw -LiteralPath (Join-Path $KitRoot "platform.txt")).Trim() -ne "windows-x86_64") {
    throw "This kit is not the Windows x86-64 payload. Use the kit matching this machine."
}
$MachineArch = if ($env:PROCESSOR_ARCHITEW6432) {
    $env:PROCESSOR_ARCHITEW6432
} else {
    $env:PROCESSOR_ARCHITECTURE
}
if ($MachineArch -ne "AMD64") {
    throw "This kit supports Windows x86-64 (AMD64); this machine reports '$MachineArch'."
}

$ChecksumFile = Join-Path $KitRoot "checksums.sha256"
if (-not (Test-Path -LiteralPath $ChecksumFile -PathType Leaf)) {
    throw "Kit checksum inventory is missing: $ChecksumFile"
}
foreach ($Line in Get-Content -LiteralPath $ChecksumFile) {
    if ($Line -notmatch '^([0-9a-f]{64})  (.+)$') {
        throw "Invalid checksum inventory line: $Line"
    }
    $Relative = $Matches[2]
    if ([IO.Path]::IsPathRooted($Relative) -or $Relative.Split('/').Contains('..')) {
        throw "Unsafe path in checksum inventory: $Relative"
    }
    $Payload = Join-Path $KitRoot ($Relative -replace '/', '\')
    if (-not (Test-Path -LiteralPath $Payload -PathType Leaf)) {
        throw "Kit payload is incomplete; missing: $Relative"
    }
    $Actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $Payload).Hash.ToLowerInvariant()
    if ($Actual -ne $Matches[1]) {
        throw "Kit payload checksum mismatch: $Relative"
    }
}

$Uv = Join-Path $KitRoot "bin\uv.exe"
$RuntimeRoot = Join-Path $KitRoot ".runtime"
$Python = Join-Path $RuntimeRoot "python\python.exe"
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    if (Test-Path -LiteralPath $RuntimeRoot) {
        throw "Cached runtime is incomplete: $RuntimeRoot. Remove only that generated folder and rerun."
    }
    $Archives = @(Get-ChildItem -LiteralPath (Join-Path $KitRoot "runtime") -Filter "cpython-*.tar.gz" -File)
    if ($Archives.Count -ne 1) {
        throw "Expected exactly one cached CPython archive; found $($Archives.Count)."
    }
    New-Item -ItemType Directory -Path $RuntimeRoot | Out-Null
    & tar.exe -xzf $Archives[0].FullName -C $RuntimeRoot
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $Python -PathType Leaf)) {
        throw "Cached CPython extraction failed. Remove only '$RuntimeRoot' before retrying."
    }
}
$PythonVersion = (& $Python -c "import platform; print(platform.python_version())").Trim()
if ($PythonVersion -ne "3.12.13") {
    throw "Cached runtime version mismatch: expected 3.12.13, got $PythonVersion."
}

$Venv = Join-Path $KitRoot ".venv"
$VenvPython = Join-Path $Venv "Scripts\python.exe"
if (-not (Test-Path -LiteralPath $VenvPython -PathType Leaf)) {
    if (Test-Path -LiteralPath $Venv) {
        throw "Kit environment is incomplete: $Venv. Remove only that generated folder and rerun."
    }
    & $Uv venv $Venv --python $Python --no-python-downloads
    if ($LASTEXITCODE -ne 0) { throw "uv could not create the kit environment." }
}

$Requirements = Join-Path $KitRoot "requirements-runtime.txt"
$Wheelhouse = Join-Path $KitRoot "wheelhouse"
if ($Mode -eq "offline") {
    & $Uv pip install --python $VenvPython --exact --require-hashes --no-index `
        --find-links $Wheelhouse --offline --no-python-downloads --requirements $Requirements
} else {
    & $Uv pip install --python $VenvPython --exact --require-hashes `
        --no-python-downloads --requirements $Requirements
}
if ($LASTEXITCODE -ne 0) { throw "Pinned runtime dependency installation failed in $Mode mode." }

$ProjectWheels = @(Get-ChildItem -LiteralPath $Wheelhouse -Filter "echelle_spectra-*.whl" -File)
if ($ProjectWheels.Count -ne 1) {
    throw "Expected exactly one echelle_spectra wheel; found $($ProjectWheels.Count)."
}
$CompanionWheels = @(Get-ChildItem -LiteralPath $Wheelhouse -Filter "spectrocube-*.whl" -File)
if ($CompanionWheels.Count -ne 1) {
    throw "Expected exactly one pinned spectrocube wheel; found $($CompanionWheels.Count)."
}
$SifParserWheels = @(Get-ChildItem -LiteralPath $Wheelhouse -Filter "sif_parser-*.whl" -File)
if ($SifParserWheels.Count -ne 1) {
    throw "Expected exactly one pinned sif_parser wheel; found $($SifParserWheels.Count)."
}
& $Uv pip install --python $VenvPython --reinstall --no-deps --no-index --offline `
    $CompanionWheels[0].FullName $SifParserWheels[0].FullName $ProjectWheels[0].FullName
if ($LASTEXITCODE -ne 0) { throw "Local project/companion wheel installation failed." }

& (Join-Path $Venv "Scripts\echelle.exe") status
if ($LASTEXITCODE -ne 0) { throw "Installation finished, but 'echelle status' failed." }
