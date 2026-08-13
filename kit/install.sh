#!/bin/sh
set -eu

mode=${1:-offline}
case "$mode" in
  offline|online) ;;
  *) printf '%s\n' "usage: ./install.sh [offline|online]" >&2; exit 2 ;;
esac

kit_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
UV_CACHE_DIR="$kit_root/.cache/uv"
export UV_CACHE_DIR
platform=$(tr -d '\r\n' < "$kit_root/platform.txt")
machine_arch=$(uname -m)
case "$platform:$machine_arch" in
  macos-aarch64:arm64|macos-x86_64:x86_64) ;;
  macos-aarch64:*)
    printf '%s\n' "This kit supports Apple Silicon macOS; this machine reports $machine_arch." >&2
    exit 2
    ;;
  macos-x86_64:*)
    printf '%s\n' "This kit supports Intel macOS; this machine reports $machine_arch." >&2
    exit 2
    ;;
  *)
    printf '%s\n' "This is not a supported macOS kit payload: $platform." >&2
    exit 2
    ;;
esac

(cd "$kit_root" && shasum -a 256 -c checksums.sha256)
chmod +x "$kit_root/bin/uv"

runtime_root="$kit_root/.runtime"
python="$runtime_root/python/bin/python3"
if [ ! -f "$python" ]; then
  if [ -e "$runtime_root" ]; then
    printf '%s\n' "Cached runtime is incomplete: $runtime_root. Remove only that generated folder and rerun." >&2
    exit 2
  fi
  set -- "$kit_root"/runtime/cpython-*.tar.gz
  if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
    printf '%s\n' "Expected exactly one cached CPython archive." >&2
    exit 2
  fi
  mkdir "$runtime_root"
  tar -xzf "$1" -C "$runtime_root"
  if [ ! -f "$python" ]; then
    printf '%s\n' "Cached CPython extraction failed. Remove only $runtime_root before retrying." >&2
    exit 2
  fi
fi
python_version=$($python -c 'import platform; print(platform.python_version())')
if [ "$python_version" != "3.12.13" ]; then
  printf '%s\n' "Cached runtime version mismatch: expected 3.12.13, got $python_version." >&2
  exit 2
fi

uv="$kit_root/bin/uv"
venv="$kit_root/.venv"
venv_python="$venv/bin/python"
if [ ! -f "$venv_python" ]; then
  if [ -e "$venv" ]; then
    printf '%s\n' "Kit environment is incomplete: $venv. Remove only that generated folder and rerun." >&2
    exit 2
  fi
  "$uv" venv "$venv" --python "$python" --no-python-downloads
fi

if [ "$mode" = offline ]; then
  "$uv" pip install --python "$venv_python" --exact --require-hashes --no-index \
    --find-links "$kit_root/wheelhouse" --offline --no-python-downloads \
    --requirements "$kit_root/requirements-runtime.txt"
else
  "$uv" pip install --python "$venv_python" --exact --require-hashes \
    --no-python-downloads --requirements "$kit_root/requirements-runtime.txt"
fi

set -- "$kit_root"/wheelhouse/echelle_spectra-*.whl
if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  printf '%s\n' "Expected exactly one echelle_spectra wheel." >&2
  exit 2
fi
project_wheel=$1
set -- "$kit_root"/wheelhouse/spectrocube-*.whl
if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  printf '%s\n' "Expected exactly one pinned spectrocube wheel." >&2
  exit 2
fi
spectrocube_wheel=$1
set -- "$kit_root"/wheelhouse/sif_parser-*.whl
if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  printf '%s\n' "Expected exactly one pinned sif_parser wheel." >&2
  exit 2
fi
"$uv" pip install --python "$venv_python" --reinstall --no-deps --no-index --offline \
  "$spectrocube_wheel" "$1" "$project_wheel"
"$venv/bin/echelle" status
