"""Platform-independent tests for Packet 6 kit contracts and refusal paths."""

from __future__ import annotations

import hashlib
import re
import subprocess
import zipfile
from pathlib import Path

import pytest
import tomllib

from scripts.nifs_kit import (
    Asset,
    KitError,
    assert_external_path,
    copy_posix_script,
    copy_wheelhouse,
    extract_uv,
    fetch_asset,
    inspect_project_wheel,
    load_manifest,
    validate_wheelhouse,
    verify_checksums,
    write_checksums,
)

ROOT = Path(__file__).parents[1]
MANIFEST = ROOT / "kit" / "manifest.toml"


def _wheel(path: Path, name: str, version: str) -> Path:
    normalized = name.replace("-", "_")
    with zipfile.ZipFile(path, "w") as bundle:
        bundle.writestr(
            f"{normalized}-{version}.dist-info/METADATA",
            f"Metadata-Version: 2.4\nName: {name}\nVersion: {version}\n",
        )
    return path


def _requirement_blocks(text: str) -> list[str]:
    blocks: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        if line and not line.startswith(" "):
            if current:
                blocks.append("\n".join(current))
            current = [line]
        elif current:
            current.append(line)
    if current:
        blocks.append("\n".join(current))
    return blocks


def test_manifest_pins_supported_matrix_and_project_identity() -> None:
    manifest = load_manifest(MANIFEST)
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert manifest.application_version == pyproject["project"]["version"] == "1.6.0"
    assert manifest.toolchain["python"] == (ROOT / ".python-version").read_text().strip()
    assert set(manifest.platforms) == {
        "windows-x86_64",
        "macos-aarch64",
        "macos-x86_64",
    }
    assert manifest.toolchain == {
        "python": "3.12.13",
        "python_tag": "cp312",
        "uv": "0.11.31",
        "pip": "25.2",
        "build": "1.5.0",
        "setuptools": "84.0.0",
        "wheel": "0.48.0",
        "twine": "7.0.0",
        "source_date_epoch": 1786665600,
    }
    assert {
        key: platform.minimum_os_version
        for key, platform in manifest.platforms.items()
    } == {
        "windows-x86_64": "10.0",
        "macos-aarch64": "14.0",
        "macos-x86_64": "13.0",
    }
    assert "macosx_14_0_arm64" in manifest.platforms["macos-aarch64"].pip_platforms
    for platform in manifest.platforms.values():
        for asset in (platform.uv, platform.python):
            assert asset.url.startswith("https://github.com/astral-sh/")
            assert asset.size > 20_000_000
            assert re.fullmatch(r"[0-9a-f]{64}", asset.sha256)


def test_companion_source_is_the_same_exact_commit_as_uv_lock_source() -> None:
    manifest = load_manifest(MANIFEST)
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    companion = next(item for item in manifest.companions if item.name == "spectrocube")
    sif_parser = next(item for item in manifest.companions if item.name == "sif_parser")
    source = pyproject["tool"]["uv"]["sources"]["spectrocube"]
    assert companion.name == "spectrocube"
    assert companion.version == "0.2.0"
    assert pyproject["project"]["optional-dependencies"]["spectrocube"] == [
        "spectrocube>=0.2.0"
    ]
    assert companion.source_url == source["git"]
    assert companion.source_ref == source["rev"]
    assert sif_parser.version == "0.3.6"
    assert sif_parser.source_ref == (
        "sha256:5faf1e156b1aef7835968bc3df651ad828c4d42f1c2c60a5bd26784519c2a111"
    )
    assert sif_parser.source_size == 18722


def test_runtime_exports_are_exact_hashed_and_free_of_local_sources() -> None:
    files = [ROOT / "kit" / "requirements-runtime.txt"] + sorted(
        (ROOT / "kit" / "requirements").glob("*.txt")
    )
    names_by_target: dict[str, set[str]] = {}
    for path in files:
        blocks = _requirement_blocks(path.read_text(encoding="utf-8"))
        assert blocks
        names = set()
        for block in blocks:
            match = re.match(r"^([A-Za-z0-9_.-]+)==([^ ;\\]+)", block)
            assert match is not None, block
            names.add(match.group(1).lower().replace("_", "-"))
            assert "--hash=sha256:" in block
            assert "file:" not in block
            assert "git+" not in block
            assert " @ " not in block
        names_by_target[path.stem] = names
        assert "spectrocube" not in names
        assert "sif-parser" not in names
        assert "echelle-spectra" not in names
    expected = {
        "lmfit",
        "matplotlib",
        "numpy",
        "pandas",
        "pyqt5",
        "pyqtgraph",
        "scipy",
        "peakutils",
        "plotly",
        "nbformat",
        "netcdf4",
        "xarray",
    }
    for names in names_by_target.values():
        assert expected <= names
        assert "tzdata" in names


def test_lock_is_limited_to_the_declared_kit_runtime_and_platforms() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    environments = pyproject["tool"]["uv"]["environments"]
    required = pyproject["tool"]["uv"]["required-environments"]
    assert all("python_full_version == '3.12.13'" in marker for marker in environments)
    assert {"win32", "darwin"} == {
        "win32" if "win32" in marker else "darwin" for marker in environments
    }
    assert any("AMD64" in marker for marker in required)
    assert any("arm64" in marker for marker in required)
    assert any("x86_64" in marker for marker in required)
    lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))
    assert lock["requires-python"] == ">=3.9"
    locked_names = {package["name"] for package in lock["package"]}
    assert {"echelle-spectra", "spectrocube", "pyqt5", "netcdf4"} <= locked_names
    spectrocube = next(package for package in lock["package"] if package["name"] == "spectrocube")
    assert spectrocube["version"] == "0.2.0"
    assert spectrocube["source"]["git"].endswith(
        "0b02ac96e386c7121ca2c30f6d36a76518e4e83a"
    )


def test_manifest_refuses_schema_url_and_checksum_drift(tmp_path: Path) -> None:
    original = MANIFEST.read_text(encoding="utf-8")
    cases = (
        original.replace("schema = 1", "schema = 2", 1),
        original.replace("https://github.com/astral-sh/uv", "file:///tmp/uv", 1),
        original.replace(
            "410c2fd3126ff621c9450a21cfc200002c7540dc48d130069a8f619cdb0a811b",
            "not-a-digest",
            1,
        ),
    )
    for index, text in enumerate(cases):
        path = tmp_path / f"bad-{index}.toml"
        path.write_text(text, encoding="utf-8")
        with pytest.raises(KitError):
            load_manifest(path)


def test_fetch_asset_offline_refuses_missing_and_corrupt_cache(tmp_path: Path) -> None:
    payload = b"verified input"
    asset = Asset(
        filename="input.bin",
        url="https://example.invalid/input.bin",
        size=len(payload),
        sha256=hashlib.sha256(payload).hexdigest(),
    )
    with pytest.raises(KitError, match="offline asset cache is incomplete"):
        fetch_asset(asset, tmp_path, offline=True)
    (tmp_path / asset.filename).write_bytes(b"corrupt")
    with pytest.raises(KitError, match="size mismatch"):
        fetch_asset(asset, tmp_path, offline=True)


def test_extract_uv_selects_only_the_named_binary(tmp_path: Path) -> None:
    archive = tmp_path / "uv.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("uv-target/uv.exe", b"binary")
        bundle.writestr("uv-target/uvx.exe", b"not selected")
    destination = tmp_path / "path with spaces" / "uv.exe"
    extract_uv(archive, "uv.exe", destination)
    assert destination.read_bytes() == b"binary"


def test_project_wheel_identity_and_malformed_wheel_refusal(tmp_path: Path) -> None:
    wheel = _wheel(tmp_path / "echelle_spectra-1.6.0-py3-none-any.whl", "echelle_spectra", "1.6.0")
    assert inspect_project_wheel(wheel) == ("echelle_spectra", "1.6.0")
    malformed = tmp_path / "malformed.whl"
    malformed.write_bytes(b"not a zip")
    with pytest.raises(KitError, match="cannot inspect"):
        inspect_project_wheel(malformed)


def test_external_path_guard_and_wheelhouse_refusal(tmp_path: Path) -> None:
    with pytest.raises(KitError, match="outside the repository"):
        assert_external_path(ROOT / "generated", ROOT, "kit destination")
    assert assert_external_path(tmp_path / "path with spaces", ROOT, "kit destination")
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.mkdir()
    (source / "readme.txt").write_text("not a wheel", encoding="utf-8")
    with pytest.raises(KitError, match="non-wheel"):
        copy_wheelhouse(source, destination)


def test_posix_launchers_are_copied_with_lf_shebangs(tmp_path: Path) -> None:
    source = tmp_path / "source.sh"
    destination = tmp_path / "installed.sh"
    source.write_bytes(b"#!/bin/sh\r\nset -eu\r\nprintf 'ok\\n'\r\n")

    copy_posix_script(source, destination)

    assert destination.read_bytes() == b"#!/bin/sh\nset -eu\nprintf 'ok\\n'\n"

    source.write_bytes(b"#!/bin/sh\nset -eu\rprintf bad\n")
    with pytest.raises(KitError, match="bare carriage return"):
        copy_posix_script(source, destination)


def test_wheelhouse_validation_scratch_is_removed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    manifest = load_manifest(MANIFEST)
    platform = manifest.platforms["windows-x86_64"]
    scratch = tmp_path / "scratch"

    def fake_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        target = scratch / "wheelhouse-validation-windows-x86_64"
        target.mkdir(parents=True)
        (target / "artifact.whl").write_bytes(b"temporary")
        return subprocess.CompletedProcess([], 0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    validate_wheelhouse(
        uv="uv",
        manifest=manifest,
        platform=platform,
        requirements=tmp_path / "requirements.txt",
        wheelhouse=tmp_path / "wheelhouse",
        scratch=scratch,
        cache=tmp_path / "cache",
    )
    assert not (scratch / "wheelhouse-validation-windows-x86_64").exists()


def test_checksum_inventory_detects_missing_mismatch_extra_and_traversal(tmp_path: Path) -> None:
    root = tmp_path / "kit with spaces"
    root.mkdir()
    payload = root / "payload.txt"
    payload.write_text("good", encoding="utf-8")
    write_checksums(root)
    verify_checksums(root)
    payload.write_text("changed", encoding="utf-8")
    with pytest.raises(KitError, match="checksum mismatch"):
        verify_checksums(root)
    payload.write_text("good", encoding="utf-8")
    (root / "extra.txt").write_text("extra", encoding="utf-8")
    with pytest.raises(KitError, match="unchecksummed"):
        verify_checksums(root)
    (root / "extra.txt").unlink()
    (root / "checksums.sha256").write_text(f"{'0' * 64}  ../escape\n", encoding="utf-8")
    with pytest.raises(KitError, match="unsafe checksum path"):
        verify_checksums(root)


def test_three_command_readme_and_installers_expose_offline_refusals() -> None:
    readme = (ROOT / "README-KIT.md").read_text(encoding="utf-8")
    operator_guide = (ROOT / "docs" / "operator-cheat-sheet.md").read_text(
        encoding="utf-8"
    )
    windows = re.search(r"## Windows x86-64: three commands.*?```powershell\n(.*?)```", readme, re.S)
    macos = re.search(r"## macOS: three commands.*?```bash\n(.*?)```", readme, re.S)
    assert windows is not None and len(windows.group(1).strip().splitlines()) == 3
    assert macos is not None and len(macos.group(1).strip().splitlines()) == 3
    assert "OPERATOR-CHEAT-SHEET.md" in readme
    assert "Lab is optional development convenience" in operator_guide
    assert "Portable kit on Windows" in operator_guide
    assert "Portable kit on macOS" in operator_guide
    powershell = (ROOT / "kit" / "install.ps1").read_text(encoding="utf-8")
    shell = (ROOT / "kit" / "install.sh").read_text(encoding="utf-8")
    for text in (powershell, shell):
        assert "offline" in text
        assert "online" in text
        assert "no-index" in text
        assert "no-python-downloads" in text
        assert "checksum" in text.lower()
        assert "3.12.13" in text
        assert "UV_CACHE_DIR" in text
    assert 'machine_os=$(uname -s)' in shell
    assert '"$machine_os" != "Darwin"' in shell
    assert (
        'python_version=$("$python" -c '
        "'import platform; print(platform.python_version())')"
    ) in shell
    assert re.search(r"\$\(\$[A-Za-z_]", shell) is None
    for relative in ("kit/install.sh", "kit/echelle"):
        payload = (ROOT / relative).read_bytes()
        assert payload.startswith(b"#!/bin/sh\n")
        assert b"\r" not in payload
    attributes = (ROOT / ".gitattributes").read_text(encoding="utf-8")
    assert "kit/install.sh text eol=lf" in attributes
    assert "kit/echelle text eol=lf" in attributes


def test_all_existing_console_entry_points_remain_declared() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert set(pyproject["project"]["scripts"]) == {
        "echelle",
        "echelle_spectra",
        "echelle-calib",
        "echelle-snapshot",
        "echelle-pattern",
        "echelle-align",
        "echelle-validate-lines",
        "echelle-wavelength-qc",
        "echelle-nist-overlay",
        "echelle-spectrocube",
    }
