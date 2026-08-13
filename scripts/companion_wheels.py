"""Build the two manifest-pinned universal companion wheels for a NIFS kit."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

from scripts.nifs_kit import (
    Asset,
    KitError,
    assert_external_path,
    inspect_project_wheel,
    load_manifest,
    sha256_file,
    verify_asset,
)
from scripts.reproducible_build import _check_tool_versions, normalize_wheel


def _clean_source(path: Path, label: str) -> Path:
    source = path.resolve()
    if not any((source / name).is_file() for name in ("pyproject.toml", "setup.py")):
        raise KitError(f"{label} source is missing pyproject.toml or setup.py: {source}")
    if (source / ".git").exists():
        raise KitError(f"{label} source must be an exported copy without .git")
    return source


def _build_wheel(
    *, source: Path, destination: Path, python: str, epoch: int, expected: tuple[str, str]
) -> Path:
    raw = destination / "raw"
    raw.mkdir(parents=True)
    env = os.environ.copy()
    env.update({"PYTHONHASHSEED": "0", "SOURCE_DATE_EPOCH": str(epoch), "TZ": "UTC"})
    try:
        subprocess.run(
            [python, "-m", "build", "--no-isolation", "--wheel", "--outdir", str(raw)],
            cwd=source,
            env=env,
            check=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise KitError(f"companion wheel build failed for {expected[0]}: {exc}") from exc
    wheels = list(raw.glob("*.whl"))
    if len(wheels) != 1:
        raise KitError(f"expected one {expected[0]} wheel; found {len(wheels)}")
    identity = inspect_project_wheel(wheels[0])
    if identity != expected:
        raise KitError(f"companion identity mismatch: expected {expected}, got {identity}")
    output = destination.parent / wheels[0].name
    normalize_wheel(wheels[0], output, epoch)
    return output


def build_companion_wheels(
    *,
    manifest_path: Path,
    spectrocube_source: Path,
    sif_parser_sdist: Path,
    destination: Path,
    python: str,
) -> dict[str, str]:
    """Build normalized companion wheels after verifying both immutable inputs."""

    manifest = load_manifest(manifest_path)
    repository = manifest.path.parent.parent
    destination = assert_external_path(destination, repository, "companion destination")
    if destination.exists() and any(destination.iterdir()):
        raise KitError(f"companion destination is not empty: {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    _check_tool_versions(manifest.toolchain)
    specifications = {item.name: item for item in manifest.companions}
    sif = specifications["sif_parser"]
    verify_asset(
        sif_parser_sdist,
        Asset(
            filename=sif_parser_sdist.name,
            url=sif.source_url,
            size=int(sif.source_size or 0),
            sha256=sif.source_ref.removeprefix("sha256:"),
        ),
    )
    epoch = int(manifest.toolchain["source_date_epoch"])
    with tempfile.TemporaryDirectory(prefix="echelle-companions-", dir=destination.parent) as temp:
        scratch = Path(temp)
        spectrocube = specifications["spectrocube"]
        first_scratch = _build_wheel(
            source=_clean_source(spectrocube_source, "SpectroCube"),
            destination=scratch / "spectrocube",
            python=python,
            epoch=epoch,
            expected=(spectrocube.name, spectrocube.version),
        )
        with tarfile.open(sif_parser_sdist, "r:gz") as archive:
            archive.extractall(scratch / "sif-parser-source", filter="data")
        roots = [
            path
            for path in (scratch / "sif-parser-source").iterdir()
            if path.is_dir()
            and any((path / name).is_file() for name in ("pyproject.toml", "setup.py"))
        ]
        if len(roots) != 1:
            raise KitError(f"expected one sif_parser source root; found {len(roots)}")
        second_scratch = _build_wheel(
            source=_clean_source(roots[0], "sif_parser"),
            destination=scratch / "sif-parser",
            python=python,
            epoch=epoch,
            expected=(sif.name, sif.version),
        )
        outputs = [destination / first_scratch.name, destination / second_scratch.name]
        for source, output in zip((first_scratch, second_scratch), outputs, strict=True):
            shutil.copy2(source, output)
    digests = {path.name: sha256_file(path) for path in outputs}
    (destination / "SHA256SUMS").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sorted(digests.items())),
        encoding="utf-8",
        newline="\n",
    )
    return digests


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("kit/manifest.toml"))
    parser.add_argument("--spectrocube-source", required=True, type=Path)
    parser.add_argument("--sif-parser-sdist", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--python", default=sys.executable)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        digests = build_companion_wheels(
            manifest_path=args.manifest,
            spectrocube_source=args.spectrocube_source,
            sif_parser_sdist=args.sif_parser_sdist,
            destination=args.destination,
            python=args.python,
        )
    except KitError as exc:
        print(f"companion build error: {exc}", file=sys.stderr)
        return 2
    for name, digest in sorted(digests.items()):
        print(f"{digest}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
