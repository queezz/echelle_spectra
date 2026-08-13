"""Build and normalize wheel/sdist archives under a checked byte contract."""

from __future__ import annotations

import argparse
import copy
import gzip
import importlib.metadata
import os
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile
from pathlib import Path

from scripts.nifs_kit import KitError, load_manifest, sha256_file


def _zip_timestamp(epoch: int) -> tuple[int, int, int, int, int, int]:
    value = time.gmtime(max(epoch, 315532800))
    return value.tm_year, value.tm_mon, value.tm_mday, value.tm_hour, value.tm_min, value.tm_sec


def normalize_wheel(source: Path, destination: Path, epoch: int) -> None:
    """Rewrite ZIP metadata, ordering, permissions, and compression deterministically."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    timestamp = _zip_timestamp(epoch)
    with zipfile.ZipFile(source) as incoming, zipfile.ZipFile(
        destination,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
        strict_timestamps=True,
    ) as outgoing:
        for original in sorted(incoming.infolist(), key=lambda item: item.filename):
            info = zipfile.ZipInfo(original.filename, date_time=timestamp)
            info.create_system = 3
            info.compress_type = zipfile.ZIP_DEFLATED
            info.flag_bits = original.flag_bits & 0x800
            mode = 0o755 if original.is_dir() else 0o644
            info.external_attr = (mode & 0xFFFF) << 16
            info.comment = b""
            info.extra = b""
            outgoing.writestr(info, incoming.read(original.filename))


def normalize_sdist(source: Path, destination: Path, epoch: int) -> None:
    """Rewrite gzip and tar metadata under one stable POSIX archive contract."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(source, "r:gz") as incoming, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=epoch, compresslevel=9) as zipped:
            with tarfile.open(fileobj=zipped, mode="w", format=tarfile.PAX_FORMAT) as outgoing:
                for original in sorted(incoming.getmembers(), key=lambda item: item.name):
                    member = copy.copy(original)
                    member.uid = 0
                    member.gid = 0
                    member.uname = ""
                    member.gname = ""
                    member.mtime = epoch
                    member.pax_headers = {}
                    if member.isdir():
                        member.mode = 0o755
                    elif member.isfile():
                        member.mode = 0o644
                    content = incoming.extractfile(original) if original.isfile() else None
                    outgoing.addfile(member, content)


def _check_tool_versions(toolchain: dict[str, object]) -> None:
    for distribution in ("build", "setuptools", "wheel"):
        expected = str(toolchain[distribution])
        try:
            actual = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError as exc:
            raise KitError(f"release tool is not installed: {distribution}=={expected}") from exc
        if actual != expected:
            raise KitError(
                f"release tool mismatch: expected {distribution}=={expected}, got {actual}"
            )


def _validate_clean_source(source: Path) -> Path:
    source = source.resolve()
    required = (source / "pyproject.toml", source / "README.md", source / "src")
    missing = [str(path.name) for path in required if not path.exists()]
    if missing:
        raise KitError(f"source copy is incomplete; missing: {', '.join(missing)}")
    if (source / ".git").exists():
        raise KitError(
            "source must be an exported clean copy without .git, not the live checkout"
        )
    return source


def _validate_destination(source: Path, destination: Path) -> Path:
    destination = destination.resolve()
    if destination == source or source in destination.parents:
        raise KitError("distribution destination must be outside the clean source copy")
    if destination.exists() and any(destination.iterdir()):
        raise KitError(f"distribution destination is not empty: {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def build_distributions(
    *, source: Path, destination: Path, manifest_path: Path, python: str
) -> dict[str, str]:
    """Build once and emit normalized byte-stable wheel and sdist archives."""

    source = _validate_clean_source(source)
    destination = _validate_destination(source, destination)
    manifest = load_manifest(manifest_path)
    _check_tool_versions(manifest.toolchain)
    epoch = int(manifest.toolchain["source_date_epoch"])
    env = os.environ.copy()
    env.update(
        {
            "PYTHONHASHSEED": "0",
            "SOURCE_DATE_EPOCH": str(epoch),
            "TZ": "UTC",
        }
    )
    with tempfile.TemporaryDirectory(prefix="echelle-dist-build-", dir=destination.parent) as temp:
        raw = Path(temp) / "raw"
        raw.mkdir()
        try:
            subprocess.run(
                [python, "-m", "build", "--no-isolation", "--wheel", "--sdist", "--outdir", str(raw)],
                cwd=source,
                env=env,
                check=True,
            )
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            raise KitError(f"distribution build failed: {exc}") from exc
        wheels = sorted(raw.glob("*.whl"))
        sdists = sorted(raw.glob("*.tar.gz"))
        if len(wheels) != 1 or len(sdists) != 1:
            raise KitError(
                f"expected one wheel and one sdist; found {len(wheels)} wheel(s), "
                f"{len(sdists)} sdist(s)"
            )
        normalize_wheel(wheels[0], destination / wheels[0].name, epoch)
        normalize_sdist(sdists[0], destination / sdists[0].name, epoch)
    artifacts = sorted(path for path in destination.iterdir() if path.suffix in {".whl", ".gz"})
    digests = {path.name: sha256_file(path) for path in artifacts}
    (destination / "SHA256SUMS").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sorted(digests.items())),
        encoding="utf-8",
        newline="\n",
    )
    return digests


def compare_distributions(first: Path, second: Path) -> dict[str, str]:
    """Require identical artifact names and SHA-256 digests across two builds."""

    def inventory(root: Path) -> dict[str, str]:
        if not root.is_dir():
            raise KitError(f"distribution directory does not exist: {root}")
        return {
            path.name: sha256_file(path)
            for path in sorted(root.iterdir())
            if path.suffix in {".whl", ".gz"}
        }

    left = inventory(first)
    right = inventory(second)
    if not left or set(left) != set(right):
        raise KitError(
            f"artifact names differ: first={sorted(left)}, second={sorted(right)}"
        )
    mismatches = [name for name in left if left[name] != right[name]]
    if mismatches:
        details = ", ".join(
            f"{name} ({left[name]} != {right[name]})" for name in mismatches
        )
        raise KitError(f"distribution bytes are not reproducible: {details}")
    return left


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="Build normalized wheel and sdist archives.")
    build.add_argument("--source", required=True, type=Path)
    build.add_argument("--destination", required=True, type=Path)
    build.add_argument("--manifest", type=Path, default=Path("kit/manifest.toml"))
    build.add_argument("--python", default=sys.executable)
    compare = commands.add_parser("compare", help="Compare two independent build directories.")
    compare.add_argument("first", type=Path)
    compare.add_argument("second", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "build":
            digests = build_distributions(
                source=args.source,
                destination=args.destination,
                manifest_path=args.manifest,
                python=args.python,
            )
        else:
            digests = compare_distributions(args.first, args.second)
    except KitError as exc:
        print(f"build error: {exc}", file=sys.stderr)
        return 2
    for name, digest in sorted(digests.items()):
        print(f"{digest}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
