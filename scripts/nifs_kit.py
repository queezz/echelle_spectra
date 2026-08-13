"""Assemble and validate a pinned portable NIFS kit outside the checkout.

The generated payload is deliberately not a Python package artifact. It is a
USB-friendly folder containing one target's uv binary, CPython archive,
dependency wheelhouse, project wheel, installers, and a complete checksum
inventory. This module uses only the standard library; pip is invoked through
the pinned host uv when a wheelhouse must be fetched.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.parse
import urllib.request
import uuid
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

import tomllib

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
VERSION_RE = re.compile(r"^[0-9]+(?:\.[0-9]+){1,2}$")
CHECKSUMS_NAME = "checksums.sha256"


class KitError(RuntimeError):
    """An actionable portable-kit validation or assembly failure."""


@dataclass(frozen=True)
class Asset:
    """One immutable upstream input."""

    filename: str
    url: str
    size: int
    sha256: str


@dataclass(frozen=True)
class PlatformSpec:
    """Narrow adapter configuration for one supported target."""

    key: str
    os_name: str
    arch: str
    minimum_os_version: str
    uv_python_platform: str
    pip_platforms: tuple[str, ...]
    uv: Asset
    uv_binary: str
    python: Asset
    python_executable: str
    installer: str
    launcher: str


@dataclass(frozen=True)
class CompanionSpec:
    """A required unpublished companion wheel pinned to immutable source."""

    name: str
    version: str
    wheel_prefix: str
    source_kind: str
    source_url: str
    source_ref: str
    source_size: int | None


@dataclass(frozen=True)
class KitManifest:
    """Validated declarative kit contract."""

    path: Path
    application_name: str
    application_version: str
    wheel_prefix: str
    runtime_extra: str
    companions: tuple[CompanionSpec, ...]
    toolchain: dict[str, Any]
    platforms: dict[str, PlatformSpec]


def _required(mapping: dict[str, Any], key: str, kind: type, context: str) -> Any:
    value = mapping.get(key)
    if not isinstance(value, kind):
        raise KitError(f"{context}.{key} must be {kind.__name__}")
    return value


def _validate_asset(
    table: dict[str, Any], prefix: str, context: str
) -> Asset:
    filename = _required(table, f"{prefix}_archive", str, context)
    url = _required(table, f"{prefix}_url", str, context)
    size = _required(table, f"{prefix}_size", int, context)
    digest = _required(table, f"{prefix}_sha256", str, context)
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or not parsed.netloc:
        raise KitError(f"{context}.{prefix}_url must be an absolute HTTPS URL")
    if PurePosixPath(filename).name != filename or "\\" in filename:
        raise KitError(f"{context}.{prefix}_archive must be a plain filename")
    if urllib.parse.unquote(PurePosixPath(parsed.path).name) != filename:
        raise KitError(f"{context}.{prefix}_url does not name {filename}")
    if size <= 0:
        raise KitError(f"{context}.{prefix}_size must be positive")
    if not SHA256_RE.fullmatch(digest):
        raise KitError(f"{context}.{prefix}_sha256 must be 64 lowercase hex characters")
    return Asset(filename=filename, url=url, size=size, sha256=digest)


def load_manifest(path: Path) -> KitManifest:  # noqa: C901
    """Load and strictly validate the authoritative TOML manifest."""

    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise KitError(f"cannot read kit manifest {path}: {exc}") from exc
    if data.get("schema") != 1:
        raise KitError("kit manifest schema must be 1")
    application = _required(data, "application", dict, "manifest")
    toolchain = _required(data, "toolchain", dict, "manifest")
    platform_tables = _required(data, "platform", dict, "manifest")
    companion_tables = _required(data, "companion", list, "manifest")
    name = _required(application, "name", str, "application")
    version = _required(application, "version", str, "application")
    wheel_prefix = _required(application, "wheel_prefix", str, "application")
    runtime_extra = _required(application, "runtime_extra", str, "application")
    if name != "echelle_spectra":
        raise KitError("application.name must be echelle_spectra")
    if not VERSION_RE.fullmatch(version):
        raise KitError("application.version must be an exact X.Y.Z version")
    if wheel_prefix != "echelle_spectra-":
        raise KitError("application.wheel_prefix must be echelle_spectra-")
    companions: list[CompanionSpec] = []
    for index, value in enumerate(companion_tables):
        context = f"companion[{index}]"
        if not isinstance(value, dict):
            raise KitError(f"{context} must be a table")
        companion = CompanionSpec(
            name=_required(value, "name", str, context),
            version=_required(value, "version", str, context),
            wheel_prefix=_required(value, "wheel_prefix", str, context),
            source_kind=_required(value, "source_kind", str, context),
            source_url=_required(value, "source_url", str, context),
            source_ref=_required(value, "source_ref", str, context),
            source_size=value.get("source_size"),
        )
        expected_companions = {"spectrocube": "0.1.0", "sif_parser": "0.3.6"}
        if expected_companions.get(companion.name) != companion.version:
            raise KitError("companions must be spectrocube 0.1.0 and sif_parser 0.3.6")
        if companion.wheel_prefix != f"{companion.name}-":
            raise KitError(f"{context}.wheel_prefix must be {companion.name}-")
        parsed_source = urllib.parse.urlparse(companion.source_url)
        if parsed_source.scheme != "https" or not parsed_source.netloc:
            raise KitError(f"{context}.source_url must be an absolute HTTPS URL")
        if companion.source_kind == "git":
            if not re.fullmatch(r"[0-9a-f]{40}", companion.source_ref):
                raise KitError(f"{context}.source_ref must be a full lowercase Git commit")
            if companion.source_size is not None:
                raise KitError(f"{context}.source_size is only valid for an sdist")
        elif companion.source_kind == "sdist":
            if not re.fullmatch(r"sha256:[0-9a-f]{64}", companion.source_ref):
                raise KitError(f"{context}.source_ref must be a sha256 digest")
            if not isinstance(companion.source_size, int) or companion.source_size <= 0:
                raise KitError(f"{context}.source_size must be positive for an sdist")
        else:
            raise KitError(f"{context}.source_kind must be git or sdist")
        companions.append(companion)
    if {item.name for item in companions} != {"spectrocube", "sif_parser"}:
        raise KitError("manifest must contain exactly the two required companion wheels")
    for key in (
        "python",
        "uv",
        "pip",
        "build",
        "setuptools",
        "wheel",
        "twine",
    ):
        value = _required(toolchain, key, str, "toolchain")
        if not VERSION_RE.fullmatch(value):
            raise KitError(f"toolchain.{key} must be an exact X.Y.Z version")
    source_epoch = _required(toolchain, "source_date_epoch", int, "toolchain")
    if source_epoch < 315532800:
        raise KitError("toolchain.source_date_epoch must be at or after 1980-01-01")

    expected = {"windows-x86_64", "macos-aarch64", "macos-x86_64"}
    if set(platform_tables) != expected:
        raise KitError(f"platform keys must be exactly: {', '.join(sorted(expected))}")
    platforms: dict[str, PlatformSpec] = {}
    for key, value in platform_tables.items():
        if not isinstance(value, dict):
            raise KitError(f"platform.{key} must be a table")
        context = f"platform.{key}"
        pip_platforms = _required(value, "pip_platforms", list, context)
        if not pip_platforms or not all(isinstance(item, str) and item for item in pip_platforms):
            raise KitError(f"{context}.pip_platforms must be a non-empty string list")
        platforms[key] = PlatformSpec(
            key=key,
            os_name=_required(value, "os", str, context),
            arch=_required(value, "arch", str, context),
            minimum_os_version=_required(value, "minimum_os_version", str, context),
            uv_python_platform=_required(value, "uv_python_platform", str, context),
            pip_platforms=tuple(pip_platforms),
            uv=_validate_asset(value, "uv", context),
            uv_binary=_required(value, "uv_binary", str, context),
            python=_validate_asset(value, "python", context),
            python_executable=_required(value, "python_executable", str, context),
            installer=_required(value, "installer", str, context),
            launcher=_required(value, "launcher", str, context),
        )
        if not re.fullmatch(r"[0-9]+(?:\.[0-9]+){1,2}", platforms[key].minimum_os_version):
            raise KitError(f"{context}.minimum_os_version must be numeric")
    return KitManifest(
        path=path.resolve(),
        application_name=name,
        application_version=version,
        wheel_prefix=wheel_prefix,
        runtime_extra=runtime_extra,
        companions=tuple(companions),
        toolchain=dict(toolchain),
        platforms=platforms,
    )


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_asset(path: Path, asset: Asset) -> None:
    """Refuse a missing, wrong-sized, or checksum-mismatched asset."""

    if not path.is_file():
        raise KitError(f"required asset is missing: {path}")
    actual_size = path.stat().st_size
    if actual_size != asset.size:
        raise KitError(
            f"asset size mismatch for {path.name}: expected {asset.size}, got {actual_size}"
        )
    actual_digest = sha256_file(path)
    if actual_digest != asset.sha256:
        raise KitError(
            f"asset checksum mismatch for {path.name}: "
            f"expected {asset.sha256}, got {actual_digest}"
        )


def fetch_asset(asset: Asset, cache: Path, *, offline: bool) -> Path:
    """Fetch once into an explicit external cache, verifying before every use."""

    cache.mkdir(parents=True, exist_ok=True)
    target = cache / asset.filename
    if target.exists():
        verify_asset(target, asset)
        return target
    if offline:
        raise KitError(f"offline asset cache is incomplete: {target.name} is missing")
    temporary = cache / f".{asset.filename}.partial-{uuid.uuid4().hex}"
    try:
        with urllib.request.urlopen(asset.url, timeout=120) as response, temporary.open("wb") as out:
            shutil.copyfileobj(response, out, length=1024 * 1024)
        verify_asset(temporary, asset)
        os.replace(temporary, target)
    except Exception as exc:
        temporary.unlink(missing_ok=True)
        if isinstance(exc, KitError):
            raise
        raise KitError(f"failed to download {asset.url}: {exc}") from exc
    return target


def extract_uv(archive: Path, binary_name: str, destination: Path) -> None:
    """Extract only the uv executable, never archive paths."""

    matches: list[tuple[str, bytes]] = []
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as bundle:
            for name in bundle.namelist():
                if PurePosixPath(name).name == binary_name and not name.endswith("/"):
                    matches.append((name, bundle.read(name)))
    else:
        try:
            with tarfile.open(archive, "r:*") as bundle:
                for member in bundle.getmembers():
                    if PurePosixPath(member.name).name != binary_name or not member.isfile():
                        continue
                    handle = bundle.extractfile(member)
                    if handle is None:
                        raise KitError(f"cannot read {member.name} from {archive.name}")
                    matches.append((member.name, handle.read()))
        except tarfile.TarError as exc:
            raise KitError(f"unsupported uv archive {archive.name}: {exc}") from exc
    if len(matches) != 1:
        names = ", ".join(name for name, _ in matches) or "none"
        raise KitError(
            f"expected one {binary_name} in {archive.name}; found {len(matches)} ({names})"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(matches[0][1])
    destination.chmod(0o755)


def inspect_project_wheel(path: Path) -> tuple[str, str]:
    """Return normalized Name and Version from a wheel's METADATA."""

    if not path.is_file() or path.suffix != ".whl":
        raise KitError(f"project artifact must be one wheel file: {path}")
    try:
        with zipfile.ZipFile(path) as bundle:
            metadata_names = [
                name for name in bundle.namelist() if name.endswith(".dist-info/METADATA")
            ]
            if len(metadata_names) != 1:
                raise KitError(
                    f"wheel {path.name} must contain exactly one .dist-info/METADATA"
                )
            metadata = bundle.read(metadata_names[0]).decode("utf-8")
    except (OSError, UnicodeDecodeError, zipfile.BadZipFile) as exc:
        raise KitError(f"cannot inspect project wheel {path}: {exc}") from exc
    fields: dict[str, str] = {}
    for line in metadata.splitlines():
        if ": " in line:
            key, value = line.split(": ", 1)
            fields.setdefault(key, value)
    name = re.sub(r"[-_.]+", "_", fields.get("Name", "")).lower()
    version = fields.get("Version", "")
    if not name or not version:
        raise KitError(f"wheel {path.name} metadata is missing Name or Version")
    return name, version


def assert_external_path(path: Path, repository: Path, label: str) -> Path:
    """Refuse generated payload/cache paths inside the synced repository."""

    resolved = path.resolve()
    root = repository.resolve()
    if resolved == root or root in resolved.parents:
        raise KitError(f"{label} must be outside the repository: {resolved}")
    return resolved


def _run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    try:
        subprocess.run(command, cwd=cwd, env=env, check=True)
    except FileNotFoundError as exc:
        raise KitError(f"required command is unavailable: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise KitError(f"command failed with exit code {exc.returncode}: {' '.join(command)}") from exc


def check_uv_version(uv: str, expected: str, *, cwd: Path) -> None:
    """Prevent a silently different resolver/downloader from assembling the kit."""

    try:
        result = subprocess.run(
            [uv, "--version"], cwd=cwd, check=True, capture_output=True, text=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise KitError(f"cannot execute host uv at {uv}") from exc
    match = re.search(r"\buv\s+([0-9]+(?:\.[0-9]+){2})\b", result.stdout)
    if match is None or match.group(1) != expected:
        actual = result.stdout.strip() or "unknown"
        raise KitError(f"host uv must be {expected}; got {actual}")


def download_wheelhouse(
    *,
    uv: str,
    manifest: KitManifest,
    platform: PlatformSpec,
    requirements: Path,
    destination: Path,
    cache: Path,
) -> None:
    """Download only compatible binary wheels with pip's hash enforcement."""

    command = [
        uv,
        "tool",
        "run",
        "--from",
        f"pip=={manifest.toolchain['pip']}",
        "pip",
        "download",
        "--require-hashes",
        "--only-binary=:all:",
        "--implementation",
        "cp",
        "--python-version",
        manifest.toolchain["python"],
        "--abi",
        manifest.toolchain["python_tag"],
        "--abi",
        "abi3",
        "--abi",
        "none",
        "--dest",
        str(destination),
        "--requirement",
        str(requirements),
    ]
    for tag in platform.pip_platforms:
        command.extend(["--platform", tag])
    env = os.environ.copy()
    env["UV_CACHE_DIR"] = str(cache)
    _run(command, cwd=manifest.path.parent.parent, env=env)


def copy_wheelhouse(source: Path, destination: Path) -> None:
    """Reuse a previously fetched wheelhouse without accepting other file types."""

    if not source.is_dir():
        raise KitError(f"wheelhouse source is not a directory: {source}")
    wheels = sorted(source.glob("*.whl"))
    unexpected = sorted(path.name for path in source.iterdir() if path.is_file() and path.suffix != ".whl")
    if unexpected:
        raise KitError(f"wheelhouse source contains non-wheel files: {', '.join(unexpected)}")
    if not wheels:
        raise KitError(f"wheelhouse source contains no wheels: {source}")
    destination.mkdir(parents=True, exist_ok=True)
    for wheel in wheels:
        shutil.copy2(wheel, destination / wheel.name)


def validate_wheelhouse(
    *,
    uv: str,
    manifest: KitManifest,
    platform: PlatformSpec,
    requirements: Path,
    wheelhouse: Path,
    scratch: Path,
    cache: Path,
) -> None:
    """Prove pinned target wheels resolve from the wheelhouse without a network."""

    target = scratch / f"wheelhouse-validation-{platform.key}"
    if target.exists():
        shutil.rmtree(target)
    command = [
        uv,
        "tool",
        "run",
        "--offline",
        "--from",
        f"pip=={manifest.toolchain['pip']}",
        "pip",
        "download",
        "--dest",
        str(target),
        "--implementation",
        "cp",
        "--python-version",
        manifest.toolchain["python"],
        "--abi",
        manifest.toolchain["python_tag"],
        "--abi",
        "abi3",
        "--abi",
        "none",
        "--only-binary=:all:",
        "--require-hashes",
        "--no-index",
        "--find-links",
        str(wheelhouse),
        "--requirement",
        str(requirements),
    ]
    for tag in platform.pip_platforms:
        command.extend(["--platform", tag])
    env = os.environ.copy()
    env["UV_CACHE_DIR"] = str(cache)
    try:
        _run(command, cwd=manifest.path.parent.parent, env=env)
    finally:
        if target.exists():
            shutil.rmtree(target)


def write_checksums(root: Path) -> Path:
    """Write a sorted, relative, machine-path-free checksum inventory."""

    checksum_path = root / CHECKSUMS_NAME
    files = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.name != CHECKSUMS_NAME
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(root).as_posix()}" for path in files]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return checksum_path


def verify_checksums(root: Path, *, refuse_extras: bool = True) -> None:  # noqa: C901
    """Verify the generated payload inventory and optionally reject extra files."""

    checksum_path = root / CHECKSUMS_NAME
    if not checksum_path.is_file():
        raise KitError(f"kit checksum inventory is missing: {checksum_path}")
    listed: set[Path] = set()
    for number, line in enumerate(checksum_path.read_text(encoding="utf-8").splitlines(), 1):
        if "  " not in line:
            raise KitError(f"invalid checksum line {number}: expected SHA256, two spaces, path")
        digest, relative_text = line.split("  ", 1)
        relative = PurePosixPath(relative_text)
        if not SHA256_RE.fullmatch(digest):
            raise KitError(f"invalid checksum digest on line {number}")
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise KitError(f"unsafe checksum path on line {number}: {relative_text}")
        local = Path(*relative.parts)
        if local in listed:
            raise KitError(f"duplicate checksum path on line {number}: {relative_text}")
        listed.add(local)
        path = root / local
        if not path.is_file():
            raise KitError(f"checksummed payload is missing: {relative_text}")
        actual = sha256_file(path)
        if actual != digest:
            raise KitError(
                f"payload checksum mismatch for {relative_text}: expected {digest}, got {actual}"
            )
    if refuse_extras:
        actual_files = {
            path.relative_to(root)
            for path in root.rglob("*")
            if path.is_file() and path.name != CHECKSUMS_NAME
        }
        extras = sorted(path.as_posix() for path in actual_files - listed)
        if extras:
            raise KitError(f"kit contains unchecksummed files: {', '.join(extras)}")


def assemble_kit(  # noqa: C901
    *,
    manifest: KitManifest,
    platform_key: str,
    destination: Path,
    cache: Path,
    project_wheel: Path,
    companion_wheels: tuple[Path, ...],
    host_uv: str,
    wheelhouse_source: Path | None,
    offline_assets: bool,
) -> Path:
    """Assemble one target atomically into an explicit external destination."""

    repository = manifest.path.parent.parent
    destination = assert_external_path(destination, repository, "kit destination")
    cache = assert_external_path(cache, repository, "asset/cache directory")
    if platform_key not in manifest.platforms:
        raise KitError(
            f"unsupported platform {platform_key!r}; choose: "
            f"{', '.join(sorted(manifest.platforms))}"
        )
    if destination.exists():
        try:
            verify_checksums(destination)
        except KitError as exc:
            raise KitError(
                f"destination already exists but is not a complete valid kit: {destination}; {exc}"
            ) from exc
        return destination
    check_uv_version(host_uv, manifest.toolchain["uv"], cwd=repository)
    name, version = inspect_project_wheel(project_wheel)
    if name != manifest.application_name or version != manifest.application_version:
        raise KitError(
            f"project wheel identity mismatch: expected {manifest.application_name} "
            f"{manifest.application_version}, got {name} {version}"
        )
    if len(companion_wheels) != len(manifest.companions):
        raise KitError(
            f"expected {len(manifest.companions)} companion wheel(s); "
            f"got {len(companion_wheels)}"
        )
    companion_inputs: dict[str, tuple[Path, str]] = {}
    for wheel_path in companion_wheels:
        companion_name, companion_version = inspect_project_wheel(wheel_path)
        if companion_name in companion_inputs:
            raise KitError(f"duplicate companion wheel identity: {companion_name}")
        companion_inputs[companion_name] = (wheel_path, companion_version)
    for specification in manifest.companions:
        wheel_input = companion_inputs.get(specification.name)
        if wheel_input is None or wheel_input[1] != specification.version:
            actual = "missing" if wheel_input is None else wheel_input[1]
            raise KitError(
                f"companion wheel identity mismatch: expected {specification.name} "
                f"{specification.version} from {specification.source_ref}, got "
                f"{specification.name} {actual}"
            )
    platform = manifest.platforms[platform_key]
    stage = destination.with_name(f".{destination.name}.partial-{uuid.uuid4().hex}")
    stage.mkdir(parents=True)
    try:
        uv_archive = fetch_asset(platform.uv, cache / "downloads", offline=offline_assets)
        python_archive = fetch_asset(platform.python, cache / "downloads", offline=offline_assets)
        extract_uv(uv_archive, platform.uv_binary, stage / "bin" / platform.uv_binary)
        runtime = stage / "runtime"
        runtime.mkdir()
        shutil.copy2(python_archive, runtime / platform.python.filename)

        requirements_source = (
            manifest.path.parent / "requirements" / f"{platform.key}.txt"
        )
        if not requirements_source.is_file():
            raise KitError(f"runtime requirements lock is missing: {requirements_source}")
        requirements = stage / "requirements-runtime.txt"
        shutil.copy2(requirements_source, requirements)
        wheelhouse = stage / "wheelhouse"
        wheelhouse.mkdir()
        if wheelhouse_source is None:
            if offline_assets:
                raise KitError("offline assembly requires --wheelhouse-source")
            download_wheelhouse(
                uv=host_uv,
                manifest=manifest,
                platform=platform,
                requirements=requirements,
                destination=wheelhouse,
                cache=cache / "uv",
            )
        else:
            source = assert_external_path(wheelhouse_source, repository, "wheelhouse source")
            copy_wheelhouse(source, wheelhouse)
        shutil.copy2(project_wheel, wheelhouse / project_wheel.name)
        for companion_wheel, _ in companion_inputs.values():
            shutil.copy2(companion_wheel, wheelhouse / companion_wheel.name)
        validate_wheelhouse(
            uv=host_uv,
            manifest=manifest,
            platform=platform,
            requirements=requirements,
            wheelhouse=wheelhouse,
            scratch=cache / "validation",
            cache=cache / "uv",
        )

        shutil.copy2(repository / "README-KIT.md", stage / "README-KIT.md")
        shutil.copy2(manifest.path, stage / "kit-manifest.toml")
        shutil.copy2(manifest.path.parent / platform.installer, stage / platform.installer)
        shutil.copy2(manifest.path.parent / platform.launcher, stage / platform.launcher)
        (stage / "platform.txt").write_text(platform.key + "\n", encoding="utf-8", newline="\n")
        write_checksums(stage)
        verify_checksums(stage)
        stage.rename(destination)
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        raise
    return destination


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Assemble one verified portable NIFS kit outside the repository."
    )
    parser.add_argument("--manifest", type=Path, default=Path("kit/manifest.toml"))
    parser.add_argument("--platform", required=True, dest="platform_key")
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--project-wheel", required=True, type=Path)
    parser.add_argument(
        "--companion-wheel",
        required=True,
        action="append",
        type=Path,
        dest="companion_wheels",
        help="Pinned unpublished companion wheel; repeat once per manifest companion.",
    )
    parser.add_argument("--host-uv", default="uv")
    parser.add_argument("--wheelhouse-source", type=Path)
    parser.add_argument(
        "--offline-assets",
        action="store_true",
        help="Forbid downloads and require verified assets plus --wheelhouse-source.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry for repository-local kit assembly."""

    args = _parser().parse_args(argv)
    try:
        result = assemble_kit(
            manifest=load_manifest(args.manifest),
            platform_key=args.platform_key,
            destination=args.destination,
            cache=args.cache,
            project_wheel=args.project_wheel,
            companion_wheels=tuple(args.companion_wheels),
            host_uv=args.host_uv,
            wheelhouse_source=args.wheelhouse_source,
            offline_assets=args.offline_assets,
        )
    except KitError as exc:
        print(f"kit error: {exc}", file=sys.stderr)
        return 2
    print(f"Verified kit: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
