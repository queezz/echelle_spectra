"""Calibration snapshots: immutable, digested campaign artifacts.

A snapshot directory is identified by its folder name and contains one
``snapshot.toml`` binder.  The binder names every calibration input by role and
records a SHA-256 digest, so downstream cubes and run receipts can refer to one
stable snapshot id instead of a loose collection of filenames.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

try:  # pragma: no cover - selected by the running Python version
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib


SNAPSHOT_SCHEMA = "echelle-snapshot/v1"
SNAPSHOT_ID_RE = re.compile(r"^[0-9]{8}_[a-z0-9]+(?:-[a-z0-9]+)?$")

ROLE_FILENAMES = {
    "pattern": "pattern.txt",
    "wavelength": "wavelength.txt",
    "sphere": "sphere.sif",
    "sphere_background": "sphere_bg.sif",
    "integral": "integral.txt",
}
REQUIRED_ROLES = frozenset(ROLE_FILENAMES)


class SnapshotError(ValueError):
    """Base error for an invalid snapshot request or artifact."""


class SnapshotValidationError(SnapshotError):
    """A snapshot failed one or more integrity checks."""

    def __init__(self, errors: Sequence[str]):
        self.errors = tuple(errors)
        super().__init__("; ".join(self.errors))


@dataclass(frozen=True)
class Artifact:
    """One file recorded in a snapshot manifest."""

    role: str
    path: str
    sha256: str
    size_bytes: int
    source_name: str = ""
    label: str = ""


@dataclass(frozen=True)
class Snapshot:
    """A parsed and validated snapshot."""

    root: Path
    manifest: Mapping[str, Any]
    artifacts: tuple[Artifact, ...]

    @property
    def snapshot_id(self) -> str:
        return str(self.manifest["id"])

    @property
    def detector(self) -> str:
        return str(self.manifest["detector"])

    @property
    def lamps(self) -> tuple[str, ...]:
        return tuple(str(item) for item in self.manifest.get("lamps", []))

    def artifact_for_role(self, role: str) -> Artifact:
        """Return the unique required artifact for *role*."""

        matches = [artifact for artifact in self.artifacts if artifact.role == role]
        if len(matches) != 1:
            raise SnapshotValidationError(
                [f"snapshot {self.snapshot_id!r} does not contain exactly one {role!r} artifact"]
            )
        return matches[0]

    def calibration_files(self) -> dict[str, str]:
        """Return the established ``Calibrations`` filename vocabulary."""

        return {
            "orders": self.artifact_for_role("pattern").path,
            "wavelength": self.artifact_for_role("wavelength").path,
            "sphr": self.artifact_for_role("sphere").path,
            "bkgr": self.artifact_for_role("sphere_background").path,
            "integral": self.artifact_for_role("integral").path,
        }

    def provenance_attrs(self) -> dict[str, str]:
        """Return portable, complete NetCDF-string provenance attributes."""

        manifest_sha256, _ = file_digest(self.root / "snapshot.toml")
        artifact_records = [
            {
                "role": artifact.role,
                "path": artifact.path,
                "sha256": artifact.sha256,
                "size_bytes": artifact.size_bytes,
                **({"source_name": artifact.source_name} if artifact.source_name else {}),
                **({"label": artifact.label} if artifact.label else {}),
            }
            for artifact in self.artifacts
        ]
        return {
            "snapshot_id": self.snapshot_id,
            "snapshot_manifest_sha256": manifest_sha256,
            "snapshot_manifest_json": json.dumps(
                self.manifest,
                default=str,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ),
            "calibration_file_digests_json": json.dumps(
                artifact_records,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ),
        }


def validate_snapshot_id(snapshot_id: str) -> str:
    """Return a normalized id or raise with the public naming convention."""

    normalized = str(snapshot_id).strip().lower()
    if not SNAPSHOT_ID_RE.fullmatch(normalized):
        raise SnapshotError(
            "snapshot id must be YYYYMMDD_<detector>[-rev], for example "
            "20250926_cmos or 20250926_cmos-r1"
        )
    return normalized


def file_digest(path: Path) -> tuple[str, int]:
    """Return SHA-256 and byte size for *path*."""

    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _toml_string(value: object) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _toml_array(values: Iterable[object]) -> str:
    return "[" + ", ".join(_toml_string(value) for value in values) + "]"


def render_manifest(  # noqa: C901 - explicit TOML value kinds stay legible here
    *,
    snapshot_id: str,
    detector: str,
    lamps: Sequence[str],
    artifacts: Sequence[Artifact],
    notes: str = "",
    base_snapshot: str | None = None,
    validity: Mapping[str, object] | None = None,
    alignment: Mapping[str, object] | None = None,
    qc: Mapping[str, object] | None = None,
    created_utc: str | None = None,
) -> str:
    """Render a deterministic, hand-editable ``snapshot.toml``."""

    created = created_utc or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    lines = [
        "# Calibration snapshot binder generated by echelle_spectra.",
        "# This is ordinary hand-editable TOML; rerun snapshot validation after edits.",
        f"schema = {_toml_string(SNAPSHOT_SCHEMA)}",
        f"id = {_toml_string(snapshot_id)}",
        f"created_utc = {_toml_string(created)}",
        f"detector = {_toml_string(detector)}",
        f"lamps = {_toml_array(lamps)}",
    ]
    if base_snapshot:
        lines.append(f"base_snapshot = {_toml_string(base_snapshot)}")
    if notes:
        lines.append(f"notes = {_toml_string(notes)}")

    for heading, values in (
        ("validity", validity),
        ("alignment", alignment),
        ("qc", qc),
    ):
        if not values:
            continue
        lines.extend(["", f"[{heading}]"])
        for key, value in values.items():
            if value is None:
                continue
            if isinstance(value, bool):
                rendered = "true" if value else "false"
            elif isinstance(value, (int, float)):
                rendered = str(value)
            elif isinstance(value, (list, tuple)):
                rendered = _toml_array(value)
            else:
                rendered = _toml_string(value)
            lines.append(f"{key} = {rendered}")

    for artifact in artifacts:
        lines.extend(
            [
                "",
                "[[artifacts]]",
                f"role = {_toml_string(artifact.role)}",
                f"path = {_toml_string(artifact.path)}",
                f"sha256 = {_toml_string(artifact.sha256)}",
                f"size_bytes = {artifact.size_bytes}",
            ]
        )
        if artifact.source_name:
            lines.append(f"source_name = {_toml_string(artifact.source_name)}")
        if artifact.label:
            lines.append(f"label = {_toml_string(artifact.label)}")
    return "\n".join(lines) + "\n"


def _safe_relative_path(raw: object) -> Path | None:
    text = str(raw)
    posix = PurePosixPath(text)
    windows = PureWindowsPath(text)
    if not text or posix.is_absolute() or windows.is_absolute() or windows.drive:
        return None
    if ".." in posix.parts or ".." in windows.parts:
        return None
    return Path(*posix.parts)


def _parse_artifact(item: object, index: int, errors: list[str]) -> Artifact | None:
    if not isinstance(item, dict):
        errors.append(f"artifact {index} must be a TOML table")
        return None
    role = str(item.get("role", "")).strip()
    path = str(item.get("path", "")).strip()
    sha256 = str(item.get("sha256", "")).strip().lower()
    try:
        size_bytes = int(item.get("size_bytes", -1))
    except (TypeError, ValueError):
        size_bytes = -1
    if not role:
        errors.append(f"artifact {index} has no role")
    if _safe_relative_path(path) is None:
        errors.append(f"artifact {index} path must stay inside the snapshot: {path!r}")
    if not re.fullmatch(r"[0-9a-f]{64}", sha256):
        errors.append(f"artifact {index} has an invalid SHA-256 digest")
    if size_bytes < 0:
        errors.append(f"artifact {index} has an invalid size_bytes")
    if not role or _safe_relative_path(path) is None or size_bytes < 0:
        return None
    return Artifact(
        role=role,
        path=path,
        sha256=sha256,
        size_bytes=size_bytes,
        source_name=str(item.get("source_name", "")),
        label=str(item.get("label", "")),
    )


def load_snapshot(  # noqa: C901 - one pass accumulates every validation failure
    path: str | Path, *, verify_files: bool = True
) -> Snapshot:
    """Load and validate a snapshot directory or its manifest path."""

    supplied = Path(path)
    manifest_path = supplied if supplied.is_file() else supplied / "snapshot.toml"
    root = manifest_path.parent
    errors: list[str] = []
    if not manifest_path.is_file():
        raise SnapshotValidationError([f"snapshot manifest not found: {manifest_path}"])
    try:
        with manifest_path.open("rb") as stream:
            manifest = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise SnapshotValidationError([f"cannot read snapshot manifest: {exc}"]) from exc

    if manifest.get("schema") != SNAPSHOT_SCHEMA:
        errors.append(f"schema must be {SNAPSHOT_SCHEMA!r}")
    raw_id = str(manifest.get("id", ""))
    try:
        snapshot_id = validate_snapshot_id(raw_id)
    except SnapshotError as exc:
        errors.append(str(exc))
        snapshot_id = raw_id
    if root.name != snapshot_id:
        errors.append(
            f"snapshot id {snapshot_id!r} must match its folder name {root.name!r}"
        )
    detector = str(manifest.get("detector", "")).strip().lower()
    if not detector:
        errors.append("detector is required")
    elif snapshot_id and "_" in snapshot_id:
        id_detector = snapshot_id.split("_", 1)[1].split("-", 1)[0]
        if detector != id_detector:
            errors.append(
                f"detector {detector!r} does not match snapshot id detector {id_detector!r}"
            )
    lamps = manifest.get("lamps")
    if (
        not isinstance(lamps, list)
        or not lamps
        or not all(str(item).strip() for item in lamps)
    ):
        errors.append("lamps must be a non-empty TOML array")

    base_snapshot = manifest.get("base_snapshot")
    if base_snapshot:
        try:
            normalized_base = validate_snapshot_id(str(base_snapshot))
            if normalized_base == snapshot_id:
                errors.append("base_snapshot cannot refer to the snapshot itself")
        except SnapshotError as exc:
            errors.append(f"invalid base_snapshot: {exc}")

    raw_artifacts = manifest.get("artifacts", [])
    if not isinstance(raw_artifacts, list):
        errors.append("artifacts must be an array of TOML tables")
        raw_artifacts = []
    artifacts = tuple(
        artifact
        for index, item in enumerate(raw_artifacts, start=1)
        if (artifact := _parse_artifact(item, index, errors)) is not None
    )
    roles = [artifact.role for artifact in artifacts]
    lamp_names = (
        {str(item).strip().casefold() for item in lamps}
        if isinstance(lamps, list)
        else set()
    )
    for artifact in artifacts:
        if artifact.role != "lamp":
            continue
        if not artifact.label:
            errors.append(f"lamp artifact has no species label: {artifact.path}")
        elif artifact.label.casefold() not in lamp_names:
            errors.append(
                f"lamp artifact label {artifact.label!r} is not listed in the top-level lamps array"
            )
    for role in sorted(REQUIRED_ROLES - set(roles)):
        errors.append(f"required artifact role is missing: {role}")
    for role in sorted(REQUIRED_ROLES):
        if roles.count(role) > 1:
            errors.append(f"required artifact role appears more than once: {role}")

    root_resolved = root.resolve()
    seen_paths: set[str] = set()
    for artifact in artifacts:
        if artifact.path in seen_paths:
            errors.append(f"artifact path appears more than once: {artifact.path}")
            continue
        seen_paths.add(artifact.path)
        relative = _safe_relative_path(artifact.path)
        if relative is None:
            continue
        artifact_path = root / relative
        try:
            resolved = artifact_path.resolve()
        except OSError as exc:
            errors.append(f"cannot resolve artifact {artifact.path}: {exc}")
            continue
        if root_resolved != resolved and root_resolved not in resolved.parents:
            errors.append(f"artifact resolves outside the snapshot: {artifact.path}")
            continue
        if not artifact_path.is_file():
            errors.append(f"artifact file not found: {artifact.path}")
            continue
        if verify_files:
            actual_sha256, actual_size = file_digest(artifact_path)
            if actual_sha256 != artifact.sha256:
                errors.append(f"artifact digest mismatch: {artifact.path}")
            if actual_size != artifact.size_bytes:
                errors.append(f"artifact size mismatch: {artifact.path}")

    if errors:
        raise SnapshotValidationError(errors)
    return Snapshot(root=root, manifest=manifest, artifacts=artifacts)


def create_snapshot(  # noqa: C901 - atomic assembly keeps all cleanup in one scope
    destination_root: str | Path,
    *,
    snapshot_id: str,
    detector: str,
    files: Mapping[str, str | Path],
    lamps: Sequence[str],
    lamp_files: Sequence[tuple[str, str | Path]] = (),
    notes: str = "",
    base_snapshot: str | None = None,
    validity: Mapping[str, object] | None = None,
    alignment: Mapping[str, object] | None = None,
    qc: Mapping[str, object] | None = None,
) -> Snapshot:
    """Assemble an immutable snapshot folder and validate it before publish.

    Construction occurs in a temporary sibling directory.  The final directory
    appears only after every copy, digest, manifest write, and validation has
    succeeded.  An existing snapshot is always refused.
    """

    normalized_id = validate_snapshot_id(snapshot_id)
    detector_name = str(detector).strip().lower()
    id_detector = normalized_id.split("_", 1)[1].split("-", 1)[0]
    if detector_name != id_detector:
        raise SnapshotError(
            f"detector {detector_name!r} must match the detector in snapshot id "
            f"{id_detector!r}"
        )
    lamp_names = tuple(str(item).strip() for item in lamps if str(item).strip())
    if not lamp_names:
        raise SnapshotError("at least one lamp name is required")
    if base_snapshot:
        base_snapshot = validate_snapshot_id(base_snapshot)
        if base_snapshot == normalized_id:
            raise SnapshotError("base_snapshot cannot refer to the snapshot itself")

    missing_roles = REQUIRED_ROLES - set(files)
    extra_roles = set(files) - REQUIRED_ROLES
    if missing_roles:
        raise SnapshotError("missing calibration files: " + ", ".join(sorted(missing_roles)))
    if extra_roles:
        raise SnapshotError("unknown calibration file roles: " + ", ".join(sorted(extra_roles)))

    source_files: dict[str, Path] = {}
    for role, raw_path in files.items():
        source = Path(raw_path)
        if not source.is_file():
            raise SnapshotError(f"{role} source is not a file: {source}")
        source_files[role] = source
    lamp_sources = [(str(label).strip(), Path(path)) for label, path in lamp_files]
    for label, source in lamp_sources:
        if not label:
            raise SnapshotError("every lamp source needs a species label")
        if label.casefold() not in {item.casefold() for item in lamp_names}:
            raise SnapshotError(f"lamp source label {label!r} is not listed in lamps")
        if not source.is_file():
            raise SnapshotError(f"lamp source is not a file: {source}")
    names = [source.name.casefold() for _, source in lamp_sources]
    if len(names) != len(set(names)):
        raise SnapshotError("lamp source filenames must be unique")

    destination_parent = Path(destination_root)
    destination_parent.mkdir(parents=True, exist_ok=True)
    destination = destination_parent / normalized_id
    if destination.exists():
        raise SnapshotError(f"snapshot already exists and will not be replaced: {destination}")

    staging_parent = Path(
        tempfile.mkdtemp(prefix=f".{normalized_id}.staging-", dir=str(destination_parent))
    )
    staging = staging_parent / normalized_id
    staging.mkdir()
    try:
        artifacts: list[Artifact] = []
        for role, filename in ROLE_FILENAMES.items():
            source = source_files[role]
            target = staging / filename
            shutil.copy2(source, target)
            sha256, size_bytes = file_digest(target)
            artifacts.append(
                Artifact(
                    role=role,
                    path=target.relative_to(staging).as_posix(),
                    sha256=sha256,
                    size_bytes=size_bytes,
                    source_name=source.name,
                )
            )
        if lamp_sources:
            lamps_dir = staging / "lamps"
            lamps_dir.mkdir()
            for label, source in lamp_sources:
                target = lamps_dir / source.name
                shutil.copy2(source, target)
                sha256, size_bytes = file_digest(target)
                artifacts.append(
                    Artifact(
                        role="lamp",
                        path=target.relative_to(staging).as_posix(),
                        sha256=sha256,
                        size_bytes=size_bytes,
                        source_name=source.name,
                        label=label,
                    )
                )

        manifest_text = render_manifest(
            snapshot_id=normalized_id,
            detector=detector_name,
            lamps=lamp_names,
            artifacts=artifacts,
            notes=notes,
            base_snapshot=base_snapshot,
            validity=validity,
            alignment=alignment,
            qc=qc,
        )
        (staging / "snapshot.toml").write_text(manifest_text, encoding="utf-8")
        load_snapshot(staging)
        os.replace(staging, destination)
        shutil.rmtree(staging_parent, ignore_errors=True)
    except Exception:
        shutil.rmtree(staging_parent, ignore_errors=True)
        raise
    return load_snapshot(destination)
