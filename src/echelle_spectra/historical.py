"""Thin immutable manifests over historical calibration filenames.

A binder states what a past calibration actually used: the artifact filenames as
they were named at the time, and — where the file travels with this package —
its recorded size and SHA-256.  Artifacts too large to package, such as a
campaign folder's sphere pair, are named without an identity and are captured
when the operator supplies the folder holding them.  The binder is never a copy
of the calibration; :func:`import_historical_snapshot` is what turns one into a
registrable ``echelle-snapshot/v1`` folder.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from datetime import date, datetime, timezone
from importlib.resources import files
from pathlib import Path, PurePosixPath
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

from .snapshot import Snapshot, SnapshotError, create_snapshot, file_digest

HISTORICAL_SCHEMA = "echelle-historical-calibration/v1"

#: Roles a reconstructable snapshot needs before it can be imported at all.
RECONSTRUCTABLE_ROLES = frozenset(
    {"pattern", "wavelength", "sphere", "sphere_background", "integral"}
)

#: The bench writes its solved alignment beside the calibration artifacts under
#: this name; an imported binder that carries one does the same.
ALIGNMENT_SETTINGS_NAME = "alignment.toml"


class HistoricalError(ValueError):
    """A binder is malformed, or its artifacts are not reachable."""


@dataclass(frozen=True)
class HistoricalArtifact:
    role: str
    path: str
    sha256: str
    size_bytes: int
    source_name: str
    label: str = ""
    resolved_path: Path | None = None
    identity_recorded: bool = True

    @property
    def available(self) -> bool:
        """Whether this artifact was found in one of the searched roots."""

        return self.resolved_path is not None


@dataclass(frozen=True)
class HistoricalManifest:
    manifest_path: Path
    snapshot_id: str
    detector: str
    artifacts: tuple[HistoricalArtifact, ...]
    metadata: dict[str, Any]

    @property
    def lamps(self) -> tuple[str, ...]:
        declared = [str(item).strip() for item in self.metadata.get("lamps", [])]
        return tuple(item for item in declared if item)

    @property
    def acquired_date(self) -> str:
        return str(self.metadata.get("acquired_date", ""))

    @property
    def missing(self) -> tuple[HistoricalArtifact, ...]:
        return tuple(item for item in self.artifacts if not item.available)

    def artifact_for_role(self, role: str) -> HistoricalArtifact | None:
        return next((item for item in self.artifacts if item.role == role), None)


def packaged_artifact_root() -> Path:
    """Return the calibration-file root that ships inside this package."""

    return Path(files("echelle_spectra.resources").joinpath("calibration_files"))


def _artifact_roots(extra: tuple[str | Path, ...]) -> tuple[Path, ...]:
    return (packaged_artifact_root(), *(Path(item) for item in extra))


def _relative(item: dict[str, Any]) -> PurePosixPath:
    relative = PurePosixPath(str(item.get("path", "")))
    if not str(relative) or relative.is_absolute() or ".." in relative.parts:
        raise HistoricalError(f"historical artifact path is not portable: {relative}")
    return relative


def _resolve_artifact(
    relative: PurePosixPath, roots: tuple[Path, ...]
) -> Path | None:
    for root in roots:
        candidate = root.joinpath(*relative.parts)
        if candidate.is_file():
            return candidate
    return None


def _load_artifact(
    item: dict[str, Any], roots: tuple[Path, ...]
) -> HistoricalArtifact:
    relative = _relative(item)
    source_name = str(item.get("source_name", ""))
    if source_name != relative.name:
        raise HistoricalError(f"historical source name was not preserved: {relative}")
    declared_sha256 = str(item.get("sha256", "")).strip().lower()
    resolved = _resolve_artifact(relative, roots)
    if resolved is None:
        # The binder still states what was used; only the identity is deferred.
        return HistoricalArtifact(
            role=str(item["role"]),
            path=relative.as_posix(),
            sha256=declared_sha256,
            size_bytes=int(item.get("size_bytes", 0)),
            source_name=source_name,
            label=str(item.get("label", "")),
            resolved_path=None,
            identity_recorded=bool(declared_sha256),
        )
    digest, size = file_digest(resolved)
    if declared_sha256:
        if digest != declared_sha256 or size != int(item.get("size_bytes", -1)):
            raise HistoricalError(f"historical artifact identity mismatch: {relative}")
    return HistoricalArtifact(
        role=str(item["role"]),
        path=relative.as_posix(),
        sha256=digest,
        size_bytes=size,
        source_name=source_name,
        label=str(item.get("label", "")),
        resolved_path=resolved,
        identity_recorded=bool(declared_sha256),
    )


def load_historical_manifest(
    path: str | Path, *, artifact_roots: tuple[str | Path, ...] = ()
) -> HistoricalManifest:
    """Validate a thin binder without copying or renaming its source artifacts.

    Artifacts are looked for in the packaged calibration files first and then in
    each supplied root, in order.  An artifact that no root holds is reported as
    unavailable rather than refused here: listing the historical epochs must not
    require every campaign folder to be mounted.  The import refuses instead.
    """

    manifest_path = Path(path)
    try:
        with manifest_path.open("rb") as stream:
            payload = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise HistoricalError(f"cannot read historical manifest {manifest_path}: {exc}") from exc
    if payload.get("schema") != HISTORICAL_SCHEMA:
        raise HistoricalError(f"historical manifest schema must be {HISTORICAL_SCHEMA!r}")
    roots = _artifact_roots(artifact_roots)
    artifacts = tuple(_load_artifact(item, roots) for item in payload.get("artifacts", []))
    roles = {artifact.role for artifact in artifacts}
    if not RECONSTRUCTABLE_ROLES <= roles:
        raise HistoricalError("historical manifest lacks reconstructable snapshot roles")
    snapshot_id = str(payload.get("snapshot_id", ""))
    if not re.fullmatch(r"\d{8}_[a-z0-9]+", snapshot_id):
        raise HistoricalError(f"invalid historical snapshot identity: {snapshot_id!r}")
    return HistoricalManifest(
        manifest_path=manifest_path,
        snapshot_id=snapshot_id,
        detector=str(payload.get("detector", "")),
        artifacts=artifacts,
        metadata=payload,
    )


def bundled_historical_manifests(
    *, artifact_roots: tuple[str | Path, ...] = ()
) -> tuple[HistoricalManifest, ...]:
    root = Path(files("echelle_spectra.resources").joinpath("historical_manifests"))
    return tuple(
        load_historical_manifest(path, artifact_roots=artifact_roots)
        for path in sorted(root.glob("*.toml"))
    )


def find_bundled_manifest(
    identity: str, *, artifact_roots: tuple[str | Path, ...] = ()
) -> HistoricalManifest:
    """Return one bundled binder by snapshot id, manifest name, or file path."""

    candidate = Path(identity)
    if candidate.is_file():
        return load_historical_manifest(candidate, artifact_roots=artifact_roots)
    bundled = bundled_historical_manifests(artifact_roots=artifact_roots)
    for manifest in bundled:
        if identity in {manifest.snapshot_id, manifest.manifest_path.stem}:
            return manifest
    known = ", ".join(
        f"{manifest.snapshot_id} ({manifest.manifest_path.stem})" for manifest in bundled
    )
    raise HistoricalError(
        f"no bundled historical binder is named {identity!r}; known binders: {known}"
    )


def _alignment_table(settings_path: Path) -> dict[str, object]:
    """Read a bench alignment settings file into a snapshot ``[alignment]`` table."""

    from .tools.calibration_alignment import load_alignment_settings

    try:
        settings = load_alignment_settings(settings_path)
    except (KeyError, OSError, TypeError, ValueError):
        return {}
    return {
        "dx_px": float(settings.transform.dx_px),
        "dy_px": float(settings.transform.dy_px),
        "rotation_deg": float(settings.transform.theta_deg),
        "rms_px": float(settings.rms_px),
    }


def _validity(
    manifest: HistoricalManifest,
    *,
    valid_from: str | None,
    valid_to: str | None,
    shot_from: int | None,
    shot_to: int | None,
) -> dict[str, object]:
    start = valid_from or manifest.acquired_date
    if not start and shot_from is None and shot_to is None:
        raise HistoricalError(
            f"binder {manifest.snapshot_id} records no acquired_date; state the epoch "
            "start with --valid-from YYYY-MM-DD"
        )
    epoch: dict[str, object] = {}
    for name, value in (("date_from", start), ("date_to", valid_to)):
        if not value:
            continue
        try:
            epoch[name] = date.fromisoformat(str(value)).isoformat()
        except ValueError as exc:
            raise HistoricalError(f"--{name.replace('_', '-')} must be ISO YYYY-MM-DD") from exc
    if shot_from is not None:
        epoch["shot_from"] = int(shot_from)
    if shot_to is not None:
        epoch["shot_to"] = int(shot_to)
    return epoch


def _require_reachable(
    manifest: HistoricalManifest, artifact_roots: tuple[str | Path, ...]
) -> None:
    missing = manifest.missing
    if not missing:
        return
    searched = ", ".join(str(root) for root in _artifact_roots(artifact_roots))
    names = ", ".join(f"{item.role}={item.path}" for item in missing)
    raise HistoricalError(
        f"binder {manifest.snapshot_id} names artifact(s) no searched root holds: {names}. "
        f"Searched: {searched}. Supply the campaign folder that holds them with "
        "--artifact-root DIR (repeat the flag for several folders)."
    )


def import_historical_snapshot(
    identity: str,
    destination_root: str | Path,
    *,
    valid_from: str | None = None,
    valid_to: str | None = None,
    shot_from: int | None = None,
    shot_to: int | None = None,
    artifact_roots: tuple[str | Path, ...] = (),
) -> Snapshot:
    """Convert one historical binder into a registrable snapshot folder.

    Every artifact the binder names must be reachable in the packaged
    calibration files or in a supplied ``--artifact-root``; an unreachable one
    is refused by name, because a snapshot that quietly substituted a different
    epoch's file would be worse than no snapshot at all.  The created snapshot
    records where it came from in an ``[imported]`` table, so a reader can tell
    a retroactively imported epoch from one the bench measured live.
    """

    manifest = find_bundled_manifest(identity, artifact_roots=artifact_roots)
    _require_reachable(manifest, artifact_roots)
    files_by_role = {
        role: manifest.artifact_for_role(role).resolved_path
        for role in sorted(RECONSTRUCTABLE_ROLES)
    }
    lamp_files = [
        (artifact.label, artifact.resolved_path)
        for artifact in manifest.artifacts
        if artifact.role == "lamp"
    ]
    lamps = manifest.lamps or tuple(sorted({label for label, _ in lamp_files if label}))
    if not lamps:
        raise HistoricalError(
            f"binder {manifest.snapshot_id} declares no lamps; add a lamps array to the binder"
        )
    unlabelled = [
        artifact.path
        for artifact in manifest.artifacts
        if artifact.role == "lamp" and not artifact.label
    ]
    if unlabelled:
        raise HistoricalError(
            "historical lamp artifact(s) carry no species label: " + ", ".join(unlabelled)
        )
    settings = manifest.artifact_for_role("alignment_settings")
    binder_sha256, _ = file_digest(manifest.manifest_path)
    try:
        snapshot = create_snapshot(
            destination_root,
            snapshot_id=manifest.snapshot_id,
            detector=manifest.detector,
            files=files_by_role,
            lamps=lamps,
            lamp_files=lamp_files,
            notes=str(manifest.metadata.get("notes", "")),
            validity=_validity(
                manifest,
                valid_from=valid_from,
                valid_to=valid_to,
                shot_from=shot_from,
                shot_to=shot_to,
            ),
            alignment=(
                _alignment_table(settings.resolved_path)
                if settings is not None and settings.resolved_path is not None
                else None
            ),
            imported={
                "from_schema": HISTORICAL_SCHEMA,
                "binder": manifest.manifest_path.name,
                "binder_sha256": binder_sha256,
                "imported_utc": datetime.now(timezone.utc)
                .replace(microsecond=0)
                .isoformat(),
            },
        )
    except SnapshotError as exc:
        raise HistoricalError(f"historical import failed: {exc}") from exc
    if settings is not None and settings.resolved_path is not None:
        shutil.copy2(settings.resolved_path, snapshot.root / ALIGNMENT_SETTINGS_NAME)
    return snapshot
