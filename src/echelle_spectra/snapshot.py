"""Calibration snapshots: immutable, digested campaign artifacts.

A snapshot directory is identified by its folder name and contains one
``snapshot.toml`` binder.  The binder names every calibration input by role and
records a SHA-256 digest, so downstream cubes and run receipts can refer to one
stable snapshot id instead of a loose collection of filenames.

An artifact reaches the binder in one of two ways, and the binder says which.
A *copied* artifact lives inside the snapshot folder at a relative path that
cannot leave it — that is what the computed files (the pattern, the corrected
wavelength table, the sphere's spectral reference) always are.  A *referenced*
artifact keeps living where it was measured, and the binder records the path
back to it beside the same digest and size.  Raw detector frames are referenced:
the calibration folder already holds the lamp and sphere SIFs, and a snapshot
sitting two levels below them has nothing to gain from a second copy of 380 MB.

The digest stays the identity authority either way.  The path only says where
the bytes live, and validation re-reads them there.
"""

from __future__ import annotations

import dataclasses
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

#: An artifact whose bytes were copied into the snapshot folder.  This is what a
#: binder means when it states no ``kind`` at all, which is what every snapshot
#: written before references existed does.
ARTIFACT_COPIED = "copied"
#: An artifact the snapshot points back at instead of copying.  ``path`` is then
#: POSIX-relative to the snapshot root (``../../sphere-0.1s-x3.sif``) when the
#: source shares the calibration folder, and absolute when it does not.
ARTIFACT_REFERENCED = "referenced"
ARTIFACT_KINDS = (ARTIFACT_COPIED, ARTIFACT_REFERENCED)

#: The raw detector frames a thin snapshot references rather than copies.  The
#: computed files are not here on purpose: the snapshot owns those outright.
RAW_SOURCE_ROLES = frozenset({"sphere", "sphere_background", "lamp"})


class SnapshotError(ValueError):
    """Base error for an invalid snapshot request or artifact."""


class SnapshotValidationError(SnapshotError):
    """A snapshot failed one or more integrity checks."""

    def __init__(self, errors: Sequence[str]):
        self.errors = tuple(errors)
        super().__init__("; ".join(self.errors))


@dataclass(frozen=True)
class Artifact:
    """One file recorded in a snapshot manifest.

    ``path`` is read according to ``kind``: inside the snapshot for a copied
    artifact, and back out to the measured file for a referenced one.
    ``resolved_path`` is filled in by :func:`load_snapshot`, which is the moment
    the manifest text becomes one absolute location on this machine.
    """

    role: str
    path: str
    sha256: str
    size_bytes: int
    source_name: str = ""
    label: str = ""
    kind: str = ARTIFACT_COPIED
    resolved_path: Path | None = None

    @property
    def is_reference(self) -> bool:
        """Whether the bytes live outside the snapshot folder."""

        return self.kind == ARTIFACT_REFERENCED


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

    def path_for(self, artifact: Artifact) -> Path:
        """Return where *artifact*'s bytes actually are on this machine."""

        if artifact.resolved_path is not None:
            return artifact.resolved_path
        if artifact.is_reference:
            return reference_target(self.root, artifact.path)
        return self.root / artifact.path

    def source_path(self, role: str) -> Path:
        """Return the readable location of the unique artifact for *role*."""

        return self.path_for(self.artifact_for_role(role))

    def _calibration_name(self, role: str) -> str:
        """Name one role the way ``build_calibration`` wants to receive it.

        A copied artifact keeps its bare role filename, which the loader
        resolves inside the snapshot folder and therefore cannot be shadowed by
        the working directory.  A referenced one is handed over already resolved,
        because no folder-relative name could reach it.
        """

        artifact = self.artifact_for_role(role)
        if artifact.is_reference:
            return str(self.path_for(artifact))
        return artifact.path

    def calibration_files(self) -> dict[str, str]:
        """Return the established ``Calibrations`` filename vocabulary."""

        return {
            "orders": self._calibration_name("pattern"),
            "wavelength": self._calibration_name("wavelength"),
            "sphr": self._calibration_name("sphere"),
            "bkgr": self._calibration_name("sphere_background"),
            "integral": self._calibration_name("integral"),
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
                # Stated only when it is not the historical default, so a cube
                # written from an all-copied snapshot carries the same
                # provenance string it always did.
                **({"kind": artifact.kind} if artifact.is_reference else {}),
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


def _is_absolute_text(text: str) -> bool:
    """Whether *text* names an absolute path in either path flavour."""

    windows = PureWindowsPath(text)
    return PurePosixPath(text).is_absolute() or windows.is_absolute() or bool(windows.drive)


def reference_target(snapshot_root: str | Path, recorded: str) -> Path:
    """Resolve a referenced artifact's recorded path to one absolute location.

    Relative references are read against the snapshot folder rather than the
    working directory, so a snapshot resolves to the same files whatever
    directory the process was started in.
    """

    text = str(recorded)
    if _is_absolute_text(text):
        return Path(text).resolve()
    parts = PurePosixPath(text).parts
    return Path(snapshot_root).joinpath(*parts).resolve()


def reference_path(snapshot_root: str | Path, source: str | Path) -> str:
    """Return the manifest text pointing from a snapshot back at *source*.

    A source inside the calibration folder that holds the snapshot is recorded
    relative to the snapshot root — ``../../sphere-0.1s-x3.sif`` for the ordinary
    bench layout — so the whole calibration folder can be moved or copied to
    another machine without a single edit.  A source from anywhere else is
    recorded absolutely, because a relative path across unrelated trees would
    only pretend to be portable.
    """

    root = Path(snapshot_root).resolve()
    target = Path(source).resolve()
    parents = list(root.parents)
    # The calibration folder is the snapshot root's grandparent in the bench
    # layout, <calibration folder>/calibrations/<id>; a shallower root simply
    # offers less tree to stay inside.
    tree = parents[1] if len(parents) >= 2 else root.parent
    try:
        target.relative_to(tree)
        relative = os.path.relpath(target, root)
    except ValueError:
        return target.as_posix()
    return PurePosixPath(*Path(relative).parts).as_posix()


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
    imported: Mapping[str, object] | None = None,
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
        # Where a retroactively converted epoch came from. Live bench snapshots
        # simply have no [imported] table.
        ("imported", imported),
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
        lines.extend(["", "[[artifacts]]", f"role = {_toml_string(artifact.role)}"])
        if artifact.is_reference:
            # Said out loud, because a reader must never have to guess from the
            # shape of a path whether the bytes are in here or out there.
            lines.append(f"kind = {_toml_string(artifact.kind)}")
        lines.extend(
            [
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
    # No key means the artifact was copied in, which is what every binder
    # written before references existed says by saying nothing.
    kind = str(item.get("kind", ARTIFACT_COPIED)).strip().lower()
    try:
        size_bytes = int(item.get("size_bytes", -1))
    except (TypeError, ValueError):
        size_bytes = -1
    if not role:
        errors.append(f"artifact {index} has no role")
    if kind not in ARTIFACT_KINDS:
        errors.append(
            f"artifact {index} has an unknown kind {kind!r}; expected one of "
            + ", ".join(ARTIFACT_KINDS)
        )
        return None
    path_ok = bool(path) if kind == ARTIFACT_REFERENCED else _safe_relative_path(path) is not None
    if not path_ok:
        errors.append(
            f"artifact {index} has no path"
            if kind == ARTIFACT_REFERENCED
            else f"artifact {index} path must stay inside the snapshot: {path!r}"
        )
    if not re.fullmatch(r"[0-9a-f]{64}", sha256):
        errors.append(f"artifact {index} has an invalid SHA-256 digest")
    if size_bytes < 0:
        errors.append(f"artifact {index} has an invalid size_bytes")
    if not role or not path_ok or size_bytes < 0:
        return None
    return Artifact(
        role=role,
        path=path,
        sha256=sha256,
        size_bytes=size_bytes,
        source_name=str(item.get("source_name", "")),
        label=str(item.get("label", "")),
        kind=kind,
    )


def _verify_artifact(  # noqa: C901 - one artifact, every way it can be wrong
    artifact: Artifact,
    *,
    root: Path,
    reference_root: Path,
    verify_files: bool,
    verify_references: bool = True,
    errors: list[str],
) -> Artifact:
    """Locate one artifact, check it, and return it carrying that location."""

    if artifact.is_reference:
        if not verify_references:
            # Not even a stat.  A referenced frame is 380 MB of raw detector
            # sitting wherever it was measured, which on this instrument is
            # routinely a share that may not be mounted; a light reading of the
            # binder must not go looking for it, let alone hash it.  The path
            # stays unresolved and :meth:`Snapshot.path_for` works it out on
            # demand, as it always has.
            return artifact
        target = reference_target(reference_root, artifact.path)
        if not target.is_file():
            # The absolute path, because "not found: ../../sphere.sif" tells an
            # operator nothing about which folder was actually looked in.
            errors.append(f"referenced {artifact.role} source not found: {target}")
            return artifact
        if verify_files:
            actual_sha256, actual_size = file_digest(target)
            if actual_sha256 != artifact.sha256:
                errors.append(f"referenced {artifact.role} source digest mismatch: {target}")
            if actual_size != artifact.size_bytes:
                errors.append(f"referenced {artifact.role} source size mismatch: {target}")
        return dataclasses.replace(artifact, resolved_path=target)

    relative = _safe_relative_path(artifact.path)
    if relative is None:
        return artifact
    artifact_path = root / relative
    try:
        resolved = artifact_path.resolve()
    except OSError as exc:
        errors.append(f"cannot resolve artifact {artifact.path}: {exc}")
        return artifact
    root_resolved = root.resolve()
    if root_resolved != resolved and root_resolved not in resolved.parents:
        errors.append(f"artifact resolves outside the snapshot: {artifact.path}")
        return artifact
    if not artifact_path.is_file():
        errors.append(f"artifact file not found: {artifact.path}")
        return artifact
    if verify_files:
        actual_sha256, actual_size = file_digest(artifact_path)
        if actual_sha256 != artifact.sha256:
            errors.append(f"artifact digest mismatch: {artifact.path}")
        if actual_size != artifact.size_bytes:
            errors.append(f"artifact size mismatch: {artifact.path}")
    return dataclasses.replace(artifact, resolved_path=resolved)


def load_snapshot(  # noqa: C901 - one pass accumulates every validation failure
    path: str | Path,
    *,
    verify_files: bool = True,
    verify_references: bool = True,
    reference_root: str | Path | None = None,
) -> Snapshot:
    """Load and validate a snapshot directory or its manifest path.

    ``reference_root`` states where referenced sources should be resolved from
    when the folder is not yet standing in its final place — which is exactly
    the situation :func:`create_snapshot` validates from, one staging directory
    away from where the snapshot is about to live.

    ``verify_references=False`` is the *light* check: schema, roles, and the
    digests of the artifacts the snapshot owns — all of them kilobyte-sized
    files inside the folder — with the referenced raw frames left entirely
    alone.  It is what a reader asking "is this calibration done?" needs, and
    it answers in milliseconds whether or not the frames are reachable.
    """

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

    resolution_root = Path(reference_root) if reference_root is not None else root
    seen_paths: set[str] = set()
    located: list[Artifact] = []
    for artifact in artifacts:
        if artifact.path in seen_paths:
            errors.append(f"artifact path appears more than once: {artifact.path}")
            located.append(artifact)
            continue
        seen_paths.add(artifact.path)
        located.append(
            _verify_artifact(
                artifact,
                root=root,
                reference_root=resolution_root,
                verify_files=verify_files,
                verify_references=verify_references,
                errors=errors,
            )
        )
    artifacts = tuple(located)

    if errors:
        raise SnapshotValidationError(errors)
    return Snapshot(root=root, manifest=manifest, artifacts=artifacts)


@dataclass(frozen=True)
class SnapshotReading:
    """One saved snapshot folder, read the way somebody glancing at it would.

    A snapshot that has been saved says so only in its own files, and an
    operator who has to open ``snapshot.toml`` and compare dates before he can
    tell whether this morning's work landed has not been told anything (owner,
    2026-08-18: "I'm not sure if that's done or not. No clear indications for
    me").  This is the one object that answers it: what the calibration is
    called, when it was saved, which lamps it used, how tight the alignment
    came out, what the sphere comparison said, and whether the folder still
    verifies against its own digests.

    It never raises.  A binder that cannot be parsed, or that fails validation,
    is a reading with ``errors`` — because "this one is broken, here is the
    file" is an answer and an exception is not.
    """

    root: Path
    snapshot_id: str
    created_utc: str = ""
    lamps: tuple[str, ...] = ()
    alignment_rms_px: float | None = None
    sphere_comparison: str = ""
    errors: tuple[str, ...] = ()

    @property
    def valid(self) -> bool:
        """Whether the light check found nothing wrong."""

        return not self.errors

    @property
    def saved_local(self) -> str:
        """``created_utc`` on the clock the operator is reading it by.

        The binder records UTC, which is right for a file and wrong for the
        question "did I save this this morning?".  A stamp that cannot be
        parsed is handed back as it was written rather than dropped.
        """

        text = self.created_utc.strip()
        if not text:
            return ""
        try:
            moment = datetime.fromisoformat(text.replace("Z", "+00:00"))
        except ValueError:
            return text
        if moment.tzinfo is None:
            moment = moment.replace(tzinfo=timezone.utc)
        return moment.astimezone().strftime("%Y-%m-%d %H:%M")

    def summary(self) -> str:
        """One plain sentence-fragment naming this calibration and its verdict.

        Every part is stated only when the binder has it, so an older or
        hand-written snapshot reads short rather than reading a row of dashes.
        The verdict is always last, and always there.
        """

        parts: list[str] = []
        if self.saved_local:
            parts.append(f"saved {self.saved_local}")
        if self.lamps:
            parts.append("/".join(self.lamps))
        if self.alignment_rms_px is not None:
            parts.append(f"RMS {self.alignment_rms_px:.2f} px")
        if self.sphere_comparison:
            parts.append(f"sphere comparison {self.sphere_comparison.replace('-', ' ')}")
        parts.append(
            "validated"
            if self.valid
            else f"DOES NOT VALIDATE — {self.errors[0]}"
        )
        return f"{self.snapshot_id} — " + ", ".join(parts)


def _manifest_table(manifest: Mapping[str, Any], heading: str) -> Mapping[str, Any]:
    """One optional TOML table, or an empty one — never a TypeError."""

    table = manifest.get(heading)
    return table if isinstance(table, dict) else {}


def read_snapshot_folder(path: str | Path) -> SnapshotReading:
    """Read one saved snapshot folder cheaply enough to say so at first paint.

    The verdict comes from :func:`load_snapshot`'s light check — the same
    schema, role and digest checks ``echelle snapshot validate`` runs, minus
    the referenced raw frames — so nothing here reads more than a few kilobytes
    inside the folder itself.
    """

    root = Path(path)
    try:
        with (root / "snapshot.toml").open("rb") as stream:
            manifest: Mapping[str, Any] = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError):
        # Unparseable: the validator below says so in full, and the displayable
        # facts simply are not there to display.
        manifest = {}
    try:
        load_snapshot(root, verify_references=False)
    except SnapshotValidationError as exc:
        errors = tuple(exc.errors)
    else:
        errors = ()

    raw_lamps = manifest.get("lamps")
    rms = _manifest_table(manifest, "alignment").get("rms_px")
    comparison = _manifest_table(manifest, "qc").get("sphere_comparison")
    return SnapshotReading(
        root=root,
        # The folder name is the fallback because a binder whose id is missing
        # is exactly the binder whose folder name has to carry the identity.
        snapshot_id=str(manifest.get("id", "") or root.name),
        created_utc=str(manifest.get("created_utc", "")),
        lamps=tuple(
            str(item).strip()
            for item in (raw_lamps if isinstance(raw_lamps, list) else ())
            if str(item).strip()
        ),
        alignment_rms_px=(
            float(rms) if isinstance(rms, (int, float)) and not isinstance(rms, bool) else None
        ),
        sphere_comparison=str(comparison or "").strip(),
        errors=errors,
    )


def saved_snapshots_in(output_root: str | Path) -> tuple[SnapshotReading, ...]:
    """Every snapshot saved directly inside *output_root*, newest first.

    One level deep and no deeper: a snapshot lives at ``<output root>/<id>``,
    the settings bundles live beside them in a folder that holds no binder, and
    the raw frames the binders point at can be anywhere at all — including a
    share whose recursive walk would take the bench off the air.  Anything that
    is not a directory holding a ``snapshot.toml`` is simply not a snapshot.
    """

    root = Path(output_root)
    try:
        entries = sorted(root.iterdir())
    except OSError:
        return ()
    readings: list[SnapshotReading] = []
    for entry in entries:
        try:
            holds_binder = entry.is_dir() and (entry / "snapshot.toml").is_file()
        except OSError:  # pragma: no cover - an entry that vanished mid-scan
            continue
        if holds_binder:
            readings.append(read_snapshot_folder(entry))
    # ``created_utc`` is ISO-8601, so it sorts as text; the id breaks ties
    # between two snapshots written inside the same second.
    readings.sort(key=lambda reading: (reading.created_utc, reading.snapshot_id), reverse=True)
    return tuple(readings)


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
    imported: Mapping[str, object] | None = None,
    reference_raw: bool = False,
) -> Snapshot:
    """Assemble an immutable snapshot folder and validate it before publish.

    Construction occurs in a temporary sibling directory.  The final directory
    appears only after every copy, digest, manifest write, and validation has
    succeeded.  An existing snapshot is always refused.

    With ``reference_raw`` the raw detector frames — the sphere pair and every
    lamp — are digested where they are and recorded as references instead of
    being copied in.  The computed files are still copied, so the snapshot owns
    everything it computed and points at everything it merely read.
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
    if not reference_raw and len(names) != len(set(names)):
        # Only copies collide: two referenced lamps keep their own folders.
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

        def referenced(role: str, source: Path, label: str = "") -> Artifact:
            """Digest *source* where it lies and point the manifest at it."""

            sha256, size_bytes = file_digest(source)
            return Artifact(
                role=role,
                # Relative to the published folder, never to the staging one it
                # is assembled in: the recorded path must survive the rename.
                path=reference_path(destination, source),
                sha256=sha256,
                size_bytes=size_bytes,
                source_name=source.name,
                label=label,
                kind=ARTIFACT_REFERENCED,
            )

        def copied(role: str, source: Path, target: Path, label: str = "") -> Artifact:
            # Not copy2: copy2's copystat forwards BSD file flags, and a source
            # inside Dropbox carries UF_TRACKED, which an exFAT destination (a
            # field drive) answers with EINVAL -- an errno copystat does not
            # forgive, crashing the save after every number was computed.  The
            # bytes and the mtime are what the manifest cares about; flags a
            # drive cannot hold are dropped, not fatal.
            shutil.copyfile(source, target)
            try:
                shutil.copystat(source, target)
            except OSError:
                stat = source.stat()
                os.utime(target, (stat.st_atime, stat.st_mtime))
            sha256, size_bytes = file_digest(target)
            return Artifact(
                role=role,
                path=target.relative_to(staging).as_posix(),
                sha256=sha256,
                size_bytes=size_bytes,
                source_name=source.name,
                label=label,
            )

        for role, filename in ROLE_FILENAMES.items():
            source = source_files[role]
            if reference_raw and role in RAW_SOURCE_ROLES:
                artifacts.append(referenced(role, source))
            else:
                artifacts.append(copied(role, source, staging / filename))
        if lamp_sources and reference_raw:
            artifacts.extend(
                referenced("lamp", source, label=label) for label, source in lamp_sources
            )
        elif lamp_sources:
            lamps_dir = staging / "lamps"
            lamps_dir.mkdir()
            for label, source in lamp_sources:
                artifacts.append(
                    copied("lamp", source, lamps_dir / source.name, label=label)
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
            imported=imported,
        )
        (staging / "snapshot.toml").write_text(manifest_text, encoding="utf-8")
        load_snapshot(staging, reference_root=destination)
        os.replace(staging, destination)
        shutil.rmtree(staging_parent, ignore_errors=True)
    except Exception:
        shutil.rmtree(staging_parent, ignore_errors=True)
        raise
    return load_snapshot(destination)
