"""Thin immutable manifests over historical calibration filenames."""

from __future__ import annotations

import re
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path, PurePosixPath
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib

from .snapshot import file_digest

HISTORICAL_SCHEMA = "echelle-historical-calibration/v1"


@dataclass(frozen=True)
class HistoricalArtifact:
    role: str
    path: str
    sha256: str
    size_bytes: int
    source_name: str


@dataclass(frozen=True)
class HistoricalManifest:
    manifest_path: Path
    snapshot_id: str
    detector: str
    artifacts: tuple[HistoricalArtifact, ...]
    metadata: dict[str, Any]


def load_historical_manifest(path: str | Path) -> HistoricalManifest:
    """Validate a thin binder without copying or renaming its source artifacts."""

    manifest_path = Path(path)
    with manifest_path.open("rb") as stream:
        payload = tomllib.load(stream)
    if payload.get("schema") != HISTORICAL_SCHEMA:
        raise ValueError(f"historical manifest schema must be {HISTORICAL_SCHEMA!r}")
    resource_root = Path(files("echelle_spectra.resources").joinpath("calibration_files"))
    artifacts = []
    for item in payload.get("artifacts", []):
        relative = PurePosixPath(str(item.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"historical artifact path is not portable: {relative}")
        artifact_path = resource_root.joinpath(*relative.parts)
        digest, size = file_digest(artifact_path)
        if digest != item.get("sha256") or size != int(item.get("size_bytes", -1)):
            raise ValueError(f"historical artifact identity mismatch: {relative}")
        if artifact_path.name != item.get("source_name"):
            raise ValueError(f"historical source name was not preserved: {relative}")
        artifacts.append(
            HistoricalArtifact(
                role=str(item["role"]),
                path=relative.as_posix(),
                sha256=digest,
                size_bytes=size,
                source_name=artifact_path.name,
            )
        )
    required = {"pattern", "wavelength", "sphere", "sphere_background", "integral"}
    roles = {artifact.role for artifact in artifacts}
    if not required <= roles:
        raise ValueError("historical manifest lacks reconstructable snapshot roles")
    snapshot_id = str(payload.get("snapshot_id", ""))
    if not re.fullmatch(r"\d{8}_[a-z0-9]+", snapshot_id):
        raise ValueError(f"invalid historical snapshot identity: {snapshot_id!r}")
    return HistoricalManifest(
        manifest_path=manifest_path,
        snapshot_id=snapshot_id,
        detector=str(payload.get("detector", "")),
        artifacts=tuple(artifacts),
        metadata=payload,
    )


def bundled_historical_manifests() -> tuple[HistoricalManifest, ...]:
    root = Path(files("echelle_spectra.resources").joinpath("historical_manifests"))
    return tuple(load_historical_manifest(path) for path in sorted(root.glob("*.toml")))
