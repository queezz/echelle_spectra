"""Science-line drift sampling, verdicts, and immutable refinements."""

from __future__ import annotations

import json
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .campaign_run import sha256_file
from .snapshot import REQUIRED_ROLES, Snapshot, create_snapshot, load_snapshot
from .tools.line_catalog import load_line_table

DRIFT_SCHEMA = "echelle-drift-evidence/v1"
ALIGNMENT_THRESHOLD_NM = 0.03
ALIGNMENT_MAX_RESIDUAL_NM = 0.06
REPAIR_LIMIT_NM = 0.25
SHIFT_CONSISTENCY_NM = 0.04
MINIMUM_CENTROIDS = 2
MINIMUM_SNR = 4.0


class DriftError(ValueError):
    """Drift evidence or refinement is unsafe or incomplete."""


def select_sample_paths(
    paths: list[str | Path], *, every: int = 1, shots: set[str] | None = None
) -> list[Path]:
    """Select interval samples plus explicitly named shots."""

    if every < 1:
        raise DriftError("sample interval must be at least 1")
    ordered = sorted(Path(path) for path in paths)
    chosen = {path for index, path in enumerate(ordered) if index % every == 0}
    if shots:
        for path in ordered:
            if path.stem in shots or any(shot in path.stem for shot in shots):
                chosen.add(path)
    return sorted(chosen)


def _spectrum(ds) -> tuple[np.ndarray, np.ndarray]:
    wavelength = np.asarray(ds.coords["wavelength"].values, dtype=float)
    intensity = ds["intensity"]
    reduce_dims = [dim for dim in intensity.dims if dim != "wavelength"]
    if reduce_dims:
        intensity = intensity.median(dim=reduce_dims, skipna=True)
    return wavelength, np.asarray(intensity.values, dtype=float)


def centroid_evidence(
    wavelength: np.ndarray,
    intensity: np.ndarray,
    *,
    expected_nm: float,
    window_nm: float = 0.4,
) -> dict[str, Any]:
    """Measure a baseline-subtracted centroid with explicit sufficiency evidence."""

    selected = np.abs(wavelength - expected_nm) <= window_nm
    x = wavelength[selected]
    y = intensity[selected]
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 5:
        return {"status": "insufficient-data", "reason": "fewer than five samples in window"}
    edge_count = max(1, x.size // 5)
    baseline = float(np.median(np.r_[y[:edge_count], y[-edge_count:]]))
    signal = np.clip(y - baseline, 0.0, None)
    noise = float(np.std(np.r_[y[:edge_count], y[-edge_count:]] - baseline))
    peak = float(np.max(signal))
    snr = peak / max(noise, np.finfo(float).eps)
    if not np.isfinite(snr) or snr < MINIMUM_SNR or signal.sum() <= 0:
        return {
            "status": "insufficient-data",
            "reason": f"peak SNR {snr:.3g} is below {MINIMUM_SNR:g}",
            "snr": snr,
        }
    centroid = float(np.sum(x * signal) / np.sum(signal))
    return {
        "status": "measured",
        "centroid_nm": centroid,
        "residual_nm": centroid - expected_nm,
        "snr": snr,
        "window_nm": window_nm,
    }


def verdict_from_evidence(evidence: list[dict[str, Any]]) -> tuple[str, dict[str, float]]:
    residuals = np.asarray(
        [item["residual_nm"] for item in evidence if item.get("status") == "measured"],
        dtype=float,
    )
    if residuals.size < MINIMUM_CENTROIDS:
        return "insufficient-data", {}
    median = float(np.median(residuals))
    maximum = float(np.max(np.abs(residuals)))
    spread = float(np.max(np.abs(residuals - median)))
    summary = {
        "median_shift_nm": median,
        "maximum_absolute_residual_nm": maximum,
        "maximum_shift_deviation_nm": spread,
    }
    if abs(median) <= ALIGNMENT_THRESHOLD_NM and maximum <= ALIGNMENT_MAX_RESIDUAL_NM:
        return "aligned", summary
    if abs(median) <= REPAIR_LIMIT_NM and spread <= SHIFT_CONSISTENCY_NM:
        return "shifted", summary
    return "misaligned-beyond-repair", summary


def audit_cubes(
    cube_paths: list[str | Path],
    *,
    every: int = 1,
    shots: set[str] | None = None,
    families: tuple[str, ...] = ("balmer", "fulcher"),
) -> dict[str, Any]:
    """Audit sampled saved cubes and return an immutable JSON-ready verdict."""

    import xarray as xr

    selected = select_sample_paths(cube_paths, every=every, shots=shots)
    per_line: list[dict[str, Any]] = []
    snapshot_ids: set[str] = set()
    sampled = []
    for path in selected:
        with xr.open_dataset(path) as opened:
            ds = opened.load()
        snapshot_ids.add(str(ds.attrs.get("snapshot_id", "unassigned")))
        wavelength, intensity = _spectrum(ds)
        shot_number = str(ds.attrs.get("shot_number", path.stem))
        sampled.append(
            {
                "cube": path.name,
                "sha256": sha256_file(path),
                "shot_number": shot_number,
            }
        )
        minimum, maximum = float(np.min(wavelength)), float(np.max(wavelength))
        for family in families:
            for line in load_line_table(family):
                if not minimum <= line.wavelength_nm <= maximum:
                    continue
                result = centroid_evidence(
                    wavelength, intensity, expected_nm=line.wavelength_nm
                )
                per_line.append(
                    {
                        "cube": path.name,
                        "shot_number": shot_number,
                        "family": family,
                        "line": line.label,
                        "expected_nm": line.wavelength_nm,
                        "source_reference": line.source_reference,
                        **result,
                    }
                )
    verdict, summary = verdict_from_evidence(per_line)
    repair = ""
    if verdict == "shifted":
        repair = (
            "echelle drift refine DRIFT_EVIDENCE.json --calibrations calibrations "
            f"--accept-shift {summary['median_shift_nm']:.9g}"
        )
    return {
        "schema": DRIFT_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample_rule": {"every": every, "shots": sorted(shots or set())},
        "snapshot_ids": sorted(snapshot_ids),
        "sampled_cubes": sampled,
        "thresholds_nm": {
            "aligned_median": ALIGNMENT_THRESHOLD_NM,
            "aligned_maximum": ALIGNMENT_MAX_RESIDUAL_NM,
            "repair_limit": REPAIR_LIMIT_NM,
            "shift_consistency": SHIFT_CONSISTENCY_NM,
            "minimum_centroids": MINIMUM_CENTROIDS,
            "minimum_snr": MINIMUM_SNR,
        },
        "verdict": verdict,
        "summary": summary,
        "lines": per_line,
        "repair_command": repair,
    }


def write_drift_evidence(path: str | Path, payload: dict[str, Any]) -> Path:
    destination = Path(path)
    if destination.exists():
        raise FileExistsError(f"drift evidence is immutable and already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return destination


def _shift_wavelength_table(source: Path, destination: Path, correction_nm: float) -> None:
    output = []
    for line in source.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            output.append(line)
            continue
        fields = stripped.split()
        if len(fields) < 5:
            output.append(line)
            continue
        try:
            fields[4] = f"{float(fields[4]) + correction_nm:.9f}"
        except ValueError:
            output.append(line)
        else:
            output.append("\t".join(fields))
    destination.write_text("\n".join(output) + "\n", encoding="utf-8", newline="\n")


def _next_refinement_id(base: Snapshot, root: Path) -> str:
    stem = re.sub(r"-r\d+$", "", base.snapshot_id)
    revision = 1
    while (root / f"{stem}-r{revision}").exists():
        revision += 1
    return f"{stem}-r{revision}"


def create_refinement_snapshot(
    evidence_path: str | Path,
    *,
    calibrations_root: str | Path,
    accepted_shift_nm: float,
) -> tuple[Snapshot, Path]:
    """Accept a shifted verdict and emit a new immutable ``-rN`` snapshot."""

    evidence_file = Path(evidence_path)
    accepted_path = evidence_file.with_suffix(".accepted.json")
    if accepted_path.exists():
        raise FileExistsError(f"accepted drift record already exists: {accepted_path}")
    evidence = json.loads(evidence_file.read_text(encoding="utf-8"))
    if evidence.get("schema") != DRIFT_SCHEMA or evidence.get("verdict") != "shifted":
        raise DriftError("only a shifted sampled verdict can create a refinement snapshot")
    measured = float(evidence.get("summary", {}).get("median_shift_nm", np.nan))
    if not np.isfinite(measured) or not np.isclose(
        accepted_shift_nm, measured, rtol=0.0, atol=1e-9
    ):
        raise DriftError("--accept-shift must exactly acknowledge the sampled median shift")
    snapshot_ids = list(evidence.get("snapshot_ids", []))
    if len(snapshot_ids) != 1:
        raise DriftError("refinement evidence must refer to exactly one base snapshot")
    root = Path(calibrations_root)
    base = load_snapshot(root / snapshot_ids[0])
    refined_id = _next_refinement_id(base, root)
    with tempfile.TemporaryDirectory(prefix="echelle-refinement-") as temporary:
        adjusted = Path(temporary) / "wavelength.txt"
        _shift_wavelength_table(
            base.root / base.artifact_for_role("wavelength").path,
            adjusted,
            correction_nm=-measured,
        )
        files = {
            role: (
                adjusted
                if role == "wavelength"
                else base.root / base.artifact_for_role(role).path
            )
            for role in REQUIRED_ROLES
        }
        lamps = [
            (artifact.label, base.root / artifact.path)
            for artifact in base.artifacts
            if artifact.role == "lamp"
        ]
        snapshot = create_snapshot(
            root,
            snapshot_id=refined_id,
            detector=base.detector,
            files=files,
            lamps=base.lamps,
            lamp_files=lamps,
            notes=f"Accepted science-line refinement from {evidence_file.name}.",
            base_snapshot=base.snapshot_id,
            validity=base.manifest.get("validity"),
            alignment={
                "method": "sampled Balmer/Fulcher rigid wavelength shift",
                "observed_shift_nm": measured,
                "applied_correction_nm": -measured,
                "evidence_sha256": sha256_file(evidence_file),
            },
            qc={"sampled_verdict": "shifted", "accepted": True},
        )
    accepted = {
        **evidence,
        "accepted": True,
        "accepted_snapshot_id": snapshot.snapshot_id,
        # Accepting a shift condemns the audited snapshot: only the refinement
        # carries the corrected wavelength table, so only the refinement is
        # authorized for bulk processing.
        "snapshot_ids": [snapshot.snapshot_id],
        "base_snapshot_ids": sorted(snapshot_ids),
        "base_evidence_sha256": sha256_file(evidence_file),
    }
    write_drift_evidence(accepted_path, accepted)
    return snapshot, accepted_path


def _require_authorized_snapshots(selected: set[str], authorized: set[str]) -> None:
    if not selected <= authorized:
        missing = ", ".join(sorted(selected - authorized))
        raise DriftError(f"drift evidence does not cover selected snapshot(s): {missing}")


def require_sampled_verdict(path: str | Path, snapshot_ids: set[str]) -> dict[str, Any]:
    """Validate the bulk-processing prerequisite without treating insufficiency as aligned."""

    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != DRIFT_SCHEMA:
        raise DriftError("unsupported drift evidence schema")
    verdict = payload.get("verdict")
    if verdict == "aligned":
        _require_authorized_snapshots(
            snapshot_ids, set(map(str, payload.get("snapshot_ids", [])))
        )
        return payload
    if verdict == "shifted" and payload.get("accepted"):
        refined = str(payload.get("accepted_snapshot_id") or "")
        recorded = set(map(str, payload.get("snapshot_ids", [])))
        # Records written before the refinement narrowed its own authorization
        # list both ids; everything that is not the refinement is superseded.
        superseded = set(map(str, payload.get("base_snapshot_ids") or [])) or (
            recorded - {refined}
        )
        condemned = sorted(snapshot_ids & superseded)
        if condemned:
            raise DriftError(
                f"registry still selects {', '.join(condemned)}; the accepted correction "
                f"produced {refined} — update the registry validity to point at the "
                "refined snapshot"
            )
        _require_authorized_snapshots(snapshot_ids, {refined} if refined else set())
        return payload
    if verdict == "insufficient-data":
        raise DriftError("bulk processing refused: sampled verdict is insufficient-data")
    raise DriftError(f"bulk processing refused: sampled verdict is {verdict}")
