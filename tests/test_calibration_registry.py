from __future__ import annotations

from datetime import date
from pathlib import Path

import pytest

from echelle_spectra.calibration_registry import (
    REGISTRY_SCHEMA,
    CalibrationRegistryError,
    CalibrationSourceIdentity,
    load_calibration_registry,
    source_identity_from_path,
)
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot


def _snapshot(
    root: Path,
    snapshot_id: str,
    *,
    validity: dict[str, object],
) -> Path:
    sources = root.parent / f"{snapshot_id}-sources"
    sources.mkdir()
    files: dict[str, Path] = {}
    for role in ROLE_FILENAMES:
        source = sources / f"{role}.dat"
        source.write_text(f"{snapshot_id}/{role}\n", encoding="utf-8")
        files[role] = source
    return create_snapshot(
        root,
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity=validity,
    ).root


def _registry(path: Path, *snapshot_ids: str) -> Path:
    lines = [f'schema = "{REGISTRY_SCHEMA}"']
    for snapshot_id in snapshot_ids:
        lines.extend(["", "[[epochs]]", f'snapshot_id = "{snapshot_id}"'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _identity(*, shot: int | None = None, day: str | None = None) -> CalibrationSourceIdentity:
    return CalibrationSourceIdentity(
        path=Path("fixture.SIF"),
        shot_number=shot,
        acquisition_date=date.fromisoformat(day) if day else None,
    )


def test_shot_boundaries_are_inclusive_and_select_both_sides(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(
        snapshots,
        "20240101_cmos",
        validity={"shot_from": 100, "shot_to": 199},
    )
    _snapshot(
        snapshots,
        "20250101_cmos",
        validity={"shot_from": 200, "shot_to": 299},
    )
    registry = load_calibration_registry(
        _registry(tmp_path / "calibration_registry.toml", "20240101_cmos", "20250101_cmos"),
        snapshots_root=snapshots,
    )

    assert registry.resolve(_identity(shot=199)).snapshot_id == "20240101_cmos"
    assert registry.resolve(_identity(shot=200)).snapshot_id == "20250101_cmos"
    assert registry.resolve(_identity(shot=201)).snapshot_id == "20250101_cmos"


def test_date_boundaries_are_inclusive(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(
        snapshots,
        "20240101_cmos",
        validity={"date_from": "2024-01-01", "date_to": "2024-12-31"},
    )
    _snapshot(
        snapshots,
        "20250101_cmos",
        validity={"date_from": "2025-01-01", "date_to": "2025-12-31"},
    )
    registry = load_calibration_registry(
        _registry(tmp_path / "calibration_registry.toml", "20250101_cmos", "20240101_cmos"),
        snapshots_root=snapshots,
    )

    assert registry.resolve(_identity(day="2024-12-30")).snapshot_id == "20240101_cmos"
    assert registry.resolve(_identity(day="2024-12-31")).snapshot_id == "20240101_cmos"
    assert registry.resolve(_identity(day="2025-01-01")).snapshot_id == "20250101_cmos"


def test_shot_and_date_bounds_are_an_explicit_and_match(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(
        snapshots,
        "20250101_cmos",
        validity={
            "shot_from": 200,
            "shot_to": 299,
            "date_from": "2025-01-01",
            "date_to": "2025-12-31",
        },
    )
    registry = load_calibration_registry(
        _registry(tmp_path / "registry.toml", "20250101_cmos"),
        snapshots_root=snapshots,
    )

    assert registry.resolve(_identity(shot=250, day="2025-06-01")).snapshot_id == "20250101_cmos"
    with pytest.raises(CalibrationRegistryError, match="no calibration epoch matches"):
        registry.resolve(_identity(shot=199, day="2025-06-01"))
    with pytest.raises(CalibrationRegistryError, match="missing required.*acquisition date"):
        registry.resolve(_identity(shot=250))


def test_registry_preserves_declared_order_without_mapping_iteration(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(snapshots, "20240101_cmos", validity={"shot_from": 1, "shot_to": 9})
    _snapshot(snapshots, "20250101_cmos", validity={"shot_from": 10, "shot_to": 19})
    registry = load_calibration_registry(
        _registry(tmp_path / "registry.toml", "20250101_cmos", "20240101_cmos"),
        snapshots_root=snapshots,
    )
    assert [epoch.snapshot_id for epoch in registry.epochs] == [
        "20250101_cmos",
        "20240101_cmos",
    ]
    assert [epoch.position for epoch in registry.epochs] == [1, 2]


def test_overlapping_inclusive_boundaries_are_rejected(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(snapshots, "20240101_cmos", validity={"shot_from": 1, "shot_to": 10})
    _snapshot(snapshots, "20250101_cmos", validity={"shot_from": 10, "shot_to": 20})
    path = _registry(tmp_path / "registry.toml", "20240101_cmos", "20250101_cmos")
    with pytest.raises(CalibrationRegistryError, match="overlap is ambiguous.*inclusive"):
        load_calibration_registry(path, snapshots_root=snapshots)


def test_no_match_and_missing_identity_fail_clearly(tmp_path: Path) -> None:
    snapshots = tmp_path / "calibrations"
    _snapshot(snapshots, "20240101_cmos", validity={"shot_from": 100, "shot_to": 199})
    registry = load_calibration_registry(
        _registry(tmp_path / "registry.toml", "20240101_cmos"), snapshots_root=snapshots
    )
    with pytest.raises(CalibrationRegistryError, match="no calibration epoch matches"):
        registry.resolve(_identity(shot=99))
    with pytest.raises(CalibrationRegistryError, match="missing required.*shot number"):
        registry.resolve(_identity())


@pytest.mark.parametrize("defect", ["missing", "invalid", "digest"])
def test_missing_invalid_or_digest_mismatched_snapshot_is_rejected(
    tmp_path: Path, defect: str
) -> None:
    snapshots = tmp_path / "calibrations"
    snapshot_id = "20240101_cmos"
    if defect != "missing":
        snapshot = _snapshot(
            snapshots,
            snapshot_id,
            validity={"shot_from": 1, "shot_to": 9},
        )
        if defect == "invalid":
            (snapshot / "snapshot.toml").write_text("schema = 'broken'\n", encoding="utf-8")
        else:
            (snapshot / "pattern.txt").write_text("tampered\n", encoding="utf-8")
    path = _registry(tmp_path / "registry.toml", snapshot_id)
    with pytest.raises(CalibrationRegistryError, match="snapshot.*invalid"):
        load_calibration_registry(path, snapshots_root=snapshots)


def test_path_parser_supports_lhd_shots_and_dates_without_confusing_them() -> None:
    shot = source_identity_from_path(Path("NIFS") / "193778_Echelle.SIF")
    assert shot.shot_number == 193778
    assert shot.acquisition_date is None

    dated = source_identity_from_path(Path("2025-09-26") / "lamp.SIF")
    assert dated.shot_number is None
    assert dated.acquisition_date == date(2025, 9, 26)

    both = source_identity_from_path(Path("2025-09-26") / "shot_193778.SIF")
    assert both.shot_number == 193778
    assert both.acquisition_date == date(2025, 9, 26)
