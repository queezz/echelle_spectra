"""Packet F6 — historical epochs reach the registry, and the date scan is bounded.

The import is exercised through the shipped ``echelle snapshot import-historical``
verb and the packaged binders, so the three real epochs the campaign has to reach
are the ones under test rather than a synthetic stand-in.
"""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import pytest

from echelle_spectra import cli
from echelle_spectra.calibration_registry import (
    load_calibration_registry,
    source_identity_from_path,
)
from echelle_spectra.historical import (
    HistoricalError,
    bundled_historical_manifests,
    import_historical_snapshot,
)
from echelle_spectra.snapshot import load_snapshot

# The 2025 campaign folder's sphere pair and Neon lamp are far too large to
# package, so the binder names them without an identity and the import captures
# it from whatever root the operator supplies.
SUPPLIED_2025_ARTIFACTS = (
    "sphere-0.1s-x3.sif",
    "sphere-0.1s-x3-bg.sif",
    "Ne-0.02s-x3-bright-lines.sif",
    "Ne-0.02s-x3-bright-lines_bg.sif",
)

EPOCHS = (
    ("20190529_cmos", "2019-05-29", "2024-03-04"),
    ("20240305_cmos", "2024-03-05", "2025-09-25"),
    ("20250926_cmos", "2025-09-26", None),
)


def _campaign_folder(tmp_path: Path) -> Path:
    """Stand in for local/20250926_calib, which never travels with the package."""

    folder = tmp_path / "20250926_calib"
    folder.mkdir()
    for name in SUPPLIED_2025_ARTIFACTS:
        (folder / name).write_bytes(f"campaign artifact {name}".encode())
    return folder


def _registry(path: Path, *snapshot_ids: str) -> Path:
    lines = ['schema = "echelle-calibration-registry/v1"']
    for snapshot_id in snapshot_ids:
        lines.extend(["", "[[epochs]]", f'snapshot_id = "{snapshot_id}"'])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Historical binders reach the registry
# ---------------------------------------------------------------------------


def test_the_2025_binder_names_the_campaign_sphere_pair_it_actually_used() -> None:
    """The binder states what was used, not the nearest file that happened to ship."""

    manifest = next(
        item for item in bundled_historical_manifests() if item.snapshot_id == "20250926_cmos"
    )
    sphere = manifest.artifact_for_role("sphere")
    background = manifest.artifact_for_role("sphere_background")
    assert (sphere.source_name, background.source_name) == (
        "sphere-0.1s-x3.sif",
        "sphere-0.1s-x3-bg.sif",
    )
    # No packaged root holds them, so their identity is honestly absent rather
    # than borrowed from the 2024 pair.
    assert not sphere.available and not background.available
    assert not sphere.identity_recorded
    assert {item.role for item in manifest.missing} == {"sphere", "sphere_background", "lamp"}


def test_every_bundled_epoch_converts_registers_and_resolves_a_dated_source(
    tmp_path: Path,
) -> None:
    calibrations = tmp_path / "calibrations"
    roots = (str(_campaign_folder(tmp_path)),)
    for snapshot_id, valid_from, valid_to in EPOCHS:
        snapshot = import_historical_snapshot(
            snapshot_id,
            calibrations,
            valid_from=valid_from,
            valid_to=valid_to,
            artifact_roots=roots,
        )
        assert snapshot.snapshot_id == snapshot_id
        assert snapshot.manifest["validity"]["date_from"] == valid_from
        # An imported epoch says so; a bench-measured one has no such table.
        assert snapshot.manifest["imported"]["binder_sha256"]
        # The snapshot, not the binder, is the identity authority from here on:
        # every artifact is digested and re-verified by the loader.
        load_snapshot(snapshot.root)

    registry = load_calibration_registry(
        _registry(tmp_path / "calibration_registry.toml", *(item[0] for item in EPOCHS)),
        snapshots_root=calibrations,
    )
    assert [epoch.snapshot_id for epoch in registry.epochs] == [item[0] for item in EPOCHS]

    shots = tmp_path / "nifs" / "2024-03-05"
    shots.mkdir(parents=True)
    source = shots / "193778_Echelle.SIF"
    source.write_bytes(b"raw")
    assert registry.resolve_source(source, root=shots).snapshot_id == "20240305_cmos"
    assert registry.resolve_source(source, root=shots).snapshot.detector == "cmos"


def test_the_2025_import_refuses_and_names_the_root_to_supply(tmp_path: Path) -> None:
    with pytest.raises(HistoricalError) as refusal:
        import_historical_snapshot("20250926_cmos", tmp_path / "calibrations")
    message = str(refusal.value)
    assert "sphere=sphere-0.1s-x3.sif" in message
    assert "--artifact-root DIR" in message
    assert "Searched:" in message
    assert not (tmp_path / "calibrations" / "20250926_cmos").exists()


def test_a_supplied_root_holding_a_different_file_cannot_masquerade(tmp_path: Path) -> None:
    """A recorded identity is verified wherever the artifact is found."""

    binder = tmp_path / "binder.toml"
    binder.write_text(
        "\n".join(
            [
                'schema = "echelle-historical-calibration/v1"',
                'snapshot_id = "20250926_cmos"',
                'detector = "cmos"',
                'acquired_date = "2025-09-26"',
                'lamps = ["Ne"]',
                *(
                    line
                    for role, name in (
                        ("pattern", "pattern.txt"),
                        ("wavelength", "wavelength.txt"),
                        ("sphere_background", "sphere-0.1s-x3-bg.sif"),
                        ("integral", "integral.txt"),
                    )
                    for line in (
                        "",
                        "[[artifacts]]",
                        f'role = "{role}"',
                        f'path = "{name}"',
                        f'source_name = "{name}"',
                    )
                ),
                "",
                "[[artifacts]]",
                'role = "sphere"',
                'path = "sphere-0.1s-x3.sif"',
                'source_name = "sphere-0.1s-x3.sif"',
                "size_bytes = 11",
                f'sha256 = "{"0" * 64}"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    supplied = tmp_path / "supplied"
    supplied.mkdir()
    for name in (
        "pattern.txt",
        "wavelength.txt",
        "sphere-0.1s-x3.sif",
        "sphere-0.1s-x3-bg.sif",
        "integral.txt",
    ):
        (supplied / name).write_text("substitute\n", encoding="utf-8")

    with pytest.raises(HistoricalError, match="identity mismatch"):
        import_historical_snapshot(
            str(binder), tmp_path / "calibrations", artifact_roots=(str(supplied),)
        )


def test_import_historical_is_reachable_through_the_umbrella_command(
    tmp_path: Path, capsys
) -> None:
    code = cli.main(
        [
            "snapshot",
            "import-historical",
            "20240305_cmos",
            "--calibrations",
            str(tmp_path / "calibrations"),
            "--valid-from",
            "2024-03-05",
            "--shot-from",
            "150000",
        ]
    )
    printed = capsys.readouterr().out
    assert code == 0
    assert "Imported historical calibration 20240305_cmos" in printed
    assert "shot_from=150000" in printed
    snapshot = load_snapshot(tmp_path / "calibrations" / "20240305_cmos")
    assert snapshot.manifest["validity"] == {
        "date_from": "2024-03-05",
        "shot_from": 150000,
    }
    assert snapshot.lamps == ("ThAr",)


def test_an_unknown_binder_name_lists_the_binders_that_do_exist(tmp_path: Path, capsys) -> None:
    code = cli.main(
        ["snapshot", "import-historical", "20260101_cmos", "--calibrations", str(tmp_path)]
    )
    assert code == 1
    assert "20250926_cmos" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# The date scan is bounded by the source root the operator named
# ---------------------------------------------------------------------------


def test_a_dated_volume_above_the_declared_root_supplies_no_acquisition_date(
    tmp_path: Path,
) -> None:
    root = tmp_path / "2019-05-29-archive-usb" / "nifs" / "shots"
    root.mkdir(parents=True)
    source = root / "193778_Echelle.SIF"
    source.write_bytes(b"raw")

    # Unbounded, the volume's own name dates the shot — the review's failure.
    assert source_identity_from_path(source).acquisition_date == date(2019, 5, 29)

    bounded = source_identity_from_path(source, root=root)
    assert bounded.acquisition_date is None
    assert bounded.shot_number == 193778


def test_the_declared_roots_own_name_still_dates_the_sources_below_it(tmp_path: Path) -> None:
    root = tmp_path / "2019-05-29-archive-usb" / "2024-03-05"
    root.mkdir(parents=True)
    source = root / "lamp.SIF"
    source.write_bytes(b"raw")

    identity = source_identity_from_path(source, root=root)
    assert identity.acquisition_date == date(2024, 3, 5)


def test_the_same_file_resolves_identically_named_absolutely_or_relatively(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "2024-03-05"
    root.mkdir()
    source = root / "193778_Echelle.SIF"
    source.write_bytes(b"raw")

    absolute = source_identity_from_path(source, root=root)
    monkeypatch.chdir(tmp_path)
    relative = source_identity_from_path(
        Path("2024-03-05") / "193778_Echelle.SIF", root=Path("2024-03-05")
    )
    assert relative.acquisition_date == absolute.acquisition_date == date(2024, 3, 5)
    assert relative.shot_number == absolute.shot_number == 193778

    # A working directory that is itself dated changes nothing either.
    dated_cwd = tmp_path / "2019-05-29"
    dated_cwd.mkdir()
    monkeypatch.chdir(dated_cwd)
    from_elsewhere = source_identity_from_path(
        Path(os.path.relpath(source, dated_cwd)), root=root
    )
    assert from_elsewhere.acquisition_date == date(2024, 3, 5)


def test_a_single_dated_file_is_gated_and_exported_under_the_same_epoch(
    tmp_path: Path,
) -> None:
    """The gate and the export must bound the date scan identically."""

    import sys
    from unittest.mock import MagicMock, patch

    import xarray as xr

    from echelle_spectra.drift import DRIFT_SCHEMA, write_drift_evidence
    from echelle_spectra.spectrocube_cli import ExportResult
    from echelle_spectra.spectrocube_cli import main as process_main

    calibrations = tmp_path / "calibrations"
    roots = (str(_campaign_folder(tmp_path)),)
    for snapshot_id, valid_from, valid_to in EPOCHS:
        import_historical_snapshot(
            snapshot_id,
            calibrations,
            valid_from=valid_from,
            valid_to=valid_to,
            artifact_roots=roots,
        )
    registry = _registry(
        tmp_path / "calibration_registry.toml", *(item[0] for item in EPOCHS)
    )
    evidence = write_drift_evidence(
        tmp_path / "drift.json",
        {
            "schema": DRIFT_SCHEMA,
            "verdict": "aligned",
            "snapshot_ids": ["20240305_cmos"],
            "summary": {"median_shift_px": 0.0},
        },
    )

    # The file carries no shot and no date of its own; only the folder the
    # operator named dates it.
    source = tmp_path / "2024-03-05" / "lamp.SIF"
    source.parent.mkdir()
    source.write_bytes(b"raw")

    def export(sif: Path, nc_out: Path, **kwargs) -> ExportResult:
        xr.Dataset(attrs={"snapshot_id": kwargs["snapshot"].snapshot_id}).to_netcdf(nc_out)
        return ExportResult("exported")

    output = tmp_path / "lamp.nc"
    with (
        patch.dict(sys.modules, {"spectrocube": MagicMock()}),
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch("echelle_spectra.spectrocube_cli._export_one", side_effect=export),
        pytest.raises(SystemExit) as result,
    ):
        process_main(
            [
                str(source),
                "-o",
                str(output),
                "--registry",
                str(registry),
                "--calibrations",
                str(calibrations),
                "--drift-verdict",
                str(evidence),
            ]
        )
    assert int(result.value.code) == 0
    with xr.open_dataset(output) as cube:
        assert cube.attrs["snapshot_id"] == "20240305_cmos"


def test_shot_digit_widths_are_read_exactly_as_documented() -> None:
    # Five to eight digits are LHD shots ...
    assert source_identity_from_path(Path("12345.SIF")).shot_number == 12345
    assert source_identity_from_path(Path("12345678_Echelle.SIF")).shot_number == 12345678
    # ... unless the eight digits are a calendar date, which the reachable guard
    # now actually rejects instead of being unreachable behind a seven-digit cap.
    dated = source_identity_from_path(Path("20240305_Echelle.SIF"))
    assert dated.shot_number is None
    assert dated.acquisition_date == date(2024, 3, 5)
    # Fewer than five leading digits is not an LHD shot and is not read as one;
    # an explicit token still is.
    assert source_identity_from_path(Path("1234.SIF")).shot_number is None
    assert source_identity_from_path(Path("shot_42.SIF")).shot_number == 42
