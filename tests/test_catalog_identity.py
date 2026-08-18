"""Packet F6 — a drive keeps one identity, and merges never lose fresh rows.

The processing paths run through the shipped ``echelle process`` entry point
with only the SIF reader and the exporter replaced, so drive identity, the
per-drive catalog, and the optional auto-merge are exercised where an operator
meets them.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from echelle_spectra.campaign_run import (
    DRIVE_ID_NAME,
    GATE_SAMPLE,
    GATE_UNGATED,
    GATE_UNRECORDED,
    GATE_VERDICT,
    ReceiptLookupError,
    RunReceipt,
    default_volume_label,
    ensure_drive_identity,
    read_drive_identity,
)
from echelle_spectra.campaign_tools_cli import catalog_main
from echelle_spectra.catalog import (
    CATALOG_NAME,
    build_drive_catalog,
    discover_drive_id,
    load_catalog,
    merge_catalogs,
    refresh_catalog_row,
    source_catalog_path,
)
from echelle_spectra.drift import (
    DATE_ATTRIBUTES,
    DRIFT_SCHEMA,
    cube_date,
    write_drift_evidence,
)
from echelle_spectra.recalibration import recalibrate_cube
from echelle_spectra.snapshot import ROLE_FILENAMES, create_snapshot
from echelle_spectra.spectrocube_cli import ExportResult, main


@pytest.fixture()
def fake_spectrocube():
    with patch.dict(sys.modules, {"spectrocube": MagicMock()}):
        yield


def _cube_writer(paths: list[Path]):
    """Write a catalogable cube carrying the attributes the CLI selected for it."""

    def export(sif: Path, nc_out: Path, **kwargs) -> ExportResult:
        paths.append(nc_out)
        xr.Dataset(
            {"intensity": (("wavelength",), np.arange(4.0))},
            coords={"wavelength": np.arange(410.0, 414.0)},
            attrs={"shot_number": sif.stem, **dict(kwargs["extra_attrs"])},
        ).to_netcdf(nc_out)
        return ExportResult("exported")

    return export


def _sparing_writer():
    """A fake exporter that keeps the real one's promise not to rewrite a cube.

    The plain writer above always writes, which would let a second run claim
    every cube the first one published.  The trip loop's second pass meets cubes
    a sampling run already wrote, and the real exporter skips them, so a test
    about who published what has to skip them too.
    """

    write = _cube_writer([])

    def export(sif: Path, nc_out: Path, **kwargs) -> ExportResult:
        if nc_out.exists() and not kwargs.get("overwrite"):
            return ExportResult("skipped", "output exists; use --overwrite to replace it")
        return write(sif, nc_out, **kwargs)

    return export


def _process(argv: list[str], *, export=None) -> int:
    with (
        patch("echelle_spectra.tools.loader.build_calibration", return_value=object()),
        patch(
            "echelle_spectra.spectrocube_cli._export_one",
            side_effect=export or _cube_writer([]),
        ),
        pytest.raises(SystemExit) as result,
    ):
        main(argv)
    return int(result.value.code)


def _drive(root: Path, count: int = 2) -> Path:
    shots = root / "shots"
    shots.mkdir(parents=True)
    for index in range(count):
        (shots / f"20000{index}_Echelle.SIF").write_text(f"raw {index}", encoding="utf-8")
    return shots


def _run_argv(drive: Path, tmp_path: Path, *extra: str) -> list[str]:
    return [
        str(drive / "shots"),
        "-o",
        str(drive / "cubes"),
        "--runs-dir",
        str(tmp_path / "runs" / drive.name),
        *extra,
    ]


# ---------------------------------------------------------------------------
# One drive, one identity
# ---------------------------------------------------------------------------


def test_first_processing_writes_the_drive_id_and_the_catalog_keys_on_it(
    tmp_path: Path, fake_spectrocube
) -> None:
    drive = tmp_path / "drive-a"
    _drive(drive)
    assert _process(_run_argv(drive, tmp_path, "--volume-label", "NIFS-A")) == 0

    marker = drive / "shots" / DRIVE_ID_NAME
    assert marker.is_file()
    identity = read_drive_identity(drive / "shots")
    assert identity is not None and len(identity.drive_id) == 32
    assert identity.label == "NIFS-A"

    catalog = load_catalog(drive / "cubes" / CATALOG_NAME)
    assert catalog["drive_id"] == identity.drive_id
    assert catalog["volume_label"] == "NIFS-A"

    # A second run reuses the announced identity rather than minting a new one.
    before = marker.read_text(encoding="utf-8")
    assert _process(_run_argv(drive, tmp_path, "--volume-label", "NIFS-A")) == 0
    assert marker.read_text(encoding="utf-8") == before


def test_a_reconnected_drive_under_another_root_keeps_one_catalog_identity(
    tmp_path: Path, fake_spectrocube
) -> None:
    first = tmp_path / "E-drive"
    _drive(first)
    assert _process(_run_argv(first, tmp_path, "--volume-label", "NIFS-A")) == 0
    original = load_catalog(first / "cubes" / CATALOG_NAME)

    # The same disk returns under a different letter or mount point.
    second = tmp_path / "Volumes" / "NIFS-A"
    second.parent.mkdir()
    shutil.copytree(first, second)
    shutil.rmtree(second / "cubes")
    assert _process(_run_argv(second, tmp_path, "--volume-label", "NIFS-A")) == 0
    reconnected = load_catalog(second / "cubes" / CATALOG_NAME)
    assert reconnected["drive_id"] == original["drive_id"]

    merged = load_catalog(
        merge_catalogs(
            [first / "cubes" / CATALOG_NAME, second / "cubes" / CATALOG_NAME],
            tmp_path / "all-years.json",
        )
    )
    assert len(merged["sources"]) == 1
    assert merged["sources"][0]["drive_id"] == original["drive_id"]


def test_a_root_that_cannot_store_the_id_falls_back_and_says_so_in_the_receipt(
    tmp_path: Path, fake_spectrocube, monkeypatch: pytest.MonkeyPatch, capsys
) -> None:
    """Simulate a read-only mount by refusing exactly the drive-id write."""

    original_write_text = Path.write_text

    def refuse(self: Path, *args, **kwargs):
        if DRIVE_ID_NAME in self.name:
            raise OSError(30, "Read-only file system")
        return original_write_text(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", refuse)
    drive = tmp_path / "archive"
    _drive(drive)
    assert _process(_run_argv(drive, tmp_path, "--volume-label", "ARCHIVE-2019")) == 0
    assert "not writable" in capsys.readouterr().err

    monkeypatch.undo()
    assert not (drive / "shots" / DRIVE_ID_NAME).exists()
    receipt = next(RunReceipt.load(path.parent) for path in (tmp_path / "runs").rglob("run.toml"))
    assert receipt.drive_id == ""
    assert "volume label alone" in receipt.drive_warning
    assert receipt.volume_label == "ARCHIVE-2019"
    catalog = load_catalog(drive / "cubes" / CATALOG_NAME)
    assert catalog["drive_id"] == ""
    assert "volume label alone" in catalog["run"]["drive_warning"]
    # Unidentified drives still merge, keyed on the only identity they have.
    merged = load_catalog(
        merge_catalogs([drive / "cubes" / CATALOG_NAME], tmp_path / "all-years.json")
    )
    assert merged["sources"][0]["volume_label"] == "ARCHIVE-2019"


def test_the_fallback_label_never_degrades_to_an_empty_or_anchor_only_name(
    tmp_path: Path,
) -> None:
    labelled = tmp_path / "NIFS-A" / "shots"
    labelled.mkdir(parents=True)
    label = default_volume_label(labelled)
    assert label
    # On Windows the anchor is the drive letter; on POSIX it names nothing, so
    # the named folder is used instead of a silent blank.
    assert label == labelled.anchor.rstrip("\\/") or label == "shots"
    assert default_volume_label(Path(labelled.anchor or "/")) != ""

    identity = ensure_drive_identity(labelled)
    assert identity.label == label
    assert identity.stored is True


# ---------------------------------------------------------------------------
# Merge integrity
# ---------------------------------------------------------------------------


def _drive_catalog(root: Path, *, label: str, drive_id: str, generated_at: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    path = root / CATALOG_NAME
    path.write_text(
        json.dumps(
            {
                "schema": "echelle-drive-catalog/v1",
                "generated_at": generated_at,
                "volume_label": label,
                "drive_id": drive_id,
                "cube_root": ".",
                "run": {"state": "completed", "gate": GATE_UNGATED},
                "cubes": [{"path": f"{generated_at[:10]}.nc", "snapshot_id": "20250101_cmos"}],
                "errors": [],
            }
        ),
        encoding="utf-8",
    )
    return path


def test_a_stale_index_merged_after_a_fresh_catalog_cannot_revert_its_rows(
    tmp_path: Path,
) -> None:
    drive = tmp_path / "drive"
    stale_catalog = _drive_catalog(
        drive, label="NIFS-A", drive_id="deadbeef", generated_at="2026-01-01T00:00:00.000+00:00"
    )
    stale_index = merge_catalogs([stale_catalog], tmp_path / "stale-index.json")
    fresh_catalog = _drive_catalog(
        drive, label="NIFS-A", drive_id="deadbeef", generated_at="2026-08-14T00:00:00.000+00:00"
    )

    # The review's exact footgun: the stale all-years index is merged last.
    merged = load_catalog(
        merge_catalogs([fresh_catalog, stale_index], tmp_path / "all-years.json")
    )
    assert len(merged["sources"]) == 1
    source = merged["sources"][0]
    assert source["generated_at"] == "2026-08-14T00:00:00.000+00:00"
    assert [row["path"] for row in source["cubes"]] == ["2026-08-14.nc"]

    # Order genuinely does not decide it: the other order agrees.
    other = load_catalog(
        merge_catalogs([stale_index, fresh_catalog], tmp_path / "other-order.json")
    )
    assert other["sources"][0]["cubes"] == source["cubes"]


def test_rows_store_the_catalog_path_relative_to_a_relocatable_drive_root(
    tmp_path: Path,
) -> None:
    drive = tmp_path / "E-drive" / "cubes"
    catalog = _drive_catalog(
        drive, label="NIFS-A", drive_id="deadbeef", generated_at="2026-08-14T00:00:00.000+00:00"
    )
    index = merge_catalogs([catalog], tmp_path / "all-years.json")
    source = load_catalog(index)["sources"][0]
    # Portable across operating systems: no drive letter, no absolute prefix.
    assert source["catalog_path"] == CATALOG_NAME
    assert source["available"] is True

    # The drive comes back at a different mount point.
    relocated = tmp_path / "Volumes" / "NIFS-A" / "cubes"
    relocated.parent.mkdir(parents=True)
    shutil.move(str(drive), str(relocated))
    assert source_catalog_path(load_catalog(index)["sources"][0]).is_file() is False

    remerged = load_catalog(
        merge_catalogs([index, relocated / CATALOG_NAME], tmp_path / "all-years.json")
    )
    assert len(remerged["sources"]) == 1
    assert remerged["sources"][0]["available"] is True
    assert remerged["sources"][0]["catalog_path"] == CATALOG_NAME
    assert Path(remerged["sources"][0]["drive_root"]) == relocated.resolve()


# ---------------------------------------------------------------------------
# Auto-merge and the gate a catalog row remembers
# ---------------------------------------------------------------------------


def test_a_completed_run_auto_merges_only_when_a_central_index_is_configured(
    tmp_path: Path, fake_spectrocube
) -> None:
    quiet = tmp_path / "quiet"
    _drive(quiet)
    assert _process(_run_argv(quiet, tmp_path, "--volume-label", "QUIET")) == 0
    assert not (tmp_path / "all-years.json").exists()

    loud = tmp_path / "loud"
    _drive(loud)
    index = tmp_path / "all-years.json"
    assert (
        _process(_run_argv(loud, tmp_path, "--volume-label", "LOUD", "--central-index", str(index)))
        == 0
    )
    merged = load_catalog(index)
    assert [source["volume_label"] for source in merged["sources"]] == ["LOUD"]
    assert len(merged["sources"][0]["cubes"]) == 2

    # The manual flow still folds the quiet drive in beside it.
    both = load_catalog(merge_catalogs([index, quiet / "cubes" / CATALOG_NAME], index))
    assert sorted(source["volume_label"] for source in both["sources"]) == ["LOUD", "QUIET"]


def test_catalog_rows_record_the_gate_that_authorized_the_run(
    tmp_path: Path, fake_spectrocube
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    drive = tmp_path / "drive"
    _drive(drive)
    assert (
        _process(
            _run_argv(
                drive,
                tmp_path,
                "--registry",
                str(registry),
                "--calibrations",
                str(snapshots),
                "--sample",
                "1",
            )
        )
        == 0
    )
    catalog = load_catalog(drive / "cubes" / CATALOG_NAME)
    assert catalog["run"]["gate"] == GATE_SAMPLE
    assert catalog["run"]["sample"] is True
    assert [row["gate"] for row in catalog["cubes"]] == [GATE_SAMPLE]

    # A cube this run never published keeps the unrecorded word instead of
    # borrowing the sample's authorization.
    xr.Dataset(
        {"intensity": (("wavelength",), np.arange(2.0))},
        coords={"wavelength": np.arange(410.0, 412.0)},
    ).to_netcdf(drive / "cubes" / "stranger.nc")
    receipt_dir = next((tmp_path / "runs").rglob("run.toml")).parent
    rebuilt = load_catalog(
        build_drive_catalog(
            drive / "cubes", volume_label="NIFS-A", receipt_dir=receipt_dir
        )
    )
    gates = {row["path"]: row["gate"] for row in rebuilt["cubes"]}
    assert gates["stranger.nc"] == GATE_UNRECORDED
    assert set(gates.values()) == {GATE_SAMPLE, GATE_UNRECORDED}


def _epoch_registry(tmp_path: Path) -> tuple[Path, Path]:
    snapshots = tmp_path / "calibrations"
    sources = tmp_path / "snapshot-sources"
    sources.mkdir()
    files: dict[str, Path] = {}
    for role in ROLE_FILENAMES:
        source = sources / f"{role}.dat"
        source.write_text(f"20250101_cmos/{role}\n", encoding="utf-8")
        files[role] = source
    create_snapshot(
        snapshots,
        snapshot_id="20250101_cmos",
        detector="cmos",
        files=files,
        lamps=("ThAr",),
        validity={"shot_from": 100000, "shot_to": 299999},
    )
    registry = tmp_path / "calibration_registry.toml"
    registry.write_text(
        'schema = "echelle-calibration-registry/v1"\n\n'
        '[[epochs]]\nsnapshot_id = "20250101_cmos"\n',
        encoding="utf-8",
    )
    return registry, snapshots


# ---------------------------------------------------------------------------
# The input/output seam: the rehearsal sequence, typed with the ROOT paths
#
# The drive announces itself in the tree that is *read* and the cubes are
# written somewhere else entirely, so a catalog that searched around the cubes
# found no drive id at all; and `--receipt-dir` was given the runs root, which
# holds no run.toml of its own, so the run was silently not identified either.
# The first real end-to-end run therefore produced a catalog with an empty
# drive_id, `run: null`, and every cube "unrecorded (pre-gate receipt)" -- while
# `echelle status`, reading that very tree, reported the run completed.
# ---------------------------------------------------------------------------


def test_the_rehearsal_sequence_with_root_paths_fully_identifies_the_catalog(
    tmp_path: Path, fake_spectrocube
) -> None:
    registry, snapshots = _epoch_registry(tmp_path)
    drive = tmp_path / "NIFS"
    _drive(drive, count=3)
    runs = tmp_path / "runs"  # the --runs-dir ROOT, exactly as the loop prints it
    cubes = drive / "cubes"
    gated = [
        str(drive / "shots"),
        "-o",
        str(cubes),
        "--runs-dir",
        str(runs),
        "--registry",
        str(registry),
        "--calibrations",
        str(snapshots),
        "--volume-label",
        "NIFS-A",
    ]

    # Steps 4-5: the unverified first sample. Then 6-7: the audited verdict
    # spends that sample and takes the rest of the drive.
    assert _process([*gated, "--sample", "1"], export=_sparing_writer()) == 0
    evidence = write_drift_evidence(
        tmp_path / "drift-evidence-001.json",
        {"schema": DRIFT_SCHEMA, "verdict": "aligned", "snapshot_ids": ["20250101_cmos"]},
    )
    assert _process([*gated, "--drift-verdict", str(evidence)], export=_sparing_writer()) == 0

    # The old empty-identity outcome, reproduced before the fix is asked for.
    assert (drive / "shots" / DRIVE_ID_NAME).is_file()
    assert discover_drive_id(cubes) == ""
    assert not (runs / "run.toml").exists()

    # Step 8, typed with the roots.
    catalog = load_catalog(
        build_drive_catalog(
            cubes,
            volume_label="NIFS-A",
            receipt_dir=runs,
            output=cubes / "catalog.json",
        )
    )

    announced = read_drive_identity(drive / "shots")
    assert announced is not None
    assert catalog["drive_id"] == announced.drive_id != ""
    assert catalog["run"] is not None
    assert catalog["run"]["state"] == "completed"
    assert catalog["run"]["gate"] == GATE_VERDICT
    # Nothing is left unattributed: the sampled cube is identified by the
    # sample's own receipt rather than condemned because the later run only
    # skipped it.
    assert sorted(row["gate"] for row in catalog["cubes"]) == [
        GATE_SAMPLE,
        GATE_VERDICT,
        GATE_VERDICT,
    ]
    assert GATE_UNRECORDED not in {row["gate"] for row in catalog["cubes"]}

    # The receipts themselves say which drive they belong to, not only which
    # folder: two drives ending in the same folder name stay tellable apart.
    assert all(
        "nifs-a-shots" in manifest.parent.name for manifest in runs.rglob("run.toml")
    )


def test_a_receipt_folder_that_identifies_nothing_is_refused_by_name(
    tmp_path: Path,
) -> None:
    cubes = tmp_path / "cubes"
    cubes.mkdir()
    empty = tmp_path / "not-runs"
    empty.mkdir()

    with pytest.raises(ReceiptLookupError) as absent:
        build_drive_catalog(cubes, volume_label="NIFS-A", receipt_dir=tmp_path / "gone")
    assert "--runs-dir" in str(absent.value)

    with pytest.raises(ReceiptLookupError) as barren:
        build_drive_catalog(cubes, volume_label="NIFS-A", receipt_dir=empty)
    assert "no run receipt under" in str(barren.value)
    assert "--runs-dir" in str(barren.value)


def test_a_runs_root_holding_another_drives_receipts_is_refused_in_one_line(
    tmp_path: Path, fake_spectrocube, capsys
) -> None:
    drive = tmp_path / "drive-a"
    _drive(drive)
    runs = tmp_path / "runs"
    assert (
        _process(
            [
                str(drive / "shots"),
                "-o",
                str(drive / "cubes"),
                "--runs-dir",
                str(runs),
                "--volume-label",
                "NIFS-A",
            ]
        )
        == 0
    )
    capsys.readouterr()

    stranger = tmp_path / "elsewhere"
    stranger.mkdir()
    code = catalog_main(
        [
            "build",
            str(stranger),
            "--volume-label",
            "NIFS-B",
            "--receipt-dir",
            str(runs),
        ]
    )

    assert code == 1
    refusal = capsys.readouterr().err
    assert refusal.startswith("ERROR: ")
    assert refusal.count("\n") == 1
    assert "none of which wrote cubes to" in refusal
    assert "(from --receipt-dir)" in refusal
    assert not (stranger / CATALOG_NAME).exists()


# ---------------------------------------------------------------------------
# Recalibration refreshes what the catalog claims
# ---------------------------------------------------------------------------


def _recal_snapshot(root: Path, snapshot_id: str, *, shift: float = 0.0):
    sources = root / f"sources-{snapshot_id}"
    sources.mkdir(parents=True)
    files = {}
    for role, name in {
        "pattern": "pattern.txt",
        "sphere": "sphere.sif",
        "sphere_background": "sphere_bg.sif",
        "integral": "integral.txt",
    }.items():
        path = sources / name
        path.write_text("same" if role == "pattern" else "fixture", encoding="utf-8")
        files[role] = path
    wavelength = sources / "wavelength.txt"
    wavelength.write_text(
        "# fixture\n# order from to center wavelength\n"
        + "\n".join(
            f"0 {index} {index + 1} {index} {410.0 + index + shift} line" for index in range(5)
        )
        + "\n",
        encoding="utf-8",
    )
    files["wavelength"] = wavelength
    return create_snapshot(
        root / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=["ThAr"],
        validity={"shot_from": 1, "shot_to": 999999},
    )


def _provenance_cube(snapshot) -> xr.Dataset:
    wavelength = np.arange(410.0, 415.0)
    factor = np.full(wavelength.size, 2.0)
    ds = xr.Dataset(
        {"intensity": (("frame", "wavelength"), (np.arange(1.0, 6.0) * factor)[None, :])},
        coords={
            "frame": [0],
            "wavelength": wavelength,
            "detector_pixel": ("wavelength", np.arange(wavelength.size, dtype=float)),
            "echelle_order": ("wavelength", np.zeros(wavelength.size, dtype=np.int64)),
        },
        attrs={
            "spectrocube_version": "0.2.0",
            "instrument_id": "echelle",
            "calibration_type": "absolute",
            "intensity_units": "W/m2/nm/sr",
            "wavelength_medium": "air",
            "shot_number": "12345",
            "created_at": "2025-01-01T00:00:00+00:00",
            **snapshot.provenance_attrs(),
        },
    )
    ds.coords["wavelength"].attrs.update(units="nm", medium="air")
    ds.coords["detector_pixel"].attrs.update(
        units="pixel", detector_axis="column", reference_frame="raw_detector", index_origin=0
    )
    ds["applied_absolute_calibration_factor"] = (
        ("wavelength",),
        factor,
        {
            "units": "W/m2/nm/sr per (counts/s)",
            "source_units": "counts/s",
            "absolute_kind": "wmsr",
            "application": (
                "stored_intensity = source_signal * applied_absolute_calibration_factor"
            ),
        },
    )
    return ds


def test_recal_refreshes_the_adjacent_catalog_row_and_bumps_the_drive(
    tmp_path: Path,
) -> None:
    pytest.importorskip("spectrocube", reason="spectrocube package not installed")
    old = _recal_snapshot(tmp_path, "20240101_cmos")
    new = _recal_snapshot(tmp_path, "20250101_cmos", shift=0.1)
    cubes = tmp_path / "cubes"
    cubes.mkdir()
    source = cubes / "12345.nc"
    _provenance_cube(old).to_netcdf(source)

    catalog_path = build_drive_catalog(cubes, volume_label="NIFS-A", drive_id="deadbeef")
    before = load_catalog(catalog_path)
    original = next(row for row in before["cubes"] if row["path"] == "12345.nc")
    index = merge_catalogs([catalog_path], tmp_path / "all-years.json")

    output = cubes / "12345-r1.nc"
    recalibrate_cube(source, output, new_snapshot_path=new.root, update_factor=False)

    after = load_catalog(catalog_path)
    revised = next(row for row in after["cubes"] if row["path"] == "12345-r1.nc")
    assert revised["snapshot_id"] == "20250101_cmos"
    assert revised["sha256"] != original["sha256"]
    # The untouched original row is preserved beside the revision.
    assert {row["path"] for row in after["cubes"]} == {"12345.nc", "12345-r1.nc"}
    # The drive is now newer than what the index remembers, so a later merge
    # propagates the correction instead of the index winning on list order.
    assert after["generated_at"] > before["generated_at"]
    propagated = load_catalog(merge_catalogs([index, catalog_path], index))
    assert {row["path"] for row in propagated["sources"][0]["cubes"]} == {
        "12345.nc",
        "12345-r1.nc",
    }


def test_recal_outside_any_catalog_is_silent_rather_than_inventing_one(
    tmp_path: Path,
) -> None:
    cubes = tmp_path / "loose"
    cubes.mkdir()
    (cubes / "12345.nc").write_bytes(b"not a cube")
    assert refresh_catalog_row(cubes / "12345.nc") is None
    assert not (cubes / CATALOG_NAME).exists()


# ---------------------------------------------------------------------------
# t_start makes date selection exact
# ---------------------------------------------------------------------------


def test_export_writes_t_start_from_the_sif_acquisition_time() -> None:
    pytest.importorskip("spectrocube", reason="spectrocube package not installed")
    from echelle_spectra.tools.spectrocube_export import to_spectrocube

    class _Spectrum:
        pass

    spectrum = _Spectrum()
    spectrum.wavelength = np.linspace(400.0, 700.0, 8)
    spectrum.counts = np.arange(8.0)[None, :]
    # The Unix timestamp shape an Andor header records in ExperimentTime.
    spectrum.info = {"NumberOfFrames": 1, "ExposureTime": 0.5, "ExperimentTime": 1540951489}
    spectrum.fpth = "/data/193778_Echelle.sif"

    cube = to_spectrocube(spectrum, units="counts", squeeze_single_frame=False)
    # Stored in UTC: a cube travels between machines, and a naive local time
    # would silently mean a different instant on each of them.
    assert cube.ds.attrs["t_start"] == "2018-10-31T02:04:49+00:00"

    spectrum.info.pop("ExperimentTime")
    assert "t_start" not in to_spectrocube(spectrum, units="counts").ds.attrs


def test_drift_date_selection_prefers_t_start_over_the_source_filename(
    tmp_path: Path,
) -> None:
    from echelle_spectra.drift import _filter_by_date

    assert DATE_ATTRIBUTES[0] == "t_start"
    attrs = {
        "t_start": "2025-09-26T01:02:03+00:00",
        "source_file": "/archive/2019-05-29/193778_Echelle.sif",
        "created_at": "2030-01-01T00:00:00+00:00",
    }
    assert cube_date(attrs) == (__import__("datetime").date(2025, 9, 26), "t_start")

    for name, cube_attrs in (
        ("dated.nc", attrs),
        ("older.nc", {"source_file": "/archive/2019-05-29/x.sif"}),
    ):
        xr.Dataset(attrs=cube_attrs).to_netcdf(tmp_path / name)
    selected = _filter_by_date(
        sorted(tmp_path.glob("*.nc")), date_from="2025-01-01", date_to="2025-12-31"
    )
    assert [path.name for path in selected] == ["dated.nc"]
