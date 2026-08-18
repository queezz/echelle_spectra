"""A snapshot that points at the light instead of copying it.

The calibration folder holds the lamp and sphere frames, and it holds them
completely.  A snapshot written beside them keeps everything it computed and
records where the raw frames are, with the digest that says they are still the
same frames.  These tests pin that: no raw bytes inside, the reference resolved
from the snapshot folder rather than the working directory, a clear refusal
naming the absolute path when a referenced frame goes missing or changes — and
the older all-copied snapshots still loading exactly as they did.
"""

from __future__ import annotations

import json
import os
from datetime import date
from pathlib import Path

import pytest

from echelle_spectra.calibration_registry import (
    CalibrationSourceIdentity,
    load_calibration_registry,
)
from echelle_spectra.snapshot import (
    ARTIFACT_COPIED,
    ARTIFACT_REFERENCED,
    Snapshot,
    SnapshotValidationError,
    create_snapshot,
    load_snapshot,
    read_snapshot_folder,
    saved_snapshots_in,
)

RAW_SPHERE = b"sphere pixels" * 16
RAW_SPHERE_BACKGROUND = b"sphere background pixels" * 16
RAW_LAMP = b"thar pixels" * 16


@pytest.fixture
def calibration_folder(tmp_path: Path) -> Path:
    """A calibration folder shaped the way the bench leaves one behind."""

    (tmp_path / "sphere-0.1s-x3.sif").write_bytes(RAW_SPHERE)
    (tmp_path / "sphere-0.1s-x3-bg.sif").write_bytes(RAW_SPHERE_BACKGROUND)
    (tmp_path / "ThAr-0.3s-x3.sif").write_bytes(RAW_LAMP)
    for name in ("pattern.txt", "wavelength.txt", "integral.txt"):
        (tmp_path / name).write_text(f"{name} rows\n", encoding="utf-8")
    return tmp_path


def _sources(folder: Path) -> dict[str, Path]:
    return {
        "pattern": folder / "pattern.txt",
        "wavelength": folder / "wavelength.txt",
        "integral": folder / "integral.txt",
        "sphere": folder / "sphere-0.1s-x3.sif",
        "sphere_background": folder / "sphere-0.1s-x3-bg.sif",
    }


def _thin(
    folder: Path,
    *,
    snapshot_id: str = "20260813_cmos",
    lamps: tuple[str, ...] = ("ThAr",),
    **extra,
):
    return create_snapshot(
        folder / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=_sources(folder),
        lamps=lamps,
        lamp_files=[("ThAr", folder / "ThAr-0.3s-x3.sif")],
        reference_raw=True,
        **extra,
    )


def _restamp(snapshot: Snapshot, created_utc: str) -> None:
    """Rewrite one binder's save moment; the artifact digests are untouched."""

    manifest = snapshot.root / "snapshot.toml"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    manifest.write_text(
        "\n".join(
            f'created_utc = "{created_utc}"' if line.startswith("created_utc = ") else line
            for line in lines
        )
        + "\n",
        encoding="utf-8",
    )


def _fat(folder: Path, *, snapshot_id: str = "20260813_cmos"):
    return create_snapshot(
        folder / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=_sources(folder),
        lamps=("ThAr",),
        lamp_files=[("ThAr", folder / "ThAr-0.3s-x3.sif")],
    )


# ---------------------------------------------------------------------------
# What a thin snapshot is
# ---------------------------------------------------------------------------


def test_a_thin_snapshot_holds_what_it_computed_and_points_at_the_rest(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)

    inside = {path.name for path in snapshot.root.rglob("*")}
    assert inside == {
        "snapshot.toml",
        "pattern.txt",
        "wavelength.txt",
        "integral.txt",
    }
    # Not one raw byte, and no lamps/ folder of copies to hold them in.
    assert not any(name.lower().endswith(".sif") for name in inside)

    referenced = {
        artifact.role: artifact for artifact in snapshot.artifacts if artifact.is_reference
    }
    assert set(referenced) == {"sphere", "sphere_background", "lamp"}
    assert referenced["sphere"].path == "../../sphere-0.1s-x3.sif"
    assert referenced["lamp"].path == "../../ThAr-0.3s-x3.sif"
    assert referenced["lamp"].label == "ThAr"
    assert referenced["sphere"].size_bytes == len(RAW_SPHERE)
    assert snapshot.source_path("sphere") == calibration_folder / "sphere-0.1s-x3.sif"
    assert snapshot.source_path("sphere").read_bytes() == RAW_SPHERE

    # The computed files are still the snapshot's own, copied and contained.
    copied = {artifact.role for artifact in snapshot.artifacts if not artifact.is_reference}
    assert copied == {"pattern", "wavelength", "integral"}


def test_the_manifest_says_which_entries_are_references_rather_than_hinting(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)
    text = (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")

    assert text.count(f'kind = "{ARTIFACT_REFERENCED}"') == 3
    # A copied entry is the unstated default, exactly as every binder written
    # before references existed states it.
    assert f'kind = "{ARTIFACT_COPIED}"' not in text
    # And nothing machine-specific leaked into the reference paths.
    assert str(calibration_folder) not in text


def test_a_source_outside_the_calibration_folder_is_referenced_absolutely(
    calibration_folder: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    """A relative path across unrelated trees would only pretend to be portable."""

    elsewhere = tmp_path_factory.mktemp("archive")
    lamp = elsewhere / "ThAr-archive.sif"
    lamp.write_bytes(RAW_LAMP)

    snapshot = create_snapshot(
        calibration_folder / "calibrations",
        snapshot_id="20260813_cmos",
        detector="cmos",
        files=_sources(calibration_folder),
        lamps=("ThAr",),
        lamp_files=[("ThAr", lamp)],
        reference_raw=True,
    )

    artifact = next(item for item in snapshot.artifacts if item.role == "lamp")
    assert Path(artifact.path).is_absolute()
    assert snapshot.path_for(artifact) == lamp.resolve()


# ---------------------------------------------------------------------------
# What validation says when a reference stops being true
# ---------------------------------------------------------------------------


def test_a_missing_referenced_frame_is_refused_by_its_absolute_path(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)
    missing = calibration_folder / "sphere-0.1s-x3.sif"
    missing.unlink()

    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root)

    assert caught.value.errors == (
        f"referenced sphere source not found: {missing.resolve()}",
    )


def test_a_changed_referenced_frame_is_refused_by_its_absolute_path(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)
    changed = calibration_folder / "ThAr-0.3s-x3.sif"
    changed.write_bytes(RAW_LAMP + b"extra")

    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root)

    resolved = changed.resolve()
    assert f"referenced lamp source digest mismatch: {resolved}" in caught.value.errors
    assert f"referenced lamp source size mismatch: {resolved}" in caught.value.errors


def test_an_unknown_kind_is_refused_rather_than_guessed_at(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)
    manifest = snapshot.root / "snapshot.toml"
    manifest.write_text(
        manifest.read_text(encoding="utf-8").replace(
            f'kind = "{ARTIFACT_REFERENCED}"', 'kind = "linked"', 1
        ),
        encoding="utf-8",
    )

    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root, verify_files=False)

    assert any("unknown kind 'linked'" in error for error in caught.value.errors)


# ---------------------------------------------------------------------------
# Reading a saved folder cheaply: "is my calibration done?"
# ---------------------------------------------------------------------------


def test_the_light_reading_verifies_its_own_files_and_leaves_the_frames_alone(
    calibration_folder: Path,
) -> None:
    """Owner, 2026-08-18: "I'm not sure if that's done or not."

    The answer has to arrive without hashing 380 MB of raw detector that may
    not even be mounted, so the light check digests the snapshot's own
    kilobyte-sized files and never touches a referenced source.
    """

    snapshot = _thin(
        calibration_folder,
        lamps=("Hg", "Ne", "ThAr", "Xe"),
        alignment={"rms_px": 0.53},
        qc={"sphere_comparison": "ready"},
    )
    # As if the share the frames live on were not mounted this morning.
    (calibration_folder / "sphere-0.1s-x3.sif").unlink()
    (calibration_folder / "ThAr-0.3s-x3.sif").unlink()

    reading = read_snapshot_folder(snapshot.root)

    assert reading.valid and reading.errors == ()
    assert reading.snapshot_id == "20260813_cmos"
    assert reading.lamps == ("Hg", "Ne", "ThAr", "Xe")
    assert reading.alignment_rms_px == pytest.approx(0.53)
    assert reading.sphere_comparison == "ready"
    assert reading.saved_local  # the binder's UTC stamp, on this machine's clock
    summary = reading.summary()
    assert summary.startswith("20260813_cmos — saved ")
    assert "Hg/Ne/ThAr/Xe" in summary
    assert "RMS 0.53 px" in summary
    assert "sphere comparison ready" in summary
    assert summary.endswith("validated")

    # The full check is unchanged, and still says the frames are gone.
    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root)
    assert any("sphere source not found" in error for error in caught.value.errors)


def test_the_light_reading_names_the_file_whose_digest_moved(
    calibration_folder: Path,
) -> None:
    snapshot = _thin(calibration_folder)
    (snapshot.root / "pattern.txt").write_text("edited by hand\n", encoding="utf-8")

    reading = read_snapshot_folder(snapshot.root)

    assert not reading.valid
    assert "artifact digest mismatch: pattern.txt" in reading.errors
    assert "DOES NOT VALIDATE — artifact digest mismatch: pattern.txt" in reading.summary()


def test_saved_snapshots_are_listed_newest_first_and_nothing_else_is_walked(
    calibration_folder: Path,
) -> None:
    """One level deep: the settings bundles are not snapshots, and neither is
    anything under the folders the binders point back at."""

    first = _thin(calibration_folder, snapshot_id="20260813_cmos")
    second = _thin(calibration_folder, snapshot_id="20260814_cmos")
    # Newest is the one saved most recently, which is not the same fact as the
    # one whose id carries the later acquisition date — a re-run of an old
    # campaign is exactly the case that separates them.
    _restamp(first, "2026-08-20T10:46:00+00:00")
    _restamp(second, "2026-08-14T09:00:00+00:00")
    root = calibration_folder / "calibrations"
    (root / "configs").mkdir()
    (root / "configs" / "campaign.toml").write_text("not a snapshot\n", encoding="utf-8")
    (root / "loose-note.txt").write_text("not a snapshot either\n", encoding="utf-8")

    readings = saved_snapshots_in(root)

    assert [reading.snapshot_id for reading in readings] == [
        "20260813_cmos",
        "20260814_cmos",
    ]
    assert all(reading.valid for reading in readings)
    assert saved_snapshots_in(calibration_folder / "nowhere") == ()


# ---------------------------------------------------------------------------
# The older, all-copied snapshots
# ---------------------------------------------------------------------------


def test_an_all_copied_snapshot_still_validates_exactly_as_before(
    calibration_folder: Path,
) -> None:
    snapshot = _fat(calibration_folder)

    assert {path.name for path in snapshot.root.rglob("*")} == {
        "snapshot.toml",
        "pattern.txt",
        "wavelength.txt",
        "integral.txt",
        "sphere.sif",
        "sphere_bg.sif",
        "lamps",
        "ThAr-0.3s-x3.sif",
    }
    assert "kind" not in (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")
    assert all(not artifact.is_reference for artifact in snapshot.artifacts)
    assert load_snapshot(snapshot.root).calibration_files() == {
        "orders": "pattern.txt",
        "wavelength": "wavelength.txt",
        "sphr": "sphere.sif",
        "bkgr": "sphere_bg.sif",
        "integral": "integral.txt",
    }
    # Removing the copy is still the plain "artifact file not found", by the
    # relative name that is all a contained artifact ever had.
    (snapshot.root / "sphere.sif").unlink()
    with pytest.raises(SnapshotValidationError) as caught:
        load_snapshot(snapshot.root)
    assert "artifact file not found: sphere.sif" in caught.value.errors


# ---------------------------------------------------------------------------
# What the readers make of a mixed manifest
# ---------------------------------------------------------------------------


def test_a_mixed_manifest_round_trips_through_every_reader(
    calibration_folder: Path,
) -> None:
    snapshot = load_snapshot(_thin(calibration_folder).root)

    files = snapshot.calibration_files()
    # Computed files keep their bare role names, which the loader resolves
    # inside the snapshot folder; referenced frames arrive already resolved.
    assert files["orders"] == "pattern.txt"
    assert Path(files["sphr"]).is_absolute()
    assert Path(files["sphr"]).read_bytes() == RAW_SPHERE
    assert Path(files["bkgr"]).read_bytes() == RAW_SPHERE_BACKGROUND
    assert (snapshot.root / files["wavelength"]).is_file()

    records = json.loads(snapshot.provenance_attrs()["calibration_file_digests_json"])
    by_role = {record["role"]: record for record in records}
    assert by_role["sphere"]["kind"] == ARTIFACT_REFERENCED
    assert by_role["sphere"]["path"] == "../../sphere-0.1s-x3.sif"
    assert by_role["sphere"]["source_name"] == "sphere-0.1s-x3.sif"
    # A copied entry says nothing about kind, so a cube written from an
    # all-copied snapshot carries the provenance string it always did.
    assert "kind" not in by_role["pattern"]


def test_a_relative_reference_is_read_from_the_snapshot_not_the_shell(
    calibration_folder: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    """The CWD trap: a same-named decoy underfoot must not be read instead."""

    snapshot_root = _thin(calibration_folder).root
    decoy_home = tmp_path_factory.mktemp("elsewhere")
    decoy_parent = decoy_home / "calibrations"
    decoy_parent.mkdir()
    (decoy_home / "sphere-0.1s-x3.sif").write_bytes(b"a completely different sphere")

    previous = Path.cwd()
    os.chdir(decoy_parent)
    try:
        loaded = load_snapshot(snapshot_root)
        assert loaded.source_path("sphere").read_bytes() == RAW_SPHERE
        assert loaded.source_path("sphere") == (
            calibration_folder / "sphere-0.1s-x3.sif"
        ).resolve()
    finally:
        os.chdir(previous)


def test_a_snapshot_opened_by_a_relative_path_still_finds_its_frames(
    calibration_folder: Path,
) -> None:
    snapshot_id = _thin(calibration_folder).snapshot_id
    previous = Path.cwd()
    os.chdir(calibration_folder)
    try:
        loaded = load_snapshot(Path("calibrations") / snapshot_id)
        assert loaded.source_path("sphere").read_bytes() == RAW_SPHERE
    finally:
        os.chdir(previous)


def test_a_snapshot_derived_from_a_thin_one_points_at_the_same_frames(
    calibration_folder: Path,
) -> None:
    """What a refinement does: re-solve one table, inherit every reference.

    ``create_refinement_snapshot`` builds its new folder from the base's
    resolved sources and carries the base's thinness across, so accepting a
    drift refinement never quietly duplicates the raw frames.
    """

    base = _thin(calibration_folder)
    derived = create_snapshot(
        calibration_folder / "calibrations",
        snapshot_id="20260813_cmos-r1",
        detector="cmos",
        files={
            role: base.source_path(role)
            for role in ("pattern", "wavelength", "integral", "sphere", "sphere_background")
        },
        lamps=base.lamps,
        lamp_files=[
            (artifact.label, base.path_for(artifact))
            for artifact in base.artifacts
            if artifact.role == "lamp"
        ],
        base_snapshot=base.snapshot_id,
        reference_raw=any(artifact.is_reference for artifact in base.artifacts),
    )

    assert not any(path.suffix.lower() == ".sif" for path in derived.root.rglob("*"))
    assert derived.artifact_for_role("sphere").path == "../../sphere-0.1s-x3.sif"
    assert derived.source_path("sphere") == base.source_path("sphere")


def test_show_and_validate_say_where_the_referenced_frames_are(
    calibration_folder: Path, capsys
) -> None:
    from echelle_spectra import snapshot_cli

    root = _thin(calibration_folder).root

    assert snapshot_cli.main(["validate", str(root)]) == 0
    validated = capsys.readouterr().out
    assert "3 referenced source(s) verified where they live" in validated

    assert snapshot_cli.main(["show", str(root)]) == 0
    shown = capsys.readouterr().out
    assert f"reference: sphere -> {calibration_folder / 'sphere-0.1s-x3.sif'}" in shown
    assert f"reference: lamp/ThAr -> {calibration_folder / 'ThAr-0.3s-x3.sif'}" in shown


# ---------------------------------------------------------------------------
# End to end, the way a run reaches one
# ---------------------------------------------------------------------------


def test_the_registry_resolves_a_thin_snapshot_down_to_readable_frames(
    calibration_folder: Path,
) -> None:
    _thin(calibration_folder, validity={"date_from": "2026-08-13"})
    registry_path = calibration_folder / "calibration_registry.toml"
    registry_path.write_text(
        'schema = "echelle-calibration-registry/v1"\n'
        "\n[[epochs]]\n"
        'snapshot_id = "20260813_cmos"\n',
        encoding="utf-8",
    )

    registry = load_calibration_registry(
        registry_path, snapshots_root=calibration_folder / "calibrations"
    )
    epoch = registry.resolve(
        CalibrationSourceIdentity(
            Path("201234_Echelle.SIF"), acquisition_date=date(2026, 9, 1)
        )
    )

    assert epoch.snapshot_id == "20260813_cmos"
    files = epoch.snapshot.calibration_files()
    resolved = {
        key: Path(value) if Path(value).is_absolute() else epoch.snapshot.root / value
        for key, value in files.items()
    }
    assert all(path.is_file() for path in resolved.values())
    assert resolved["sphr"].read_bytes() == RAW_SPHERE
