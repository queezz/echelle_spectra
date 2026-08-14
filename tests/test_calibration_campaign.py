"""UI-independent tests for Packet 5 calibration campaign memory."""

from __future__ import annotations

import hashlib
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import pytest

from echelle_spectra.calibration_bench import (
    BenchFrame,
    CalibrationBenchSession,
    FrameLoader,
    StableSifWatcher,
)
from echelle_spectra.calibration_campaign import (
    AbsoluteCalibrationResult,
    CalibrationCampaignSession,
    ChecklistState,
    ComparisonState,
    ExposureState,
    MeasurementRole,
    SaveState,
    TomlState,
    catalog_lines_for_order,
    compute_absolute_calibration,
    evaluate_exposure,
    suggest_file_roles,
)
from echelle_spectra.calibration_registry import (
    CalibrationSourceIdentity,
    load_calibration_registry,
)
from echelle_spectra.snapshot import SnapshotError, load_snapshot
from echelle_spectra.tools.calibration_alignment import (
    CalibrationTableLine,
    load_alignment_settings,
    load_wavelength_table,
)

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10
    import tomli as tomllib


def _sources(tmp_path: Path) -> dict[str, Path]:
    answer = {}
    for name in (
        "pattern.txt",
        "wavelength.txt",
        "integral.txt",
        "sphere.sif",
        "sphere_bg.sif",
        "previous_sphere.sif",
        "previous_sphere_bg.sif",
        "thar.sif",
        "thar_bg.sif",
    ):
        path = tmp_path / name
        path.write_bytes((name + "\n").encode())
        answer[name] = path
    return answer


def _aligned_session(tmp_path: Path, shift_px: float = 1.0) -> CalibrationBenchSession:
    columns = 100
    pattern = np.column_stack(
        [np.full(columns, 12, dtype=float), np.full(columns, 30, dtype=float)]
    )
    lines = (
        CalibrationTableLine(0, 20, 30, 25, 600.0, "ThI", "ok"),
        CalibrationTableLine(1, 50, 60, 55, 610.0, "ArI", "ok"),
    )
    x = np.arange(columns, dtype=float)
    spectra = (
        20 + 800 * np.exp(-0.5 * ((x - 25 - shift_px) / 1.5) ** 2),
        20 + 700 * np.exp(-0.5 * ((x - 55 - shift_px) / 1.5) ** 2),
    )
    images = np.zeros((1, 44, columns), dtype=float)
    images[0, 12, :] = spectra[0]
    images[0, 30, :] = spectra[1]
    session = CalibrationBenchSession(pattern, lines, minimum_snr=3.0)
    session.accept_frame(
        BenchFrame(
            tmp_path / "thar.sif",
            images,
            images[0],
            spectra,
            {"ExposureTime": 0.1},
        )
    )
    assert session.fit_anchor_at(0, 25 + shift_px).accepted
    assert session.fit_anchor_at(1, 55 + shift_px).accepted
    return session


_CURATED_TABLE = """\
# Curated fixture wavelength lookup table
# order from to center wavelength band
0\t0020.000\t0030.000\t00025.0000\t00600.00000\tThI  # ok
0\t0065.000\t0075.000\t00070.0000\t00604.50000\tArI  # ok
1\t0050.000\t0060.000\t00055.0000\t00610.00000\tArI  # ok
1\t0080.000\t0090.000\t00085.0000\t00614.25000\tThI  # ?
"""


def _curated_sources(tmp_path: Path) -> dict[str, Path]:
    """Give the shared fixture a wavelength table that really has rows."""

    sources = _sources(tmp_path)
    sources["wavelength.txt"].write_text(_CURATED_TABLE, encoding="utf-8")
    return sources


def _saved_snapshot(
    tmp_path: Path,
    sources: dict[str, Path],
    alignment: CalibrationBenchSession,
    *,
    snapshot_id: str = "20250813_cmos",
    validity=None,
):
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", snapshot_id, alignment)
    snapshot = campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        alignment=alignment,
        validity=validity,
    )
    return campaign, snapshot


def _campaign(tmp_path: Path, sources: dict[str, Path]) -> CalibrationCampaignSession:
    return CalibrationCampaignSession(
        pattern_source=sources["pattern.txt"],
        wavelength_source=sources["wavelength.txt"],
        integral_source=sources["integral.txt"],
        required_lamps=("Th-Ar",),
        previous_sphere=sources["previous_sphere.sif"],
        previous_sphere_background=sources["previous_sphere_bg.sif"],
    )


def _calculator(**values) -> AbsoluteCalibrationResult:
    scale = 1.1 if "previous" not in values["sphere"].name else 1.0
    wavelength = np.linspace(400.0, 700.0, 100)
    return AbsoluteCalibrationResult(wavelength, np.full(100, scale))


def _classify_complete(
    campaign: CalibrationCampaignSession,
    sources: dict[str, Path],
    frame: BenchFrame,
) -> None:
    campaign.classify_file(
        sources["sphere.sif"], MeasurementRole.SPHERE, frame=frame
    )
    campaign.classify_file(
        sources["sphere_bg.sif"], MeasurementRole.SPHERE_BACKGROUND, frame=frame
    )
    campaign.classify_file(
        sources["thar.sif"], MeasurementRole.LAMP, lamp_family="ThAr", frame=frame
    )
    campaign.classify_file(
        sources["thar_bg.sif"],
        MeasurementRole.LAMP_BACKGROUND,
        lamp_family="ThAr",
        frame=frame,
    )


def test_filename_help_never_ticks_an_ambiguous_role(tmp_path):
    sources = _sources(tmp_path)
    campaign = _campaign(tmp_path, sources)
    alignment = _aligned_session(tmp_path)

    suggestion = campaign.observe_file(tmp_path / "background.sif")

    assert not suggestion.is_unambiguous
    assert set(suggestion.roles) == {
        MeasurementRole.SPHERE_BACKGROUND,
        MeasurementRole.LAMP_BACKGROUND,
    }
    checklist = {item.key: item for item in campaign.checklist(alignment)}
    assert checklist["sphere-background"].state is ChecklistState.WAITING
    assert checklist["lamp-ThAr-background"].state is ChecklistState.WAITING


def test_role_suggestions_are_help_not_automatic_classification(tmp_path):
    source = tmp_path / "sphere-0.1s-bg.sif"
    source.touch()
    suggestion = suggest_file_roles(source)
    assert suggestion.roles == (MeasurementRole.SPHERE_BACKGROUND,)


@pytest.mark.parametrize(
    ("peak", "expected"),
    [
        (65535.0, ExposureState.SATURATED),
        (1000.0, ExposureState.DIM),
        (30000.0, ExposureState.GOOD),
    ],
)
def test_exposure_guidance_names_next_action(tmp_path, peak, expected):
    images = np.zeros((1, 8, 8), dtype=float)
    images[0, 3, 3] = peak
    frame = BenchFrame(
        tmp_path / "lamp.sif",
        images,
        images[0],
        (np.zeros(8),),
        {"ExposureTime": 0.2},
    )
    result = evaluate_exposure(frame)
    assert result.state is expected
    assert result.next_action
    if expected is ExposureState.SATURATED:
        assert "Lower exposure" in result.next_action
    elif expected is ExposureState.DIM:
        assert "Increase exposure" in result.next_action
    else:
        assert "Continue" in result.next_action


def test_shared_catalog_line_help_maps_packaged_thar_rows():
    rows = (
        CalibrationTableLine(8, 0, 10, 100, 578.0, "ThI", "ok"),
        CalibrationTableLine(8, 20, 30, 2000, 640.0, "ArI", "ok"),
    )
    help_lines = catalog_lines_for_order(rows, 8, "ThAr")
    assert help_lines
    assert {item.line.family for item in help_lines} == {"thar"}
    assert all(100 <= item.detector_pixel <= 2000 for item in help_lines)
    assert all(item.line.source_reference for item in help_lines)


def test_checklist_self_ticks_from_explicit_measurements_and_results(tmp_path):
    sources = _sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)

    before = {item.key: item for item in campaign.checklist(alignment)}
    assert before["sphere"].state is ChecklistState.DONE
    assert before["lamp-ThAr-signal"].state is ChecklistState.DONE
    assert before["sphere-comparison"].state is ChecklistState.WAITING
    assert before["alignment"].state is ChecklistState.DONE

    comparison = campaign.compute_sphere_comparison(_calculator)
    assert comparison.state is ComparisonState.READY
    assert comparison.median_ratio == pytest.approx(1.1)
    after = {item.key: item for item in campaign.checklist(alignment)}
    assert after["sphere-comparison"].state is ChecklistState.DONE


def test_missing_previous_pair_is_own_insufficient_data_state(tmp_path):
    sources = _sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = CalibrationCampaignSession(
        pattern_source=sources["pattern.txt"],
        wavelength_source=sources["wavelength.txt"],
        integral_source=sources["integral.txt"],
    )
    _classify_complete(campaign, sources, alignment.frame)

    comparison = campaign.compute_sphere_comparison(_calculator)

    assert comparison.state is ComparisonState.INSUFFICIENT_DATA
    assert comparison.candidate is not None
    assert "unavailable" in comparison.reason


def test_tomls_are_commented_parseable_and_machine_path_free(tmp_path):
    sources = _sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)

    paths = campaign.write_tomls(
        tmp_path / "configs", "20250813_cmos", alignment
    )

    assert campaign.toml_state is TomlState.GENERATED
    assert set(paths) == {"campaign", "alignment", "export"}
    for path in paths.values():
        text = path.read_text(encoding="utf-8")
        assert text.startswith("#")
        assert str(tmp_path) not in text
        with path.open("rb") as stream:
            assert tomllib.load(stream)


def test_complete_rehearsal_saves_and_validates_through_snapshot_api(tmp_path):
    sources = _sources(tmp_path)
    original_digests = {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in sources.items()
    }
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", "20250813_cmos", alignment)

    snapshot = campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20250813_cmos",
        detector="CMOS",
        alignment=alignment,
        notes="2025 fixture rehearsal",
        base_snapshot="20240305_cmos",
    )

    assert campaign.save_state is SaveState.VALIDATED
    assert load_snapshot(snapshot.root).snapshot_id == "20250813_cmos"
    lamp_artifacts = [item for item in snapshot.artifacts if item.role == "lamp"]
    assert len(lamp_artifacts) == 2
    assert {item.label for item in lamp_artifacts} == {"ThAr"}
    assert snapshot.manifest["qc"]["sphere_comparison"] == "ready"
    assert original_digests == {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in sources.items()
    }


def test_existing_snapshot_refusal_is_recoverable_with_new_identity(tmp_path):
    sources = _sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", "20250813_cmos", alignment)
    campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20250813_cmos",
        detector="cmos",
        alignment=alignment,
    )

    second = _campaign(tmp_path, sources)
    _classify_complete(second, sources, alignment.frame)
    second.compute_sphere_comparison(_calculator)
    second.write_tomls(tmp_path / "retry-configs", "20250813_cmos", alignment)
    with pytest.raises(SnapshotError, match="will not be replaced"):
        second.save_snapshot(
            tmp_path / "calibrations",
            snapshot_id="20250813_cmos",
            detector="cmos",
            alignment=alignment,
        )
    assert second.save_state is SaveState.FAILED

    second.write_tomls(tmp_path / "configs", "20250813_cmos-r1", alignment)
    recovered = second.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20250813_cmos-r1",
        detector="cmos",
        alignment=alignment,
    )
    assert recovered.snapshot_id == "20250813_cmos-r1"
    assert second.save_state is SaveState.VALIDATED


def test_snapshot_identity_must_match_generated_tomls(tmp_path):
    sources = _sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", "20250813_cmos", alignment)

    assert not campaign.ready_for_snapshot("20250813_cmos-r1", alignment)
    with pytest.raises(SnapshotError, match="generated TOMLs target"):
        campaign.save_snapshot(
            tmp_path / "calibrations",
            snapshot_id="20250813_cmos-r1",
            detector="cmos",
            alignment=alignment,
        )


def test_saved_wavelength_table_is_the_measured_one_not_the_base(tmp_path):
    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path, shift_px=1.0)
    base = sources["wavelength.txt"]
    base_digest = hashlib.sha256(base.read_bytes()).hexdigest()

    campaign, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    saved = snapshot.root / "wavelength.txt"
    assert hashlib.sha256(saved.read_bytes()).hexdigest() != base_digest
    dx_px = alignment.transform.dx_px
    assert dx_px == pytest.approx(1.0, abs=0.02)
    base_rows = load_wavelength_table(base)
    saved_rows = load_wavelength_table(saved)
    assert len(saved_rows) == len(base_rows) == 4
    for original, moved in zip(base_rows, saved_rows):
        assert moved.center_pixel == pytest.approx(original.center_pixel + dx_px, abs=5e-4)
        assert moved.pixel_from == pytest.approx(original.pixel_from + dx_px, abs=5e-4)
        assert moved.wavelength_nm == original.wavelength_nm
        assert moved.species == original.species
    header = saved.read_text(encoding="utf-8")
    assert "Base wavelength file: wavelength.txt" in header
    assert "Base pattern file: pattern.txt" in header
    assert "Alignment dataset: 20250813_cmos" in header
    assert "rigid detector transform" in header
    assert "Fitted lines: 2" in header
    assert campaign.wavelength_correction.applied
    assert snapshot.manifest["alignment"]["wavelength_correction_applied"] is True
    assert hashlib.sha256(base.read_bytes()).hexdigest() == base_digest


def test_transform_that_moves_nothing_copies_the_base_table_byte_for_byte(tmp_path):
    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path, shift_px=0.0)
    base = sources["wavelength.txt"]

    campaign, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    saved = snapshot.root / "wavelength.txt"
    assert saved.read_bytes() == base.read_bytes()
    assert not campaign.wavelength_correction.applied
    assert "no line" in campaign.wavelength_correction.reason
    assert snapshot.manifest["alignment"]["wavelength_correction_applied"] is False
    detail = {item.key: item for item in campaign.checklist(alignment)}["snapshot"].detail
    assert "copied unchanged" in detail


def test_bench_snapshot_is_registrable_without_hand_editing(tmp_path):
    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)

    _campaign_session, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    assert snapshot.manifest["validity"] == {"date_from": "2025-08-13"}
    assert load_snapshot(snapshot.root).snapshot_id == "20250813_cmos"
    registry_path = tmp_path / "calibration-registry.toml"
    registry_path.write_text(
        'schema = "echelle-calibration-registry/v1"\n'
        "\n[[epochs]]\n"
        'snapshot_id = "20250813_cmos"\n',
        encoding="utf-8",
    )
    registry = load_calibration_registry(
        registry_path, snapshots_root=tmp_path / "calibrations"
    )
    epoch = registry.resolve(
        CalibrationSourceIdentity(Path("20250901_lamp.sif"), acquisition_date=date(2025, 9, 1))
    )
    assert epoch.snapshot_id == "20250813_cmos"
    assert epoch.date_from == date(2025, 8, 13)
    assert epoch.date_to is None


def test_snapshot_folder_carries_round_tripping_alignment_settings(tmp_path):
    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)

    _campaign_session, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    settings_path = snapshot.root / "alignment.toml"
    assert settings_path.is_file()
    settings = load_alignment_settings(settings_path)
    assert settings.instrument_id == "cmos"
    assert settings.base_wavelength_file == "wavelength.txt"
    assert settings.base_pattern_file == "pattern.txt"
    assert settings.alignment_dataset_id == "20250813_cmos"
    assert settings.alignment_lamp == "ThAr"
    assert settings.signal_file == "thar.sif"
    assert settings.background_file == "thar_bg.sif"
    assert settings.output_wavelength_file == "wavelength.txt"
    assert settings.n_lines == 2
    assert settings.transform.dx_px == pytest.approx(alignment.transform.dx_px)
    assert settings.transform.dy_px == pytest.approx(alignment.transform.dy_px)
    assert settings.transform.theta_rad == pytest.approx(alignment.transform.theta_rad)
    assert settings.rms_px == pytest.approx(alignment.rms_px)


def test_explicit_epoch_start_overrides_the_snapshot_identity_date(tmp_path):
    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)

    _campaign_session, snapshot = _saved_snapshot(
        tmp_path, sources, alignment, validity={"date_from": "2025-08-01"}
    )

    assert snapshot.manifest["validity"] == {"date_from": "2025-08-01"}
    assert "2025-08-01" in (snapshot.root / "wavelength.txt").read_text(encoding="utf-8")


def test_packaged_2025_watch_to_validated_snapshot_rehearsal(tmp_path):
    resources = Path(__file__).parents[1] / "src" / "echelle_spectra" / "resources"
    calibration_dir = resources / "calibration_files"
    pattern_path = calibration_dir / "pattern_CMOS_20250926.txt"
    wavelength_path = (
        calibration_dir
        / "alignments"
        / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    )
    sphere_path = calibration_dir / "sphere_cmos_20240305.sif"
    sphere_background_path = calibration_dir / "sphere_cmos_20240305_bkg.sif"
    integral_path = calibration_dir / "integrating_sphere.txt"
    lamp_source = calibration_dir / "ThAr-0.3s-x3_20240305.sif"
    source_digest = hashlib.sha256(lamp_source.read_bytes()).hexdigest()

    watch_dir = tmp_path / "watch"
    watch_dir.mkdir()
    watched_lamp = watch_dir / "ThAr_2025_fixture.sif"
    shutil.copy2(lamp_source, watched_lamp)
    watched_lamp.touch()
    watcher = StableSifWatcher(watch_dir, required_unchanged_polls=2)
    modified_ns = watched_lamp.stat().st_mtime_ns
    watcher.poll(now_ns=modified_ns + 2_000_000_000)
    ready = watcher.poll(now_ns=modified_ns + 2_000_000_000)
    assert ready.ready_path == watched_lamp
    synthetic_background = tmp_path / "ThAr_bg_synthetic_fixture.sif"
    synthetic_background.write_bytes(b"synthetic fixture: no packaged 2025 lamp background")

    pattern = np.loadtxt(pattern_path, dtype=int)
    lines = load_wavelength_table(wavelength_path)
    alignment = CalibrationBenchSession(pattern, lines, minimum_snr=3.0)
    assert alignment.load_file(ready.ready_path, FrameLoader(pattern))
    first = lines[0]
    second = next(line for line in lines if line.order_idx != first.order_idx)
    assert alignment.fit_anchor_at(first.order_idx, first.center_pixel).accepted
    assert alignment.fit_anchor_at(second.order_idx, second.center_pixel).accepted

    campaign = CalibrationCampaignSession(
        pattern_source=pattern_path,
        wavelength_source=wavelength_path,
        integral_source=integral_path,
        required_lamps=("ThAr",),
        previous_sphere=sphere_path,
        previous_sphere_background=sphere_background_path,
    )
    campaign.classify_file(sphere_path, MeasurementRole.SPHERE)
    campaign.classify_file(sphere_background_path, MeasurementRole.SPHERE_BACKGROUND)
    campaign.classify_file(watched_lamp, MeasurementRole.LAMP, lamp_family="ThAr")
    campaign.classify_file(
        synthetic_background,
        MeasurementRole.LAMP_BACKGROUND,
        lamp_family="ThAr",
    )
    comparison = campaign.compute_sphere_comparison(compute_absolute_calibration)
    assert comparison.state is ComparisonState.READY
    assert comparison.median_ratio == pytest.approx(1.0)
    campaign.write_tomls(
        tmp_path / "configs", "20250926_cmos-fixture", alignment
    )
    snapshot = campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20250926_cmos-fixture",
        detector="cmos",
        alignment=alignment,
        notes=(
            "Software rehearsal: accepted 2025 geometry/wavelength resources, "
            "historical packaged ThAr/sphere SIFs, synthetic lamp background placeholder."
        ),
        base_snapshot="20240305_cmos",
    )

    assert load_snapshot(snapshot.root).snapshot_id == "20250926_cmos-fixture"
    assert campaign.save_state is SaveState.VALIDATED
    assert hashlib.sha256(lamp_source.read_bytes()).hexdigest() == source_digest
