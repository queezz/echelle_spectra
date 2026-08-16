"""UI-independent tests for calibration campaign memory and exposure triage."""

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
    ReferenceState,
    SaveState,
    TomlState,
    catalog_lines_for_order,
    compute_absolute_calibration,
    evaluate_exposure,
    expected_lines_for_order,
    lamp_reference_set,
    measure_saturation_clusters,
    normalize_lamp_name,
    suggest_file_roles,
    triage_exposure,
    triage_for_role,
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
        suggested_lamps=("Th-Ar",),
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


#: The owner's real 2025-09-26 calibration folder: Ne bright/dim pairs and a
#: sphere pair, with the names the acquisition software actually wrote.
_REAL_2025_NAMES = (
    "Ne-0.02s-x3-bright-lines.sif",
    "Ne-0.02s-x3-bright-lines_bg.sif",
    "Ne-0.1s-x3-dimm-lines.sif",
    "Ne-0.1s-x3-dimm-lines-bg.sif",
    "sphere-0.1s-x3.sif",
    "sphere-0.1s-x3-bg.sif",
)


def _real_2025_folder(tmp_path: Path) -> dict[str, Path]:
    """Mirror the real Ne-only campaign folder, references included."""

    folder = tmp_path / "20250926_calib"
    folder.mkdir()
    sources = {}
    for name in _REAL_2025_NAMES:
        path = folder / name
        path.write_bytes((name + "\n").encode())
        sources[name] = path
    for name, text in (
        ("pattern_CMOS_20250926.txt", "pattern\n"),
        ("wavelength.txt", _CURATED_TABLE),
        ("integrating_sphere.txt", "integral\n"),
    ):
        path = folder / name
        path.write_text(text, encoding="utf-8")
        sources[name] = path
    return sources


def _ne_campaign(sources: dict[str, Path]) -> CalibrationCampaignSession:
    return CalibrationCampaignSession(
        pattern_source=sources["pattern_CMOS_20250926.txt"],
        wavelength_source=sources["wavelength.txt"],
        integral_source=sources["integrating_sphere.txt"],
    )


def _assign_real_2025_roles(campaign: CalibrationCampaignSession) -> None:
    """Assign every role by hand, exactly as the operator would click them."""

    roles = {
        "sphere-0.1s-x3.sif": (MeasurementRole.SPHERE, ""),
        "sphere-0.1s-x3-bg.sif": (MeasurementRole.SPHERE_BACKGROUND, ""),
        "Ne-0.02s-x3-bright-lines.sif": (MeasurementRole.LAMP, "Ne"),
        "Ne-0.02s-x3-bright-lines_bg.sif": (MeasurementRole.LAMP_BACKGROUND, "Ne"),
        "Ne-0.1s-x3-dimm-lines.sif": (MeasurementRole.LAMP, "Ne"),
        "Ne-0.1s-x3-dimm-lines-bg.sif": (MeasurementRole.LAMP_BACKGROUND, "Ne"),
    }
    for name, (role, lamp) in roles.items():
        source = campaign.pattern_source.parent / name
        campaign.classify_file(source, role, lamp_family=lamp)


def test_real_2025_ne_folder_reaches_a_registrable_snapshot_without_thar(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    alignment = _aligned_session(tmp_path)
    _assign_real_2025_roles(campaign)
    assert campaign.assigned_lamps == ("Ne",)

    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", "20250926_cmos", alignment)
    snapshot = campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20250926_cmos",
        detector="cmos",
        alignment=alignment,
    )

    assert campaign.save_state is SaveState.VALIDATED
    checklist = campaign.checklist(alignment)
    assert all(
        item.state is ChecklistState.DONE for item in checklist if item.blocking
    ), [item.key for item in checklist if item.blocking and item.state is not ChecklistState.DONE]
    keys = {item.key for item in checklist}
    assert "lamp-Ne-signal" in keys and "lamp-Ne-background" in keys
    assert not any("ThAr" in key for key in keys)
    assert snapshot.lamps == ("Ne",)
    registry_path = tmp_path / "calibration-registry.toml"
    registry_path.write_text(
        'schema = "echelle-calibration-registry/v1"\n'
        "\n[[epochs]]\n"
        'snapshot_id = "20250926_cmos"\n',
        encoding="utf-8",
    )
    registry = load_calibration_registry(
        registry_path, snapshots_root=tmp_path / "calibrations"
    )
    epoch = registry.resolve(
        CalibrationSourceIdentity(
            Path("20251001_lamp.sif"), acquisition_date=date(2025, 10, 1)
        )
    )
    assert epoch.snapshot_id == "20250926_cmos"
    manifest_text = (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")
    campaign_toml = (tmp_path / "configs" / "20250926_cmos" / "campaign.toml").read_text(
        encoding="utf-8"
    )
    assert "ThAr" not in manifest_text
    assert "ThAr" not in campaign_toml


def test_a_thar_pair_adds_its_own_steps_to_the_same_checklist(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    alignment = _aligned_session(tmp_path)
    _assign_real_2025_roles(campaign)
    ne_only = {item.key for item in campaign.checklist(alignment)}

    for name, role in (
        ("ThAr-0.3s-x3.sif", MeasurementRole.LAMP),
        ("ThAr-0.3s-x3-bg.sif", MeasurementRole.LAMP_BACKGROUND),
    ):
        path = sources["wavelength.txt"].parent / name
        path.write_bytes((name + "\n").encode())
        campaign.classify_file(path, role, lamp_family="ThAr")

    with_thar = {item.key for item in campaign.checklist(alignment)}

    assert campaign.assigned_lamps == ("Ne", "ThAr")
    assert with_thar - ne_only == {"lamp-ThAr-signal", "lamp-ThAr-background"}


def test_the_previous_campaigns_lamps_are_advice_and_never_block(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    alignment = _aligned_session(tmp_path)
    _assign_real_2025_roles(campaign)

    suggestion = {item.key: item for item in campaign.checklist(alignment)}[
        "lamp-suggestions"
    ]

    assert suggestion.state is ChecklistState.SUGGESTION
    assert not suggestion.blocking
    assert "last time: Ne" in suggestion.detail
    assert "consider ThAr" in suggestion.detail
    assert campaign.ready_for_snapshot("20250926_cmos", alignment) is False
    campaign.compute_sphere_comparison(_calculator)
    campaign.write_tomls(tmp_path / "configs", "20250926_cmos", alignment)
    assert campaign.ready_for_snapshot("20250926_cmos", alignment) is True


def test_every_unfinished_step_names_what_unblocks_it(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    alignment = _aligned_session(tmp_path)

    items = {item.key: item for item in campaign.checklist(alignment)}

    assert all(
        item.unblocked_by
        for item in items.values()
        if item.blocking and item.state is ChecklistState.WAITING
    )
    assert "sphere" in items["sphere-comparison"].unblocked_by
    assert "no lamp is needed" in items["sphere-comparison"].unblocked_by
    assert "any lamp name works" in items["lamp-any"].unblocked_by
    campaign.classify_file(
        sources["sphere-0.1s-x3.sif"], MeasurementRole.SPHERE
    )
    campaign.classify_file(
        sources["sphere-0.1s-x3-bg.sif"], MeasurementRole.SPHERE_BACKGROUND
    )
    after = {item.key: item for item in campaign.checklist(alignment)}
    assert "Compute factors" in after["sphere-comparison"].unblocked_by


@pytest.mark.parametrize("lamp", ["Kr", "Xe-2", "ne", "TH-AR"])
def test_any_lamp_name_is_accepted_and_normalized_only_when_known(tmp_path, lamp):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    source = sources["Ne-0.02s-x3-bright-lines.sif"]

    record = campaign.classify_file(source, MeasurementRole.LAMP, lamp_family=lamp)

    assert record.lamp_family == normalize_lamp_name(lamp)
    assert campaign.assigned_lamps == (record.lamp_family,)
    assert f"lamp-{record.lamp_family}-signal" in {
        item.key for item in campaign.checklist(_aligned_session(tmp_path))
    }


def test_a_meaninglessly_named_file_still_takes_any_role(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    nameless = sources["wavelength.txt"].parent / "IMG_0042.sif"
    nameless.write_bytes(b"nameless\n")

    suggestion = campaign.observe_file(nameless)
    record = campaign.classify_file(
        nameless, MeasurementRole.LAMP, lamp_family="Kr"
    )

    assert not suggestion.is_unambiguous
    assert "assign one by hand" in suggestion.reason
    assert record.role is MeasurementRole.LAMP
    assert record.lamp_family == "Kr"


def test_filename_help_prefills_the_real_2025_names(tmp_path):
    sources = _real_2025_folder(tmp_path)

    bright = suggest_file_roles(sources["Ne-0.02s-x3-bright-lines.sif"])
    bright_bg = suggest_file_roles(sources["Ne-0.02s-x3-bright-lines_bg.sif"])
    dim_bg = suggest_file_roles(sources["Ne-0.1s-x3-dimm-lines-bg.sif"])
    sphere_bg = suggest_file_roles(sources["sphere-0.1s-x3-bg.sif"])

    assert bright.roles == (MeasurementRole.LAMP,) and bright.lamp_name == "Ne"
    assert bright_bg.roles == (MeasurementRole.LAMP_BACKGROUND,)
    assert bright_bg.lamp_name == "Ne"
    assert dim_bg.roles == (MeasurementRole.LAMP_BACKGROUND,)
    assert sphere_bg.roles == (MeasurementRole.SPHERE_BACKGROUND,)


def test_loaded_frames_are_triaged_before_any_role_exists(tmp_path):
    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    frame = _triage_frame(tmp_path, "whatever.sif", _cosmic_singles)

    record = campaign.record_frame(frame)

    assert record.triage.state is ExposureState.GOOD
    assert record.triage.saturation.anomalous_pixels == 3
    assert campaign.loaded[frame.path] is record
    assert not campaign.measurements
    items = {item.key: item for item in campaign.checklist(_aligned_session(tmp_path))}
    assert items["files"].state is ChecklistState.DONE
    assert "1 file(s) triaged" in items["files"].detail


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
    # No lamp was assigned, so the procedure asks for a lamp without naming one.
    assert checklist["lamp-any"].state is ChecklistState.WAITING


def test_role_suggestions_are_help_not_automatic_classification(tmp_path):
    source = tmp_path / "sphere-0.1s-bg.sif"
    source.touch()
    suggestion = suggest_file_roles(source)
    assert suggestion.roles == (MeasurementRole.SPHERE_BACKGROUND,)


def test_a_suggested_role_is_reported_as_unconfirmed_not_as_assigned(tmp_path):
    """F14 item 1: the whole real folder pre-fills, and none of it is assigned.

    This is what the owner saw: every Role control read correctly while the
    campaign held nothing, so the procedure said "no file carries this role
    yet" and the factor computation failed. The campaign now names that
    in-between state instead of leaving the surface to imply an assignment.
    """

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    for name in _REAL_2025_NAMES:
        campaign.record_frame(
            _triage_frame(sources[name].parent, name, _bright_line)
        )

    pending = dict(campaign.unconfirmed_suggestions())
    assert not campaign.measurements
    assert pending[sources["sphere-0.1s-x3.sif"]] is MeasurementRole.SPHERE
    assert (
        pending[sources["sphere-0.1s-x3-bg.sif"]] is MeasurementRole.SPHERE_BACKGROUND
    )

    checklist = {item.key: item for item in campaign.checklist(_aligned_session(tmp_path))}
    assert checklist["sphere"].state is ChecklistState.WAITING
    assert "only suggested by its filename" in checklist["sphere"].detail
    assert "sphere-0.1s-x3.sif" in checklist["sphere"].detail
    assert "Confirm suggested roles" in checklist["sphere"].unblocked_by
    assert "carry no role yet" in checklist["files"].detail

    comparison = campaign.compute_sphere_comparison(lambda **_values: None)
    assert comparison.state is ComparisonState.FAILED
    assert "only shows a filename suggestion" in comparison.reason

    confirmed = campaign.confirm_suggested_roles()

    assert {record.path.name for record in confirmed} == set(_REAL_2025_NAMES)
    assert campaign.unconfirmed_suggestions() == ()
    assert campaign.measurements[sources["sphere-0.1s-x3.sif"]].role is (
        MeasurementRole.SPHERE
    )
    assert campaign.assigned_lamps == ("Ne",)
    checklist = {item.key: item for item in campaign.checklist(_aligned_session(tmp_path))}
    assert checklist["sphere"].state is ChecklistState.DONE
    assert checklist["sphere"].detail == "sphere-0.1s-x3.sif"


def test_saturation_never_fails_a_lamp_frame_but_still_fails_a_sphere(tmp_path):
    """F14 item 2, in the owner's words: of course the dim series saturates."""

    triage = triage_exposure(_triage_frame(tmp_path, "Ne-dim.sif", _saturated_cluster))
    assert triage.state is ExposureState.SATURATED

    for role in (MeasurementRole.LAMP, MeasurementRole.LAMP_BACKGROUND):
        verdict = triage_for_role(triage, role)
        assert verdict.state is ExposureState.SATURATED
        assert not verdict.blocking
        assert verdict.is_usable
        assert "fit unsaturated lines only" in verdict.headline
        assert "cluster(s)" in verdict.headline
        assert "FAILED" not in verdict.headline.upper()
        assert "lower the exposure" not in verdict.advice.casefold()

    for role in (MeasurementRole.SPHERE, MeasurementRole.SPHERE_BACKGROUND):
        verdict = triage_for_role(triage, role)
        assert verdict.blocking
        assert not verdict.is_usable
        assert "SATURATED" in verdict.headline

    # A frame with no role yet keeps the front-door verdict untouched: triage
    # is the front door and needs nothing but a file.
    unassigned = triage_for_role(triage, None)
    assert unassigned.blocking
    assert unassigned.headline == triage.headline

    healthy = triage_exposure(_triage_frame(tmp_path, "Ne.sif", _healthy))
    assert not triage_for_role(healthy, MeasurementRole.LAMP).blocking
    assert not triage_for_role(healthy, MeasurementRole.SPHERE).blocking


def test_a_classified_lamp_frame_is_never_told_to_lower_its_exposure(tmp_path):
    """The next-action line follows the same law as the verdict."""

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    folder = sources["Ne-0.1s-x3-dimm-lines.sif"].parent
    for name in ("Ne-0.1s-x3-dimm-lines.sif", "sphere-0.1s-x3.sif"):
        campaign.record_frame(_triage_frame(folder, name, _saturated_cluster))

    lamp = campaign.classify_file(
        sources["Ne-0.1s-x3-dimm-lines.sif"], MeasurementRole.LAMP, lamp_family="Ne"
    )
    sphere = campaign.classify_file(
        sources["sphere-0.1s-x3.sif"], MeasurementRole.SPHERE
    )

    assert lamp.exposure is not None and sphere.exposure is not None
    assert "expected" in lamp.exposure.next_action
    assert "unsaturated lines only" in lamp.exposure.next_action
    assert "Lower exposure" not in lamp.exposure.next_action
    assert "Lower exposure" in sphere.exposure.next_action

    checklist = {item.key: item for item in campaign.checklist(_aligned_session(tmp_path))}
    # One frame is saturated on purpose; only the sphere frame counts against us.
    assert "1 saturated (1 saturated on purpose, on lamp frames)" in (
        checklist["files"].detail
    )


def _triage_frame(tmp_path: Path, name: str, painter) -> BenchFrame:
    """Build a synthetic detector frame with a realistic background floor."""

    rng = np.random.default_rng(11)
    images = rng.normal(300.0, 12.0, size=(1, 64, 64))
    painter(images)
    return BenchFrame(
        tmp_path / name,
        images,
        images[0],
        (np.zeros(64),),
        {"ExposureTime": 0.2},
    )


def _bright_line(images: np.ndarray, level: float = 45000.0) -> None:
    images[0, 30:34, 20:24] = level


def _healthy(images: np.ndarray) -> None:
    _bright_line(images)


def _saturated_cluster(images: np.ndarray) -> None:
    _bright_line(images)
    images[0, 10:13, 40] = 65535.0


def _cosmic_singles(images: np.ndarray) -> None:
    _bright_line(images)
    for row, column in ((5, 5), (50, 7), (20, 60)):
        images[0, row, column] = 65535.0


def _too_dim(images: np.ndarray) -> None:
    _bright_line(images, level=700.0)


@pytest.mark.parametrize(
    ("name", "painter", "expected", "clusters", "anomalies"),
    [
        ("healthy.sif", _healthy, ExposureState.GOOD, 0, 0),
        ("saturated.sif", _saturated_cluster, ExposureState.SATURATED, 1, 0),
        ("cosmic.sif", _cosmic_singles, ExposureState.GOOD, 0, 3),
        ("dim.sif", _too_dim, ExposureState.DIM, 0, 0),
    ],
)
def test_triage_judges_clusters_not_lone_full_scale_pixels(
    tmp_path, name, painter, expected, clusters, anomalies
):
    triage = triage_exposure(_triage_frame(tmp_path, name, painter))

    assert triage.state is expected
    assert triage.saturation.cluster_count == clusters
    assert triage.saturation.anomalous_pixels == anomalies
    assert triage.headline
    assert triage.histogram.counts.sum() == 64 * 64
    assert triage.top_histogram.edges[-1] >= triage.full_scale
    assert 0.0 < triage.headroom_fraction <= 1.0
    if expected is ExposureState.GOOD:
        # The lone full-scale pixels are anomalies, so headroom comes from the
        # brightest real pixel and the frame stays usable.
        assert triage.saturation.clean_peak_counts == pytest.approx(45000.0)
        assert triage.is_usable
    if anomalies:
        assert "anomalies" in triage.headline
        assert "not saturation" in triage.headline


def test_a_lone_spike_is_never_mistaken_for_signal(tmp_path):
    dark = np.zeros((1, 8, 8), dtype=float)
    dark[0, 2, 2] = 65535.0
    frame = BenchFrame(tmp_path / "spike.sif", dark, dark[0], (np.zeros(8),), {})

    triage = triage_exposure(frame)

    # The frame is dark apart from one cosmic hit, so it is dim, not saturated,
    # and its headroom comes from the real pixels rather than from the spike.
    assert triage.state is ExposureState.DIM
    assert triage.saturation.anomalous_pixels == 1
    assert triage.saturation.cluster_count == 0
    assert triage.peak_counts == pytest.approx(0.0)


def test_a_frame_of_nothing_but_a_spike_is_still_judged(tmp_path):
    only = np.full((1, 1, 1), 65535.0)
    frame = BenchFrame(tmp_path / "one.sif", only, only[0], (np.zeros(1),), {})

    triage = triage_exposure(frame)

    assert triage.state is not ExposureState.NO_DATA
    assert triage.peak_counts == pytest.approx(65535.0)
    assert triage.saturation.anomalous_pixels == 1


def test_saturation_clustering_is_per_frame_four_connectivity(tmp_path):
    images = np.zeros((2, 8, 8), dtype=float)
    # The same pixel saturated in both frames is one repeated hot pixel, and
    # two diagonally touching pixels are two separate cosmic hits.
    images[0, 4, 4] = 65535.0
    images[1, 4, 4] = 65535.0
    images[0, 1, 1] = 65535.0
    images[0, 2, 2] = 65535.0

    clusters = measure_saturation_clusters(images)

    assert clusters.cluster_count == 0
    assert clusters.anomalous_pixels == 4
    assert not clusters.is_saturated

    images[1, 6, 6] = 65535.0
    images[1, 6, 7] = 65535.0
    joined = measure_saturation_clusters(images)
    assert joined.cluster_count == 1
    assert joined.cluster_pixels == 2
    assert joined.largest_cluster_pixels == 2


@pytest.mark.parametrize(
    ("painter", "expected", "phrase"),
    [
        (_saturated_cluster, ExposureState.SATURATED, "Lower exposure"),
        (_too_dim, ExposureState.DIM, "Increase exposure"),
        (_healthy, ExposureState.GOOD, "Continue"),
        (_cosmic_singles, ExposureState.GOOD, "Continue"),
    ],
)
def test_exposure_guidance_names_next_action(tmp_path, painter, expected, phrase):
    result = evaluate_exposure(_triage_frame(tmp_path, "lamp.sif", painter))

    assert result.state is expected
    assert phrase in result.next_action
    if painter is _cosmic_singles:
        assert result.anomalous_pixels == 3
        assert result.saturated_pixels == 0


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
        suggested_lamps=("ThAr",),
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


def test_the_lamp_pair_names_the_signal_and_its_own_background(tmp_path):
    """F16 item 5: line fitting needs the signal, and the pair to subtract."""

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    assert campaign.lamp_pair("Ne") == (None, None)
    assert campaign.lamp_pair("") == (None, None)

    _assign_real_2025_roles(campaign)
    signal, background = campaign.lamp_pair("Ne")

    assert signal is not None and signal.name.endswith("dimm-lines.sif")
    assert background is not None and background.name.endswith("dimm-lines-bg.sif")
    # A lamp nobody measured has no pair, and says so without raising.
    assert campaign.lamp_pair("ThAr") == (None, None)


# ----------------------------------------------------------------------
# Packet F17 item 2 — one line source for sticks, rows, and counts
# ----------------------------------------------------------------------


def _packaged_wavelength_rows():
    resources = Path(__file__).parents[1] / "src" / "echelle_spectra" / "resources"
    return load_wavelength_table(
        resources
        / "calibration_files"
        / "alignments"
        / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    )


def test_expected_lines_are_the_lamps_own_rows_not_a_second_catalog():
    """The owner's own two counter-examples, pinned against the real table.

    Live on his 2025 Ne data the bench labelled NeI 640.225 in order 7 while
    the expected-lines table did not carry it, and labelled three NeI lines in
    order 6 under a panel reading "0 expected Ne lines in this order". The
    cause was two sources: sticks from the curated wavelength table, rows from
    the packaged NIST Ne cache interpolated onto the order. That cache stops at
    638.3 nm and starts at 580.4 nm, so it can neither see 640.225 nor reach
    order 6 at all.
    """

    rows = _packaged_wavelength_rows()
    reference = lamp_reference_set("Ne", rows)
    assert reference.state is ReferenceState.MATCHED

    order_seven = expected_lines_for_order(reference, 7)
    assert [f"{line.wavelength_nm:.3f}" for line in order_seven] != []
    assert any(
        line.species == "NeI" and abs(line.wavelength_nm - 640.225) < 0.001
        for line in order_seven
    ), [line.label for line in order_seven]

    # Order 6 carries neon and now says so; the packaged cache holds none of it.
    order_six = expected_lines_for_order(reference, 6)
    assert len(order_six) == 3
    assert {line.species for line in order_six} == {"NeI"}
    assert catalog_lines_for_order(rows, 6, "Ne") == ()

    # Every listed line is an anchorable curated row, in detector order.
    for line in (*order_six, *order_seven):
        assert line.row in reference.lines
        assert line.detector_pixel == line.row.center_pixel
        assert line.label == f"{line.species} {line.wavelength_nm:.3f}"
    pixels = [line.detector_pixel for line in order_seven]
    assert pixels == sorted(pixels)

    # No other element ever leaks in, whichever order is asked for.
    for order_idx in range(12):
        listed = expected_lines_for_order(reference, order_idx)
        assert {line.species for line in listed} <= {"NeI", "NeII"}


def test_expected_lines_annotate_from_the_catalog_without_obeying_it():
    """The packaged cache may enrich a row; it may never add or remove one."""

    rows = _packaged_wavelength_rows()
    reference = lamp_reference_set("Ne", rows)
    listed = expected_lines_for_order(reference, 9)
    assert listed
    # Order 9 sits inside the packaged cache's span, so intensities appear...
    assert any(line.relative_intensity is not None for line in listed)
    assert all(
        line.catalog is None or abs(line.catalog.wavelength_nm - line.wavelength_nm) <= 0.05
        for line in listed
    )
    # ...and the count is still exactly the curated rows of that order.
    assert len(listed) == len(
        [row for row in reference.lines if row.order_idx == 9]
    )


def test_an_unreferenceable_lamp_lists_nothing_and_says_why():
    """A lamp with no catalog contributes no expected lines, honestly."""

    rows = _packaged_wavelength_rows()
    unknown = lamp_reference_set("Kr", rows)
    assert unknown.state is ReferenceState.NO_CATALOG
    assert expected_lines_for_order(unknown, 9) == ()
    assert "no line catalog for Kr" in unknown.message

    # With no lamp assigned at all the whole table is the fallback, so the
    # sticks the bench already draws are exactly what the panel lists.
    fallback = expected_lines_for_order(None, 9, fallback_lines=rows)
    assert len(fallback) == len([row for row in rows if row.order_idx == 9])

    # The strongest-order hint reads the same list it is a hint about.
    neon = lamp_reference_set("Ne", rows)
    assert neon.best_order == max(
        {row.order_idx for row in neon.lines},
        key=lambda order: len([r for r in neon.lines if r.order_idx == order]),
    )


def test_the_snapshot_records_which_vetted_line_set_anchored_it(tmp_path):
    """F19 rider: RMS says how self-consistent a fit was, never whose it is.

    A later reader deciding whether to believe a snapshot needs to know which
    lines were trusted and on whose authority, so the manifest carries the
    vetted set by name — and carries its absence by name too, rather than
    leaving a blank that reads like an oversight.
    """

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    _campaign_obj, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    manifest = tomllib.loads(
        (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")
    )

    # The fixture table is nobody's vetted set, and says so plainly.
    assert manifest["alignment"]["vetted_set"] == ""
    assert manifest["alignment"]["vetted_lineage"] == ["wavelength.txt"]


def test_a_snapshot_anchored_on_the_bh_paper_set_names_it(tmp_path):
    """The same manifest field, when the lineage really does reach the paper.

    The lineage is a table's own header rather than its filename, so the
    fixture declares the base it was derived from exactly as a real aligned
    table does, and inherits the vetting the same way.
    """

    sources = _curated_sources(tmp_path)
    sources["wavelength.txt"].write_text(
        "# Adjusted wavelength calibration lookup table\n"
        "# Base wavelength file: Th_wavelength_CMOS_20240305.txt\n"
        + _CURATED_TABLE,
        encoding="utf-8",
    )
    alignment = _aligned_session(tmp_path)
    _campaign_obj, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    manifest = tomllib.loads(
        (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")
    )

    assert manifest["alignment"]["vetted_set"] == "BH paper"
    assert manifest["alignment"]["vetted_set_source"] == (
        "Th_wavelength_CMOS_20240305.txt"
    )
    assert manifest["alignment"]["vetted_lineage"] == [
        "wavelength.txt",
        "Th_wavelength_CMOS_20240305.txt",
    ]
