"""UI-independent tests for calibration campaign memory and exposure triage."""

from __future__ import annotations

import hashlib
import os
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import pytest

from echelle_spectra import calibration_campaign as campaign_module
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
from echelle_spectra.tools.line_catalog import load_line_table
from echelle_spectra.tools.nist_lamp_calibration import (
    lamp_species,
    normalize_species_key,
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
    # The advice names the control the bench actually carries: the Compute
    # factors button is gone, and the strip's verb for this step is this one.
    assert "Measure sensitivity" in after["sphere-comparison"].unblocked_by


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


def test_a_folder_whose_names_leave_nothing_to_ask_assigns_itself(tmp_path):
    """F21 item 1: confirming what nobody doubts is a dance, not a decision.

    The owner's folders always read the same way ("the suggestions are always
    like so"), and every role in them had to be confirmed one row at a time.
    A set of filenames that says exactly one thing is not a question.
    """

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    for name in _REAL_2025_NAMES:
        campaign.record_frame(_triage_frame(sources[name].parent, name, _bright_line))

    pending = dict(campaign.unanimous_suggestions())
    assert set(pending) == {sources[name] for name in _REAL_2025_NAMES}
    # The bright/dim pair is not a clash: a lamp family is shot twice on
    # purpose, and only the roles one file can hold are exclusive.
    assert pending[sources["Ne-0.02s-x3-bright-lines.sif"]] is MeasurementRole.LAMP
    assert pending[sources["Ne-0.1s-x3-dimm-lines.sif"]] is MeasurementRole.LAMP

    applied = campaign.apply_unanimous_suggestions()

    assert {record.path.name for record in applied} == set(_REAL_2025_NAMES)
    assert campaign.unconfirmed_suggestions() == ()
    assert campaign.unanimous_suggestions() == ()
    assert campaign.assigned_lamps == ("Ne",)
    assert campaign.measurements[sources["sphere-0.1s-x3.sif"]].role is (
        MeasurementRole.SPHERE
    )
    checklist = {item.key: item for item in campaign.checklist(_aligned_session(tmp_path))}
    assert checklist["sphere"].state is ChecklistState.DONE


def test_two_files_claiming_one_role_are_never_applied_unasked(tmp_path):
    """A role only one file can hold, proposed by two, is the operator's call."""

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    second_sphere = sources["sphere-0.1s-x3.sif"].parent / "sphere-0.2s-x3.sif"
    second_sphere.write_bytes(b"sphere again\n")
    for name in (*_REAL_2025_NAMES, second_sphere.name):
        campaign.record_frame(
            _triage_frame(sources["sphere-0.1s-x3.sif"].parent, name, _bright_line)
        )

    assert campaign.unanimous_suggestions() == ()
    assert campaign.apply_unanimous_suggestions() == ()
    # Nothing was applied at all — not even the files nobody was in doubt about.
    assert not campaign.measurements
    # And the confirm flow is exactly where it was, for every one of them.
    assert len(campaign.unconfirmed_suggestions()) == len(_REAL_2025_NAMES) + 1


def test_one_unreadable_name_leaves_the_whole_drop_to_the_operator(tmp_path):
    """Half a folder guessed is worse than a folder asked about."""

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    nameless = sources["sphere-0.1s-x3.sif"].parent / "IMG_0042.sif"
    nameless.write_bytes(b"nameless\n")
    for name in (*_REAL_2025_NAMES, nameless.name):
        campaign.record_frame(
            _triage_frame(sources["sphere-0.1s-x3.sif"].parent, name, _bright_line)
        )

    assert campaign.unanimous_suggestions() == ()
    assert not campaign.apply_unanimous_suggestions()

    # A lamp with no readable lamp name is the same kind of doubt: which lamp?
    campaign = _ne_campaign(sources)
    lamp = sources["sphere-0.1s-x3.sif"].parent / "lamp-0.1s.sif"
    lamp.write_bytes(b"lamp\n")
    campaign.record_frame(_triage_frame(lamp.parent, lamp.name, _bright_line))
    assert campaign.observed[lamp].is_unambiguous
    assert campaign.unanimous_suggestions() == ()


def test_a_file_the_operator_unassigned_is_never_reassigned_by_a_later_drop(tmp_path):
    """An applied role is an ordinary assignment, and so is taking one off."""

    sources = _real_2025_folder(tmp_path)
    campaign = _ne_campaign(sources)
    for name in _REAL_2025_NAMES:
        campaign.record_frame(_triage_frame(sources[name].parent, name, _bright_line))
    campaign.apply_unanimous_suggestions()
    sphere = sources["sphere-0.1s-x3.sif"]
    assert campaign.remove_classification(sphere)
    later = sphere.parent / "Ne-0.3s-x3-lines.sif"
    later.write_bytes(b"one more lamp exposure\n")
    campaign.record_frame(_triage_frame(later.parent, later.name, _bright_line))

    applied = campaign.apply_unanimous_suggestions(declined=[sphere])

    # The new file lands; the one taken off by hand stays off.
    assert [record.path.name for record in applied] == [later.name]
    assert sphere not in campaign.measurements
    assert campaign.unanimous_suggestions(declined=[sphere]) == ()


def test_saturation_never_fails_a_lamp_frame_but_still_fails_a_sphere(tmp_path):
    """F14 item 2, in the owner's words: of course the dim series saturates."""

    triage = triage_exposure(_triage_frame(tmp_path, "Ne-dim.sif", _saturated_cluster))
    assert triage.state is ExposureState.SATURATED

    # The lamp SIGNAL is the frame this law is about: it is shot to saturate
    # its strong lines, so clustered saturation there is information.
    verdict = triage_for_role(triage, MeasurementRole.LAMP)
    assert verdict.state is ExposureState.SATURATED
    assert not verdict.blocking
    assert verdict.is_usable
    assert "fit unsaturated lines only" in verdict.headline
    assert "cluster(s)" in verdict.headline
    assert "FAILED" not in verdict.headline.upper()
    assert "lower the exposure" not in verdict.advice.casefold()

    # A saturated BACKGROUND is a different animal and always was: a frame
    # shot with the light off cannot legitimately reach full scale, so it is a
    # fault rather than expected saturation, whichever pair it belongs to.
    for role in (
        MeasurementRole.SPHERE,
        MeasurementRole.SPHERE_BACKGROUND,
        MeasurementRole.LAMP_BACKGROUND,
    ):
        verdict = triage_for_role(triage, role)
        assert verdict.blocking
        assert not verdict.is_usable
        assert "saturated" in verdict.headline.casefold()

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


def test_the_bench_saves_a_snapshot_that_points_at_the_frames_it_measured(tmp_path):
    """The calibration folder holds the light; the snapshot holds the arithmetic.

    A bench save used to leave a second copy of every sphere and lamp SIF one
    folder away from the originals.  It records where they are and what they
    hash to instead, so the folder that already held them completely is still
    the only place they live.
    """

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)

    _campaign_obj, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    inside = {path.name for path in snapshot.root.rglob("*")}
    assert not any(name.lower().endswith(".sif") for name in inside)
    assert "lamps" not in inside
    assert inside == {
        "snapshot.toml",
        "pattern.txt",
        "wavelength.txt",
        "integral.txt",
        "alignment.toml",
    }

    referenced = {
        artifact.path: artifact
        for artifact in snapshot.artifacts
        if artifact.is_reference
    }
    assert set(referenced) == {
        "../../sphere.sif",
        "../../sphere_bg.sif",
        "../../thar.sif",
        "../../thar_bg.sif",
    }
    for path, artifact in referenced.items():
        source = sources[Path(path).name]
        assert snapshot.path_for(artifact) == source.resolve()
        assert artifact.sha256 == hashlib.sha256(source.read_bytes()).hexdigest()
        assert artifact.size_bytes == source.stat().st_size

    # And the whole folder still passes the same validator, unchanged.
    assert load_snapshot(snapshot.root).snapshot_id == "20250813_cmos"


def test_the_generated_export_config_names_the_sphere_where_it_really_is(tmp_path):
    """The config bundle points back out too, or it would name a copy nobody made."""

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)

    paths = campaign.write_tomls(
        tmp_path / "calibrations" / "configs",
        "20250813_cmos",
        alignment,
        snapshot_root=tmp_path / "calibrations" / "20250813_cmos",
    )

    with paths["export"].open("rb") as stream:
        calibration = tomllib.load(stream)["calibration"]
    assert calibration["sphere"] == "../../sphere.sif"
    assert calibration["sphere_background"] == "../../sphere_bg.sif"
    # The computed files are the snapshot's own and are still named plainly.
    assert calibration["order_pattern"] == "pattern.txt"


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
    the packaged NIST Ne cache interpolated onto the order.

    The cache has since been widened from the 580.4--638.3 nm Fulcher window to
    the whole 380--810 nm instrument range, so it now reaches both. That closes
    the gap that made the two sources visibly disagree and leaves the real
    point standing on its own: the panel lists the lamp's curated rows, and the
    catalog having plenty to say about an order changes nothing about which
    rows are listed.
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

    # Order 6 carries neon and says so. The widened cache now has plenty of
    # neon to offer this order — and the panel still lists three rows, because
    # the curated table is what it lists and the catalog only annotates.
    order_six = expected_lines_for_order(reference, 6)
    assert len(order_six) == 3
    assert {line.species for line in order_six} == {"NeI"}
    assert len(catalog_lines_for_order(rows, 6, "Ne")) > len(order_six)

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


def test_the_snapshot_records_whether_science_lines_ever_agreed(tmp_path):
    """F19 second rider: RMS and validation are different questions.

    ``rms_px`` is the anchors agreeing with each other in pixels.  Whether the
    solution agrees with the lines physics knows is the question the BH paper
    was actually held to, and a manifest that carries only the first would let
    an unvalidated snapshot read as a validated one.
    """

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    _campaign_obj, snapshot = _saved_snapshot(tmp_path, sources, alignment)

    manifest = tomllib.loads(
        (snapshot.root / "snapshot.toml").read_text(encoding="utf-8")
    )
    alignment_block = manifest["alignment"]

    # The fixture is a thorium/argon table on a lamp that emits no hydrogen
    # light, so the honest answer is that the question could not be asked —
    # recorded as its own state rather than as a missing field.
    assert alignment_block["science_validation"] == "no-frame"
    assert "Balmer or Fulcher" in alignment_block["science_validation_note"]
    assert "science_residual_rms_nm" not in alignment_block
    # And the self-consistency number is still there beside it, unchanged.
    assert isinstance(alignment_block["rms_px"], float)


def test_existing_settings_can_be_deliberately_rewritten(tmp_path):
    """F21 item 11: "already exists" named the path but offered no way out.

    Refusing to clobber stays the default. The second, explicit press replaces
    the bundle — and stages and parses the whole new one first, so a failure
    part way through leaves the old files exactly where they were.
    """

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    configs = tmp_path / "configs"

    first = campaign.write_tomls(configs, "20250813_cmos", alignment)
    original = first["campaign"].read_text(encoding="utf-8")

    with pytest.raises(SnapshotError) as refused:
        campaign.write_tomls(configs, "20250813_cmos", alignment)
    assert "already exists" in str(refused.value)
    assert first["campaign"].read_text(encoding="utf-8") == original

    again = campaign.write_tomls(
        configs, "20250813_cmos", alignment, overwrite=True
    )

    assert campaign.toml_state is TomlState.GENERATED
    assert again["campaign"].exists()
    assert again["campaign"].read_text(encoding="utf-8") == original
    # No staging debris is left beside the bundle it replaced.
    assert [path.name for path in configs.iterdir()] == ["20250813_cmos"]


def _overwrite_ready(tmp_path: Path):
    """A campaign with one published bundle, ready to be rewritten."""

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)
    configs = tmp_path / "configs"
    published = campaign.write_tomls(configs, "20250813_cmos", alignment)
    return campaign, alignment, configs, published


def test_an_undoable_failed_overwrite_puts_the_old_bundle_back(tmp_path, monkeypatch):
    campaign, alignment, configs, published = _overwrite_ready(tmp_path)
    original = published["campaign"].read_text(encoding="utf-8")
    real_replace = os.replace
    moves = []

    def fail_the_publish(source, destination):
        # The first move parks the old bundle; the second publishes the new one.
        moves.append(destination)
        if len(moves) == 2:
            raise OSError("publish denied")
        return real_replace(source, destination)

    monkeypatch.setattr(campaign_module.os, "replace", fail_the_publish)

    with pytest.raises(OSError, match="publish denied"):
        campaign.write_tomls(configs, "20250813_cmos", alignment, overwrite=True)
    monkeypatch.undo()

    assert campaign.toml_state is TomlState.FAILED
    assert published["campaign"].read_text(encoding="utf-8") == original
    # Both the staging tree and the empty parking directory are gone.
    assert [path.name for path in configs.iterdir()] == ["20250813_cmos"]


def test_a_failed_overwrite_never_deletes_the_bundle_it_moved_aside(
    tmp_path, monkeypatch
):
    """The moved-aside bundle is the only copy of the old files.

    Parked inside the staging tree, it was deleted by the very cleanup that
    runs when publishing fails — so the one path that promised the old files
    survived was the path that removed them.
    """

    campaign, alignment, configs, published = _overwrite_ready(tmp_path)
    original = published["campaign"].read_text(encoding="utf-8")
    destination = published["campaign"].parent
    real_replace = os.replace

    def refuse_the_destination(source, destination_path):
        # Publishing fails, and putting the old bundle back fails with it.
        if Path(destination_path) == destination:
            raise OSError("destination is locked")
        return real_replace(source, destination_path)

    monkeypatch.setattr(campaign_module.os, "replace", refuse_the_destination)

    with pytest.raises(SnapshotError) as failure:
        campaign.write_tomls(configs, "20250813_cmos", alignment, overwrite=True)
    monkeypatch.undo()

    assert campaign.toml_state is TomlState.FAILED
    rescued = Path(str(failure.value).rsplit("kept at ", 1)[1].strip())
    assert rescued.is_absolute()
    assert rescued.is_dir()
    assert (rescued / "campaign.toml").read_text(encoding="utf-8") == original
    # And nothing on any path deleted them: every generated file of the old
    # bundle is still readable somewhere under the configuration root.
    survivors = {
        path.name
        for path in configs.rglob("*.toml")
        if path.parent.name == "20250813_cmos"
    }
    assert survivors == {"campaign.toml", "alignment.toml", "export.toml"}


def _peaked_frame(tmp_path, name, peak, *, width: int = 1):
    """A frame whose brightest raw pixel is *peak* counts.

    ``width`` widens that peak into a connected cluster, which is what a
    deliberately saturated frame carries: a lone full-scale pixel is read as a
    cosmic ray and never becomes the frame's peak.
    """

    columns = 60
    images = np.zeros((1, 20, columns), dtype=float)
    images[0, 10, 30 : 30 + width] = float(peak)
    spectra = (np.full(columns, float(peak) / 10.0),)
    return BenchFrame(tmp_path / name, images, images[0], spectra, {"ExposureTime": 0.1})


class TestABackgroundIsJudgedAsABackground:
    """The owner, deadpan, on a lamp background reading DIM with advice to
    expose it for twenty seconds: "Background is dim. Wow, genius."

    A background is shot with the light off. Dark is the whole of its job, so
    DIM is its correct state and raising its exposure would spoil the pair.
    What it can get wrong is the opposite, and that is what gets loud.
    """

    def _campaign(self, tmp_path):
        sources = _sources(tmp_path)
        return _campaign(tmp_path, sources), sources

    def test_a_dim_background_praises_itself_and_never_says_expose_longer(
        self, tmp_path
    ):
        triage = triage_exposure(_peaked_frame(tmp_path, "bg.sif", 200))
        assert triage.state is ExposureState.DIM

        for role in (
            MeasurementRole.LAMP_BACKGROUND,
            MeasurementRole.SPHERE_BACKGROUND,
        ):
            verdict = triage_for_role(triage, role, partner_peak=40000.0)

            assert verdict.label == "background"
            # Nothing to say and nothing to do: the bench does not explain
            # darkness to the person who shot the dark frame, and it offers no
            # action panel saying there is no action (owner: "I know that darks
            # are dark. Did you write it for yourself?").
            assert verdict.headline == ""
            assert verdict.next_action == ""
            assert not verdict.blocking
            assert verdict.is_usable
            assert "increase" not in verdict.advice.casefold()

    def test_a_background_approaching_its_signal_gets_loud(self, tmp_path):
        # Half the signal's counts: a leak, an open shutter, or the wrong file.
        triage = triage_exposure(_peaked_frame(tmp_path, "bg.sif", 20000))

        verdict = triage_for_role(
            triage, MeasurementRole.LAMP_BACKGROUND, partner_peak=40000.0
        )

        assert verdict.blocking
        assert not verdict.is_usable
        assert "50% of its lamp signal" in verdict.headline
        assert "not dark" in verdict.headline
        assert "light leak" in verdict.advice
        assert "reshoot" in verdict.next_action.casefold()

    def test_a_background_well_under_its_signal_stays_quiet(self, tmp_path):
        # A tenth of the signal is an ordinary, healthy background.
        triage = triage_exposure(_peaked_frame(tmp_path, "bg.sif", 4000))

        verdict = triage_for_role(
            triage, MeasurementRole.LAMP_BACKGROUND, partner_peak=40000.0
        )

        assert not verdict.blocking
        assert verdict.headline == ""

    def test_the_signal_frame_it_is_compared_against_is_its_own_partner(
        self, tmp_path
    ):
        campaign, sources = self._campaign(tmp_path)
        bright = _peaked_frame(tmp_path, "Ne-bright.sif", 40000)
        dark = _peaked_frame(tmp_path, "Ne-bg.sif", 300)
        for frame in (bright, dark):
            frame.path.write_bytes(b"sif\n")
        campaign.record_frame(bright)
        campaign.record_frame(dark)
        campaign.classify_file(
            bright.path, MeasurementRole.LAMP, lamp_family="Ne", frame=bright
        )
        campaign.classify_file(
            dark.path, MeasurementRole.LAMP_BACKGROUND, lamp_family="Ne", frame=dark
        )

        assert campaign.partner_peak(dark.path, MeasurementRole.LAMP_BACKGROUND) == (
            pytest.approx(40000.0)
        )
        # A lamp signal of another element is not this background's partner.
        assert campaign.partner_peak(bright.path, MeasurementRole.LAMP) is None
        verdict = campaign.role_triage(dark.path)
        assert verdict.label == "background"
        assert not verdict.blocking

    def test_a_background_is_judged_against_the_faintest_frame_it_serves(
        self, tmp_path
    ):
        """A lamp family is shot bright and dim; the same background serves both.

        Judged against the bright frame — deliberately at full scale — the leak
        threshold sits an order of magnitude above the dim series, and a
        background carrying most of the dim signal it is subtracted from passes
        as dark.
        """

        campaign, _sources = self._campaign(tmp_path)
        frames = {
            "bright": _peaked_frame(tmp_path, "Ne-bright.sif", 65535, width=4),
            "dim": _peaked_frame(tmp_path, "Ne-dim.sif", 3000),
            "background": _peaked_frame(tmp_path, "Ne-bg.sif", 2500),
        }
        roles = {
            "bright": MeasurementRole.LAMP,
            "dim": MeasurementRole.LAMP,
            "background": MeasurementRole.LAMP_BACKGROUND,
        }
        for key, frame in frames.items():
            frame.path.write_bytes(b"sif\n")
            campaign.record_frame(frame)
            # Spelled as the operator would type it; the family is normalized
            # on the way in and compared normalized on the way out.
            campaign.classify_file(
                frame.path, roles[key], lamp_family="ne", frame=frame
            )

        background = frames["background"].path
        assert campaign.partner_peak(
            background, MeasurementRole.LAMP_BACKGROUND
        ) == pytest.approx(3000.0)
        verdict = campaign.role_triage(background)
        assert verdict.blocking
        assert verdict.is_background
        assert "83% of its lamp signal" in verdict.headline

    def test_without_a_partner_a_background_is_still_read_as_one(self, tmp_path):
        triage = triage_exposure(_peaked_frame(tmp_path, "bg.sif", 200))

        verdict = triage_for_role(triage, MeasurementRole.SPHERE_BACKGROUND)

        assert verdict.label == "background"
        assert not verdict.blocking


# ---------------------------------------------------------------------------
# Two sphere pairs of the SAME physical response, shot at different exposures.
#
# Spheres are not usually shot at a shared exposure time, so a factor path that
# forgets to divide each sphere's counts by *its own* ExposureTime would report
# the exposure ratio as if it were a change in the lamp.  The fixture below is
# that whole claim in one frame set: identical count *rates*, exposures 0.1 s
# and 0.3 s, counts scaled to match -- the comparison must read 1.0.
# ---------------------------------------------------------------------------

_SPHERE_DIMW = 60
_SPHERE_DIMO = 55
_SPHERE_DV = 8  # Calibrations default order half-width
_SPHERE_ROWS = (20, 40)
_PREVIOUS_EXPOSURE_S = 0.1
_CANDIDATE_EXPOSURE_S = 0.3


def _sphere_count_rates() -> tuple[np.ndarray, np.ndarray]:
    """Two blaze-like humps that cross inside the orders' wavelength overlap."""

    columns = np.arange(_SPHERE_DIMW, dtype=float)
    lower = 1000.0 - 0.4 * (columns - 30.0) ** 2
    upper = 1000.0 - 0.4 * (columns - 22.0) ** 2
    return lower, upper


def _paint_order(image: np.ndarray, row: int, values: np.ndarray) -> None:
    rows = np.arange(row - _SPHERE_DV, row + _SPHERE_DV + 1)
    image[np.ix_(rows, np.arange(_SPHERE_DIMW))] = values / (2 * _SPHERE_DV + 1)


def _rate_frames() -> tuple[np.ndarray, np.ndarray]:
    """Sphere and background count *rates* (counts per second), not counts."""

    lower_rate, upper_rate = _sphere_count_rates()
    sphere = np.zeros((_SPHERE_DIMO, _SPHERE_DIMW))
    _paint_order(sphere, _SPHERE_ROWS[0], lower_rate)
    _paint_order(sphere, _SPHERE_ROWS[1], upper_rate)

    background = np.zeros((_SPHERE_DIMO, _SPHERE_DIMW))
    _paint_order(background, _SPHERE_ROWS[0], np.full(_SPHERE_DIMW, 60.0))
    _paint_order(background, _SPHERE_ROWS[1], np.full(_SPHERE_DIMW, 60.0))
    return sphere, background


@pytest.fixture
def sphere_pairs_at_two_exposures(monkeypatch: pytest.MonkeyPatch):
    """Serve both sphere pairs from one count-rate model at two exposures."""

    from echelle_spectra.tools import echelle as echelle_module

    sphere_rate, background_rate = _rate_frames()

    def read_image(fpth, spec="black", crop=(0, -1), exptime=1):
        name = Path(str(fpth)).name.casefold()
        exposure = _CANDIDATE_EXPOSURE_S if "0.3s" in name else _PREVIOUS_EXPOSURE_S
        rate = background_rate if "-bg" in name else sphere_rate
        images = (rate * exposure)[np.newaxis]
        info = {
            "NumberOfFrames": 1,
            "xbin": 1,
            "ybin": 1,
            "size": np.array([_SPHERE_DIMW, _SPHERE_DIMO]),
            "ExposureTime": exposure,
            "CycleTime": 1.0,
        }
        return images.copy(), info

    monkeypatch.setattr(echelle_module, "read_image", read_image)
    return read_image


def _sphere_fixture_sources(tmp_path: Path) -> dict[str, Path]:
    root = tmp_path / "campaign"
    root.mkdir(parents=True, exist_ok=True)

    pattern = root / "pattern.txt"
    np.savetxt(
        pattern,
        np.column_stack(
            [
                np.full(_SPHERE_DIMW, _SPHERE_ROWS[0], dtype=int),
                np.full(_SPHERE_DIMW, _SPHERE_ROWS[1], dtype=int),
            ]
        ),
        fmt="%d",
    )

    rows = [(27, pixel, 436.0 - 0.05 * pixel) for pixel in (2, 15, 30, 45, 58)]
    rows += [(28, pixel, 434.0 - 0.05 * pixel) for pixel in (2, 15, 30, 45, 58)]
    wavelength = root / "wavelength.txt"
    wavelength.write_text(
        "# synthetic lamp identification\n"
        "# order from to center wavelength\n"
        + "".join(
            f"{order:d} {pixel - 1:d} {pixel + 1:d} {pixel:d} {value:.6f}\n"
            for order, pixel, value in rows
        ),
        encoding="utf-8",
    )

    integral = root / "integral.txt"
    micrometres = np.linspace(0.425, 0.440, 9)
    np.savetxt(
        integral,
        np.column_stack([micrometres, np.full(micrometres.size, 1.0e-2)]),
        fmt="%.8f",
    )

    names = {
        "sphere": "sphere-0.3s.sif",
        "sphere_background": "sphere-0.3s-bg.sif",
        "previous_sphere": "sphere-0.1s.sif",
        "previous_sphere_background": "sphere-0.1s-bg.sif",
    }
    answer = {"pattern": pattern, "wavelength": wavelength, "integral": integral}
    for role, name in names.items():
        path = root / name
        path.write_bytes(b"synthetic sphere frame\n")
        answer[role] = path
    return answer


class TestSphereFactorsCarryEachSpheresOwnExposure:
    """A reported sphere ratio must survive a pair that differs only in exposure."""

    def _campaign(self, sources: dict[str, Path]) -> CalibrationCampaignSession:
        campaign = CalibrationCampaignSession(
            pattern_source=sources["pattern"],
            wavelength_source=sources["wavelength"],
            integral_source=sources["integral"],
            previous_sphere=sources["previous_sphere"],
            previous_sphere_background=sources["previous_sphere_background"],
        )
        campaign.classify_file(sources["sphere"], MeasurementRole.SPHERE)
        campaign.classify_file(
            sources["sphere_background"], MeasurementRole.SPHERE_BACKGROUND
        )
        return campaign

    def test_same_response_at_0_3_s_and_0_1_s_compares_at_one(
        self, tmp_path, sphere_pairs_at_two_exposures
    ):
        sources = _sphere_fixture_sources(tmp_path)
        campaign = self._campaign(sources)

        # The real engine, not a stub: this is the calculator the bench uses.
        comparison = campaign.compute_sphere_comparison()

        assert comparison.state is ComparisonState.READY, comparison.reason
        assert comparison.sample_count >= 20
        # Anything but 1.0 here is the 3x exposure ratio leaking into a number
        # the operator would read as a change in the integrating sphere.
        assert comparison.median_ratio == pytest.approx(1.0, rel=1e-12)
        assert comparison.p05_ratio == pytest.approx(1.0, rel=1e-12)
        assert comparison.p95_ratio == pytest.approx(1.0, rel=1e-12)

    def test_the_two_factor_curves_agree_column_by_column(
        self, tmp_path, sphere_pairs_at_two_exposures
    ):
        """The median agreeing is not enough; every shared column must agree."""

        sources = _sphere_fixture_sources(tmp_path)
        campaign = self._campaign(sources)

        comparison = campaign.compute_sphere_comparison()

        assert comparison.candidate is not None and comparison.previous is not None
        # Same response, three times the counts: without each file's own
        # exposure the candidate factor would come out three times small.
        np.testing.assert_allclose(
            comparison.candidate.factors_wmsr,
            comparison.previous.factors_wmsr,
            rtol=1e-12,
        )

    def test_a_shot_is_calibrated_by_its_own_exposure_too(
        self, tmp_path, sphere_pairs_at_two_exposures
    ):
        """The applied-factor path divides by the shot's exposure, not the sphere's."""

        from echelle_spectra.tools.loader import build_calibration, load_spectrum

        sources = _sphere_fixture_sources(tmp_path)
        calibration = build_calibration(
            sources["pattern"].parent,
            "CMOS",
            calibration_files={
                "orders": str(sources["pattern"]),
                "wavelength": str(sources["wavelength"]),
                "sphr": str(sources["sphere"]),
                "bkgr": str(sources["sphere_background"]),
                "integral": str(sources["integral"]),
            },
        )

        # The same light, recorded once for 0.1 s and once for 0.3 s.
        short = load_spectrum(sources["previous_sphere"], calibration=calibration)
        long_ = load_spectrum(sources["sphere"], calibration=calibration)

        assert short.info["ExposureTime"] == pytest.approx(_PREVIOUS_EXPOSURE_S)
        assert long_.info["ExposureTime"] == pytest.approx(_CANDIDATE_EXPOSURE_S)
        short_counts = np.asarray(short.counts, dtype=float)
        long_counts = np.asarray(long_.counts, dtype=float)
        assert not np.allclose(short_counts, long_counts)
        np.testing.assert_allclose(
            np.asarray(long_.wmsr, dtype=float),
            np.asarray(short.wmsr, dtype=float),
            rtol=1e-12,
        )


#: The owner's real 2019 folder, in its own dialect: "IS" is how he says
#: integrating sphere, and the lamp was shot signal-only because in 2019 nobody
#: shot lamp backgrounds.  Neither fact can be changed now — the frames are
#: seven years old — so the bench has to read both.
_REAL_2019_NAMES = (
    "IS-1s.sif",
    "IS_bg.sif",
    "Ne_1s_10fr.sif",
)


def _real_2019_folder(tmp_path: Path) -> dict[str, Path]:
    folder = tmp_path / "20190314"
    folder.mkdir()
    sources = {}
    for name in _REAL_2019_NAMES:
        path = folder / name
        path.write_bytes((name + "\n").encode())
        sources[name] = path
    for name, text in (
        ("pattern.txt", "pattern\n"),
        ("wavelength.txt", _CURATED_TABLE),
        ("integral.txt", "integral\n"),
    ):
        path = folder / name
        path.write_text(text, encoding="utf-8")
        sources[name] = path
    return sources


def _2019_campaign(sources: dict[str, Path]) -> CalibrationCampaignSession:
    return CalibrationCampaignSession(
        pattern_source=sources["pattern.txt"],
        wavelength_source=sources["wavelength.txt"],
        integral_source=sources["integral.txt"],
        suggested_lamps=("Ne",),
    )


def _signal_only_campaign(tmp_path: Path):
    """The owner's live state: sphere pair, one lamp signal, no lamp background."""

    sources = _real_2019_folder(tmp_path)
    campaign = _2019_campaign(sources)
    campaign.classify_file(sources["IS-1s.sif"], MeasurementRole.SPHERE)
    campaign.classify_file(sources["IS_bg.sif"], MeasurementRole.SPHERE_BACKGROUND)
    campaign.classify_file(
        sources["Ne_1s_10fr.sif"], MeasurementRole.LAMP, lamp_family="Ne"
    )
    campaign.compute_sphere_comparison(_calculator)
    return campaign, sources


# ---------------------------------------------------------------------------
# Xenon is a real lamp
# ---------------------------------------------------------------------------


def test_xenon_no_longer_walls_at_the_lamp_preset_registry():
    """Regression pin for the wall: the three refusals a Xe lamp used to hit.

    Assigning the lamp always worked — ``normalize_lamp_name`` accepts any
    path-safe name — and everything downstream of it refused: no ``xe`` preset,
    no ``XeI``/``XeII`` species, no ``xe`` catalog family.  The operator could
    name the lamp and then had nothing to fit against.
    """

    assert normalize_lamp_name("Xe") == "Xe"
    assert normalize_lamp_name("xenon") == "Xe"
    assert lamp_species(["xe"]) == ("XeI", "XeII")
    assert normalize_species_key("Xe I") == "XeI"
    assert normalize_species_key("XeII") == "XeII"
    assert campaign_module.catalog_family_for_lamp("Xe") == "xe"
    assert load_line_table("xe")


def test_an_assigned_xenon_lamp_gets_expected_lines_and_a_fit_reference(tmp_path):
    """The panel and the click-to-fit reach the new cache through the ordinary path."""

    sources = _curated_sources(tmp_path)
    rows = load_wavelength_table(sources["wavelength.txt"])
    order = rows[0].order_idx

    listed = catalog_lines_for_order(rows, order, "Xe")
    assert listed, "an assigned Xe lamp draws no expected lines"
    assert {entry.line.family for entry in listed} == {"xe"}
    assert all(entry.line.species.startswith("Xe ") for entry in listed)
    assert all(entry.line.relative_intensity is not None for entry in listed)
    # Same door as any other lamp: nothing about Xe is special-cased.
    assert catalog_lines_for_order(rows, order, "xenon") == listed

    reference = lamp_reference_set("Xe", rows)
    assert reference.lamp == "Xe"
    assert reference.catalog_family == "xe"
    assert reference.catalog_label == "Xe"
    assert reference.species == ("XeI", "XeII")


def test_a_xenon_lamp_says_the_vetted_set_carries_no_xenon_rows(tmp_path):
    """F19's pedigree statement, not an error: nobody has vetted a xenon line.

    Before this packet the answer was "no line catalog for Xe", which blamed
    the package for a gap that belongs to the wavelength table.  The catalog
    exists now; what is missing is anybody's signature on a xenon measurement,
    and that is what the bench says.
    """

    sources = _curated_sources(tmp_path)
    rows = load_wavelength_table(sources["wavelength.txt"])
    reference = lamp_reference_set("Xe", rows)

    assert reference.state is ReferenceState.NO_ROWS
    assert not reference.is_referenceable
    assert "packaged Xe catalog" in reference.message
    assert "XeI, XeII" in reference.message
    assert reference.lines == ()
    # And it contributes nothing to the panel rather than raising there.
    assert expected_lines_for_order(reference, rows[0].order_idx) == ()


# ---------------------------------------------------------------------------
# "IS" means integrating sphere
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, role",
    [
        ("IS-1s.sif", MeasurementRole.SPHERE),
        ("IS_1s.sif", MeasurementRole.SPHERE),
        ("is-0.5s-x3.sif", MeasurementRole.SPHERE),
        ("IS.sif", MeasurementRole.SPHERE),
        ("IS-bg.sif", MeasurementRole.SPHERE_BACKGROUND),
        ("IS_bg.sif", MeasurementRole.SPHERE_BACKGROUND),
        ("is_background.sif", MeasurementRole.SPHERE_BACKGROUND),
    ],
)
def test_the_is_prefix_prefills_the_sphere_roles(name, role):
    """The owner's shorthand joins the suggestion vocabulary, as a pre-fill only."""

    suggestion = suggest_file_roles(name)
    assert suggestion.roles == (role,)
    assert suggestion.is_unambiguous
    assert "integrating-sphere" in suggestion.reason
    assert "confirm" in suggestion.reason


@pytest.mark.parametrize(
    "name", ["isotope-scan.sif", "island.sif", "issue-42.sif", "Ne-is-bright.sif"]
)
def test_a_word_that_merely_starts_with_is_is_not_a_sphere_frame(name):
    """Only the head of the name, and only where no letter follows it."""

    suggestion = suggest_file_roles(name)
    sphere_named = (
        suggestion.is_unambiguous and suggestion.roles[0] is MeasurementRole.SPHERE
    )
    assert not sphere_named


def test_a_folder_of_is_xe_and_ne_names_assigns_itself(tmp_path):
    """The unanimous rule reads the new vocabulary like any other suggestion."""

    sources = _curated_sources(tmp_path)
    campaign = _campaign(tmp_path, sources)
    names = (
        "IS-1s.sif",
        "IS-1s-bg.sif",
        "Xe-2s.sif",
        "Xe-2s_bg.sif",
        "Ne-0.5s.sif",
        "Ne-0.5s_bg.sif",
    )
    for name in names:
        (tmp_path / name).write_bytes((name + "\n").encode())
        campaign.record_frame(_triage_frame(tmp_path, name, _bright_line))

    applied = {
        record.path.name: record for record in campaign.apply_unanimous_suggestions()
    }

    assert set(applied) == set(names)
    assert applied["IS-1s.sif"].role is MeasurementRole.SPHERE
    assert applied["IS-1s-bg.sif"].role is MeasurementRole.SPHERE_BACKGROUND
    assert applied["Xe-2s.sif"].role is MeasurementRole.LAMP
    assert applied["Xe-2s.sif"].lamp_family == "Xe"
    assert applied["Xe-2s_bg.sif"].role is MeasurementRole.LAMP_BACKGROUND
    assert applied["Xe-2s_bg.sif"].lamp_family == "Xe"
    assert applied["Ne-0.5s.sif"].lamp_family == "Ne"
    assert campaign.assigned_lamps == ("Ne", "Xe")


# ---------------------------------------------------------------------------
# A signal-only lamp folder can be saved
# ---------------------------------------------------------------------------


def test_a_signal_only_lamp_campaign_saves_and_records_the_absent_background(tmp_path):
    """The 2019 folder has no lamp background and never will; the save proceeds.

    Demanding a complete lamp pair meant a historical folder could solve an
    alignment and then never write it down.  The sphere pair is still required
    — the absolute factor is computed from the difference of those two frames —
    but a lamp background is an improvement to a fit, not a licence to record
    the fit that was made without it.
    """

    alignment = _aligned_session(tmp_path)
    campaign, _named = _signal_only_campaign(tmp_path)

    assert campaign.lamps_without_background() == ("Ne",)

    paths = campaign.write_tomls(tmp_path / "configs", "20190314_cmos", alignment)

    assert campaign.toml_state is TomlState.GENERATED
    assert campaign.save_state is SaveState.READY
    assert campaign.ready_for_snapshot("20190314_cmos", alignment)

    text = paths["campaign"].read_text(encoding="utf-8")
    assert "background_shot = false" in text
    assert "no Ne lamp background frame exists in this campaign" in text
    assert "unsubtracted" in text
    payload = tomllib.loads(text)
    lamp = [row for row in payload["measurements"] if row["role"] == "lamp"]
    assert len(lamp) == 1
    assert lamp[0]["lamp_family"] == "Ne"
    assert lamp[0]["background_shot"] is False

    # And the snapshot the settings unlock is reachable, not merely composed.
    snapshot = campaign.save_snapshot(
        tmp_path / "calibrations",
        snapshot_id="20190314_cmos",
        detector="cmos",
        alignment=alignment,
    )
    assert snapshot.snapshot_id == "20190314_cmos"
    assert campaign.save_state is SaveState.VALIDATED


def test_a_complete_lamp_pair_still_records_its_background(tmp_path):
    """The pair-complete path is unchanged, and now says so in the record."""

    sources = _curated_sources(tmp_path)
    alignment = _aligned_session(tmp_path)
    campaign = _campaign(tmp_path, sources)
    _classify_complete(campaign, sources, alignment.frame)
    campaign.compute_sphere_comparison(_calculator)

    paths = campaign.write_tomls(tmp_path / "configs", "20250813_cmos", alignment)

    assert campaign.lamps_without_background() == ()
    payload = tomllib.loads(paths["campaign"].read_text(encoding="utf-8"))
    lamp = [row for row in payload["measurements"] if row["role"] == "lamp"]
    assert lamp and all(row["background_shot"] is True for row in lamp)
    assert all("background_note" not in row for row in lamp)


def test_a_missing_lamp_background_is_a_note_and_never_the_next_step(tmp_path):
    """The checklist states the fact and blocks nothing with it."""

    alignment = _aligned_session(tmp_path)
    campaign, _named = _signal_only_campaign(tmp_path)

    checklist = {item.key: item for item in campaign.checklist(alignment)}
    row = checklist["lamp-Ne-background"]
    assert row.state is ChecklistState.SUGGESTION
    assert not row.blocking
    assert "signal only" in row.detail
    assert "no background shot" in row.detail
    assert checklist["lamp-Ne-signal"].state is ChecklistState.DONE
    # Nothing blocking is left pointing at a frame that does not exist.
    blocking = [
        item
        for item in campaign.checklist(alignment)
        if item.blocking and item.state is not ChecklistState.DONE
    ]
    assert all("background" not in item.label for item in blocking)


def test_a_refused_save_leaves_no_folders_standing(tmp_path):
    """Nothing is created on disk until the bundle composes.

    An empty ``calibrations/configs/`` beside a bench that had just refused was
    the only trace the refusal left anywhere, and it read as a half-written
    save rather than as no save at all.
    """

    sources = _real_2019_folder(tmp_path)
    campaign = _2019_campaign(sources)
    campaign.classify_file(sources["IS-1s.sif"], MeasurementRole.SPHERE)
    alignment = _aligned_session(tmp_path)
    root = tmp_path / "calibrations" / "configs"

    with pytest.raises(SnapshotError) as refusal:
        campaign.write_tomls(root, "20190314_cmos", alignment)

    assert "sphere pair" in str(refusal.value)
    assert "one lamp signal" in str(refusal.value)
    assert not root.exists()
    assert not root.parent.exists()
    assert campaign.toml_state is TomlState.FAILED
