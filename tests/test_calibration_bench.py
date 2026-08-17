"""UI-independent tests for the live calibration bench workflow."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from echelle_spectra.calibration_bench import (
    AlignmentState,
    Anchor,
    BenchFrame,
    CalibrationBenchSession,
    FileLoadState,
    FrameLoader,
    StableFileState,
    StableSifWatcher,
    band_center_offsets,
)
from echelle_spectra.tools.calibration_alignment import (
    CalibrationTableLine,
    DetectorWindowSaturation,
    LineCentroidFit,
    load_wavelength_table,
)


def _pattern(columns: int = 120) -> np.ndarray:
    pattern = np.empty((columns, 3), dtype=float)
    pattern[:, 0] = 12
    pattern[:, 1] = 30
    pattern[:, 2] = 48
    return pattern


def _line(order: int, center: float, wavelength: float) -> CalibrationTableLine:
    return CalibrationTableLine(
        order,
        center - 5,
        center + 5,
        center,
        wavelength,
        "NeI",
        "ok",
    )


def _frame(
    tmp_path: Path,
    *,
    centers: tuple[float, float, float] = (32.0, 72.0, 92.0),
    saturated: tuple[int, int] | None = None,
) -> BenchFrame:
    columns = 120
    x = np.arange(columns, dtype=float)
    spectra = tuple(
        10.0 + 500.0 * np.exp(-0.5 * ((x - center) / 1.6) ** 2)
        for center in centers
    )
    images = np.zeros((2, 64, columns), dtype=float)
    pattern = _pattern(columns)
    for order, spectrum in enumerate(spectra):
        for column, value in enumerate(spectrum):
            images[:, int(pattern[column, order]), column] = value
    if saturated is not None:
        order, column = saturated
        images[0, int(pattern[column, order]), column] = 65535
    return BenchFrame(
        path=tmp_path / "fixture.sif",
        images=images,
        detector_image=images.mean(axis=0),
        order_spectra=spectra,
        metadata={"NumberOfFrames": 2},
    )


def _anchor(line: CalibrationTableLine, measured: float) -> Anchor:
    fit = LineCentroidFit(line, measured, 1.5, 100.0, 5.0, 20.0, True)
    saturation = DetectorWindowSaturation(line, 100.0, 20, 0, 0.0)
    return Anchor(line, fit, saturation)


def test_watcher_requires_repeated_unchanged_observations(tmp_path):
    watcher = StableSifWatcher(tmp_path, required_unchanged_polls=2)
    assert watcher.poll().state is StableFileState.EMPTY

    source = tmp_path / "lamp.SIF"
    source.write_bytes(b"incomplete")
    first = watcher.poll()
    second = watcher.poll()
    third = watcher.poll()

    assert first.state is StableFileState.CHANGING
    assert first.unchanged_polls == 1
    assert second.state is StableFileState.STABLE
    assert second.ready_path == source
    assert third.state is StableFileState.ALREADY_EMITTED
    assert third.ready_path is None


def test_watcher_restarts_stability_after_growth_and_prefers_newest(tmp_path):
    first_path = tmp_path / "first.sif"
    first_path.write_bytes(b"a")
    watcher = StableSifWatcher(tmp_path, required_unchanged_polls=2)
    watcher.poll()
    assert watcher.poll().ready_path == first_path

    first_path.write_bytes(b"still-growing")
    assert watcher.poll().state is StableFileState.CHANGING
    assert watcher.poll().ready_path == first_path

    newest = tmp_path / "second.sif"
    newest.write_bytes(b"complete")
    newest.touch()
    candidate = watcher.poll()
    ready = watcher.poll()
    assert candidate.path == newest
    assert candidate.state is StableFileState.CHANGING
    assert ready.ready_path == newest


def test_watcher_reports_missing_folder_without_raising(tmp_path):
    watcher = StableSifWatcher(tmp_path / "missing")
    result = watcher.poll()
    assert result.state is StableFileState.FAILED
    assert result.reason


def test_watcher_enforces_minimum_age_after_file_stops_changing(tmp_path):
    source = tmp_path / "young.sif"
    source.write_bytes(b"complete")
    modified_ns = source.stat().st_mtime_ns
    watcher = StableSifWatcher(
        tmp_path,
        required_unchanged_polls=2,
        minimum_age_s=2.0,
    )

    watcher.poll(now_ns=modified_ns + 500_000_000)
    too_young = watcher.poll(now_ns=modified_ns + 1_000_000_000)
    ready = watcher.poll(now_ns=modified_ns + 2_000_000_000)

    assert too_young.state is StableFileState.CHANGING
    assert ready.ready_path == source


def test_anchor_mutations_drive_alignment_states(tmp_path):
    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0))
    session = CalibrationBenchSession(_pattern(), lines, minimum_snr=3.0)
    session.accept_frame(_frame(tmp_path))
    assert session.alignment_state is AlignmentState.EMPTY

    first = session.fit_anchor_at(0, 32.0)
    assert first.accepted
    assert session.alignment_state is AlignmentState.COLLECTING
    assert len(session.anchors) == 1

    second = session.fit_anchor_at(1, 72.0)
    assert second.accepted
    assert session.alignment_state is AlignmentState.ALIGNED
    assert session.transform is not None
    np.testing.assert_allclose(session.transform.dx_px, 2.0, atol=0.05)
    assert session.rms_px is not None
    assert session.rms_px < 0.05
    assert len(session.residuals) == 2

    assert session.remove_anchor(second.anchor.key)
    assert session.alignment_state is AlignmentState.COLLECTING
    assert session.transform is None
    session.clear_anchors()
    assert session.alignment_state is AlignmentState.EMPTY


def test_replacing_anchor_recomputes_rigid_fit(tmp_path):
    first_line = _line(0, 30.0, 600.0)
    second_line = _line(1, 70.0, 610.0)
    session = CalibrationBenchSession(_pattern(), (first_line, second_line))
    session.accept_frame(_frame(tmp_path))
    session.upsert_anchor(_anchor(first_line, 31.0))
    session.upsert_anchor(_anchor(second_line, 71.0))
    assert session.transform is not None
    assert session.transform.dx_px == pytest.approx(1.0)

    session.upsert_anchor(_anchor(second_line, 73.0))
    assert len(session.anchors) == 2
    assert session.transform is not None
    assert session.rms_px is not None
    assert session.rms_px > 0


def test_saturated_and_low_signal_clicks_are_refused_then_recover(tmp_path):
    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0))
    session = CalibrationBenchSession(_pattern(), lines, minimum_snr=3.0)
    session.accept_frame(_frame(tmp_path, saturated=(0, 32)))
    saturated = session.fit_anchor_at(0, 32.0)
    assert not saturated.accepted
    # The rejection names the line rather than ordering a reshoot: a dim-series
    # lamp frame is *meant* to saturate its strong lines (F14 item 2).
    assert "saturated on the raw detector" in saturated.reason
    assert "600.000 nm" in saturated.reason
    assert "lower exposure" not in saturated.reason.casefold()
    assert not session.anchors

    dim = _frame(tmp_path, centers=(115.0, 115.0, 115.0))
    session.accept_frame(dim)
    failed = session.fit_anchor_at(0, 30.0)
    assert not failed.accepted
    assert not session.anchors

    session.accept_frame(_frame(tmp_path))
    assert session.fit_anchor_at(0, 32.0).accepted
    assert session.fit_anchor_at(1, 72.0).accepted
    assert session.alignment_state is AlignmentState.ALIGNED


def _dim_series_frame(tmp_path: Path, frames: int = 3) -> BenchFrame:
    """A dim-series lamp acquisition: strong line clipped, weak lines clean.

    This is the exposure the owner shoots on purpose — the bright/dim pair
    exists so the strong lines saturate while the weak ones finally emerge.
    """

    columns = 120
    x = np.arange(columns, dtype=float)
    pattern = _pattern(columns)
    strong = 10.0 + 90000.0 * np.exp(-0.5 * ((x - 32.0) / 1.6) ** 2)
    weak_a = 10.0 + 900.0 * np.exp(-0.5 * ((x - 72.0) / 1.6) ** 2)
    weak_b = 10.0 + 700.0 * np.exp(-0.5 * ((x - 92.0) / 1.6) ** 2)
    spectra = (np.minimum(strong, 65535.0), weak_a, weak_b)
    images = np.zeros((frames, 64, columns), dtype=float)
    for order, spectrum in enumerate(spectra):
        for column, value in enumerate(spectrum):
            images[:, int(pattern[column, order]), column] = value
    per_frame = tuple(tuple(np.array(row, dtype=float) for row in spectra) for _ in range(frames))
    return BenchFrame(
        path=tmp_path / "Ne-0.1s-x3-dimm-lines.sif",
        images=images,
        detector_image=images.mean(axis=0),
        order_spectra=spectra,
        metadata={"NumberOfFrames": frames},
        frame_order_spectra=per_frame,
    )


def test_per_anchor_saturation_rejects_only_the_saturated_lines(tmp_path):
    """F14 item 2: a dim-series frame still yields its unsaturated anchors.

    The frame as a whole is saturated, which is the point of shooting it. The
    fit is guarded line by line, so the clipped strong line is refused while
    the weak lines it was shot for fit and solve the alignment.
    """

    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0), _line(2, 90.0, 620.0))
    session = CalibrationBenchSession(_pattern(), lines, minimum_snr=3.0)
    session.accept_frame(_dim_series_frame(tmp_path))

    clipped = session.fit_anchor_at(0, 32.0)
    assert not clipped.accepted
    assert "saturated on the raw detector" in clipped.reason
    assert not session.anchors

    assert session.fit_anchor_at(1, 72.0).accepted
    assert session.fit_anchor_at(2, 92.0).accepted
    assert session.alignment_state is AlignmentState.ALIGNED
    assert len(session.anchors) == 2


def test_a_three_frame_acquisition_fits_the_mean_or_one_frame(tmp_path):
    """F14 item 4: the owner shoots three frames; both readings are offered."""

    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0), _line(2, 90.0, 620.0))
    session = CalibrationBenchSession(_pattern(), lines, minimum_snr=3.0)
    session.accept_frame(_dim_series_frame(tmp_path, frames=3))

    assert session.frame_count == 3
    assert session.selected_frame is None
    assert session.frame_choice_label() == "mean of 3 frame(s)"
    assert session.active_images().shape[0] == 3

    assert session.fit_anchor_at(1, 72.0).accepted
    assert session.fit_anchor_at(2, 92.0).accepted
    assert session.alignment_state is AlignmentState.ALIGNED

    # A centroid belongs to the spectrum it was measured on, so switching the
    # reading drops the anchors instead of mixing two of them.
    session.set_selected_frame(1)
    assert session.frame_choice_label() == "frame 2 of 3"
    assert session.active_images().shape[0] == 1
    assert not session.anchors
    assert session.alignment_state is AlignmentState.EMPTY

    assert session.fit_anchor_at(1, 72.0).accepted
    np.testing.assert_allclose(
        session.active_order_spectra()[1], session.frame.frame_order_spectra[1][1]
    )
    with pytest.raises(IndexError):
        session.set_selected_frame(7)

    session.set_selected_frame(None)
    assert session.selected_frame is None
    assert not session.anchors


def test_alignment_failure_keeps_anchors_and_recovers_after_mutation(tmp_path):
    valid_a = _line(0, 30.0, 600.0)
    valid_b = _line(1, 70.0, 610.0)
    invalid = _line(9, 50.0, 620.0)
    session = CalibrationBenchSession(_pattern(), (valid_a, valid_b, invalid))
    session.accept_frame(_frame(tmp_path))
    session.upsert_anchor(_anchor(valid_a, 32.0))
    session.upsert_anchor(_anchor(invalid, 52.0))
    assert session.alignment_state is AlignmentState.FAILED
    assert "outside pattern" in session.last_error
    assert len(session.anchors) == 2

    session.remove_anchor((invalid.order_idx, invalid.center_pixel, invalid.wavelength_nm))
    session.upsert_anchor(_anchor(valid_b, 72.0))
    assert session.alignment_state is AlignmentState.ALIGNED
    assert not session.last_error


def test_file_failure_preserves_last_good_frame_and_next_file_recovers(tmp_path):
    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0))
    session = CalibrationBenchSession(_pattern(), lines)
    good = _frame(tmp_path)
    session.accept_frame(good)
    session.upsert_anchor(_anchor(lines[0], 32.0))

    def fail_loader(_path):
        raise OSError("truncated")

    assert not session.load_file(tmp_path / "bad.sif", fail_loader)
    assert session.file_state is FileLoadState.FAILED
    assert session.frame is good
    assert len(session.anchors) == 1

    recovered = BenchFrame(
        tmp_path / "recovered.sif",
        good.images,
        good.detector_image,
        good.order_spectra,
        good.metadata,
    )
    assert session.load_file(recovered.path, lambda _path: recovered)
    assert session.file_state is FileLoadState.LOADED
    assert session.frame is recovered
    assert session.alignment_state is AlignmentState.EMPTY
    assert not session.anchors
    assert not session.last_error


def test_packaged_fixture_loads_with_accepted_2025_pattern():
    package = Path(__file__).parents[1] / "src" / "echelle_spectra" / "resources"
    calibration_dir = package / "calibration_files"
    pattern = np.loadtxt(calibration_dir / "pattern_CMOS_20250926.txt", dtype=int)
    loader = FrameLoader(pattern)

    frame = loader(calibration_dir / "ThAr-0.3s-x3_20240305.sif")

    assert frame.detector_image.shape == (2160, 2560)
    assert len(frame.order_spectra) == pattern.shape[1] == 29
    assert all(spectrum.size > 2000 for spectrum in frame.order_spectra)
    assert np.isfinite(frame.detector_image).all()

    lines = load_wavelength_table(
        calibration_dir
        / "alignments"
        / "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    )
    session = CalibrationBenchSession(pattern, lines, minimum_snr=3.0)
    session.accept_frame(frame)
    first = lines[0]
    second = next(line for line in lines if line.order_idx != first.order_idx)
    assert session.fit_anchor_at(first.order_idx, first.center_pixel).accepted
    assert session.fit_anchor_at(second.order_idx, second.center_pixel).accepted
    assert session.alignment_state is AlignmentState.ALIGNED
    assert session.rms_px is not None and np.isfinite(session.rms_px)


def _background_frame(tmp_path: Path, level: float = 40.0) -> BenchFrame:
    """A lineless hump: what a lamp background actually looks like."""

    columns = 120
    x = np.arange(columns, dtype=float)
    hump = level + 0.02 * (x - 60.0) ** 2
    spectra = tuple(np.array(hump) for _ in range(3))
    images = np.tile(hump, (2, 64, 1))
    return BenchFrame(
        tmp_path / "Ne_bg.sif",
        images,
        images.mean(axis=0),
        spectra,
        {"ExposureTime": 0.1},
    )


def test_the_fit_measures_the_signal_minus_its_assigned_background(tmp_path):
    """F16 item 5: what echelle-align measures is signal minus background."""

    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))
    session.accept_frame(_frame(tmp_path))
    raw = session.active_order_spectra()[0].copy()
    assert session.background_path is None
    assert session.background_label() == "no background subtracted"

    background = _background_frame(tmp_path)
    session.use_background_frame(background)

    assert session.background_path == background.path
    assert "Ne_bg.sif" in session.background_label()
    corrected = session.active_order_spectra()[0]
    assert np.allclose(corrected, raw - background.order_spectra[0])
    # The peak survives the subtraction; only the hump under it goes.
    assert corrected.argmax() == raw.argmax()

    # Dropping the pairing restores the raw signal.
    session.use_background_frame(None)
    assert np.allclose(session.active_order_spectra()[0], raw)


def test_changing_the_background_drops_anchors_measured_without_it(tmp_path):
    """A centroid belongs to the spectrum it was measured on."""

    session = CalibrationBenchSession(
        _pattern(), (_line(0, 32.0, 585.2), _line(1, 72.0, 640.1))
    )
    session.accept_frame(_frame(tmp_path))
    assert session.fit_anchor_at(0, 32.0).accepted
    assert session.fit_anchor_at(1, 72.0).accepted
    assert session.alignment_state is AlignmentState.ALIGNED

    session.use_background_frame(_background_frame(tmp_path))
    assert not session.anchors
    assert session.alignment_state is AlignmentState.EMPTY

    # Re-declaring the same background is not a change and costs nothing.
    assert session.fit_anchor_at(0, 32.0).accepted
    session.use_background_frame(_background_frame(tmp_path))
    assert session.anchors


def test_a_new_acquisition_forgets_the_previous_pairing(tmp_path):
    """A new frame may be another lamp; its background is re-established."""

    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))
    session.accept_frame(_frame(tmp_path))
    session.use_background_frame(_background_frame(tmp_path))
    assert session.background_path is not None

    session.accept_frame(_frame(tmp_path))
    assert session.background_path is None


def test_a_background_of_the_wrong_shape_is_refused_not_applied(tmp_path):
    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))
    session.accept_frame(_frame(tmp_path))
    wrong = _background_frame(tmp_path)
    wrong = BenchFrame(
        wrong.path,
        wrong.images,
        wrong.detector_image,
        wrong.order_spectra[:2],
        wrong.metadata,
    )
    with pytest.raises(ValueError, match="order counts"):
        session.use_background_frame(wrong)
    assert session.background_path is None


def test_the_fit_lands_where_the_lamps_own_lines_are(tmp_path):
    """F16 item 7: a lamp lives in a few orders out of twenty-nine."""

    from echelle_spectra.calibration_campaign import lamp_reference_set

    lines = (
        _line(0, 32.0, 585.2),
        _line(1, 40.0, 590.1),
        _line(1, 72.0, 594.4),
        _line(1, 92.0, 597.5),
        _line(2, 50.0, 600.2),
    )
    session = CalibrationBenchSession(_pattern(), lines)
    session.accept_frame(_frame(tmp_path))
    # With no lamp assigned there is nothing to land on and nothing is claimed.
    assert session.best_reference_order() is None

    session.use_lamp_reference(lamp_reference_set("Ne", lines))
    assert session.best_reference_order() == 1

    # A lamp whose catalog this table never carries lands nowhere, quietly.
    session.use_lamp_reference(lamp_reference_set("Hg", lines))
    assert session.best_reference_order() is None


# --- The solved correction moves the rows a click is measured against -------
#
# The 2019 folder is what named this: its transform solved dx +18.02 px, the
# sticks stayed on the 2024 table's pixels, and every click after the solve
# opened its +/-18 px window with the real peak outside it.  The fixture shifts
# by 20 px for the same reason the screenshot did it at 18.02: past the window
# radius, the fit has only the line's flank to work with.


def _shifted_session(tmp_path: Path, shift_px: float = 20.0) -> CalibrationBenchSession:
    """A frame whose lines all sit *shift_px* right of the curated table."""

    lines = (_line(0, 30.0, 600.0), _line(1, 70.0, 610.0), _line(2, 90.0, 620.0))
    session = CalibrationBenchSession(_pattern(), lines, minimum_snr=3.0)
    session.accept_frame(
        _frame(tmp_path, centers=(30.0 + shift_px, 70.0 + shift_px, 90.0 + shift_px))
    )
    return session


def _solve_two_anchors(session: CalibrationBenchSession, shift_px: float = 20.0) -> None:
    """Bootstrap the way an operator does: click the peaks that are visible."""

    assert session.fit_anchor_at(0, 30.0 + shift_px).accepted
    assert session.fit_anchor_at(1, 70.0 + shift_px).accepted
    assert session.transform is not None
    assert session.transform.dx_px == pytest.approx(shift_px, abs=0.05)


def test_before_a_solve_the_rows_sit_where_the_table_put_them(tmp_path):
    """First anchors are found against base pixels; nothing else is known yet."""

    session = _shifted_session(tmp_path)

    assert session.transform is None
    assert [line.center_pixel for line in session.lines_for_order(2)] == [90.0]
    assert session.corrected_rows(session.lines) == session.lines

    # And the click is still what centres the window, which is the only reason
    # a table 20 px out can be anchored by hand at all.
    result = session.fit_anchor_at(2, 110.0)
    assert result.accepted
    assert result.anchor.fit.center_pixel == pytest.approx(110.0, abs=0.05)
    assert result.anchor.line.center_pixel == 90.0


def test_a_solved_transform_moves_the_expected_rows_onto_the_frame(tmp_path):
    """The displayed row follows the correction; its identity does not."""

    session = _shifted_session(tmp_path)
    _solve_two_anchors(session)

    shown = session.lines_for_order(2)
    assert [line.center_pixel for line in shown] == [pytest.approx(110.0, abs=0.05)]
    # The interval travels with the centre, so the row stays one whole line.
    assert shown[0].pixel_from == pytest.approx(105.0, abs=0.05)
    assert shown[0].pixel_to == pytest.approx(115.0, abs=0.05)
    assert shown[0].wavelength_nm == 620.0

    # The curated table itself is untouched: the correction is a view of it.
    assert [line.center_pixel for line in session.lines] == [30.0, 70.0, 90.0]


def test_a_click_on_the_corrected_row_fits_the_peak_that_is_there(tmp_path):
    session = _shifted_session(tmp_path)
    _solve_two_anchors(session)

    result = session.fit_anchor_at(2, 110.0)

    assert result.accepted
    assert result.anchor.fit.center_pixel == pytest.approx(110.0, abs=0.05)
    assert result.anchor.line.wavelength_nm == 620.0


def test_a_click_on_the_base_row_still_matches_and_still_finds_the_peak(tmp_path):
    """The defect itself: 20 px is inside the 30 px match radius, so the click
    was always accepted — and then measured through a window centred 20 px off,
    which held the line's flank and not its apex.  That fit succeeded and
    returned 108, two pixels of pure error into the solve.  The window now
    opens on the corrected pixel instead."""

    session = _shifted_session(tmp_path)
    _solve_two_anchors(session)

    result = session.fit_anchor_at(2, 90.0)

    assert result.accepted
    assert result.anchor.fit.center_pixel == pytest.approx(110.0, abs=0.05)
    # The anchor is the curated row, so the next solve still measures the whole
    # correction rather than re-solving against its own output.
    assert result.anchor.line.center_pixel == 90.0
    assert result.anchor.key == (2, 90.0, 620.0)


def test_the_correction_does_not_fold_into_itself_across_solves(tmp_path):
    """Re-solving with corrected rows on screen must not shrink the shift."""

    session = _shifted_session(tmp_path)
    _solve_two_anchors(session)

    assert session.fit_anchor_at(2, 110.0).accepted
    assert session.transform.dx_px == pytest.approx(20.0, abs=0.05)
    assert session.rms_px < 0.05

    # Removing an anchor and placing it again lands on the same numbers.
    session.remove_anchor((2, 90.0, 620.0))
    assert session.fit_anchor_at(2, 90.0).accepted
    assert session.transform.dx_px == pytest.approx(20.0, abs=0.05)


def test_a_click_beyond_the_match_radius_of_the_corrected_row_is_refused(tmp_path):
    """The radius is tested against the position the operator can see."""

    session = _shifted_session(tmp_path)
    _solve_two_anchors(session)

    # 70 px is 20 from the base row (inside the old radius) and 40 from the
    # corrected one, which is where the row is drawn and where it now counts.
    refused = session.fit_anchor_at(2, 70.0)

    assert not refused.accepted
    assert "not near a known calibration row" in refused.reason


# ---------------------------------------------------------------------------
# The vertical blind spot: a wrong pattern SOURCE.
#
# Every anchor this bench fits is a dispersion measurement — the centroid
# supplies the column and the detector row that goes with it is read out of the
# reference pattern, never off the frame — so a pattern belonging to another era
# passes the fit, the residuals and the save without a mark.  The owner's 2024
# folder calibrated against the packaged 2025 pattern is exactly that, and its
# snapshot saved green.  The sphere frame is the one frame that lights every
# order, so it is the one frame that can answer the question the fit cannot ask.
# ---------------------------------------------------------------------------

_SPHERE_COLUMNS = 640
_SPHERE_ROWS = 260
#: Four order bands, sixty rows apart, sloping the way real traces do.
_SPHERE_ORDER_ROWS = (40.0, 100.0, 160.0, 220.0)
_SPHERE_SLOPE = 0.01
#: Wider than the extraction half width, as a real order band is.
_SPHERE_BAND_SIGMA = 5.0
_SPHERE_LEVEL = 20000.0
_SPHERE_FLOOR = 40.0


def _sphere_traces(shift_rows: float = 0.0) -> np.ndarray:
    columns = np.arange(_SPHERE_COLUMNS, dtype=float)
    return np.column_stack(
        [row + _SPHERE_SLOPE * columns + shift_rows for row in _SPHERE_ORDER_ROWS]
    )


def _sphere_pattern(shift_rows: float = 0.0) -> np.ndarray:
    """The integer pattern file, as one is really written."""

    return np.rint(_sphere_traces(shift_rows)).astype(int)


def _sphere_image(
    shift_rows: float = 0.0, *, dark_orders: tuple[int, ...] = ()
) -> np.ndarray:
    """A continuum frame whose bands sit *shift_rows* off the base pattern."""

    rows = np.arange(_SPHERE_ROWS, dtype=float)[:, None]
    traces = _sphere_traces(shift_rows)
    image = np.full((_SPHERE_ROWS, _SPHERE_COLUMNS), _SPHERE_FLOOR)
    for order_idx in range(traces.shape[1]):
        if order_idx in dark_orders:
            continue
        centers = traces[:, order_idx][None, :]
        image = image + _SPHERE_LEVEL * np.exp(
            -0.5 * ((rows - centers) / _SPHERE_BAND_SIGMA) ** 2
        )
    return image


def test_band_centres_read_the_row_shift_the_fit_cannot_see():
    """A sphere three and a half rows off its pattern says so, per order."""

    reading = band_center_offsets(_sphere_image(3.5), _sphere_pattern())

    assert reading.measured
    assert reading.median_offset_rows == pytest.approx(3.5, abs=0.3)
    assert reading.order_count == len(_SPHERE_ORDER_ROWS)
    for order in reading.orders:
        assert order.offset_rows == pytest.approx(3.5, abs=0.3), order
    assert reading.exceeds()
    assert "rows above the chosen pattern" in reading.summary()
    assert "median of 4 order(s)" in reading.summary()


def test_a_sphere_below_its_pattern_is_reported_as_below():
    """The sign is words, not a convention the reader has to look up."""

    reading = band_center_offsets(_sphere_image(-3.5), _sphere_pattern())

    assert reading.median_offset_rows == pytest.approx(-3.5, abs=0.3)
    assert "rows below the chosen pattern" in reading.summary()


def test_a_sphere_on_its_own_pattern_reads_as_no_shift():
    """The matched era is quiet: near zero, and under the threshold."""

    reading = band_center_offsets(_sphere_image(), _sphere_pattern())

    assert reading.measured
    assert reading.median_offset_rows == pytest.approx(0.0, abs=0.2)
    assert not reading.exceeds()


def test_a_saturated_band_column_is_left_out_of_the_reading():
    """Saturation is excluded the way triage counts it: clusters, not spikes."""

    image = _sphere_image(3.5)
    saturated_columns = range(0, _SPHERE_COLUMNS, 64)
    for column in saturated_columns:
        row = int(round(_sphere_traces(3.5)[column, 1]))
        image[row - 1 : row + 2, column] = 65535.0
    # One lone full-scale pixel, in an order that keeps every column.
    image[int(round(_sphere_traces(3.5)[0, 0])), 0] = 65535.0

    reading = band_center_offsets(image, _sphere_pattern())

    saturated = reading.orders[1]
    assert saturated.offset_rows is None
    assert saturated.columns_excluded == len(list(saturated_columns))
    assert reading.columns_excluded == saturated.columns_excluded
    # The lone spike neither excluded its column nor moved the answer.
    assert reading.orders[0].columns_excluded == 0
    assert reading.orders[0].offset_rows == pytest.approx(3.5, abs=0.3)
    assert reading.median_offset_rows == pytest.approx(3.5, abs=0.3)


def test_an_order_carrying_no_light_is_unmeasured_never_zero():
    """A dark order states its absence and never averages into the median."""

    reading = band_center_offsets(
        _sphere_image(3.5, dark_orders=(2,)), _sphere_pattern()
    )

    dark = reading.orders[2]
    assert dark.offset_rows is None
    assert "no usable light" in dark.reason
    assert reading.order_count == len(_SPHERE_ORDER_ROWS) - 1
    assert reading.median_offset_rows == pytest.approx(3.5, abs=0.3)


def test_a_frame_with_no_light_at_all_says_so_instead_of_zero():
    reading = band_center_offsets(
        np.full((_SPHERE_ROWS, _SPHERE_COLUMNS), 12.0), _sphere_pattern()
    )

    assert not reading.measured
    assert reading.median_offset_rows is None
    assert "no order carried enough light" in reading.reason
    assert not reading.exceeds()


def test_a_pattern_from_another_detector_is_refused_by_its_own_width():
    reading = band_center_offsets(_sphere_image(), _sphere_pattern()[:100])

    assert not reading.measured
    assert "100 column(s) wide" in reading.reason
    assert f"frame is {_SPHERE_COLUMNS}" in reading.reason


def test_a_frame_stack_is_averaged_before_it_is_read():
    """Three exposures of one acquisition are one reading, as they are shot."""

    stack = np.stack([_sphere_image(3.5)] * 3)

    reading = band_center_offsets(stack, _sphere_pattern())

    assert reading.median_offset_rows == pytest.approx(3.5, abs=0.3)


def test_the_reading_window_never_reaches_into_the_neighbouring_order():
    """Traces closer than the default window shrink it rather than mix orders."""

    columns = np.arange(_SPHERE_COLUMNS, dtype=float)
    crowded = np.rint(
        np.column_stack([40 + 0.0 * columns, 58 + 0.0 * columns])
    ).astype(int)
    rows = np.arange(_SPHERE_ROWS, dtype=float)[:, None]
    image = np.full((_SPHERE_ROWS, _SPHERE_COLUMNS), _SPHERE_FLOOR)
    for center in (43.5, 61.5):
        image = image + _SPHERE_LEVEL * np.exp(-0.5 * ((rows - center) / 2.0) ** 2)

    reading = band_center_offsets(image, crowded)

    assert reading.measured
    for order in reading.orders:
        assert order.offset_rows == pytest.approx(3.5, abs=0.3), order


#: The owner's own light, when the working folder travels with the checkout.
#: The packaged 2025 pattern and the 2025 sphere are the matched era, so this
#: reading is the one that has to come out near zero — a guard that cried wolf
#: on the era it was written for would be turned off by the second evening.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_LIVE_SPHERE = _REPO_ROOT / "local" / "20250926_calib" / "sphere-0.1s-x3.sif"
_PACKAGED_2025_PATTERN = (
    _REPO_ROOT
    / "src"
    / "echelle_spectra"
    / "resources"
    / "calibration_files"
    / "pattern_CMOS_20250926.txt"
)


@pytest.mark.skipif(
    not _LIVE_SPHERE.is_file(), reason="local/20250926_calib is not in this checkout"
)
def test_the_2025_sphere_sits_on_the_packaged_2025_pattern():
    """The matched era reads quiet, through the bench's own reader."""

    pattern = np.loadtxt(_PACKAGED_2025_PATTERN, dtype=int)
    frame = FrameLoader(pattern)(_LIVE_SPHERE)

    reading = band_center_offsets(frame.detector_image, pattern)

    assert reading.measured
    assert reading.order_count == pattern.shape[1]
    # Loosely pinned: this is a measurement of real optics, not a fixture.
    assert abs(reading.median_offset_rows) < 1.0
    assert not reading.exceeds()


# ---------------------------------------------------------------------------
# The pattern is not a constant.
#
# It was a constructor argument for as long as the only way to change it was to
# close the bench.  The bench can extract one from its own sphere now, so both
# the state machine and the reader take one mid-session — and the anchors that
# were fitted on the old geometry go, because an anchor's detector row is read
# out of the pattern it was placed against.
# ---------------------------------------------------------------------------


def test_adopting_a_pattern_drops_the_anchors_fitted_on_the_old_geometry(tmp_path):
    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))
    session.accept_frame(_frame(tmp_path))
    session.upsert_anchor(_anchor(_line(0, 32.0, 585.2), 32.4))
    session.upsert_anchor(_anchor(_line(1, 72.0, 588.1), 72.5))
    assert session.transform is not None

    moved = _pattern() + 5.0
    cleared = session.adopt_pattern(moved)

    assert cleared == 2
    assert not session.anchors
    assert session.transform is None
    assert session.rms_px is None
    assert session.residuals == ()
    assert session.alignment_state is AlignmentState.EMPTY
    np.testing.assert_allclose(session.pattern, moved)
    # The frame itself is kept: it is re-extracted by the caller, not lost.
    assert session.frame is not None
    assert session.file_state is FileLoadState.LOADED


def test_a_pattern_with_another_order_count_releases_the_open_frame(tmp_path):
    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))
    session.accept_frame(_frame(tmp_path))
    session.set_selected_order(2)

    session.adopt_pattern(np.full((120, 2), 12.0))

    assert session.frame is None
    assert session.file_state is FileLoadState.WAITING
    assert session.alignment_state is AlignmentState.WAITING_FOR_FRAME
    assert session.selected_order == 1


def test_a_pattern_that_is_not_a_table_is_refused(tmp_path):
    session = CalibrationBenchSession(_pattern(), (_line(0, 32.0, 585.2),))

    with pytest.raises(ValueError):
        session.adopt_pattern(np.arange(10.0))

    assert session.pattern.shape == (120, 3)


def test_the_frame_loader_reads_later_frames_off_the_adopted_pattern():
    loader = FrameLoader(_pattern().astype(int))

    loader.adopt_pattern((_pattern() + 7).astype(int))

    np.testing.assert_array_equal(loader.pattern, (_pattern() + 7).astype(int))
    with pytest.raises(ValueError):
        loader.adopt_pattern(np.arange(5))
