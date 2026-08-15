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
