"""Headless GUI tests for known-line overlay toggles and zoom behavior."""

from __future__ import annotations

import colorsys
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pyqtgraph as pg
import pytest
from PyQt5 import QtWidgets

from echelle_spectra.resources.window_layout import Ui_MainWindow
from echelle_spectra.tools.image_line_overlay import _brighter
from echelle_spectra.tools.line_catalog import LINE_FAMILIES, load_line_table
from echelle_spectra.tools.line_overlay import (
    LINE_OVERLAY_STYLES,
    SPECTRUM_CURVE_COLORS,
    LineOverlayManager,
    choose_label_lines,
    select_overlay_lines,
)

#: The span the 20250926 CMOS wavelength solution reaches, which is the view
#: the operator opens on and the one the detector overlay always asks for.
DETECTOR_SPAN_NM = (401.5, 801.9)

#: Bright Ne I lines the owner found unlabeled on
#: ``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif``.  Every one of them is
#: in the packaged NIST cache and every one is below the broad-view strength
#: floor of 0.08 — 667.83 at 0.05 and 671.70 at 0.007 are genuinely weak in the
#: database and unmissable in a real neon lamp — and every one is an OK-marked
#: row of the curated 20240305 wavelength table.
UNLABELED_FIELD_LINES = (594.4834, 609.6163, 653.2882, 667.8276, 671.7043)


@pytest.fixture(scope="module")
def qt_app():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def test_gui_exposes_separate_overlay_toggles_off_by_default(qt_app):
    window = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(window)
    assert list(ui.line_overlay_checks) == ["balmer", "fulcher", "thar", "ne", "hg"]
    assert all(not checkbox.isChecked() for checkbox in ui.line_overlay_checks.values())
    # The two toggles that sit beside them start off as well.
    assert not ui.order_trace_check.isChecked()
    assert not ui.cursor_link_check.isChecked()
    window.close()


def test_manager_toggle_zoom_and_disabled_state(qt_app):
    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    plot.setXRange(600.0, 616.0, padding=0)
    manager = LineOverlayManager(max_labels=8)
    manager.register_plot("spectrum", plot)

    assert manager.rendered_lines("spectrum", "fulcher") == ()
    manager.set_family_visible("fulcher", True)
    qt_app.processEvents()
    broad = manager.rendered_lines("spectrum", "fulcher")
    assert 1 < len(broad) <= 8
    assert all(600.0 <= line.wavelength_nm <= 616.0 for line in broad)

    plot.setXRange(400.0, 680.0, padding=0)
    qt_app.processEvents()
    manager.refresh()
    assert len(manager.rendered_lines("spectrum", "fulcher")) <= 4

    plot.setXRange(601.7, 603.3, padding=0)
    qt_app.processEvents()
    manager.refresh()
    assert [line.label for line in manager.rendered_lines("spectrum", "fulcher")] == [
        "Q1(0-0)",
        "Q2(0-0)",
        "Q3(0-0)",
    ]

    manager.set_family_visible("fulcher", False)
    assert manager.rendered_lines("spectrum", "fulcher") == ()
    assert not manager.is_family_visible("fulcher")
    widget.close()


def test_atomic_zoom_reveals_weaker_lines_without_overcrowding():
    """The strength threshold still rules every row the curated table is silent on."""

    broad = select_overlay_lines("thar", 578.0, 640.0, max_labels=10)
    narrow = select_overlay_lines("thar", 618.0, 619.0, max_labels=20)
    assert broad
    assert narrow
    uncurated = [line for line in broad if not line.curated]
    assert uncurated
    assert min(line.relative_intensity for line in uncurated) >= 0.08
    assert min(
        line.relative_intensity
        for line in narrow
        if not line.curated and line.relative_intensity is not None
    ) < 0.08


def test_a_vetted_line_never_loses_a_selection_contest_to_a_database_number():
    """The union rule, on the five lines the owner found unmarked in the field.

    Each is below the broad-view strength floor, so the threshold alone drops
    every one of them; each is an OK-marked row of the curated table, so the
    union keeps every one of them.  Curation buys standing, not a free pass:
    the rows the threshold selects on its own are all still there.
    """

    selected = select_overlay_lines("ne", *DETECTOR_SPAN_NM, max_labels=200)
    by_wavelength = {round(line.wavelength_nm, 4): line for line in selected}

    for wavelength in UNLABELED_FIELD_LINES:
        found = [
            line for line in selected if abs(line.wavelength_nm - wavelength) < 0.01
        ]
        assert found, f"{wavelength} nm is still unselected in the full-span view"
        assert found[0].curated
        assert found[0].relative_intensity < 0.08, (
            f"{wavelength} nm would have been selected without its curation, "
            "so it pins nothing"
        )

    strong = [
        line
        for line in load_line_table("ne")
        if DETECTOR_SPAN_NM[0] <= line.wavelength_nm <= DETECTOR_SPAN_NM[1]
        and not line.curated
        and line.relative_intensity is not None
        and line.relative_intensity >= 0.08
    ]
    assert strong
    for line in strong:
        assert round(line.wavelength_nm, 4) in by_wavelength, line.label


def test_a_broad_view_of_a_neon_lamp_is_not_a_view_of_neon_ions():
    """The whole detector, and the marks land where the lamp puts light.

    Pinned on the 401.5--801.9 nm span of the 20250926 CMOS solution, which is
    what the 2-D detector overlay asks for.  Before the cache was widened this
    returned 33 rows spanning 580--638 nm, 16 of them Ne II.

    Adding the curated union must not undo that: the curated table carries no
    OK-marked Ne II row at all, so the ion count is the threshold's alone and
    stays where it was.  Measured on the owner's bright frame, exactly one Ne
    box lands on dark detector for an ionic reason, and it is this one.
    """

    selected = select_overlay_lines("ne", *DETECTOR_SPAN_NM, max_labels=200)
    assert selected
    ions = [line for line in selected if line.species == "Ne II"]
    assert len(ions) == 1, [line.label for line in ions]
    assert not any(line.curated for line in ions)

    # The strong red neon lines the old cache stopped short of, including the
    # brightest blob on the owner's frame at 640.2248 nm.
    for wavelength in (640.2248, 650.6528, 692.9467, 703.2413, 724.5167):
        assert any(
            abs(line.wavelength_nm - wavelength) < 0.01 for line in selected
        ), f"{wavelength} nm missing from the broad Ne selection"


def _synthetic_spectrum(peaks, *, low=580.0, high=680.0, samples=20001):
    """A flat trace with a Gaussian at each ``(wavelength, height)`` given."""

    wavelengths = np.linspace(low, high, samples)
    values = np.full(samples, 10.0)
    for wavelength, height in peaks:
        values += height * np.exp(-0.5 * ((wavelengths - wavelength) / 0.02) ** 2)
    return wavelengths, values


def test_a_crowded_view_labels_the_lines_the_measured_spectrum_is_loudest_at():
    """Placement follows the data on screen, not the database.

    Four candidates and room for two.  The two the cached strengths would pick
    are deliberately the two the synthetic frame shows nothing at, so a ranking
    that ignored the measurement could not accidentally pass.
    """

    lines = select_overlay_lines("ne", 580.0, 680.0, max_labels=400)
    ranked = sorted(lines, key=lambda line: -(line.relative_intensity or 0.0))
    loud_one, loud_two = ranked[-1], ranked[-2]
    quiet_one, quiet_two = ranked[0], ranked[1]
    candidates = (loud_one, loud_two, quiet_one, quiet_two)

    spectrum = _synthetic_spectrum(
        ((loud_one.wavelength_nm, 9000.0), (loud_two.wavelength_nm, 6000.0))
    )

    chosen = choose_label_lines(candidates, 2, spectrum=spectrum)
    assert {line.wavelength_nm for line in chosen} == {
        loud_one.wavelength_nm,
        loud_two.wavelength_nm,
    }
    assert [line.wavelength_nm for line in chosen] == sorted(
        line.wavelength_nm for line in chosen
    )

    # Without a measurement there is nothing else to go on, so the cached
    # strengths rank them — which is what used to happen always.
    blind = choose_label_lines(candidates, 2)
    assert {line.wavelength_nm for line in blind} == {
        quiet_one.wavelength_nm,
        quiet_two.wavelength_nm,
    }


def test_every_selected_line_keeps_its_stick_when_the_labels_run_out(qt_app):
    """The cap rations text, never markers — a stick is how zoom is invited."""

    widget = pg.PlotWidget()
    plot = widget.getPlotItem()
    plot.setXRange(*DETECTOR_SPAN_NM, padding=0)
    manager = LineOverlayManager(max_labels=6)
    manager.register_plot("calibrated", plot, labels=True)
    loudest = 671.7043
    manager.set_measured_spectrum(
        "calibrated", *_synthetic_spectrum(((loudest, 9000.0),), low=400.0, high=810.0)
    )
    manager.set_family_visible("ne", True)
    qt_app.processEvents()
    manager.refresh()

    drawn = manager.rendered_lines("calibrated", "ne")
    labelled = manager.labelled_lines("calibrated", "ne")
    assert len(labelled) < len(drawn)
    assert len(labelled) <= 6
    assert set(labelled) <= set(drawn)
    # Every field line keeps a stick, and the one the frame is loud at is named.
    for wavelength in UNLABELED_FIELD_LINES:
        assert any(abs(line.wavelength_nm - wavelength) < 0.01 for line in drawn)
    assert any(abs(line.wavelength_nm - loudest) < 0.01 for line in labelled)
    # One item per drawn line, and only the labelled ones carry text.
    items = manager._items["calibrated"]["ne"]
    assert len(items) == len(drawn)
    # pyqtgraph only gives an InfiniteLine a ``label`` attribute when it was
    # built with text, so its absence is exactly "this stick carries none".
    assert len(
        [item for item in items if getattr(item, "label", None) is not None]
    ) == len(labelled)
    widget.close()


def test_no_family_colour_hides_in_a_spectrum_curve():
    """Overlay hues must contrast with both curves they are drawn over.

    The Ne family shipped a green stick over the green counts curve, which is
    a marker that annotates by disappearing.  Red and yellow-green bracket the
    warm arc, so the palette lives in the cool one; this pins that it stays
    there and that no two families collapse into each other either.
    """

    def hue(color: str) -> float:
        red, green, blue, _alpha = pg.mkColor(color).getRgb()
        return colorsys.rgb_to_hls(red / 255, green / 255, blue / 255)[0] * 360.0

    def separation(first: float, second: float) -> float:
        gap = abs(first - second) % 360.0
        return min(gap, 360.0 - gap)

    curves = {name: hue(color) for name, color in SPECTRUM_CURVE_COLORS.items()}
    families = {name: hue(style.color) for name, style in LINE_OVERLAY_STYLES.items()}
    assert set(families) == set(LINE_FAMILIES)

    for family, family_hue in families.items():
        for curve, curve_hue in curves.items():
            assert separation(family_hue, curve_hue) >= 40.0, (
                f"{family} at {family_hue:.0f} deg is too close to the "
                f"{curve} curve at {curve_hue:.0f} deg"
            )

    names = sorted(families)
    for index, first in enumerate(names):
        for second in names[index + 1 :]:
            assert separation(families[first], families[second]) >= 25.0, (
                f"{first} and {second} are not far enough apart in hue"
            )


def test_the_detector_boxes_use_the_same_family_hue_as_the_sticks():
    """One family, one hue, both views — the boxes only lift it toward white."""

    for family, style in LINE_OVERLAY_STYLES.items():
        lifted = _brighter(style.color)
        base = pg.mkColor(style.color).getRgb()[:3]
        lit = pg.mkColor(lifted).getRgb()[:3]
        assert all(after >= before for before, after in zip(base, lit)), family
        assert lit != base, family
        assert (
            colorsys.rgb_to_hls(*(value / 255 for value in lit))[2] > 0.15
        ), f"{family} washes out to grey on the detector image"
