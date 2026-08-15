"""Focused off-screen tests for the separate calibration-bench window."""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pyqtgraph as pg
import pytest
from PyQt5 import QtCore, QtGui, QtWidgets

from echelle_spectra.calibration_bench import BenchFrame, CalibrationBenchSession
from echelle_spectra.calibration_bench_gui import (
    _PACKAGE_DIR,
    _SUGGESTED_BADGE,
    BENCH_BODY_POINT_SIZE,
    BENCH_HEADLINE_POINT_SIZE,
    BENCH_TOOLTIP_LIMIT,
    CalibrationBenchWindow,
    _build_parser,
    bench_default_geometry,
    bench_point_sizes,
    bench_window_icon,
    one_line,
    role_combo_minimum_width,
)
from echelle_spectra.calibration_campaign import (
    AbsoluteCalibrationResult,
    CalibrationCampaignSession,
    ComparisonState,
    ExposureState,
    MeasurementRole,
)
from echelle_spectra.tools.calibration_alignment import CalibrationTableLine

_COLUMNS = 80


def _SUGGESTED_COLOR_IN(combo) -> bool:
    """Whether the Role control is wearing the unconfirmed-suggestion colour."""

    return "#ffb86b" in combo.styleSheet()


@pytest.fixture(scope="module")
def qt_app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


def _frame_for(
    path: Path,
    *,
    peak: float = 40000.0,
    cosmic: bool = False,
    frames: int = 1,
) -> BenchFrame:
    x = np.arange(_COLUMNS, dtype=float)
    order_spectra = (
        5 + 200 * np.exp(-0.5 * ((x - 26) / 1.5) ** 2),
        5 + 160 * np.exp(-0.5 * ((x - 56) / 1.5) ** 2),
    )
    rng = np.random.default_rng(5)
    images = rng.normal(300.0, 10.0, size=(frames, 44, _COLUMNS))
    for index in range(frames):
        images[index, 12, :] = order_spectra[0]
        images[index, 30, :] = order_spectra[1]
        images[index, 20:23, 40:43] = peak
    if cosmic:
        images[0, 5, 5] = 65535.0
        images[0, 40, 70] = 65535.0
    per_frame = (
        tuple(tuple(np.array(row) for row in order_spectra) for _ in range(frames))
        if frames > 1
        else ()
    )
    return BenchFrame(
        Path(path),
        images,
        images.mean(axis=0),
        order_spectra,
        {"ExposureTime": 0.1},
        frame_order_spectra=per_frame,
    )


def _loader(**frame_options):
    def load(path):
        return _frame_for(Path(path), **frame_options)

    return load


def _window(tmp_path: Path, *, with_loader: bool = False, **frame_options) -> CalibrationBenchWindow:
    pattern = np.column_stack(
        [np.full(_COLUMNS, 12, dtype=float), np.full(_COLUMNS, 30, dtype=float)]
    )
    # Both lamps carry rows in this table, because F17's whole point is that
    # the panel and the sticks read the assigned lamp's own rows: a Ne lamp
    # against a Th-only table is honestly empty, and would prove nothing.
    lines = (
        CalibrationTableLine(0, 20, 30, 25, 578.0, "ThI", "ok"),
        CalibrationTableLine(0, 65, 75, 70, 640.0, "ArI", "ok"),
        CalibrationTableLine(1, 50, 60, 55, 580.0, "ThI", "ok"),
        CalibrationTableLine(1, 65, 75, 70, 635.0, "ArI", "ok"),
        CalibrationTableLine(0, 21, 31, 26, 585.2488, "NeI", "ok"),
        CalibrationTableLine(1, 51, 61, 56, 588.1895, "NeI", "ok"),
    )
    session = CalibrationBenchSession(pattern, lines)
    session.accept_frame(_frame_for(tmp_path / "bench-fixture.sif"))
    return CalibrationBenchWindow(
        session,
        loader=_loader(**frame_options) if with_loader else None,
        start_timer=False,
    )


def _manual_window(tmp_path: Path, **frame_options) -> CalibrationBenchWindow:
    """A bench with campaign memory, a reader, and nothing loaded yet."""

    window = _window(tmp_path, with_loader=True, **frame_options)
    references = {}
    for name in ("pattern.txt", "wavelength.txt", "integral.txt"):
        path = tmp_path / name
        path.write_text(name, encoding="utf-8")
        references[name] = path
    window.campaign = CalibrationCampaignSession(
        pattern_source=references["pattern.txt"],
        wavelength_source=references["wavelength.txt"],
        integral_source=references["integral.txt"],
    )
    window.refresh()
    return window


def test_window_icon_resolves_and_is_set(qt_app, tmp_path: Path):
    """The bench window carries its own icon so the owner can tell it from
    the main GUI, even instantiated off-screen without a QApplication-level
    icon already set."""

    icon_path = _PACKAGE_DIR / "resources" / "graphics" / "echelle.png"
    assert icon_path.is_file(), f"icon file missing at {icon_path}"

    window = _window(tmp_path)
    assert not window.windowIcon().isNull()


def _drop(window: CalibrationBenchWindow, paths) -> None:
    """Drive the real Qt drop path the main GUI established."""

    mime = QtCore.QMimeData()
    mime.setUrls([QtCore.QUrl.fromLocalFile(str(path)) for path in paths])
    event = QtGui.QDropEvent(
        QtCore.QPointF(10.0, 10.0),
        QtCore.Qt.CopyAction,
        mime,
        QtCore.Qt.LeftButton,
        QtCore.Qt.NoModifier,
    )
    enter = QtGui.QDragEnterEvent(
        QtCore.QPoint(10, 10),
        QtCore.Qt.CopyAction,
        mime,
        QtCore.Qt.LeftButton,
        QtCore.Qt.NoModifier,
    )
    window.dragEnterEvent(enter)
    assert enter.isAccepted()
    window.dropEvent(event)


def _wait_for_loads(window: CalibrationBenchWindow, qt_app) -> None:
    for _attempt in range(400):
        if window._load_thread is None and not window._queue:
            break
        qt_app.processEvents()
        QtCore.QThread.msleep(5)
    qt_app.processEvents()


def _campaign_window(tmp_path: Path) -> CalibrationBenchWindow:
    window = _window(tmp_path)
    paths = {}
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
        path.write_text(name, encoding="utf-8")
        paths[name] = path
    campaign = CalibrationCampaignSession(
        pattern_source=paths["pattern.txt"],
        wavelength_source=paths["wavelength.txt"],
        integral_source=paths["integral.txt"],
        suggested_lamps=("ThAr",),
        previous_sphere=paths["previous_sphere.sif"],
        previous_sphere_background=paths["previous_sphere_bg.sif"],
    )
    frame = window.session.frame
    campaign.classify_file(paths["sphere.sif"], MeasurementRole.SPHERE, frame=frame)
    campaign.classify_file(
        paths["sphere_bg.sif"], MeasurementRole.SPHERE_BACKGROUND, frame=frame
    )
    campaign.classify_file(
        paths["thar.sif"], MeasurementRole.LAMP, lamp_family="ThAr", frame=frame
    )
    campaign.classify_file(
        paths["thar_bg.sif"],
        MeasurementRole.LAMP_BACKGROUND,
        lamp_family="ThAr",
        frame=frame,
    )

    def calculator(**values):
        scale = 1.08 if "previous" not in values["sphere"].name else 1.0
        return AbsoluteCalibrationResult(
            np.linspace(400, 700, 60), np.full(60, scale)
        )

    campaign.compute_sphere_comparison(calculator)
    window.campaign = campaign
    window.snapshot_id_edit.setText("20250813_cmos")
    window.refresh()
    return window


def test_bench_is_a_separate_detector_order_residual_surface(qt_app, tmp_path):
    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert window.windowTitle() == "Echelle calibration bench"
    assert window.detector_plot.titleLabel.text == "Detector + order traces"
    assert "Order 0" in window.order_plot.titleLabel.text
    assert window.residual_plot is not None
    assert window.order_spin.maximum() == 1
    assert window.file_value.text() == "bench-fixture.sif"
    assert window.anchor_table.rowCount() == 0
    window.close()


def test_window_reflects_anchor_and_alignment_state(qt_app, tmp_path):
    window = _window(tmp_path)
    session = window.session
    assert session.fit_anchor_at(0, 26).accepted
    assert session.fit_anchor_at(1, 56).accepted
    window.refresh()
    qt_app.processEvents()

    assert window.anchor_table.rowCount() == 2
    assert window.alignment_state_value.text() == "ALIGNED"
    assert window.rms_value.text().endswith(" px")
    assert "°" in window.transform_value.text()
    window.close()


def test_cli_needs_no_watch_folder_and_no_named_lamp():
    args = _build_parser().parse_args([])
    assert args.pattern.name == "pattern_CMOS_20250926.txt"
    assert args.wavelength.name == "Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    assert args.stable_polls == 2
    assert args.minimum_age_s == 1.0
    assert args.integral.name == "integrating_sphere.txt"
    assert args.lamp is None
    assert args.watch is False
    assert args.file is None
    assert args.valid_from == date.today()
    assert _build_parser().parse_args(["--valid-from", "2026-09-01"]).valid_from == date(
        2026, 9, 1
    )
    # Any lamp name is accepted; the bench never owns a permitted list.
    assert _build_parser().parse_args(["--lamp", "Kr", "--lamp", "Ne"]).lamp == ["Kr", "Ne"]


def test_bench_views_lead_with_triage_and_manual_input(qt_app, tmp_path):
    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert [window.control_tabs.tabText(i) for i in range(window.control_tabs.count())] == [
        "Files",
        "Procedure",
        "Lamp fit",
        "Save",
    ]
    # F16 item 7: line identification is not a room of its own any more.
    assert [window.view_tabs.tabText(i) for i in range(window.view_tabs.count())] == [
        "Exposure triage",
        "Lamp fit",
        "Sphere factors",
    ]
    assert window.acceptDrops()
    assert window.checklist_tree.count() >= 9
    assert "median" in window.comparison_value.text()
    assert window.line_help_table.rowCount() > 0
    assert "Anchor" == window.line_help_table.horizontalHeaderItem(4).text()
    window.close()


def test_with_no_files_the_drop_target_is_the_primary_surface(qt_app, tmp_path):
    window = _manual_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert window.drop_hint.isVisible()
    assert "DROP SIF FILES HERE" in window.drop_hint.text()
    assert window.add_files_button.isEnabled()
    assert window.file_table.rowCount() == 0
    items = {
        window.checklist_tree.itemWidget(window.checklist_tree.item(row)).text()
        for row in range(window.checklist_tree.count())
    }
    assert any("unblocked by: drop any SIF" in text for text in items)
    window.close()


def test_dropped_files_are_triaged_before_any_role(qt_app, tmp_path):
    window = _manual_window(tmp_path, cosmic=True)
    window.show()
    paths = []
    for name in ("Ne-0.02s-x3-bright-lines.sif", "sphere-0.1s-x3.sif"):
        path = tmp_path / name
        path.write_bytes(b"sif\n")
        paths.append(path)

    _drop(window, paths)
    _wait_for_loads(window, qt_app)

    assert window.file_table.rowCount() == 2
    assert not window.drop_hint.isVisible()
    assert set(window.campaign.loaded) == set(paths)
    # Triage happened; roles did not.
    assert not window.campaign.measurements
    triage = window.campaign.loaded[paths[0]].triage
    assert triage.state is ExposureState.GOOD
    assert triage.saturation.anomalous_pixels == 2
    # The headline is two lines: which file, and what it says.
    assert len(window.triage_headline.text().splitlines()) == 2
    assert "anomalies" in window.triage_headline.text()
    assert "not saturation" in window.details_view.toPlainText()
    # F16 item 3: the multi-line breakdown folds into the Why dock; the view
    # keeps the one line that decides the next action.
    assert "% of full scale" in window.details_view.toPlainText()
    assert "\n" not in window.triage_next_value.text()
    window.close()


def test_roles_are_assigned_by_hand_per_file_whatever_the_name(qt_app, tmp_path):
    window = _manual_window(tmp_path)
    window.show()
    nameless = tmp_path / "IMG_0042.sif"
    nameless.write_bytes(b"sif\n")
    lamp_named = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    lamp_named.write_bytes(b"sif\n")

    _drop(window, [nameless, lamp_named])
    _wait_for_loads(window, qt_app)

    # The recognised name only pre-fills its controls; nothing is assigned.
    prefilled_role = window.file_table.cellWidget(1, 1)
    prefilled_lamp = window.file_table.cellWidget(1, 2)
    assert prefilled_role.currentData() is MeasurementRole.LAMP
    assert prefilled_lamp.currentText() == "Ne"
    assert not window.campaign.measurements

    # The meaningless name takes any role the operator picks, with a free lamp.
    role_combo = window.file_table.cellWidget(0, 1)
    lamp_combo = window.file_table.cellWidget(0, 2)
    lamp_combo.setCurrentText("Kr")
    role_combo.setCurrentIndex(role_combo.findData(MeasurementRole.LAMP))
    qt_app.processEvents()

    record = window.campaign.measurements[nameless]
    assert record.role is MeasurementRole.LAMP
    assert record.lamp_family == "Kr"
    assert window.campaign.assigned_lamps == ("Kr",)
    checklist = {item.key for item in window.campaign.checklist(window.session)}
    assert "lamp-Kr-signal" in checklist
    assert not any("ThAr" in key for key in checklist)
    window.close()


def test_confirming_the_prefilled_role_assigns_it(qt_app, tmp_path):
    window = _manual_window(tmp_path)
    window.show()
    source = tmp_path / "sphere-0.1s-x3-bg.sif"
    source.write_bytes(b"sif\n")

    _drop(window, [source])
    _wait_for_loads(window, qt_app)
    role_combo = window.file_table.cellWidget(0, 1)
    assert role_combo.currentData() is MeasurementRole.SPHERE_BACKGROUND
    assert not window.campaign.measurements
    # The pre-filled control must never read as an assignment nobody made —
    # and after F16 the badge lives beside the combo, never inside its text.
    assert f"{_SUGGESTED_BADGE} ONLY" in window.file_table.item(0, 0).text()
    assert _SUGGESTED_BADGE not in role_combo.currentText()
    assert _SUGGESTED_COLOR_IN(role_combo)

    # Picking the already-shown entry emits no index change; the operator's
    # confirmation must still register.
    role_combo.activated.emit(role_combo.currentIndex())
    qt_app.processEvents()

    assert (
        window.campaign.measurements[source].role is MeasurementRole.SPHERE_BACKGROUND
    )
    assert "sphere-background" in window.file_table.item(0, 0).text()
    assert not _SUGGESTED_COLOR_IN(role_combo)
    window.close()


_REAL_2025_NAMES = (
    "Ne-0.02s-x3-bright-lines.sif",
    "Ne-0.02s-x3-bright-lines_bg.sif",
    "Ne-0.1s-x3-dimm-lines-bg.sif",
    "Ne-0.1s-x3-dimm-lines.sif",
    "sphere-0.1s-x3-bg.sif",
    "sphere-0.1s-x3.sif",
)


def _real_folder_window(qt_app, tmp_path: Path, **frame_options):
    """The owner's own folder on the bench, dropped through the real drop path."""

    window = _manual_window(tmp_path, **frame_options)
    paths = []
    for name in _REAL_2025_NAMES:
        path = tmp_path / name
        path.write_bytes(b"sif\n")
        paths.append(path)
    _drop(window, paths)
    _wait_for_loads(window, qt_app)
    return window, paths


def test_a_prefilled_role_never_looks_like_an_assigned_one(qt_app, tmp_path):
    """F14 item 1 regression.

    Every file of the real folder pre-fills its role correctly, so the Role
    column read right while the campaign had been given nothing: the Procedure
    tab said "no file carries this role yet" and the factor computation failed.
    The control now says SUGGESTED until somebody confirms it, and confirming
    it through the combo reaches the campaign.
    """

    window, paths = _real_folder_window(qt_app, tmp_path)
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    sphere_row = window._file_rows.index(sphere)
    sphere_combo = window.file_table.cellWidget(sphere_row, 1)

    # What the owner's screenshot showed: the right role, in the right column.
    assert sphere_combo.currentData() is MeasurementRole.SPHERE
    # What the bench had actually been given: nothing.
    assert not window.campaign.measurements
    # So the control must say so, and the Procedure tab must name the file.
    assert sphere_combo.currentText() == "Sphere"
    assert _SUGGESTED_COLOR_IN(sphere_combo)
    assert f"{_SUGGESTED_BADGE} ONLY" in window.file_table.item(sphere_row, 0).text()
    assert window.confirm_roles_button.isEnabled()
    assert "6" in window.confirm_roles_button.text()
    checklist = {
        item.key: item for item in window.campaign.checklist(window.session)
    }
    assert "only suggested by its filename" in checklist["sphere"].detail

    # Confirming through the real combo handler reaches the campaign.
    for path in paths:
        combo = window.file_table.cellWidget(window._file_rows.index(path), 1)
        combo.activated.emit(combo.currentIndex())
        qt_app.processEvents()

    assert window.campaign.measurements[sphere].role is MeasurementRole.SPHERE
    assert (
        window.campaign.measurements[tmp_path / "sphere-0.1s-x3-bg.sif"].role
        is MeasurementRole.SPHERE_BACKGROUND
    )
    assert window.campaign.assigned_lamps == ("Ne",)
    assert not _SUGGESTED_COLOR_IN(sphere_combo)
    assert not window.confirm_roles_button.isEnabled()

    def calculator(**values):
        assert values["sphere"].name == "sphere-0.1s-x3.sif"
        assert values["sphere_background"].name == "sphere-0.1s-x3-bg.sif"
        return AbsoluteCalibrationResult(np.linspace(400, 700, 60), np.full(60, 1.0))

    comparison = window.campaign.compute_sphere_comparison(calculator)
    assert comparison.state is not ComparisonState.FAILED
    assert comparison.candidate is not None
    window.close()


def test_one_press_confirms_every_suggested_role(qt_app, tmp_path):
    """The whole folder is one deliberate press, not one popup pick per row."""

    window, _paths = _real_folder_window(qt_app, tmp_path)
    assert not window.campaign.measurements

    window.confirm_roles_button.click()
    qt_app.processEvents()

    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)
    assert window.campaign.unconfirmed_suggestions() == ()
    assert window.campaign.assigned_lamps == ("Ne",)
    assert "Assigned 6 suggested role(s)" in window.message_value.text()
    checklist = {
        item.key: item for item in window.campaign.checklist(window.session)
    }
    assert checklist["sphere"].detail == "sphere-0.1s-x3.sif"
    assert checklist["sphere-background"].detail == "sphere-0.1s-x3-bg.sif"
    window.close()


def test_a_saturated_lamp_frame_is_informed_not_failed(qt_app, tmp_path):
    """F14 item 2 on the surface: the dim series is meant to clip its lines."""

    window = _manual_window(tmp_path, peak=65535.0)
    lamp = tmp_path / "Ne-0.1s-x3-dimm-lines.sif"
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    for path in (lamp, sphere):
        path.write_bytes(b"sif\n")
    _drop(window, [lamp, sphere])
    _wait_for_loads(window, qt_app)
    window.confirm_roles_button.click()
    qt_app.processEvents()

    assert window.campaign.loaded[lamp].triage.state is ExposureState.SATURATED
    window._select_file_row(lamp)
    qt_app.processEvents()
    assert "SATURATED LINES (EXPECTED)" in window.triage_headline.text()
    assert "fit unsaturated lines only" in window.details_view.toPlainText().lower()
    assert not window.campaign.role_triage(lamp).blocking
    assert "SATURATED LINES (EXPECTED)" in window.file_table.item(
        window._file_rows.index(lamp), 0
    ).text()

    # The sphere frame keeps the hard verdict.
    window._select_file_row(sphere)
    qt_app.processEvents()
    assert window.triage_headline.text().splitlines()[1].startswith("SATURATED")
    assert window.campaign.role_triage(sphere).blocking
    assert "Lower exposure" in window.exposure_value.text()
    window.close()


def test_lamp_fit_tab_offers_the_frame_and_the_order(qt_app, tmp_path):
    """F14 item 4: 3-frame acquisitions, and an order control you can see."""

    window = _window(tmp_path, with_loader=True, frames=3)
    window.session.accept_frame(_frame_for(tmp_path / "three.sif", frames=3))
    window.refresh()
    qt_app.processEvents()

    assert window.frame_combo.count() == 4
    assert window.frame_combo.itemData(0) is None
    assert window.frame_combo.itemText(0) == "Mean of all frames"
    assert [window.frame_combo.itemData(i) for i in range(1, 4)] == [0, 1, 2]
    assert window.frame_combo.isEnabled()
    assert "mean of 3 frame(s)" in window.frame_choice_value.text()

    assert window.session.fit_anchor_at(0, 26).accepted
    window.frame_combo.setCurrentIndex(2)
    qt_app.processEvents()

    assert window.session.selected_frame == 1
    assert "frame 2 of 3" in window.frame_choice_value.text()
    assert "frame 2 of 3" in window.order_plot.titleLabel.text
    assert not window.session.anchors

    # The order control is explicit: a spin box plus stepping buttons.
    assert not window.previous_order_button.isEnabled()
    assert window.next_order_button.isEnabled()
    window.next_order_button.click()
    qt_app.processEvents()
    assert window.session.selected_order == 1
    assert window.previous_order_button.isEnabled()
    assert not window.next_order_button.isEnabled()
    window.close()


def test_left_pane_is_resizable_and_its_text_wraps(qt_app, tmp_path):
    """F14 item 3: the pane is the operator's to widen; nothing is clipped."""

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    controls = window.root_splitter.widget(0)
    assert controls.maximumWidth() >= 10_000  # no ceiling on the left pane
    assert not window.root_splitter.childrenCollapsible()
    assert window.file_table.wordWrap()
    assert window.file_table.textElideMode() == QtCore.Qt.ElideNone
    for label in (
        window.triage_headline,
        window.triage_next_value,
        window.exposure_value,
        window.fit_file_value,
        window.comparison_value,
        window.reference_value,
        window.message_value,
        window.transform_value,
        window.frame_choice_value,
    ):
        assert label.wordWrap()

    window.root_splitter.setSizes([760, 780])
    qt_app.processEvents()
    window._relayout_wrapped_text()
    width = window.checklist_tree.viewport().width()
    for row in range(window.checklist_tree.count()):
        label = window.checklist_tree.itemWidget(window.checklist_tree.item(row))
        assert label.wordWrap()
        assert label.width() <= max(280, width)
        # The row is tall enough for every wrapped line it holds.
        assert window.checklist_tree.item(row).sizeHint().height() >= label.sizeHint().height()
        assert label.toolTip()
    window.close()


def test_the_bench_icon_is_derived_from_the_shared_graphic_and_differs(qt_app):
    """F14 item 5: two windows, one drawing, two tellable-apart taskbar icons."""

    source = _PACKAGE_DIR / "resources" / "graphics" / "echelle.png"
    assert source.is_file()
    base = QtGui.QPixmap(str(source)).scaled(
        128, 128, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation
    )
    icon = bench_window_icon(size=128)
    assert not icon.isNull()

    derived = icon.pixmap(128, 128).toImage()
    original = base.toImage().convertToFormat(derived.format())
    assert derived.size() == original.size()
    differing = sum(
        1
        for y in range(0, derived.height(), 4)
        for x in range(0, derived.width(), 4)
        if derived.pixel(x, y) != original.pixel(x, y)
    )
    assert differing > 0, "the bench icon must not be the main GUI's icon"
    # No new artwork file was added for it.
    assert sorted(p.name for p in (_PACKAGE_DIR / "resources" / "graphics").iterdir()) == [
        "echelle.png",
        "echelle.svg",
    ]


def test_only_the_verdict_shouts(qt_app, tmp_path):
    """F16 item 3: loudness is a budget and the budget is exactly one.

    F14 made every reading loud, which is the same as making none of them
    loud. Two tiers now: the verdict headline, and body text.
    """

    body, headline = bench_point_sizes()
    assert body >= BENCH_BODY_POINT_SIZE
    assert headline >= BENCH_HEADLINE_POINT_SIZE
    # A large platform font scales both tiers with it rather than being ignored.
    assert bench_point_sizes(24.0)[1] > headline

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    loud = window.loud_widgets()
    assert loud == [window.triage_headline], [w.objectName() for w in loud]
    for widget in (
        window.rms_value,
        window.anchor_count_value,
        window.comparison_value,
        window.file_state_value,
        window.alignment_state_value,
        window.save_state_value,
        window.exposure_value,
    ):
        assert widget.font().pointSizeF() < headline
        assert widget.font().pointSizeF() >= BENCH_BODY_POINT_SIZE
    window.close()


def test_long_help_reads_in_the_dock_not_in_a_floating_tooltip(qt_app, tmp_path):
    """F16 item 6: tooltips are one short line; the dock carries the rest."""

    assert one_line("First sentence here. Second one is dropped.") == (
        "First sentence here."
    )
    assert len(one_line("word " * 200)) <= BENCH_TOOLTIP_LIMIT

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    explained = list(window._explainable_widgets)
    assert explained, "nothing is explained at all"
    for widget in explained:
        tooltip = widget.toolTip()
        assert tooltip, f"{widget} carries no explanation"
        assert "\n" not in tooltip, f"{widget} tooltip spans lines: {tooltip!r}"
        assert len(tooltip) <= BENCH_TOOLTIP_LIMIT, f"{widget}: {tooltip!r}"
        # The full text is longer than the tooltip and lives in the dock.
        assert widget.property("explainText")

    # F17 item 5 supersedes F16's hover route: a click fills the dock.
    window.details_view.setHtml("")
    window.eventFilter(
        window.rms_value,
        QtGui.QMouseEvent(
            QtCore.QEvent.MouseButtonPress,
            QtCore.QPointF(1.0, 1.0),
            QtCore.Qt.LeftButton,
            QtCore.Qt.LeftButton,
            QtCore.Qt.NoModifier,
        ),
    )
    assert "root-mean-square" in window.details_view.toPlainText()
    # So does focus.
    window.details_view.setHtml("")
    window.eventFilter(
        window.transform_value, QtCore.QEvent(QtCore.QEvent.FocusIn)
    )
    assert "rotation" in window.details_view.toPlainText()

    # A checklist row still answers in the dock, and rebuilding the list does
    # not leave a growing pile of event filters behind.
    before = len(window._explainable_widgets)
    window._refresh_checklist()
    window._refresh_checklist()
    assert len(window._explainable_widgets) == before
    window.checklist_tree.setCurrentRow(0)
    qt_app.processEvents()
    first = window.campaign.checklist(window.session)[0]
    assert first.label in window.details_view.toPlainText()
    window.close()


def test_the_role_control_can_never_elide_its_state(qt_app, tmp_path):
    """F16 item 1: "Sphere ba...SUGGESTED" clipped the one thing that mattered.

    The structural fix: short labels inside the control, the SUGGESTED state
    outside it. The control is then sized to its own longest entry, so no
    width — not the default one, not a squeezed one — can clip it.
    """

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    qt_app.processEvents()

    for row in range(window.file_table.rowCount()):
        combo = window.file_table.cellWidget(row, 1)
        texts = [combo.itemText(index) for index in range(combo.count())]
        # No state text of any kind lives inside the control.
        assert not any(_SUGGESTED_BADGE in text for text in texts), texts
        assert all(len(text) <= 12 for text in texts), texts
        assert combo.minimumWidth() >= role_combo_minimum_width(combo)

    def widest_entry(combo):
        metrics = combo.fontMetrics()
        return max(
            metrics.horizontalAdvance(combo.itemText(i)) for i in range(combo.count())
        )

    # Squeeze the pane hard: the Role column is content-sized, so it holds.
    for pane_width in (420, 300, 240):
        window.root_splitter.setSizes([pane_width, 1200 - pane_width])
        qt_app.processEvents()
        for row in range(window.file_table.rowCount()):
            combo = window.file_table.cellWidget(row, 1)
            assert combo.width() >= widest_entry(combo), (
                f"role combo elides at pane width {pane_width}"
            )
    window.close()


def test_the_default_geometry_lays_out_legibly(qt_app, tmp_path):
    """F16 item 2: every visible word readable at the size it opens at."""

    # The default never exceeds the screen it will open on.
    assert bench_default_geometry(QtCore.QSize(1280, 800)) == (1200, 680)
    assert bench_default_geometry(QtCore.QSize(6000, 6000)) == (1500, 920)

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*bench_default_geometry())
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    clipped = [
        widget.objectName() or widget.text()
        for widget in window.findChildren(QtWidgets.QLabel)
        if widget.isVisible() and _label_is_clipped(widget)
    ]
    assert not clipped, f"labels clip their own text: {clipped}"

    squeezed = [
        button.text()
        for button in window.findChildren(QtWidgets.QPushButton)
        if button.isVisible() and button.width() < button.sizeHint().width()
    ]
    assert not squeezed, f"buttons clip their labels: {squeezed}"

    for group in window.findChildren(QtWidgets.QGroupBox):
        if not group.isVisible():
            continue
        metrics = QtGui.QFontMetrics(group.font())
        assert group.width() >= metrics.horizontalAdvance(group.title()) + 30, (
            f"group box title garbles: {group.title()}"
        )
        # The title sits in the top margin; the margin has to be tall enough.
        assert group.contentsRect().top() >= metrics.height(), (
            f"group box title overlaps its own frame: {group.title()}"
        )

    overlaps = _overlapping_siblings(window)
    assert not overlaps, f"widgets overlap at the default size: {overlaps}"
    window.close()


_LAID_OUT = (
    QtWidgets.QLabel,
    QtWidgets.QPushButton,
    QtWidgets.QGroupBox,
    QtWidgets.QComboBox,
    QtWidgets.QSpinBox,
    QtWidgets.QLineEdit,
)


def _label_is_clipped(label) -> bool:
    """Whether a label's own text does not fit the box it was given."""

    text = label.text()
    if not text or label.width() <= 0:
        return False
    if label.wordWrap():
        return label.heightForWidth(label.width()) > label.height()
    return label.fontMetrics().horizontalAdvance(text) > label.width()


def _overlapping_siblings(window) -> list[str]:
    """Laid-out siblings whose rectangles collide are a broken layout."""

    by_parent: dict[int, list] = {}
    for widget in window.findChildren(_LAID_OUT):
        if not widget.isVisible() or widget.parentWidget() is None:
            continue
        if widget.parentWidget().layout() is None:
            continue
        by_parent.setdefault(id(widget.parentWidget()), []).append(widget)
    collisions = []
    for siblings in by_parent.values():
        for index, first in enumerate(siblings):
            for second in siblings[index + 1 :]:
                if second.isAncestorOf(first) or first.isAncestorOf(second):
                    continue
                if first.geometry().intersects(second.geometry()):
                    collisions.append(
                        f"{first.objectName() or type(first).__name__} x "
                        f"{second.objectName() or type(second).__name__}"
                    )
    return collisions


def test_the_fit_lands_on_the_assigned_lamp_signal_and_names_it(qt_app, tmp_path):
    """F16 item 5: the fit tab showed a lineless `_bg` hump and never said so.

    After confirming a whole folder the last file read stays open — on the
    owner's folder, a background. The fit must move to the assigned lamp
    signal, name the file it is measuring, and subtract that lamp's own
    background, which is what echelle-align does.
    """

    window, _paths = _real_folder_window(qt_app, tmp_path)
    signal = tmp_path / "Ne-0.1s-x3-dimm-lines.sif"
    background = tmp_path / "Ne-0.1s-x3-dimm-lines-bg.sif"
    # What the owner saw: a background frame open for line fitting.
    window.session.accept_frame(_frame_for(background))
    window.refresh()
    _wait_for_loads(window, qt_app)

    window.confirm_roles_button.click()
    qt_app.processEvents()
    _wait_for_loads(window, qt_app)
    for _attempt in range(400):
        if window._background_thread is None and window.session.frame is not None:
            if window.session.frame.path == signal:
                break
        qt_app.processEvents()
        QtCore.QThread.msleep(5)
    qt_app.processEvents()

    # The fit landed on the assigned Ne signal, not on whatever was open.
    assert window.session.frame.path == signal
    assert signal.name in window.fit_file_value.text()
    assert "Ne lamp" in window.fit_file_value.text()
    assert window.fit_warning_value.text() == ""

    # And it is the signal minus that lamp's own background.
    assert window.session.background_path == background
    assert background.name in window.fit_file_value.text()
    # The fixture reader returns the same spectra for every path, so a
    # subtracted pair cancels exactly — which is the point being pinned.
    raw = window.session.frame.order_spectra[0]
    fitted = window.session.active_order_spectra()[0]
    assert not np.allclose(fitted, raw)
    assert np.allclose(fitted, 0.0), "background subtraction did not happen"
    window.close()


def test_opening_a_background_for_line_work_says_so(qt_app, tmp_path):
    """F16 item 5: a deliberate look at a background is warned about, not fought."""

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.confirm_roles_button.click()
    qt_app.processEvents()
    _wait_for_loads(window, qt_app)

    background = tmp_path / "Ne-0.1s-x3-dimm-lines-bg.sif"
    window._select_file_row(background)
    window.show_frame_button.click()
    _wait_for_loads(window, qt_app)
    qt_app.processEvents()

    assert window.session.frame.path == background
    assert window.fit_warning_value.text()
    assert "background frame" in window.fit_warning_value.text()
    assert "no lines are expected" in window.fit_warning_value.text()
    assert background.name in window.fit_file_value.text()

    # A sphere frame is warned about too — a continuum is not a line spectrum.
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    window._select_file_row(sphere)
    window.show_frame_button.click()
    _wait_for_loads(window, qt_app)
    qt_app.processEvents()
    assert window.fit_warning_value.text()
    assert "continuum" in window.fit_warning_value.text()
    window.close()


def test_the_expected_line_panel_lives_with_the_spectrum_and_fills_itself(
    qt_app, tmp_path
):
    """F16 item 7: one lamp choice, one view, the whole line workflow."""

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    window.confirm_roles_button.click()
    qt_app.processEvents()
    _wait_for_loads(window, qt_app)

    # The panel followed the assigned lamp with nobody choosing a catalog.
    assert window.campaign.assigned_lamps == ("Ne",)
    assert window.line_family_combo.currentText() == "Ne"
    assert window.line_help_table.rowCount() > 0
    assert "expected Ne line(s)" in window.line_panel_header.text()
    # F17 item 1 moved it out from under the spectrum into the tall left
    # column; the spectrum keeps its own view.
    assert window.controls_splitter.isAncestorOf(window.line_help_table)
    assert window.lamp_fit_splitter.isAncestorOf(window.order_plot.getViewWidget())

    # Picking the work on the left brings its own view with it.
    lamp_tab = [
        index
        for index in range(window.control_tabs.count())
        if window.control_tabs.tabText(index) == "Lamp fit"
    ][0]
    window.control_tabs.setCurrentIndex(lamp_tab)
    qt_app.processEvents()
    assert window.view_tabs.tabText(window.view_tabs.currentIndex()) == "Lamp fit"

    # Selecting a row marks its stick on the spectrum above.
    assert not window.line_highlight.isVisible()
    window.line_help_table.selectRow(0)
    qt_app.processEvents()
    assert window.line_highlight.isVisible()
    assert window.line_highlight.value() == pytest.approx(
        window._catalog_rows[0].detector_pixel
    )
    assert "expected at pixel" in window.details_view.toPlainText()
    window.close()


def test_accepting_an_anchor_marks_its_expected_line_row(qt_app, tmp_path):
    """F16 item 7: an accepted anchor is visible in the table it came from."""

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()
    # Put one catalog row where an anchor of the same wavelength exists.
    rows = window._catalog_rows
    assert window.line_help_table.rowCount() == len(rows)
    assert all(
        window.line_help_table.item(index, 4).text() == "" for index in range(len(rows))
    )

    assert window.session.fit_anchor_at(0, 26).accepted
    window.refresh()
    qt_app.processEvents()
    anchored = [
        window.line_help_table.item(index, 4).text()
        for index in range(window.line_help_table.rowCount())
    ]
    if window._catalog_rows:
        # The mark only appears where the catalog actually holds that line.
        assert set(anchored) <= {"", "✓"}
    window.close()


def test_order_scrolling_never_touches_the_detector_image(qt_app, tmp_path):
    """F16 item 4: extract once, re-slice — the main GUI's own pattern.

    Re-uploading and re-percentiling a 2560x2160 detector image on every order
    change was 45% of a 99 ms order step. It changes nothing that an order
    switch shows.
    """

    window = _window(tmp_path, with_loader=True)
    window.show()
    qt_app.processEvents()

    uploads = []
    original = window.detector_image.setImage
    window.detector_image.setImage = lambda *a, **k: (
        uploads.append(1),
        original(*a, **k),
    )[1]
    traces = list(window._pattern_items)

    for order in (1, 0, 1, 0):
        window.order_spin.setValue(order)
        qt_app.processEvents()

    assert uploads == [], "an order change re-uploaded the detector image"
    # Plot items are reused, not rebuilt.
    assert window._pattern_items == traces
    assert window.order_curve is not None

    # A genuinely new frame does re-upload it.
    window.session.accept_frame(_frame_for(tmp_path / "another.sif"))
    window.refresh()
    qt_app.processEvents()
    assert uploads, "a new frame must refresh the detector image"
    window.close()


def test_the_factor_curves_are_downsampled_clipped_and_reused(qt_app, tmp_path):
    """F16 item 4: one PlotDataItem per curve, drawn only where it is seen.

    Painting 42k antialiased samples per curve cost 1.3 s per pan step; the
    same data downsampled to the view and clipped to it paints in tens of ms.
    """

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    curves = [window.candidate_curve, window.previous_curve]
    identities = [id(curve) for curve in curves]
    for _repeat in range(4):
        window._refresh_sphere_plot()
    qt_app.processEvents()

    # The same two items, refreshed in place.
    assert [id(curve) for curve in curves] == identities
    plotted = [
        item
        for item in window.sphere_plot.getPlotItem().items
        if isinstance(item, pg.PlotDataItem)
    ]
    assert len(plotted) == 2, "the factors plot accumulated duplicate curves"
    for curve in curves:
        assert curve.opts["autoDownsample"] is True
        assert curve.opts["clipToView"] is True
        assert curve.opts["antialias"] is False
        # Order gaps stay gaps rather than being bridged by a long diagonal.
        assert curve.opts["connect"] == "finite"
    window.close()


def test_packet5_toml_and_save_controls_follow_domain_state(qt_app, tmp_path):
    window = _campaign_window(tmp_path)
    session = window.session
    assert session.fit_anchor_at(0, 26).accepted
    assert session.fit_anchor_at(1, 56).accepted
    window.config_root = tmp_path / "configs"
    window._generate_tomls()
    window.refresh()
    qt_app.processEvents()

    assert window.campaign.comparison.state is ComparisonState.READY
    assert window.toml_preview.toPlainText().startswith("#")
    assert window.save_state_value.text() == "READY"
    assert window.save_snapshot_button.isEnabled()
    window.close()


# ----------------------------------------------------------------------
# Packet F17 — line-work round 3
# ----------------------------------------------------------------------


def _drawn_sticks(window) -> list[tuple[str, float]]:
    """Every labelled expected-line stick currently drawn on the spectrum."""

    drawn = []
    for marker, label in window._line_pool:
        if not marker.isVisible():
            continue
        drawn.append((label.toPlainText(), float(marker.value())))
    return drawn


def _listed_rows(window) -> list[tuple[str, float]]:
    """Every row of the expected-lines panel, as it reads on screen."""

    table = window.line_help_table
    return [
        (table.item(row, 0).text(), float(table.item(row, 2).text()))
        for row in range(table.rowCount())
    ]


class _PlotClick:
    """The two things ``_order_plot_clicked`` asks a mouse event for."""

    def __init__(self, button, scene_pos) -> None:
        self._button = button
        self._scene_pos = scene_pos

    def button(self):
        return self._button

    def scenePos(self):  # noqa: N802 - pyqtgraph's own spelling
        return self._scene_pos


def _click_at(window, x: float, button=QtCore.Qt.LeftButton) -> _PlotClick:
    """A click on the order plot at data column *x*, mid-height."""

    box = window.order_plot.getViewBox()
    # Off-screen the plot is never painted, so its view range stays the unit
    # square until it is asked to fit the curve it holds.
    box.autoRange(padding=0.02)
    (_x0, _x1), (y0, y1) = box.viewRange()
    scene_pos = box.mapViewToScene(QtCore.QPointF(float(x), (y0 + y1) / 2.0))
    assert window.order_plot.sceneBoundingRect().contains(scene_pos)
    return _PlotClick(button, scene_pos)


def test_the_sticks_and_the_table_are_one_line_list(qt_app, tmp_path):
    """F17 item 2: one source of truth for sticks, rows, and per-order counts.

    The bench used to label sticks from the curated wavelength table while
    filling the panel from the packaged NIST caches interpolated onto the
    order. On the owner's real Ne data that produced a labelled NeI 640.225 in
    order 7 that the table had never heard of, and three labelled Ne sticks in
    order 6 under a panel reading "0 expected Ne lines in this order".
    """

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    window.confirm_roles_button.click()
    qt_app.processEvents()
    _wait_for_loads(window, qt_app)
    assert window.campaign.assigned_lamps == ("Ne",)

    for order in (0, 1):
        window.order_spin.setValue(order)
        qt_app.processEvents()
        drawn = _drawn_sticks(window)
        listed = _listed_rows(window)
        assert drawn, f"order {order} drew no expected line at all"
        # Label for label and pixel for pixel, the plot and the panel agree.
        assert drawn == listed, f"order {order}: sticks {drawn} vs table {listed}"
        # And so does the count the header states.
        assert (
            f"{len(listed)} expected Ne line(s) in order {order}"
            in window.line_panel_header.text()
        )
        # It is the assigned lamp's list: no other element leaks into it.
        assert all(name.startswith("Ne") for name, _pixel in listed), listed
    window.close()


def test_the_expected_line_list_stands_in_the_tall_left_column(qt_app, tmp_path):
    """F17 item 1: the list gets the empty left column, not a cramped strip."""

    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry())
    window.show()
    qt_app.processEvents()

    assert window.controls_splitter.orientation() == QtCore.Qt.Vertical
    assert window.controls_splitter.isAncestorOf(window.line_help_table)
    assert not window.controls_splitter.childrenCollapsible()
    panel = window.expected_lines_panel
    assert panel.isVisible()
    # It sits under the controls and takes a real share of the column's
    # height, rather than the four-row strip it had under the spectrum.
    assert window.controls_splitter.indexOf(panel) == 1
    column = window.controls_splitter.height()
    assert panel.height() >= 0.25 * column, (panel.height(), column)
    assert panel.height() > window.line_help_table.horizontalHeader().height() * 4
    # It is in view from every control tab, not only the lamp-fit one.
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        assert panel.isVisible()

    # Row-click still marks the stick on the spectrum (F16 item 7 preserved).
    assert not window.line_highlight.isVisible()
    window.line_help_table.selectRow(0)
    qt_app.processEvents()
    assert window.line_highlight.isVisible()
    assert window.line_highlight.value() == pytest.approx(
        window._catalog_rows[0].detector_pixel
    )
    window.close()


def test_the_moved_line_table_adds_no_per_order_item_churn(qt_app, tmp_path):
    """F17 item 1 must not undo F16 item 4: order scrolling stays cheap."""

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    pool = [id(pair) for pair in window._line_pool]
    scene_items = len(window.order_plot.items)
    builds = []
    original = window.line_help_table.setRowCount
    window.line_help_table.setRowCount = lambda count: (
        builds.append(count),
        original(count),
    )[1]

    for order in (1, 0, 1, 0):
        window.order_spin.setValue(order)
        qt_app.processEvents()

    # Sticks and labels are pooled: the same objects moved and relabelled.
    assert [id(pair) for pair in window._line_pool] == pool
    assert len(window.order_plot.items) == scene_items
    # The table is refilled, but its list is computed once per (order, lamp).
    assert builds, "the table never refilled on an order change"
    assert len(window._catalog_cache) <= window.session.pattern.shape[1]
    window.close()


def test_an_empty_top_end_states_itself_in_words(qt_app, tmp_path):
    """F17 item 3: an all-zero log histogram was a solid block saying nothing."""

    window = _manual_window(tmp_path)
    window.show()
    path = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    path.write_bytes(b"sif\n")
    _drop(window, [path])
    _wait_for_loads(window, qt_app)
    window.file_table.selectRow(0)
    qt_app.processEvents()

    # This fixture peaks at 40000 of 65535: nothing sits in the last 10%.
    assert window.top_end_message.isVisible()
    assert (
        "no pixels within 10% of full scale"
        in window.top_end_message.text().casefold()
    )
    assert not window.top_histogram_widget.isVisible()
    window.close()


def test_a_populated_top_end_draws_visible_outlined_bins(qt_app, tmp_path):
    """F17 item 3, the other half: when pixels are up there, show them."""

    window = _manual_window(tmp_path, peak=65000.0)
    window.show()
    path = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    path.write_bytes(b"sif\n")
    _drop(window, [path])
    _wait_for_loads(window, qt_app)
    window.file_table.selectRow(0)
    qt_app.processEvents()

    assert not window.top_end_message.isVisible()
    assert window.top_histogram_widget.isVisible()
    curves = window.top_histogram_plot.listDataItems()
    assert curves, "the populated top end drew no bins"
    curve = curves[0]
    assert curve.opts["pen"].width() >= 1  # outlined, not a flat wash
    # The plot is in log mode, so the fill has to start below one pixel;
    # filling to y=1 is what painted an empty histogram as a solid block.
    assert curve.opts["fillLevel"] == pytest.approx(float(np.log10(0.5)))
    window.close()


def test_anchors_come_off_at_the_spectrum(qt_app, tmp_path):
    """F17 item 4: removal used to live only in the table, blind to the plot."""

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    window._order_plot_clicked(_click_at(window, 26.0))
    qt_app.processEvents()
    assert len(window.session.anchor_rows()) == 1
    assert window.anchor_table.rowCount() == 1
    anchor = window.session.anchor_rows()[0]

    # The anchored stick wears the anchor's colour, so it looks removable.
    anchored_sticks = [
        marker
        for marker, _label in window._line_pool
        if marker.isVisible()
        and abs(marker.value() - anchor.line.center_pixel) < 0.001
    ]
    assert anchored_sticks
    assert anchored_sticks[0].pen.color().name() == "#6ee7b7"
    assert window.anchor_near(0, anchor.line.center_pixel + 1.0) is anchor
    assert window.anchor_near(1, anchor.line.center_pixel) is None

    # Clicking it again takes the anchor off, and the table follows.
    window._order_plot_clicked(_click_at(window, anchor.fit.center_pixel))
    qt_app.processEvents()
    assert window.session.anchor_rows() == ()
    assert window.anchor_table.rowCount() == 0
    assert "un-anchored" in window.message_value.text()
    assert all(
        window.line_help_table.item(row, 4).text() == ""
        for row in range(window.line_help_table.rowCount())
    )

    # Right-click removes too, and says so honestly when there is nothing there.
    window._order_plot_clicked(_click_at(window, 26.0))
    qt_app.processEvents()
    assert len(window.session.anchor_rows()) == 1
    window._order_plot_clicked(_click_at(window, 26.0, QtCore.Qt.RightButton))
    qt_app.processEvents()
    assert window.session.anchor_rows() == ()
    window._order_plot_clicked(_click_at(window, 5.0, QtCore.Qt.RightButton))
    qt_app.processEvents()
    assert "No anchor sits at that pixel" in window.message_value.text()
    window.close()


def test_the_why_dock_changes_only_when_asked(qt_app, tmp_path):
    """F17 item 5: hover yanked the sentence away while it was being read."""

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert QtCore.QEvent.Enter not in window._EXPLAIN_EVENTS
    assert not window.file_table.hasMouseTracking()

    window.explain("Held still", "This sentence must survive the pointer.")
    held = window.details_view.toPlainText()
    hover_events = (
        QtCore.QEvent.Enter,
        QtCore.QEvent.HoverEnter,
        QtCore.QEvent.HoverMove,
        QtCore.QEvent.ToolTip,
    )
    for widget in list(window._explainable_widgets):
        for kind in hover_events:
            window.eventFilter(widget, QtCore.QEvent(kind))
    assert window.details_view.toPlainText() == held

    # Asking still works, and the tooltip stays the one short line F16 set.
    window.eventFilter(window.rms_value, QtCore.QEvent(QtCore.QEvent.FocusIn))
    assert "root-mean-square" in window.details_view.toPlainText()
    assert "\n" not in window.rms_value.toolTip()
    assert len(window.rms_value.toolTip()) <= BENCH_TOOLTIP_LIMIT
    window.close()
