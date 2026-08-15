"""Focused off-screen tests for the separate calibration-bench window."""

from __future__ import annotations

import os
from datetime import date
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pytest
from PyQt5 import QtCore, QtGui, QtWidgets

from echelle_spectra.calibration_bench import BenchFrame, CalibrationBenchSession
from echelle_spectra.calibration_bench_gui import (
    _PACKAGE_DIR,
    _SUGGESTED_SUFFIX,
    BENCH_BODY_POINT_SIZE,
    BENCH_HEADLINE_POINT_SIZE,
    BENCH_READING_POINT_SIZE,
    CalibrationBenchWindow,
    _build_parser,
    bench_point_sizes,
    bench_window_icon,
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
    lines = (
        CalibrationTableLine(0, 20, 30, 25, 578.0, "ThI", "ok"),
        CalibrationTableLine(0, 65, 75, 70, 640.0, "ArI", "ok"),
        CalibrationTableLine(1, 50, 60, 55, 580.0, "ThI", "ok"),
        CalibrationTableLine(1, 65, 75, 70, 635.0, "ArI", "ok"),
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
    assert [window.view_tabs.tabText(i) for i in range(window.view_tabs.count())] == [
        "Exposure triage",
        "Lamp alignment",
        "Line identification",
        "Sphere factors",
    ]
    assert window.acceptDrops()
    assert window.checklist_tree.count() >= 9
    assert "median" in window.comparison_value.text()
    assert window.line_help_table.rowCount() > 0
    assert "Packaged source" == window.line_help_table.horizontalHeaderItem(4).text()
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
    assert "anomalies" in window.triage_headline.text()
    assert "not saturation" in window.triage_headline.text()
    assert "% of full scale" in window.triage_detail.text()
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
    # The pre-filled control must never read as an assignment nobody made.
    assert "SUGGESTED ONLY" in window.file_table.item(0, 0).text()
    assert role_combo.currentText().endswith(_SUGGESTED_SUFFIX)

    # Picking the already-shown entry emits no index change; the operator's
    # confirmation must still register.
    role_combo.activated.emit(role_combo.currentIndex())
    qt_app.processEvents()

    assert (
        window.campaign.measurements[source].role is MeasurementRole.SPHERE_BACKGROUND
    )
    assert "sphere-background" in window.file_table.item(0, 0).text()
    assert not role_combo.currentText().endswith(_SUGGESTED_SUFFIX)
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
    assert sphere_combo.currentText().startswith("Sphere signal")
    assert sphere_combo.currentText().endswith(_SUGGESTED_SUFFIX)
    assert "SUGGESTED ONLY" in window.file_table.item(sphere_row, 0).text()
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
    assert not sphere_combo.currentText().endswith(_SUGGESTED_SUFFIX)
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
    assert "FIT UNSATURATED LINES ONLY" in window.triage_headline.text()
    assert "SATURATED —" not in window.triage_headline.text()
    assert "SATURATED LINES (EXPECTED)" in window.file_table.item(
        window._file_rows.index(lamp), 0
    ).text()

    # The sphere frame keeps the hard verdict.
    window._select_file_row(sphere)
    qt_app.processEvents()
    assert "SATURATED —" in window.triage_headline.text()
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
        window.triage_detail,
        window.exposure_value,
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


def test_readings_are_sized_to_be_read_at_a_glance(qt_app, tmp_path):
    """F14 item 6: text is data here; the numbers are loud and explained."""

    body, reading, headline = bench_point_sizes()
    assert body >= BENCH_BODY_POINT_SIZE
    assert reading >= BENCH_READING_POINT_SIZE
    assert headline >= BENCH_HEADLINE_POINT_SIZE
    # A large platform font scales every tier with it rather than being ignored.
    assert bench_point_sizes(24.0)[2] > headline

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert window.triage_headline.font().pointSizeF() >= BENCH_HEADLINE_POINT_SIZE
    for widget in (
        window.rms_value,
        window.anchor_count_value,
        window.comparison_value,
        window.file_state_value,
        window.alignment_state_value,
        window.save_state_value,
    ):
        assert widget.font().pointSizeF() >= BENCH_READING_POINT_SIZE
        assert widget.font().bold()
    for widget in (
        window.comparison_value,
        window.rms_value,
        window.anchor_count_value,
        window.transform_value,
        window.frame_choice_value,
        window.confirm_roles_button,
    ):
        assert widget.toolTip(), f"{widget} carries no explanation"

    # The whys live behind a click, in one dock, not as body prose.
    assert window.details_dock is not None
    window.checklist_tree.setCurrentRow(0)
    qt_app.processEvents()
    first = window.campaign.checklist(window.session)[0]
    assert first.label in window.details_view.toPlainText()
    window.explain("Fit RMS", window.rms_value.toolTip())
    assert "root-mean-square" in window.details_view.toPlainText()
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
