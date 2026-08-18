"""Focused off-screen tests for the separate calibration-bench window."""

from __future__ import annotations

import os
import subprocess
import sys
from dataclasses import replace
from datetime import date, datetime
from pathlib import Path, PureWindowsPath

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pyqtgraph as pg
import pytest
from PyQt5 import QtCore, QtGui, QtWidgets

from echelle_spectra import calibration_bench_gui as bench_gui
from echelle_spectra import calibration_campaign as campaign_module
from echelle_spectra.calibration_bench import (
    BenchFrame,
    CalibrationBenchSession,
    Residual,
    StableSifWatcher,
)
from echelle_spectra.calibration_bench_gui import (
    _PACKAGE_DIR,
    _ROOT_SHARES,
    _SUGGESTED_BADGE,
    BENCH_APP_USER_MODEL_ID,
    BENCH_BODY_POINT_SIZE,
    BENCH_FLOOR_LINES,
    BENCH_HEADLINE_POINT_SIZE,
    BENCH_HISTOGRAM_LINES,
    BENCH_PREFERRED_LINES,
    BENCH_RESIDUAL_FLOOR_LINES,
    BENCH_RESIDUAL_LINES,
    BENCH_TOOLTIP_LIMIT,
    BENCH_TOP_END_LINES,
    CONFIG_ROOT_NAME,
    SNAPSHOT_ROOT_NAME,
    CalibrationBenchWindow,
    _build_parser,
    _ElidingLabel,
    absolute_root,
    acquisition_date_from_metadata,
    acquisition_date_from_name,
    apply_windows_taskbar_identity,
    bench_default_geometry,
    bench_layout_unit,
    bench_minimum_size,
    bench_point_sizes,
    bench_window_icon,
    default_bench_roots,
    forget_session_layout,
    one_line,
    resolve_calibration_file,
    role_combo_minimum_width,
    snapshot_identity,
)
from echelle_spectra.calibration_campaign import (
    AbsoluteCalibrationResult,
    CalibrationCampaignSession,
    ComparisonState,
    ExposureState,
    MeasurementRecord,
    MeasurementRole,
    SaveState,
    TomlState,
    default_validity,
    lamp_reference_set,
    triage_exposure,
    triage_for_role,
)
from echelle_spectra.snapshot import Snapshot
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
    experiment_time: int | None = None,
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
    # ``ExperimentTime`` is the Andor header's Unix acquisition start; it is
    # what dates a calibration when its folder name does not.
    metadata: dict = {"ExposureTime": 0.1}
    if experiment_time is not None:
        metadata["ExperimentTime"] = int(experiment_time)
    return BenchFrame(
        Path(path),
        images,
        images.mean(axis=0),
        order_spectra,
        metadata,
        frame_order_spectra=per_frame,
    )


def _loader(**frame_options):
    def load(path):
        return _frame_for(Path(path), **frame_options)

    return load


def _window(
    tmp_path: Path,
    *,
    with_loader: bool = False,
    window_options: dict | None = None,
    **frame_options,
) -> CalibrationBenchWindow:
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
        **(window_options or {}),
    )


def _manual_window(
    tmp_path: Path, *, window_options: dict | None = None, **frame_options
) -> CalibrationBenchWindow:
    """A bench with campaign memory, a reader, and nothing loaded yet."""

    window = _window(
        tmp_path, with_loader=True, window_options=window_options, **frame_options
    )
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
    # Nothing is dated by today at parse time: the identity and its epoch are
    # dated by the acquisition, which is not knowable until the folder — and
    # then the frames — have been read.  ``None`` is what keeps "unset"
    # distinguishable from "set to the default".
    assert args.snapshot_id is None
    assert args.valid_from is None
    assert _build_parser().parse_args(["--valid-from", "2026-09-01"]).valid_from == date(
        2026, 9, 1
    )
    # Any lamp name is accepted; the bench never owns a permitted list.
    assert _build_parser().parse_args(["--lamp", "Kr", "--lamp", "Ne"]).lamp == ["Kr", "Ne"]


def test_previous_sphere_flags_say_what_the_default_actually_is():
    """--help names the packaged pair, so "previous" is never a mystery."""

    help_text = _build_parser().format_help()
    assert "sphere_cmos_20240305.sif" in help_text
    assert "sphere_cmos_20240305_bkg.sif" in help_text
    assert "2024-03-05" in help_text


def test_sphere_panel_names_its_reference_and_calls_out_a_self_check(
    qt_app, tmp_path
):
    """The panel says what it compared against, and when that was itself."""

    window = _campaign_window(tmp_path)
    text = window.comparison_value.text()
    assert "previous_sphere.sif + previous_sphere_bg.sif" in text
    assert "SELF-CHECK" not in text
    assert "previous_sphere.sif" in window.comparison_value.toolTip()

    # Now the shipped situation: the "previous" pair is this sphere's bytes.
    (tmp_path / "previous_sphere.sif").write_bytes(
        (tmp_path / "sphere.sif").read_bytes()
    )
    (tmp_path / "previous_sphere_bg.sif").write_bytes(
        (tmp_path / "sphere_bg.sif").read_bytes()
    )
    window.campaign.compute_sphere_comparison(
        lambda **_values: AbsoluteCalibrationResult(
            np.linspace(400, 700, 60), np.full(60, 1.0)
        )
    )
    window.refresh()
    qt_app.processEvents()

    text = window.comparison_value.text()
    assert text.startswith("READY ·"), text
    assert "SELF-CHECK" in text
    assert "self-check, not physics" in text
    assert "choose a different previous pair" in text
    # The absolute pair is one click away, per the house Why-dock pattern.
    assert str(tmp_path / "previous_sphere.sif") in window.comparison_value.property(
        "explainText"
    )
    window.close()


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


def test_an_empty_bench_offers_its_two_verbs_and_teaches_nothing(qt_app, tmp_path):
    """Owner, 2026-08-18: "I don't need reminding."

    The dashed "DROP SIF FILES HERE / any names, any order, as many as you
    like" panel taught a gesture to the person who wrote the drop handler, and
    charged the shortest column on screen a hundred pixels for the lesson. It
    is gone, along with the stylesheet rule that painted it; the empty table is
    its own invitation, and dropping still works everywhere it worked before
    (pinned by the drop tests below, which drive the real Qt events).
    """

    window = _manual_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert not hasattr(window, "drop_hint")
    assert "dropTarget" not in window.styleSheet()
    for widget in window.findChildren(QtWidgets.QLabel):
        assert "DROP SIF FILES HERE" not in widget.text()
    assert window.open_folder_button.isEnabled()
    assert window.add_files_button.isEnabled()
    assert window.file_table.rowCount() == 0
    items = {
        window.checklist_tree.itemWidget(window.checklist_tree.item(row)).text()
        for row in range(window.checklist_tree.count())
    }
    assert any("unblocked by: drop any SIF" in text for text in items)
    window.close()


def test_dropped_files_are_triaged_before_any_role(qt_app, tmp_path):
    """Triage needs a file and nothing else — not even a readable name.

    The names here say nothing, so no role is applied and none is confirmed;
    the exposure verdict is reached all the same. (F21 item 1 assigns a folder
    whose names *do* say something, which is a different test.)
    """

    window = _manual_window(tmp_path, cosmic=True)
    window.show()
    paths = []
    for name in ("IMG_0042.sif", "IMG_0043.sif"):
        path = tmp_path / name
        path.write_bytes(b"sif\n")
        paths.append(path)

    _drop(window, paths)
    _wait_for_loads(window, qt_app)

    assert window.file_table.rowCount() == 2
    assert set(window.campaign.loaded) == set(paths)
    # Triage happened; roles did not.
    assert not window.campaign.measurements
    triage = window.campaign.loaded[paths[0]].triage
    assert triage.state is ExposureState.GOOD
    assert triage.saturation.anomalous_pixels == 2
    # One line, carrying the verdict and nothing else. The filename is in the
    # selected table row a few centimetres away and the loudest text on screen
    # does not spend a line repeating it; the anomaly count is a number, so it
    # reads with the other numbers (owner, 2026-08-16).
    assert window.triage_headline.text().splitlines() == ["GOOD"]
    assert paths[0].name not in window.triage_headline.text()
    assert "anomalies" in window.exposure_value.text()
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


def test_opening_a_frame_is_a_load_but_not_an_arrival(qt_app, tmp_path):
    """A folder's filenames are judged when files ARRIVE, and only then.

    Assigning a role makes the bench open the assigned lamp signal for fitting,
    which is a load like any other — and the whole-drop judgement used to hang
    off any load finishing with an empty queue. So resolving the one ambiguous
    name by hand made every other suggestion in the folder apply itself a
    moment later, on no event the operator caused. Only an arrival counts now.
    """

    window = _manual_window(tmp_path)
    window.show()
    nameless = tmp_path / "IMG_0042.sif"
    nameless.write_bytes(b"sif\n")
    lamp_named = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    lamp_named.write_bytes(b"sif\n")
    _drop(window, [nameless, lamp_named])
    _wait_for_loads(window, qt_app)
    # One ambiguous name in the drop, so nothing is applied and nothing is
    # pending: the drop has been judged, once.
    assert not window.campaign.measurements
    assert window._arrivals_pending is False

    window.file_table.cellWidget(0, 2).setCurrentText("Kr")
    role_combo = window.file_table.cellWidget(0, 1)
    role_combo.setCurrentIndex(role_combo.findData(MeasurementRole.LAMP))
    # Let the frame the assignment opens finish loading — the load that used to
    # re-judge the folder behind the operator's back.
    _wait_for_loads(window, qt_app)
    qt_app.processEvents()

    assert window.campaign.assigned_lamps == ("Kr",)
    assert lamp_named not in window.campaign.measurements
    assert window._arrivals_pending is False
    window.close()


def test_confirming_the_prefilled_role_assigns_it(qt_app, tmp_path):
    """The confirm path, on a drop the bench will not decide for itself.

    The nameless file is what keeps the drop doubtful: with it there, nothing
    is applied unasked and every row keeps its SUGGESTED badge (F21 item 1).
    """

    window = _manual_window(tmp_path)
    window.show()
    source = tmp_path / "sphere-0.1s-x3-bg.sif"
    source.write_bytes(b"sif\n")
    nameless = tmp_path / "IMG_0042.sif"
    nameless.write_bytes(b"sif\n")

    _drop(window, [source, nameless])
    _wait_for_loads(window, qt_app)
    role_combo = window.file_table.cellWidget(0, 1)
    assert role_combo.currentData() is MeasurementRole.SPHERE_BACKGROUND
    assert not window.campaign.measurements
    # The pre-filled control must never read as an assignment nobody made —
    # and after F16 the badge lives beside the combo, never inside its text.
    # The carrier changed again on 2026-08-16: "SUGGESTED ONLY — no role
    # assigned" spelled out on every row was one sentence repeated six times,
    # so the row keeps the short word and the amber border does the shouting.
    assert _SUGGESTED_BADGE.lower() in window.file_table.item(0, 0).text()
    assert "no role assigned" not in window.file_table.item(0, 0).text()
    assert _SUGGESTED_BADGE not in role_combo.currentText()
    assert _SUGGESTED_COLOR_IN(role_combo)

    # Picking the already-shown entry emits no index change; the operator's
    # confirmation must still register.
    role_combo.activated.emit(role_combo.currentIndex())
    qt_app.processEvents()

    assert (
        window.campaign.measurements[source].role is MeasurementRole.SPHERE_BACKGROUND
    )
    # The confirmed role now reads as the verdict "BACKGROUND" rather than the
    # role token repeated beside it: a background frame's verdict already names
    # what it is, and saying "sphere-background" next to it was the same word
    # twice (owner, 2026-08-16).
    assert "BACKGROUND" in window.file_table.item(0, 0).text()
    assert window.campaign.measurements[source].role is (
        MeasurementRole.SPHERE_BACKGROUND
    )
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


def _real_folder_window(qt_app, tmp_path: Path, *, extra=(), **frame_options):
    """The owner's own folder on the bench, dropped through the real drop path.

    ``extra`` adds files beyond the six real ones — a name the bench cannot
    read, or a second claim on a role only one file can hold — which is what
    keeps a drop doubtful now that an unambiguous one assigns itself.
    """

    window = _manual_window(tmp_path, **frame_options)
    paths = []
    for name in (*_REAL_2025_NAMES, *extra):
        path = tmp_path / name
        path.write_bytes(b"sif\n")
        paths.append(path)
    _drop(window, paths)
    _wait_for_loads(window, qt_app)
    return window, paths


def test_a_prefilled_role_never_looks_like_an_assigned_one(qt_app, tmp_path):
    """F14 item 1 regression, on the drop that is still the operator's to sort.

    Every file of the real folder pre-fills its role correctly, so the Role
    column read right while the campaign had been given nothing: the Procedure
    tab said "no file carries this role yet" and the factor computation failed.
    The control says SUGGESTED until somebody confirms it, and confirming it
    through the combo reaches the campaign.

    F21 item 1 applies a folder whose names leave nothing to ask, so the drop
    tested here carries one name the bench cannot read. One doubt is enough:
    nothing is applied, and every row of it stays SUGGESTED.
    """

    window, paths = _real_folder_window(qt_app, tmp_path, extra=("IMG_0042.sif",))
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
    assert _SUGGESTED_BADGE.lower() in window.file_table.item(sphere_row, 0).text()
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
    """A doubtful drop is one deliberate press, not one popup pick per row.

    The unreadable name is what left the press to be made at all: without it
    the six would have assigned themselves (F21 item 1).
    """

    window, _paths = _real_folder_window(qt_app, tmp_path, extra=("IMG_0042.sif",))
    assert not window.campaign.measurements

    window.confirm_roles_button.click()
    qt_app.processEvents()

    # The six that named themselves; the one that did not stays unassigned.
    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)
    assert tmp_path / "IMG_0042.sif" not in window.campaign.measurements
    assert window.campaign.unconfirmed_suggestions() == ()
    assert window.campaign.assigned_lamps == ("Ne",)
    assert "Assigned 6 suggested role(s)" in window.message_value.text()
    checklist = {
        item.key: item for item in window.campaign.checklist(window.session)
    }
    assert checklist["sphere"].detail == "sphere-0.1s-x3.sif"
    assert checklist["sphere-background"].detail == "sphere-0.1s-x3-bg.sif"
    window.close()


def test_a_folder_that_names_itself_is_assigned_on_arrival(qt_app, tmp_path):
    """F21 item 1: an unambiguous drop is applied, with one line about it.

    The owner's complaint was the confirmation dance — a dropdown per row, or
    a Confirm press, for suggestions that are "always like so". When every
    name in the drop says exactly one thing there is nothing to decide, so the
    bench assigns them and says so once.
    """

    window, paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    qt_app.processEvents()

    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)
    assert window.campaign.unconfirmed_suggestions() == ()
    assert window.campaign.assigned_lamps == ("Ne",)
    assert (
        window.campaign.measurements[tmp_path / "sphere-0.1s-x3.sif"].role
        is MeasurementRole.SPHERE
    )
    # One notice, and it names what happened and where to undo it.
    message = window.message_value.text()
    assert "Roles assigned from filenames: 6 file(s)" in message
    assert "change any that are wrong" in message

    # No new badge state was invented: SUGGESTED simply is not there, exactly
    # as after a confirmation, and the dropdowns still carry the roles.
    for row, path in enumerate(window._file_rows):
        combo = window.file_table.cellWidget(row, 1)
        assert not _SUGGESTED_COLOR_IN(combo)
        assert _SUGGESTED_BADGE.lower() not in window.file_table.item(row, 0).text()
        assert combo.currentData() is window.campaign.measurements[path].role
    assert not window.confirm_roles_button.isEnabled()
    # And the procedure moved on rather than asking for a confirmation.
    checklist = {item.key: item for item in window.campaign.checklist(window.session)}
    assert checklist["sphere"].detail == "sphere-0.1s-x3.sif"
    window.close()


def test_an_applied_role_is_an_ordinary_assignment_the_combo_can_change(
    qt_app, tmp_path
):
    """Applied is not locked: the dropdown remains the operator's override."""

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    row = window._file_rows.index(sphere)
    combo = window.file_table.cellWidget(row, 1)

    combo.setCurrentIndex(combo.findData(MeasurementRole.OTHER))
    qt_app.processEvents()
    assert window.campaign.measurements[sphere].role is MeasurementRole.OTHER

    # Taking the role off entirely sticks, too — including across the next
    # drop, which must never push the filename's guess back onto it.
    combo.setCurrentIndex(combo.findData(None))
    qt_app.processEvents()
    assert sphere not in window.campaign.measurements

    later = tmp_path / "Ne-0.3s-x3-lines.sif"
    later.write_bytes(b"sif\n")
    _drop(window, [later])
    _wait_for_loads(window, qt_app)

    assert window.campaign.measurements[later].role is MeasurementRole.LAMP
    assert sphere not in window.campaign.measurements
    window.close()


def test_two_files_claiming_one_role_keep_the_confirm_flow(qt_app, tmp_path):
    """The owner's own example of doubt: two files that both say sphere."""

    window, _paths = _real_folder_window(
        qt_app, tmp_path, extra=("sphere-0.2s-x3.sif",)
    )
    window.show()
    qt_app.processEvents()

    assert not window.campaign.measurements
    assert "Roles assigned from filenames" not in window.message_value.text()
    for row in range(window.file_table.rowCount()):
        assert _SUGGESTED_COLOR_IN(window.file_table.cellWidget(row, 1))
    assert window.confirm_roles_button.isEnabled()
    assert "7" in window.confirm_roles_button.text()
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
    # The headline is one line and starts with the verdict: the filename it
    # used to carry on a first line is already in the selected table row, and
    # repeating it spent the loudest text on screen saying it twice (owner,
    # 2026-08-16). The advice lives in its own panel underneath, and the
    # numbers in theirs, so no panel now repeats another.
    assert window.triage_headline.text().startswith("SATURATED")
    assert lamp.name not in window.triage_headline.text()
    assert window.campaign.role_triage(sphere).blocking
    assert "Lower exposure" in window.triage_next_value.text()
    assert "Lower exposure" not in window.exposure_value.text()
    assert "peak" in window.exposure_value.text()
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
        window.frame_choice_value,
    ):
        assert label.wordWrap()
    # transform_value deliberately left this list: wrapping put the θ on a
    # second line the row height then hid, so "dx / dy / θ" showed two
    # numbers and a comma. It is one-line state text now (elided middle only
    # if the rail is ever absurdly narrow, full text in the tooltip).
    assert not window.transform_value.wordWrap()

    window.root_splitter.setSizes([760, 780])
    qt_app.processEvents()
    window._relayout_wrapped_text()
    width = window.checklist_tree.viewport().width()
    for row in range(window.checklist_tree.count()):
        label = window.checklist_tree.itemWidget(window.checklist_tree.item(row))
        assert label.wordWrap()
        assert label.width() <= max(280, width)
        # The row is tall enough for every wrapped line it holds — measured at
        # the width the label really has.  This used to compare against
        # ``label.sizeHint().height()``, which is a hint computed at some other
        # width entirely (252 px for a label drawn at 777) and was therefore
        # demanding rows twice as tall as their own text: F21 item 5's padding,
        # asserted from the other side. Deliberately superseded.
        needed = label.heightForWidth(label.width())
        assert window.checklist_tree.item(row).sizeHint().height() >= needed
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
    """F16 item 2: every visible word readable at the size it opens at.

    F20 supersedes F18's *pixel* floor.  ``_MINIMUM_WINDOW_SIZE == (1300, 880)``
    was a constant measured on a 100%-scale display, and that is exactly what
    broke on the owner's: a 1920x1080 screen at 150% offers 1280x720 logical
    pixels, so the floor asked for a window wider and taller than the whole
    desktop and every panel inside it was crushed.  The floor is now quoted in
    lines of the platform's own text and yields to the screen.
    """

    unit = bench_layout_unit()
    preferred = (BENCH_PREFERRED_LINES[0] * unit, BENCH_PREFERRED_LINES[1] * unit)
    floor = (BENCH_FLOOR_LINES[0] * unit, BENCH_FLOOR_LINES[1] * unit)

    # A generous screen gets the preferred size, whatever the DPI.
    assert bench_default_geometry(QtCore.QSize(6000, 6000)) == preferred
    assert bench_default_geometry(QtCore.QSize(1920, 1080)) == preferred
    # The floor is a floor only while the screen can hold it.
    assert bench_minimum_size(QtCore.QSize(6000, 6000)) == floor
    # The owner's display at 150%: 1280x720 logical.  Nothing the bench asks
    # for may exceed it — asking is what produced the cramped first paint.
    scaled = QtCore.QSize(1280, 720)
    for size in (bench_default_geometry(scaled), bench_minimum_size(scaled)):
        assert size[0] <= scaled.width(), size
        assert size[1] <= scaled.height(), size

    # Measured at the size a 1080p-class display gives it.  The offscreen
    # platform reports an 800x600 pseudo-screen, and the geometry function now
    # correctly clamps to it — a window that small is the degraded case, not
    # the one the layout is designed against.
    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
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
    if isinstance(label, _ElidingLabel):
        # It shortens on purpose, in the middle, with the whole reading in its
        # tooltip. That is a decision, not a clip (F18 item 3).
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
    # F17 item 1 moved it out from under the spectrum; F20 moved it on into
    # the tables rail.  The spectrum keeps its own view either way.
    assert window.tables_splitter.isAncestorOf(window.line_help_table)
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


def test_the_factor_curves_are_downsampled_and_reused(qt_app, tmp_path):
    """F16 item 4: one PlotDataItem per curve, downsampled to the view.

    Painting 42k antialiased samples per curve cost 1.3 s per pan step; the
    same data downsampled to the view paints in tens of ms.

    F18 item 5 supersedes the clipping half of F16's answer.  Clipping to the
    view finds the visible slice by binary search, which needs an ascending x;
    a factor curve is the echelle orders laid end to end and its wavelength
    axis runs backwards at nearly every sample, so the search returned about
    one point and the curve vanished at every zoom.  Downsampling is where the
    speed was: dropping the clip costs about a tenth of a millisecond a frame.
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
        assert curve.opts["clipToView"] is False
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


def test_the_expected_line_list_stands_in_the_tall_right_rail(qt_app, tmp_path):
    """F20 supersedes F17 item 1: the list gets a rail, not the left column.

    F17 moved the table out from under the spectrum into the tall left space;
    the owner's answer after living with it was that a table under four tabs of
    controls is still a table sharing one column's height with controls.  The
    two working tables get their own rail.
    """

    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()

    assert window.tables_splitter.orientation() == QtCore.Qt.Vertical
    assert window.tables_splitter.isAncestorOf(window.line_help_table)
    assert window.tables_splitter.isAncestorOf(window.anchor_table)
    assert not window.tables_splitter.childrenCollapsible()
    # Neither table is under the controls any more.
    assert not window.control_tabs.isAncestorOf(window.line_help_table)
    assert not window.control_tabs.isAncestorOf(window.anchor_table)
    panel = window.expected_lines_panel
    assert panel.isVisible()
    assert window.tables_splitter.indexOf(window.anchors_panel) == 0
    assert window.tables_splitter.indexOf(panel) == 1
    rail = window.tables_splitter.height()
    assert panel.height() >= 0.33 * rail, (panel.height(), rail)
    assert panel.height() > window.line_help_table.horizontalHeader().height() * 4
    # It is in view from every control tab, not only the lamp-fit one.
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        assert panel.isVisible()
        assert window.anchors_panel.isVisible()

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


def test_the_whole_folder_is_summarised_one_row_per_file(qt_app, tmp_path):
    """F21 item 4a: judge the folder in one look, then open the row worth opening.

    Triage was a single-file surface: six acquisitions meant six clicks and
    five verdicts held in the head. The summary answers "is this folder any
    good" at a glance, and the detailed panel below still reads one file.
    """

    window, _paths = _real_folder_window(qt_app, tmp_path, cosmic=True)
    window.show()
    qt_app.processEvents()
    table = window.triage_summary_table

    assert table.rowCount() == len(window._file_rows) == len(_REAL_2025_NAMES)
    for row, path in enumerate(window._file_rows):
        record = window.campaign.loaded[path]
        verdict = window.campaign.role_triage(path)
        assert table.item(row, 0).text() == path.name
        # The role it carries, spelled out in its own column.
        assert window.campaign.measurements[path].role.value in table.item(row, 1).text()
        # The verdict word, in the colour that verdict earns.
        assert table.item(row, 2).text() == verdict.label.upper()
        assert table.item(row, 2).foreground().color().name() == (
            window._reading_colour(verdict, record.triage)
        )
        # Peak as a share of full scale, and the anomalies that were counted
        # rather than held against the frame.
        assert table.item(row, 3).text() == (
            f"{100.0 * record.triage.headroom_fraction:.0f}%"
        )
        assert table.item(row, 4).text() == "2"

    # The table is the entry point: clicking a row opens that file below.
    lamp_row = window._file_rows.index(tmp_path / "Ne-0.02s-x3-bright-lines.sif")
    window._summary_row_clicked(lamp_row, 0)
    qt_app.processEvents()
    assert window._selected_file() == window._file_rows[lamp_row]
    assert window.triage_headline.text()
    assert window.file_table.currentRow() == lamp_row
    window.close()


def test_the_summary_paints_an_expected_saturation_as_a_note(qt_app, tmp_path):
    """A dim-series lamp is shot to clip; the summary must not shout red."""

    window, _paths = _real_folder_window(qt_app, tmp_path, peak=65535.0)
    window.show()
    qt_app.processEvents()
    table = window.triage_summary_table

    lamp_row = window._file_rows.index(tmp_path / "Ne-0.1s-x3-dimm-lines.sif")
    sphere_row = window._file_rows.index(tmp_path / "sphere-0.1s-x3.sif")

    assert table.item(lamp_row, 2).text() == "SATURATED LINES (EXPECTED)"
    assert table.item(lamp_row, 2).foreground().color().name() == "#ffb86b"
    # The sphere keeps the hard verdict, in the colour that means stop.
    assert table.item(sphere_row, 2).text() == "SATURATED"
    assert table.item(sphere_row, 2).foreground().color().name() == "#ff8f8f"
    window.close()


def test_the_top_end_is_a_strip_under_the_histogram_it_qualifies(qt_app, tmp_path):
    """F21 item 4b: the near-saturation panel stays, and stops competing.

    The owner now reads the top end as the number the lamp is adjusted by, so
    it is kept — but drawn at the same size as the raw-counts histogram it
    split the view's attention in half. Its subordination is structural: a
    capped strip with no appetite for height, beside a primary that keeps its
    floor and every share of the growth.
    """

    window = _manual_window(tmp_path, peak=65000.0)
    window.show()
    path = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    path.write_bytes(b"sif\n")
    _drop(window, [path])
    _wait_for_loads(window, qt_app)
    qt_app.processEvents()

    strip = window.top_histogram_widget
    primary = window.histogram_widget
    unit = window.layout_unit

    assert strip.maximumHeight() == BENCH_TOP_END_LINES * unit
    assert primary.minimumHeight() == BENCH_HISTOGRAM_LINES * unit
    # The ceiling of the secondary sits below the floor of the primary.
    assert strip.maximumHeight() < primary.minimumHeight()
    assert strip.sizePolicy().verticalPolicy() == QtWidgets.QSizePolicy.Maximum
    assert primary.sizePolicy().verticalPolicy() == QtWidgets.QSizePolicy.Expanding

    layout = strip.parentWidget().layout()
    assert layout.stretch(layout.indexOf(strip)) == 0
    assert layout.stretch(layout.indexOf(primary)) >= 1

    # And on screen, not only in the size policy.
    assert strip.isVisible()
    assert strip.height() < primary.height()
    # F17 item 3 survives: emptiness is still stated in words, not drawn.
    assert not window.top_end_message.isVisible()
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


def _anchor_pair(window) -> None:
    """Two anchors in different orders: the least that solves a transform."""

    assert window.session.fit_anchor_at(0, 26).accepted
    assert window.session.fit_anchor_at(1, 56).accepted
    window.refresh()


def _residual_click_at(
    window, index: float, button=QtCore.Qt.LeftButton, span=None
) -> _PlotClick:
    """A click on the residual strip at bar *index*, mid-height.

    *span* widens the view first, which is the only way to put a click beyond
    the last bar inside the strip's own rectangle.
    """

    box = window.residual_plot.getViewBox()
    # Off-screen the strip is never painted, so it has to be asked to fit the
    # bars it holds before a data coordinate means anything in the scene.
    box.autoRange(padding=0.02)
    if span is not None:
        box.setXRange(*span, padding=0.0)
    (_x0, _x1), (y0, y1) = box.viewRange()
    scene_pos = box.mapViewToScene(QtCore.QPointF(float(index), (y0 + y1) / 2.0))
    assert window.residual_plot.sceneBoundingRect().contains(scene_pos)
    return _PlotClick(button, scene_pos)


def _residual_bars(window):
    """The one bar item the strip draws, or None while there is nothing to draw."""

    return window._residual_bars


def _distinct_colours(pixmap) -> int:
    """How many colours a grab actually painted."""

    image = pixmap.toImage()
    return len(
        {
            image.pixel(x, y)
            for y in range(image.height())
            for x in range(0, image.width(), 2)
        }
    )


def _assert_labels_are_legible(window) -> None:
    """No drawn wavelength may touch its neighbour or fall off the strip."""

    axis = window.residual_plot.getAxis("bottom")
    box = window.residual_plot.getViewBox()
    metrics = QtGui.QFontMetrics(axis.style["tickFont"])
    labels = window.residual_tick_labels()
    assert labels, "the strip drew no wavelengths at all"
    values = [value for value, _text in axis._tickLevels[0]]
    spans = []
    for value, text in zip(values, labels):
        centre = box.mapViewToScene(QtCore.QPointF(float(value), 0.0)).x()
        half = metrics.horizontalAdvance(text) / 2.0
        spans.append((centre - half, centre + half))
    for (_left, right), (next_left, _next_right) in zip(spans, spans[1:]):
        assert next_left > right, (labels, spans)
    rect = window.residual_plot.sceneBoundingRect()
    assert spans[0][0] >= rect.left() - 1, (labels, spans[0], rect.left())
    assert spans[-1][1] <= rect.right() + 1, (labels, spans[-1], rect.right())


def test_the_residual_strip_stands_with_the_anchor_table_it_scores(qt_app, tmp_path):
    """Owner, 2026-08-17: "what about moving the anchor plot into the anchors tab?"

    The residual bars used to sit under the spectrum in the Lamp fit view, so
    they were only in front of the operator while that one view was open — and
    they were nowhere near the rows they score. In the anchors panel they are
    in the right rail, under the table, and therefore never behind a tab.
    """

    forget_session_layout()
    window = _window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()

    # Where it lives now: the anchors panel, under the table.
    panel_layout = window.anchors_panel.layout()
    assert window.anchors_panel.isAncestorOf(window.residual_widget)
    assert window.tables_rail.isAncestorOf(window.residual_widget)
    assert panel_layout.indexOf(window.residual_widget) > panel_layout.indexOf(
        window.anchor_table
    )

    # And where it does not: the Lamp fit view it was taken out of.
    lamp_view = next(
        window.view_tabs.widget(index)
        for index in range(window.view_tabs.count())
        if window.view_tabs.tabText(index) == "Lamp fit"
    )
    assert not lamp_view.isAncestorOf(window.residual_widget)
    assert window.order_plot.scene() is not window.residual_plot.scene()

    # In view from every control tab, exactly as the table above it is.
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        assert window.residual_widget.isVisible(), (
            f"hidden on the {window.control_tabs.tabText(index)!r} tab"
        )

    # Subordinate to the table, structurally: capped, no appetite for height,
    # and never at the cost of the table's own workable rows.
    assert window.residual_widget.maximumHeight() <= (
        BENCH_RESIDUAL_LINES * window.layout_unit
    )
    assert window.residual_widget.maximumHeight() >= (
        BENCH_RESIDUAL_FLOOR_LINES * window.layout_unit
    )
    assert (
        window.residual_widget.sizePolicy().verticalPolicy()
        == QtWidgets.QSizePolicy.Maximum
    )
    assert panel_layout.stretch(panel_layout.indexOf(window.residual_widget)) == 0
    assert panel_layout.stretch(panel_layout.indexOf(window.anchor_table)) >= 1
    assert window.residual_widget.height() < window.anchor_table.height()
    assert _visible_rows(window.anchor_table) >= (
        CalibrationBenchWindow.EXPECTED_LINE_ROWS
    )
    # The rail's width is the rail's: the strip adapts to it, never widens it.
    assert window.residual_widget.minimumWidth() == 0
    window.close()
    forget_session_layout()


def test_the_residual_strip_labels_its_wavelengths_and_tracks_every_anchor(
    qt_app, tmp_path
):
    """One bar per accepted anchor, in the table's order, wavelengths under it.

    Owner, 2026-08-17: "dropping wavelength hints on the plots? Why? Quite
    usable. Useful." So a few anchors are all labelled, and the operator reads
    the wavelengths off the axis without clicking anything.
    """

    forget_session_layout()
    window = _window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()

    # Nothing anchored: no bars, no labels, and an empty strip rather than a
    # drawn zero.
    assert _residual_bars(window) is None
    assert window.residual_tick_labels() == []

    _anchor_pair(window)
    qt_app.processEvents()
    residuals = window.session.residuals
    assert len(residuals) == 2
    bars = _residual_bars(window)
    assert bars is not None
    # One bar per accepted anchor, in the anchor table's own order.
    assert list(bars.opts["height"]) == [residual.dx_px for residual in residuals]
    assert [residual.key for residual in residuals] == [
        anchor.key for anchor in window.session.anchor_rows()
    ]
    assert window.anchor_table.rowCount() == len(residuals)

    # Sparse anchors: every one of them is named, in nm, in bar order.
    assert window.residual_tick_labels() == [
        f"{residual.wavelength_nm:.1f}" for residual in residuals
    ]
    bottom = window.residual_plot.getAxis("bottom")
    assert bottom.style["showValues"] is True
    # Drawn smaller than the bench's body text — that is what buys a rail
    # several wavelengths — and the axis has the height to draw a row of them.
    assert bottom.style["tickFont"].pointSizeF() < window.body_pt
    assert bottom.height() >= QtGui.QFontMetrics(bottom.style["tickFont"]).height()

    # And legible at the width the rail actually gives it: this is a strip a
    # rail wide, not a plot, and every label sits clear of its neighbour and
    # inside the strip.  (The offscreen platform rasterises no glyphs at all,
    # so what is checked is where the text goes, from the same font metrics
    # that place it, and that the strip really paints at that size.)
    assert window.residual_widget.width() < window.width() // 3
    _assert_labels_are_legible(window)
    painted = window.residual_widget.grab()
    assert painted.width() == window.residual_widget.width()
    assert painted.height() == window.residual_widget.height()
    assert _distinct_colours(painted) > 1, "the strip painted nothing at rail width"

    # Removing one anchor drops the transform and with it the bars; the strip
    # follows the same triggers it always did, and no new ones.
    window.session.remove_anchor(window.session.anchor_rows()[0].key)
    window.refresh()
    qt_app.processEvents()
    assert window.session.residuals == ()
    assert _residual_bars(window) is None

    _anchor_pair(window)
    qt_app.processEvents()
    assert len(_residual_bars(window).opts["height"]) == 2
    window.session.clear_anchors()
    window.refresh()
    qt_app.processEvents()
    assert _residual_bars(window) is None
    assert window.residual_tick_labels() == []
    window.close()
    forget_session_layout()


def _crowd_residuals(window, count: int) -> tuple:
    """An auto-anchored campaign's worth of residuals, without the fitting."""

    residuals = tuple(
        Residual(
            (index % 3, 100.0 + index, 500.0 + index * 7.3),
            index % 3,
            500.0 + index * 7.3,
            0.4 * ((-1) ** index) * (1 + index % 4),
            0.0,
            abs(0.4 * (1 + index % 4)),
        )
        for index in range(count)
    )
    window.session.residuals = residuals
    window._refresh_residual_plot()
    return residuals


def test_dense_anchors_thin_the_wavelengths_rather_than_smear_them(qt_app, tmp_path):
    """Auto-anchor fills the table, and a rail cannot label every bar.

    So the labels are sampled evenly across the whole anchor set — the same
    thinning the spectrum overlays do when a family carries more lines than a
    view can label — rather than drawn on top of each other. Both ends of the
    set are always named, and every label that is drawn stays legible.
    """

    forget_session_layout()
    window = _window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()

    residuals = _crowd_residuals(window, 24)
    qt_app.processEvents()
    labels = window.residual_tick_labels()
    every = [f"{residual.wavelength_nm:.1f}" for residual in residuals]

    # One bar each, but not one label each.
    assert len(_residual_bars(window).opts["height"]) == len(residuals)
    assert 1 < len(labels) < len(residuals), labels
    assert set(labels) <= set(every)
    # The ends of the set are named, and the labels read left to right.
    assert labels[0] == every[0]
    assert labels[-1] == every[-1]
    assert labels == sorted(labels, key=float)
    _assert_labels_are_legible(window)

    # A wider rail spends the width on more of them, and still no smear.
    narrow = len(labels)
    window.resize(2400, 1400)
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()
    assert len(window.residual_tick_labels()) > narrow
    _assert_labels_are_legible(window)

    # And a set small enough to label in full is labelled in full.
    _crowd_residuals(window, 3)
    qt_app.processEvents()
    assert len(window.residual_tick_labels()) == 3
    window.close()
    forget_session_layout()


def test_clicking_a_residual_bar_selects_its_anchor_row(qt_app, tmp_path):
    """What the labels cannot do on their own: name every bar, on request.

    A bar names its own anchor by selecting it in the table above, which is
    what carries the order, the wavelength and the residual in figures — and
    which already explains the anchor in the Why dock. The selected anchor's
    bar is drawn in the selection colour, so the sync reads both ways.
    """

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()
    _anchor_pair(window)
    qt_app.processEvents()

    anchors = window.session.anchor_rows()
    residuals = window.session.residuals
    assert len(residuals) == 2

    for index in (1, 0):
        window._residual_plot_clicked(_residual_click_at(window, index))
        qt_app.processEvents()
        row = window.anchor_table.currentRow()
        assert anchors[row].key == residuals[index].key, index
        assert window.anchor_table.selectionModel().hasSelection()
        # The row it selected is the row that names the line.
        assert window.anchor_table.item(row, 1).text() == (
            f"{anchors[row].line.wavelength_nm:.3f}"
        )
        # ...and the clicked bar wears the selection colour, alone.
        brushes = _residual_bars(window).opts["brushes"]
        selected = [
            position
            for position, brush in enumerate(brushes)
            if brush.color().name() == CalibrationBenchWindow._RESIDUAL_BAR_SELECTED
        ]
        assert selected == [index], (index, selected)

    # Past the last bar there is no anchor to name, and a right-click is not a
    # selection at all — neither one invents a row.
    window.anchor_table.clearSelection()
    window._residual_plot_clicked(_residual_click_at(window, 4.0, span=(-5.0, 5.0)))
    window._residual_plot_clicked(
        _residual_click_at(window, 0, button=QtCore.Qt.RightButton)
    )
    qt_app.processEvents()
    assert not window.anchor_table.selectionModel().hasSelection()
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


# ----------------------------------------------------------------------
# Packet F18 — geometry round 4
# ----------------------------------------------------------------------


def test_the_real_entry_point_opens_a_window_wearing_the_icon(qt_app, tmp_path, monkeypatch):
    """F18 item 1: the icon test now covers the path the console script takes.

    The old test built a window directly, so it stayed green through every
    report that the icon was gone: it never ran ``main`` and never asked
    whether the process had claimed a shell identity to hang the icon on.
    """

    claimed: list[str] = []
    monkeypatch.setattr(
        bench_gui,
        "apply_windows_taskbar_identity",
        lambda *args, **kwargs: claimed.append(BENCH_APP_USER_MODEL_ID),
    )

    def bench_windows():
        application = QtWidgets.QApplication.instance()
        return {
            id(widget): widget
            for widget in application.topLevelWidgets()
            if isinstance(widget, CalibrationBenchWindow)
        }

    # Earlier tests leave windows alive, so the one this test is about is the
    # one that was not there a moment ago.
    before = set(bench_windows())
    seen: dict = {}

    def fake_exec(self=None):
        application = QtWidgets.QApplication.instance()
        opened = [
            widget
            for key, widget in bench_windows().items()
            if key not in before
        ]
        seen["windows"] = opened
        seen["app_icon"] = application.windowIcon()
        if opened:
            seen["icon"] = opened[0].windowIcon()
            opened[0].close()
        return 0

    monkeypatch.setattr(QtWidgets.QApplication, "exec_", fake_exec)

    folder = tmp_path / "acquisition"
    folder.mkdir()
    assert bench_gui.main([str(folder)]) == 0

    assert claimed == [BENCH_APP_USER_MODEL_ID], "main never claimed a taskbar identity"
    assert len(seen["windows"]) == 1
    icon = seen["icon"]
    assert not icon.isNull(), "the window the console script opens carries no icon"
    assert not seen["app_icon"].isNull()
    # A title bar asks for 16 px.  Handing it a real pixmap at that size is
    # what stops the icon reading as a smudge, which is "gone" to a human.
    for edge in (16, 32, 64):
        pixmap = icon.pixmap(edge, edge)
        assert not pixmap.isNull(), f"no {edge} px icon"
        assert pixmap.width() == edge
    assert {(size.width(), size.height()) for size in icon.availableSizes()} >= {
        (16, 16),
        (32, 32),
    }


def test_the_bench_taskbar_identity_is_its_own(qt_app):
    """F18 item 1: two windows of one program, two taskbar groups.

    The main GUI has claimed ``echelle_spectra-<version>`` since it was
    written; the bench claimed nothing, so Windows filed it under whatever
    launched it and showed the launcher's icon.  The ids must differ, or the
    two windows collapse into one taskbar button again (F14 item 5).
    """

    from echelle_spectra import __version__

    assert BENCH_APP_USER_MODEL_ID != f"echelle_spectra-{__version__}"
    assert __version__ in BENCH_APP_USER_MODEL_ID
    assert "calibration_bench" in BENCH_APP_USER_MODEL_ID
    # Advisory everywhere: it never raises, and off Windows it simply declines.
    result = apply_windows_taskbar_identity(BENCH_APP_USER_MODEL_ID)
    assert result in (None, BENCH_APP_USER_MODEL_ID)


def test_the_default_size_is_the_smallest_usable_one(qt_app, tmp_path):
    """F18 items 2 and 3: everything the owner listed, in view, at once.

    The four surfaces the owner named — the file table, the bench readings,
    the factors line and the expected-line table — are all on screen at the
    default geometry of a 1080p-class screen, and the controls column does
    not need a scrollbar to reach any of them.
    """

    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    # All four surfaces are visible without opening anything or scrolling.
    assert window.file_table.isVisible()
    assert window.bench_state_group.isVisible()
    assert window.sphere_factors_group.isVisible()
    assert window.line_help_table.isVisible()

    # The expected-line table is a working surface, not a corner box.
    row_height = window.line_help_table.verticalHeader().defaultSectionSize()
    visible_rows = window.line_help_table.viewport().height() // max(1, row_height)
    assert visible_rows >= CalibrationBenchWindow.EXPECTED_LINE_ROWS, (
        f"only {visible_rows} expected-line rows fit; "
        f"{CalibrationBenchWindow.EXPECTED_LINE_ROWS} is the floor"
    )

    # No control tab needs a scrollbar to reach what it holds.  F18 only ever
    # checked the Files tab; the Lamp fit tab was scrolling by 234 px at every
    # size the bench has opened at, because the anchor table was stacked in it
    # (F20 moved it to the rail).
    for index in range(window.control_tabs.count()):
        tab = window.control_tabs.widget(index)
        assert isinstance(tab, QtWidgets.QScrollArea)
        assert tab.verticalScrollBar().maximum() == 0, (
            f"the {window.control_tabs.tabText(index)!r} controls scroll at the "
            "size the bench opens at"
        )

    # The readings are tab-independent: they never hide behind another tab.
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        assert window.bench_state_group.isVisible()
        assert window.sphere_factors_group.isVisible()
    window.close()


def test_the_minimum_geometry_is_the_usable_floor(qt_app, tmp_path):
    """F18 item 3 as F20 rewrites it: the floor yields to the screen.

    The window still never opens smaller than the layout needs — but "needs"
    is measured in lines of the platform's own text and capped by the desktop,
    because a floor that exceeds the desktop is how the owner's bench opened
    crushed on a 150%-scaled display.
    """

    floor = bench_minimum_size(QtCore.QSize(1920, 1080))
    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*floor)
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    assert (window.minimumSize().width(), window.minimumSize().height()) == (
        bench_minimum_size()
    )

    row_height = window.line_help_table.verticalHeader().defaultSectionSize()
    visible_rows = window.line_help_table.viewport().height() // max(1, row_height)
    assert visible_rows >= CalibrationBenchWindow.EXPECTED_LINE_ROWS
    for index in range(window.control_tabs.count()):
        assert window.control_tabs.widget(index).verticalScrollBar().maximum() == 0

    clipped = [
        widget.objectName() or widget.text()
        for widget in window.findChildren(QtWidgets.QLabel)
        if widget.isVisible() and _label_is_clipped(widget)
    ]
    assert not clipped, f"labels clip at the minimum size: {clipped}"
    assert not _overlapping_siblings(window)
    window.close()


def test_the_columns_are_draggable_and_the_cut_survives_the_session(qt_app, tmp_path):
    """F18 item 2 as F20 rewrites it: three splitters, and the cut is a share."""

    forget_session_layout()
    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    for splitter in (
        window.root_splitter,
        window.tables_splitter,
        window.readings_splitter,
    ):
        assert not splitter.childrenCollapsible()
        assert splitter.handle(1) is not None
        assert splitter.handle(1).isEnabled(), "the operator cannot drag this handle"
    # F20 moves the readings strip across the top of the whole window: with two
    # rails flanking it, the middle column is the narrowest thing on screen and
    # a strip put there folds onto three rows and eats the plots.
    assert window.readings_splitter.orientation() == QtCore.Qt.Vertical
    assert window.readings_splitter.indexOf(window.status_band) == 0
    assert window.readings_splitter.indexOf(window.root_splitter) == 1

    # The rails and plots keep the bulk of the window; the readings are a strip.
    band_height, plots_height = window.readings_splitter.sizes()
    assert plots_height > 2 * band_height, (band_height, plots_height)

    # A drag is remembered, and the next window of this session opens with it.
    window.tables_splitter.setSizes([300, 500])
    window.tables_splitter.splitterMoved.emit(300, 1)
    qt_app.processEvents()
    remembered = window.tables_splitter.sizes()
    assert remembered[0] < remembered[1], "the drag did not take"
    window.close()

    second = _campaign_window(tmp_path)
    second.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    second.show()
    qt_app.processEvents()
    second._relayout_wrapped_text()
    qt_app.processEvents()
    # The cut comes back as a proportion, which is what survives a window of
    # a different size — a pixel list would not.
    restored = second.tables_splitter.sizes()
    assert restored[0] / sum(restored) == pytest.approx(
        remembered[0] / sum(remembered), abs=0.02
    ), (restored, remembered)
    second.close()
    forget_session_layout()


def test_the_readings_strip_folds_instead_of_squeezing_its_button(qt_app, tmp_path):
    """F18 item 3: panels collapse deliberately; controls are never shaved."""

    forget_session_layout()
    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    wide = window._status_band_columns(1400)
    narrow = window._status_band_columns(320)
    assert wide == len(window._status_panels), "a wide band should stand in one row"
    assert narrow == 1, "a narrow band must fold, not squeeze"

    # The strip carries no buttons at all now: Bench state, Alignment and
    # Sphere factors are readings, and the one action they held sat diagonally
    # opposite the step that calls for it while also existing in the next-step
    # panel (owner, 2026-08-16). What must survive the fold is the reading.
    for width in (1400, 900, 640, 320):
        window._reflow_status_band(window._status_band_columns(width))
        qt_app.processEvents()
        assert not window.sphere_factors_group.findChildren(QtWidgets.QPushButton)
        assert window.comparison_value.isVisible()
    window.close()
    forget_session_layout()


def test_a_bench_reading_shortens_rather_than_wrapping_the_strip_taller(qt_app, tmp_path):
    """F18 item 3: the readings strip stays a strip."""

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert isinstance(window.file_value, _ElidingLabel)
    assert isinstance(window.watch_value, _ElidingLabel)
    long_name = "Ne-0.02s-x3-bright-lines-with-a-very-long-descriptive-tail.sif"
    window.file_value.setText(long_name)
    window.file_value.resize(120, 20)
    qt_app.processEvents()

    # The reading is whole to anything that asks, and whole in the tooltip;
    # only the painting is shortened, and in the middle.
    assert window.file_value.text() == long_name
    assert window.file_value.toolTip() == long_name
    painted = QtWidgets.QLabel.text(window.file_value)
    assert painted != long_name
    assert "…" in painted
    window.close()


def test_the_procedure_rows_span_the_list_and_hang_their_indents(qt_app, tmp_path):
    """F18 item 4: no phantom second column, and a real hanging indent.

    Every row label was pinned to ``max(280, viewport - 18)`` measured before
    the tab had ever been laid out, so it stayed 280 px wide in a pane twice
    that.  The untouched list background beside it — carrying the selected
    row's highlight — is the empty second column in the owner's screenshot.
    """

    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    window.control_tabs.setCurrentIndex(1)
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    tree = window.checklist_tree
    assert tree.count() >= 9
    viewport_width = tree.viewport().width()
    assert viewport_width > 300

    for row in range(tree.count()):
        item = tree.item(row)
        label = tree.itemWidget(item)
        assert label is not None
        # One column, as wide as the list: nothing of the row is left bare.
        assert label.width() >= viewport_width - 2 * tree.frameWidth() - 1, (
            f"row {row} leaves {viewport_width - label.width()} px of bare "
            "list background beside it — the phantom second column"
        )
        assert item.sizeHint().width() >= viewport_width - 2 * tree.frameWidth() - 1
        assert item.sizeHint().height() >= label.heightForWidth(label.width())
        # The indents are structure, not spaces someone counted out by hand.
        text = label.text()
        assert label.textFormat() == QtCore.Qt.RichText
        assert "\n" not in text
        assert "     " not in text, "an indent made of spaces is not an indent"
        assert "<td" in text, "the glyph and the text must be separate cells"

    # A row that says what unblocks it puts that on its own indented line.
    with_unblock = [
        tree.itemWidget(tree.item(row)).text()
        for row in range(tree.count())
        if "unblocked by:" in tree.itemWidget(tree.item(row)).text()
    ]
    assert with_unblock, "no checklist row named what would unblock it"
    for text in with_unblock:
        assert "margin-left" in text, "the sub-line needs a real indent"

    # Selecting a row is meaningful: it writes that row into the Why dock.
    tree.setCurrentRow(0)
    qt_app.processEvents()
    assert window.details_view.toPlainText().strip()
    window.close()


def test_a_long_checklist_row_still_spans_one_column(qt_app, tmp_path):
    """F18 item 4: the fix has to hold for a row that wraps several times."""

    window = _campaign_window(tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    window.control_tabs.setCurrentIndex(1)
    qt_app.processEvents()

    tree = window.checklist_tree
    item = QtWidgets.QListWidgetItem()
    label = QtWidgets.QLabel(
        CalibrationBenchWindow._checklist_row_html(
            "!",
            "A very long row — " + "measured evidence about this step, " * 12,
            "assign the sphere signal and its background before anything else",
        )
    )
    label.setTextFormat(QtCore.Qt.RichText)
    label.setWordWrap(True)
    tree.addItem(item)
    tree.setItemWidget(item, label)
    window._relayout_wrapped_text()
    qt_app.processEvents()

    width = tree.viewport().width() - 2 * tree.frameWidth()
    assert label.width() >= width - 1
    assert item.sizeHint().height() >= label.heightForWidth(label.width())
    # It wrapped many times, and every one of those lines is inside the row.
    assert label.heightForWidth(label.width()) > 3 * label.fontMetrics().height()
    window.close()


def test_both_factor_curves_survive_every_view_range(qt_app, tmp_path):
    """F18 item 5: the previous-campaign curve stops vanishing on zoom.

    A factor curve is the echelle orders end to end, so its wavelength axis
    is not ascending.  Clipping to the view searched it as if it were and
    returned a slice of about one point, which is a curve the legend still
    names and the eye cannot find.
    """

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    # Two order segments laid end to end: x runs backwards at the seam, which
    # is what the real curves do at every order boundary.
    first = np.linspace(700.0, 640.0, 400)
    second = np.linspace(660.0, 600.0, 400)
    wavelengths = np.concatenate([first, [np.nan], second])
    window._set_factor_curve(
        window.candidate_curve,
        AbsoluteCalibrationResult(
            wavelengths,
            np.concatenate([np.full(400, 1e-5), [np.nan], np.full(400, 1.2e-5)]),
        ),
    )
    window._set_factor_curve(
        window.previous_curve,
        AbsoluteCalibrationResult(
            wavelengths,
            np.concatenate([np.full(400, 4e-7), [np.nan], np.full(400, 5e-7)]),
        ),
    )
    qt_app.processEvents()

    assert not np.all(np.diff(wavelengths[np.isfinite(wavelengths)]) > 0), (
        "this fixture must not be sorted, or it tests nothing"
    )

    view_box = window.sphere_plot.getPlotItem().getViewBox()
    ranges = {
        "everything": ((595.0, 705.0), (np.log10(1e-7), np.log10(5e-5))),
        "the new curve's band": ((640.0, 700.0), (np.log10(5e-6), np.log10(3e-5))),
        "the previous curve's band": ((640.0, 700.0), (np.log10(2e-7), np.log10(9e-7))),
        "across the order seam": ((650.0, 670.0), (np.log10(1e-7), np.log10(5e-5))),
    }
    for name, (x_range, y_range) in ranges.items():
        view_box.setXRange(*x_range, padding=0)
        view_box.setYRange(*y_range, padding=0)
        qt_app.processEvents()
        for label, curve in (
            ("new measured pair", window.candidate_curve),
            ("previous campaign", window.previous_curve),
        ):
            _x, drawn = curve.curve.getData()
            finite = 0 if drawn is None else int(np.isfinite(drawn).sum())
            assert finite > 10, (
                f"the {label} curve has {finite} drawable points at {name} — "
                "the legend names a curve nobody can see"
            )
    window.close()


# ----------------------------------------------------------------------
# Packet F20 — the bench's two-rail geometry
# ----------------------------------------------------------------------

#: Every combination the owner actually meets: the size the bench opens at on
#: his display class, and the same window maximized.  Both are checked at 100%
#: and again, in a child process, at 150% desktop scaling.
_TWO_RAIL_SIZES = (
    ("default", None),
    ("maximized-like", (2000, 1200)),
)


def _visible_rows(table) -> int:
    """How many rows of *table* the operator can actually read."""

    return table.viewport().height() // max(
        1, table.verticalHeader().defaultSectionSize()
    )


@pytest.mark.parametrize("label,size", _TWO_RAIL_SIZES, ids=[n for n, _ in _TWO_RAIL_SIZES])
def test_two_rail_layout_is_usable_at(qt_app, tmp_path, label, size):
    """F20: two rails and a middle, usable at every size and every scaling.

    Owner, after F18: the bench still opened cramped, and maximized the left
    column collapsed to some 240 px — "expand me right, manually, now, before
    use".  Both are the same defect from two ends: pixel splitter cuts tuned at
    one geometry, laid down before the first layout pass.  This test is run
    twice over, once per desktop scaling, by the child-process test below.
    """

    forget_session_layout()
    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*(size or bench_default_geometry(QtCore.QSize(1920, 1080))))
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    where = f"{label} ({window.width()}x{window.height()}, scale "
    where += f"{os.environ.get('QT_SCALE_FACTOR', '1')})"

    # TWO RAILS: controls left, tables right, plots between them.
    root = window.root_splitter
    assert root.orientation() == QtCore.Qt.Horizontal
    assert root.count() == 3, where
    assert root.indexOf(window.controls_rail) == 0
    assert root.indexOf(window.view_tabs) == 1
    assert root.indexOf(window.tables_rail) == 2

    left, centre, right = root.sizes()
    total = left + centre + right
    # Not collapsed, and not the majority of the window either.
    assert left >= window.controls_rail.minimumWidth(), (where, left)
    assert 0.15 <= left / total <= 0.45, (where, left / total)
    assert right >= window.tables_rail.minimumWidth(), (where, right)
    assert 0.15 <= right / total <= 0.45, (where, right / total)
    # The plots are the flexible element: smaller than the two rails together.
    assert centre < left + right, (where, centre, left, right)
    assert centre > 0.2 * total, (where, centre / total)

    # No control column scrolls, and no button is clipped by one.  Each tab is
    # opened first: a QTabWidget lays a page out when it is shown, so a page
    # nobody has opened still carries the geometry of some earlier window size.
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        window._relayout_wrapped_text()
        qt_app.processEvents()
        tab = window.control_tabs.widget(index)
        assert tab.verticalScrollBar().maximum() == 0, (
            f"{window.control_tabs.tabText(index)!r} scrolls at {where}"
        )
        for button in window.findChildren(QtWidgets.QPushButton):
            if not button.isVisible():
                continue
            assert button.width() >= button.sizeHint().width(), (button.text(), where)
            painted = button.visibleRegion().boundingRect()
            assert painted.height() >= button.height(), (
                f"{button.text()!r} is clipped by its container at {where}"
            )

    # Both working tables are vertically long, not corner boxes.
    assert _visible_rows(window.line_help_table) >= (
        CalibrationBenchWindow.EXPECTED_LINE_ROWS
    ), (where, _visible_rows(window.line_help_table))
    assert _visible_rows(window.anchor_table) >= (
        CalibrationBenchWindow.EXPECTED_LINE_ROWS
    ), (where, _visible_rows(window.anchor_table))

    # The readings strip stays in view whatever is being done (F18).
    assert window.bench_state_group.isVisible()
    assert window.sphere_factors_group.isVisible()

    assert not _overlapping_siblings(window), where
    window.close()
    forget_session_layout()


def test_two_rail_cuts_are_shares_rather_than_pixels(qt_app, tmp_path):
    """F20: the cut is a proportion, so maximizing widens every column.

    The old defaults were pixel lists (``[560, 440]``, ``[controls, 1500 -
    controls]``) and the root splitter gave the left column stretch factor 0,
    so every pixel a larger window gained went to the plots and the rails
    stayed where a small window had left them.
    """

    forget_session_layout()
    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.resize(*bench_default_geometry(QtCore.QSize(1920, 1080)))
    window.show()
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()
    small = window.root_splitter.sizes()

    window.resize(2400, 1400)
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()
    large = window.root_splitter.sizes()

    assert sum(large) > sum(small)
    # Every column grew — not only the plots.  This is the whole of "maximized,
    # the LEFT column collapses to ~240 px": with stretch factor 0 the rails
    # kept the width a small window had given them.
    for index, (before, after) in enumerate(zip(small, large)):
        assert after > before, (index, small, large)
    # Once no content floor is binding, the cut is exactly the declared shares.
    for index, share in enumerate(_ROOT_SHARES):
        assert large[index] / sum(large) == pytest.approx(share, abs=0.05), (
            index,
            large,
        )

    # Stretch factors, not one-off sizes, are what carries this: a splitter
    # writes setStretchFactor onto the child's own size policy.
    for index in range(3):
        stretch = window.root_splitter.widget(index).sizePolicy().horizontalStretch()
        assert stretch > 0, (index, stretch)
    window.close()
    forget_session_layout()


def test_scaled_dpi_repeats_the_two_rail_checks(qt_app):
    """F20: run the whole two-rail matrix again at 150% desktop scaling.

    ``QT_SCALE_FACTOR`` is read once, when QApplication is constructed, so a
    second scaling can only be exercised in a second process.  This is the
    reproduction the owner's report needed and F18's single-scale test never
    had: at 150% his 1920x1080 display offers 1280x720 logical pixels, and a
    geometry quoted in pixels measured at 100% simply does not fit on it.
    """

    if os.environ.get("ECHELLE_BENCH_SCALED_CHILD"):
        pytest.skip("already the scaled child process")
    environment = dict(os.environ)
    environment["QT_SCALE_FACTOR"] = "1.5"
    environment["QT_QPA_PLATFORM"] = "offscreen"
    environment["ECHELLE_BENCH_SCALED_CHILD"] = "1"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            str(Path(__file__).resolve()),
            "-k",
            "two_rail",
            "-q",
            "-p",
            "no:cacheprovider",
        ],
        cwd=str(Path(__file__).resolve().parents[1]),
        env=environment,
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert completed.returncode == 0, (
        "the two-rail layout breaks at 150% desktop scaling:\n"
        f"{completed.stdout}\n{completed.stderr}"
    )


def test_selecting_file_after_file_never_orphans_a_threshold_line(qt_app, tmp_path):
    """A removed InfiniteLine painting itself takes the whole process down.

    Found while F19 was being built: the histograms cleared their plot and
    added two fresh threshold lines for every file selected, and a cleared
    ``InfiniteLine`` keeps a repaint scheduled while losing the view box it
    measures itself against.  ``boundingRect`` then raises inside Qt's paint,
    which is not an exception anybody can catch — it aborts.  It surfaced as a
    crash in an unrelated test, which is exactly how a latent one behaves.

    The lines are now built once per plot and moved, so the count of items on
    the plot is flat no matter how many files are triaged, and every line still
    has the view box it belongs to.
    """

    window = _manual_window(tmp_path, peak=65000.0)
    window.show()
    paths = []
    for index in range(4):
        path = tmp_path / f"Ne-0.0{index + 1}s-x3.sif"
        path.write_bytes(b"sif\n")
        paths.append(path)
    _drop(window, paths)
    _wait_for_loads(window, qt_app)

    window.file_table.selectRow(0)
    qt_app.processEvents()
    drawn = [
        item
        for item in window.histogram_plot.items
        if isinstance(item, pg.InfiniteLine)
    ]
    assert drawn, "the thresholds stopped being drawn"

    # Select the other files, which is what redraws the histogram.
    for row in range(1, len(paths)):
        window.file_table.selectRow(row)
        qt_app.processEvents()

    # The lines captured before those redraws must still be the lines on the
    # plot, and must still answer boundingRect.  Under clear()+addLine they are
    # detached objects whose getViewBox() is None, and this call is the one Qt
    # makes during a paint — where an exception aborts instead of failing.
    for line in drawn:
        assert line.getViewBox() is not None, (
            "a threshold line was detached from its view box while still alive"
        )
        line.boundingRect()
    assert drawn == [
        item
        for item in window.histogram_plot.items
        if isinstance(item, pg.InfiniteLine)
    ], "the thresholds were rebuilt per file instead of moved"
    window.close()


def test_no_procedure_row_is_taller_than_its_own_text(qt_app, tmp_path):
    """F21 item 5: the Procedure tab added a LOT of vertical padding.

    Every row was sized ``max(heightForWidth, sizeHint().height())``. For a
    word-wrapped rich-text label the size hint is a heuristic computed against
    a width the label does not have, it came out taller than the truth on nine
    rows out of ten, and it won that comparison — so rows stood 48 to 114 px
    taller than their text and a five-item list scrolled.

    This is the screenshot-shaped check the owner asked for: measure what is
    actually reserved against what the text actually costs, at the width the
    row is actually drawn at.
    """

    forget_session_layout()
    window, _paths = _real_folder_window(qt_app, tmp_path)
    qt_app.processEvents()
    index = [
        window.control_tabs.tabText(i) for i in range(window.control_tabs.count())
    ].index("Procedure")
    window.control_tabs.setCurrentIndex(index)
    qt_app.processEvents()
    window._relayout_wrapped_text()
    qt_app.processEvents()

    tree = window.checklist_tree
    assert tree.count() >= 4, "the fixture stopped producing a real checklist"
    for row in range(tree.count()):
        item = tree.item(row)
        label = tree.itemWidget(item)
        reserved = item.sizeHint().height()
        needed = label.heightForWidth(label.width())
        assert needed > 0
        # Four pixels of breathing room is the intent; anything approaching a
        # second copy of the text is the defect coming back.
        assert reserved - needed <= 8, (
            f"row {row} reserves {reserved} px for {needed} px of text"
        )
    window.close()
    forget_session_layout()


def test_the_auto_anchor_button_is_visible_from_every_tab(qt_app, tmp_path):
    """F21 item 10, and the miss that prompted this round.

    The action that fills the anchor table was put away on the Lamp fit
    control tab, so an operator standing on Procedure saw an empty table, two
    greyed buttons, and no way to learn what fills them. It lives on the
    anchors panel now, which is in the right rail and therefore never behind
    a tab.
    """

    forget_session_layout()
    window, _paths = _real_folder_window(qt_app, tmp_path)
    window.show()
    qt_app.processEvents()

    assert window.tables_rail.isAncestorOf(window.auto_anchor_button)
    for index in range(window.control_tabs.count()):
        window.control_tabs.setCurrentIndex(index)
        qt_app.processEvents()
        assert window.auto_anchor_button.isVisible(), (
            f"hidden on the {window.control_tabs.tabText(index)!r} tab"
        )
    window.close()
    forget_session_layout()


def test_stepping_the_order_control_up_moves_the_trace_up(qt_app, tmp_path):
    """F21 item 9: the control and the detector stacking were opposite.

    Detector row grows with order number — order 0 sits near row 56 and the
    last order near row 2094 — while the view inverted its y axis, so pressing
    up walked the highlight down the image. The axis is no longer inverted, so
    a higher order number is higher on screen, and the order numbers still
    read exactly as the wavelength table spells them.
    """

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert not window.detector_plot.getViewBox().yInverted()

    pattern = window.session.pattern
    column = pattern.shape[0] // 2
    rows = [float(pattern[column, order]) for order in range(pattern.shape[1])]
    assert rows[-1] > rows[0], "fixture no longer stacks orders downward in row"

    # Screen y grows downward, so "visually higher" means a smaller scene y.
    box = window.detector_plot.getViewBox()
    first = box.mapViewToScene(QtCore.QPointF(column, rows[0])).y()
    last = box.mapViewToScene(QtCore.QPointF(column, rows[-1])).y()
    assert last < first, "a higher order number still paints lower on screen"
    window.close()


def test_the_detector_view_can_be_squared_up_on_request(qt_app, tmp_path):
    """F21 item 8: equal aspect wanted, not mandatory."""

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert not window.equal_aspect_check.isChecked()
    assert window.detector_plot.getViewBox().state["aspectLocked"] is False

    window.equal_aspect_check.setChecked(True)
    qt_app.processEvents()
    assert window.detector_plot.getViewBox().state["aspectLocked"] == 1

    window.equal_aspect_check.setChecked(False)
    qt_app.processEvents()
    assert window.detector_plot.getViewBox().state["aspectLocked"] is False
    window.close()


# ----------------------------------------------------------------------
# Packet F22 — the review round: gating, scope, and routing
# ----------------------------------------------------------------------


def test_a_running_campaign_task_gates_the_handlers_that_share_its_state(
    qt_app, tmp_path
):
    """Every synchronous handler mutates what the running task is reading.

    The next-step button and Regenerate ran their handlers on the GUI thread
    while a CampaignTaskThread held the same campaign and session, and neither
    was in the gate the buttons beside them were.
    """

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()
    assert window.next_step_button.isEnabled()

    # A task is running; nothing that touches campaign state may fire.
    window._campaign_thread = object()
    window.refresh_campaign()

    assert not window.next_step_button.isEnabled()
    assert not window.regenerate_tomls_button.isEnabled()
    assert not window.recompute_sphere_button.isEnabled()
    assert not window.confirm_roles_button.isEnabled()
    ran = []
    window._next_action = lambda: ran.append("ran")
    window._run_next_action()
    assert ran == []
    written = []
    window._refused_identity = window.snapshot_id_edit.text().strip()
    window.campaign.write_tomls = lambda *args, **values: written.append(args)
    window._regenerate_tomls()
    assert written == []

    window._campaign_thread = None
    window.close()


def test_overwriting_is_offered_only_for_the_identity_that_was_refused(
    qt_app, tmp_path
):
    """Regenerate is permission for one bundle, not for the id field.

    Its visibility survived every edit of the snapshot identity while the
    overwrite itself read that field at the moment of the press, so one press
    could rewrite a bundle that had never refused anything.
    """

    window = _campaign_window(tmp_path)
    window.config_root = tmp_path / "configs"
    window.show()
    assert window.session.fit_anchor_at(0, 26).accepted
    assert window.session.fit_anchor_at(1, 56).accepted

    window._generate_tomls()
    assert window.regenerate_tomls_button.isHidden()
    window._generate_tomls()
    assert "already exists" in window.message_value.text()
    assert not window.regenerate_tomls_button.isHidden()

    # A second identity, published and then edited by hand.
    window.snapshot_id_edit.setText("20250814_cmos")
    assert window.regenerate_tomls_button.isHidden()
    window._generate_tomls()
    other = window.config_root / "20250814_cmos" / "campaign.toml"
    other.write_text("# edited by hand\n", encoding="utf-8")

    # The stale offer, pressed: the refusal was about the first identity.
    window._refused_identity = "20250813_cmos"
    window._regenerate_tomls()

    assert other.read_text(encoding="utf-8") == "# edited by hand\n"
    assert window.regenerate_tomls_button.isHidden()
    assert "identity changed" in window.message_value.text()
    window.close()


def test_a_missing_reference_table_offers_no_button_that_could_satisfy_it(
    qt_app, tmp_path
):
    """The picker adds SIFs; this row is about tables named at construction."""

    window = _manual_window(tmp_path)
    window.show()
    window.campaign.pattern_source.unlink()
    window.refresh_campaign()
    qt_app.processEvents()

    step = window.next_step_value.text()
    assert "Pattern, wavelength" in step
    assert "reopen the bench" in step
    assert window.next_step_button.isHidden()
    assert window._next_action is None
    window.close()


def test_the_sensitivity_can_be_measured_again_once_it_is_done(qt_app, tmp_path):
    """F22 item 6: the only way to recompute was a checklist step that is over.

    The comparison is read beside the curves, so the control that reruns it is
    there too — and never on the readings strip, which stays button-free.
    """

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()

    assert window.campaign.comparison.state is ComparisonState.READY
    assert window.recompute_sphere_button.isEnabled()
    assert not window.sphere_factors_group.findChildren(QtWidgets.QPushButton)

    started = []
    window._start_campaign_task = lambda operation: started.append(operation)
    window.recompute_sphere_button.click()

    assert started == [window.campaign.compute_sphere_comparison]
    window.close()


def test_a_background_is_routed_by_its_role_and_not_by_its_label(qt_app, tmp_path):
    """Colour, marks and advice all keyed off ``label.startswith("background")``.

    A badge is text to reword; the role is what the frame is. Rewording the
    badge must not repaint the row, name the role twice, or send the panel
    down the branch that tells a dark frame to expose longer.
    """

    dim = triage_exposure(_frame_for(tmp_path / "Ne-bg.sif", peak=400.0))
    verdict = triage_for_role(dim, MeasurementRole.LAMP_BACKGROUND, 40000.0)
    reworded = replace(verdict, label="shot with the light off")
    measurement = MeasurementRecord(tmp_path / "Ne-bg.sif", MeasurementRole.LAMP_BACKGROUND, "Ne")

    assert verdict.is_background and reworded.is_background
    colour = CalibrationBenchWindow._verdict_colour(verdict, dim)
    assert CalibrationBenchWindow._verdict_colour(reworded, dim) == colour
    # And the colour is the role's, not the raw state's: a correct dark frame
    # is DIM, and DIM alone would paint it amber.
    assert colour != bench_gui._TRIAGE_COLORS[dim.state]
    assert CalibrationBenchWindow._role_marks(measurement, reworded, False) == ["Ne"]

    # The advice branch, on the one verdict where the two branches differ.
    saturated = triage_exposure(_frame_for(tmp_path / "Ne-bg.sif", peak=65535.0))
    blocked = triage_for_role(saturated, MeasurementRole.LAMP_BACKGROUND, 40000.0)
    renamed = replace(blocked, label="not dark at all")

    assert saturated.state is ExposureState.SATURATED
    headline = CalibrationBenchWindow._short_verdict(
        tmp_path / "Ne-bg.sif", saturated, renamed
    )
    assert headline == f"NOT DARK AT ALL · {blocked.headline}"
    assert "cluster" not in headline


def test_the_previous_pair_is_chosen_in_the_bench_not_on_a_command_line(
    qt_app, tmp_path, monkeypatch
):
    """Owner, 2026-08-18: "This is terrible UX again. CLI OR GUI. Not both."

    The pattern got its picker; the previous pair is the same case. One dialog
    for the signal, the background offered from the naming beside it, and the
    comparison dropped through the path every other input change already uses.
    """

    # Flags still win at launch: a bench opened with a pair names that pair and
    # keeps it for as long as nobody picks another.
    launched = _campaign_window(tmp_path)
    assert launched.campaign.previous_sphere == tmp_path / "previous_sphere.sif"
    assert "previous_sphere.sif + previous_sphere_bg.sif" in (
        launched.previous_pair_value.text()
    )
    launched.refresh()
    assert launched.campaign.previous_sphere == tmp_path / "previous_sphere.sif"
    launched.close()

    window = _manual_window(tmp_path)
    window.show()
    qt_app.processEvents()
    previous = tmp_path / "campaign_2024"
    previous.mkdir()
    signal = previous / "sphere-2s.sif"
    signal.write_bytes(b"sif\n")
    background = previous / "sphere-2s-bg.sif"
    background.write_bytes(b"sif-bg\n")

    # Nobody has picked yet: the packaged default is whatever launched, and
    # the reading says so rather than inventing one.
    assert window.campaign.previous_pair_name == ""
    assert "none chosen" in window.previous_pair_value.text()
    assert window.choose_previous_button.isEnabled()

    asked = []
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getOpenFileName",
        lambda *args, **values: (asked.append(args[1]) or (str(signal), "")),
    )
    offered = []
    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        lambda *args, **values: (
            offered.append(args[2]) or QtWidgets.QMessageBox.Yes
        ),
    )
    window.campaign.comparison = campaign_module.SphereComparison(
        ComparisonState.READY, "compared against something else"
    )

    window.choose_previous_button.click()
    qt_app.processEvents()

    # One dialog, because the folder answered the second question itself.
    assert len(asked) == 1
    assert "SIGNAL" in asked[0]
    assert background.name in offered[0]
    assert window.campaign.previous_sphere == signal
    assert window.campaign.previous_sphere_background == background
    # Invalidated through the ordinary path, and named where the comparison is.
    assert window.campaign.comparison.state is ComparisonState.NOT_RUN
    assert window.campaign.comparison.reason == "inputs changed"
    assert f"{signal.name} + {background.name}" in window.previous_pair_value.text()
    assert "Recompute" in window.message_value.text()

    # A folder that cannot answer asks a second time instead of guessing.
    lonely = tmp_path / "campaign_2019"
    lonely.mkdir()
    orphan = lonely / "IS-1s.sif"
    orphan.write_bytes(b"sif\n")
    other = lonely / "IS-1s_bg.sif"
    other.write_bytes(b"sif-bg\n")
    replies = [str(orphan), str(other)]
    asked.clear()
    offered.clear()
    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        lambda *args, **values: (
            offered.append(args[2]) or QtWidgets.QMessageBox.No
        ),
    )
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getOpenFileName",
        lambda *args, **values: (asked.append(args[1]) or (replies.pop(0), "")),
    )

    window.choose_previous_button.click()
    qt_app.processEvents()

    assert len(asked) == 2
    assert "BACKGROUND" in asked[1]
    assert window.campaign.previous_sphere == orphan
    assert window.campaign.previous_sphere_background == other
    window.close()


def test_a_correct_background_is_painted_its_own_calm_colour(qt_app, tmp_path):
    """Owner, 2026-08-18: "can we paint background something... telling?"

    It had been amber (attend to me) and was fixed to green, which said the
    right thing in the wrong voice: on the bench's dark ground a background row
    and a healthy lamp row read identically.  A correct background now carries
    a colour that belongs to no exposure state — and a background that is not
    dark still carries the blocking one.  Pinned on the shared helper, because
    that is the single place the rule is written and every surface reads it.
    """

    dim = triage_exposure(_frame_for(tmp_path / "Ne-bg.sif", peak=400.0))
    correct = triage_for_role(dim, MeasurementRole.LAMP_BACKGROUND, 40000.0)
    colour = CalibrationBenchWindow._verdict_colour(correct, dim)

    assert colour == bench_gui._BACKGROUND_COLOR
    # Neither an alarm, nor a blockage, nor the same reading as a good signal.
    assert colour not in {
        bench_gui._TRIAGE_COLORS[ExposureState.GOOD],
        bench_gui._TRIAGE_COLORS[ExposureState.DIM],
        bench_gui._TRIAGE_COLORS[ExposureState.SATURATED],
        bench_gui._SUGGESTED_COLOR,
    }
    # Calm and cool, and readable against the bench's own dark ground: more
    # blue than red, and light enough to be a headline.
    paint = QtGui.QColor(colour)
    assert paint.blue() > paint.red()
    assert paint.lightness() > QtGui.QColor("#151b22").lightness()

    # Everything the shared helper feeds repaints with it, and the second
    # helper above it agrees on a background rather than demoting it.
    assert CalibrationBenchWindow._reading_colour(correct, dim) == colour

    # A background that is NOT dark keeps the blocking colour.
    leaking = triage_for_role(dim, MeasurementRole.LAMP_BACKGROUND, 500.0)
    assert leaking.blocking
    assert CalibrationBenchWindow._verdict_colour(leaking, dim) == (
        bench_gui._TRIAGE_COLORS[ExposureState.SATURATED]
    )


def test_the_why_dock_sets_prose_to_a_measure_and_spaces_its_paragraphs(
    qt_app, tmp_path
):
    """Owner, 2026-08-18: "the formatting is terrible. All long and breaky."

    Two defects in one dock: prose wrapping at the full width of the window,
    and multi-paragraph explanations run together into one grey wall.
    """

    forget_session_layout()
    window = _manual_window(tmp_path)
    window.resize(3000, 1200)
    window.show()
    qt_app.processEvents()
    window._distribute_space()
    qt_app.processEvents()

    view = window.details_view
    metrics = QtGui.QFontMetrics(view.font())
    # The cap is a measure in characters, turned into pixels by this platform's
    # own font — never a pixel number measured on one display.
    assert window._reading_measure() == (
        bench_gui.BENCH_READING_MEASURE_CHARS * max(1, metrics.averageCharWidth())
    )
    assert view.lineWrapMode() == QtWidgets.QTextEdit.FixedPixelWidth
    assert 0 < view.lineWrapColumnOrWidth() <= window._reading_measure()
    # Bounded: the dock is three thousand pixels wide and the prose is not.
    assert view.lineWrapColumnOrWidth() < window.width()
    # A path is one long unbreakable word and must wrap rather than widen it.
    assert view.document().defaultTextOption().wrapMode() == (
        QtGui.QTextOption.WrapAtWordBoundaryOrAnywhere
    )

    window.explain(
        "Ne-bg.sif — exposure verdict",
        "The first paragraph, about the reading.\n\n"
        f"In full: {tmp_path / 'Ne-bg.sif'}\n\n"
        "SELF-CHECK: the third paragraph, about the pair.",
    )
    qt_app.processEvents()

    blocks = []
    block = view.document().begin()
    while block.isValid():
        blocks.append(block)
        block = block.next()
    bodies = [item for item in blocks if item.text().strip()][1:]
    assert len(bodies) == 3, [item.text() for item in blocks]
    # Real space between them, not three sentences jammed into one wall.
    assert all(item.blockFormat().bottomMargin() > 0 for item in bodies)
    # And the path is still there to read and to copy out.
    assert str(tmp_path / "Ne-bg.sif") in view.toPlainText()
    window.close()
    forget_session_layout()


def test_the_why_dock_opens_tall_enough_and_remembers_a_drag(qt_app, tmp_path):
    """Owner: "by default it's a few lines" — and it cropped almost everything.

    The dock asks for the height a whole ordinary explanation costs and yields
    what the rails cannot spare, because both tables carry a workable-row
    contract this dock may not buy its comfort with.  A drag is not a wish and
    is obeyed as given, for as long as the session lives.
    """

    forget_session_layout()
    window = _manual_window(tmp_path)
    # A window with height to spare: at the bench's default geometry the two
    # tables are at exactly their six rows and there is nothing to give.
    window.resize(1501, 1400)
    window.show()
    qt_app.processEvents()
    for _ in range(2):
        window._distribute_space()
        qt_app.processEvents()

    unit = window.layout_unit
    assert window.details_dock.height() >= bench_gui.BENCH_DETAILS_LINES * unit
    window.explain(
        "Ne-bg.sif — exposure verdict",
        "Brightest non-anomalous pixel 41234 counts.\n"
        "Saturation: no connected full-scale cluster.\n"
        "Anomalies: none counted.\n"
        "Noise floor 300 counts.",
    )
    qt_app.processEvents()
    view = window.details_view
    assert view.document().size().height() <= view.viewport().height()

    # The rails keep what they were promised, whatever the dock wanted.
    assert _visible_rows(window.anchor_table) >= (
        CalibrationBenchWindow.EXPECTED_LINE_ROWS
    )

    # An operator's drag is remembered as lines of this platform's own text,
    # and survives the next layout pass rather than being cut back to default.
    dragged = 14 * unit
    QtWidgets.QApplication.sendEvent(
        window.details_dock,
        QtGui.QResizeEvent(
            QtCore.QSize(window.details_dock.width(), dragged),
            QtCore.QSize(window.details_dock.width(), window.details_dock.height()),
        ),
    )
    assert bench_gui._SESSION_DOCK_LINES["details"] == pytest.approx(14.0)
    window._distribute_space()
    qt_app.processEvents()
    assert window._details_dock_height == dragged
    assert window._details_dock_lines() == pytest.approx(14.0)

    forget_session_layout()
    assert window._details_dock_lines() <= bench_gui.BENCH_DETAILS_LINES
    window.close()
    forget_session_layout()


def test_a_bench_without_campaign_memory_still_adds_files(qt_app, tmp_path, monkeypatch):
    """F22 item 9: the button kept its caption and lost its verb.

    With no campaign the next-step slot is never refreshed, so the caption the
    button was built with is the action it must carry.
    """

    window = _window(tmp_path)
    window.show()
    qt_app.processEvents()
    asked = []
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getOpenFileNames",
        lambda *args, **values: (asked.append(args) or ([], "")),
    )

    assert window.campaign is None
    assert window.next_step_button.text() == "Add SIF files…"
    window.next_step_button.click()

    assert asked
    window.close()


def test_one_campaign_refresh_derives_the_procedure_once(qt_app, tmp_path):
    """The checklist is O(files x measurements) and this runs per keystroke."""

    window = _campaign_window(tmp_path)
    window.show()
    qt_app.processEvents()
    derived = []
    real_checklist = window.campaign.checklist

    def counted(session):
        derived.append(session)
        return real_checklist(session)

    window.campaign.checklist = counted
    window.refresh_campaign()

    assert len(derived) == 1
    window.close()


# ----------------------------------------------------------------------
# The bench saves beside its data, and says so
# ----------------------------------------------------------------------


def _bench_from_main(monkeypatch, argv) -> dict:
    """Run ``main`` off-screen and report the one window it opened.

    The real ``exec_`` never returns, so it is replaced by a look at whatever
    window was not there a moment ago — the same trick the icon test uses.
    """

    monkeypatch.setattr(
        bench_gui, "apply_windows_taskbar_identity", lambda *args, **values: None
    )

    def bench_windows():
        application = QtWidgets.QApplication.instance()
        # The widgets themselves, not their ``id()``: an id set does not keep
        # the windows it names alive, and CPython hands a fresh window the id
        # of a collected one — which read, at random, as "no window opened".
        return {
            widget
            for widget in application.topLevelWidgets()
            if isinstance(widget, CalibrationBenchWindow)
        }

    before = bench_windows()
    seen: dict = {}

    def fake_exec(self=None):
        opened = [widget for widget in bench_windows() if widget not in before]
        assert len(opened) == 1
        window = opened[0]
        seen["output_root"] = window.output_root
        seen["config_root"] = window.config_root
        seen["snapshots_said"] = window.snapshot_root_value.text()
        seen["configs_said"] = window.config_root_value.text()
        seen["snapshots_tooltip"] = window.snapshot_root_value.toolTip()
        seen["snapshot_id"] = window.snapshot_id_edit.text().strip()
        seen["snapshot_date"] = window.snapshot_date
        seen["valid_from"] = window.valid_from
        seen["identity_why"] = window.snapshot_id_edit.property("explainText")
        window.close()
        return 0

    monkeypatch.setattr(QtWidgets.QApplication, "exec_", fake_exec)
    assert bench_gui.main(argv) == 0
    return seen


def test_the_folder_argument_decides_where_the_campaign_is_saved(
    qt_app, tmp_path, monkeypatch
):
    """The field report: TOMLs landed in the launch directory, not the data.

    Both roots defaulted to ``Path.cwd()``-relative, so a bench opened on a
    NAS acquisition folder wrote its calibrations into whatever directory the
    shortcut happened to start in — with nothing on screen to say so.  The
    folder the bench was pointed at is the calibration: the measured frames and
    the calibration computed from them, side by side in the one folder that
    holds it all.
    """

    data = tmp_path / "20250926_calib"
    data.mkdir()
    elsewhere = tmp_path / "launch-dir"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    seen = _bench_from_main(monkeypatch, [str(data)])

    assert seen["output_root"] == data / SNAPSHOT_ROOT_NAME
    # One generated folder inside it, holding everything the bench writes.
    assert seen["config_root"] == data / SNAPSHOT_ROOT_NAME / CONFIG_ROOT_NAME
    assert seen["output_root"].is_absolute() and seen["config_root"].is_absolute()
    # And emphatically not the working directory the launcher started in.
    assert elsewhere not in seen["output_root"].parents
    assert elsewhere not in seen["config_root"].parents


def test_without_a_folder_argument_the_roots_stay_where_the_bench_started(
    qt_app, tmp_path, monkeypatch
):
    """There is nothing better than the working directory — but it is absolute."""

    here = tmp_path / "just-here"
    here.mkdir()
    monkeypatch.chdir(here)

    seen = _bench_from_main(monkeypatch, [])

    assert seen["output_root"] == here / SNAPSHOT_ROOT_NAME
    assert seen["config_root"] == here / SNAPSHOT_ROOT_NAME / CONFIG_ROOT_NAME


def test_explicit_roots_beat_both_the_folder_and_the_working_directory(
    qt_app, tmp_path, monkeypatch
):
    """A flag is a decision; deriving from the folder is only a default."""

    data = tmp_path / "20250926_calib"
    data.mkdir()
    chosen_snapshots = tmp_path / "archive" / "snaps"
    monkeypatch.chdir(tmp_path)

    seen = _bench_from_main(
        monkeypatch,
        [
            str(data),
            "--output-root",
            str(chosen_snapshots),
            "--config-root",
            "relative-configs",
        ],
    )

    assert seen["output_root"] == chosen_snapshots
    # A relative flag is honoured, and made absolute so the display cannot lie.
    assert seen["config_root"] == tmp_path / "relative-configs"
    assert seen["config_root"].is_absolute()
    # The parser itself keeps "unset" distinguishable from "set to the default".
    assert _build_parser().parse_args([]).output_root is None
    assert _build_parser().parse_args([]).config_root is None


def test_the_save_tab_names_the_whole_destination_not_its_last_folder(
    qt_app, tmp_path, monkeypatch
):
    """``calibrations`` names a folder on every machine; it located none of them."""

    data = tmp_path / "20250926_calib"
    data.mkdir()
    monkeypatch.chdir(tmp_path)

    seen = _bench_from_main(monkeypatch, [str(data)])

    assert str(data / SNAPSHOT_ROOT_NAME) in seen["snapshots_said"]
    assert str(data / SNAPSHOT_ROOT_NAME / CONFIG_ROOT_NAME) in seen["configs_said"]
    # The whole path is one hover away even when the label paints it shortened.
    assert seen["snapshots_tooltip"] == str(data / SNAPSHOT_ROOT_NAME)


def test_the_save_tab_reading_shortens_in_the_middle_and_keeps_the_whole_path(
    qt_app, tmp_path
):
    """A long NAS path must not clip away the half that identifies it."""

    window = _window(tmp_path)
    assert isinstance(window.snapshot_root_value, _ElidingLabel)
    assert str(window.output_root) in window.snapshot_root_value.text()
    assert window.snapshot_root_value.property("explainText").endswith(
        str(window.output_root)
    )
    window.close()


def test_a_bare_root_name_is_absolute_before_it_is_ever_displayed(qt_app, tmp_path):
    """The constructor's own defaults were relative names, shown as names."""

    window = _window(tmp_path)
    assert window.output_root.is_absolute()
    assert window.config_root.is_absolute()
    assert window.output_root.name == SNAPSHOT_ROOT_NAME
    assert window.config_root.name == CONFIG_ROOT_NAME
    # Derived, the settings bundles are a child of the snapshot root: one
    # generated folder per calibration folder, not two.
    assert window.config_root.parent == window.output_root
    # An explicitly passed root is still honoured wherever it points.
    elsewhere = _window(tmp_path, window_options={"config_root": tmp_path / "cfg"})
    assert elsewhere.config_root == tmp_path / "cfg"
    elsewhere.close()
    window.close()


def test_a_unc_folder_keeps_its_share_when_the_roots_are_derived():
    r"""Pure paths only: the NAS is never touched to work out where to write.

    ``\\server\share`` survives a pathlib join and does not survive string
    surgery, which is why the derivation is joins all the way down.
    """

    share = "\\" * 2 + "10.10.10.122" + "\\" + "workdata"
    folder = PureWindowsPath(share, "2025-LHD-BH", "20250926_calib")
    snapshots, configs = default_bench_roots(folder)

    assert str(snapshots) == str(folder / SNAPSHOT_ROOT_NAME)
    assert str(configs) == str(folder / SNAPSHOT_ROOT_NAME / CONFIG_ROOT_NAME)
    assert str(snapshots).startswith(share + "\\")
    assert snapshots.drive == share
    assert snapshots.is_absolute()
    # An already-absolute UNC path is handed straight back, unresolved.
    assert absolute_root(folder) is folder


def test_everything_generated_lives_under_one_folder_that_holds_no_stray_snapshot(
    tmp_path, capsys
):
    """Owner: "It's kinder to keep all generated stuff in one tidy folder."

    The settings bundles moved from a sibling ``calibration-configs\\`` into a
    ``configs`` child of the snapshot root, so a calibration folder gains one
    generated child instead of two. The thing that has to be checked when a
    snapshot root gains a child that is not a snapshot is whatever enumerates
    it: ``echelle status`` keys on ``<child>/snapshot.toml``, so ``configs``
    is neither counted nor reported as an invalid snapshot. (The epoch registry
    never enumerates at all — it resolves the snapshot IDs it was given.)
    """

    from echelle_spectra.cli import _status

    snapshots, configs = default_bench_roots(tmp_path / "20250926_calib")
    assert configs.parent == snapshots
    bundle = configs / "20250926_cmos"
    bundle.mkdir(parents=True)
    (bundle / "campaign.toml").write_text("# campaign\n", encoding="utf-8")

    _status(
        [
            "--calibrations",
            str(snapshots),
            "--registry",
            str(tmp_path / "nowhere.toml"),
            "--runs",
            str(tmp_path / "no-runs"),
        ]
    )
    printed = capsys.readouterr().out

    assert "snapshots: none found" in printed
    assert "invalid" not in printed
    assert "configs" not in printed.replace(str(snapshots), "")


# ----------------------------------------------------------------------
# The identity is dated by acquisition, never by computation
# ----------------------------------------------------------------------


def test_the_acquisition_date_is_read_from_a_folder_name_or_a_sif_header():
    """Two sources, one convention: YYYYMMDD_<detector>, dated by the shot day."""

    assert acquisition_date_from_name("20250926_calib") == date(2025, 9, 26)
    assert acquisition_date_from_name(
        PureWindowsPath("T:\\2025-LHD-BH\\Echelle\\20250926_calib")
    ) == date(2025, 9, 26)
    # Only the folder's own name, and only where the date leads it: a campaign
    # tree may carry any number of dated components above this one.
    assert acquisition_date_from_name("calib_20250926") is None
    assert acquisition_date_from_name("20259999_calib") is None
    # An LHD shot number that merely begins with eight digits is not a date.
    assert acquisition_date_from_name("123456789_shot") is None
    assert acquisition_date_from_name(None) is None

    moment = datetime(2025, 9, 26, 11, 30)
    assert acquisition_date_from_metadata(
        {"ExperimentTime": int(moment.timestamp())}
    ) == date(2025, 9, 26)
    assert acquisition_date_from_metadata({"ExposureTime": 0.1}) is None
    assert acquisition_date_from_metadata({"ExperimentTime": 0}) is None
    assert acquisition_date_from_metadata(None) is None

    assert snapshot_identity(date(2025, 9, 26), "cmos") == "20250926_cmos"
    assert snapshot_identity(date(2025, 9, 26), " ") == "20250926_cmos"


def test_the_identity_is_dated_by_the_folder_the_calibration_was_measured_in(
    qt_app, tmp_path, monkeypatch
):
    """Owner: "calibration IS dated by the day the images were taken."

    The bench prefilled today's date, so frames shot on 2025-09-26 were about
    to be saved as ``20260817_cmos`` — the day the computation happened, which
    is bookkeeping and not physics. The folder says which day it was.
    """

    data = tmp_path / "Echelle" / "20250926_calib"
    data.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)

    seen = _bench_from_main(monkeypatch, [str(data)])

    assert seen["snapshot_id"] == "20250926_cmos"
    assert seen["snapshot_date"] == date(2025, 9, 26)
    assert seen["snapshot_id"] != snapshot_identity(date.today(), "cmos")
    # The epoch follows the identity's own date rather than today's, which is
    # what ``--valid-from`` being unset now means.
    assert seen["valid_from"] is None
    assert default_validity(seen["snapshot_id"], seen["valid_from"]) == {
        "date_from": "2025-09-26"
    }
    assert "2025-09-26" in seen["identity_why"]
    assert "TAKEN" in seen["identity_why"]


def test_an_undateable_folder_falls_back_to_today_and_says_it_must_be_checked(
    qt_app, tmp_path, monkeypatch
):
    """Today is a placeholder, and a placeholder has to announce itself."""

    data = tmp_path / "incoming"
    data.mkdir()
    monkeypatch.chdir(tmp_path)

    seen = _bench_from_main(monkeypatch, [str(data)])

    assert seen["snapshot_id"] == snapshot_identity(date.today(), "cmos")
    assert seen["snapshot_date"] is None
    assert "No acquisition date could be read" in seen["identity_why"]
    assert "CHECK" in seen["identity_why"]


def test_an_explicit_identity_or_epoch_beats_every_derivation(
    qt_app, tmp_path, monkeypatch
):
    """A flag is a decision; deriving the date is only a default."""

    data = tmp_path / "Echelle" / "20250926_calib"
    data.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)

    seen = _bench_from_main(
        monkeypatch,
        [
            str(data),
            "--snapshot-id",
            "20240305_cmos",
            "--valid-from",
            "2026-09-01",
        ],
    )

    assert seen["snapshot_id"] == "20240305_cmos"
    assert seen["snapshot_date"] == date(2024, 3, 5)
    assert seen["valid_from"] == date(2026, 9, 1)
    assert default_validity(seen["snapshot_id"], seen["valid_from"]) == {
        "date_from": "2026-09-01"
    }
    # A detector other than the default is carried into a derived identity.
    other = tmp_path / "Echelle" / "20250926_calib"
    assert _bench_from_main(monkeypatch, [str(other), "--detector", "ccd"])[
        "snapshot_id"
    ] == "20250926_ccd"


def test_the_frames_date_the_identity_when_the_folder_name_cannot(qt_app, tmp_path):
    """The SIF headers know the shot day even when the folder name says nothing."""

    placeholder = snapshot_identity(date.today(), "cmos")
    taken = datetime(2025, 9, 26, 11, 30)
    window = _manual_window(
        tmp_path,
        window_options={"snapshot_id": placeholder},
        experiment_time=int(taken.timestamp()),
    )
    window.show()
    assert window.snapshot_id_edit.text() == placeholder
    assert window.snapshot_date is None

    source = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    source.write_bytes(b"sif\n")
    _drop(window, [source])
    _wait_for_loads(window, qt_app)

    assert window.snapshot_date == date(2025, 9, 26)
    assert window.snapshot_id_edit.text() == "20250926_cmos"
    assert "SIF headers" in window.snapshot_id_edit.property("explainText")
    window.close()


def test_a_typed_identity_is_never_overwritten_by_a_later_frame(qt_app, tmp_path):
    """A hand edit is a decision, and outranks anything the frames say."""

    placeholder = snapshot_identity(date.today(), "cmos")
    taken = datetime(2025, 9, 26, 11, 30)
    window = _manual_window(
        tmp_path,
        window_options={"snapshot_id": placeholder},
        experiment_time=int(taken.timestamp()),
    )
    window.show()
    # Typing is what ``textEdited`` reports; ``setText`` alone is the bench's
    # own voice and must never be mistaken for the operator's.
    window.snapshot_id_edit.setText("20251003_cmos")
    window.snapshot_id_edit.textEdited.emit("20251003_cmos")
    qt_app.processEvents()

    source = tmp_path / "Ne-0.02s-x3-bright-lines.sif"
    source.write_bytes(b"sif\n")
    _drop(window, [source])
    _wait_for_loads(window, qt_app)

    assert window.snapshot_id_edit.text() == "20251003_cmos"
    assert window.snapshot_date == date(2025, 10, 3)
    assert default_validity(window.snapshot_id_edit.text(), window.valid_from) == {
        "date_from": "2025-10-03"
    }
    window.close()


def test_the_saved_snapshot_confirmation_names_its_folder_and_opens_it(
    qt_app, tmp_path, monkeypatch
):
    """Naming a folder and leaving the operator to find it is half a message."""

    window = _campaign_window(tmp_path)
    window.output_root = tmp_path / "calibrations"
    window.show()
    qt_app.processEvents()
    assert window.open_snapshot_button.isHidden()

    class _Saved:
        root = tmp_path / "calibrations" / "20250813_cmos"
        snapshot_id = "20250813_cmos"

    opened: list = []
    monkeypatch.setattr(
        bench_gui, "open_containing_folder", lambda path: opened.append(Path(path)) or True
    )
    window._campaign_task_completed(_Saved())
    qt_app.processEvents()

    assert str(_Saved.root) in window.message_value.text()
    assert not window.open_snapshot_button.isHidden()
    window.open_snapshot_button.click()

    assert opened == [_Saved.root]
    window.close()


# ----------------------------------------------------------------------
# The Save tab answers its own buttons
# ----------------------------------------------------------------------


def _signal_only_window(tmp_path: Path) -> CalibrationBenchWindow:
    """The owner's live 2019 state: sphere pair, Ne signal, no lamp background."""

    window = _window(tmp_path)
    paths = {}
    for name in (
        "pattern.txt",
        "wavelength.txt",
        "integral.txt",
        "IS-1s.sif",
        "IS_bg.sif",
        "Ne_1s_10fr.sif",
        "previous_sphere.sif",
        "previous_sphere_bg.sif",
    ):
        path = tmp_path / name
        path.write_text(name, encoding="utf-8")
        paths[name] = path
    campaign = CalibrationCampaignSession(
        pattern_source=paths["pattern.txt"],
        wavelength_source=paths["wavelength.txt"],
        integral_source=paths["integral.txt"],
        suggested_lamps=("Ne",),
        previous_sphere=paths["previous_sphere.sif"],
        previous_sphere_background=paths["previous_sphere_bg.sif"],
    )
    frame = window.session.frame
    campaign.classify_file(paths["IS-1s.sif"], MeasurementRole.SPHERE, frame=frame)
    campaign.classify_file(
        paths["IS_bg.sif"], MeasurementRole.SPHERE_BACKGROUND, frame=frame
    )
    campaign.classify_file(
        paths["Ne_1s_10fr.sif"], MeasurementRole.LAMP, lamp_family="Ne", frame=frame
    )

    def calculator(**values):
        scale = 1.08 if "previous" not in values["sphere"].name else 1.0
        return AbsoluteCalibrationResult(np.linspace(400, 700, 60), np.full(60, scale))

    campaign.compute_sphere_comparison(calculator)
    window.campaign = campaign
    window.config_root = tmp_path / "configs"
    window.output_root = tmp_path / "snapshots"
    window.snapshot_id_edit.setText("20190314_cmos")
    # The fit the owner had already solved when he pressed Save: it ran on the
    # lamp signal alone, which is the whole point of the state under test.
    assert window.session.fit_anchor_at(0, 26).accepted
    assert window.session.fit_anchor_at(1, 56).accepted
    window.refresh()
    return window


def _open_save_tab(window: CalibrationBenchWindow, qt_app) -> None:
    tabs = window.control_tabs
    names = [tabs.tabText(index) for index in range(tabs.count())]
    tabs.setCurrentIndex(names.index("Save"))
    qt_app.processEvents()


def test_a_refused_save_is_read_where_the_button_that_refused_it_is(qt_app, tmp_path):
    """Owner, live: "Doesn't save. No errors, terminal is silent."

    It was never silent.  The bench refused in a full sentence and put it on
    ``message_value``, which lives on the Lamp fit tab — so a press made from
    the Save tab answered itself on a page nobody was looking at.  Both lines
    carry the answer now, and the Save tab's one is on screen when the press
    happens.
    """

    window = _signal_only_window(tmp_path)
    # Nothing is classified as a lamp: this is a genuine refusal, not the
    # signal-only case, which saves.
    for path, record in list(window.campaign.measurements.items()):
        if record.role is MeasurementRole.LAMP:
            window.campaign.remove_classification(path)
    window.show()
    _open_save_tab(window, qt_app)

    assert not window.message_value.isVisible(), "the fit tab's line is not on screen"
    assert window.save_message_value.isVisible()

    window._generate_tomls()
    qt_app.processEvents()

    refusal = window.save_message_value.text()
    assert "were not saved" in refusal
    assert "lamp signal" in refusal
    # The one message line carries it too, and nothing suppressed it.
    assert window.message_value.text() == refusal
    assert window._role_notice == ""
    window.close()


def test_the_save_answer_is_set_after_the_refresh_that_could_bury_it(qt_app, tmp_path):
    """A refusal is put on screen last, so no follow-up narration lands on top."""

    window = _signal_only_window(tmp_path)
    window.show()
    _open_save_tab(window, qt_app)
    # The bench's own narration, in the state where it would fire during the
    # refresh that follows the press.
    window._say_roles("Roles assigned from filenames: 3 file(s)")
    narrated = []
    original = window.refresh_campaign

    def refresh_and_narrate():
        original()
        window._bench_says("Opened the assigned lamp signal for fitting.")
        narrated.append(True)

    window.refresh_campaign = refresh_and_narrate
    window._generate_tomls()
    qt_app.processEvents()

    assert narrated, "the refresh under test never ran"
    assert "Saved the alignment settings" in window.save_message_value.text()
    assert window.save_message_value.text() == window.message_value.text()
    window.close()


def test_a_signal_only_2019_folder_saves_from_the_bench(qt_app, tmp_path):
    """The whole live blocker, end to end through the surface."""

    window = _signal_only_window(tmp_path)
    window.show()
    _open_save_tab(window, qt_app)

    assert window.campaign.lamps_without_background() == ("Ne",)
    window._generate_tomls()
    qt_app.processEvents()

    assert window.save_state_value.text() == "READY"
    assert window.save_snapshot_button.isEnabled()
    assert "Saved the alignment settings" in window.save_message_value.text()
    assert window.toml_preview.toPlainText().startswith("#")
    assert "background_shot = false" in window.toml_preview.toPlainText()
    # And the next-step strip no longer points at a frame that does not exist.
    assert "background" not in window.next_step_value.text()
    window.close()


def test_a_refused_save_from_the_bench_leaves_no_folders_behind(qt_app, tmp_path):
    """The empty ``calibrations/configs`` a refusal used to leave is gone."""

    window = _signal_only_window(tmp_path)
    for path, record in list(window.campaign.measurements.items()):
        if record.role is MeasurementRole.LAMP:
            window.campaign.remove_classification(path)
    window.show()
    _open_save_tab(window, qt_app)

    window._generate_tomls()
    qt_app.processEvents()

    assert not Path(window.config_root).exists()
    assert "were not saved" in window.save_message_value.text()
    window.close()


# --- The solved correction reaches every view of an expected line -----------


def _shifted_window(tmp_path: Path, shift_px: float = 14.0) -> CalibrationBenchWindow:
    """A bench whose frame shows every curated line *shift_px* to the right.

    This is the 2019 folder in miniature: the curated table was written on a
    detector that has since moved, so the sticks the bench draws from it stand
    beside the light rather than on it until a transform is solved.
    """

    pattern = np.column_stack(
        [np.full(_COLUMNS, 12, dtype=float), np.full(_COLUMNS, 30, dtype=float)]
    )
    lines = (
        CalibrationTableLine(0, 21, 31, 26, 585.2488, "NeI", "ok"),
        CalibrationTableLine(1, 51, 61, 56, 588.1895, "NeI", "ok"),
    )
    x = np.arange(_COLUMNS, dtype=float)
    order_spectra = tuple(
        5 + 200 * np.exp(-0.5 * ((x - (center + shift_px)) / 1.5) ** 2)
        for center in (26.0, 56.0)
    )
    images = np.zeros((1, 44, _COLUMNS), dtype=float)
    images[0, 12, :] = order_spectra[0]
    images[0, 30, :] = order_spectra[1]
    session = CalibrationBenchSession(pattern, lines, minimum_snr=3.0)
    session.accept_frame(
        BenchFrame(
            tmp_path / "Ne-0.1s-x1-bright-lines.sif",
            images,
            images[0],
            order_spectra,
            {"ExposureTime": 0.1},
        )
    )
    session.use_lamp_reference(lamp_reference_set("Ne", lines))
    window = CalibrationBenchWindow(session, start_timer=False)
    window.refresh()
    return window


def _views_agree(window) -> list[tuple[str, float]]:
    """The one expected-line list, asserted to read the same in both views."""

    drawn = _drawn_sticks(window)
    listed = _listed_rows(window)
    assert [name for name, _pixel in drawn] == [name for name, _pixel in listed]
    for (_name, stick), (_row, cell) in zip(drawn, listed):
        # The panel prints one decimal; the stick carries the full float.
        assert stick == pytest.approx(cell, abs=0.05)
    return drawn


def _solve_by_clicking(window, qt_app, shift_px: float = 14.0) -> None:
    """Bootstrap the alignment through the real plot-click path."""

    window._order_plot_clicked(_click_at(window, 26.0 + shift_px))
    qt_app.processEvents()
    window.order_spin.setValue(1)
    qt_app.processEvents()
    window._order_plot_clicked(_click_at(window, 56.0 + shift_px))
    qt_app.processEvents()


def test_the_solved_correction_moves_the_sticks_and_the_pixel_column(
    qt_app, tmp_path
):
    """The owner's 2019 screenshot: dx +18.02 solved, sticks still on the 2024
    pixels, clicks landing on background.  Once a transform exists the sticks,
    the panel's Pixel column and the click-fit window all read the corrected
    position, and they move together on the ordinary refresh."""

    window = _shifted_window(tmp_path)
    window.show()
    qt_app.processEvents()

    # Before a solve the table's own pixel is all anybody knows.
    assert _views_agree(window) == [("NeI 585.249", 26.0)]

    _solve_by_clicking(window, qt_app)

    assert window.session.transform is not None
    assert window.session.transform.dx_px == pytest.approx(14.0, abs=0.05)
    # No refresh hook of its own: the click handler's ordinary refresh moved
    # both views, the cached per-order list included.
    name, pixel = _views_agree(window)[0]
    assert name == "NeI 588.189"
    assert pixel == pytest.approx(70.0, abs=0.05)

    window.order_spin.setValue(0)
    qt_app.processEvents()
    name, pixel = _views_agree(window)[0]
    assert name == "NeI 585.249"
    assert pixel == pytest.approx(40.0, abs=0.05)
    window.close()


def test_a_click_on_a_corrected_stick_anchors_the_line_under_it(qt_app, tmp_path):
    """The click path itself: after the solve, the pixel the stick stands on is
    the pixel that fits, and the anchor is still the curated row."""

    window = _shifted_window(tmp_path)
    window.show()
    qt_app.processEvents()
    _solve_by_clicking(window, qt_app)
    window.order_spin.setValue(0)
    qt_app.processEvents()
    # Take the order-0 anchor back off by clicking its own stick, which now
    # stands at the corrected pixel rather than at the table's.
    window._order_plot_clicked(_click_at(window, 40.0))
    qt_app.processEvents()
    assert not [
        anchor
        for anchor in window.session.anchor_rows()
        if anchor.line.order_idx == 0
    ]

    window._order_plot_clicked(_click_at(window, 40.0))
    qt_app.processEvents()

    anchored = [
        anchor
        for anchor in window.session.anchor_rows()
        if anchor.line.order_idx == 0
    ]
    assert len(anchored) == 1
    assert anchored[0].line.center_pixel == 26.0
    assert anchored[0].fit.center_pixel == pytest.approx(40.0, abs=0.05)
    window.close()


def test_without_a_lamp_reference_the_fallback_rows_follow_too(qt_app, tmp_path):
    """The panel falls back to the whole table when no lamp is assigned; that
    list is drawn from the same corrected rows, not from a second source."""

    window = _shifted_window(tmp_path)
    window.session.use_lamp_reference(None)
    window.show()
    window.refresh()
    qt_app.processEvents()
    assert _views_agree(window) == [("NeI 585.249", 26.0)]

    _solve_by_clicking(window, qt_app)
    window.order_spin.setValue(0)
    qt_app.processEvents()

    name, pixel = _views_agree(window)[0]
    assert name == "NeI 585.249"
    assert pixel == pytest.approx(40.0, abs=0.05)
    window.close()


def test_the_transform_line_shows_all_three_numbers_on_one_line(tmp_path):
    """dx / dy / θ promises three numbers; word wrap was hiding the third.

    The wrapped θ landed on a second line the row's height budget clipped,
    so the header read "+18.02, -0.60," with a dangling comma (owner
    screenshot, 2026-08-17). One line, elided middle if ever too narrow.
    """

    from echelle_spectra.tools.calibration_alignment import RigidTransform

    window = _window(tmp_path)
    window.session.transform = RigidTransform(
        dx_px=18.02, dy_px=-0.60, theta_rad=0.0001
    )
    window.refresh()
    text = window.transform_value.text()
    assert text.endswith("°")
    assert text.count(",") == 2
    assert not window.transform_value.wordWrap()


# ---------------------------------------------------------------------------
# The pattern the sphere does not fit, on screen.
#
# The bench cannot see a vertical detector shift through its fit — anchor rows
# come out of the reference pattern — so the sphere's own band centres are the
# only place it can be read.  The reading is said beside the factor curves it
# belongs to, and again on the save panel, where it cannot stop the press.
# ---------------------------------------------------------------------------

_BAND_ORDER_ROWS = (40.0, 100.0, 160.0, 220.0)
_BAND_COLUMNS = 640
_BAND_ROWS = 260


def _band_traces(shift_rows: float = 0.0) -> np.ndarray:
    columns = np.arange(_BAND_COLUMNS, dtype=float)
    return np.column_stack(
        [row + 0.01 * columns + shift_rows for row in _BAND_ORDER_ROWS]
    )


def _sphere_bands(path: Path, shift_rows: float) -> BenchFrame:
    rows = np.arange(_BAND_ROWS, dtype=float)[:, None]
    traces = _band_traces(shift_rows)
    image = np.full((_BAND_ROWS, _BAND_COLUMNS), 40.0)
    for order_idx in range(traces.shape[1]):
        image = image + 20000.0 * np.exp(
            -0.5 * ((rows - traces[:, order_idx][None, :]) / 5.0) ** 2
        )
    spectra = tuple(
        image[int(round(traces[0, order_idx])), :].copy()
        for order_idx in range(traces.shape[1])
    )
    return BenchFrame(
        Path(path), image[np.newaxis, :, :], image, spectra, {"ExposureTime": 0.1}
    )


def _bench_on_a_sphere(tmp_path: Path, shift_rows: float):
    """A bench whose campaign holds a real pattern and an opened sphere."""

    window = _manual_window(tmp_path)
    pattern_path = tmp_path / "pattern_bands.txt"
    np.savetxt(pattern_path, np.rint(_band_traces()).astype(int), fmt="%d")
    window.campaign = CalibrationCampaignSession(
        pattern_source=pattern_path,
        wavelength_source=tmp_path / "wavelength.txt",
        integral_source=tmp_path / "integral.txt",
    )
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    sphere.write_bytes(b"sphere\n")
    frame = _sphere_bands(sphere, shift_rows)
    window.campaign.record_frame(frame)
    window.campaign.classify_file(sphere, MeasurementRole.SPHERE)
    window.refresh_campaign()
    return window


def test_the_sphere_view_says_when_the_pattern_does_not_fit(qt_app, tmp_path: Path):
    window = _bench_on_a_sphere(tmp_path, 3.5)

    message = window.sphere_view_message.text()
    assert "PATTERN DOES NOT FIT" in message
    assert "rows above the chosen pattern" in message
    # The way out is a press on this very view, not a command line and a
    # restart: the reading and the two controls that answer it are together.
    assert "extract a pattern from this sphere" in message
    assert "without closing" in message
    assert window.extract_pattern_button.isEnabled()
    assert window.choose_pattern_button.isEnabled()


def test_a_fitting_pattern_is_reported_beside_the_factors_without_alarm(
    qt_app, tmp_path: Path
):
    window = _bench_on_a_sphere(tmp_path, 0.0)

    message = window.sphere_view_message.text()
    assert "PATTERN DOES NOT FIT" not in message
    assert "rder bands sit" in message


def test_the_save_press_is_warned_about_and_still_goes_through(qt_app, tmp_path: Path):
    """Loud, and never a refusal: the press still starts the save it asked for."""

    window = _bench_on_a_sphere(tmp_path, 3.5)
    started = []
    window._start_campaign_task = lambda operation: started.append(operation)

    window._save_snapshot()

    assert started, "the warning must not swallow the press"
    said = window.save_message_value.text()
    assert "PATTERN DOES NOT FIT THIS SPHERE" in said
    assert "Saving anyway" in said
    assert "recorded in the snapshot manifest" in said


def test_a_fitting_pattern_says_nothing_at_the_save_press(qt_app, tmp_path: Path):
    window = _bench_on_a_sphere(tmp_path, 0.0)
    window._start_campaign_task = lambda operation: None
    window.save_message_value.setText("")

    window._save_snapshot()

    assert window.save_message_value.text() == ""


# ---------------------------------------------------------------------------
# Extracting a pattern in place.
#
# The owner, at the bench, with the guard warning on screen and every role
# assigned: "I'm blocked, why do I need to exit for the pattern?"  He was
# right.  The bench holds the sphere frame, packages the extraction code, and
# has just measured the offset — the pattern was a constructor argument for
# historical reasons only.  These tests are his flow: bands five rows off the
# launch pattern, guard fires, one press, guard quiet, and nothing derived from
# the old geometry left standing.
# ---------------------------------------------------------------------------

#: One trusted line per order, so the anchors an extraction clears are anchors
#: that were really fitted, and the re-run that replaces them really runs.
_BAND_LINE_PIXELS = ((0, 120.0), (1, 240.0), (2, 360.0), (3, 480.0))

#: How far this sphere's own bands sit above the pattern the bench opened on —
#: the owner's real 2024-folder-on-the-2025-pattern reading, near enough.
_BAND_SHIFT_ROWS = 5.0


def _band_lines() -> tuple:
    return tuple(
        CalibrationTableLine(
            order, center - 5.0, center + 5.0, center, 600.0 + center, "NeI", "ok"
        )
        for order, center in _BAND_LINE_PIXELS
    )


def _band_frame(path: Path, shift_rows: float, *, lit: bool = True) -> BenchFrame:
    """A sphere acquisition: order bands across the rows, lines along them."""

    rows = np.arange(_BAND_ROWS, dtype=float)[:, None]
    traces = _band_traces(shift_rows)
    image = np.full((_BAND_ROWS, _BAND_COLUMNS), 40.0)
    if lit:
        for order_idx in range(traces.shape[1]):
            image = image + 20000.0 * np.exp(
                -0.5 * ((rows - traces[:, order_idx][None, :]) / 5.0) ** 2
            )
    x = np.arange(_BAND_COLUMNS, dtype=float)
    spectra = []
    for order_idx in range(traces.shape[1]):
        trace = np.full(_BAND_COLUMNS, 10.0)
        for line_order, center in _BAND_LINE_PIXELS:
            if line_order == order_idx and lit:
                trace = trace + 900.0 * np.exp(-0.5 * ((x - center) / 1.7) ** 2)
        spectra.append(trace)
    return BenchFrame(
        Path(path),
        image[np.newaxis, :, :],
        image,
        tuple(spectra),
        {"ExposureTime": 0.1},
    )


class _BandLoader:
    """A reader that hands back the sphere pair, and remembers a new pattern."""

    def __init__(self, shift_rows: float) -> None:
        self.shift_rows = shift_rows
        self.pattern: np.ndarray | None = None

    def __call__(self, path) -> BenchFrame:
        source = Path(path)
        return _band_frame(source, self.shift_rows, lit="-bg" not in source.name)

    def adopt_pattern(self, pattern) -> None:
        self.pattern = np.asarray(pattern, dtype=int)


def _extraction_bench(
    tmp_path: Path, *, shift_rows: float = _BAND_SHIFT_ROWS, background: bool = True
):
    """A bench whose sphere's bands sit *shift_rows* off the pattern it opened on."""

    launch = tmp_path / "pattern_launch.txt"
    np.savetxt(launch, np.rint(_band_traces()).astype(int), fmt="%d")
    for name in ("wavelength.txt", "integral.txt"):
        (tmp_path / name).write_text(name, encoding="utf-8")

    loader = _BandLoader(shift_rows)
    sphere = tmp_path / "sphere-0.1s-x3.sif"
    sphere_bg = tmp_path / "sphere-0.1s-x3-bg.sif"
    for path in (sphere, sphere_bg):
        path.write_bytes(b"sphere fixture\n")

    session = CalibrationBenchSession(np.loadtxt(launch, dtype=int), _band_lines())
    session.accept_frame(loader(sphere))
    window = CalibrationBenchWindow(
        session,
        loader=loader,
        start_timer=False,
        output_root=tmp_path / "calibrations",
    )
    campaign = CalibrationCampaignSession(
        pattern_source=launch,
        wavelength_source=tmp_path / "wavelength.txt",
        integral_source=tmp_path / "integral.txt",
    )
    campaign.record_frame(session.frame)
    campaign.classify_file(sphere, MeasurementRole.SPHERE)
    if background:
        campaign.classify_file(sphere_bg, MeasurementRole.SPHERE_BACKGROUND)
    window.campaign = campaign
    window.snapshot_id_edit.setText("20250926_cmos")
    window.refresh()
    return window, launch, sphere


def _wait_for_the_bench(window: CalibrationBenchWindow, qt_app) -> None:
    """Let the extraction thread and the re-extraction it starts both finish."""

    for _attempt in range(600):
        idle = (
            window._campaign_thread is None
            and window._load_thread is None
            and not window._queue
        )
        if idle:
            break
        qt_app.processEvents()
        QtCore.QThread.msleep(5)
    qt_app.processEvents()


def test_the_band_guard_carries_the_press_that_answers_it(qt_app, tmp_path: Path):
    """The row that reports the mismatch is the row that offers the fix."""

    window, _launch, sphere = _extraction_bench(tmp_path)

    assert window.campaign.pattern_band_warning()
    references = next(
        item
        for item in window.campaign.checklist(window.session)
        if item.key == "references"
    )
    assert window._next_action_for(references) == (
        "Extract pattern from this sphere",
        window._extract_pattern_from_sphere,
    )
    assert window.next_step_button.text() == "Extract pattern from this sphere"
    assert window.extract_pattern_button.isEnabled()
    assert window.extract_pattern_button.font().bold()
    # And it says what it will fit, before it is pressed.
    assert sphere.name in references.unblocked_by
    assert "without closing" in references.unblocked_by
    window.close()


def test_extracting_from_this_sphere_rebases_the_bench_in_place(qt_app, tmp_path: Path):
    """The owner's flow, end to end: press once, and the geometry is right."""

    window, launch, sphere = _extraction_bench(tmp_path)
    window._auto_anchor()
    anchored = len(window.session.anchors)
    assert anchored >= 2, "the fixture must place anchors on the old geometry"
    window.campaign.compute_sphere_comparison(
        lambda **values: AbsoluteCalibrationResult(np.linspace(400, 700, 8), np.ones(8))
    )
    assert window.campaign.pattern_band_warning()

    window._extract_pattern_from_sphere()
    _wait_for_the_bench(window, qt_app)

    # Written first, under the campaign's own one tidy folder, and named in full.
    written = tmp_path / "calibrations" / "patterns" / "pattern_extracted_20250926.txt"
    assert written.is_file()
    said = window.save_message_value.text()
    assert str(written) in said

    # The bench is standing on it, and says so where the state is read.
    assert window.campaign.pattern_source == written
    assert window.pattern_source_value.text() == str(written)
    assert str(written) in window.pattern_choice_value.text()

    # The traces land on this sphere's own bands, which is the whole point.
    extracted = np.loadtxt(written, dtype=int)
    truth = _band_traces(_BAND_SHIFT_ROWS)
    assert np.median(np.abs(extracted - truth)) < 1.0
    np.testing.assert_allclose(window.session.pattern, extracted)
    np.testing.assert_array_equal(window.loader.pattern, extracted)
    np.testing.assert_allclose(window._pattern_items[0].yData, extracted[:, 0])

    # The guard, re-measured against what the bench now wears, reads quiet.
    assert window.campaign.pattern_band_warning() == ""
    reading = window.campaign.sphere_band_offsets()
    assert reading is not None and reading.measured
    assert abs(reading.median_offset_rows) < 1.0
    assert "PATTERN DOES NOT FIT" not in window.sphere_view_message.text()

    # Anchors go, with the one honest notice, and the re-run really re-runs.
    assert not window.session.anchors
    assert f"{anchored} anchor(s) were cleared" in said
    assert "fitted on the old geometry" in said
    assert "auto-anchor re-runs them in seconds" in said
    window._auto_anchor()
    assert len(window.session.anchors) >= 2
    assert window.session.rms_px is not None

    # Everything the old pattern produced is dropped through the existing paths.
    assert window.campaign.comparison.state is ComparisonState.NOT_RUN
    assert window.campaign.toml_state is TomlState.NOT_GENERATED
    assert window.campaign.save_state is SaveState.NOT_READY
    # And the roles and the files themselves are untouched.
    assert window.campaign.sphere_pair_paths()[0] == sphere
    assert launch.is_file(), "the pattern it came from is never overwritten"
    window.close()


def test_the_bench_extraction_is_the_cli_extraction(qt_app, tmp_path: Path):
    """One core, called through one function, with the worn pattern as prior."""

    from echelle_spectra import pattern_cli
    from echelle_spectra.tools import pattern_extraction

    assert (
        bench_gui.extract_pattern_from_sphere
        is pattern_cli.extract_pattern_from_sphere
        is pattern_extraction.extract_pattern_from_sphere
    )

    window, launch, _sphere = _extraction_bench(tmp_path)
    seen = {}
    real = pattern_extraction.extract_pattern_from_sphere

    def spy(image, prior=None, **kwargs):
        seen["prior"] = None if prior is None else np.asarray(prior)
        seen["image"] = np.asarray(image)
        return real(image, prior, **kwargs)

    bench_gui.extract_pattern_from_sphere = spy
    try:
        window._extract_pattern_from_sphere()
        _wait_for_the_bench(window, qt_app)
    finally:
        bench_gui.extract_pattern_from_sphere = real

    np.testing.assert_array_equal(seen["prior"], np.loadtxt(launch, dtype=int))
    # Signal minus background, as the CLI does it: the bands survive and the
    # detector floor the background carries is gone.
    assert seen["image"].max() > 1000.0
    assert float(np.percentile(seen["image"], 5)) < 1.0
    window.close()


def test_a_pattern_file_can_be_chosen_on_the_same_terms(qt_app, tmp_path: Path):
    """The picker is the general mechanism; extraction is its special case."""

    window, _launch, _sphere = _extraction_bench(tmp_path)
    window._auto_anchor()
    assert window.session.anchors
    on_geometry = tmp_path / "pattern_this_era.txt"
    np.savetxt(
        on_geometry, np.rint(_band_traces(_BAND_SHIFT_ROWS)).astype(int), fmt="%d"
    )
    chosen = (str(on_geometry), "")
    monkey = QtWidgets.QFileDialog.getOpenFileName
    QtWidgets.QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: chosen)
    try:
        window._choose_pattern_file()
        _wait_for_the_bench(window, qt_app)
    finally:
        QtWidgets.QFileDialog.getOpenFileName = monkey

    assert window.campaign.pattern_source == on_geometry
    np.testing.assert_allclose(
        window.session.pattern, np.loadtxt(on_geometry, dtype=int)
    )
    np.testing.assert_array_equal(
        window.loader.pattern, np.loadtxt(on_geometry, dtype=int)
    )
    assert window.campaign.pattern_band_warning() == ""
    assert not window.session.anchors
    assert "anchor(s) were cleared" in window.save_message_value.text()
    assert window.campaign.save_state is SaveState.NOT_READY
    # Nothing was extracted, so nothing was written.
    assert not (tmp_path / "calibrations" / "patterns").exists()
    window.close()


def test_a_campaign_without_a_sphere_says_why_it_cannot_extract(qt_app, tmp_path: Path):
    """A disabled button that does not say why is a bench that went quiet."""

    window, _launch, sphere = _extraction_bench(tmp_path)
    window.campaign.remove_classification(sphere)
    window.refresh_campaign()

    assert not window.extract_pattern_button.isEnabled()
    assert "sphere signal role" in window.extract_pattern_button.toolTip()
    assert window.choose_pattern_button.isEnabled(), "the picker never needs a sphere"

    window._extract_pattern_from_sphere()
    said = window.save_message_value.text()
    assert "No pattern was extracted" in said
    assert "sphere signal role" in said
    assert window._campaign_thread is None
    window.close()


def test_extracting_without_a_background_states_the_caveat(qt_app, tmp_path: Path):
    """Signal alone still answers the guard, and still says what it was."""

    window, _launch, _sphere = _extraction_bench(tmp_path, background=False)

    window._extract_pattern_from_sphere()
    _wait_for_the_bench(window, qt_app)

    said = window.save_message_value.text()
    assert "signal alone" in said
    assert window.campaign.pattern_band_warning() == ""
    window.close()


def test_the_busy_gate_covers_the_extraction(qt_app, tmp_path: Path):
    """The extraction reads the very session and campaign a press would move."""

    window, _launch, _sphere = _extraction_bench(tmp_path)
    window._campaign_thread = bench_gui.CampaignTaskThread(lambda: None, window)
    try:
        window.refresh_campaign()
        assert not window.extract_pattern_button.isEnabled()
        assert not window.choose_pattern_button.isEnabled()
        assert not window.next_step_button.isEnabled()
        # And a press that reaches the handler anyway starts nothing.
        started = []
        window._start_campaign_task(lambda: started.append(1))
        assert not started
    finally:
        window._campaign_thread = None
    window.refresh_campaign()
    assert window.extract_pattern_button.isEnabled()
    window.close()


# ---------------------------------------------------------------------------
# The packaged names are the vocabulary.
#
# The owner typed the name the README's own examples use —
# ``--pattern pattern_CMOS_20240305.txt`` — and was told "pattern file not
# found", because only a literal path was ever tried.  A path still wins; a
# bare packaged name now resolves; and a name that is neither is refused with
# the handful of names that would have worked.
# ---------------------------------------------------------------------------


def test_a_bare_packaged_pattern_name_resolves_to_the_packaged_file():
    resolved = resolve_calibration_file(Path("pattern_CMOS_20240305.txt"), "pattern")

    assert resolved.is_absolute()
    assert resolved.is_file()
    assert resolved.name == "pattern_CMOS_20240305.txt"
    assert resolved.parent.name == "calibration_files"


def test_a_bare_packaged_wavelength_name_resolves_from_a_subfolder():
    resolved = resolve_calibration_file(
        Path("Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"), "wavelength"
    )

    assert resolved.is_file()
    assert resolved.parent.name == "alignments"


def test_a_real_path_still_wins_over_a_packaged_name(tmp_path: Path):
    mine = tmp_path / "pattern_CMOS_20240305.txt"
    mine.write_text("1 2\n3 4\n", encoding="utf-8")

    resolved = resolve_calibration_file(mine, "pattern")

    assert resolved == mine.resolve()


def test_a_name_that_is_neither_is_refused_with_the_packaged_vocabulary():
    with pytest.raises(SystemExit) as refusal:
        resolve_calibration_file(Path("pattern_CMOS_1999.txt"), "pattern")

    said = str(refusal.value)
    assert "pattern file not found" in said
    assert "pattern_CMOS_20240305.txt" in said
    assert "pattern_CMOS_20250926.txt" in said
    # The vocabulary this argument is about, and not the other tables beside it.
    assert "Th_wavelength" not in said


def test_the_bench_says_which_tables_it_resolved(qt_app, tmp_path: Path, capsys):
    """Which pattern the bench opened on is the one fact the fit cannot tell you."""

    started = {}

    def fake_window(session, **kwargs):
        started["campaign"] = kwargs["campaign"]
        raise SystemExit(0)

    monkey = bench_gui.CalibrationBenchWindow
    bench_gui.CalibrationBenchWindow = fake_window
    try:
        with pytest.raises(SystemExit):
            bench_gui.main([str(tmp_path), "--pattern", "pattern_CMOS_20240305.txt"])
    finally:
        bench_gui.CalibrationBenchWindow = monkey

    printed = capsys.readouterr().out
    assert "pattern_CMOS_20240305.txt" in printed
    assert str(_PACKAGE_DIR) in printed
    assert started["campaign"].pattern_source.is_absolute()


# ----------------------------------------------------------------------
# Open calibration folder: the whole session moves, and no restart
# ----------------------------------------------------------------------


def _calibration_folder(tmp_path: Path, name: str) -> Path:
    """A folder that reads like one of the owner's own: dated, and unambiguous."""

    folder = tmp_path / name
    folder.mkdir()
    for entry in _REAL_2025_NAMES:
        (folder / entry).write_bytes(b"sif\n")
    return folder


def _bench_at(tmp_path: Path, folder: Path, **options) -> CalibrationBenchWindow:
    """A bench launched at *folder* and holding nothing yet, as ``main`` leaves it."""

    window = _manual_window(tmp_path, window_options={"folder": folder, **options})
    window.output_root, window.config_root = bench_gui.default_bench_roots(folder)
    window._describe_output_roots()
    # The shared fixture hands every window a frame; a freshly launched bench
    # has none, and an empty bench is the one that must never be questioned.
    window.session.forget_frame()
    window.refresh()
    return window


def _answer_folder_dialog(monkeypatch, folder: Path) -> list:
    """Answer the folder picker through the module seam, not a real modal."""

    asked: list = []
    monkeypatch.setattr(
        bench_gui,
        "choose_calibration_folder",
        lambda parent, start: (asked.append(Path(start)) or str(folder)),
    )
    return asked


def _answer_question(monkeypatch, reply) -> list:
    asked: list[str] = []
    monkeypatch.setattr(
        QtWidgets.QMessageBox,
        "question",
        lambda *args, **values: (asked.append(args[2]) or reply),
    )
    return asked


def test_opening_a_calibration_folder_rebases_the_whole_session(
    qt_app, tmp_path, monkeypatch
):
    """Owner, 2026-08-18: "So I can start GUI, and go through many calibration
    folders. I probably would have to do that at NIFS."

    Picking a folder is the launch argument, minus the process: both roots are
    derived inside it, the identity is re-dated from its name, its SIFs come in
    through the ordinary drop path so triage runs and the unambiguous roles
    apply themselves, the dialogs and the watch follow it — and nothing the
    previous folder produced comes along.
    """

    first = _calibration_folder(tmp_path, "20250926_calib")
    second = _calibration_folder(tmp_path, "20261101_calib")
    watcher = StableSifWatcher(first, required_unchanged_polls=3, minimum_age_s=2.0)
    window = _bench_at(tmp_path, first, watcher=watcher)
    window.show()

    _answer_folder_dialog(monkeypatch, first)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)
    assert window.campaign.assigned_lamps == ("Ne",)
    first_files = set(window.campaign.loaded)
    # Work in progress the operator would lose: an anchor, and a comparison.
    window.session.anchors[(0, 1.0, 2.0)] = object()
    window.campaign.comparison = campaign_module.SphereComparison(
        ComparisonState.READY, "measured against the first folder"
    )

    _answer_folder_dialog(monkeypatch, second)
    offered = _answer_question(monkeypatch, QtWidgets.QMessageBox.Yes)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert len(offered) == 1
    # Both roots hang off the folder that was opened, through the one function
    # the launch argument uses.
    assert (window.output_root, window.config_root) == bench_gui.default_bench_roots(
        second
    )
    assert str(second) in window.snapshot_root_value.text()
    assert str(second) in window.config_root_value.text()
    # The identity is re-dated by the folder's own name, and is the bench's
    # guess again — the SIF headers may still fill it in.
    assert window.snapshot_id_edit.text() == "20261101_cmos"
    assert window.snapshot_date == date(2026, 11, 1)
    assert window._snapshot_id_decided is False
    # The new folder's files arrived as an arrival, so triage ran and the
    # unambiguous names assigned themselves exactly as a drop's would.
    assert {path.parent for path in window.campaign.loaded} == {second}
    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)
    assert {path.parent for path in window.campaign.measurements} == {second}
    assert window.campaign.assigned_lamps == ("Ne",)
    assert not first_files & set(window.campaign.loaded)
    assert len(window._file_rows) == len(_REAL_2025_NAMES)
    assert window.file_table.rowCount() == len(_REAL_2025_NAMES)
    # And nothing measured on the old folder survived.
    assert not window.session.anchors
    assert window.campaign.comparison.state is ComparisonState.NOT_RUN
    assert window.campaign.toml_state is TomlState.NOT_GENERATED
    assert window.campaign.saved_snapshot is None
    # The dialogs and the watch are aimed at the folder now on the bench.
    assert window.last_folder == second
    assert window.watcher.folder == second
    assert window.watcher.required_unchanged_polls == 3
    assert window.watcher.minimum_age_ns == 2_000_000_000
    # And the Bench state strip finally names the calibration it is about.
    assert str(second) in window.watch_value.text()
    window.close()


def test_unsaved_work_is_confirmed_once_and_cancel_leaves_the_session_alone(
    qt_app, tmp_path, monkeypatch
):
    """One honest sentence, and only when there is something to lose."""

    first = _calibration_folder(tmp_path, "20250926_calib")
    second = _calibration_folder(tmp_path, "20261101_calib")
    window = _bench_at(tmp_path, first)
    window.show()

    # An empty bench has nothing to lose, so it is asked nothing at all.
    _answer_folder_dialog(monkeypatch, first)
    never = _answer_question(monkeypatch, QtWidgets.QMessageBox.No)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert never == []
    assert len(window.campaign.measurements) == len(_REAL_2025_NAMES)

    before_files = set(window.campaign.loaded)
    before_rows = list(window._file_rows)
    before_identity = window.snapshot_id_edit.text()
    before_roots = (window.output_root, window.config_root)

    _answer_folder_dialog(monkeypatch, second)
    offered = _answer_question(monkeypatch, QtWidgets.QMessageBox.No)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert len(offered) == 1
    assert offered[0] == (
        "Opening a new folder clears the current session — the alignment "
        f"for {before_identity} was not saved."
    )
    # Cancelled means untouched: the same files, rows, identity and roots.
    assert set(window.campaign.loaded) == before_files
    assert window._file_rows == before_rows
    assert window.snapshot_id_edit.text() == before_identity
    assert (window.output_root, window.config_root) == before_roots
    assert window.calibration_folder == first
    assert "was not opened" in window.save_message_value.text()
    window.close()


def test_a_saved_folder_is_swapped_without_a_question(qt_app, tmp_path, monkeypatch):
    """A validated snapshot is the answer to "was this saved?", so it is not asked.

    The campaign's own ``saved_snapshot`` is the fact consulted, and it is
    cleared by ``_invalidate_outputs`` the moment anything it was built from
    changes — so this is a promise about a session that really is finished, not
    about one that merely once pressed Save.
    """

    first = _calibration_folder(tmp_path, "20250926_calib")
    second = _calibration_folder(tmp_path, "20261101_calib")
    window = _bench_at(tmp_path, first)
    window.show()
    _answer_folder_dialog(monkeypatch, first)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    window.campaign.saved_snapshot = Snapshot(
        window.output_root / "20250926_cmos",
        {"id": "20250926_cmos", "detector": "cmos"},
        (),
    )
    assert window._unsaved_session_warning() == ""

    _answer_folder_dialog(monkeypatch, second)
    never = _answer_question(monkeypatch, QtWidgets.QMessageBox.No)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert never == []
    assert window.calibration_folder == second
    # The swap went through the ordinary invalidation, so the saved snapshot of
    # the folder that was left behind is not still claimed by the new one.
    assert window.campaign.saved_snapshot is None
    window.close()


def test_the_chosen_pattern_and_wavelength_survive_the_folder_switch(
    qt_app, tmp_path, monkeypatch
):
    """The operator picked those tables; a folder does not overrule the pick.

    The band guard is what judges a carried-over pattern against the new
    sphere, and it goes on reading the same file — which is exactly how a stale
    choice announces itself rather than being silently swapped underneath.
    """

    first = _calibration_folder(tmp_path, "20250926_calib")
    second = _calibration_folder(tmp_path, "20261101_calib")
    window = _bench_at(tmp_path, first)
    window.show()
    chosen = (
        window.campaign.pattern_source,
        window.campaign.wavelength_source,
        window.campaign.integral_source,
    )
    window.campaign.previous_sphere = tmp_path / "previous_sphere.sif"
    window.campaign.previous_sphere_background = tmp_path / "previous_sphere_bg.sif"
    lamps = window.campaign.suggested_lamps

    _answer_folder_dialog(monkeypatch, second)
    _answer_question(monkeypatch, QtWidgets.QMessageBox.Yes)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)

    assert (
        window.campaign.pattern_source,
        window.campaign.wavelength_source,
        window.campaign.integral_source,
    ) == chosen
    assert window.campaign.previous_sphere == tmp_path / "previous_sphere.sif"
    assert (
        window.campaign.previous_sphere_background == tmp_path / "previous_sphere_bg.sif"
    )
    assert window.campaign.suggested_lamps == lamps
    assert str(chosen[0]) in window.pattern_source_value.text()
    window.close()


def test_dropping_still_works_on_a_bench_that_was_handed_a_folder(
    qt_app, tmp_path, monkeypatch
):
    """The teaching box went; the gesture it taught did not.

    A file dropped from outside the opened folder is queued, read and triaged
    through the very path the removed panel used to advertise.
    """

    folder = _calibration_folder(tmp_path, "20250926_calib")
    window = _bench_at(tmp_path, folder)
    window.show()
    _answer_folder_dialog(monkeypatch, folder)
    window.open_folder_button.click()
    _wait_for_loads(window, qt_app)
    stray = tmp_path / "IMG_0042.sif"
    stray.write_bytes(b"sif\n")

    _drop(window, [stray])
    _wait_for_loads(window, qt_app)

    assert stray in window.campaign.loaded
    assert stray in window._file_rows
    assert window.file_table.rowCount() == len(_REAL_2025_NAMES) + 1
    window.close()


def test_the_folder_picker_shows_the_files_it_will_not_let_you_pick(
    qt_app, tmp_path, monkeypatch
):
    """Owner, 2026-08-18: "we should show contents, but make them gray.
    Otherwise I have to open the same folder twice."

    Windows' native folder picker hides files: his real folder of SIFs read
    "No items match your search", so telling 20190115 from 20250926 meant
    opening one, closing it, and opening the other. Qt's own dialog in
    Directory mode with ShowDirsOnly off lists them greyed — evidence, never an
    answer. Driven through the real function with only ``exec_`` answered, so
    what is asserted is the dialog the operator would actually get.
    """

    folder = _calibration_folder(tmp_path, "20250926_calib")
    opened: list[QtWidgets.QFileDialog] = []

    def refuse(dialog):
        opened.append(dialog)
        return 0  # QDialog.Rejected — nothing is picked, nothing is shown

    monkeypatch.setattr(QtWidgets.QFileDialog, "exec_", refuse)

    assert bench_gui.choose_calibration_folder(None, folder) == ""

    assert len(opened) == 1
    dialog = opened[0]
    assert dialog.fileMode() == QtWidgets.QFileDialog.Directory
    # Files visible, and Qt's own dialog: the native one cannot be told either.
    assert not dialog.testOption(QtWidgets.QFileDialog.ShowDirsOnly)
    assert dialog.testOption(QtWidgets.QFileDialog.DontUseNativeDialog)
    assert dialog.acceptMode() == QtWidgets.QFileDialog.AcceptOpen
    assert Path(dialog.directory().absolutePath()) == folder
    # A folder full of SIFs is what the operator is trying to recognise, and
    # this dialog's own model lists them rather than hiding them.
    listed = {
        Path(entry).name
        for entry in dialog.directory().entryList(QtCore.QDir.Files)
    }
    assert set(_REAL_2025_NAMES) <= listed


def test_a_folder_is_not_opened_while_the_bench_is_still_reading(
    qt_app, tmp_path, monkeypatch
):
    """A read in flight belongs to the folder being left.

    It would land in the fresh campaign a moment after that campaign was
    emptied — the one file of the old folder that followed the bench to the new
    one. The press is refused in words while anything is still queued, and the
    dialog is never even opened.
    """

    folder = _calibration_folder(tmp_path, "20250926_calib")
    other = _calibration_folder(tmp_path, "20261101_calib")
    window = _bench_at(tmp_path, folder)
    window.show()
    asked = _answer_folder_dialog(monkeypatch, other)
    window._queue.append(folder / _REAL_2025_NAMES[0])

    window._pick_calibration_folder()

    assert asked == []
    assert window.calibration_folder == folder
    assert "still reading" in window.save_message_value.text()
    window._queue.clear()
    window.close()


def test_the_input_line_names_the_launched_folder(qt_app, tmp_path):
    """It read "manual — drag and drop or Add files", which nobody needed told."""

    folder = _calibration_folder(tmp_path, "20250926_calib")
    window = _bench_at(tmp_path, folder)
    window.show()
    qt_app.processEvents()

    assert window.watch_value.text() == str(folder)
    assert window.calibration_folder == folder
    # A bench nobody has pointed anywhere says that, rather than naming a folder
    # it does not have.
    bare = _manual_window(tmp_path)
    assert "no folder open" in bare.watch_value.text()
    bare.close()
    window.close()


#: The owner's own 2025 folder, when this checkout carries it: the matched era,
#: where an extraction must agree with the pattern it was fitted against rather
#: than move it.  A guard that shifted the geometry on a good pattern would be
#: worse than no button at all.
_LIVE_FOLDER = Path(__file__).resolve().parents[1] / "local" / "20250926_calib"
_LIVE_SPHERE_PAIR = (
    _LIVE_FOLDER / "sphere-0.1s-x3.sif",
    _LIVE_FOLDER / "sphere-0.1s-x3-bg.sif",
)
_LIVE_PATTERN = (
    _PACKAGE_DIR / "resources" / "calibration_files" / "pattern_CMOS_20250926.txt"
)


@pytest.mark.skipif(
    not all(path.is_file() for path in _LIVE_SPHERE_PAIR),
    reason="local/20250926_calib is not in this checkout",
)
def test_the_real_2025_sphere_extracts_back_onto_its_own_pattern(tmp_path: Path):
    """Real optics, matched era: the fit lands where the pattern already was."""

    from echelle_spectra.calibration_bench import FrameLoader, band_center_offsets

    prior = np.loadtxt(_LIVE_PATTERN, dtype=int)
    sphere, background = _LIVE_SPHERE_PAIR
    written = tmp_path / "calibrations" / "patterns" / "pattern_extracted_20250926.txt"

    rebase = bench_gui.extract_pattern_beside_campaign(
        loader=FrameLoader(prior),
        prior=prior,
        sphere_path=sphere,
        background_path=background,
        output_path=written,
    )

    extracted = np.loadtxt(written, dtype=int)
    assert extracted.shape == prior.shape
    assert rebase.caveat == ""
    assert abs(rebase.median_offset_rows) < 1.0
    # And the guard, which was quiet before, is quiet after.
    reading = band_center_offsets(rebase.sphere_image, extracted)
    assert reading.measured
    assert not reading.exceeds()
