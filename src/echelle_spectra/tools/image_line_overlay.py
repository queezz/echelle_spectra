"""Known-line markers on the raw 2-D detector image.

The spectrum overlays in :mod:`~echelle_spectra.tools.line_overlay` answer
"where is this line in wavelength".  On the detector image the same question is
geometric: a loaded calibration already carries both halves of the map — the
per-order wavelength solution sampled at every column, and the order pattern
giving that order's centre row at every column — so inverting the solution
turns a catalog wavelength into a ``(column, row)`` on the sensor.

Three conventions this module keeps rather than re-derives:

* ``detector_pixel`` is the **pre-flip raw column**.  ``Spectrum`` flips its
  stitched 1-D arrays when ``direction < 0`` (wavelength falling with column on
  the black Echelle) so the spectrum plots read left-to-right in wavelength.
  The image is never flipped, and neither is ``clbr.order_wavel``, which is
  indexed by raw column.  Inverting it therefore lands on the raw column
  directly — applying the flip here would mirror every mark.
* In an echelle overlap the same wavelength is **physically exposed twice**, on
  two adjacent orders, and the image shows both blobs.  So the image marks both:
  ``clbr.order_borders`` — the mask that stitches the 1-D spectrum — decides
  which of the two is the *primary* mark, and the twin is drawn beside it in a
  secondary pen rather than hidden.  Hiding it, as this module first shipped,
  left a real blob unmarked and made the mark on the stitch owner read as
  misplaced.
* The mark is a **box, not a bracket**.  The first shipping drew two ticks just
  outside the extraction band, which on this instrument's fat PSF is inside the
  blob: measured on ``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif`` a Ne
  line is ~6 px FWHM across columns but ~26 px FWHM down rows, while the
  extraction band is only ``2*dv+1 = 17`` px.  Ticks at ±8..±12 rows were being
  painted over the brightest pixels of the very blob they marked.  The box now
  spans the order's band out to the neighbouring traces, so its edges land in
  the dark gutter between orders where they can be seen.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore

from .line_catalog import LINE_FAMILIES, SpectralLine
from .line_overlay import LINE_OVERLAY_STYLES, select_overlay_lines

__all__ = [
    "CURSOR_LINK_COLOR",
    "CURSOR_LINK_RATE_HZ",
    "DEFAULT_LINE_WIDTH_PX",
    "CursorLink",
    "DetectorGeometry",
    "DetectorLineMark",
    "DetectorLineOverlay",
    "OrderTraceOverlay",
]

#: A detector pixel spans ``[index, index + 1)`` in image coordinates, so a mark
#: meant for the middle of a pixel is drawn half a pixel up and to the right of
#: its index.  Marks themselves carry plain index space, which is what the
#: wavelength solution, the pattern, and ``detector_pixel`` all speak.
_PIXEL_CENTER = 0.5

#: Rough line width in columns, the box's full width.  Measured on the owner's
#: bright Ne frame ``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif``: over 26
#: unsaturated Ne blobs the column FWHM was 6.0 px (median, 6.0 px at the 90th
#: percentile), so ~2x FWHM is 12 px.  Wide enough that a box wraps its blob with
#: visible margin, narrow enough that neighbouring lines 20 px apart keep their
#: own boxes.
DEFAULT_LINE_WIDTH_PX = 12.0

#: Order traces on the detector image, matching the calibration bench's own
#: detector view (``calibration_bench_gui._refresh_pattern_traces``) so the two
#: windows draw the same thing the same way.
ORDER_TRACE_COLOR = "#7b91a4"
ORDER_TRACE_WIDTH = 0.9


def _brighter(color: str, amount: float = 0.45) -> str:
    """Blend a family colour toward white so it reads over saturated pixels."""

    red, green, blue, _alpha = pg.mkColor(color).getRgb()
    lift = lambda value: int(round(value + (255 - value) * amount))
    return "#{:02x}{:02x}{:02x}".format(lift(red), lift(green), lift(blue))


@dataclass(frozen=True)
class DetectorLineMark:
    """Where one catalog line is expected to land on the sensor.

    ``primary`` is ``True`` for the order the stitched 1-D spectrum takes this
    wavelength from, and ``False`` for the twin the overlap exposes on the
    neighbouring order.
    """

    family: str
    label: str
    wavelength_nm: float
    order: int
    order_index: int
    column: float
    row: float
    half_height: float
    half_width: float = DEFAULT_LINE_WIDTH_PX / 2.0
    primary: bool = True


class DetectorGeometry:
    """The λ → (column, row-band) map a loaded calibration already carries."""

    def __init__(
        self,
        pattern,
        order_wavel,
        *,
        order_ids=None,
        order_borders=None,
        half_height: float = 8.0,
        line_width_px: float = DEFAULT_LINE_WIDTH_PX,
    ):
        pattern = np.asarray(pattern, dtype=float)
        order_wavel = np.asarray(order_wavel, dtype=float)
        if pattern.ndim != 2 or order_wavel.ndim != 2:
            raise ValueError("pattern and wavelength solution must both be 2-D")
        if (pattern.shape[1], pattern.shape[0]) != order_wavel.shape:
            raise ValueError(
                "pattern {}x{} and wavelength solution {}x{} disagree about the "
                "detector".format(*pattern.shape, *order_wavel.shape)
            )
        self.pattern = pattern
        self.order_wavel = order_wavel
        self.columns = int(pattern.shape[0])
        self.order_count = int(pattern.shape[1])
        self.half_height = float(half_height)
        self.half_width = float(line_width_px) / 2.0

        ids = np.asarray(order_ids if order_ids is not None else range(self.order_count))
        if ids.size != self.order_count:
            ids = np.arange(self.order_count)
        self.order_ids = tuple(int(value) for value in ids)

        # Two inverted solutions per order, wavelengths ascending beside the raw
        # columns they came from so a falling order inverts like a rising one.
        # ``_carried`` is every column the order actually exposes — that is what
        # the image shows, overlap twins included.  ``_owned`` is the same
        # solution trimmed by the stitch mask, and answers only "which of the
        # two orders does the 1-D spectrum read this wavelength from".
        self._grid = np.arange(self.columns, dtype=float)
        self._carried = self._invert(None)
        self._owned = self._carried if order_borders is None else self._invert(order_borders)

    def _invert(self, order_borders):
        solutions: list[tuple[np.ndarray, np.ndarray] | None] = []
        for index in range(self.order_count):
            wavelengths = self.order_wavel[index]
            usable = np.isfinite(wavelengths)
            if order_borders is not None:
                usable = usable & np.asarray(order_borders[index], dtype=bool)
            columns = self._grid[usable]
            values = wavelengths[usable]
            if values.size < 2:
                solutions.append(None)
                continue
            ascending = np.argsort(values)
            solutions.append((values[ascending], columns[ascending]))
        return solutions

    @classmethod
    def from_calibration(
        cls, calibration, *, line_width_px: float = DEFAULT_LINE_WIDTH_PX
    ) -> "DetectorGeometry | None":
        """Read a started :class:`~echelle_spectra.tools.echelle.Calibrations`.

        Returns ``None`` for anything that cannot answer the question — a
        calibration that never finished loading, or one whose tables disagree —
        so a caller can say so instead of failing the image display.
        """

        if calibration is None:
            return None
        pattern = getattr(calibration, "pattern", None)
        order_wavel = getattr(calibration, "order_wavel", None)
        if pattern is None or order_wavel is None:
            return None
        try:
            return cls(
                pattern,
                order_wavel,
                order_ids=getattr(calibration, "order_ids", None),
                order_borders=getattr(calibration, "order_borders", None),
                half_height=float(getattr(calibration, "dv", 8)),
                line_width_px=line_width_px,
            )
        except (ValueError, TypeError, IndexError):
            return None

    def wavelength_span(self) -> tuple[float, float] | None:
        """The nanometre range the whole detector can show, or ``None``."""

        spans = [solution[0] for solution in self._carried if solution is not None]
        if not spans:
            return None
        return (
            float(min(values[0] for values in spans)),
            float(max(values[-1] for values in spans)),
        )

    def column_for(self, order_index: int, wavelength_nm: float) -> float | None:
        """Invert one order's wavelength solution; ``None`` when out of range."""

        return self._column_from(self._carried[order_index], wavelength_nm)

    @staticmethod
    def _column_from(solution, wavelength_nm: float) -> float | None:
        if solution is None:
            return None
        values, columns = solution
        target = float(wavelength_nm)
        if target < values[0] or target > values[-1]:
            return None
        return float(np.interp(target, values, columns))

    def row_at(self, order_index: int, column: float) -> float:
        """The order's centre row at a (possibly fractional) column."""

        return float(np.interp(column, self._grid, self.pattern[:, order_index]))

    def band_half_height(self, order_index: int, column: float) -> float:
        """Half the order's band at one column: the reach to its neighbour.

        The extraction band is only ``±dv`` rows, but the light is not — a line
        on this instrument runs ~26 px FWHM down the rows against a 17-px
        extraction band.  Taking half the distance to the nearer neighbouring
        trace gives a box that wraps the whole blob and stops in the dark gutter
        between orders, and never reaches into the neighbour's own band.  It is
        floored at ``±dv`` so the box always contains the rows the spectrum was
        actually extracted from, and falls back to it when there is no
        neighbour to measure against.
        """

        center = self.row_at(order_index, column)
        gaps = [
            abs(self.row_at(index, column) - center)
            for index in (order_index - 1, order_index + 1)
            if 0 <= index < self.order_count
        ]
        if not gaps:
            return self.half_height
        return max(self.half_height, min(gaps) / 2.0)

    def _mark(self, line: SpectralLine, index: int, column: float, primary: bool):
        return DetectorLineMark(
            family=line.family,
            label=line.label,
            wavelength_nm=float(line.wavelength_nm),
            order=self.order_ids[index],
            order_index=index,
            column=column,
            row=self.row_at(index, column),
            half_height=self.band_half_height(index, column),
            half_width=self.half_width,
            primary=primary,
        )

    def owner_index(self, wavelength_nm: float) -> int | None:
        """Which order the stitched 1-D spectrum reads this wavelength from."""

        for index in range(self.order_count):
            if self._column_from(self._owned[index], wavelength_nm) is not None:
                return index
        return None

    def marks_for_line(self, line: SpectralLine) -> tuple[DetectorLineMark, ...]:
        """Place one catalog line on **every** order that exposes it.

        An echelle overlap puts the same line on two adjacent orders, and the
        detector shows both blobs.  The one the stitch mask owns comes back
        ``primary``; its twin comes back as a duplicate.  When the trimmed mask
        owns none of them — a wavelength that falls in the seam itself — the
        first carrying order is named primary so a line is never all-duplicate.
        """

        owner = self.owner_index(line.wavelength_nm)
        placed = []
        for index in range(self.order_count):
            column = self.column_for(index, line.wavelength_nm)
            if column is None:
                continue
            placed.append((index, column))
        if not placed:
            return ()
        if owner is None or owner not in [index for index, _ in placed]:
            owner = placed[0][0]
        return tuple(
            self._mark(line, index, column, index == owner) for index, column in placed
        )

    def mark_for(self, line: SpectralLine) -> DetectorLineMark | None:
        """The primary mark for one line, or ``None`` when no order carries it."""

        for mark in self.marks_for_line(line):
            if mark.primary:
                return mark
        return None

    def marks_for(self, lines) -> tuple[DetectorLineMark, ...]:
        """Place every catalog line that lands on the sensor, twins included."""

        return tuple(mark for line in lines for mark in self.marks_for_line(line))

    # -- the same map read forwards, for a cursor rather than for a catalog --

    def wavelength_at(self, order_index: int, column: float) -> float | None:
        """One order's sampled solution read at a (possibly fractional) column.

        ``None`` off the end of the part of the order the solution covers — a
        partial order runs out of sensor and its solution is NaN-padded there,
        and a padded column carries no wavelength to report.
        """

        values = self.order_wavel[order_index]
        usable = np.isfinite(values)
        if int(usable.sum()) < 2:
            return None
        columns = self._grid[usable]
        if column < columns[0] or column > columns[-1]:
            return None
        return float(np.interp(column, columns, values[usable]))

    def order_at(self, column: float, row: float) -> int | None:
        """Which order's band a detector row falls in at one column.

        The band is the same reach the line boxes use — out to half the
        distance to the neighbouring trace, floored at ``dv`` — so a cursor
        anywhere on a blob names the order that blob belongs to.  A cursor in
        the dark gutter beyond the outermost trace belongs to no order, and the
        honest answer there is ``None`` rather than the nearest one.
        """

        if not 0.0 <= column <= self.columns - 1:
            return None
        nearest: tuple[float, int] | None = None
        for index in range(self.order_count):
            distance = abs(float(row) - self.row_at(index, column))
            if distance > self.band_half_height(index, column):
                continue
            if nearest is None or distance < nearest[0]:
                nearest = (distance, index)
        return None if nearest is None else nearest[1]

    def image_position(self, wavelength_nm: float):
        """Where the stitched spectrum reads one wavelength from the sensor.

        Returns ``(order_index, column, row)`` for the order the stitch mask
        owns — the primary mark's order — or ``None`` when no order carries it.
        """

        index = self.owner_index(wavelength_nm)
        if index is None:
            return None
        column = self.column_for(index, wavelength_nm)
        if column is None:
            return None
        return index, column, self.row_at(index, column)


class DetectorLineOverlay:
    """Pooled per-family line boxes over the 2-D detector image.

    At most two :class:`pyqtgraph.PlotDataItem` per family — one for the marks
    the stitched spectrum owns, one for the twins an order overlap exposes on
    the neighbouring order — each carrying every box of that kind as
    disconnected segments.  So a family costs two items and two draw calls
    however many lines it places.  Nothing is created until a family is first
    switched on, and nothing recomputes on pan or zoom: the marks depend on the
    calibration and the toggles, and on nothing else.

    ``line_width_px`` is handed to any geometry this overlay builds itself from
    a calibration.  A :class:`DetectorGeometry` passed in ready-made already
    carries its own, and keeps it — the box's extent belongs to the geometry
    that measured it, not to the widget drawing it.
    """

    def __init__(
        self,
        plot,
        *,
        max_marks: int = 200,
        line_width_px: float = DEFAULT_LINE_WIDTH_PX,
    ):
        self._plot = plot
        self._geometry: DetectorGeometry | None = None
        self._visible = {family: False for family in LINE_FAMILIES}
        self._items: dict[tuple[str, bool], pg.PlotDataItem] = {}
        self._marks: dict[str, tuple[DetectorLineMark, ...]] = {
            family: () for family in LINE_FAMILIES
        }
        self.max_marks = int(max_marks)
        self.line_width_px = float(line_width_px)

    @property
    def geometry(self) -> DetectorGeometry | None:
        return self._geometry

    def set_geometry(self, source) -> bool:
        """Adopt the geometry of a freshly loaded image, or drop a stale one.

        Accepts a :class:`DetectorGeometry`, a started ``Calibrations``, or
        ``None``.  Returns whether the overlay can now place marks.
        """

        if source is None or isinstance(source, DetectorGeometry):
            self._geometry = source
        else:
            self._geometry = DetectorGeometry.from_calibration(
                source, line_width_px=self.line_width_px
            )
        self.refresh()
        return self._geometry is not None

    def set_family_visible(self, family: str, visible: bool) -> int:
        """Show or hide one family; returns how many **lines** are now drawn.

        A line doubled across an order overlap is one line with two boxes, so
        the count the status bar reports is the number of primary marks.  Ask
        :meth:`duplicate_count` for the twins.
        """

        key = family.strip().lower()
        if key not in self._visible:
            known = ", ".join(LINE_FAMILIES)
            raise ValueError(f"unknown line family {family!r}; known families: {known}")
        self._visible[key] = bool(visible)
        self.refresh(key)
        return self.line_count(key)

    def is_family_visible(self, family: str) -> bool:
        return self._visible[family]

    def marks(self, family: str) -> tuple[DetectorLineMark, ...]:
        """Every box currently drawn for one family, primarily for QA."""

        return self._marks[family]

    def primary_marks(self, family: str) -> tuple[DetectorLineMark, ...]:
        return tuple(mark for mark in self._marks[family] if mark.primary)

    def duplicate_marks(self, family: str) -> tuple[DetectorLineMark, ...]:
        return tuple(mark for mark in self._marks[family] if not mark.primary)

    def line_count(self, family: str) -> int:
        """How many catalog lines this family currently marks."""

        return len(self.primary_marks(family))

    def duplicate_count(self, family: str) -> int:
        """How many of those lines are doubled onto a neighbouring order."""

        return len(self.duplicate_marks(family))

    def item(self, family: str) -> pg.PlotDataItem | None:
        """The pooled primary item, or ``None`` while it costs nothing."""

        return self._items.get((family, True))

    def duplicate_item(self, family: str) -> pg.PlotDataItem | None:
        """The pooled overlap-twin item, or ``None`` while it costs nothing."""

        return self._items.get((family, False))

    def refresh(self, family: str | None = None) -> None:
        """Redraw after a toggle or a calibration change — never on a view change."""

        for key in (family,) if family is not None else LINE_FAMILIES:
            self._refresh_family(key)

    def _refresh_family(self, family: str) -> None:
        marks = self._compute(family) if self._visible[family] else ()
        self._marks[family] = marks
        for primary in (True, False):
            subset = [mark for mark in marks if mark.primary is primary]
            item = self._items.get((family, primary))
            if not subset:
                if item is not None:
                    item.setData(x=[], y=[])
                    item.setVisible(False)
                continue
            if item is None:
                item = self._create_item(family, primary)
            columns, rows = self._segments(subset)
            item.setData(x=columns, y=rows, connect="pairs")
            item.setVisible(True)

    def _compute(self, family: str) -> tuple[DetectorLineMark, ...]:
        if self._geometry is None:
            return ()
        span = self._geometry.wavelength_span()
        if span is None:
            return ()
        lines = select_overlay_lines(family, span[0], span[1], max_labels=self.max_marks)
        return self._geometry.marks_for(lines)

    def _create_item(self, family: str, primary: bool) -> pg.PlotDataItem:
        style = LINE_OVERLAY_STYLES[family]
        if primary:
            pen = pg.mkPen(_brighter(style.color), width=2.0)
        else:
            # Same hue, unmistakably secondary: the twin an overlap exposes is
            # real light, but the stitched spectrum is not reading it from here.
            pen = pg.mkPen(style.color, width=2.0, style=QtCore.Qt.DashLine)
            pen.setDashPattern([float(value) for value in style.dash])
        item = pg.PlotDataItem(pen=pen, connect="pairs", antialias=False)
        item.setZValue(20 if primary else 19)
        self._plot.addItem(item, ignoreBounds=True)
        self._items[(family, primary)] = item
        return item

    def _segments(self, marks) -> tuple[np.ndarray, np.ndarray]:
        """One open rectangle per line: four ``connect="pairs"`` edges.

        Width is a rough line width; height is the order's band at that column.
        The box is not filled — the point is to ring the blob, never to hide it.
        """

        columns = np.empty(len(marks) * 8, dtype=float)
        rows = np.empty(len(marks) * 8, dtype=float)
        for index, mark in enumerate(marks):
            left = mark.column + _PIXEL_CENTER - mark.half_width
            right = mark.column + _PIXEL_CENTER + mark.half_width
            low = mark.row + _PIXEL_CENTER - mark.half_height
            high = mark.row + _PIXEL_CENTER + mark.half_height
            columns[8 * index : 8 * index + 8] = (
                left, right,  # bottom edge
                right, right,  # right edge
                right, left,  # top edge
                left, left,  # left edge
            )
            rows[8 * index : 8 * index + 8] = (
                low, low,
                low, high,
                high, high,
                high, low,
            )
        return columns, rows


#: The cursor link belongs to no line family, so it takes none of their hues:
#: a neutral near-white that reads over both the dark detector and the plots.
CURSOR_LINK_COLOR = "#e6ecf2"

#: How often the pointer is allowed to move the markers.  Mouse moves arrive
#: far faster than a plot can usefully be redrawn, so they are proxied down to
#: this; ~30 Hz is smooth to the hand and cheap to the scene.
CURSOR_LINK_RATE_HZ = 30.0

#: Half the arms of the crosshair drawn on the image, in detector pixels.
CURSOR_CROSS_HALF_PX = 16.0


class CursorLink:
    """One pointer, marked in both views — while it is switched on.

    The detector image and the spectrum plots show the same light in two
    coordinate systems, and the loaded calibration already maps between them.
    This reads that map at the pointer: a cursor over the image names the order
    its row falls in and the wavelength its column carries, and marks that
    wavelength on the spectra; a cursor over a spectrum marks where the
    stitched trace read that wavelength from the sensor.

    Off is genuinely free.  Nothing is connected and no item exists until
    :meth:`set_enabled` is called with ``True``: mouse movement over the image
    is the highest-rate event in the window, and a viewer that is not using the
    link must not pay a slot call for every one of them.  Switched on, the
    moves are rate-limited through :class:`pyqtgraph.SignalProxy`, and each
    view owns exactly one marker item that is moved and hidden — never rebuilt
    per move.
    """

    def __init__(
        self,
        image_plot,
        spectrum_plots,
        *,
        readout=None,
        rate_hz: float = CURSOR_LINK_RATE_HZ,
        color: str = CURSOR_LINK_COLOR,
    ):
        self._image_plot = image_plot
        self._spectrum_plots = dict(spectrum_plots)
        self._readout = readout
        self._rate_hz = float(rate_hz)
        self._color = color
        self._geometry: DetectorGeometry | None = None
        self._enabled = False
        self._proxies: list[pg.SignalProxy] = []
        self._image_marker: pg.PlotDataItem | None = None
        self._spectrum_markers: dict[str, pg.InfiniteLine] = {}
        self._label = ""

    # ------------------------------------------------------------- state ---

    @property
    def geometry(self) -> DetectorGeometry | None:
        return self._geometry

    @property
    def is_enabled(self) -> bool:
        return self._enabled

    @property
    def label(self) -> str:
        """What the link last had to say, or ``""`` when it is saying nothing."""

        return self._label

    def proxy_count(self) -> int:
        """How many mouse-move proxies are connected; zero while switched off."""

        return len(self._proxies)

    def set_geometry(self, source) -> bool:
        """Adopt a freshly loaded calibration, or drop a stale one."""

        if source is None or isinstance(source, DetectorGeometry):
            self._geometry = source
        else:
            self._geometry = DetectorGeometry.from_calibration(source)
        if self._geometry is None:
            self.clear()
        return self._geometry is not None

    def set_enabled(self, enabled: bool) -> bool:
        """Switch the link on or off; returns whether it can map anything yet."""

        wanted = bool(enabled)
        if wanted == self._enabled:
            return self._geometry is not None
        self._enabled = wanted
        if wanted:
            self._connect()
        else:
            self._disconnect()
            self.clear()
        return self._geometry is not None

    # ------------------------------------------------------------ wiring ---

    def _scenes(self):
        """Every distinct scene the linked plots live in, in a stable order.

        The two spectrum plots share one ``GraphicsLayoutWidget`` and therefore
        one scene, so connecting per plot would deliver every move twice.
        """

        scenes = []
        for plot in (self._image_plot, *self._spectrum_plots.values()):
            scene = plot.scene()
            if scene is not None and not any(scene is known for known in scenes):
                scenes.append(scene)
        return scenes

    def _connect(self) -> None:
        for scene in self._scenes():
            self._proxies.append(
                pg.SignalProxy(
                    scene.sigMouseMoved,
                    rateLimit=self._rate_hz,
                    slot=self._pointer_moved,
                )
            )

    def _disconnect(self) -> None:
        for proxy in self._proxies:
            proxy.disconnect()
        self._proxies = []

    def _pointer_moved(self, event) -> None:
        position = event[0]
        if self._image_plot.sceneBoundingRect().contains(position):
            point = self._image_plot.getViewBox().mapSceneToView(position)
            self.show_for_image_point(
                point.x() - _PIXEL_CENTER, point.y() - _PIXEL_CENTER
            )
            return
        for name, plot in self._spectrum_plots.items():
            if plot.sceneBoundingRect().contains(position):
                point = plot.getViewBox().mapSceneToView(position)
                self.show_for_wavelength(point.x(), source=name)
                return
        self.clear()

    # ------------------------------------------------------------ marking ---

    def show_for_image_point(self, column: float, row: float):
        """Mark the wavelength a detector position carries; ``None`` if none.

        Returns ``(order, wavelength_nm)`` using the order's own *identifier*,
        which is what the operator reads on the calibration bench, not its
        index in the pattern.
        """

        if not self._enabled or self._geometry is None:
            return None
        index = self._geometry.order_at(column, row)
        if index is None:
            self.clear()
            return None
        wavelength = self._geometry.wavelength_at(index, column)
        if wavelength is None:
            self.clear()
            return None
        order = self._geometry.order_ids[index]
        self._mark_spectra(wavelength)
        self._hide_image_marker()
        self._announce(order, wavelength)
        return order, float(wavelength)

    def show_for_wavelength(self, wavelength_nm: float, *, source: str | None = None):
        """Mark where the stitched spectrum read one wavelength from the sensor.

        Returns ``(order, column, row)``, or ``None`` when no order owns it.
        """

        if not self._enabled or self._geometry is None:
            return None
        placed = self._geometry.image_position(wavelength_nm)
        if placed is None:
            self.clear()
            return None
        index, column, row = placed
        order = self._geometry.order_ids[index]
        self._mark_image(column, row)
        self._mark_spectra(wavelength_nm, skip=source)
        self._announce(order, wavelength_nm)
        return order, float(column), float(row)

    def clear(self) -> None:
        """Take both marks down — the pointer has left, or has nothing to say."""

        self._hide_image_marker()
        for marker in self._spectrum_markers.values():
            marker.setVisible(False)
        if self._label:
            self._label = ""
            if self._readout is not None:
                self._readout("")

    def image_marker(self) -> pg.PlotDataItem | None:
        """The pooled crosshair, or ``None`` while the link costs nothing."""

        return self._image_marker

    def spectrum_marker(self, name: str) -> pg.InfiniteLine | None:
        """One plot's pooled marker, or ``None`` while the link costs nothing."""

        return self._spectrum_markers.get(name)

    def _announce(self, order: int, wavelength_nm: float) -> None:
        self._label = f"order {order} · {wavelength_nm:.2f} nm"
        if self._readout is not None:
            self._readout(self._label)

    def _pen(self):
        return pg.mkPen(self._color, width=1.0)

    def _mark_spectra(self, wavelength_nm: float, *, skip: str | None = None) -> None:
        for name, plot in self._spectrum_plots.items():
            if name == skip:
                # The pointer itself is the mark on the plot it is over, and a
                # stale line left under it would be the only thing lying.
                stale = self._spectrum_markers.get(name)
                if stale is not None:
                    stale.setVisible(False)
                continue
            marker = self._spectrum_markers.get(name)
            if marker is None:
                marker = pg.InfiniteLine(angle=90, movable=False, pen=self._pen())
                marker.setZValue(30)
                plot.addItem(marker, ignoreBounds=True)
                self._spectrum_markers[name] = marker
            marker.setPos(float(wavelength_nm))
            marker.setVisible(True)

    def _mark_image(self, column: float, row: float) -> None:
        if self._image_marker is None:
            self._image_marker = pg.PlotDataItem(
                pen=self._pen(), connect="pairs", antialias=False
            )
            self._image_marker.setZValue(30)
            self._image_plot.addItem(self._image_marker, ignoreBounds=True)
        x = column + _PIXEL_CENTER
        y = row + _PIXEL_CENTER
        self._image_marker.setData(
            x=[x - CURSOR_CROSS_HALF_PX, x + CURSOR_CROSS_HALF_PX, x, x],
            y=[y, y, y - CURSOR_CROSS_HALF_PX, y + CURSOR_CROSS_HALF_PX],
            connect="pairs",
        )
        self._image_marker.setVisible(True)

    def _hide_image_marker(self) -> None:
        if self._image_marker is not None:
            self._image_marker.setVisible(False)


class OrderTraceOverlay:
    """The order pattern drawn over the detector image, as one pooled curve.

    The calibration bench shows the same traces on its own detector view; this
    is that view's answer for the main window, where the operator is checking
    that the frame in front of him sits on the calibration he loaded.  Every
    order lives in a single :class:`pyqtgraph.PlotDataItem` joined by ``NaN``
    gaps and ``connect="finite"``, so the whole pattern costs one item and one
    draw call.  Off by default, nothing built until it is switched on, and
    nothing recomputed on pan or zoom.
    """

    def __init__(
        self,
        plot,
        *,
        color: str = ORDER_TRACE_COLOR,
        width: float = ORDER_TRACE_WIDTH,
    ):
        self._plot = plot
        self._geometry: DetectorGeometry | None = None
        self._item: pg.PlotDataItem | None = None
        self._visible = False
        self.color = color
        self.width = float(width)

    @property
    def geometry(self) -> DetectorGeometry | None:
        return self._geometry

    @property
    def is_visible(self) -> bool:
        return self._visible

    def item(self) -> pg.PlotDataItem | None:
        """The pooled curve, or ``None`` while the overlay costs nothing."""

        return self._item

    def order_count(self) -> int:
        """How many traces are currently drawn."""

        if not self._visible or self._geometry is None:
            return 0
        return self._geometry.order_count

    def set_geometry(self, source) -> bool:
        """Adopt a freshly loaded calibration's pattern, or drop a stale one."""

        if source is None or isinstance(source, DetectorGeometry):
            self._geometry = source
        else:
            self._geometry = DetectorGeometry.from_calibration(source)
        self.refresh()
        return self._geometry is not None

    def set_visible(self, visible: bool) -> int:
        """Show or hide the traces; returns how many orders are now drawn."""

        self._visible = bool(visible)
        self.refresh()
        return self.order_count()

    def refresh(self) -> None:
        columns, rows = self._curve()
        if columns is None:
            if self._item is not None:
                self._item.setData(x=[], y=[])
                self._item.setVisible(False)
            return
        if self._item is None:
            self._item = pg.PlotDataItem(
                pen=pg.mkPen(self.color, width=self.width),
                connect="finite",
                antialias=False,
            )
            self._item.setZValue(10)
            self._plot.addItem(self._item, ignoreBounds=True)
        self._item.setData(x=columns, y=rows, connect="finite")
        self._item.setVisible(True)

    def _curve(self):
        """Every trace end to end, separated by a ``NaN`` the pen does not cross."""

        if not self._visible or self._geometry is None:
            return None, None
        pattern = self._geometry.pattern
        span = np.arange(pattern.shape[0], dtype=float) + _PIXEL_CENTER
        gap = np.array([np.nan])
        columns = np.concatenate(
            [value for _ in range(pattern.shape[1]) for value in (span, gap)]
        )
        rows = np.concatenate(
            [
                value
                for index in range(pattern.shape[1])
                for value in (pattern[:, index] + _PIXEL_CENTER, gap)
            ]
        )
        return columns, rows
