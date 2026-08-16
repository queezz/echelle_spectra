"""Known-line markers on the raw 2-D detector image.

The spectrum overlays in :mod:`~echelle_spectra.tools.line_overlay` answer
"where is this line in wavelength".  On the detector image the same question is
geometric: a loaded calibration already carries both halves of the map — the
per-order wavelength solution sampled at every column, and the order pattern
giving that order's centre row at every column — so inverting the solution
turns a catalog wavelength into a ``(column, row)`` on the sensor.

Two conventions this module keeps rather than re-derives:

* ``detector_pixel`` is the **pre-flip raw column**.  ``Spectrum`` flips its
  stitched 1-D arrays when ``direction < 0`` (wavelength falling with column on
  the black Echelle) so the spectrum plots read left-to-right in wavelength.
  The image is never flipped, and neither is ``clbr.order_wavel``, which is
  indexed by raw column.  Inverting it therefore lands on the raw column
  directly — applying the flip here would mirror every mark.
* Which order owns a wavelength is decided by ``clbr.order_borders``, the same
  mask that stitches the 1-D spectrum, so a line in the overlap of two orders
  is marked on the order the spectrum plot actually shows it from.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg

from .line_catalog import LINE_FAMILIES, SpectralLine
from .line_overlay import LINE_OVERLAY_STYLES, select_overlay_lines

__all__ = [
    "DetectorGeometry",
    "DetectorLineMark",
    "DetectorLineOverlay",
]

#: A detector pixel spans ``[index, index + 1)`` in image coordinates, so a mark
#: meant for the middle of a pixel is drawn half a pixel up and to the right of
#: its index.  Marks themselves carry plain index space, which is what the
#: wavelength solution, the pattern, and ``detector_pixel`` all speak.
_PIXEL_CENTER = 0.5


@dataclass(frozen=True)
class DetectorLineMark:
    """Where one catalog line is expected to land on the sensor."""

    family: str
    label: str
    wavelength_nm: float
    order: int
    order_index: int
    column: float
    row: float
    half_height: float


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

        ids = np.asarray(order_ids if order_ids is not None else range(self.order_count))
        if ids.size != self.order_count:
            ids = np.arange(self.order_count)
        self.order_ids = tuple(int(value) for value in ids)

        # One inverted solution per order: wavelengths ascending beside the raw
        # columns they came from, so a monotonically falling order inverts with
        # the same interpolation as a rising one.
        self._grid = np.arange(self.columns, dtype=float)
        self._solutions: list[tuple[np.ndarray, np.ndarray] | None] = []
        for index in range(self.order_count):
            wavelengths = self.order_wavel[index]
            usable = np.isfinite(wavelengths)
            if order_borders is not None:
                usable = usable & np.asarray(order_borders[index], dtype=bool)
            columns = self._grid[usable]
            values = wavelengths[usable]
            if values.size < 2:
                self._solutions.append(None)
                continue
            ascending = np.argsort(values)
            self._solutions.append((values[ascending], columns[ascending]))

    @classmethod
    def from_calibration(cls, calibration) -> "DetectorGeometry | None":
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
            )
        except (ValueError, TypeError, IndexError):
            return None

    def wavelength_span(self) -> tuple[float, float] | None:
        """The nanometre range the whole detector can show, or ``None``."""

        spans = [solution[0] for solution in self._solutions if solution is not None]
        if not spans:
            return None
        return (
            float(min(values[0] for values in spans)),
            float(max(values[-1] for values in spans)),
        )

    def column_for(self, order_index: int, wavelength_nm: float) -> float | None:
        """Invert one order's wavelength solution; ``None`` when out of range."""

        solution = self._solutions[order_index]
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

    def mark_for(self, line: SpectralLine) -> DetectorLineMark | None:
        """Place one catalog line, or ``None`` when no order carries it."""

        for index in range(self.order_count):
            column = self.column_for(index, line.wavelength_nm)
            if column is None:
                continue
            return DetectorLineMark(
                family=line.family,
                label=line.label,
                wavelength_nm=float(line.wavelength_nm),
                order=self.order_ids[index],
                order_index=index,
                column=column,
                row=self.row_at(index, column),
                half_height=self.half_height,
            )
        return None

    def marks_for(self, lines) -> tuple[DetectorLineMark, ...]:
        """Place every catalog line that lands on the sensor, once each."""

        placed = (self.mark_for(line) for line in lines)
        return tuple(mark for mark in placed if mark is not None)


class DetectorLineOverlay:
    """Pooled per-family tick marks over the 2-D detector image.

    One :class:`pyqtgraph.PlotDataItem` per family carries every mark of that
    family as disconnected segments, so a family costs one item and one draw
    call however many lines it places.  Nothing is created until a family is
    first switched on, and nothing recomputes on pan or zoom: the marks depend
    on the calibration and the toggles, and on nothing else.
    """

    def __init__(self, plot, *, max_marks: int = 200, tick_px: float = 4.0):
        self._plot = plot
        self._geometry: DetectorGeometry | None = None
        self._visible = {family: False for family in LINE_FAMILIES}
        self._items: dict[str, pg.PlotDataItem] = {}
        self._marks: dict[str, tuple[DetectorLineMark, ...]] = {
            family: () for family in LINE_FAMILIES
        }
        self.max_marks = int(max_marks)
        self.tick_px = float(tick_px)

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
            self._geometry = DetectorGeometry.from_calibration(source)
        self.refresh()
        return self._geometry is not None

    def set_family_visible(self, family: str, visible: bool) -> int:
        """Show or hide one family; returns how many marks are now drawn."""

        key = family.strip().lower()
        if key not in self._visible:
            known = ", ".join(LINE_FAMILIES)
            raise ValueError(f"unknown line family {family!r}; known families: {known}")
        self._visible[key] = bool(visible)
        self.refresh(key)
        return len(self._marks[key])

    def is_family_visible(self, family: str) -> bool:
        return self._visible[family]

    def marks(self, family: str) -> tuple[DetectorLineMark, ...]:
        """The marks currently drawn for one family, primarily for QA."""

        return self._marks[family]

    def item(self, family: str) -> pg.PlotDataItem | None:
        """The pooled item for one family, or ``None`` while it costs nothing."""

        return self._items.get(family)

    def refresh(self, family: str | None = None) -> None:
        """Redraw after a toggle or a calibration change — never on a view change."""

        for key in (family,) if family is not None else LINE_FAMILIES:
            self._refresh_family(key)

    def _refresh_family(self, family: str) -> None:
        marks = self._compute(family) if self._visible[family] else ()
        self._marks[family] = marks
        item = self._items.get(family)
        if not marks:
            if item is not None:
                item.setData(x=[], y=[])
                item.setVisible(False)
            return
        if item is None:
            item = self._create_item(family)
        columns, rows = self._segments(marks)
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

    def _create_item(self, family: str) -> pg.PlotDataItem:
        style = LINE_OVERLAY_STYLES[family]
        item = pg.PlotDataItem(pen=pg.mkPen(style.color, width=1.2), connect="pairs")
        item.setZValue(20)
        self._plot.addItem(item, ignoreBounds=True)
        self._items[family] = item
        return item

    def _segments(self, marks) -> tuple[np.ndarray, np.ndarray]:
        """Bracket each line with two short ticks, leaving the order band clear."""

        columns = np.repeat([mark.column + _PIXEL_CENTER for mark in marks], 4)
        rows = np.empty(len(marks) * 4, dtype=float)
        for index, mark in enumerate(marks):
            center = mark.row + _PIXEL_CENTER
            low = center - mark.half_height
            high = center + mark.half_height
            rows[4 * index : 4 * index + 4] = (
                low - self.tick_px,
                low,
                high,
                high + self.tick_px,
            )
        return columns, rows
