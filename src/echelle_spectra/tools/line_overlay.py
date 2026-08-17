"""Range-aware pyqtgraph overlays backed by the shared line catalog."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import pyqtgraph as pg

from .line_catalog import LINE_FAMILIES, SpectralLine, load_line_table

__all__ = [
    "LINE_OVERLAY_STYLES",
    "SPECTRUM_CURVE_COLORS",
    "LineOverlayManager",
    "OverlayStyle",
    "select_overlay_lines",
]


#: The two spectrum curves an overlay is drawn on top of, by plot name.
#:
#: The raw counts trace is yellow-green and the calibrated trace is pure red,
#: and a marker is useless when it is the colour of the curve it annotates —
#: the Ne family shipped green sticks over the green counts curve and vanished
#: into it.  The palette below is answerable to these two values, so they live
#: here beside it rather than only in the GUI that sets the pens.
SPECTRUM_CURVE_COLORS: Mapping[str, str] = {
    "counts": "#6ac600",
    "calibrated": "#ff0000",
}


@dataclass(frozen=True)
class OverlayStyle:
    """Visual identity for one line family."""

    color: str
    dash: tuple[int, int]


#: One hue per family, used by the 1-D sticks and labels and, lifted toward
#: white, by the 2-D detector boxes — so a family is the same colour in both
#: views and the operator learns five colours rather than ten.
#:
#: Every hue sits in the cool arc.  The two curve hues are red at 0 deg and
#: yellow-green at 88 deg, and they *bracket* the warm colours: orange, amber,
#: and yellow all fall between the two curves and read as a dim version of one
#: of them on the dark plot background.  So the palette starts past green at
#: 157 deg and runs to magenta at 314 deg, which leaves every family at least
#: 45 deg from either curve.
#:
#: Within that arc the three lamp families are placed as far apart as it
#: allows — 157 / 228 / 314 deg — because they are the ones an operator turns
#: on together while identifying a lamp.  The hydrogen families fill the gaps,
#: and no two families are closer than 33 deg.  Dash patterns stay distinct as
#: a second channel, but colour does the work.
LINE_OVERLAY_STYLES: Mapping[str, OverlayStyle] = {
    "balmer": OverlayStyle("#b377f8", (6, 3)),  # violet, 268 deg
    "fulcher": OverlayStyle("#2bcdee", (3, 2)),  # cyan, 190 deg
    "thar": OverlayStyle("#30e8a2", (5, 2)),  # mint, 157 deg
    "ne": OverlayStyle("#5677fb", (4, 2)),  # blue, 228 deg
    "hg": OverlayStyle("#f368d2", (2, 2)),  # magenta, 314 deg
}


def _atomic_strength_threshold(width_nm: float) -> float:
    """How strong a cached lamp row must be to earn room in a view this wide.

    These are floors on
    :attr:`~.line_catalog.SpectralLine.relative_intensity`, which is a
    lamp-context strength spanning five decades across the packaged caches, so
    they are not evenly spaced.  They were set against the owner's bright Ne
    frame ``local/20250926_calib/Ne-0.02s-x3-bright-lines.sif``: over the whole
    401--802 nm detector the widest floor keeps 36 Ne rows of which 21 sit on a
    measurable blob and exactly one is an ion, where a floor of 0.02 would keep
    253 rows, 228 of them on dark detector and 170 of those Ne II.  0.08 rather
    than 0.10 because it costs Ne nothing — the same 36 rows — and buys Hg its
    576.96 nm line back.

    The narrowest view keeps everything.  Zooming that far in *is* the request
    to see the weak rows, ionized stages included.
    """

    if width_nm <= 1.0:
        return 0.0
    if width_nm <= 5.0:
        return 0.002
    if width_nm <= 20.0:
        return 0.02
    return 0.08


def _evenly_spaced(lines: list[SpectralLine], count: int) -> list[SpectralLine]:
    if count <= 1:
        return lines[:1]
    last = len(lines) - 1
    indices = {round(index * last / (count - 1)) for index in range(count)}
    return [lines[index] for index in sorted(indices)]


def select_overlay_lines(
    family: str,
    minimum_nm: float,
    maximum_nm: float,
    *,
    max_labels: int = 14,
) -> tuple[SpectralLine, ...]:
    """Choose legible labels for the current horizontal view.

    Atomic tables reveal progressively weaker cached NIST rows as the user
    zooms.  Dense equal-weight molecular tables are sampled across a broad view
    and become complete in a narrow view.
    """

    if max_labels < 1:
        return ()
    low, high = sorted((float(minimum_nm), float(maximum_nm)))
    width = high - low
    in_view = [
        line
        for line in load_line_table(family)
        if low <= line.wavelength_nm <= high
    ]
    if family in {"thar", "ne", "hg"}:
        threshold = _atomic_strength_threshold(width)
        in_view = [
            line
            for line in in_view
            if line.relative_intensity is not None
            and line.relative_intensity >= threshold
        ]
    if len(in_view) <= max_labels:
        return tuple(in_view)
    if any(line.relative_intensity is not None for line in in_view):
        strongest = sorted(
            in_view,
            key=lambda line: (
                -(line.relative_intensity or 0.0),
                line.wavelength_nm,
                line.label,
            ),
        )[:max_labels]
        return tuple(sorted(strongest, key=lambda line: line.wavelength_nm))
    return tuple(_evenly_spaced(in_view, max_labels))


class LineOverlayManager:
    """Own labeled line items for one or more wavelength plots."""

    def __init__(self, *, max_labels: int = 14):
        self.max_labels = max_labels
        self._plots: dict[str, object] = {}
        self._labels: dict[str, bool] = {}
        self._visible = {family: False for family in LINE_FAMILIES}
        self._items: dict[str, dict[str, list[pg.InfiniteLine]]] = {}
        self._rendered: dict[str, dict[str, tuple[SpectralLine, ...]]] = {}

    def register_plot(self, name: str, plot, *, labels: bool = True) -> None:
        """Attach a wavelength plot and update it whenever its x range changes."""

        if name in self._plots:
            raise ValueError(f"overlay plot {name!r} is already registered")
        self._plots[name] = plot
        self._labels[name] = labels
        self._items[name] = {family: [] for family in LINE_FAMILIES}
        self._rendered[name] = {family: () for family in LINE_FAMILIES}
        plot.getViewBox().sigXRangeChanged.connect(
            lambda *_args, plot_name=name: self.refresh(plot_name)
        )

    def set_family_visible(self, family: str, visible: bool) -> None:
        """Show or hide a family on every registered plot."""

        key = family.strip().lower()
        if key not in self._visible:
            known = ", ".join(LINE_FAMILIES)
            raise ValueError(f"unknown line family {family!r}; known families: {known}")
        self._visible[key] = bool(visible)
        for name in self._plots:
            self.refresh(name, family=key)

    def is_family_visible(self, family: str) -> bool:
        return self._visible[family]

    def rendered_lines(self, plot_name: str, family: str) -> tuple[SpectralLine, ...]:
        """Return the catalog records currently drawn, primarily for QA."""

        return self._rendered[plot_name][family]

    def refresh(self, plot_name: str | None = None, *, family: str | None = None) -> None:
        """Rebuild visible markers after data or zoom changes."""

        names = (plot_name,) if plot_name is not None else tuple(self._plots)
        families = (family,) if family is not None else LINE_FAMILIES
        for name in names:
            plot = self._plots[name]
            low, high = plot.getViewBox().viewRange()[0]
            for key in families:
                for item in self._items[name][key]:
                    plot.removeItem(item)
                self._items[name][key] = []
                self._rendered[name][key] = ()
                if not self._visible[key]:
                    continue
                table = load_line_table(key)
                family_low = max(low, table[0].wavelength_nm)
                family_high = min(high, table[-1].wavelength_nm)
                family_width = max(0.0, family_high - family_low)
                view_width = max(float(high - low), 1e-12)
                occupied_pixels = float(plot.getViewBox().width()) * family_width / view_width
                spaced_limit = max(1, int(occupied_pixels / 52.0))
                selected = select_overlay_lines(
                    key,
                    low,
                    high,
                    max_labels=min(self.max_labels, spaced_limit),
                )
                style = LINE_OVERLAY_STYLES[key]
                pen = pg.mkPen(style.color, width=1.2, dash=style.dash)
                for index, line in enumerate(selected):
                    label = line.label if self._labels[name] else None
                    label_options = None
                    if label is not None:
                        label_options = {
                            "position": 0.96 - 0.12 * (index % 3),
                            "color": style.color,
                            "fill": (0, 0, 0, 155),
                        }
                    item = pg.InfiniteLine(
                        pos=line.wavelength_nm,
                        angle=90,
                        movable=False,
                        pen=pen,
                        label=label,
                        labelOpts=label_options,
                    )
                    item.setZValue(20)
                    plot.addItem(item, ignoreBounds=True)
                    self._items[name][key].append(item)
                self._rendered[name][key] = selected
