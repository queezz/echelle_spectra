"""Range-aware pyqtgraph overlays backed by the shared line catalog."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import pyqtgraph as pg

from .line_catalog import LINE_FAMILIES, SpectralLine, load_line_table

__all__ = [
    "LINE_OVERLAY_STYLES",
    "LineOverlayManager",
    "OverlayStyle",
    "select_overlay_lines",
]


@dataclass(frozen=True)
class OverlayStyle:
    """Visual identity for one line family."""

    color: str
    dash: tuple[int, int]


LINE_OVERLAY_STYLES: Mapping[str, OverlayStyle] = {
    "balmer": OverlayStyle("#f778ba", (6, 3)),
    "fulcher": OverlayStyle("#36c9d0", (3, 2)),
    "thar": OverlayStyle("#ff9f43", (5, 2)),
    "ne": OverlayStyle("#70d66b", (4, 2)),
    "hg": OverlayStyle("#ffd166", (2, 2)),
}


def _atomic_strength_threshold(width_nm: float) -> float:
    if width_nm <= 1.0:
        return 0.0
    if width_nm <= 5.0:
        return 0.05
    if width_nm <= 20.0:
        return 0.12
    return 0.25


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
