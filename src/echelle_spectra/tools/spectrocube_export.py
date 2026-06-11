"""Export calibrated Echelle spectra to the SpectroCube format.

This module provides a thin bridge between ``echelle_spectra``'s
:class:`~echelle_spectra.tools.echelle.Spectrum` object and the neutral
SpectroCube container defined by the ``spectrocube`` package.

Design notes
------------
* The conversion lives here (in ``echelle_spectra``) because only this
  package knows the meaning of its internal data products.
* ``spectrocube`` is an **optional** dependency; an :exc:`ImportError` with a
  helpful message is raised at call-time if it is absent.
* The exporter never re-implements file I/O; it always delegates to
  :meth:`spectrocube.SpectroCube.save` / :meth:`spectrocube.SpectroCube.load`.
"""

from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from spectrocube import SpectroCube

__all__ = ["to_spectrocube", "export_spectrocube"]

# ---------------------------------------------------------------------------
# Unit / calibration-type metadata table
# ---------------------------------------------------------------------------

#: Maps the *units* keyword used by :class:`~echelle_spectra.tools.echelle.Spectrum`
#: to SpectroCube metadata fields and the matching attribute name on Spectrum.
_UNIT_INFO: dict[str, dict[str, str]] = {
    "counts": {
        "calibration_type": "counts",
        "intensity_units": "counts",
        "spectrum_attr": "counts",
    },
    "wm": {
        "calibration_type": "absolute",
        "intensity_units": "W/m2/nm",
        "spectrum_attr": "wm",
    },
    "wmsr": {
        "calibration_type": "absolute",
        "intensity_units": "W/m2/nm/sr",
        "spectrum_attr": "wmsr",
    },
    "phmsr": {
        "calibration_type": "absolute",
        "intensity_units": "ph/s/nm/sr",
        "spectrum_attr": "phmsr",
    },
}

_DEFAULT_INSTRUMENT_ID = "echelle"

# Fraction of axis points allowed to be non-monotonic before we refuse to
# co-sort and instead raise an error.  Echelle order stitching typically
# produces a handful of seam reversals, well below 1% of the axis.
_MAX_NONMONOTONIC_FRACTION = 0.02

# Absolute floor: always allow up to this many reversals regardless of axis
# size.  An echelle spectrum has ~30 orders so well under this in practice.
_MAX_NONMONOTONIC_ABSOLUTE = 50


def _json_attr(value: object) -> str:
    """Serialize structured metadata into a stable NetCDF string attribute."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _resolve_calibration_file(folder: object, filename: object) -> Path | None:
    if not filename:
        return None
    path = Path(str(filename))
    if not path.is_absolute() and folder is not None:
        path = Path(str(folder)) / path
    return path


def _file_digest(path: Path | None) -> dict[str, object] | None:
    if path is None or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(path),
        "size_bytes": int(path.stat().st_size),
        "sha256": digest.hexdigest(),
    }


# ---------------------------------------------------------------------------
# Wavelength-axis normalization
# ---------------------------------------------------------------------------


def _normalize_wavelength_axis(
    wavelength: np.ndarray,
    intensity_2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(wavelength, intensity_2d)`` with a strictly increasing wavelength axis.

    SpectroCube requires a strictly monotonically increasing wavelength
    coordinate.  ``Spectrum.wavelength`` from ``echelle_spectra`` may be:

    * Strictly increasing — pass through unchanged.
    * Strictly decreasing — reverse both arrays (cheap, exact mapping).
    * Globally increasing/decreasing with a small number of local seam
      reversals at order boundaries.  In this case the paired
      ``(wavelength, intensity)`` data is co-sorted by wavelength.  This is
      safe because each wavelength stays paired with its own intensity slice.
    * Genuinely scrambled — :exc:`ValueError`.

    Duplicate wavelength values are also rejected with a clear error.

    Parameters
    ----------
    wavelength : np.ndarray
        1-D wavelength axis (nm).
    intensity_2d : np.ndarray
        2-D ``(frame, wavelength)`` intensity array.  The wavelength axis is
        the **last** dimension.

    Returns
    -------
    (np.ndarray, np.ndarray)
        Normalized wavelength and intensity arrays.

    Raises
    ------
    ValueError
        If the axis cannot be cleanly converted to strictly increasing.
    """
    if wavelength.ndim != 1:
        raise ValueError(f"wavelength must be 1-D, got ndim={wavelength.ndim}")
    if intensity_2d.shape[-1] != wavelength.size:
        raise ValueError(
            f"intensity last axis ({intensity_2d.shape[-1]}) does not match "
            f"wavelength size ({wavelength.size})"
        )

    dw = np.diff(wavelength)
    n_neg = int(np.sum(dw < 0))
    n_pos = int(np.sum(dw > 0))
    n_zero = int(np.sum(dw == 0))
    n_total = dw.size

    # Strictly increasing — nothing to do.
    if n_neg == 0 and n_zero == 0:
        return wavelength, intensity_2d

    # Strictly decreasing — exact reversal preserves pixel mapping.
    if n_pos == 0 and n_zero == 0:
        return wavelength[::-1].copy(), intensity_2d[..., ::-1].copy()

    # Mostly increasing (or decreasing) with a few seam reversals — co-sort.
    # We require: a clear majority direction, and the number of "wrong-way"
    # steps is small (either <= _MAX_NONMONOTONIC_FRACTION of the axis, or
    # <= _MAX_NONMONOTONIC_ABSOLUTE in absolute count — whichever is more
    # permissive, so tiny test arrays and real ~3000-point axes both work).
    minor = min(n_pos, n_neg)
    if (
        minor <= _MAX_NONMONOTONIC_ABSOLUTE
        or minor / max(n_total, 1) <= _MAX_NONMONOTONIC_FRACTION
    ):
        order = np.argsort(wavelength, kind="stable")
        wl_sorted = wavelength[order]
        intensity_sorted = intensity_2d[..., order]

        # Reject duplicates: SpectroCube needs *strictly* increasing.
        if np.any(np.diff(wl_sorted) == 0):
            dup_count = int(np.sum(np.diff(wl_sorted) == 0))
            raise ValueError(
                f"Wavelength axis contains {dup_count} duplicate value(s); "
                "cannot export to SpectroCube which requires a strictly "
                "increasing wavelength coordinate."
            )
        return wl_sorted, intensity_sorted

    # Too many reversals — refuse to guess at the right ordering.
    raise ValueError(
        "Spectrum wavelength axis is not monotonic and has too many "
        f"reversals ({minor}/{n_total} = {minor / n_total:.1%}) to safely "
        "co-sort. Check the calibration order-stitching result.\n"
        f"  wavelength[:5] = {wavelength[:5]}\n"
        f"  wavelength[-5:] = {wavelength[-5:]}"
    )


def _drop_nonfinite_columns(
    wavelength: np.ndarray,
    intensity_2d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Remove wavelength columns that cannot be plotted or serialized cleanly."""
    valid = np.isfinite(wavelength) & np.all(np.isfinite(intensity_2d), axis=0)
    dropped = int(valid.size - np.count_nonzero(valid))
    if dropped == 0:
        return wavelength, intensity_2d, 0
    return wavelength[valid], intensity_2d[:, valid], dropped


def _crop_wavelength_axis(
    wavelength: np.ndarray,
    intensity_2d: np.ndarray,
    *,
    min_nm: float | None = None,
    max_nm: float | None = None,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Crop wavelength columns to the requested inclusive wavelength range."""
    if min_nm is None and max_nm is None:
        return wavelength, intensity_2d, 0

    keep = np.ones(wavelength.shape, dtype=bool)
    if min_nm is not None:
        keep &= wavelength >= float(min_nm)
    if max_nm is not None:
        keep &= wavelength <= float(max_nm)
    dropped = int(keep.size - np.count_nonzero(keep))
    if not np.any(keep):
        raise ValueError(
            "Wavelength crop removed all columns. "
            f"Requested min_nm={min_nm!r}, max_nm={max_nm!r}; "
            f"available range is {wavelength[0]:.6g}–{wavelength[-1]:.6g} nm."
        )
    if dropped == 0:
        return wavelength, intensity_2d, 0
    return wavelength[keep], intensity_2d[:, keep], dropped


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def to_spectrocube(
    spectrum,
    *,
    units: str = "counts",
    instrument_id: str | None = None,
    wavelength_medium: str = "air",
    calibration_source: str | None = "integrating sphere (echelle_spectra)",
    squeeze_single_frame: bool = True,
    drop_nonfinite_columns: bool = True,
    wavelength_min_nm: float | None = None,
    wavelength_max_nm: float | None = None,
    **extra_attrs,
) -> "SpectroCube":
    """Convert a calibrated :class:`~echelle_spectra.tools.echelle.Spectrum`
    to a :class:`spectrocube.SpectroCube` in memory.

    Parameters
    ----------
    spectrum : echelle_spectra.tools.echelle.Spectrum
        A fully-calibrated ``Spectrum`` object.  At minimum the ``wavelength``
        and ``counts`` attributes must be present.
    units : str
        Which intensity quantity to export.  One of:

        * ``"counts"``  – raw detector counts (default)
        * ``"wm"``      – absolute spectral irradiance, W m⁻² nm⁻¹
        * ``"wmsr"``    – absolute spectral radiance, W m⁻² sr⁻¹ nm⁻¹
        * ``"phmsr"``   – absolute photon radiance, ph s⁻¹ m⁻² sr⁻¹ nm⁻¹

    instrument_id : str, optional
        Short identifier stored in SpectroCube metadata.  Defaults to
        ``"echelle"`` or the ``spectrometer`` field in ``spectrum.info`` when
        present.
    wavelength_medium : str
        ``"air"`` (default) or ``"vacuum"``.
    calibration_source : str, optional
        Human-readable description of the absolute flux standard.  Only
        written to metadata when *units* is not ``"counts"``.  Defaults to
        ``"integrating sphere (echelle_spectra)"``.
    squeeze_single_frame : bool
        When *True* and the spectrum contains exactly one frame, store
        intensity as a 1-D ``(wavelength,)`` array.  When *False*, always
        keep the 2-D ``(frame, wavelength)`` shape.  Default: ``True``.
    drop_nonfinite_columns : bool
        Drop wavelength columns where the wavelength or any frame intensity is
        non-finite.  This is useful for absolute calibration columns where a
        zero sphere response would otherwise create infinities.  Dropped column
        counts are recorded in metadata.  Default: ``True``.
    wavelength_min_nm, wavelength_max_nm : float, optional
        Inclusive wavelength crop bounds in nm.  When supplied, the exported
        wavelength axis and intensity are cropped and the original range is
        recorded in metadata.
    **extra_attrs
        Any additional key/value pairs are forwarded to
        :meth:`~spectrocube.SpectroCube.from_arrays` as global NetCDF
        attributes (e.g. ``notes="my experiment"``).

    Returns
    -------
    spectrocube.SpectroCube
        A new SpectroCube object ready to validate and/or save.

    Raises
    ------
    ImportError
        If the ``spectrocube`` package is not installed.
    ValueError
        If *units* is not one of the recognised keys, or if the wavelength
        axis is non-monotonic (neither purely increasing nor purely
        decreasing).

    Examples
    --------
    >>> sc = to_spectrocube(my_spectrum, units="wm", instrument_id="BlackEchelle")
    >>> sc.save("shot42.nc")
    """
    try:
        from spectrocube import SpectroCube
    except ImportError as exc:
        raise ImportError(
            "The 'spectrocube' package is required for SpectroCube export.\n"
            "Install it with:\n"
            "    pip install spectrocube\n"
            "or, for local development:\n"
            "    pip install -e /path/to/2026-spectrocube"
        ) from exc

    if units not in _UNIT_INFO:
        raise ValueError(
            f"Unknown units {units!r}.  Choose one of: {sorted(_UNIT_INFO)}"
        )

    unit_meta = _UNIT_INFO[units]

    # --- wavelength ---
    wavelength = np.asarray(spectrum.wavelength, dtype=float)

    # --- intensity ---
    intensity_2d = np.asarray(
        getattr(spectrum, unit_meta["spectrum_attr"]), dtype=float
    )
    if intensity_2d.ndim == 1:
        intensity_2d = intensity_2d[np.newaxis, :]

    wavelength, intensity_2d = _normalize_wavelength_axis(wavelength, intensity_2d)
    original_wavelength_min_nm = float(wavelength[0])
    original_wavelength_max_nm = float(wavelength[-1])
    original_wavelength_points = int(wavelength.size)

    dropped_nonfinite_columns = 0
    if drop_nonfinite_columns:
        wavelength, intensity_2d, dropped_nonfinite_columns = _drop_nonfinite_columns(
            wavelength,
            intensity_2d,
        )
    wavelength, intensity_2d, dropped_wavelength_crop_columns = _crop_wavelength_axis(
        wavelength,
        intensity_2d,
        min_nm=wavelength_min_nm,
        max_nm=wavelength_max_nm,
    )

    n_frames = intensity_2d.shape[0]
    if squeeze_single_frame and n_frames == 1:
        intensity = intensity_2d[0]  # 1D → (wavelength,)
        dims: tuple[str, ...] | None = None
        coords: dict[str, np.ndarray] | None = None
    else:
        intensity = intensity_2d  # 2D → (frame, wavelength)
        dims = ("frame", "wavelength")
        coords = {"frame": np.arange(n_frames, dtype=float)}

    # --- instrument_id ---
    if instrument_id is None:
        info = getattr(spectrum, "info", {})
        instrument_id = (
            info.get("spectrometer")
            or _DEFAULT_INSTRUMENT_ID
        )

    # --- build metadata attrs ---
    attrs: dict[str, object] = {
        "source_package": "echelle_spectra",
        "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }

    info = getattr(spectrum, "info", {})
    if "ExposureTime" in info:
        attrs["exposure_s"] = float(info["ExposureTime"])
    if "CycleTime" in info:
        attrs["frame_interval_s"] = float(info["CycleTime"])

    fpth = getattr(spectrum, "fpth", None)
    if fpth is not None:
        attrs["source_file"] = str(fpth)

    shotnumber = getattr(spectrum, "shotnumber", None)
    if shotnumber is not None:
        attrs["shot_number"] = str(shotnumber)

    bg_frames = info.get("BackgroundFrames", [])
    if bg_frames:
        attrs["background_frames"] = str(bg_frames)

    if dropped_nonfinite_columns:
        attrs["dropped_nonfinite_wavelength_columns"] = dropped_nonfinite_columns
    if dropped_wavelength_crop_columns:
        if wavelength_min_nm is not None:
            attrs["wavelength_crop_min_nm"] = float(wavelength_min_nm)
        if wavelength_max_nm is not None:
            attrs["wavelength_crop_max_nm"] = float(wavelength_max_nm)
        attrs["original_wavelength_min_nm"] = original_wavelength_min_nm
        attrs["original_wavelength_max_nm"] = original_wavelength_max_nm
        attrs["original_wavelength_points"] = original_wavelength_points
        attrs["dropped_wavelength_crop_columns"] = dropped_wavelength_crop_columns

    calibration_folder = getattr(spectrum, "calibration_folder", None)
    if calibration_folder is not None:
        attrs["calibration_folder"] = str(calibration_folder)

    calibration_files = getattr(spectrum, "calibration_files", {})
    if calibration_files:
        if calibration_files.get("orders"):
            attrs["calibration_order_pattern_file"] = str(calibration_files["orders"])
        if calibration_files.get("wavelength"):
            attrs["wavelength_calibration_file"] = str(calibration_files["wavelength"])
        if calibration_files.get("integral"):
            attrs["absolute_calibration_integral_file"] = str(calibration_files["integral"])
        file_digests = {
            key: digest
            for key in ("orders", "wavelength", "integral")
            if (
                digest := _file_digest(
                    _resolve_calibration_file(
                        calibration_folder,
                        calibration_files.get(key),
                    )
                )
            )
        }
        if file_digests:
            attrs["calibration_file_digests_json"] = _json_attr(file_digests)

    for attr_name in (
        "calibration_order_count",
        "calibration_detector_width_px",
        "calibration_order_half_width_px",
    ):
        value = getattr(spectrum, attr_name, None)
        if value is not None:
            attrs[attr_name] = int(value)

    order_border_pixel_ranges = getattr(spectrum, "order_border_pixel_ranges", None)
    if order_border_pixel_ranges:
        attrs["order_border_pixel_ranges_json"] = _json_attr(order_border_pixel_ranges)

    order_wavelength_ranges_nm = getattr(spectrum, "order_wavelength_ranges_nm", None)
    if order_wavelength_ranges_nm:
        attrs["order_wavelength_ranges_nm_json"] = _json_attr(order_wavelength_ranges_nm)

    if unit_meta["calibration_type"] == "absolute" and calibration_source is not None:
        attrs["calibration_source"] = calibration_source

    attrs.update(extra_attrs)

    return SpectroCube.from_arrays(
        wavelength=wavelength,
        intensity=intensity,
        instrument_id=str(instrument_id),
        calibration_type=unit_meta["calibration_type"],
        intensity_units=unit_meta["intensity_units"],
        wavelength_medium=wavelength_medium,
        dims=dims,
        coords=coords,
        **attrs,
    )


def export_spectrocube(
    spectrum,
    path: str,
    *,
    units: str = "counts",
    instrument_id: str | None = None,
    wavelength_medium: str = "air",
    calibration_source: str | None = "integrating sphere (echelle_spectra)",
    squeeze_single_frame: bool = True,
    drop_nonfinite_columns: bool = True,
    wavelength_min_nm: float | None = None,
    wavelength_max_nm: float | None = None,
    validate: bool = True,
    **extra_attrs,
) -> "SpectroCube":
    """Convert a calibrated spectrum and save it to a SpectroCube NetCDF file.

    This is a convenience wrapper around :func:`to_spectrocube` followed by
    :meth:`~spectrocube.SpectroCube.save`.

    Parameters
    ----------
    spectrum : echelle_spectra.tools.echelle.Spectrum
        Calibrated spectrum to export.
    path : str
        Output file path.  Should end in ``.nc``.
    units, instrument_id, wavelength_medium, calibration_source,
    squeeze_single_frame, drop_nonfinite_columns, wavelength_min_nm,
    wavelength_max_nm, **extra_attrs
        Forwarded to :func:`to_spectrocube` — see its docstring for details.
    validate : bool
        Passed to :meth:`~spectrocube.SpectroCube.save`.  Default: ``True``.

    Returns
    -------
    spectrocube.SpectroCube
        The SpectroCube object that was saved to *path*.

    Examples
    --------
    >>> sc = export_spectrocube(my_spectrum, "output/shot42.nc", units="wm")
    >>> # Reload later:
    >>> from spectrocube import SpectroCube
    >>> sc2 = SpectroCube.load("output/shot42.nc")
    """
    sc = to_spectrocube(
        spectrum,
        units=units,
        instrument_id=instrument_id,
        wavelength_medium=wavelength_medium,
        calibration_source=calibration_source,
        squeeze_single_frame=squeeze_single_frame,
        drop_nonfinite_columns=drop_nonfinite_columns,
        wavelength_min_nm=wavelength_min_nm,
        wavelength_max_nm=wavelength_max_nm,
        **extra_attrs,
    )
    sc.save(path, validate=validate)
    return sc
