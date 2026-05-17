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
    squeeze_single_frame, **extra_attrs
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
        **extra_attrs,
    )
    sc.save(path, validate=validate)
    return sc
