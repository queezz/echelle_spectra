"""Helpers for loading SIF files into Spectrum objects without the GUI.

This module exists so the CLI and tests can load and extract spectra without
importing any Qt/GUI machinery.  The loading logic mirrors what the GUI does
in ``LoadCalibrationsThread`` + ``LoadImageThread``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .echelle import Spectrum

_PKG_DIR = Path(__file__).parent.parent
_DEFAULT_CALIBRATION_DIR = _PKG_DIR / "resources" / "calibration_files"

# Calibration file sets.  CMOS defaults track the accepted 20250926/LHD
# wavelength-validation gate; historical tables remain available by explicitly
# passing a custom filename map to Calibrations.
DEFAULT_CAMERA_FILENAMES: dict[str, dict[str, str]] = {
    "CMOS": {
        "orders": "pattern_CMOS_20250926.txt",
        "wavelength": "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt",
        "sphr": "sphere_cmos_20240305.sif",
        "bkgr": "sphere_cmos_20240305_bkg.sif",
        "integral": "integrating_sphere.txt",
    },
    "CCD": {
        "orders": "pattern.txt",
        "wavelength": "Th_wavelength.txt",
        "sphr": "absolute_20170613_b8_0.2_v2.sif",
        "bkgr": "absolute_20170613_b8_0.2_bkg.sif",
        "integral": "integrating_sphere.txt",
    },
}

_CAMERA_FILENAMES = DEFAULT_CAMERA_FILENAMES


def _normalize_calibration_file_override(
    value: str, calibration_folder: str | Path | None = None
) -> str:
    """Resolve one filename override without letting the process CWD outrank it.

    Snapshot-sourced maps name their files by role (``pattern.txt``,
    ``wavelength.txt``, ``sphere.sif``, ...), so a same-named file in whatever
    directory the process happens to run from could silently replace a
    calibration input.  A name that exists inside *calibration_folder* is
    therefore always kept relative and resolved against that folder.  Only names
    the calibration folder does not carry fall back to working-directory
    resolution, which the shipped export config still relies on for its
    repo-relative table paths.
    """
    path = Path(value)
    if path.is_absolute():
        return str(path)
    if calibration_folder is not None and (Path(calibration_folder) / path).exists():
        return str(value)
    if path.exists():
        return str(path.resolve())
    return str(value)


def build_calibration(
    calibration_folder: str | Path | None = None,
    camera: str = "CMOS",
    calibration_files: dict[str, str] | None = None,
):
    """Load and return a ``Calibrations`` object ready for use.

    This is separated from :func:`load_spectrum` so a single calibration can
    be reused across many SIF files in batch conversion (avoids reloading
    heavy sphere images on every file).

    Parameters
    ----------
    calibration_folder : str or Path, optional
        Directory containing calibration files.  Defaults to the bundled
        ``resources/calibration_files/`` inside the installed package.
    camera : str
        Which calibration file set to use — ``"CMOS"`` (default) or ``"CCD"``.
    calibration_files : dict, optional
        Filename overrides for the selected camera calibration set.  Values may
        be relative to *calibration_folder* or absolute paths.  A relative name
        that *calibration_folder* carries always resolves inside that folder, so
        snapshot-sourced role filenames cannot be shadowed by the working
        directory.

    Returns
    -------
    echelle_spectra.tools.echelle.Calibrations
        Fully initialised calibration object (``start()`` already called).
    """
    from .echelle import Calibrations

    camera = camera.upper()
    if camera not in DEFAULT_CAMERA_FILENAMES:
        raise ValueError(f"Unknown camera {camera!r}. Choose one of: {sorted(DEFAULT_CAMERA_FILENAMES)}")

    cal_dir = (
        Path(calibration_folder) if calibration_folder is not None else _DEFAULT_CALIBRATION_DIR
    )
    if not cal_dir.is_dir():
        raise FileNotFoundError(
            f"Calibration folder not found: {cal_dir}\n"
            "Pass --calibration-dir to point to the calibration files, or ensure the\n"
            "bundled resources/calibration_files/ directory is present."
        )

    filenames = dict(DEFAULT_CAMERA_FILENAMES[camera])
    if calibration_files:
        filenames.update(
            {
                k: _normalize_calibration_file_override(str(v), cal_dir)
                for k, v in calibration_files.items()
                if v
            }
        )

    clbr = Calibrations(folder=str(cal_dir), filenames=filenames)
    clbr.start()
    return clbr


def load_spectrum(
    sif_path: str | Path,
    calibration_folder: str | Path | None = None,
    camera: str = "CMOS",
    calibration: object | None = None,
    calibration_files: dict[str, str] | None = None,
) -> "Spectrum":
    """Load and return a calibrated ``Spectrum`` from a SIF file.

    Parameters
    ----------
    sif_path : str or Path
        Path to the ``.sif`` input file.
    calibration_folder : str or Path, optional
        Directory containing calibration files.  Ignored when *calibration*
        is supplied.  Defaults to the bundled calibration files.
    camera : str
        ``"CMOS"`` (default) or ``"CCD"``.  Ignored when *calibration* is
        supplied.
    calibration : Calibrations, optional
        Pre-built calibration object.  Pass this when converting many files
        in a loop to avoid reloading the (expensive) sphere SIF for every
        input file.
    calibration_files : dict, optional
        Filename overrides passed to :func:`build_calibration` when
        *calibration* is not supplied.

    Returns
    -------
    echelle_spectra.tools.echelle.Spectrum
        Fully extracted and calibrated spectrum.

    Raises
    ------
    FileNotFoundError
        If *sif_path* does not exist.
    ValueError
        If *camera* is not recognised.

    Notes
    -----
    Loading a calibration for the first time involves reading the integrating
    sphere ``.sif`` files from the calibration folder and can take a few
    seconds.  Reuse the *calibration* argument when converting many files.
    """
    from .echelle import EchelleImage, Spectrum

    sif_path = Path(sif_path)
    if not sif_path.is_file():
        raise FileNotFoundError(f"SIF file not found: {sif_path}")

    clbr = calibration if calibration is not None else build_calibration(
        calibration_folder, camera, calibration_files=calibration_files
    )

    em = EchelleImage(str(sif_path), clbr=clbr)
    em.calculate_order_spectra()
    em.correct_order_shapes()
    em.calculate_spectra()

    sp = Spectrum(em)
    sp.fpth = sif_path
    return sp
