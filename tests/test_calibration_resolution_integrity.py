"""Packet F2 — calibration resolution integrity on the real extraction path.

Every test here drives the shipped ``Calibrations``/``EchelleImage``/``Spectrum``
code.  Only detector pixels are synthetic: ``read_image`` is replaced so the
fixture does not need real SIF files, and everything downstream of it — the
cutting masks, the wavelength solution with its partial-order NaN pad, order
borders, absolute calibration, the negative-dispersion flip, SpectroCube export,
and snapshot recalibration — is the production implementation.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from echelle_spectra.recalibration import (
    RecalibrationError,
    _factor_kind,
    factor_from_snapshot,
    recalibrate_dataset,
)
from echelle_spectra.snapshot import create_snapshot
from echelle_spectra.tools import echelle as echelle_module
from echelle_spectra.tools.loader import build_calibration, load_spectrum
from echelle_spectra.tools.spectrocube_export import to_spectrocube

# Synthetic detector geometry.  The upper order climbs off the top of the sensor
# part way across the frame, which is the shipped 29-order CMOS situation that
# leaves a NaN pad in the wavelength solution.
DIMW = 40
DIMO = 55
DV = 8
ORDER_IDS = (27, 28)
PARTIAL_ORDER_FROM_PX = 28
LINE_PIXEL = 15
LINE_WAVELENGTH_NM = 436.0 - 0.1 * LINE_PIXEL
NEGATIVE_FACTOR_PIXEL = 20
EXPOSURE_S = 0.5
WAVELENGTH_SHIFT_NM = 0.05


def _order_centers() -> tuple[np.ndarray, np.ndarray]:
    columns = np.arange(DIMW)
    lower = np.full(DIMW, 20, dtype=int)
    upper = 40 + (columns * 0.25).astype(int)
    return lower, upper


def _sphere_orders() -> tuple[np.ndarray, np.ndarray]:
    """Two blaze-like humps that cross inside the orders' wavelength overlap."""

    columns = np.arange(DIMW, dtype=float)
    lower = 1000.0 - 2.0 * (columns - 20.0) ** 2
    upper = 1000.0 - 2.0 * (columns - 14.0) ** 2
    # One noise-negative sphere-minus-background column.
    lower[NEGATIVE_FACTOR_PIXEL] = 5.0
    return lower, upper


def _paint(image: np.ndarray, centers: np.ndarray, values: np.ndarray) -> None:
    for column, (center, value) in enumerate(zip(centers, values)):
        rows = np.arange(center - DV, center + DV + 1)
        inside = (rows >= 0) & (rows < DIMO)
        image[rows[inside], column] = value / (2 * DV + 1)


def _detector_frames() -> dict[str, np.ndarray]:
    lower_centers, upper_centers = _order_centers()
    lower_sphere, upper_sphere = _sphere_orders()

    sphere = np.zeros((DIMO, DIMW))
    _paint(sphere, lower_centers, lower_sphere)
    _paint(sphere, upper_centers, upper_sphere)

    background = np.zeros((DIMO, DIMW))
    _paint(background, lower_centers, np.full(DIMW, 30.0))
    _paint(background, upper_centers, np.full(DIMW, 30.0))

    shot = np.zeros((DIMO, DIMW))
    continuum = np.full(DIMW, 120.0)
    continuum[LINE_PIXEL] = 5000.0
    _paint(shot, lower_centers, continuum)
    _paint(shot, upper_centers, np.full(DIMW, 90.0))
    return {
        "sphere": sphere[np.newaxis],
        "background": background[np.newaxis],
        "shot": shot[np.newaxis],
    }


@pytest.fixture
def detector(monkeypatch: pytest.MonkeyPatch) -> dict[str, np.ndarray]:
    """Serve synthetic frames wherever the pipeline would read a SIF file."""

    frames = _detector_frames()

    def read_image(fpth, spec="black", crop=(0, -1), exptime=1):
        name = Path(str(fpth)).name.casefold()
        if name.endswith("_bg.sif"):
            images = frames["background"]
        elif name.endswith("sphere.sif"):
            images = frames["sphere"]
        else:
            images = frames["shot"]
        info = {
            "NumberOfFrames": int(images.shape[0]),
            "xbin": 1,
            "ybin": 1,
            "size": np.array([DIMW, DIMO]),
            "ExposureTime": EXPOSURE_S,
            "CycleTime": 1.0,
        }
        return images.copy(), info

    monkeypatch.setattr(echelle_module, "read_image", read_image)
    return frames


def _write_sources(root: Path, *, shift: float, integral_scale: float) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    lower_centers, upper_centers = _order_centers()
    pattern = root / "pattern.txt"
    np.savetxt(pattern, np.column_stack([lower_centers, upper_centers]), fmt="%d")

    rows = [
        (ORDER_IDS[0], pixel, 436.0 - 0.1 * pixel + shift) for pixel in (2, 10, 18, 26, 34)
    ]
    rows += [(ORDER_IDS[1], pixel, 433.0 - 0.1 * pixel + shift) for pixel in (1, 8, 15, 22)]
    wavelength = root / "wavelength.txt"
    wavelength.write_text(
        "# synthetic lamp identification\n"
        "# order from to center wavelength\n"
        + "".join(
            f"{order:d} {pixel - 1:d} {pixel + 1:d} {pixel:d} {value:.6f}\n"
            for order, pixel, value in rows
        ),
        encoding="utf-8",
    )

    integral = root / "integral.txt"
    micrometres = np.linspace(0.420, 0.450, 9)
    radiance = np.full(micrometres.size, integral_scale * 1.0e-2)
    np.savetxt(integral, np.column_stack([micrometres, radiance]), fmt="%.8f")

    sphere = root / "sphere.sif"
    sphere.write_bytes(b"synthetic integrating-sphere frame")
    background = root / "sphere_bg.sif"
    background.write_bytes(b"synthetic integrating-sphere background frame")
    return {
        "pattern": pattern,
        "wavelength": wavelength,
        "integral": integral,
        "sphere": sphere,
        "sphere_background": background,
    }


def _snapshot(
    root: Path,
    snapshot_id: str,
    *,
    shift: float = 0.0,
    integral_scale: float = 1.0,
):
    files = _write_sources(
        root / f"sources-{snapshot_id}", shift=shift, integral_scale=integral_scale
    )
    return create_snapshot(
        root / "calibrations",
        snapshot_id=snapshot_id,
        detector="cmos",
        files=files,
        lamps=["ThAr"],
        validity={"shot_from": 1, "shot_to": 999999},
    )


def _shot_file(root: Path) -> Path:
    path = root / "193778_Echelle.SIF"
    path.write_bytes(b"synthetic shot frame")
    return path


def _cwd_decoys(root: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Plant calibration files with the snapshot's role names in the CWD."""

    decoy = root / "working-directory"
    decoy.mkdir()
    lower_centers, upper_centers = _order_centers()
    np.savetxt(
        decoy / "pattern.txt",
        np.column_stack([lower_centers + 3, upper_centers]),
        fmt="%d",
    )
    rows = [(ORDER_IDS[0], pixel, 440.0 - 0.1 * pixel) for pixel in (2, 10, 18, 26, 34)]
    rows += [(ORDER_IDS[1], pixel, 437.0 - 0.1 * pixel) for pixel in (1, 8, 15, 22)]
    (decoy / "wavelength.txt").write_text(
        "# decoy identification\n"
        "# order from to center wavelength\n"
        + "".join(
            f"{order:d} {pixel - 1:d} {pixel + 1:d} {pixel:d} {value:.6f}\n"
            for order, pixel, value in rows
        ),
        encoding="utf-8",
    )
    (decoy / "integral.txt").write_text("0.420 1.0\n0.430 1.0\n0.440 1.0\n0.450 1.0\n")
    monkeypatch.chdir(decoy)
    return decoy


def _export(spectrum, snapshot, **kws):
    return to_spectrocube(
        spectrum,
        snapshot=snapshot,
        units="wmsr",
        squeeze_single_frame=False,
        **kws,
    )


class TestSnapshotFilenamesIgnoreTheWorkingDirectory:
    def test_snapshot_roles_resolve_against_the_snapshot_not_the_cwd(
        self, tmp_path: Path, detector, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        _cwd_decoys(tmp_path, monkeypatch)

        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )

        lower_centers, _ = _order_centers()
        np.testing.assert_array_equal(calibration.pattern[:, 0], lower_centers)
        assert calibration.order_wavel[0][0] == pytest.approx(436.0, abs=1e-9)
        assert calibration.wavelength[np.isfinite(calibration.wavelength)].max() < 437.0

    def test_repo_relative_override_still_resolves_from_the_working_directory(
        self, tmp_path: Path, detector, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The shipped export config names its tables repo-relative while its
        # sphere frames stay relative to the calibration directory.
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        elsewhere = tmp_path / "tables"
        elsewhere.mkdir()
        table = elsewhere / "aligned.txt"
        table.write_text(
            (snapshot.root / "wavelength.txt").read_text(encoding="utf-8"), encoding="utf-8"
        )
        monkeypatch.chdir(tmp_path)

        files = dict(snapshot.calibration_files())
        files["wavelength"] = "tables/aligned.txt"
        calibration = build_calibration(snapshot.root, "CMOS", calibration_files=files)

        assert Path(calibration.filenames["wavelength"]) == table.resolve()
        assert calibration.filenames["sphr"] == "sphere.sif"
        assert calibration.order_wavel[0][0] == pytest.approx(436.0, abs=1e-9)


class _Field:
    """The attribute carrier ``_factor_kind`` reads on a saved factor."""

    def __init__(self, attrs: dict[str, str]) -> None:
        self.attrs = attrs


class TestAbsoluteKindIsDeclaredNotGuessed:
    def test_declared_kind_wins_and_units_stay_the_documented_fallback(self) -> None:
        assert _factor_kind(_Field({"absolute_kind": "phmsr", "units": "W/m2/nm"})) == "phmsr"
        assert _factor_kind(_Field({"units": "W/m2/nm/sr per (counts/s)"})) == "wmsr"
        assert _factor_kind(_Field({"units": "W/m2/nm per (counts/s)"})) == "wm"
        assert _factor_kind(_Field({"units": "ph/s/nm/sr per (counts/s)"})) == "phmsr"

    def test_unusable_units_are_refused_by_name(self) -> None:
        with pytest.raises(RecalibrationError, match="absolute_kind"):
            _factor_kind(_Field({"units": ""}))
        with pytest.raises(RecalibrationError, match="absolute_kind"):
            _factor_kind(_Field({"units": "spectral radiance"}))
        with pytest.raises(RecalibrationError, match="unknown absolute_kind"):
            _factor_kind(_Field({"absolute_kind": "photons", "units": "W/m2/nm"}))


class TestRealPipelineCube:
    """One end-to-end run: extraction, flip, export, and default recalibration."""

    def test_calibration_carries_a_partial_order_and_negative_dispersion(
        self, tmp_path: Path, detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )

        assert calibration.direction < 0
        np.testing.assert_array_equal(calibration.orders_bad_shape, [1])
        np.testing.assert_array_equal(calibration.orders_bad_froms, [PARTIAL_ORDER_FROM_PX])
        n_full = int(calibration.wavelength.size)
        n_keep = int(np.count_nonzero(np.isfinite(calibration.wavelength)))
        assert n_keep < n_full
        assert calibration.absolute["wmsr"].size == n_keep
        assert int(np.count_nonzero(calibration.absolute["wmsr"] <= 0)) == 1

    def test_aligned_fields_survive_flip_drop_and_crop(self, tmp_path: Path, detector) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)

        assert spectrum.direction < 0
        assert spectrum.wavelength[0] < spectrum.wavelength[-1]

        cube = _export(spectrum, snapshot, wavelength_min_nm=431.0)
        wavelength = np.asarray(cube.wavelength)
        pixels = np.asarray(cube.ds["detector_pixel"].values)
        orders = np.asarray(cube.ds["echelle_order"].values)
        factor = np.asarray(cube.ds["applied_absolute_calibration_factor"].values)

        assert np.all(np.diff(wavelength) > 0)
        assert wavelength.size == pixels.size == orders.size == factor.size
        assert cube.intensity.shape[-1] == wavelength.size

        peak = int(np.argmax(cube.intensity[0]))
        assert int(pixels[peak]) == LINE_PIXEL
        assert int(orders[peak]) == ORDER_IDS[0]
        assert wavelength[peak] == pytest.approx(LINE_WAVELENGTH_NM, abs=1e-9)

        # The polynomial round-trip guard inside the writer accepted this cube;
        # verify the same reconstruction the reader would perform.
        payload = json.loads(cube.ds.attrs["wavelength_polynomials_json"])
        assert {record["order"] for record in payload["orders"]} == set(np.unique(orders).tolist())
        for record in payload["orders"]:
            selected = orders == record["order"]
            np.testing.assert_allclose(
                np.polyval(record["coefficients"], pixels[selected]),
                wavelength[selected],
                rtol=0.0,
                atol=5e-10,
            )

        assert wavelength.min() >= 431.0
        assert cube.ds.attrs["dropped_wavelength_crop_columns"] > 0
        assert cube.validate().ok, str(cube.validate())

    def test_one_negative_sphere_column_is_dropped_and_counted(
        self, tmp_path: Path, detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)

        cube = _export(spectrum, snapshot)
        pixels = np.asarray(cube.ds["detector_pixel"].values)
        orders = np.asarray(cube.ds["echelle_order"].values)
        factor = np.asarray(cube.ds["applied_absolute_calibration_factor"].values)

        assert cube.ds.attrs["dropped_nonpositive_factor_columns"] == 1
        assert not np.any((orders == ORDER_IDS[0]) & (pixels == NEGATIVE_FACTOR_PIXEL))
        assert np.all(factor > 0)
        assert cube.ds["applied_absolute_calibration_factor"].attrs["absolute_kind"] == "wmsr"

        with pytest.raises(ValueError, match="strictly positive"):
            _export(spectrum, snapshot, drop_nonfinite_columns=False)

    def test_default_recalibration_moves_wavelengths_and_rescales_intensity(
        self, tmp_path: Path, detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        revision = _snapshot(
            tmp_path,
            "20250202_cmos",
            shift=WAVELENGTH_SHIFT_NM,
            integral_scale=2.0,
        )
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)
        cube = _export(spectrum, snapshot, wavelength_min_nm=431.0)

        new_factor = factor_from_snapshot(cube.ds, revision)
        old_factor = np.asarray(cube.ds["applied_absolute_calibration_factor"].values)
        np.testing.assert_allclose(new_factor, old_factor * 2.0, rtol=1e-12)

        revised, event = recalibrate_dataset(cube.ds, revision)

        np.testing.assert_allclose(
            np.asarray(revised.coords["wavelength"].values),
            np.asarray(cube.wavelength) + WAVELENGTH_SHIFT_NM,
            rtol=0.0,
            atol=1e-9,
        )
        np.testing.assert_allclose(
            np.asarray(revised["intensity"].values),
            np.asarray(cube.intensity) * 2.0,
            rtol=1e-12,
        )
        np.testing.assert_array_equal(
            np.asarray(revised.coords["detector_pixel"].values),
            np.asarray(cube.ds["detector_pixel"].values),
        )
        assert event["changes"] == ["wavelength", "absolute-factor"]
        assert event["new_snapshot_id"] == revision.snapshot_id
        assert "dropped_nonpositive_factor_columns" not in event

        from spectrocube import SpectroCube

        assert SpectroCube.from_dataset(revised).validate().ok

    def test_recalibration_drops_a_non_positive_replacement_sample(
        self, tmp_path: Path, detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        revision = _snapshot(tmp_path, "20250202_cmos", shift=WAVELENGTH_SHIFT_NM)
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)
        cube = _export(spectrum, snapshot, wavelength_min_nm=431.0)

        replacement = factor_from_snapshot(cube.ds, revision)
        replacement[3] = -replacement[3]
        revised, event = recalibrate_dataset(cube.ds, revision, new_factor=replacement)

        assert revised.sizes["wavelength"] == cube.wavelength.size - 1
        assert event["dropped_nonpositive_factor_columns"] == 1
        # The attribute counts every column this cube has lost to a non-positive
        # factor: one at export, one more here.
        assert revised.attrs["dropped_nonpositive_factor_columns"] == 2
        assert np.all(
            np.asarray(revised["applied_absolute_calibration_factor"].values) > 0
        )
        assert np.all(np.isfinite(np.asarray(revised["intensity"].values)))


#: A column where the sphere measured exactly its own background: net response
#: zero, so the absolute factor there is undefined rather than merely negative.
ZERO_RESPONSE_PIXEL = 24


@pytest.fixture
def zero_response_detector(monkeypatch: pytest.MonkeyPatch) -> dict[str, np.ndarray]:
    """The same synthetic pipeline with one dead integrating-sphere column."""

    frames = _detector_frames()
    lower_centers, _upper_centers = _order_centers()
    center = lower_centers[ZERO_RESPONSE_PIXEL]
    rows = np.arange(center - DV, center + DV + 1)
    inside = (rows >= 0) & (rows < DIMO)
    # Give one sphere column exactly the background's own value, pixel for
    # pixel, so their difference there is zero rather than merely small.
    frames["sphere"][0][rows[inside], ZERO_RESPONSE_PIXEL] = 30.0 / (2 * DV + 1)

    def read_image(fpth, spec="black", crop=(0, -1), exptime=1):
        name = Path(str(fpth)).name.casefold()
        if name.endswith("_bg.sif"):
            images = frames["background"]
        elif name.endswith("sphere.sif"):
            images = frames["sphere"]
        else:
            images = frames["shot"]
        info = {
            "NumberOfFrames": int(images.shape[0]),
            "xbin": 1,
            "ybin": 1,
            "size": np.array([DIMW, DIMO]),
            "ExposureTime": EXPOSURE_S,
            "CycleTime": 1.0,
        }
        return images.copy(), info

    monkeypatch.setattr(echelle_module, "read_image", read_image)
    return frames


class TestADeadSphereColumnIsDroppedNotWarned:
    """Packet F13 — a zero sphere response is a dropped column, not a warning."""

    def test_absolute_calibration_makes_it_nan_without_a_runtime_warning(
        self, tmp_path: Path, zero_response_detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            calibration = build_calibration(
                snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
            )

        factor = calibration.absolute["wmsr"]
        assert calibration.zero_response_columns == 1
        assert int(np.count_nonzero(~np.isfinite(factor))) == 1
        # The noise-negative column the F2 fixture carries keeps its negative
        # factor: only the undefined one became NaN.
        assert int(np.count_nonzero(factor <= 0)) == 1

    def test_the_dead_column_is_dropped_and_counted_at_export(
        self, tmp_path: Path, zero_response_detector
    ) -> None:
        snapshot = _snapshot(tmp_path, "20250101_cmos")
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)

        cube = _export(spectrum, snapshot)

        pixels = np.asarray(cube.ds["detector_pixel"].values)
        orders = np.asarray(cube.ds["echelle_order"].values)
        factor = np.asarray(cube.ds["applied_absolute_calibration_factor"].values)
        assert cube.ds.attrs["dropped_nonfinite_wavelength_columns"] >= 1
        assert not np.any((orders == ORDER_IDS[0]) & (pixels == ZERO_RESPONSE_PIXEL))
        assert not np.any((orders == ORDER_IDS[0]) & (pixels == NEGATIVE_FACTOR_PIXEL))
        assert np.all(np.isfinite(factor)) and np.all(factor > 0)
        assert cube.validate().ok, str(cube.validate())
