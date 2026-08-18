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
    _replace_factor,
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


def _write_sources(
    root: Path,
    *,
    shift: float,
    integral_scale: float,
    order_offsets: tuple[float, float] = (0.0, 0.0),
    sphere_marker: str = "",
) -> dict[str, Path]:
    root.mkdir(parents=True, exist_ok=True)
    lower_centers, upper_centers = _order_centers()
    pattern = root / "pattern.txt"
    np.savetxt(pattern, np.column_stack([lower_centers, upper_centers]), fmt="%d")

    rows = [
        (ORDER_IDS[0], pixel, 436.0 - 0.1 * pixel + shift + order_offsets[0])
        for pixel in (2, 10, 18, 26, 34)
    ]
    rows += [
        (ORDER_IDS[1], pixel, 433.0 - 0.1 * pixel + shift + order_offsets[1])
        for pixel in (1, 8, 15, 22)
    ]
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
    # A distinct marker makes this era's sphere pair a distinct artifact digest,
    # which is what a real era change always is.
    sphere.write_bytes(f"synthetic integrating-sphere frame{sphere_marker}".encode())
    background = root / "sphere_bg.sif"
    background.write_bytes(
        f"synthetic integrating-sphere background frame{sphere_marker}".encode()
    )
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
    order_offsets: tuple[float, float] = (0.0, 0.0),
    sphere_marker: str = "",
):
    files = _write_sources(
        root / f"sources-{snapshot_id}",
        shift=shift,
        integral_scale=integral_scale,
        order_offsets=order_offsets,
        sphere_marker=sphere_marker,
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


#: Packet F14 -- the two eras a cross-era revision has to join.  The cube is born
#: under the old one; the new one is a genuine era change, so it brings its own
#: lamp identification *and* its own sphere pair.
OLD_ERA_ID = "20190314_cmos"
NEW_ERA_ID = "20250926_cmos"
#: The new era lands its upper order this much redder than the lower one.  That
#: relative move -- not the common shift, which slides both orders together --
#: is what walks the order border, and the border decides which detector columns
#: a snapshot publishes an absolute factor for.
NEW_ERA_ORDER_OFFSET_NM = 0.2
#: A column where the new era's sphere measured exactly its own background.  The
#: sphere never shone there, so no era change can produce a factor for it.
NEW_ERA_DEAD_SPHERE_PIXEL = 12


def _era_frames() -> dict[str, dict[str, np.ndarray]]:
    """Detector frames per era: the new one has one dead sphere column."""

    old = _detector_frames()
    new = _detector_frames()
    lower_centers, _ = _order_centers()
    center = lower_centers[NEW_ERA_DEAD_SPHERE_PIXEL]
    rows = np.arange(center - DV, center + DV + 1)
    inside = (rows >= 0) & (rows < DIMO)
    new["sphere"][0][rows[inside], NEW_ERA_DEAD_SPHERE_PIXEL] = 30.0 / (2 * DV + 1)
    return {OLD_ERA_ID: old, NEW_ERA_ID: new}


@pytest.fixture
def two_era_detector(monkeypatch: pytest.MonkeyPatch) -> dict[str, dict[str, np.ndarray]]:
    """Serve each snapshot the sphere its own era measured.

    A snapshot copies its artifacts into a folder named for its id, so the path
    the pipeline reads names the era without the fixture having to thread it
    through ``build_calibration``.
    """

    eras = _era_frames()

    def read_image(fpth, spec="black", crop=(0, -1), exptime=1):
        path = Path(str(fpth))
        frames = eras.get(path.parent.name, eras[OLD_ERA_ID])
        name = path.name.casefold()
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
    return eras


def _published_samples(snapshot) -> set[tuple[int, int]]:
    """The (order, pixel) pairs a snapshot publishes an absolute factor for.

    This is the sample set the pre-F14 ``factor_from_snapshot`` demanded a cube
    be a subset of: the snapshot's order borders intersected with its
    finite-wavelength keep mask.
    """

    calibration = build_calibration(
        snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
    )
    keep = np.isfinite(np.asarray(calibration.wavelength, dtype=float))
    pixels = np.broadcast_to(np.arange(calibration.DIMW), calibration.order_borders.shape)[
        calibration.order_borders
    ]
    orders = np.broadcast_to(
        np.asarray(calibration.order_ids)[:, None], calibration.order_borders.shape
    )[calibration.order_borders]
    return set(zip(orders[keep].tolist(), pixels[keep].tolist()))


def _old_era_cube(tmp_path: Path):
    snapshot = _snapshot(tmp_path, OLD_ERA_ID, sphere_marker=" 2019")
    calibration = build_calibration(
        snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
    )
    spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)
    return snapshot, _export(spectrum, snapshot)


def _write_alignment(snapshot, *, rms_px: float) -> None:
    """Give a snapshot the alignment record the wavelength accuracy reads."""

    (snapshot.root / "alignment.toml").write_text(
        'instrument_id = "cmos"\n'
        'base_wavelength_file = "wavelength.txt"\n'
        "n_lines = 9\n"
        f"rms_px = {rms_px!r}\n"
        "[transform]\ndx_px = 0.0\ndy_px = 0.0\ntheta_rad = 0.0\n",
        encoding="utf-8",
    )


def _new_era_snapshot(tmp_path: Path):
    return _snapshot(
        tmp_path,
        NEW_ERA_ID,
        shift=WAVELENGTH_SHIFT_NM,
        order_offsets=(0.0, NEW_ERA_ORDER_OFFSET_NM),
        sphere_marker=" 2025",
    )


class TestCrossEraRecalibration:
    """Packet F14 -- a cube born in one era, recalibrated onto another."""

    def test_a_moved_border_defeats_the_matching_samples_contract(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        # The refusal this packet exists to remove: the new era publishes a
        # different set of (order, pixel) samples, so a contract that demands
        # every cube sample appear in it fails for exactly the era change the
        # revision is for.
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        cube_samples = set(
            zip(
                np.asarray(cube.ds["echelle_order"].values, dtype=int).tolist(),
                np.asarray(cube.ds["detector_pixel"].values, dtype=int).tolist(),
            )
        )
        unmatched = sorted(cube_samples - _published_samples(new))
        assert unmatched, "fixture no longer moves the order border"
        assert (ORDER_IDS[1], 2) in unmatched

    def test_the_cube_grid_recovers_every_border_moved_sample(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        factor = factor_from_snapshot(cube.ds, new)
        orders = np.asarray(cube.ds["echelle_order"].values, dtype=int)
        pixels = np.asarray(cube.ds["detector_pixel"].values, dtype=int)
        published = _published_samples(new)
        recovered = np.asarray(
            [(int(o), int(p)) not in published for o, p in zip(orders, pixels)]
        )
        assert np.any(recovered)

        # A sample the new era's borders did not publish is still a measurement:
        # its own sphere pair minus background, through its own integral curve,
        # at that very detector column.
        calibration = build_calibration(
            new.root, "CMOS", calibration_files=new.calibration_files()
        )
        rows = {int(value): index for index, value in enumerate(calibration.order_ids)}
        expected = calibration.absolute_on_detector_grid()["wmsr"]
        for order, pixel in zip(orders[recovered], pixels[recovered]):
            if pixel == NEW_ERA_DEAD_SPHERE_PIXEL:
                continue
            assert np.isfinite(expected[rows[int(order)], int(pixel)])
        np.testing.assert_array_equal(
            factor[recovered & (pixels != NEW_ERA_DEAD_SPHERE_PIXEL)],
            [
                expected[rows[int(order)], int(pixel)]
                for order, pixel in zip(orders, pixels)
                if (int(order), int(pixel)) not in published
                and int(pixel) != NEW_ERA_DEAD_SPHERE_PIXEL
            ],
        )

    def test_a_full_cross_era_revision_applies_both_deltas(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        old, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        revised, event = recalibrate_dataset(cube.ds, new)

        assert event["changes"] == ["wavelength", "absolute-factor"]
        assert event["old_snapshot_id"] == old.snapshot_id
        assert event["new_snapshot_id"] == new.snapshot_id
        assert event["old_snapshot_manifest_json"] == old.provenance_attrs()[
            "snapshot_manifest_json"
        ]
        assert event["new_snapshot_manifest_json"] == new.provenance_attrs()[
            "snapshot_manifest_json"
        ]
        assert revised.attrs["snapshot_id"] == new.snapshot_id
        assert (
            revised.attrs["calibration_file_digests_json"]
            == new.provenance_attrs()["calibration_file_digests_json"]
        )

        # The new era's wavelength solution, evaluated on the cube's own pixels.
        payload = json.loads(revised.attrs["wavelength_polynomials_json"])
        assert payload["snapshot_id"] == new.snapshot_id
        pixels = np.asarray(revised.coords["detector_pixel"].values, dtype=float)
        orders = np.asarray(revised.coords["echelle_order"].values, dtype=int)
        for record in payload["orders"]:
            selected = orders == record["order"]
            np.testing.assert_allclose(
                np.polyval(record["coefficients"], pixels[selected]),
                np.asarray(revised.coords["wavelength"].values)[selected],
                rtol=0.0,
                atol=5e-10,
            )
        # The upper order moved by the era shift plus its own offset, the lower
        # order by the era shift alone: a real solution change, not a slide.
        for order, expected in (
            (ORDER_IDS[0], WAVELENGTH_SHIFT_NM),
            (ORDER_IDS[1], WAVELENGTH_SHIFT_NM + NEW_ERA_ORDER_OFFSET_NM),
        ):
            record = next(item for item in payload["orders"] if item["order"] == order)
            assert np.polyval(record["coefficients"], 5.0) - np.polyval(
                json.loads(cube.ds.attrs["wavelength_polynomials_json"])["orders"][
                    [item["order"] for item in payload["orders"]].index(order)
                ]["coefficients"],
                5.0,
            ) == pytest.approx(expected, abs=1e-9)

        factor = np.asarray(revised["applied_absolute_calibration_factor"].values)
        assert np.all(np.isfinite(factor)) and np.all(factor > 0)
        assert np.all(np.diff(np.asarray(revised.coords["wavelength"].values)) > 0)

        from spectrocube import SpectroCube

        assert SpectroCube.from_dataset(revised).validate().ok

    def test_extraction_geometry_is_untouched_by_the_revision(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)
        revised, _ = recalibrate_dataset(cube.ds, new)

        # Every surviving sample still names the detector column it was cut
        # from; only the wavelength and the factor attached to it changed.
        kept = set(
            zip(
                np.asarray(revised.coords["echelle_order"].values, dtype=int).tolist(),
                np.asarray(revised.coords["detector_pixel"].values, dtype=int).tolist(),
            )
        )
        original = set(
            zip(
                np.asarray(cube.ds["echelle_order"].values, dtype=int).tolist(),
                np.asarray(cube.ds["detector_pixel"].values, dtype=int).tolist(),
            )
        )
        assert kept <= original
        assert (
            revised.attrs["order_border_pixel_ranges_json"]
            == cube.ds.attrs["order_border_pixel_ranges_json"]
        )

    def test_a_dead_new_era_sphere_column_is_uncovered_dropped_and_counted(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        factor = factor_from_snapshot(cube.ds, new)
        orders = np.asarray(cube.ds["echelle_order"].values, dtype=int)
        pixels = np.asarray(cube.ds["detector_pixel"].values, dtype=int)
        dead = (orders == ORDER_IDS[0]) & (pixels == NEW_ERA_DEAD_SPHERE_PIXEL)
        assert np.count_nonzero(dead) == 1
        # The old era measured that column perfectly well; the new one did not.
        # No amount of smoothness licenses inventing a response there.
        assert np.all(np.isfinite(np.asarray(cube.ds["applied_absolute_calibration_factor"])))
        assert not np.isfinite(factor[dead][0])

        revised, event = recalibrate_dataset(cube.ds, new)

        assert event["dropped_uncovered_factor_columns"] == 1
        assert event["dropped_nonpositive_factor_columns"] == 1
        assert revised.attrs["dropped_uncovered_factor_columns"] == 1
        # The cumulative attribute already carried the column export dropped.
        assert revised.attrs["dropped_nonpositive_factor_columns"] == 2
        assert revised.sizes["wavelength"] == cube.wavelength.size - 1
        kept = set(
            zip(
                np.asarray(revised.coords["echelle_order"].values, dtype=int).tolist(),
                np.asarray(revised.coords["detector_pixel"].values, dtype=int).tolist(),
            )
        )
        assert (ORDER_IDS[0], NEW_ERA_DEAD_SPHERE_PIXEL) not in kept

    def test_wavelength_only_across_eras_names_the_full_revision(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        with pytest.raises(RecalibrationError) as raised:
            recalibrate_dataset(cube.ds, new, update_factor=False)
        message = str(raised.value)
        assert "sphere" in message and "sphere_background" in message
        assert "--wavelength-only" in message and "full recalibration" in message

    def test_wavelength_only_still_serves_a_same_sphere_target(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        # A refinement snapshot: a new lamp solution over the same sphere pair.
        refined = _snapshot(
            tmp_path, "20190401_cmos", shift=WAVELENGTH_SHIFT_NM, sphere_marker=" 2019"
        )

        revised, event = recalibrate_dataset(cube.ds, refined, update_factor=False)

        assert event["changes"] == ["wavelength"]
        np.testing.assert_array_equal(
            np.asarray(revised["applied_absolute_calibration_factor"].values),
            np.asarray(cube.ds["applied_absolute_calibration_factor"].values),
        )
        np.testing.assert_array_equal(
            np.asarray(revised["intensity"].values), np.asarray(cube.intensity)
        )
        np.testing.assert_allclose(
            np.asarray(revised.coords["wavelength"].values),
            np.asarray(cube.wavelength) + WAVELENGTH_SHIFT_NM,
            rtol=0.0,
            atol=1e-9,
        )

    def test_recalibrating_onto_its_own_snapshot_is_a_bit_exact_round_trip(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        snapshot, cube = _old_era_cube(tmp_path)

        revised, event = recalibrate_dataset(cube.ds, snapshot)

        assert event["changes"] == ["wavelength", "absolute-factor"]
        assert "dropped_nonpositive_factor_columns" not in event
        assert "dropped_uncovered_factor_columns" not in event
        for name in ("wavelength", "detector_pixel", "echelle_order"):
            np.testing.assert_array_equal(
                np.asarray(revised.coords[name].values), np.asarray(cube.ds[name].values)
            )
        for name in ("intensity", "applied_absolute_calibration_factor"):
            np.testing.assert_array_equal(
                np.asarray(revised[name].values), np.asarray(cube.ds[name].values)
            )
        assert (
            revised.attrs["order_wavelength_ranges_nm_json"]
            == cube.ds.attrs["order_wavelength_ranges_nm_json"]
        )

    def test_an_unchanged_factor_leaves_every_intensity_bit_alone(self) -> None:
        """The round trip at the scale where floating point can break it.

        The shipped cube is 2304 columns wide across 29 orders.  Dividing the
        old factor out and multiplying the same value back in loses the last bit
        on roughly one sample in ten at that width, so the round trip has to be
        an identity by construction rather than by luck: an unchanged factor
        gives a ratio of exactly one.
        """

        import xarray as xr

        generator = np.random.default_rng(20250926)
        factor = generator.random(29 * 2304) * 3.0e-9 + 1.0e-12
        intensity = generator.random(29 * 2304) * 1.0e-3
        dataset = xr.Dataset(
            {"intensity": (("frame", "wavelength"), intensity[None, :])},
            coords={"frame": [0], "wavelength": np.arange(factor.size, dtype=float)},
        )
        dataset["applied_absolute_calibration_factor"] = (
            ("wavelength",),
            factor,
            {"absolute_kind": "wmsr"},
        )
        # Dividing out and multiplying back in would move these samples.
        assert not np.array_equal((intensity / factor) * factor, intensity)

        revised, dropped, uncovered = _replace_factor(dataset, dataset, None, factor.copy())

        assert (dropped, uncovered) == (0, 0)
        np.testing.assert_array_equal(
            np.asarray(revised["intensity"].values)[0], intensity
        )


def _decoded(value):
    """Read a JSON attribute string back into data, or return nothing."""

    stripped = value.strip()
    if stripped[:1] not in "[{":
        return None
    try:
        return json.loads(stripped)
    except ValueError:
        return None


def _numbers(value):
    """Every number reachable inside an attribute value, JSON strings included."""

    if isinstance(value, bool):
        return
    if isinstance(value, (int, float, np.integer, np.floating)):
        yield float(value)
    elif isinstance(value, str):
        decoded = _decoded(value)
        if decoded is not None:
            yield from _numbers(decoded)
    elif isinstance(value, dict):
        for item in value.values():
            yield from _numbers(item)
    elif isinstance(value, (list, tuple, np.ndarray)):
        for item in value:
            yield from _numbers(item)


class TestNoAttributeOutlivesItsWavelengthSolution:
    """Packet F14 -- a revised cube must not advertise the old era's numbers."""

    #: Provenance that is *supposed* to quote the superseded solution.
    HISTORY_ATTRS = frozenset(
        {
            "recalibration_history_json",
            "previous_snapshot_manifest_json",
            "snapshot_manifest_json",
            "calibration_file_digests_json",
        }
    )

    def test_no_attribute_value_still_matches_the_old_solution(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)

        # Every wavelength the old solution derived and published as an
        # attribute.  After the revision none of them may still be quoted.
        stale = set()
        for record in json.loads(cube.ds.attrs["order_wavelength_ranges_nm_json"]):
            stale.update((float(record["min_nm"]), float(record["max_nm"])))
        for name in ("original_wavelength_min_nm", "original_wavelength_max_nm"):
            if name in cube.ds.attrs:
                stale.add(float(cube.ds.attrs[name]))
        assert len(stale) >= 2

        revised, _ = recalibrate_dataset(cube.ds, new)

        offenders = {
            name: number
            for name, value in revised.attrs.items()
            if name not in self.HISTORY_ATTRS
            for number in _numbers(value)
            if any(abs(number - old) <= 1e-12 for old in stale)
        }
        assert not offenders, f"attributes still carrying the old solution: {offenders}"

    def test_each_order_range_is_restated_from_the_new_solution(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        _, cube = _old_era_cube(tmp_path)
        new = _new_era_snapshot(tmp_path)
        revised, _ = recalibrate_dataset(cube.ds, new)

        before = {
            int(record["order"]): record
            for record in json.loads(cube.ds.attrs["order_wavelength_ranges_nm_json"])
        }
        after = {
            int(record["order"]): record
            for record in json.loads(revised.attrs["order_wavelength_ranges_nm_json"])
        }
        assert set(after) == set(before)
        for order, expected_shift in (
            (ORDER_IDS[0], WAVELENGTH_SHIFT_NM),
            (ORDER_IDS[1], WAVELENGTH_SHIFT_NM + NEW_ERA_ORDER_OFFSET_NM),
        ):
            # The order still spans the same detector pixels -- extraction
            # geometry did not move -- but every wavelength it claims did.
            assert after[order]["n_px"] == before[order]["n_px"]
            for bound in ("min_nm", "max_nm"):
                assert after[order][bound] - before[order][bound] == pytest.approx(
                    expected_shift, abs=1e-9
                )

    def test_the_pre_crop_axis_bounds_are_restated_too(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        snapshot = _snapshot(tmp_path, OLD_ERA_ID, sphere_marker=" 2019")
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)
        cube = _export(spectrum, snapshot, wavelength_min_nm=431.0)
        assert "original_wavelength_min_nm" in cube.ds.attrs

        revised, _ = recalibrate_dataset(cube.ds, _new_era_snapshot(tmp_path))

        # The bounds are rebuilt over the same border-and-pad pixel set the
        # export used, so the lower order's move is the common era shift.
        assert float(revised.attrs["original_wavelength_max_nm"]) - float(
            cube.ds.attrs["original_wavelength_max_nm"]
        ) == pytest.approx(WAVELENGTH_SHIFT_NM, abs=1e-9)
        assert float(revised.attrs["original_wavelength_min_nm"]) != pytest.approx(
            float(cube.ds.attrs["original_wavelength_min_nm"]), abs=1e-9
        )
        # A user-requested crop threshold is an input, not a derived wavelength.
        assert float(revised.attrs["wavelength_crop_min_nm"]) == 431.0

    def test_the_wavelength_accuracy_is_restated_or_withdrawn(
        self, tmp_path: Path, two_era_detector
    ) -> None:
        snapshot = _snapshot(tmp_path, OLD_ERA_ID, sphere_marker=" 2019")
        _write_alignment(snapshot, rms_px=0.4)
        calibration = build_calibration(
            snapshot.root, "CMOS", calibration_files=snapshot.calibration_files()
        )
        spectrum = load_spectrum(_shot_file(tmp_path), calibration=calibration)
        cube = _export(spectrum, snapshot)
        # 0.1 nm per pixel of dispersion, so the RMS lands as 0.04 nm.
        assert float(cube.ds.attrs["wavelength_accuracy_nm"]) == pytest.approx(0.04, abs=1e-6)

        aligned_era = _new_era_snapshot(tmp_path)
        _write_alignment(aligned_era, rms_px=0.1)
        revised, _ = recalibrate_dataset(cube.ds, aligned_era)
        assert float(revised.attrs["wavelength_accuracy_nm"]) == pytest.approx(0.01, abs=1e-6)
        assert revised.attrs["wavelength_accuracy_source"] == "snapshot alignment rms_px"

        # An era that recorded no alignment cannot make the claim at all, and the
        # cube must stop making it on that era's behalf.
        silent_era = _snapshot(
            tmp_path,
            "20250927_cmos",
            shift=WAVELENGTH_SHIFT_NM,
            order_offsets=(0.0, NEW_ERA_ORDER_OFFSET_NM),
            sphere_marker=" 2025b",
        )
        withdrawn, _ = recalibrate_dataset(cube.ds, silent_era)
        assert "wavelength_accuracy_nm" not in withdrawn.attrs
        assert "wavelength_accuracy_source" not in withdrawn.attrs


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
