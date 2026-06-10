"""CLI for synthetic line-list overlays on lamp spectra."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import _config
from .tools.nist_lamp_calibration import (
    COMMON_LAMP_PRESETS,
    COMMON_NIST_SPECIES,
    default_nist_cache_dir,
    resolve_cached_line_lists,
)
from .tools.nist_synthetic_overlay import NistOverlayConfig, run_nist_synthetic_overlay


def _default_calibration_dir() -> Path:
    return _config["base_path"] / "resources" / "calibration_files"


def _parse_species_csvs(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for item in values:
        species, sep, path = item.partition("=")
        if not sep or not species.strip() or not path.strip():
            raise argparse.ArgumentTypeError(
                "Line-list arguments must look like Species=/path/to/file.csv"
            )
        result[species.strip()] = Path(path).resolve()
    return result


def _parse_csv_values(values: list[str] | None) -> tuple[str, ...]:
    if not values:
        return ()
    items: list[str] = []
    for value in values:
        items.extend(part.strip() for part in value.split(",") if part.strip())
    return tuple(dict.fromkeys(items))


def _parse_orders(value: str) -> tuple[int, ...]:
    orders: list[int] = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            orders.extend(range(int(start), int(end) + 1))
        else:
            orders.append(int(part))
    return tuple(dict.fromkeys(orders))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-nist-overlay",
        description="Overlay lamp order spectra with cached line-list CSV sticks.",
    )
    parser.add_argument("signal", nargs="?", help="Lamp signal image, usually a .sif file.")
    parser.add_argument(
        "--background",
        default=None,
        help="Optional matching lamp background image.",
    )
    parser.add_argument(
        "--calibration-dir",
        default=str(_default_calibration_dir()),
        help="Directory containing wavelength, pattern, sphere, and integrating-sphere files.",
    )
    parser.add_argument("--wavelength", default="Th_wavelength_CMOS_20240305.txt")
    parser.add_argument("--pattern", default="pattern_CMOS_20240305.txt")
    parser.add_argument("--sphere", default="sphere_cmos_20240305.sif")
    parser.add_argument("--sphere-background", default="sphere_cmos_20240305_bkg.sif")
    parser.add_argument("--integral", default="integrating_sphere.txt")
    parser.add_argument(
        "--line-list",
        action="append",
        default=[],
        help="Species=/path/to/cached.csv. Repeat for multiple species, e.g. ThI=... ArI=...",
    )
    parser.add_argument(
        "--lamp",
        action="append",
        default=[],
        help=(
            "Common lamp preset to resolve from --line-list-dir. May be repeated or comma-separated. "
            "Known presets include: " + ", ".join(sorted(COMMON_LAMP_PRESETS))
        ),
    )
    parser.add_argument(
        "--species",
        action="append",
        default=[],
        help=(
            "NIST species to resolve from --line-list-dir, e.g. HgI,NeII,ThII. "
            "Known species include: " + ", ".join(sorted(COMMON_NIST_SPECIES))
        ),
    )
    parser.add_argument(
        "--line-list-dir",
        default=str(default_nist_cache_dir()),
        help="Directory of cached NIST ASD CSV exports used with --lamp or --species.",
    )
    parser.add_argument(
        "--list-lamps",
        action="store_true",
        help="Print common lamp presets and exit.",
    )
    parser.add_argument("--orders", default="8", help="Comma/range order list, e.g. 6-10 or 8.")
    parser.add_argument("--min-nm", type=float)
    parser.add_argument("--max-nm", type=float)
    parser.add_argument("--output-dir")
    parser.add_argument("--candidate-table-out", default=None)
    parser.add_argument("--match-tolerance-nm", type=float, default=0.05)
    parser.add_argument("--min-line-weight", type=float, default=0.18)
    parser.add_argument("--dominance-threshold", type=float, default=0.45)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.list_lamps:
        print("Common NIST lamp presets:")
        for key, preset in sorted(COMMON_LAMP_PRESETS.items()):
            species = ", ".join(preset.species)
            print(f"  {key}: {species} ({preset.description})")
        return

    if args.signal is None:
        parser.error("signal is required unless --list-lamps is used")
    if args.min_nm is None or args.max_nm is None:
        parser.error("--min-nm and --max-nm are required unless --list-lamps is used")
    if args.output_dir is None:
        parser.error("--output-dir is required unless --list-lamps is used")

    calibration_dir = Path(args.calibration_dir).resolve()
    explicit_line_lists = _parse_species_csvs(args.line_list)
    lamps = _parse_csv_values(args.lamp)
    species = _parse_csv_values(args.species)
    if lamps or species:
        nist_csvs = resolve_cached_line_lists(
            lamps=lamps,
            species=species,
            cache_dir=args.line_list_dir,
            explicit=explicit_line_lists,
        )
    else:
        nist_csvs = explicit_line_lists
    if not nist_csvs:
        parser.error(
            "provide at least one --line-list or use --lamp/--species with --line-list-dir"
        )

    config = NistOverlayConfig(
        calibration_dir=calibration_dir,
        signal_file=Path(args.signal).resolve(),
        background_file=Path(args.background).resolve() if args.background else None,
        wavelength_file=args.wavelength,
        pattern_file=args.pattern,
        sphere_file=args.sphere,
        sphere_background_file=args.sphere_background,
        integral_file=args.integral,
        nist_csvs=nist_csvs,
        output_dir=Path(args.output_dir).resolve(),
        orders=_parse_orders(args.orders),
        min_wavelength_nm=args.min_nm,
        max_wavelength_nm=args.max_nm,
        candidate_table_out=(
            Path(args.candidate_table_out).resolve() if args.candidate_table_out else None
        ),
        match_tolerance_nm=args.match_tolerance_nm,
        min_nist_weight=args.min_line_weight,
        dominance_threshold=args.dominance_threshold,
    )
    result = run_nist_synthetic_overlay(config)
    print(f"Wrote {result.output_dir}")
    print(f"NIST/list lines: {result.n_nist_lines}")
    print(f"Measured peaks: {result.n_measured_peaks}")
    print(f"Candidate anchors: {result.n_candidate_anchors}")
    if result.candidate_table is not None:
        print(f"Candidate table: {result.candidate_table}")


if __name__ == "__main__":
    main()
