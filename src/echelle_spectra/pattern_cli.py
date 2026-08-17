"""CLI for extracting Echelle order-pattern tables from sphere images."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .tools.pattern_extraction import (
    PatternExtractionConfig,
    extract_pattern_from_sphere,
    pattern_row_offsets,
    subtract_background,
)


def _parse_csv_floats(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def _parse_csv_ints(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="echelle-pattern",
        description=(
            "Extract an Echelle order-pattern table from sphere/background images. "
            "By default this previews diagnostics only; pass --output to write a pattern file."
        ),
    )
    parser.add_argument("sphere", help="Sphere/flat-field image file, usually .sif.")
    parser.add_argument("background", help="Matching background image file, usually .sif.")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="PATH",
        help="Write the fitted pattern table to PATH. Omit for preview-only mode.",
    )
    parser.add_argument(
        "--prior-pattern",
        default=None,
        metavar="PATH",
        help="Existing pattern table used for prior-guided tracing.",
    )
    parser.add_argument(
        "--reference-pattern",
        default=None,
        metavar="PATH",
        help="Pattern table to compare against. Defaults to --prior-pattern when supplied.",
    )
    parser.add_argument(
        "--expected-orders",
        type=int,
        default=29,
        help="Expected number of orders in sampled columns (default: 29).",
    )
    parser.add_argument(
        "--thresholds",
        default="0.10,0.11,0.12,0.13,0.14,0.15",
        help="Comma-separated peak thresholds to try.",
    )
    parser.add_argument(
        "--column-starts",
        default="530,680,750",
        help="Comma-separated first sampled columns to try.",
    )
    parser.add_argument(
        "--sample-step",
        type=int,
        default=150,
        help="Column step between sampled slices (default: 150).",
    )
    parser.add_argument(
        "--sample-count",
        type=int,
        default=10,
        help="Number of sampled columns per trial (default: 10).",
    )
    parser.add_argument(
        "--search-radius",
        type=int,
        default=20,
        help="Prior-guided row search radius in pixels (default: 20).",
    )
    parser.add_argument(
        "--no-prior-fit",
        action="store_true",
        help="Do not run the prior-guided fit even when --prior-pattern is supplied.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow --output to replace an existing file.",
    )
    return parser


def _load_image_pair(sphere_path: Path, background_path: Path) -> np.ndarray:
    from .tools.echelle import read_image

    sphere_images, sphere_info = read_image(str(sphere_path), spec="black")
    background_images, background_info = read_image(str(background_path), spec="black")
    image = subtract_background(sphere_images, background_images)

    print(
        "Loaded sphere:",
        sphere_images.shape,
        "frames=",
        sphere_info.get("NumberOfFrames"),
        "exposure=",
        sphere_info.get("ExposureTime"),
    )
    print(
        "Loaded background:",
        background_images.shape,
        "frames=",
        background_info.get("NumberOfFrames"),
        "exposure=",
        background_info.get("ExposureTime"),
    )
    print("Subtracted image:", image.shape)
    return image


def _print_trial_summary(trials, limit: int = 10) -> None:
    print("Top extraction trials:")
    for trial in trials[:limit]:
        status = "OK" if trial.success else "CHECK"
        print(
            f"  {status:5s} threshold={trial.threshold:.2f} "
            f"columns={trial.columns_px[0]}..{trial.columns_px[-1]} "
            f"counts={trial.peak_counts.tolist()}"
        )


def _print_reference_delta(pattern: np.ndarray, reference_path: Path) -> None:
    reference = np.loadtxt(reference_path, dtype=int)
    try:
        delta = pattern_row_offsets(pattern, reference)
    except ValueError:
        print(
            f"Reference shape differs: {reference.shape} vs {pattern.shape}",
            file=sys.stderr,
        )
        return

    print("Delta row px: pattern - reference")
    print("  median:", float(np.median(delta)))
    print("  mean:  ", float(np.mean(delta)))
    print("  std:   ", float(np.std(delta)))
    print("  min/max:", float(np.min(delta)), float(np.max(delta)))


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    sphere_path = Path(args.sphere)
    background_path = Path(args.background)
    if not sphere_path.is_file():
        print(f"ERROR: sphere file not found: {sphere_path}", file=sys.stderr)
        sys.exit(1)
    if not background_path.is_file():
        print(f"ERROR: background file not found: {background_path}", file=sys.stderr)
        sys.exit(1)

    output_path = Path(args.output) if args.output else None
    if output_path and output_path.exists() and not args.overwrite:
        print(f"ERROR: output exists; pass --overwrite: {output_path}", file=sys.stderr)
        sys.exit(1)

    config = PatternExtractionConfig(
        expected_orders=args.expected_orders,
        sample_step_px=args.sample_step,
        sample_count=args.sample_count,
    )
    thresholds = _parse_csv_floats(args.thresholds)
    column_starts = _parse_csv_ints(args.column_starts)

    image = _load_image_pair(sphere_path, background_path)
    prior_path = Path(args.prior_pattern) if args.prior_pattern else None
    prior = None if prior_path is None else np.loadtxt(prior_path, dtype=int)
    # One call, shared with the calibration bench's own in-place extraction, so
    # a pattern written here and a pattern the bench stands on are the same fit.
    try:
        extraction = extract_pattern_from_sphere(
            image,
            prior,
            config=config,
            threshold_values=thresholds,
            column_start_values=column_starts,
            search_radius_px=args.search_radius,
            use_prior=not args.no_prior_fit,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}.", file=sys.stderr)
        sys.exit(1)

    _print_trial_summary(extraction.trials)
    print("Selected trial:" if extraction.trial_succeeded else "No trial succeeded; using:")
    print("  threshold:", extraction.threshold)
    print("  sampled columns:", extraction.columns_px.tolist())
    if extraction.prior_used:
        print("Prior-guided fit:")
        print("  prior:", prior_path)
        print("  search radius:", extraction.search_radius_px)

    result = extraction
    print("Pattern shape:", result.pattern.shape)
    reference = args.reference_pattern or args.prior_pattern
    if reference:
        _print_reference_delta(result.pattern, Path(reference))

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(output_path, result.pattern, fmt="%d")
        print("Saved:", output_path)
    else:
        print("Preview only; pass --output to write a pattern table.")

    sys.exit(0)


if __name__ == "__main__":
    main()
