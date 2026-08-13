"""Tests for the echelle-validate-lines CLI."""

from __future__ import annotations

import pytest

from echelle_spectra.line_validation_cli import _build_parser, main


def test_parser_defaults_to_accepted_20250926_calibration():
    args = _build_parser().parse_args(["shot.SIF"])
    assert args.pattern == "pattern_CMOS_20250926.txt"
    assert args.wavelength == "alignments/Th_wavelength_CMOS_20240305_aligned_to_20250926.txt"
    assert args.line_set == "balmer-fulcher"
    assert args.wavelength_medium == "air"
    assert args.fulcher_table is None


def test_vacuum_convention_is_rejected_until_explicit_conversion():
    with pytest.raises(SystemExit) as exc:
        main(["shot.SIF", "--wavelength-medium", "vacuum"])
    assert exc.value.code == 1
