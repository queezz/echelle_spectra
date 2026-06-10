from __future__ import annotations

from unittest.mock import patch

import pytest

from echelle_spectra.nist_overlay_cli import _build_parser, main


def test_parser_accepts_explicit_line_lists():
    args = _build_parser().parse_args(
        [
            "lamp.sif",
            "--line-list",
            "ThI=th.csv",
            "--orders",
            "8-9",
            "--min-nm",
            "600",
            "--max-nm",
            "640",
            "--output-dir",
            "out",
        ]
    )
    assert args.line_list == ["ThI=th.csv"]
    assert args.lamp == []


def test_main_resolves_lamp_preset_from_cache_dir(tmp_path):
    for name in ("nist_th_i.csv", "nist_th_ii.csv", "nist_ar_i.csv", "nist_ar_ii.csv"):
        (tmp_path / name).write_text("obs_wl_air(nm),intens\n600.0,1\n", encoding="utf-8")

    with patch("echelle_spectra.nist_overlay_cli.run_nist_synthetic_overlay") as run_overlay:
        run_overlay.return_value.output_dir = tmp_path / "out"
        run_overlay.return_value.n_nist_lines = 4
        run_overlay.return_value.n_measured_peaks = 0
        run_overlay.return_value.n_candidate_anchors = 0
        run_overlay.return_value.candidate_table = None

        main(
            [
                "lamp.sif",
                "--lamp",
                "ThAr",
                "--line-list-dir",
                str(tmp_path),
                "--orders",
                "8",
                "--min-nm",
                "600",
                "--max-nm",
                "640",
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )

    config = run_overlay.call_args.args[0]
    assert set(config.nist_csvs) == {"ThI", "ThII", "ArI", "ArII"}


def test_main_requires_line_list_source():
    with pytest.raises(SystemExit) as exc:
        main(["lamp.sif", "--min-nm", "600", "--max-nm", "640", "--output-dir", "out"])
    assert exc.value.code == 2


def test_list_lamps_does_not_require_run_arguments(capsys):
    main(["--list-lamps"])
    assert "thar: ThI, ThII, ArI, ArII" in capsys.readouterr().out
