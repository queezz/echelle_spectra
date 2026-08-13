from pathlib import Path

from echelle_spectra.tools.config_loader import load_config


def _base(tmp_path: Path) -> Path:
    resources = tmp_path / "resources"
    resources.mkdir()
    (resources / "defaults.toml").write_text('diag_name = "default"\n', encoding="utf-8")
    return tmp_path


def test_missing_config_uses_defaults_without_writing_into_package(tmp_path: Path) -> None:
    base = _base(tmp_path)
    config = load_config(base)
    assert config["diag_name"] == "default"
    assert config["base_path"] == base
    assert not (base / "config.toml").exists()


def test_invalid_config_is_preserved_while_defaults_are_used(tmp_path: Path, capsys) -> None:
    base = _base(tmp_path)
    config_path = base / "config.toml"
    config_path.write_text("invalid = [", encoding="utf-8")
    config = load_config(base)
    assert config["diag_name"] == "default"
    assert config_path.read_text(encoding="utf-8") == "invalid = ["
    assert "using packaged defaults" in capsys.readouterr().out
