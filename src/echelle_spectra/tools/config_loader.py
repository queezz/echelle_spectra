try:
    import tomllib
except ImportError:
    import tomli as tomllib


def get_config_from_file(config_path):
    with open(config_path, "rb") as cf:
        return tomllib.load(cf)


def load_config(base_path):
    config_path = base_path / "config.toml"
    defaults_path = base_path / "resources/defaults.toml"
    try:
        config = get_config_from_file(config_path)
    except tomllib.TOMLDecodeError:
        print(f"Warning: invalid TOML in {config_path}; using packaged defaults")
        config = get_config_from_file(defaults_path)
    except OSError:
        config = get_config_from_file(defaults_path)
    config["base_path"] = base_path
    return config


if __name__ == "__main__":
    pass
