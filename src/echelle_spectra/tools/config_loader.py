import shutil
import tomli


def get_config_from_file(config_path):
    with open(config_path, "rb") as cf:
        return tomli.load(cf)


def load_config(base_path):
    try:
        config = get_config_from_file(base_path / "config.toml")
    except (tomli.TOMLDecodeError, OSError):
        print("TOMLDecodError in config_loader.py")
        print("Warning: restoring default config file")
        shutil.copy(base_path / "resources/defaults.toml", base_path / "config.toml")
        config = get_config_from_file(base_path / "config.toml")
    config["base_path"] = base_path
    return config


if __name__ == "__main__":
    pass
