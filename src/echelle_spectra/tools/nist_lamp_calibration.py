"""NIST ASD helpers for lamp-based wavelength calibration review.

This module deliberately works with cached NIST ASD exports. Network fetching
belongs in notebooks or one-off local scripts so calibration runs stay
reproducible.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "NistSpecies",
    "NistLampPreset",
    "COMMON_NIST_SPECIES",
    "COMMON_LAMP_PRESETS",
    "default_nist_cache_dir",
    "lamp_species",
    "load_nist_asd_exports",
    "normalize_species_key",
    "resolve_cached_line_lists",
]


_ROMAN_TO_STAGE = {
    "I": 1,
    "II": 2,
    "III": 3,
    "IV": 4,
}


@dataclass(frozen=True)
class NistSpecies:
    """One atomic/ionic spectrum in NIST ASD notation."""

    key: str
    element: str
    ion: str
    description: str = ""

    @property
    def nist_name(self) -> str:
        """Return NIST ASD spectrum name, for example ``"Ne II"``."""
        return f"{self.element} {self.ion}"

    @property
    def filename_tokens(self) -> tuple[str, ...]:
        """Return common cache filename tokens for this species."""
        base = f"{self.element.lower()}_{self.ion.lower()}"
        compact = f"{self.element.lower()}{self.ion.lower()}"
        dashed = f"{self.element.lower()}-{self.ion.lower()}"
        return (base, compact, dashed)


@dataclass(frozen=True)
class NistLampPreset:
    """Convenience grouping of likely NIST spectra for a lamp."""

    key: str
    species: tuple[str, ...]
    description: str = ""


def _species(element: str, ion: str, description: str = "") -> NistSpecies:
    return NistSpecies(key=f"{element}{ion}", element=element, ion=ion, description=description)


COMMON_NIST_SPECIES: dict[str, NistSpecies] = {
    item.key: item
    for item in (
        _species("Ar", "I", "neutral argon"),
        _species("Ar", "II", "singly ionized argon"),
        _species("Hg", "I", "neutral mercury"),
        _species("Hg", "II", "singly ionized mercury"),
        _species("He", "I", "neutral helium"),
        _species("He", "II", "singly ionized helium"),
        _species("H", "I", "neutral hydrogen Balmer/Paschen lines"),
        _species("H", "II", "singly ionized hydrogen"),
        _species("Ne", "I", "neutral neon"),
        _species("Ne", "II", "singly ionized neon"),
        _species("Th", "I", "neutral thorium"),
        _species("Th", "II", "singly ionized thorium"),
    )
}


COMMON_LAMP_PRESETS: dict[str, NistLampPreset] = {
    "ar": NistLampPreset("ar", ("ArI", "ArII"), "argon lamp"),
    "h": NistLampPreset("h", ("HI", "HII"), "atomic hydrogen lines"),
    "h2": NistLampPreset("h2", ("HI", "HII"), "H2 lamp atomic lines; molecular H2 is external"),
    "he": NistLampPreset("he", ("HeI", "HeII"), "helium lamp"),
    "hg": NistLampPreset("hg", ("HgI", "HgII"), "mercury lamp"),
    "ne": NistLampPreset("ne", ("NeI", "NeII"), "neon lamp"),
    "th": NistLampPreset("th", ("ThI", "ThII"), "thorium lamp"),
    "thar": NistLampPreset("thar", ("ThI", "ThII", "ArI", "ArII"), "thorium-argon lamp"),
}


def default_nist_cache_dir() -> Path:
    """Return the package-bundled cache directory for curated NIST ASD exports."""
    return Path(__file__).resolve().parents[1] / "resources" / "nist_asd_cache"


def normalize_species_key(value: str) -> str:
    """Normalize labels like ``"Th I"``, ``"th_i"``, or ``"ThI"`` to ``"ThI"``."""
    text = value.strip()
    if not text:
        raise ValueError("species label cannot be empty")

    compact = re.sub(r"[^A-Za-z]", "", text)
    for key in COMMON_NIST_SPECIES:
        if compact.lower() == key.lower():
            return key

    parts = re.findall(r"[A-Za-z]+", text)
    if len(parts) >= 2:
        element = parts[0].capitalize()
        ion = parts[1].upper()
    else:
        match = re.fullmatch(r"([A-Za-z]{1,2})(I{1,3}|IV)", text.replace("_", ""), re.I)
        if not match:
            raise ValueError(f"cannot parse NIST species label: {value!r}")
        element = match.group(1).capitalize()
        ion = match.group(2).upper()

    if ion not in _ROMAN_TO_STAGE:
        raise ValueError(f"unsupported ion stage {ion!r} in {value!r}")
    key = f"{element}{ion}"
    if key not in COMMON_NIST_SPECIES:
        raise ValueError(f"unsupported common NIST species: {key}")
    return key


def lamp_species(lamps: Iterable[str]) -> tuple[str, ...]:
    """Return de-duplicated species keys for one or more common lamp presets."""
    keys: list[str] = []
    for item in lamps:
        lamp_key = item.strip().lower()
        if not lamp_key:
            continue
        if lamp_key not in COMMON_LAMP_PRESETS:
            known = ", ".join(sorted(COMMON_LAMP_PRESETS))
            raise ValueError(f"unknown lamp preset {item!r}; known presets: {known}")
        for species in COMMON_LAMP_PRESETS[lamp_key].species:
            if species not in keys:
                keys.append(species)
    return tuple(keys)


def _parse_nist_value(value) -> float:
    if value is None or pd.isna(value):
        return float("nan")
    text = str(value).strip()
    if not text:
        return float("nan")
    text = text.replace('="', "").replace('"', "").replace("=", "").strip()
    text = re.sub(r"[^0-9Ee+\-.]", "", text)
    if text in {"", ".", "-", "+"}:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def _read_nist_export(path: Path) -> pd.DataFrame:
    """Read CSV or tab-delimited NIST ASD exports."""
    header = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    if "\t" in header and "," not in header:
        return pd.read_csv(path, sep="\t", engine="python").dropna(axis=1, how="all")
    try:
        df = pd.read_csv(path)
    except pd.errors.ParserError:
        return pd.read_csv(path, sep="\t", engine="python").dropna(axis=1, how="all")
    if len(df.columns) == 1 and "\t" in str(df.columns[0]):
        return pd.read_csv(path, sep="\t", engine="python").dropna(axis=1, how="all")
    return df


def load_nist_asd_exports(
    nist_exports: Mapping[str, Path],
    *,
    min_wavelength_nm: float,
    max_wavelength_nm: float,
) -> pd.DataFrame:
    """Load cached NIST ASD exports into normalized wavelength/weight rows."""
    frames = []
    for raw_species, raw_path in nist_exports.items():
        species = normalize_species_key(raw_species)
        path = Path(raw_path)
        df = _read_nist_export(path)
        obs = df.get("obs_wl_air(nm)", pd.Series(index=df.index, dtype=object)).map(
            _parse_nist_value
        )
        ritz = df.get("ritz_wl_air(nm)", pd.Series(index=df.index, dtype=object)).map(
            _parse_nist_value
        )
        intensity = df.get("intens", pd.Series(index=df.index, dtype=object)).map(_parse_nist_value)
        aki = df.get("Aki(s^-1)", pd.Series(index=df.index, dtype=object)).map(_parse_nist_value)
        wavelength = obs.where(np.isfinite(obs), ritz)
        weight = intensity.where(
            np.isfinite(intensity) & (intensity > 0),
            np.log10(aki.clip(lower=1.0)),
        )
        frames.append(
            pd.DataFrame(
                {
                    "species": species,
                    "wavelength_nm": wavelength,
                    "observed_nm": obs,
                    "ritz_nm": ritz,
                    "intensity": intensity,
                    "aki": aki,
                    "weight_raw": weight,
                    "source_path": str(path),
                }
            )
        )
    if not frames:
        return pd.DataFrame(
            columns=[
                "species",
                "wavelength_nm",
                "observed_nm",
                "ritz_nm",
                "intensity",
                "aki",
                "weight_raw",
                "source_path",
                "weight",
            ]
        )

    lines = pd.concat(frames, ignore_index=True)
    lines = lines[np.isfinite(lines["wavelength_nm"])].copy()
    lines = lines[
        (lines["wavelength_nm"] >= min_wavelength_nm)
        & (lines["wavelength_nm"] <= max_wavelength_nm)
    ].copy()
    lines["weight_raw"] = lines["weight_raw"].fillna(1.0).clip(lower=0.1)
    lines["weight"] = lines.groupby("species")["weight_raw"].transform(lambda x: x / x.max())
    return lines.sort_values(["wavelength_nm", "species"]).reset_index(drop=True)


def _cache_match_score(path: Path, species: NistSpecies) -> tuple[int, str]:
    name = path.name.lower()
    score = 0
    for token in species.filename_tokens:
        pattern = rf"(?<![a-z0-9]){re.escape(token)}(?![a-z0-9])"
        if re.search(pattern, name):
            score += 10
    if "nist" in name:
        score += 2
    if "normalized" in name:
        score += 1
    return (-score, name)


def resolve_cached_line_lists(
    *,
    lamps: Sequence[str] = (),
    species: Sequence[str] = (),
    cache_dir: str | Path,
    explicit: Mapping[str, Path] | None = None,
) -> dict[str, Path]:
    """Resolve common lamp/species presets to cached NIST ASD CSV files.

    Explicit mappings override cache-dir discovery. Discovery is filename based
    and expects names containing tokens such as ``nist_ne_i_580_620.csv`` or
    ``nist_thii_300_900.csv``.
    """
    resolved: dict[str, Path] = {}
    if explicit:
        for label, path in explicit.items():
            resolved[normalize_species_key(label)] = Path(path).resolve()

    wanted = list(lamp_species(lamps))
    for item in species:
        key = normalize_species_key(item)
        if key not in wanted:
            wanted.append(key)

    directory = Path(cache_dir)
    if not directory.is_dir():
        raise FileNotFoundError(f"NIST cache directory not found: {directory}")

    files = sorted(path for path in directory.rglob("*.csv") if path.is_file())
    missing: list[str] = []
    for key in wanted:
        if key in resolved:
            continue
        spec = COMMON_NIST_SPECIES[key]
        matches = [path for path in files if _cache_match_score(path, spec)[0] < 0]
        if not matches:
            missing.append(key)
            continue
        resolved[key] = sorted(matches, key=lambda path: _cache_match_score(path, spec))[
            0
        ].resolve()

    if missing:
        available = ", ".join(sorted(path.name for path in files)) or "none"
        raise FileNotFoundError(
            "missing cached NIST ASD exports for "
            + ", ".join(missing)
            + f" in {directory}; available CSV files: {available}"
        )
    return resolved
