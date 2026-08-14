# Calibration epoch registry and cube provenance

Use the epoch registry when one processing run spans more than one calibration
context. It resolves one validated immutable snapshot for each source before
the source is opened with a calibration.

## One boundary authority

Shot/date boundaries live only in each snapshot's existing `[validity]` table:

```toml
[validity]
shot_from = 193000
shot_to = 199999
date_from = "2025-09-26"
date_to = "2026-08-13"
```

All bounds are inclusive. When an epoch declares both shot and date bounds,
both predicates must match. An omitted lower or upper bound is open-ended.
Dates use canonical `YYYY-MM-DD`; shots are non-negative integers.

`calibration_registry.toml` stores only the deliberate order of snapshot IDs:

```toml
schema = "echelle-calibration-registry/v1"

[[epochs]]
snapshot_id = "20240305_cmos"

[[epochs]]
snapshot_id = "20250926_cmos"
```

The list order is preserved for display and provenance. It is not a
"first match wins" escape hatch. Every referenced snapshot is loaded with full
file-size and SHA-256 verification. The registry is refused if two inclusive
shot/date rectangles can overlap, a snapshot is missing or invalid, an artifact
digest differs, a validity bound is invalid, or an ID is repeated.

A shot-only epoch is unbounded in date and a date-only epoch is unbounded in
shot. They therefore overlap each other and cannot coexist in one registry.
Use one identity strategy consistently, or give every epoch both reviewed
dimensions; the registry never invents shot-over-date or date-over-shot
precedence.

## Source identity and selection

Before calibration loading, the command recognizes:

- an LHD-style leading 5–7 digit shot, such as `193778_Echelle.SIF`;
- an explicit `shot_193778` token;
- a path date written as `YYYYMMDD`, `YYYY-MM-DD`, or `YYYY_MM_DD`.

An eight-digit calendar date is not misread as a shot. More than one distinct
date in a path is an error. Selection evaluates every epoch, requires every
identity kind that the candidate declares, and must find exactly one match.
Missing identity, no match, and ambiguity are different errors; none falls back
to a nearby or current calibration.

Use the registry from a source checkout or portable kit:

```powershell
echelle process D:\NIFS\shots `
  -o D:\NIFS\cubes `
  --registry D:\NIFS\calibration_registry.toml `
  --calibrations D:\NIFS\calibrations `
  --units wmsr `
  --volume-label NIFS-A
```

`--registry` takes calibration authority. It cannot be combined with
`--snapshot-id`, an explicitly configured `--camera`, `--calibration-dir`, or
individual calibration-file overrides. Export settings such as units, crop,
instrument ID, and wavelength medium still follow the established precedence:
built-in defaults, then config TOML, then plan TOML for input routing, then
explicit CLI options.

The same paths can be written in a config:

```toml
[calibration]
registry = "../calibration_registry.toml"
calibrations = "../calibrations"

[export]
units = "wmsr"
```

Registry and calibration paths in a config or plan are resolved relative to
that TOML. Run `echelle status --registry ... --calibrations ...` to validate
the complete registry and print every ordered snapshot and inclusive boundary.
Every cube also carries `snapshot_id`; every receipt record carries the selected
ID even when one batch crosses epochs.

## Aligned recalibration representations

Snapshot-backed exports populate the SpectroCube 0.2 optional contract:

- `detector_pixel` is a wavelength coordinate in raw, zero-origin detector
  columns, with `detector_axis = "column"` and
  `reference_frame = "raw_detector"`. Its values are created before
  `Spectrum` applies a negative-dispersion flip. Their final positions follow
  increasing wavelength, but the numbers still mean raw detector columns.
- `echelle_order` is the non-negative wavelength-aligned order coordinate.
- `wavelength_polynomials_json` uses schema
  `spectrocube.wavelength-polynomials/v1`. Coefficients are the existing
  `numpy.polyfit` results in descending-power order and are never reversed.
  Only orders retained after partial-order removal, finite masking, crop, and
  sort are serialized. Export verifies each polynomial against its stored
  detector pixels to within `5e-10 nm`.
- `applied_absolute_calibration_factor` is present for absolute exports. For
  `wmsr`, its source is `absolute["wmsr"]`, its `source_units` are `counts/s`,
  and its units are `W/m2/nm/sr per (counts/s)`. The application string is
  exactly `stored_intensity = source_signal * applied_absolute_calibration_factor`.

One shared selection path carries wavelength, every intensity frame, detector
pixel, order ID, and factor through wavelength-direction normalization, seam
sort, the recorded finite-value mask, and the inclusive wavelength crop. The
exporter refuses misaligned shapes, non-positive retained factors, duplicate or
scrambled wavelengths, missing represented-order polynomials, or a numerical
factor/polynomial inconsistency. It does not fill or silently repair values.

For the current absolute convention:

```text
source_signal = counts / exposure_time
stored_intensity = source_signal * absolute["wmsr"]
```

Dividing a stored absolute cube by its factor therefore recovers the retained
count rate numerically.

## Complete snapshot provenance

Snapshot-backed cubes store portable provenance rather than the machine path to
the calibration folder:

| Attribute | Meaning |
| --- | --- |
| `snapshot_id` | Selected immutable calibration identity |
| `snapshot_manifest_sha256` | Digest of the exact `snapshot.toml` binder |
| `snapshot_manifest_json` | Complete binder, including base, validity, alignment, QC, lamps, notes, and artifact records |
| `calibration_file_digests_json` | Every artifact role/path/source name/size/SHA-256, including pattern, wavelength, sphere, sphere background, integral, and lamp inputs |
| `calibration_registry_schema` | Registry contract used by the writer |
| `calibration_registry_sha256` | Digest of the exact ordered registry |
| `calibration_registry_epoch_position` | One-based position of the selected entry |

Existing extraction provenance remains present: detector width, extraction
half-width, order-border raw-pixel ranges, order wavelength ranges, source,
timing, crop, and dropped-column counts.

## Limitations and recovery

- Legacy exports that do not request `--registry` retain the manual/configured
  calibration path and do not claim complete snapshot provenance.
- Counts cubes can carry detector pixels, orders, polynomials, and snapshot
  provenance, but do not carry an *applied absolute* factor because no absolute
  factor produced their stored intensity.
- Wavelength and absolute scaling can be revised from a provenance-complete
  absolute cube. Extraction geometry cannot: a changed pattern requires the raw
  SIF and a new extraction.
- On no-match or missing-identity errors, correct the source naming/evidence or
  create a new snapshot with reviewed validity. Do not widen a boundary merely
  to make processing continue.
- On a digest failure, restore the immutable artifact byte-for-byte or create a
  new snapshot ID. Do not edit the binder digest to bless changed calibration
  input.
