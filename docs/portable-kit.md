# Portable NIFS kit

The portable kit is the no-admin, no-preinstalled-Python route for campaign
machines at NIFS. Its authoritative user instructions are the repository's
[root README-KIT](https://github.com/queezz/echelle_spectra/blob/master/README-KIT.md).

The repository commits the universal resolver lock, three generated
target-specific hash-locked requirements files, pinned asset manifest,
stdlib-only assembly and
normalization tools, installers, documentation, and tests. uv executables,
CPython archives, wheelhouses, caches, project distributions, extracted
runtimes, virtual environments, and completed kits are generated only in an
explicit external destination.

## Release evidence boundary

The tooling validates payload completeness, platform identity, checksums,
offline resolvability, and byte-identical distribution output without changing
the product interfaces. The current 1.6.0 kit preserves the field-gated release
contract established on clean Windows and native macOS machines, including
network-disabled offline
rehearsals. Future kit releases must repeat those gates; local fixtures and
cross-platform wheel downloads are not substitutes.
