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

The tooling can validate payload completeness, platform identity, checksums,
offline resolvability, and byte-identical distribution output without changing
the product interfaces. Packet 6 is nevertheless a field-gated release: local
fixtures and cross-platform wheel downloads cannot substitute for execution on
a clean Windows machine and a clean Mac, and an ordinary cached install cannot
substitute for a network-disabled offline rehearsal.

The current source targets a 1.5.0 release candidate while the internal seam is
exercised. Version 1.5.0 is shipped only when all native-platform and offline
gates have durable evidence.
