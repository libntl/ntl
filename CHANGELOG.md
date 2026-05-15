# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Meson-based build system that cohabits with the legacy Perl `./configure` + Makefile.
  See `doc/build-meson.txt` for usage.
- `src/MakeDesc.cpp` accepts `-DNTL_FORCE_BPL=N` (N ∈ {32, 64}) and `-DNTL_FORCE_NO_FMA`
  to override host-derived bits-per-long and FMA detection. This enables cross-compile
  workflows that generate `mach_desc.h` on the build host for any target word width.
  Independently useful for native Makefile builds too (no change to the default path).
- GitHub Actions CI workflow `.github/workflows/meson-ci.yml`: native build on Linux,
  Intel macOS, Apple Silicon macOS, and Windows (MinGW-w64 via msys2); cross-compile
  matrix from Linux to musl, ARM, PowerPC, MinGW, FreeBSD, and RISC-V targets.
- Per-target ABI tables under `src/meson/abi-tables/` describing the platform-specific
  properties (right-shift semantics, long-double policy, RPATH style, etc.) for every
  supported triplet. New targets are added by dropping in a single INI file.
- Helper scripts under `tools/`: `sync-sources.py` regenerates the Meson source list
  from `src/mfile`; `check-sources-in-sync.py` enforces no drift in CI.

### Notes

- The auto-tuning Wizard (`TUNE=auto`) is intentionally not supported by the Meson
  build path. Users wanting Wizard-tuned parameters continue to use the legacy
  Makefile build (`./configure TUNE=auto`).
- The Windows build path uses MinGW-w64 only. MSVC support is out of scope for this
  release.
