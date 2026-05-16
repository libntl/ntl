# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

(Nothing yet for the next release after 12.0.0.)

## [12.0.0] — Unreleased

### Changed (BREAKING)

- **Removed the legacy Perl `./configure` + Makefile + Wizard build path.**
  Meson is now the sole supported build path. Calls to `cd src && ./configure ... && make`
  no longer work — the `configure` script and Makefile are gone. See
  `doc/migration-from-makefile.txt` for a side-by-side option mapping
  and worked examples (manual install, Debian-style packaging, Yggdrasil
  cross-compile). Every option that the legacy `DoConfig` accepted has
  either a Meson equivalent or a documented removal.
- **The auto-tuning Wizard is reincarnated as a Python tool.** Run
  `pip install ./tools && ntl-wizard` (or `ntl-wizard --batch` for
  CI / headless) to produce `src/meson/tune-tables/host-tuned.ini`,
  then build with `meson setup -Dtune=host build`. Full parameter
  parity with legacy `src/Wizard.cpp`. See `doc/wizard.txt`.
- **NTL version bump to 12.0.** SemVer major signals the BREAKING
  build-system change for downstream packagers.

### Removed

- `src/configure`, `src/DoConfig`, `src/mfile`, `src/cfile`, `src/Wizard`,
  `src/WizardAux`, `src/TestScript`, `src/CopyFeatures`.
- `tools/sync-sources.py`, `tools/check-sources-in-sync.py`,
  `tools/check-cfile-in-sync.py` (cohabitation-only tooling).
- `tests/meson/test_symbol_parity_native.sh` (cohabitation-only test).

### Added (feature 002)

- `tools/ntl_wizard/` — Python Wizard package (Typer + Textual TUI +
  `--batch` mode for CI). Native-only execution; refuses cross-compile
  contexts with a pointer to the static tune tables.
- `src/meson/tune-tables/{generic,x86,linux-s390x}.ini` — static tune
  tables in Meson-native INI format, ported from `DoConfig`.
- `src/meson/read-tune-table.py` — Meson-side reader for tune-table
  INIs. Strict on missing keys; forward-compat (warn) on extras.
- `-Dtune=host` and `-Dtune_artifact=PATH` Meson options for consuming
  the Wizard's output.
- `doc/build.txt` (renamed from `doc/build-meson.txt`, rewritten without
  cohabitation framing).
- `doc/migration-from-makefile.txt` — full migration guide.
- `doc/wizard.txt` — Wizard usage reference.
- `tests/ntl_wizard/` — pytest suite (38 tests).
- `tests/meson/test_no_legacy_artifacts.sh` and four other guard tests
  preventing legacy re-introduction and verifying the new tune flow.

### Added (feature 001, retained in 12.0)

- Meson-based build system. Previously cohabited with the legacy build (feature 001);
  is now the sole build path (feature 002).
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
