# `tools/` — Helper Scripts for the Meson Build

These scripts support the Meson build and its CI. They are **not** invoked by
end users; they exist to keep the Meson and Makefile build descriptions in sync
and to provide CI guard-rails.

## Scripts

| Script | Purpose | Run by |
|---|---|---|
| `sync-sources.py` | Parse `src/mfile`'s `SRC` list and emit one source path per line. With `--write`, writes to `src/meson/sources.txt`. | Maintainer (after `mfile` changes); never at build time. |
| `check-sources-in-sync.py` | Run `sync-sources.py` into a temporary file and diff against the committed `src/meson/sources.txt`. Exit non-zero on drift. | CI `lint` job. |
| `check-cfile-in-sync.py` | Verify that the `@VAR@` placeholders in `src/config.h.in` form the same set as the `@{VAR}` placeholders in `src/cfile`. Exit non-zero on drift. | CI `lint` job. |
| `sync-version.py` | Read NTL's version from upstream sources (the `version.h` `NTL_VERSION` macro) and write `version.txt` at repo root. | Maintainer (when bumping). |
| `check-commit-trailer.sh` | Verify every commit on the current branch ends with the `AI-Assisted: Claude (Spec-Driven Development, TDD methodology)` trailer. Exit non-zero if any commit is missing it. | CI `lint` job. |

## Why these scripts exist

NTL has two parallel build systems by design (the legacy Makefile and the new
Meson build cohabit per `cross-compile-roadmap.md`). Each system has its own
source list and configuration-header template. The scripts here make sure those
copies don't drift apart: when `mfile` gains a source file, `sync-sources.py`
regenerates `sources.txt` mechanically and the CI guard catches anyone who
forgets to do so.

## Conventions

- All scripts are Python 3.8+ or POSIX shell; no Perl (Perl already lives in
  `src/DoConfig`).
- Scripts are idempotent: running them twice in a row produces no diff on the
  second run.
- Scripts that mutate the source tree only do so when passed `--write`; without
  it they print what they would do (dry-run / drift-check mode).
