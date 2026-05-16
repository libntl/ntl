#!/usr/bin/env python3
"""Meson-side reader for tune-table INI files.

Used by `meson.build` when `-Dtune=host` (or `-Dtune=<static>`) is set:
parses the tune-table INI, validates the key set against the frozen
parity list exported by `tools/ntl_wizard.parameters`, emits one
`-DNTL_<KEY>=<VALUE>` flag per line on stdout for Meson to consume.

Invariants:
- forward-compat: unknown extra keys emit a warning to stderr but DO
  NOT cause a non-zero exit. Allows older Meson + newer Wizard.
- backward-compat: missing required keys cause a non-zero exit.
- version-mismatch: an artifact whose [provenance].ntl_version differs
  from the surrounding tree's version.txt is rejected (stale).
- static tune tables (no [provenance] section) are accepted unconditionally.

Exit codes:
  0  success — flags emitted on stdout.
  1  parse error, missing required key, version mismatch, etc.
"""
from __future__ import annotations

import configparser
import sys
from pathlib import Path


# Self-locate so we can import ntl_wizard.parameters from the surrounding tree.
REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS_DIR = REPO_ROOT / "tools"

# Add tools/ to sys.path so we can import ntl_wizard.parameters
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

try:
    from ntl_wizard.parameters import PARAMETERS_BY_NAME
except ImportError:
    # Fallback: hand-rolled canonical list. Kept in sync via T007 parity test.
    PARAMETERS_BY_NAME = {  # type: ignore[assignment]
        name: None for name in (
            "NTL_AVOID_BRANCHING",
            "NTL_CRT_ALTCODE",
            "NTL_CRT_ALTCODE_SMALL",
            "NTL_FFT_BIGTAB",
            "NTL_FFT_LAZYMUL",
            "NTL_GF2X_ALTCODE",
            "NTL_GF2X_ALTCODE1",
            "NTL_GF2X_NOINLINE",
            "NTL_SPMM_ULL",
            "NTL_TBL_REM",
        )
    }


def _err(message: str) -> None:
    print(f"error: {message}", file=sys.stderr)


def _warn(message: str) -> None:
    print(f"warning: {message}", file=sys.stderr)


def read_tune_table(path: Path) -> dict[str, str]:
    """Parse and validate a tune-table INI. Returns the parameter
    dict. Raises ValueError on any spec violation."""
    if not path.exists():
        raise ValueError(f"tune-table file not found: {path}")

    parser = configparser.ConfigParser()
    # ConfigParser lowercases keys by default; we need case-preserved
    # NTL_* names.
    parser.optionxform = str  # type: ignore[assignment]
    try:
        parser.read(path, encoding="utf-8")
    except configparser.Error as exc:
        raise ValueError(f"failed to parse {path}: {exc}") from exc

    if "parameters" not in parser.sections():
        raise ValueError(f"{path}: missing required [parameters] section")

    raw_params = dict(parser["parameters"])
    required = set(PARAMETERS_BY_NAME)
    given = set(raw_params)

    missing = required - given
    if missing:
        raise ValueError(
            f"{path}: tune table missing required keys: {sorted(missing)}. "
            f"This is likely a stale artifact — re-run `ntl-wizard`."
        )

    extras = given - required
    if extras:
        # Forward-compat: warn but don't fail.
        _warn(
            f"{path}: tune table contains unknown extra keys "
            f"(forward-compat, ignored): {sorted(extras)}"
        )

    # Version check (only if [provenance] is present — static tables omit it)
    if "provenance" in parser.sections():
        artifact_version = parser["provenance"].get("ntl_version", "").strip()
        tree_version = (REPO_ROOT / "version.txt").read_text(encoding="utf-8").strip()
        if artifact_version and artifact_version != tree_version:
            raise ValueError(
                f"{path}: stale tune artifact (artifact says ntl_version="
                f"{artifact_version!r}, tree version is {tree_version!r}). "
                f"Re-run `ntl-wizard` against the current tree."
            )

    # Filter to known keys only (extras are dropped after the warning).
    return {k: v for k, v in raw_params.items() if k in required}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        _err("usage: read-tune-table.py <path-to-tune-table.ini>")
        return 1

    path = Path(argv[0])
    try:
        params = read_tune_table(path)
    except ValueError as exc:
        _err(str(exc))
        return 1

    # Emit one flag per line on stdout, sorted for stable Meson output.
    for name in sorted(params):
        print(f"-D{name}={params[name]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
