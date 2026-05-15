#!/usr/bin/env python3
"""
T015: Verify the set of placeholders in src/cfile (`@{VAR}` form) matches the
set of placeholders in src/config.h.in (`@VAR@` form). Exit non-zero on drift.

Run by CI's `lint` job. Catches the case where someone edits cfile but forgets
to regenerate config.h.in (or vice versa). Reference R-008.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CFILE = _REPO_ROOT / "src" / "cfile"
_CONFIG_IN = _REPO_ROOT / "src" / "config.h.in"


_CFILE_PAT = re.compile(r"@\{([A-Za-z_][A-Za-z0-9_]*)\}")
_CONFIG_IN_PAT = re.compile(r"@([A-Za-z_][A-Za-z0-9_]*)@")


def placeholders(text: str, pattern: re.Pattern[str]) -> set[str]:
    return set(pattern.findall(text))


def main() -> int:
    if not _CFILE.is_file():
        print(f"ERROR: {_CFILE} is missing", file=sys.stderr)
        return 2
    if not _CONFIG_IN.is_file():
        print(
            f"ERROR: {_CONFIG_IN} is missing. Regenerate it from src/cfile.",
            file=sys.stderr,
        )
        return 2

    cfile_keys = placeholders(_CFILE.read_text(encoding="utf-8"), _CFILE_PAT)
    config_in_keys = placeholders(
        _CONFIG_IN.read_text(encoding="utf-8"), _CONFIG_IN_PAT
    )

    only_in_cfile = sorted(cfile_keys - config_in_keys)
    only_in_config_in = sorted(config_in_keys - cfile_keys)

    if not only_in_cfile and not only_in_config_in:
        return 0

    if only_in_cfile:
        print(
            "ERROR: placeholders present in src/cfile but missing from "
            "src/config.h.in:",
            file=sys.stderr,
        )
        for k in only_in_cfile:
            print(f"  @{{{k}}}", file=sys.stderr)
    if only_in_config_in:
        print(
            "ERROR: placeholders present in src/config.h.in but missing from "
            "src/cfile:",
            file=sys.stderr,
        )
        for k in only_in_config_in:
            print(f"  @{k}@", file=sys.stderr)
    print(
        "\nRegenerate src/config.h.in from src/cfile by replacing every "
        "`@{VAR}` with `@VAR@`.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
