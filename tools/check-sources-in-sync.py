#!/usr/bin/env python3
"""
T013: Verify src/mfile and src/meson/sources.txt agree on the library source
list. Exit non-zero on drift. Run by CI's `lint` job (FR-015 / SC-007).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SYNC = _REPO_ROOT / "tools" / "sync-sources.py"
_COMMITTED = _REPO_ROOT / "src" / "meson" / "sources.txt"


def main() -> int:
    if not _SYNC.is_file():
        print(f"ERROR: {_SYNC} is missing", file=sys.stderr)
        return 2
    if not _COMMITTED.is_file():
        print(
            f"ERROR: {_COMMITTED} is missing. Run "
            f"`python3 tools/sync-sources.py --write` to generate it.",
            file=sys.stderr,
        )
        return 2

    proc = subprocess.run(
        [sys.executable, str(_SYNC)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        print(
            f"ERROR: sync-sources.py exited {proc.returncode}\n{proc.stderr}",
            file=sys.stderr,
        )
        return 2

    expected = proc.stdout
    actual = _COMMITTED.read_text(encoding="utf-8")
    if expected == actual:
        return 0

    print(
        "ERROR: src/meson/sources.txt is out of sync with src/mfile.\n"
        "Run `python3 tools/sync-sources.py --write` and commit the result.",
        file=sys.stderr,
    )
    # Print a compact diff to help diagnosis.
    import difflib

    diff = difflib.unified_diff(
        actual.splitlines(),
        expected.splitlines(),
        fromfile="src/meson/sources.txt (committed)",
        tofile="src/mfile (current)",
        lineterm="",
    )
    for line in diff:
        print(line, file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
