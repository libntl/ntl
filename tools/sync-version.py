#!/usr/bin/env python3
"""
T016: Extract NTL's version from include/NTL/version.h and write it to
version.txt at the repo root. The Meson build's project() reads version.txt.
This decouples the Meson version from upstream NTL's release cadence — when
upstream tags a new version, running this script regenerates version.txt
mechanically.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_VERSION_H = _REPO_ROOT / "include" / "NTL" / "version.h"
_OUT_DEFAULT = _REPO_ROOT / "version.txt"


def extract_version(version_h_path: Path) -> str:
    """Return the NTL version string from `#define NTL_VERSION "...".`"""
    text = version_h_path.read_text(encoding="utf-8")
    match = re.search(r'#define\s+NTL_VERSION\s+"([^"]+)"', text)
    if not match:
        raise RuntimeError(
            f"Could not find `#define NTL_VERSION \"...\"` in {version_h_path}"
        )
    return match.group(1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version-h",
        type=Path,
        default=_VERSION_H,
        help="Path to include/NTL/version.h",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_OUT_DEFAULT,
        help="Path to write version.txt",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit non-zero if version.txt is out of sync",
    )
    args = parser.parse_args()

    version = extract_version(args.version_h)

    if args.check:
        if not args.output.is_file():
            print(f"ERROR: {args.output} is missing", file=sys.stderr)
            return 1
        existing = args.output.read_text(encoding="utf-8").strip()
        if existing != version:
            print(
                f"ERROR: version.txt is {existing!r} but version.h says "
                f"{version!r}",
                file=sys.stderr,
            )
            return 1
        return 0

    args.output.write_text(version + "\n", encoding="utf-8")
    print(f"Wrote {version} to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
