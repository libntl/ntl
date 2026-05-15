#!/usr/bin/env python3
"""
T012: Parse src/mfile and extract the NTL library source list (the `SRC`
variable). Print one source file per line. With --write, persist to
src/meson/sources.txt.

The legacy build's source list is the source of truth (FR-012 forbids
modifying mfile). This script keeps the Meson build's source list mechanically
synchronized with it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parent.parent


def extract_variable(mfile_path: Path, variable_name: str) -> list[str]:
    r"""Extract the tokens assigned to `variable_name` in a Makefile-style file.

    Recognizes the Make syntax `NAME=val val val \` with backslash line
    continuations. Stops at the first non-continued line.
    """

    text = mfile_path.read_text(encoding="utf-8")
    prefix = f"{variable_name}="
    lines = text.splitlines()
    in_var = False
    collected: list[str] = []
    for line in lines:
        if not in_var:
            stripped = line.lstrip()
            if stripped.startswith(prefix):
                in_var = True
                body = stripped[len(prefix):]
                continued = body.rstrip().endswith("\\")
                if continued:
                    body = body.rstrip()[:-1]
                collected.append(body)
                if not continued:
                    break
            continue
        # Already inside the multi-line value.
        continued = line.rstrip().endswith("\\")
        body = line.rstrip()[:-1] if continued else line
        collected.append(body)
        if not continued:
            break

    if not collected:
        raise RuntimeError(
            f"Could not locate {variable_name}= in {mfile_path}"
        )

    return " ".join(collected).split()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mfile",
        type=Path,
        default=_REPO_ROOT / "src" / "mfile",
        help="Path to src/mfile (default: repo-relative)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_REPO_ROOT / "src" / "meson" / "sources.txt",
        help="Path to write sources.txt when --write is given",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write the list to --output instead of stdout",
    )
    args = parser.parse_args()

    tokens = extract_variable(args.mfile, "SRC")
    # Defensive filter: only .cpp filenames.
    sources = sorted(t for t in tokens if t.endswith(".cpp"))
    if not sources:
        raise RuntimeError("SRC variable did not contain any .cpp files")
    body = "\n".join(sources) + "\n"

    if args.write:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(body, encoding="utf-8")
        print(f"Wrote {len(sources)} sources to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(body)
    return 0


if __name__ == "__main__":
    sys.exit(main())
