#!/usr/bin/env python3
"""Run MakeDesc in a temp directory and emit the resulting mach_desc.h to stdout.

MakeDesc.cpp writes its output to a literal file named "mach_desc.h" in its
current working directory; it does NOT write to stdout. This wrapper runs the
binary in a sandbox tempdir, reads the produced file, and writes it to stdout
so a Meson `custom_target(capture: true)` can route it to the right place.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: run-makedesc.py <path-to-MakeDesc>", file=sys.stderr)
        return 2
    makedesc = Path(sys.argv[1]).resolve()
    if not makedesc.is_file():
        print(f"ERROR: MakeDesc binary not found at {makedesc}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory(prefix="meson-makedesc-") as td:
        proc = subprocess.run(
            [str(makedesc)],
            cwd=td,
            capture_output=True,
            text=True,
        )
        # MakeDesc prints diagnostics to stderr; pass them through for the build log.
        sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            sys.stderr.write(
                f"ERROR: MakeDesc exited {proc.returncode}\n"
            )
            return proc.returncode
        produced = Path(td) / "mach_desc.h"
        if not produced.is_file():
            sys.stderr.write(
                "ERROR: MakeDesc did not produce mach_desc.h\n"
            )
            return 1
        sys.stdout.write(produced.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
