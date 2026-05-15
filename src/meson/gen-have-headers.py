#!/usr/bin/env python3
"""Emit NTL/HAVE_<feature>.h headers for every feature referenced by
include/NTL/ALL_FEATURES.h.

NTL's Makefile build runs MakeCheckFeatures, which compiles+executes a
Check<feature>.cpp probe for each feature and writes either an empty
HAVE_<feature>.h (feature absent) or a non-empty one defining
`NTL_HAVE_<feature>` (feature present).

Meson's compile-time probes (cc.compiles(), cc.has_type(), …) cover
most of these features in a cross-compile-safe way; the results are
passed to this script as `--present <feature>` arguments. Features
NOT listed via --present are emitted as empty stubs (= absent),
matching MakeCheckFeatures' fallback behavior. Features assumed
unconditionally present (COPY_TRAITS1 / CHRONO_TIME — required by
NTL_SAFE_VECTORS' constexpr trait machinery on C++11 builds) are
hardcoded.

Usage:
    gen-have-headers.py <output-directory> [--present <feature>]...
"""

from __future__ import annotations

import sys
from pathlib import Path


# Features ALL_FEATURES.h #includes (and therefore must each have a
# HAVE_<name>.h file on the include path).
ALL_FEATURES = [
    "ALIGNED_ARRAY",
    "BUILTIN_CLZL",
    "LL_TYPE",
    "SSSE3",
    "AVX",
    "PCLMUL",
    "AVX2",
    "FMA",
    "AVX512F",
    "COPY_TRAITS1",
    "COPY_TRAITS2",
    "CHRONO_TIME",
    "MACOS_TIME",
    "POSIX_TIME",
    "AES_NI",
    "KMA",
]

# Features assumed unconditionally present on any C++11-conformant build.
# COPY_TRAITS1: std::is_trivially_copyable — load-bearing for
# NTL_SAFE_VECTORS' constexpr relocatability traits.
# CHRONO_TIME: std::chrono — used by the build's GetTime5.cpp.
ALWAYS_PRESENT = {"COPY_TRAITS1", "CHRONO_TIME"}


def header_body(feature: str, present: bool) -> str:
    if not present:
        return "\n"
    return (
        f"#ifndef NTL_HAVE_{feature}\n"
        f"#define NTL_HAVE_{feature}\n"
        "#endif\n"
    )


def main() -> int:
    if len(sys.argv) < 2:
        print(
            "usage: gen-have-headers.py <outdir> [--present <feature>]...",
            file=sys.stderr,
        )
        return 2
    out_dir = Path(sys.argv[1])
    extra_present: set[str] = set()
    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == "--present" and i + 1 < len(sys.argv):
            extra_present.add(sys.argv[i + 1])
            i += 2
        else:
            print(f"unrecognized argument: {sys.argv[i]}", file=sys.stderr)
            return 2

    present_set = ALWAYS_PRESENT | extra_present

    out_dir.mkdir(parents=True, exist_ok=True)
    for feat in ALL_FEATURES:
        present = feat in present_set
        (out_dir / f"HAVE_{feat}.h").write_text(
            header_body(feat, present), encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
