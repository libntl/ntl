#!/usr/bin/env python3
"""Emit NTL/HAVE_<feature>.h headers for every feature referenced by
include/NTL/ALL_FEATURES.h.

NTL's Makefile build runs MakeCheckFeatures, which compiles+executes a
Check<feature>.cpp probe for each feature and writes either an empty
HAVE_<feature>.h (feature absent) or a non-empty one defining
`NTL_HAVE_<feature>` (feature present). The probes require executing
target binaries, which is not safe in cross mode.

For MVP we:
  - Emit a populated header (defining `NTL_HAVE_<FEATURE>`) for features
    we know are present given the spec's C++11 minimum and the standard
    library it implies. COPY_TRAITS1 (std::is_trivially_copyable) and
    CHRONO_TIME (std::chrono) are the load-bearing ones — NTL's
    NTL_SAFE_VECTORS mode is broken without COPY_TRAITS1.
  - Emit an empty stub for every other feature (= absent). NTL's source
    code degrades to portable fallback paths in that case.

A polish-phase follow-up will replace the hardcoded "always present"
list with `cc.compiles()` probes per-feature so each build gets the
optimal set for its target.

Usage:
    gen-have-headers.py <output-directory>
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

# Features assumed present on any C++11-conformant build of NTL.
# COPY_TRAITS1: std::is_trivially_copyable — load-bearing for
# NTL_SAFE_VECTORS' constexpr relocatability traits.
# CHRONO_TIME: std::chrono — used by the build's GetTime5.cpp.
PRESENT_FEATURES = {"COPY_TRAITS1", "CHRONO_TIME"}


def header_body(feature: str, present: bool) -> str:
    if not present:
        return "\n"
    return (
        f"#ifndef NTL_HAVE_{feature}\n"
        f"#define NTL_HAVE_{feature}\n"
        "#endif\n"
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: gen-have-headers.py <outdir>", file=sys.stderr)
        return 2
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    for feat in ALL_FEATURES:
        present = feat in PRESENT_FEATURES
        (out_dir / f"HAVE_{feat}.h").write_text(
            header_body(feat, present), encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
