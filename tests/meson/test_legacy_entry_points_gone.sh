#!/bin/sh
# T040 / FR-013: assert legacy entry points (./configure and make
# from the source root) are no longer functional. Guards against
# accidental re-introduction.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

# 1) ./configure must not exist anywhere a legacy user would expect it.
for candidate in ./configure src/configure; do
    if [ -x "$candidate" ] || [ -f "$candidate" ]; then
        echo "FAIL: $candidate exists; legacy build entry point should be removed" >&2
        exit 1
    fi
done

# 2) `make` from the source root or src/ should NOT succeed.
#    No top-level Makefile, no src/Makefile, and no GNUmakefile.
for d in . src; do
    for mf in Makefile makefile GNUmakefile; do
        if [ -f "$d/$mf" ]; then
            echo "FAIL: $d/$mf exists; legacy Makefile-based build should be removed" >&2
            exit 1
        fi
    done
done

# 3) Spot-check that the Meson entry point IS present.
if [ ! -f meson.build ]; then
    echo "FAIL: meson.build not at repository root" >&2
    exit 1
fi

echo "PASS: legacy entry points (./configure, make) are gone; meson.build is the canonical entry."
