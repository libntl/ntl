#!/bin/sh
# T007: src/MakeDesc.cpp must reference NTL_FORCE_NO_FMA, and running it with
# the flag must not advertise FMA. The first check is what makes this a
# meaningful TDD failure on hosts that don't have FMA anyway.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Source-level check: the FORCE_NO_FMA token must appear in MakeDesc.cpp.
if ! grep -q 'NTL_FORCE_NO_FMA' "$REPO_ROOT/src/MakeDesc.cpp"; then
    echo "FAIL: src/MakeDesc.cpp does not reference NTL_FORCE_NO_FMA" >&2
    exit 1
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

cp "$REPO_ROOT/src/MakeDesc.cpp" .
cp "$REPO_ROOT/src/MakeDescAux.cpp" .
cp -r "$REPO_ROOT/include" .

g++ -O0 -I include -DNTL_FORCE_BPL=64 -DNTL_FORCE_NO_FMA \
    -o MakeDesc MakeDesc.cpp MakeDescAux.cpp -lm 2>compile.log

./MakeDesc > mach_desc.h

# When FMA is forced off, NTL_HAVE_FMA must be 0 / absent.
if grep -q '^#define NTL_HAVE_FMA 1$' mach_desc.h ; then
    echo "FAIL: NTL_HAVE_FMA is 1 despite -DNTL_FORCE_NO_FMA" >&2
    exit 1
fi

echo "PASS: T007 NTL_FORCE_NO_FMA honored"
