#!/bin/sh
# T006: MakeDesc must honor -DNTL_FORCE_BPL=32 when run on a 64-bit host
# Test FAILS until T011 patches src/MakeDesc.cpp.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

cd "$TMPDIR"

cp "$REPO_ROOT/src/MakeDesc.cpp" .
cp "$REPO_ROOT/src/MakeDescAux.cpp" .
cp -r "$REPO_ROOT/include" .

g++ -O0 -I include -DNTL_FORCE_BPL=32 \
    -o MakeDesc MakeDesc.cpp MakeDescAux.cpp -lm 2>compile.log

./MakeDesc > mach_desc.h

if ! grep -qE '^#define NTL_BITS_PER_LONG \(32\)' mach_desc.h ; then
    echo "FAIL: NTL_BITS_PER_LONG is not 32 in mach_desc.h" >&2
    grep '^#define NTL_BITS_PER_LONG' mach_desc.h >&2 || true
    exit 1
fi

echo "PASS: T006 NTL_FORCE_BPL=32 honored"
