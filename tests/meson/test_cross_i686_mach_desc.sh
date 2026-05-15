#!/bin/sh
# T036: Cross-build for i686-linux-gnu and assert the generated mach_desc.h
# reports NTL_BITS_PER_LONG (32). Validates the FORCE_BPL plumbing on the
# cross path. SKIP when the toolchain is absent.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

if ! command -v i686-linux-gnu-g++ >/dev/null 2>&1; then
    echo "SKIP: i686-linux-gnu-g++ not installed" >&2
    exit 77
fi

TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
meson setup --cross-file=ci/cross-files/i686-linux-gnu.txt "$TMP_BUILD" >/dev/null 2>&1
meson compile -C "$TMP_BUILD" src/NTL/mach_desc.h >/dev/null 2>&1

mach=$(find "$TMP_BUILD" -name mach_desc.h | head -1)
if ! grep -qE '^#define NTL_BITS_PER_LONG \(32\)' "$mach"; then
    echo "FAIL: NTL_BITS_PER_LONG is not 32 in cross-built mach_desc.h" >&2
    grep '^#define NTL_BITS_PER_LONG' "$mach" >&2 || true
    exit 1
fi

echo "PASS: T036 cross-built mach_desc.h has NTL_BITS_PER_LONG=32"
