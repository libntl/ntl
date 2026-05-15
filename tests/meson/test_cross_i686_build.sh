#!/bin/sh
# T035: Cross-build NTL for i686-linux-gnu from an x86_64 host. Asserts the
# resulting libntl.so is ELF 32-bit Intel 80386.
#
# Requires: i686-linux-gnu-gcc/g++ toolchain (Debian/Ubuntu:
# apt-get install gcc-i686-linux-gnu g++-i686-linux-gnu).
# When the toolchain is absent, the test is treated as SKIP (exit 77).

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

if ! command -v i686-linux-gnu-g++ >/dev/null 2>&1; then
    echo "SKIP: i686-linux-gnu-g++ not installed" >&2
    exit 77
fi

TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
if ! meson setup \
        --cross-file=ci/cross-files/i686-linux-gnu.txt \
        "$TMP_BUILD" \
        > "$TMP_BUILD-setup.log" 2>&1; then
    echo "FAIL: meson setup for i686-linux-gnu:" >&2
    tail -30 "$TMP_BUILD-setup.log" >&2
    exit 1
fi

if ! meson compile -C "$TMP_BUILD" > "$TMP_BUILD-compile.log" 2>&1; then
    echo "FAIL: meson compile for i686-linux-gnu:" >&2
    tail -30 "$TMP_BUILD-compile.log" >&2
    exit 1
fi

libntl=$(find "$TMP_BUILD" -name 'libntl.so*' -type f | head -1)
if [ -z "$libntl" ]; then
    echo "FAIL: libntl.so was not produced" >&2
    exit 1
fi

if ! file "$libntl" | grep -q '80386'; then
    echo "FAIL: libntl.so is not ELF 32-bit 80386:" >&2
    file "$libntl" >&2
    exit 1
fi

echo "PASS: T035 i686-linux-gnu cross-build produces 32-bit ELF"
