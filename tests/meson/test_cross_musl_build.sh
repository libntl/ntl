#!/bin/sh
# T037: Cross-build NTL for x86_64-linux-musl from a glibc x86_64 host.
# Asserts the build completes and produces a libntl.so. SKIP when the
# toolchain is absent.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

if ! command -v x86_64-linux-musl-g++ >/dev/null 2>&1; then
    echo "SKIP: x86_64-linux-musl-g++ not installed" >&2
    exit 77
fi

TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
meson setup --cross-file=ci/cross-files/x86_64-linux-musl.txt "$TMP_BUILD" \
    > "$TMP_BUILD-setup.log" 2>&1 \
    || { echo "FAIL: meson setup:" >&2; tail -20 "$TMP_BUILD-setup.log" >&2; exit 1; }

meson compile -C "$TMP_BUILD" > "$TMP_BUILD-compile.log" 2>&1 \
    || { echo "FAIL: meson compile:" >&2; tail -20 "$TMP_BUILD-compile.log" >&2; exit 1; }

libntl=$(find "$TMP_BUILD" -name 'libntl.so*' -type f | head -1)
if [ -z "$libntl" ]; then
    echo "FAIL: libntl.so was not produced" >&2
    exit 1
fi

echo "PASS: T037 x86_64-linux-musl cross-build succeeds"
