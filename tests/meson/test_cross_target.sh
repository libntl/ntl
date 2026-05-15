#!/bin/sh
# Cross-target build test (Phases 5-8). Takes a triplet name; verifies the
# cross-build succeeds and the resulting libntl matches the expected
# architecture. Exits 77 (SKIP) when the cross-toolchain for the given
# triplet is not installed.
#
# Used to satisfy T043-T077 in tasks.md without one shell script per
# triplet — the per-triplet logic is small enough to be data-driven.
#
# Usage:
#   test_cross_target.sh <triplet>
#
# Examples:
#   test_cross_target.sh aarch64-linux-gnu
#   test_cross_target.sh x86_64-w64-mingw32

set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <triplet>" >&2
    exit 2
fi
triplet="$1"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cross_file="$REPO_ROOT/ci/cross-files/$triplet.txt"
if [ ! -f "$cross_file" ]; then
    echo "FAIL: no cross-file at $cross_file" >&2
    exit 1
fi

# Discover the compiler from the cross-file. Crude grep for `cpp = '...'`.
cxx=$(awk -F"'" '/^cpp[[:space:]]*=/ { print $2; exit }' "$cross_file")
if [ -z "$cxx" ]; then
    echo "FAIL: cross-file $cross_file has no `cpp =` line" >&2
    exit 1
fi
if ! command -v "$cxx" >/dev/null 2>&1; then
    echo "SKIP: cross-toolchain compiler '$cxx' not installed" >&2
    exit 77
fi

# Expected `file` substring per architecture.
case "$triplet" in
    i686-linux-gnu|i686-w64-mingw32)            expected_arch='80386' ;;
    x86_64-linux-gnu|x86_64-linux-musl|x86_64-apple-darwin|x86_64-w64-mingw32|x86_64-unknown-freebsd)
                                                expected_arch='x86-64' ;;
    aarch64-linux-gnu|aarch64-linux-musl|aarch64-apple-darwin)
                                                expected_arch='aarch64' ;;
    armv7l-linux-gnueabihf-musl)                expected_arch='ARM' ;;
    powerpc64le-linux-gnu)                      expected_arch='PowerPC' ;;
    riscv64-linux-gnu)                          expected_arch='RISC-V' ;;
    *)                                          expected_arch='' ;;
esac

TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
if ! meson setup --cross-file="$cross_file" "$TMP_BUILD" \
        > "$TMP_BUILD-setup.log" 2>&1; then
    echo "FAIL: meson setup for $triplet:" >&2
    tail -30 "$TMP_BUILD-setup.log" >&2
    exit 1
fi

# Build step is REQUIRED for every triplet (clarification Q4).
if ! meson compile -C "$TMP_BUILD" > "$TMP_BUILD-compile.log" 2>&1; then
    echo "FAIL: meson compile for $triplet:" >&2
    tail -30 "$TMP_BUILD-compile.log" >&2
    exit 1
fi

libntl=$(find "$TMP_BUILD" -name 'libntl*' -type f \
              ! -name '*.p' \
              \( -name 'libntl.so*' -o -name 'libntl-*.dll' -o -name 'libntl*.dylib' \) \
              | head -1)
if [ -z "$libntl" ]; then
    echo "FAIL: no libntl artifact produced for $triplet" >&2
    find "$TMP_BUILD" -name 'libntl*' >&2 || true
    exit 1
fi

if [ -n "$expected_arch" ]; then
    if ! file "$libntl" | grep -q "$expected_arch"; then
        echo "FAIL: $libntl does not match expected '$expected_arch':" >&2
        file "$libntl" >&2
        exit 1
    fi
fi

echo "PASS: $triplet cross-build produced $(basename "$libntl") ($expected_arch)"
