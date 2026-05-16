#!/bin/sh
# T027: mach_desc.h produced by the Meson build on Linux x86_64 must match
# the one produced by the Makefile build (after stripping comments and
# sorting). Validates the FORCE_BPL plumbing on the native path.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

# 1. Meson-generated mach_desc.h
cd "$REPO_ROOT"
meson setup "$TMP_BUILD/meson" >"$TMP_BUILD/setup.log" 2>&1 \
    || { echo "FAIL: meson setup:" >&2; cat "$TMP_BUILD/setup.log" >&2; exit 1; }
meson compile -C "$TMP_BUILD/meson" mach_desc.h >"$TMP_BUILD/compile.log" 2>&1 \
    || { echo "FAIL: building mach_desc.h:" >&2; cat "$TMP_BUILD/compile.log" >&2; exit 1; }
meson_mach=$(find "$TMP_BUILD/meson" -name mach_desc.h | head -1)

# 2. Makefile-generated mach_desc.h. Run MakeDesc the way DoConfig does
# (it writes to `./mach_desc.h` in its cwd, NOT to stdout).
mkdir -p "$TMP_BUILD/make-makedesc"
cp "$REPO_ROOT/src/MakeDesc.cpp" "$REPO_ROOT/src/MakeDescAux.cpp" \
    "$TMP_BUILD/make-makedesc/"
cp -r "$REPO_ROOT/include" "$TMP_BUILD/make-makedesc/"
g++ -O0 -I"$TMP_BUILD/make-makedesc/include" \
    -o "$TMP_BUILD/make-makedesc/MakeDesc" \
    "$TMP_BUILD/make-makedesc/MakeDesc.cpp" \
    "$TMP_BUILD/make-makedesc/MakeDescAux.cpp" -lm
(cd "$TMP_BUILD/make-makedesc" && ./MakeDesc 2>/dev/null)
cp "$TMP_BUILD/make-makedesc/mach_desc.h" "$TMP_BUILD/make-mach.h"

normalize() {
    # Strip C comments and blank lines, then sort.
    sed 's|//.*$||; s|/\*.*\*/||' "$1" | grep -v '^[[:space:]]*$' | sort
}

normalize "$meson_mach"        > "$TMP_BUILD/m1"
normalize "$TMP_BUILD/make-mach.h" > "$TMP_BUILD/m2"

if ! diff -q "$TMP_BUILD/m1" "$TMP_BUILD/m2" >/dev/null; then
    echo "FAIL: mach_desc.h differs between Meson and Makefile paths:" >&2
    diff -u "$TMP_BUILD/m2" "$TMP_BUILD/m1" | head -30 >&2
    exit 1
fi

echo "PASS: T027 mach_desc.h parity"
