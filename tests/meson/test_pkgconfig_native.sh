#!/bin/sh
# T029: After `meson install`, a downstream user must be able to compile and
# link a small program using `pkg-config --cflags --libs ntl`. Validates FR-006.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

cd "$REPO_ROOT"
meson setup "$TMP_BUILD/build" --prefix=/usr/local >/dev/null 2>&1 \
    || { echo "FAIL: meson setup" >&2; exit 1; }
meson compile -C "$TMP_BUILD/build" >/dev/null 2>&1 \
    || { echo "FAIL: meson compile" >&2; exit 1; }

DESTDIR="$TMP_BUILD/install"
DESTDIR="$DESTDIR" meson install -C "$TMP_BUILD/build" --quiet \
    || { echo "FAIL: meson install" >&2; exit 1; }

# Locate the installed ntl.pc.
NTL_PC=$(find "$DESTDIR" -name ntl.pc | head -1)
if [ -z "$NTL_PC" ]; then
    echo "FAIL: ntl.pc was not installed" >&2
    exit 1
fi

# Use pkg-config against the installed tree.
export PKG_CONFIG_PATH="$(dirname "$NTL_PC")"
# Some installs use ${prefix} expansions assuming the configured prefix; we
# wrote /usr/local but installed to DESTDIR. Override prefix on the command
# line to make pkg-config point at the DESTDIR-relative locations.
prefix_override="$DESTDIR/usr/local"
CFLAGS=$(pkg-config --define-variable=prefix="$prefix_override" --cflags ntl)
LIBS=$(pkg-config --define-variable=prefix="$prefix_override" --libs ntl)

cat >"$TMP_BUILD/main.cpp" <<'EOF'
#include <NTL/ZZ.h>
#include <iostream>
int main() {
    NTL::ZZ a = NTL::conv<NTL::ZZ>(2);
    NTL::ZZ b = NTL::conv<NTL::ZZ>(3);
    std::cout << (a + b) << std::endl;
    return 0;
}
EOF

g++ $CFLAGS -o "$TMP_BUILD/main" "$TMP_BUILD/main.cpp" $LIBS \
    || { echo "FAIL: sample program failed to compile/link" >&2; exit 1; }

# libdir can be lib/, lib64/, lib/x86_64-linux-gnu/, etc. — discover from
# pkg-config itself rather than guessing.
LIB_PATH=$(pkg-config --define-variable=prefix="$prefix_override" \
               --variable=libdir ntl)
LD_LIBRARY_PATH="$LIB_PATH:${LD_LIBRARY_PATH:-}" \
    "$TMP_BUILD/main" > "$TMP_BUILD/out"
if [ "$(cat "$TMP_BUILD/out")" != "5" ]; then
    echo "FAIL: sample program output was not '5': $(cat "$TMP_BUILD/out")" >&2
    exit 1
fi

echo "PASS: T029 pkg-config flow"
