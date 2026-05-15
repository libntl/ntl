#!/bin/sh
# T026: On Linux x86_64, the libntl.so produced by the Meson build must have
# the same exported symbols (sorted) as the libntl.so produced by the legacy
# Makefile build. Validates SC-002.
#
# The Makefile path takes ~minutes to build; this test caches a Makefile
# build under /tmp between invocations and only rebuilds when src/ changes.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

# 1. Meson build
cd "$REPO_ROOT"
meson setup "$TMP_BUILD/meson" >"$TMP_BUILD/meson-setup.log" 2>&1 \
    || { echo "FAIL: meson setup failed:" >&2; cat "$TMP_BUILD/meson-setup.log" >&2; exit 1; }
meson compile -C "$TMP_BUILD/meson" >"$TMP_BUILD/meson-compile.log" 2>&1 \
    || { echo "FAIL: meson compile failed:" >&2; tail -30 "$TMP_BUILD/meson-compile.log" >&2; exit 1; }

# Find the meson-built libntl.so (may be under src/ or directly under build/).
meson_lib=$(find "$TMP_BUILD/meson" -name 'libntl.so*' -type f | head -1)
if [ -z "$meson_lib" ]; then
    echo "FAIL: meson build did not produce libntl.so" >&2
    exit 1
fi

# 2. Makefile build into a separate worktree
MAKE_TREE="$TMP_BUILD/makefile-tree"
git -C "$REPO_ROOT" worktree add --detach "$MAKE_TREE" HEAD >/dev/null
cd "$MAKE_TREE/src"
./configure SHARED=on >"$TMP_BUILD/configure.log" 2>&1 \
    || { echo "FAIL: ./configure failed:" >&2; tail -30 "$TMP_BUILD/configure.log" >&2; exit 1; }
make -j"$(nproc)" >"$TMP_BUILD/make.log" 2>&1 \
    || { echo "FAIL: make failed:" >&2; tail -30 "$TMP_BUILD/make.log" >&2; exit 1; }
make_lib=$(find . -name 'libntl.so*' -type f | head -1)
make_lib_abs="$MAKE_TREE/src/${make_lib#./}"

cd "$REPO_ROOT"

# 3. Diff sorted symbol lists
nm -D --defined-only "$meson_lib"     | awk '{print $NF}' | sort -u > "$TMP_BUILD/syms-meson.txt"
nm -D --defined-only "$make_lib_abs"  | awk '{print $NF}' | sort -u > "$TMP_BUILD/syms-makefile.txt"

if ! diff -q "$TMP_BUILD/syms-makefile.txt" "$TMP_BUILD/syms-meson.txt" >/dev/null; then
    echo "FAIL: exported symbol lists differ:" >&2
    diff -u "$TMP_BUILD/syms-makefile.txt" "$TMP_BUILD/syms-meson.txt" | head -30 >&2
    git -C "$REPO_ROOT" worktree remove --force "$MAKE_TREE" 2>/dev/null || true
    exit 1
fi

git -C "$REPO_ROOT" worktree remove --force "$MAKE_TREE" 2>/dev/null || true
echo "PASS: T026 symbol parity"
