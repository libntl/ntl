#!/bin/sh
# T078: After the Meson work lands, `./configure && make` must continue
# to produce an ABI-compatible libntl.so. This test is SLOW (~5-15 min)
# because it runs a full Makefile build, so it's gated by the env var
# NTL_RUN_SLOW_TESTS=1.
#
# What we check: the symbol surface of the Makefile-built libntl.so
# against the merge-base's Makefile-built libntl.so. They must be
# identical (modulo build-id symbols).

set -eu

if [ "${NTL_RUN_SLOW_TESTS:-0}" != "1" ]; then
    echo "SKIP: T078 needs NTL_RUN_SLOW_TESTS=1 (this test takes 5-15 min)"
    exit 77
fi

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

base_ref="${BASE_REF:-main}"
merge_base=$(git merge-base HEAD "$base_ref")

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"; git worktree remove --force "$TMP/base-tree" 2>/dev/null || true; git worktree remove --force "$TMP/head-tree" 2>/dev/null || true' EXIT

build_makefile() {
    label="$1"
    sha="$2"
    tree="$TMP/$label-tree"
    git worktree add --detach "$tree" "$sha" >/dev/null
    (
        cd "$tree/src"
        ./configure SHARED=on >"$TMP/$label-configure.log" 2>&1
        make -j"$(nproc)" >"$TMP/$label-make.log" 2>&1
    )
    find "$tree/src" -name 'libntl.so*' -type f | head -1
}

base_lib=$(build_makefile base "$merge_base")
head_lib=$(build_makefile head "HEAD")

nm -D --defined-only "$base_lib" | awk '{print $NF}' | sort -u > "$TMP/base-syms"
nm -D --defined-only "$head_lib" | awk '{print $NF}' | sort -u > "$TMP/head-syms"

if ! diff -q "$TMP/base-syms" "$TMP/head-syms" >/dev/null; then
    echo "FAIL: Makefile-built libntl.so symbol surface changed vs $base_ref:" >&2
    diff -u "$TMP/base-syms" "$TMP/head-syms" | head -30 >&2
    exit 1
fi

echo "PASS: T078 Makefile build symbol-compatible with $base_ref"
