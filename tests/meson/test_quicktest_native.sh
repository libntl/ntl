#!/bin/sh
# T028: The Meson build must run QuickTest, BerlekampTest, and ZZTest
# successfully under `meson test`. Validates FR-007 on the native path.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
meson setup "$TMP_BUILD" >"$TMP_BUILD-setup.log" 2>&1 \
    || { echo "FAIL: meson setup:" >&2; cat "$TMP_BUILD-setup.log" >&2; exit 1; }
meson compile -C "$TMP_BUILD" >"$TMP_BUILD-compile.log" 2>&1 \
    || { echo "FAIL: meson compile:" >&2; tail -30 "$TMP_BUILD-compile.log" >&2; exit 1; }

failed=0
for t in QuickTest BerlekampTest ZZTest; do
    if ! meson test -C "$TMP_BUILD" "$t" >"$TMP_BUILD-test-$t.log" 2>&1; then
        echo "FAIL: $t failed:" >&2
        tail -20 "$TMP_BUILD-test-$t.log" >&2
        failed=1
    fi
done
[ "$failed" -eq 0 ] || exit 1

echo "PASS: T028 QuickTest BerlekampTest ZZTest pass"
