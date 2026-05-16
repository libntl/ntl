#!/bin/sh
# T031 — `-Dtune=<static>` must select the matching INI from
# src/meson/tune-tables/. RED until T033 ports the legacy tables.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

cd "$REPO_ROOT"

failed=0
for tune in generic x86 linux-s390x; do
    ini="src/meson/tune-tables/${tune}.ini"
    if [ ! -e "$ini" ]; then
        echo "FAIL: static tune table missing at $ini"
        failed=1
        continue
    fi
    if ! meson setup -Dtune="$tune" "$TMP_BUILD/$tune" >"$TMP_BUILD/$tune.log" 2>&1; then
        echo "FAIL: meson setup -Dtune=$tune failed"
        tail -20 "$TMP_BUILD/$tune.log"
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    exit 1
fi
echo "PASS: all static tune tables (generic, x86, linux-s390x) configure cleanly"
