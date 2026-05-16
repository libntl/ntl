#!/bin/sh
# T023: meson setup must succeed on the native host once the x86_64-linux-gnu
# ABI table exists. Fails before Phase 3 T030.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
if ! meson setup "$TMP_BUILD" 2>"$TMP_BUILD-stderr"; then
    echo "FAIL: meson setup failed:" >&2
    cat "$TMP_BUILD-stderr" >&2 || true
    exit 1
fi

if [ ! -f "$TMP_BUILD/build.ninja" ]; then
    echo "FAIL: build.ninja was not generated" >&2
    exit 1
fi

echo "PASS: T023 meson setup smoke"
