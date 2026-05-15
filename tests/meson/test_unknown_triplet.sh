#!/bin/sh
# T025: meson setup -Dabi_triplet=does-not-exist must refuse with the
# FR-013 message naming the missing ABI table file.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
if meson setup -Dabi_triplet=fake-triplet-xyz "$TMP_BUILD" \
        >"$TMP_BUILD-out" 2>&1; then
    echo "FAIL: meson setup succeeded for an unknown triplet" >&2
    exit 1
fi

if ! grep -q 'No ABI table entry' "$TMP_BUILD-out"; then
    echo "FAIL: error did not mention the missing ABI table" >&2
    cat "$TMP_BUILD-out" >&2
    exit 1
fi

echo "PASS: T025 unknown triplet rejected"
