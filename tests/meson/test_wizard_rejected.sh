#!/bin/sh
# T024: meson setup -Dtune=auto must be rejected with a clear diagnostic.
# Should pass immediately because meson.options declares tune as a `combo`
# with allowed values {generic, x86, linux-s390x}; 'auto' is rejected at
# the option-parse level (FR-014).

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)/build"
trap 'rm -rf "$(dirname "$TMP_BUILD")"' EXIT

cd "$REPO_ROOT"
if meson setup -Dtune=auto "$TMP_BUILD" >"$TMP_BUILD-out" 2>&1; then
    echo "FAIL: meson setup -Dtune=auto unexpectedly succeeded" >&2
    exit 1
fi

# Either the option-parse error names 'tune' / 'auto', or the explicit
# error() branch in meson.build does. Both are acceptable.
if ! grep -qE 'tune|auto|Wizard' "$TMP_BUILD-out"; then
    echo "FAIL: rejection message did not mention 'tune', 'auto', or 'Wizard':" >&2
    cat "$TMP_BUILD-out" >&2
    exit 1
fi

echo "PASS: T024 -Dtune=auto rejected"
