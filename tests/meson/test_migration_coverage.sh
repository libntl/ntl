#!/bin/sh
# T043 / SC-006 — every captured DoConfig option appears in
# doc/migration-from-makefile.txt (either with a Meson equivalent
# or marked "removed").

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

CAPTURED="specs/002-remove-legacy-build/captured-doconfig-options.txt"
DOC="doc/migration-from-makefile.txt"

if [ ! -f "$CAPTURED" ]; then
    echo "SKIP: $CAPTURED missing (run T044 first)"
    exit 0
fi
if [ ! -f "$DOC" ]; then
    echo "FAIL: $DOC missing (T047)"
    exit 1
fi

missing=0
while IFS=':' read -r option_name _rest; do
    option_name=$(echo "$option_name" | sed 's/[[:space:]]*$//')
    [ -z "$option_name" ] && continue
    if ! grep -q "\b$option_name\b" "$DOC"; then
        echo "FAIL: $option_name not covered in $DOC"
        missing=$((missing + 1))
    fi
done < "$CAPTURED"

if [ "$missing" -ne 0 ]; then
    echo "Migration doc missing $missing option(s) from $CAPTURED"
    exit 1
fi
echo "PASS: every captured DoConfig option is referenced in the migration doc"
