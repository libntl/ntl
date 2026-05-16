#!/bin/sh
# T029 — Asserts every legacy-build-system file enumerated in FR-001
# is gone from the source tree. RED until T037 (legacy deletion) and
# T038 (cohabitation tooling deletion) land.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

# Note: src/MakeDescAux.cpp is KEPT (defines val_int/val_uint/val_long/
# val_double/val_ldouble used by MakeDesc.cpp during host-side
# mach_desc.h generation; not part of the legacy build path).
LEGACY_PATHS="
    src/configure
    src/DoConfig
    src/Makefile
    src/mfile
    src/cfile
    src/Wizard
    src/WizardAux
    src/TestScript
    src/CopyFeatures
    src/Wizards
    tools/sync-sources.py
    tools/check-sources-in-sync.py
    tools/check-cfile-in-sync.py
    tests/meson/test_symbol_parity_native.sh
"

failed=0
for path in $LEGACY_PATHS; do
    if [ -e "$path" ]; then
        echo "FAIL: legacy artifact still present: $path" >&2
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    echo "FR-001 violation: $failed legacy-path check(s) failed." >&2
    exit 1
fi

echo "PASS: no legacy build-system artifacts found."
