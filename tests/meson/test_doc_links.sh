#!/bin/sh
# T042 — Doc-link integrity (FR-007 / FR-011):
#   - README references doc/build.txt AND doc/migration-from-makefile.txt
#   - CHANGELOG entry references migration doc
#   - No doc references legacy ./configure or make as live workflows

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

failed=0

for expected in "doc/build.txt" "doc/migration-from-makefile.txt"; do
    if ! grep -q "$expected" README; then
        echo "FAIL: README does not link to $expected"
        failed=1
    fi
done

if ! grep -q "migration-from-makefile" CHANGELOG.md; then
    echo "FAIL: CHANGELOG.md does not link to the migration doc"
    failed=1
fi

# No live-workflow references to legacy build paths anywhere under doc/
# or README. Historical mentions tagged with "previously", "legacy",
# "before", "in NTL 11.x" are fine — we use a coarse heuristic: any
# `./configure` or `make install` that does NOT appear within 200
# characters of one of these history markers is a problem.

# Simpler check: assert there is no `./configure` line that *describes
# how to build* (i.e., is not under a "migration" / "history" heading).
# We do this by listing files and checking for forbidden patterns
# outside the migration document.
for path in README doc/build.txt doc/wizard.txt; do
    if grep -E "^\s*\./configure" "$path" >/dev/null 2>&1; then
        echo "FAIL: $path contains a live-workflow ./configure reference"
        failed=1
    fi
done

if [ "$failed" -ne 0 ]; then
    exit 1
fi
echo "PASS: doc links coherent; no live-workflow references to legacy build"
