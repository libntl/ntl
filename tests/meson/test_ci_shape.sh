#!/bin/sh
# T054 — Assert .github/workflows/meson-ci.yml has the right shape for
# the post-feature-002 era: no parity job, no sync-lint jobs, but
# DOES contain the ntl-wizard-tests job.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

WF=".github/workflows/meson-ci.yml"
failed=0

if grep -qi "symbol-parity\|symbol_parity\|test_symbol_parity_native" "$WF"; then
    echo "FAIL: $WF still references symbol-parity infrastructure"
    failed=1
fi

if grep -q "sync-sources\.py\|check-sources-in-sync\|check-cfile-in-sync" "$WF"; then
    echo "FAIL: $WF still invokes legacy cohabitation sync tooling"
    failed=1
fi

if ! grep -q "ntl-wizard-tests" "$WF"; then
    echo "FAIL: $WF does not define the ntl-wizard-tests job"
    failed=1
fi

if [ "$failed" -ne 0 ]; then
    exit 1
fi
echo "PASS: CI workflow has the right post-feature-002 shape"
