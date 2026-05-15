#!/bin/sh
# T084: Verify CHANGELOG.md follows Keep a Changelog format. Runs in the
# `lint` CI job.
#
# Minimum requirements:
#  - File exists at repo root
#  - Contains the header line marking it as a changelog
#  - Has a [Unreleased] section
#  - Uses one or more of the standard category headings (Added, Changed,
#    Deprecated, Removed, Fixed, Security) under at least one version

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cl="$REPO_ROOT/CHANGELOG.md"

if [ ! -f "$cl" ]; then
    echo "FAIL: CHANGELOG.md is missing at repo root" >&2
    exit 1
fi

ok=1

if ! grep -q '^# Changelog' "$cl"; then
    echo "FAIL: CHANGELOG.md missing '# Changelog' header" >&2
    ok=0
fi

if ! grep -qE '^## \[Unreleased\]' "$cl"; then
    echo "FAIL: CHANGELOG.md missing '## [Unreleased]' section" >&2
    ok=0
fi

if ! grep -qE '^### (Added|Changed|Deprecated|Removed|Fixed|Security)' "$cl"; then
    echo "FAIL: CHANGELOG.md has no entries under a recognized category" >&2
    echo "  Allowed: Added / Changed / Deprecated / Removed / Fixed / Security" >&2
    ok=0
fi

if ! grep -q 'keepachangelog.com' "$cl" && ! grep -q 'Keep a Changelog' "$cl"; then
    echo "FAIL: CHANGELOG.md missing reference to Keep a Changelog spec" >&2
    ok=0
fi

if [ "$ok" -eq 1 ]; then
    echo "PASS: T084 CHANGELOG.md is in Keep a Changelog format"
    exit 0
fi
exit 1
