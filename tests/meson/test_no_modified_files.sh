#!/bin/sh
# T080: Cohabitation invariant. The Meson build must not touch the legacy
# Perl/Makefile path. `git diff` against the merge-base must show ZERO
# changed lines for src/{mfile,cfile,DoConfig,Makefile,Wizard*} and for
# any file under src/ other than src/MakeDesc.cpp (which gains the
# narrow FORCE_BPL flag, FR-018-compatible).
#
# This is the cheapest way to enforce FR-012 in CI; runs in the `lint` job.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

base_ref="${1:-main}"
if ! git rev-parse --verify "$base_ref" >/dev/null 2>&1; then
    echo "FAIL: base ref '$base_ref' does not exist" >&2
    exit 1
fi
merge_base=$(git merge-base HEAD "$base_ref")

# Files that MUST be unchanged by this feature.
protected="src/mfile src/cfile src/DoConfig src/Makefile"

violations=""
for f in $protected; do
    if [ ! -f "$f" ]; then
        # File doesn't exist in HEAD; nothing to check.
        continue
    fi
    if ! git diff --quiet "$merge_base" -- "$f"; then
        violations="$violations $f"
    fi
done

# Also: any Wizard* file in src/ must be unchanged.
for f in $(git ls-tree -r --name-only HEAD src/ | grep -E '^src/Wizard' || true); do
    if ! git diff --quiet "$merge_base" -- "$f"; then
        violations="$violations $f"
    fi
done

if [ -n "$violations" ]; then
    echo "FAIL: FR-012 violation — the following legacy build files were modified:" >&2
    for v in $violations; do
        echo "  $v" >&2
    done
    echo "" >&2
    echo "Diffs:" >&2
    git diff --stat "$merge_base" -- $violations >&2
    exit 1
fi

# Also: under src/, the ONLY file allowed to differ from base is
# src/MakeDesc.cpp (the FORCE_BPL/FORCE_NO_FMA patch).
unexpected_src_changes=$(git diff --name-only "$merge_base" -- src/ \
    | grep -v -E '^src/MakeDesc\.cpp$|^src/meson/|^src/NTL/|^src/config\.h\.in$|^src/meson\.build$' \
    || true)
if [ -n "$unexpected_src_changes" ]; then
    echo "FAIL: unexpected changes to existing src/ files:" >&2
    echo "$unexpected_src_changes" >&2
    echo "" >&2
    echo "Only src/MakeDesc.cpp is allowed to differ from $base_ref." >&2
    exit 1
fi

echo "PASS: T080 no protected file was modified vs $base_ref"
