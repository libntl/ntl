#!/bin/sh
# T009: check-cfile-in-sync.py must detect when cfile and config.h.in have
# diverged in placeholder set. Test FAILS until T015 implements the check.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

if [ ! -f "$REPO_ROOT/tools/check-cfile-in-sync.py" ]; then
    echo "FAIL: tools/check-cfile-in-sync.py does not exist" >&2
    exit 1
fi

FAKE_REPO="$TMPDIR/fake"
mkdir -p "$FAKE_REPO/src" "$FAKE_REPO/tools"
cp "$REPO_ROOT/src/cfile" "$FAKE_REPO/src/cfile"
cp "$REPO_ROOT/tools/check-cfile-in-sync.py" "$FAKE_REPO/tools/"

# Build a config.h.in missing one of cfile's @{VAR} placeholders, by stripping
# the entire line that contains the first @{...} occurrence.
awk '
{ gsub(/@\{([A-Za-z_][A-Za-z0-9_]*)\}/, "@\\1@"); print }
' "$REPO_ROOT/src/cfile" | sed '/@/{1d}' > "$FAKE_REPO/src/config.h.in"

if (cd "$FAKE_REPO" && python3 tools/check-cfile-in-sync.py >/dev/null 2>&1); then
    echo "FAIL: check-cfile-in-sync.py did not detect drift" >&2
    exit 1
fi

echo "PASS: T009 cfile drift detected"
