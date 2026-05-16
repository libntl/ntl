#!/bin/sh
# T008: check-sources-in-sync.py must detect drift between mfile and sources.txt.
# Test FAILS until T013 implements the check.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

if [ ! -x "$REPO_ROOT/tools/check-sources-in-sync.py" ] \
   && [ ! -f "$REPO_ROOT/tools/check-sources-in-sync.py" ]; then
    echo "FAIL: tools/check-sources-in-sync.py does not exist" >&2
    exit 1
fi

# Set up a fake repo where sources.txt is missing one entry that mfile has.
FAKE_REPO="$TMPDIR/fake"
mkdir -p "$FAKE_REPO/src/meson" "$FAKE_REPO/tools"

cp "$REPO_ROOT/src/mfile" "$FAKE_REPO/src/mfile"
cp "$REPO_ROOT/tools/sync-sources.py" "$FAKE_REPO/tools/" 2>/dev/null || {
    echo "FAIL: tools/sync-sources.py does not exist" >&2
    exit 1
}
cp "$REPO_ROOT/tools/check-sources-in-sync.py" "$FAKE_REPO/tools/"

# Create a sources.txt missing the last line so the check should fail.
python3 "$FAKE_REPO/tools/sync-sources.py" > "$FAKE_REPO/src/meson/sources.txt"
sed -i '$d' "$FAKE_REPO/src/meson/sources.txt"

if (cd "$FAKE_REPO" && python3 tools/check-sources-in-sync.py >/dev/null 2>&1); then
    echo "FAIL: check-sources-in-sync.py did not detect drift" >&2
    exit 1
fi

echo "PASS: T008 sync-sources drift detected"
