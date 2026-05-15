#!/bin/sh
# T010: pick-abi.py must refuse to load an ABI table missing a required key,
# naming the missing key in the error. Test FAILS until T018 implements the check.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

if [ ! -f "$REPO_ROOT/src/meson/pick-abi.py" ]; then
    echo "FAIL: src/meson/pick-abi.py does not exist" >&2
    exit 1
fi

# Build a minimal ABI table missing the bits_per_long key.
cat > "$TMPDIR/bad-table.ini" <<'EOF'
[properties]
arith_right_shift = 1
fma_policy = auto
long_double = target_native
x86_specializations = true
tune_table = x86
exec_mode = native
exe_wrapper =
shlib_style = elf
threading = pthread
tls_hack = false
EOF

if python3 "$REPO_ROOT/src/meson/pick-abi.py" --abi-file "$TMPDIR/bad-table.ini" \
        --triplet x86_64-linux-gnu \
        > "$TMPDIR/out" 2> "$TMPDIR/err"; then
    echo "FAIL: pick-abi.py succeeded on a table missing bits_per_long" >&2
    cat "$TMPDIR/out" >&2
    exit 1
fi

if ! grep -q 'bits_per_long' "$TMPDIR/err"; then
    echo "FAIL: error did not name the missing 'bits_per_long' key" >&2
    cat "$TMPDIR/err" >&2
    exit 1
fi

echo "PASS: T010 pick-abi rejects missing key"
