#!/bin/sh
# Run an NTL "golden-diff" test: stdin from <name>In, stdout compared to <name>Out
# via `diff -b`. Mirrors src/TestScript's behavior.
#
# Usage: run-golden-test.sh <test-program> <input-file> <expected-output-file>

set -eu

if [ "$#" -ne 3 ]; then
    echo "usage: run-golden-test.sh <prog> <input> <expected>" >&2
    exit 2
fi

prog="$1"
input="$2"
expected="$3"

tmp_out=$(mktemp)
trap 'rm -f "$tmp_out"' EXIT

if ! "$prog" < "$input" > "$tmp_out" 2>&1; then
    echo "FAIL: $prog exited non-zero. Output:" >&2
    cat "$tmp_out" >&2
    exit 1
fi

if ! diff -b "$tmp_out" "$expected"; then
    echo "FAIL: $prog output does not match $expected" >&2
    exit 1
fi

echo "PASS: $(basename "$prog")"
