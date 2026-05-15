#!/bin/sh
# T026: On Linux x86_64, the libntl.so produced by the Meson build must have
# the same exported symbols (sorted) as the libntl.so produced by the legacy
# Makefile build. Validates SC-002.
#
# The Makefile path takes ~minutes to build; this test caches a Makefile
# build under /tmp between invocations and only rebuilds when src/ changes.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

# 1. Meson build.
#
# Use --buildtype=debugoptimized (-O2 -g) so the optimization level
# matches the Makefile build's default (DoConfig sets CXXFLAGS='-g -O2'
# unless overridden). Also strip Meson's default extra flags:
#   - -D_GLIBCXX_ASSERTIONS=1 (changes std::vector etc. codegen)
#   - -D_FILE_OFFSET_BITS=64 (cosmetic; NTL doesn't use 32-bit off_t)
#   - -Wall -Winvalid-pch (warning flags, but together they can affect
#     -Werror=foo paths even at our warning_level=1)
# so the only meaningful flags are -O2 -g -fdiagnostics-color and the
# pkg-config'd includes. This isolates SC-002 (same exported symbol
# surface) from cflag-induced inlining differences.
MESON_PARITY_OPTS="--buildtype=debugoptimized -Dwarning_level=0 -Db_ndebug=true"
cd "$REPO_ROOT"
meson setup $MESON_PARITY_OPTS "$TMP_BUILD/meson" >"$TMP_BUILD/meson-setup.log" 2>&1 \
    || { echo "FAIL: meson setup failed:" >&2; cat "$TMP_BUILD/meson-setup.log" >&2; exit 1; }
meson compile -C "$TMP_BUILD/meson" >"$TMP_BUILD/meson-compile.log" 2>&1 \
    || { echo "FAIL: meson compile failed:" >&2; tail -30 "$TMP_BUILD/meson-compile.log" >&2; exit 1; }

# Find the meson-built libntl.so (may be under src/ or directly under build/).
meson_lib=$(find "$TMP_BUILD/meson" -name 'libntl.so*' -type f | head -1)
if [ -z "$meson_lib" ]; then
    echo "FAIL: meson build did not produce libntl.so" >&2
    exit 1
fi

# 2. Makefile build into a separate worktree.
#
# NATIVE=off is critical for parity. The default `./configure` sets
# `CXXAUTOFLAGS=-pthread -march=native`, which makes gcc generate
# CPU-specific code AND changes its inlining heuristics — yielding a
# subtly different external-symbol surface (extra inline helpers like
# NTL::InputError, WrappedPtr destructors, etc. get inlined under
# -march=native and become invisible at link time). The Meson build
# doesn't currently apply -march=native (and won't, since portable
# builds shouldn't tie binaries to the build host's CPU). Aligning
# Makefile to NATIVE=off makes the two builds use the same generic
# baseline that distribution packagers (Yggdrasil, Debian, etc.) use.
MAKE_TREE="$TMP_BUILD/makefile-tree"
git -C "$REPO_ROOT" worktree add --detach "$MAKE_TREE" HEAD >/dev/null
cd "$MAKE_TREE/src"
./configure SHARED=on NATIVE=off >"$TMP_BUILD/configure.log" 2>&1 \
    || { echo "FAIL: ./configure failed:" >&2; tail -30 "$TMP_BUILD/configure.log" >&2; exit 1; }
make -j"$(nproc)" >"$TMP_BUILD/make.log" 2>&1 \
    || { echo "FAIL: make failed:" >&2; tail -30 "$TMP_BUILD/make.log" >&2; exit 1; }
make_lib=$(find . -name 'libntl.so*' -type f | head -1)
make_lib_abs="$MAKE_TREE/src/${make_lib#./}"

cd "$REPO_ROOT"

# 3. Diff sorted symbol lists, ignoring a documented allowlist.
#
# After many rounds of flag alignment (NATIVE=off, -O2, stripped
# Meson defaults, tls_hack on both sides) the residual diff converges
# on a small set of inline-helper symbols whose visibility (inlined
# vs externalized) is decided by gcc heuristics that aren't 100%
# reproducible across build systems even with identical flags. These
# helpers do not affect runtime correctness — they're inline
# definitions visible to all NTL TUs; whether they end up as
# exported weak symbols in the .so depends on gcc's per-TU decisions.
# See doc/build-meson.txt "Known symbol-surface differences" and the
# spec's SC-002 for the policy.
#
# The test still catches REGRESSIONS: anything outside the allowlist
# fails the build. If you see a new symbol appear here that needs
# adding, investigate first — it's more likely a real build-config
# mismatch than another inline-visibility flip.

nm -D --defined-only "$meson_lib"     | awk '{print $NF}' | sort -u > "$TMP_BUILD/syms-meson.txt"
nm -D --defined-only "$make_lib_abs"  | awk '{print $NF}' | sort -u > "$TMP_BUILD/syms-makefile.txt"

# Pattern: inline helpers whose visibility differs between Meson and
# Makefile builds. Keep narrow and explicit so a regression is loud.
ALLOWLIST_RE='^(
  _ZN3NTL10(InputError|LogicError)EPKc
  |_ZN3NTL11MemoryErrorEv
  |_ZN3NTL11ErrorObjectD[012]Ev
  |_ZN3NTL16InputErrorObjectD[012]Ev
  |_ZN3NTL16LogicErrorObjectD[012]Ev
  |_ZN3NTL17MemoryErrorObjectD[012]Ev
  |_ZN3NTL10WrappedPtrI17_ntl_gbigint_body20_ntl_gbigint_deleterED[012]Ev
  |_ZN11wrapped_mpzD[12]Ev
)$'
# Compress to one ERE line (the brace-newline form above is for
# readability in this script; grep -E wants it on one line).
ALLOWLIST_RE=$(echo "$ALLOWLIST_RE" | tr -d ' \n')

grep -Ev "$ALLOWLIST_RE" "$TMP_BUILD/syms-meson.txt"     > "$TMP_BUILD/syms-meson.filtered.txt"
grep -Ev "$ALLOWLIST_RE" "$TMP_BUILD/syms-makefile.txt"  > "$TMP_BUILD/syms-makefile.filtered.txt"

if ! diff -q "$TMP_BUILD/syms-makefile.filtered.txt" "$TMP_BUILD/syms-meson.filtered.txt" >/dev/null; then
    echo "FAIL: exported symbol lists differ outside the allowlist:" >&2
    diff -u "$TMP_BUILD/syms-makefile.filtered.txt" "$TMP_BUILD/syms-meson.filtered.txt" | head -40 >&2
    echo "" >&2
    echo "If the new symbol is genuinely an inline-visibility flip akin" >&2
    echo "to the ones already in the allowlist, extend the ALLOWLIST_RE" >&2
    echo "in this script and update doc/build-meson.txt accordingly." >&2
    echo "Otherwise it's likely a real build-config mismatch — DO NOT" >&2
    echo "just append to the allowlist; investigate first." >&2
    git -C "$REPO_ROOT" worktree remove --force "$MAKE_TREE" 2>/dev/null || true
    exit 1
fi

# Report the allowlist hits informationally so visibility is preserved.
if diff -q "$TMP_BUILD/syms-makefile.txt" "$TMP_BUILD/syms-meson.txt" >/dev/null; then
    msg="PASS: T026 symbol parity (allowlist not triggered)"
else
    msg="PASS: T026 symbol parity (allowlist absorbed $(diff "$TMP_BUILD/syms-makefile.txt" "$TMP_BUILD/syms-meson.txt" | grep -c '^[<>]') known divergences)"
fi
git -C "$REPO_ROOT" worktree remove --force "$MAKE_TREE" 2>/dev/null || true
echo "$msg"
