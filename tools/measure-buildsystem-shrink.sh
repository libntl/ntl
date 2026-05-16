#!/bin/sh
# T058 — Measure the build-system source-line reduction from feature 002.
# Compares the legacy (HEAD~1 or first arg) tree to the current HEAD.
#
# Output: a one-page report on stdout.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

BEFORE_REF="${1:-HEAD~1}"
AFTER_REF="${2:-HEAD}"

# Files that count as "build-system source". Excludes NTL's C++ library
# code itself (src/*.cpp, src/*.h, src/lzz_*.cpp, etc.); keeps only what
# orchestrates the build.
declare_pattern() {
    cat <<'EOF'
src/configure
src/DoConfig
src/Makefile
src/mfile
src/cfile
src/Wizard
src/WizardAux
src/TestScript
src/CopyFeatures
src/MakeDescAux.cpp
meson.build
meson.options
src/meson.build
src/NTL/meson.build
src/meson/pick-abi.py
src/meson/gen-have-headers.py
src/meson/gen-gmp-aux.py
src/meson/run-makedesc.py
src/meson/run-golden-test.sh
src/meson/read-tune-table.py
tools/sync-sources.py
tools/check-sources-in-sync.py
tools/check-cfile-in-sync.py
tools/sync-version.py
tools/check-commit-trailer.sh
tools/measure-buildsystem-shrink.sh
tools/ntl_wizard/__init__.py
tools/ntl_wizard/__main__.py
tools/ntl_wizard/parameters.py
tools/ntl_wizard/platform_check.py
tools/ntl_wizard/artifacts.py
tools/ntl_wizard/session.py
tools/ntl_wizard/measure.py
tools/ntl_wizard/search.py
tools/ntl_wizard/cli.py
tools/ntl_wizard/app.py
tools/pyproject.toml
EOF
}

count_at_ref() {
    ref="$1"
    total=0
    while IFS= read -r path; do
        [ -z "$path" ] && continue
        if git show "$ref:$path" >/dev/null 2>&1; then
            lines=$(git show "$ref:$path" | wc -l)
            total=$((total + lines))
        fi
    done
    echo "$total"
}

before=$(declare_pattern | count_at_ref "$BEFORE_REF")
after=$(declare_pattern | count_at_ref "$AFTER_REF")

if [ "$before" -eq 0 ]; then
    echo "Warning: 0 lines counted at $BEFORE_REF. May be a fresh repo."
fi

delta=$((before - after))
if [ "$before" -gt 0 ]; then
    pct=$(echo "scale=1; $delta * 100 / $before" | bc 2>/dev/null || echo "?")
else
    pct="?"
fi

printf 'Build-system source line count\n'
printf '==============================\n\n'
printf '  Before (%s): %d lines\n' "$BEFORE_REF" "$before"
printf '  After  (%s): %d lines\n' "$AFTER_REF" "$after"
printf '  Delta:           %d lines (%s%%)\n' "$delta" "$pct"
printf '\n'
if [ "$pct" != "?" ]; then
    threshold=80
    pct_int=$(printf '%.0f' "$pct" 2>/dev/null || echo 0)
    if [ "$pct_int" -ge "$threshold" ]; then
        printf 'SC-004 satisfied: build-system source dropped >= %d%%.\n' "$threshold"
    else
        printf 'SC-004 NOT satisfied: %d%% reduction is below %d%% threshold.\n' "$pct_int" "$threshold"
    fi
fi
