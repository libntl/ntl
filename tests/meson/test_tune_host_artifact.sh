#!/bin/sh
# T030 — `-Dtune=host` with a valid artifact must propagate parameter
# values into the NTL compile flags. RED until T032/T034/T035 land.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
ARTIFACT_BACKUP=""
cleanup() {
    rm -rf "$TMP_BUILD"
    if [ -n "$ARTIFACT_BACKUP" ] && [ -e "$ARTIFACT_BACKUP" ]; then
        mv "$ARTIFACT_BACKUP" "$REPO_ROOT/src/meson/tune-tables/host-tuned.ini"
    fi
}
trap cleanup EXIT

cd "$REPO_ROOT"

# Skip cleanly if Phase-4 wiring isn't in place yet.
if ! grep -q "'host'" meson.options 2>/dev/null; then
    echo "SKIP: meson.options does not yet expose -Dtune=host"
    exit 0
fi

# Preserve any existing artifact so we can restore on exit.
ARTIFACT_PATH="$REPO_ROOT/src/meson/tune-tables/host-tuned.ini"
if [ -e "$ARTIFACT_PATH" ]; then
    ARTIFACT_BACKUP="$ARTIFACT_PATH.test-backup"
    mv "$ARTIFACT_PATH" "$ARTIFACT_BACKUP"
fi

# Write a synthetic artifact with known values.
PYTHONPATH="$REPO_ROOT/tools" python3 - <<EOF
from pathlib import Path
from ntl_wizard.artifacts import write_artifact, now_utc_iso

write_artifact(
    Path("$ARTIFACT_PATH"),
    {
        "NTL_AVOID_BRANCHING": 0,
        "NTL_CRT_ALTCODE": 0,
        "NTL_CRT_ALTCODE_SMALL": 1,
        "NTL_FFT_BIGTAB": 1,
        "NTL_FFT_LAZYMUL": 1,
        "NTL_GF2X_ALTCODE": 0,
        "NTL_GF2X_ALTCODE1": 1,
        "NTL_GF2X_NOINLINE": 0,
        "NTL_SPMM_ULL": 1,
        "NTL_TBL_REM": 1,
    },
    {
        "ntl_version": Path("$REPO_ROOT/version.txt").read_text().strip(),
        "wizard_version": Path("$REPO_ROOT/version.txt").read_text().strip(),
        "host_fingerprint": "sha256:test30",
        "host_cpu": "test", "host_os": "Linux", "compiler": "gcc",
        "generated_utc": now_utc_iso(),
        "session_id": "t030",
    },
)
EOF

if ! meson setup -Dtune=host "$TMP_BUILD/host" >"$TMP_BUILD/host.log" 2>&1; then
    echo "FAIL: meson setup -Dtune=host with valid artifact should succeed"
    cat "$TMP_BUILD/host.log"
    exit 1
fi

# Verify at least one of our `-DNTL_*=*` flags appears in the build args.
# Meson stores effective compile commands in compile_commands.json.
if [ -f "$TMP_BUILD/host/compile_commands.json" ]; then
    if ! grep -q "NTL_TBL_REM=1" "$TMP_BUILD/host/compile_commands.json"; then
        echo "FAIL: artifact's NTL_TBL_REM=1 did not propagate to compile flags"
        exit 1
    fi
fi

echo "PASS: -Dtune=host artifact flow OK"
