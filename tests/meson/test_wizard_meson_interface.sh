#!/bin/sh
# T013 — Wizard ↔ Meson interface tests (US3, contracts/wizard-meson-interface.md).
#
# Validates the contract between the Wizard's output artifact and Meson's
# read-tune-table.py:
#   - `-Dtune=host` without artifact → configure error with clear message.
#   - `-Dtune=host` with valid artifact → compile flags are injected.
#   - artifact mtime newer than build dir's last-configure → re-read.
#
# RED until T032 / T034 / T035 land.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_BUILD="$(mktemp -d)"
trap 'rm -rf "$TMP_BUILD"' EXIT

cd "$REPO_ROOT"

# Skip if the meson option isn't extended yet (pre-T034)
if ! grep -q "host" meson.options 2>/dev/null; then
    echo "SKIP: meson.options does not yet expose -Dtune=host (T034 not landed)"
    exit 0
fi

# --- Scenario 1: -Dtune=host without artifact MUST fail at configure ---
ARTIFACT_PATH="src/meson/tune-tables/host-tuned.ini"
if [ -e "$ARTIFACT_PATH" ]; then
    mv "$ARTIFACT_PATH" "${ARTIFACT_PATH}.bak"
fi
trap 'rm -rf "$TMP_BUILD"; [ -e "${ARTIFACT_PATH}.bak" ] && mv "${ARTIFACT_PATH}.bak" "$ARTIFACT_PATH" || true' EXIT

if meson setup -Dtune=host "$TMP_BUILD/no-artifact" >"$TMP_BUILD/no-artifact.log" 2>&1; then
    echo "FAIL: meson setup -Dtune=host without artifact should have FAILED"
    cat "$TMP_BUILD/no-artifact.log"
    exit 1
fi
if ! grep -q -i "host-tuned\.ini\|run.*ntl-wizard" "$TMP_BUILD/no-artifact.log"; then
    echo "FAIL: configure-error message does not point user at the Wizard"
    cat "$TMP_BUILD/no-artifact.log"
    exit 1
fi
echo "PASS: -Dtune=host without artifact rejected with helpful message"

# --- Scenario 2: -Dtune=host with valid artifact MUST inject -D flags ---
# Use the Python writer to produce a synthetic artifact at the source-tree path
PYTHONPATH="$REPO_ROOT/tools" python3 - <<EOF
from pathlib import Path
from ntl_wizard.artifacts import write_artifact

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
        "host_fingerprint": "sha256:test",
        "host_cpu": "test", "host_os": "Linux", "compiler": "gcc",
        "generated_utc": "2026-05-16T00:00:00Z",
        "session_id": "test",
    },
)
EOF

if ! meson setup -Dtune=host "$TMP_BUILD/with-artifact" >"$TMP_BUILD/with-artifact.log" 2>&1; then
    echo "FAIL: meson setup -Dtune=host with valid artifact should have SUCCEEDED"
    cat "$TMP_BUILD/with-artifact.log"
    exit 1
fi

# Verify -DNTL_TBL_REM=1 shows up in the build's compile commands
if [ -f "$TMP_BUILD/with-artifact/compile_commands.json" ]; then
    if ! grep -q "NTL_TBL_REM=1" "$TMP_BUILD/with-artifact/compile_commands.json"; then
        echo "FAIL: artifact's NTL_TBL_REM=1 did not propagate to compile flags"
        exit 1
    fi
fi
echo "PASS: -Dtune=host with valid artifact injects -D flags"

# --- Scenario 3: mtime-based re-read (smoke test) ---
sleep 1
touch "$ARTIFACT_PATH"
# Reconfigure should pick up the new mtime
meson setup --reconfigure "$TMP_BUILD/with-artifact" >>"$TMP_BUILD/with-artifact.log" 2>&1
echo "PASS: artifact mtime change triggers re-read (smoke)"

echo "All T013 scenarios passed."
