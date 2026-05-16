"""T027 — Minimal end-to-end integration test for ntl-wizard.

Excluded from default pytest run via the `@pytest.mark.slow` marker
(see conftest.py); CI runs it in the scheduled-weekly tier only
(research §R7 tier 3).

Asserts: invoking `ntl-wizard --batch --dry-run` against the real
NTL source tree completes successfully, prints valid stdout lines,
and exits 0. The full measure-and-write integration is left to the
non-gating weekly CI job.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


pytestmark = pytest.mark.slow


def test_wizard_dry_run_on_real_source(run_wizard, repo_root):
    """A dry-run against the actual NTL source tree at repo_root
    completes in seconds, produces a structured stdout, and exits 0."""
    code, stdout, stderr = run_wizard(
        ["--batch", "--dry-run", f"--ntl-source-dir={repo_root}"]
    )
    assert code == 0, (
        f"--batch --dry-run on real source should exit 0. "
        f"Got {code}. stderr={stderr!r}"
    )
    # Every phase's timing program should be listed as OK
    for phase_program in (
        "Poly1TimeTest.cpp",
        "Poly2TimeTest.cpp",
        "Poly3TimeTest.cpp",
        "GF2XTimeTest.cpp",
    ):
        assert phase_program in stdout, (
            f"dry-run output should mention each timing program. "
            f"Missing {phase_program}. stdout={stdout!r}"
        )


def test_wizard_status_returns_valid_json(run_wizard, tmp_cache_dir):
    """`ntl-wizard --status` returns valid JSON describing zero or
    more sessions. The fixture's empty cache dir means we expect
    zero sessions for this run."""
    code, stdout, _stderr = run_wizard(
        ["--status"],
        env_overrides={"NTL_WIZARD_CACHE_DIR": str(tmp_cache_dir)},
    )
    assert code == 0
    payload = json.loads(stdout)
    assert "sessions" in payload
    assert payload["sessions"] == []


def test_wizard_refuses_cross_target(run_wizard):
    """Asserting a wildly different target arch must produce exit 2
    and a stderr that points at the static-tune-table fallback."""
    code, _stdout, stderr = run_wizard(
        ["--batch", "--target=riscv64-linux-gnu"],
    )
    assert code == 2, f"Expected EXIT_CROSS_REFUSAL=2; got {code}"
    assert "-Dtune=" in stderr or "tune=" in stderr, (
        f"Cross-refusal should suggest static tune table. stderr={stderr!r}"
    )
