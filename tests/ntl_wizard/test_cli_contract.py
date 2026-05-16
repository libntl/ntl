"""T008 — CLI contract tests (US3).

Validates `ntl-wizard`'s CLI surface against `contracts/ntl-wizard-cli.md`:
exit codes, --version/--help behavior, no-TTY refusal, --batch flow,
--dry-run flow. RED until T022 (cli.py argparse) and T025 (__main__)
land.
"""
from __future__ import annotations

import pytest

# Exit code constants from contracts/ntl-wizard-cli.md
EXIT_OK = 0
EXIT_GENERIC_ERROR = 1
EXIT_CROSS_REFUSAL = 2
EXIT_COMPILE_FAILURE = 3
EXIT_RUNTIME_FAILURE = 4
EXIT_MEASUREMENT_NOISE = 5
EXIT_SESSION_CONFLICT = 6
EXIT_USER_INTERRUPT = 130


def test_version_flag_exits_zero(run_wizard):
    """`--version` MUST exit 0 and print to stdout."""
    code, stdout, stderr = run_wizard(["--version"])
    assert code == EXIT_OK, f"--version exit code: expected 0, got {code}; stderr={stderr!r}"
    assert stdout.strip(), "--version produced no stdout"


def test_help_flag_exits_zero(run_wizard):
    """`--help` MUST exit 0, print to stdout, and mention key flags."""
    code, stdout, stderr = run_wizard(["--help"])
    assert code == EXIT_OK
    assert "--batch" in stdout
    assert "--dry-run" in stdout
    assert "--resume" in stdout


def test_no_tty_default_mode_refuses(run_wizard):
    """Invoking the Wizard with no TTY and without --batch MUST refuse
    with a clear error suggesting --batch. (Subprocess context never
    has a controlling TTY, so this is the natural CI assertion.)"""
    code, stdout, stderr = run_wizard([])  # no --batch
    assert code != EXIT_OK, "Wizard with no TTY and no --batch should refuse"
    assert "--batch" in stderr or "--batch" in stdout, (
        f"Refusal message should suggest --batch. stdout={stdout!r}, stderr={stderr!r}"
    )


def test_dry_run_exits_zero_in_batch_mode(run_wizard, fake_ntl_source):
    """`--batch --dry-run --ntl-source-dir=<fake-src>` MUST validate
    setup and exit 0 without measuring anything."""
    code, stdout, stderr = run_wizard(
        ["--batch", "--dry-run", f"--ntl-source-dir={fake_ntl_source}"]
    )
    assert code == EXIT_OK, (
        f"--batch --dry-run exit: expected 0, got {code}; stderr={stderr!r}"
    )
    # stdout should contain at least one platform check line
    assert "platform" in stdout.lower() or "platform" in stderr.lower(), (
        "dry-run output should include the platform check result"
    )


def test_status_flag_exits_zero(run_wizard):
    """`--status` MUST exit 0 and print machine-readable JSON (per
    contract). When there is no prior session, the output is still
    valid JSON (e.g. `{"sessions": []}`)."""
    import json as _json
    code, stdout, _stderr = run_wizard(["--status"])
    assert code == EXIT_OK
    # stdout MUST be parseable JSON
    try:
        _json.loads(stdout)
    except _json.JSONDecodeError as exc:
        pytest.fail(f"--status stdout must be valid JSON; got {stdout!r}; error: {exc}")


def test_cross_target_refusal_exit_code(run_wizard):
    """Asserting a target architecture different from the build host
    MUST result in exit code 2 with the cross-vs-native error."""
    code, _stdout, stderr = run_wizard(
        ["--batch", "--target=nosuchcpu-unknown-elf"]
    )
    # exit 2 specifically, NOT a generic exit 1
    assert code == EXIT_CROSS_REFUSAL, (
        f"Expected exit code 2 (CROSS_REFUSAL); got {code}; stderr={stderr!r}"
    )
    # The error message MUST recommend the static-tune-table fallback
    assert "tune=" in stderr.lower() or "tune=" in stderr.lower(), (
        f"Cross-refusal stderr should suggest `-Dtune=...`. Got: {stderr!r}"
    )


def test_stderr_contains_error_marker_on_failure(run_wizard):
    """Per CLI contract: stderr on any error MUST contain an `error:`
    line (or, for argparse-style flag rejection, a `Usage:` + `Error:`
    pair). Trigger by passing an unknown flag."""
    code, _stdout, stderr = run_wizard(["--this-flag-does-not-exist"])
    assert code != EXIT_OK
    lowered = stderr.lower()
    assert (
        "error:" in lowered
        or "usage:" in lowered
        or "no such option" in lowered
    ), (
        f"stderr should contain an error marker (`error:`, `Usage:`, or "
        f"`no such option`). Got: {stderr!r}"
    )


def test_batch_mode_stdout_is_line_oriented(run_wizard, fake_ntl_source):
    """Batch mode stdout MUST be line-oriented and parseable. Each
    log line follows `[<utc-timestamp>] <level> <message>`. Verify
    via dry-run."""
    import re
    code, stdout, _stderr = run_wizard(
        ["--batch", "--dry-run", f"--ntl-source-dir={fake_ntl_source}"]
    )
    assert code == EXIT_OK
    timestamp_line_re = re.compile(r"^\[\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z?\] \w+ ")
    matching_lines = [
        line for line in stdout.splitlines() if timestamp_line_re.match(line)
    ]
    assert matching_lines, (
        f"Batch-mode stdout should contain at least one timestamped "
        f"log line. Got stdout={stdout!r}"
    )
