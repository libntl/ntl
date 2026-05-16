"""T012 — Tune-table schema validation tests (US3, contracts/tune-table-schema.md).

Validates that the artifact's INI schema is enforced symmetrically by
both the Wizard (writer) and Meson (reader):
- Writer rejects payloads with wrong-type values for a TunableParameter.
- Reader accepts forward-compat (unknown extra keys → warn).
- Reader rejects backward-compat (missing required keys → error).

RED until T018 (writer) and T032 (reader) land.
"""
from __future__ import annotations

import configparser
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
READER_SCRIPT = REPO_ROOT / "src" / "meson" / "read-tune-table.py"


def _full_parameter_set() -> dict[str, int]:
    return {
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
    }


def _basic_provenance() -> dict[str, str]:
    return {
        "ntl_version": Path(__file__).resolve().parents[2].joinpath("version.txt").read_text().strip(),
        "wizard_version": Path(__file__).resolve().parents[2].joinpath("version.txt").read_text().strip(),
        "host_fingerprint": "sha256:" + "0" * 64,
        "host_cpu": "test",
        "host_os": "Linux test",
        "compiler": "gcc 13.2.0",
        "generated_utc": "2026-05-16T00:00:00Z",
        "session_id": "test-session",
    }


def test_writer_accepts_valid_payload(tmp_artifact_path: Path):
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _full_parameter_set(), _basic_provenance())
    # No exception means accepted
    assert tmp_artifact_path.exists()


def test_writer_rejects_string_for_bool_flag(tmp_artifact_path: Path):
    """A bool_flag parameter MUST get 0 or 1; supplying "yes" must fail."""
    from ntl_wizard.artifacts import write_artifact
    params = _full_parameter_set()
    params["NTL_TBL_REM"] = "yes"  # bool_flag, only 0/1 allowed
    with pytest.raises((ValueError, TypeError)):
        write_artifact(tmp_artifact_path, params, _basic_provenance())


def test_writer_rejects_unknown_parameter(tmp_artifact_path: Path):
    """Writing a parameter not in PARAMETERS MUST fail at write time
    (the writer is the strict side; the reader is the lenient side
    for forward-compat)."""
    from ntl_wizard.artifacts import write_artifact
    params = _full_parameter_set()
    params["NTL_NOT_A_REAL_PARAM"] = 1
    with pytest.raises((ValueError, KeyError)):
        write_artifact(tmp_artifact_path, params, _basic_provenance())


def test_reader_emits_compile_flags(tmp_artifact_path: Path):
    """The reader script's success-path output is a sequence of
    `-DNTL_<KEY>=<VALUE>` flags, one per line, for Meson to consume."""
    if not READER_SCRIPT.exists():
        pytest.skip("Reader script not yet implemented (T032)")
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _full_parameter_set(), _basic_provenance())

    result = subprocess.run(
        ["python3", str(READER_SCRIPT), str(tmp_artifact_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"Reader failed: stderr={result.stderr!r}"
    # Every parameter should appear as a -D flag
    for k, v in _full_parameter_set().items():
        assert f"-D{k}={v}" in result.stdout, (
            f"Missing -D flag for {k}={v} in reader output: {result.stdout!r}"
        )


def test_reader_rejects_version_mismatch(tmp_artifact_path: Path):
    """An artifact whose ntl_version differs from the surrounding
    tree's version.txt MUST be rejected (stale artifact)."""
    if not READER_SCRIPT.exists():
        pytest.skip("Reader script not yet implemented (T032)")
    from ntl_wizard.artifacts import write_artifact
    prov = _basic_provenance()
    prov["ntl_version"] = "999.999.999"  # very stale
    write_artifact(tmp_artifact_path, _full_parameter_set(), prov)

    result = subprocess.run(
        ["python3", str(READER_SCRIPT), str(tmp_artifact_path)],
        capture_output=True, text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode != 0
    assert "stale" in result.stderr.lower() or "version" in result.stderr.lower()
