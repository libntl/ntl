"""T009 — Tune-artifact writer tests (US3, contracts/tune-table-schema.md).

Verifies `ntl_wizard.artifacts.write_artifact()`:
- atomic write (no partial file ever appears at target path),
- LF line endings,
- `[parameters]` section sorted by TunableParameter declaration order,
- `[provenance]` populated correctly,
- reader (`src/meson/read-tune-table.py`) rejects malformed artifacts.

RED until T018 (artifacts.py) and T032 (read-tune-table.py) land.
"""
from __future__ import annotations

import configparser
import os
import threading
from pathlib import Path

import pytest


def _minimal_parameters() -> dict[str, int | bool | str]:
    """A minimal valid parameter set keyed by name."""
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


def _minimal_provenance() -> dict[str, str]:
    return {
        "ntl_version": Path(__file__).resolve().parents[2].joinpath("version.txt").read_text().strip(),
        "wizard_version": Path(__file__).resolve().parents[2].joinpath("version.txt").read_text().strip(),
        "host_fingerprint": "sha256:" + "0" * 64,
        "host_cpu": "Generic CPU",
        "host_os": "Linux test",
        "compiler": "gcc 13.2.0",
        "generated_utc": "2026-05-16T14:32:11Z",
        "session_id": "deadbeef-test",
    }


def test_artifact_module_exists():
    try:
        from ntl_wizard import artifacts  # noqa: F401
    except ImportError as exc:
        pytest.fail(f"ntl_wizard.artifacts not importable: {exc}")


def test_write_artifact_produces_parseable_ini(tmp_artifact_path: Path):
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    assert tmp_artifact_path.exists()
    parser = configparser.ConfigParser()
    parser.read(tmp_artifact_path, encoding="utf-8")
    assert "parameters" in parser.sections()
    assert "provenance" in parser.sections()
    # Round-trip every parameter key
    for k, v in _minimal_parameters().items():
        assert parser["parameters"][k] == str(v)


def test_write_artifact_is_atomic(tmp_artifact_path: Path):
    """The writer MUST never leave a partial file at the artifact path.

    Strategy: race the writer against an inotify-style watcher and
    assert that at every observed mtime, the file at the artifact
    path is either non-existent OR fully valid (has [parameters] and
    [provenance] sections).
    """
    from ntl_wizard.artifacts import write_artifact

    observations: list[bool] = []  # True = valid INI, False = partial / missing
    stop = threading.Event()

    def watcher() -> None:
        while not stop.is_set():
            if tmp_artifact_path.exists():
                try:
                    parser = configparser.ConfigParser()
                    parser.read(tmp_artifact_path, encoding="utf-8")
                    observations.append(
                        "parameters" in parser.sections()
                        and "provenance" in parser.sections()
                    )
                except (configparser.Error, OSError):
                    observations.append(False)

    t = threading.Thread(target=watcher, daemon=True)
    t.start()
    try:
        # Write a few times to maximize race window
        for _ in range(5):
            write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    finally:
        stop.set()
        t.join(timeout=2)

    # All observations of the file at the artifact path MUST show a
    # valid INI. The atomic-rename property guarantees this.
    assert all(observations), (
        "Observed a partial file at the artifact path — "
        "write_artifact is not atomic."
    )


def test_write_artifact_uses_lf_line_endings(tmp_artifact_path: Path):
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    raw = tmp_artifact_path.read_bytes()
    assert b"\r\n" not in raw, "Artifact contains CRLF; expected LF-only."


def test_parameters_section_in_declaration_order(tmp_artifact_path: Path):
    """`[parameters]` keys MUST follow the TunableParameter declaration
    order (from ntl_wizard.parameters.PARAMETERS), NOT alphabetical.
    Keeps diffs reviewable across re-runs."""
    from ntl_wizard import parameters as pmod
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    raw = tmp_artifact_path.read_text(encoding="utf-8")

    declared_order = [p.name for p in pmod.PARAMETERS]
    # Find each key's first occurrence
    positions = {k: raw.find(k + " ") for k in declared_order if (k + " ") in raw}
    # Restrict to keys present in the artifact
    found = [k for k in declared_order if positions.get(k, -1) >= 0]
    seen_positions = [positions[k] for k in found]
    assert seen_positions == sorted(seen_positions), (
        f"[parameters] keys should appear in PARAMETERS declaration order; "
        f"got order: {found}, positions: {seen_positions}"
    )


def test_provenance_required_fields(tmp_artifact_path: Path):
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    parser = configparser.ConfigParser()
    parser.read(tmp_artifact_path, encoding="utf-8")
    p = parser["provenance"]
    for required in ("ntl_version", "wizard_version", "host_fingerprint", "generated_utc"):
        assert required in p, f"Missing required provenance key: {required}"


def test_reader_rejects_missing_key(tmp_artifact_path: Path):
    """Meson-side reader MUST refuse to consume an artifact whose
    [parameters] section is missing a required key.

    The writer is strict — it would refuse to emit such a file. So
    we hand-craft the malformed artifact here to exercise the reader's
    backward-compat error path."""
    # Write a valid artifact first to set up [provenance] correctly
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())
    # Strip one [parameters] key manually
    raw = tmp_artifact_path.read_text(encoding="utf-8")
    raw = "\n".join(
        line for line in raw.splitlines() if not line.startswith("NTL_TBL_REM ")
    )
    tmp_artifact_path.write_text(raw, encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[2]
    reader = repo_root / "src" / "meson" / "read-tune-table.py"
    if not reader.exists():
        pytest.skip("Reader script not yet implemented (T032)")

    import subprocess
    result = subprocess.run(
        ["python3", str(reader), str(tmp_artifact_path)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0, (
        f"Reader should reject artifact missing required keys. "
        f"stdout={result.stdout!r}, stderr={result.stderr!r}"
    )


def test_reader_warns_on_unknown_key(tmp_artifact_path: Path):
    """Forward-compat: an artifact with an unknown extra key in the
    [parameters] section emits a warning to stderr but does not fail.
    Allows older Meson + newer Wizard."""
    from ntl_wizard.artifacts import write_artifact
    write_artifact(tmp_artifact_path, _minimal_parameters(), _minimal_provenance())

    # Inject a stray key directly under [parameters] (right after the
    # section header). The writer is strict so we hand-edit here.
    raw = tmp_artifact_path.read_text(encoding="utf-8")
    raw = raw.replace(
        "[parameters]\n",
        "[parameters]\nNTL_FUTURE_TUNABLE = 42\n",
        1,
    )
    tmp_artifact_path.write_text(raw, encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[2]
    reader = repo_root / "src" / "meson" / "read-tune-table.py"
    if not reader.exists():
        pytest.skip("Reader script not yet implemented (T032)")

    import subprocess
    result = subprocess.run(
        ["python3", str(reader), str(tmp_artifact_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        f"Reader should ACCEPT unknown extra keys (forward-compat). "
        f"stderr={result.stderr!r}"
    )
    assert "warn" in result.stderr.lower() or "unknown" in result.stderr.lower(), (
        f"Reader should emit warning on unknown key; stderr={result.stderr!r}"
    )
