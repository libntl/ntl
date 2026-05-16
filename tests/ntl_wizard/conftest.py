"""pytest fixtures for the ntl-wizard test suite.

Pyramid-style fixtures: cheap unit-level fixtures at the top, slower
integration fixtures (real Meson sub-build) at the bottom. Tests
labelled `@pytest.mark.slow` are excluded from the default run.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import pytest


# --- Path resolution --------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS_DIR = REPO_ROOT / "tools"
SPECS_DIR = REPO_ROOT / "specs" / "002-remove-legacy-build"
CAPTURED_PARAMS_FILE = SPECS_DIR / "captured-legacy-params.txt"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the NTL repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def captured_legacy_params() -> list[str]:
    """The sorted list of legacy Wizard-tunable parameter names.

    Loaded from `specs/002-remove-legacy-build/captured-legacy-params.txt`
    (see T002). This is the ground truth for parameter-parity tests
    (T007 / FR-005a) — the test MUST fail if the Python Wizard does
    not cover exactly this set.
    """
    if not CAPTURED_PARAMS_FILE.exists():
        pytest.fail(
            f"Missing captured-legacy-params.txt at {CAPTURED_PARAMS_FILE}. "
            "Run T002 before any other Phase-2 test."
        )
    names: list[str] = []
    for raw in CAPTURED_PARAMS_FILE.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        names.append(line)
    return sorted(names)


# --- Filesystem fixtures ----------------------------------------------------

@pytest.fixture
def tmp_cache_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Ephemeral XDG-style cache dir for one test. Wizard's session
    persistence lands under here; nothing escapes the tmp tree.
    """
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setenv("NTL_WIZARD_CACHE_DIR", str(cache))
    # Also override XDG so any fallback codepath lands here too.
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    return cache


@pytest.fixture
def tmp_artifact_path(tmp_path: Path) -> Path:
    """Where a test Wizard run should write its host-tuned.ini.

    Default path inside a tmp tree so tests never touch the real
    `src/meson/tune-tables/host-tuned.ini`.
    """
    out_dir = tmp_path / "tune-tables"
    out_dir.mkdir()
    return out_dir / "host-tuned.ini"


@pytest.fixture
def fake_ntl_source(tmp_path: Path) -> Path:
    """A minimal fake NTL source tree just rich enough to satisfy the
    Wizard's "I can find my inputs" preflight checks. Real measurement
    requires the slow fixture below.
    """
    src = tmp_path / "fake-ntl"
    src.mkdir()
    (src / "version.txt").write_text("11.6.0\n")
    (src / "src").mkdir()
    for name in (
        "Poly1TimeTest.cpp",
        "Poly2TimeTest.cpp",
        "Poly3TimeTest.cpp",
        "GF2XTimeTest.cpp",
        "InitSettings.cpp",
        "DispSettings.cpp",
    ):
        (src / "src" / name).write_text(
            f"// stub for {name} — only present so the Wizard's preflight passes\n"
            "int main() { return 0; }\n"
        )
    (src / "meson.build").write_text(
        "project('ntl-fake', 'cpp', version: files('version.txt'))\n"
    )
    return src


# --- Mock measurement results -----------------------------------------------

@dataclass(frozen=True)
class FakeMeasurement:
    """Lightweight Measurement double for tests that exercise the
    parameter-search logic without running the real timing binaries.
    """
    parameter_set: dict[str, int | bool | str]
    wall_clock_seconds: float
    iteration_count: int


@pytest.fixture
def mock_measurements_poly1() -> list[FakeMeasurement]:
    """Plausible Poly1 timing results: the configuration with
    NTL_FFT_LAZYMUL=1 + NTL_SPMM_ULL=1 wins by a small margin.
    """
    return [
        FakeMeasurement({"NTL_FFT_LAZYMUL": 0, "NTL_SPMM_ULL": 0,
                         "NTL_AVOID_BRANCHING": 0}, 1.20, 1000),
        FakeMeasurement({"NTL_FFT_LAZYMUL": 1, "NTL_SPMM_ULL": 1,
                         "NTL_AVOID_BRANCHING": 0}, 0.95, 1000),
        FakeMeasurement({"NTL_FFT_LAZYMUL": 1, "NTL_SPMM_ULL": 1,
                         "NTL_AVOID_BRANCHING": 1}, 0.98, 1000),
    ]


# --- Slow / integration fixtures --------------------------------------------

@pytest.fixture(scope="session")
def have_meson() -> bool:
    """Whether the test runner has Meson + ninja available."""
    return bool(shutil.which("meson") and shutil.which("ninja"))


@pytest.fixture(scope="session")
def have_cxx_compiler() -> bool:
    """Whether the test runner has a working C++ compiler."""
    candidates = ("g++", "clang++", "c++")
    return any(shutil.which(c) for c in candidates)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip @pytest.mark.slow tests unless --run-slow is passed."""
    if config.getoption("--run-slow", default=False):
        return
    skip_slow = pytest.mark.skip(reason="slow test, pass --run-slow to enable")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run @pytest.mark.slow integration tests (real Meson sub-build).",
    )


# --- Subprocess invocation helper ------------------------------------------

@pytest.fixture
def run_wizard():
    """Invoke `ntl-wizard` (via python -m to avoid PATH issues) and
    return (returncode, stdout, stderr).
    """
    def _run(args: list[str], env_overrides: dict[str, str] | None = None,
             cwd: Path | None = None) -> tuple[int, str, str]:
        cmd = [sys.executable, "-m", "ntl_wizard", *args]
        env = os.environ.copy()
        if env_overrides:
            env.update(env_overrides)
        # Ensure the in-tree package is importable. We rely on the
        # caller pip-installing it OR on PYTHONPATH; the latter is
        # easier in CI sandboxes that may not allow pip install.
        env.setdefault("PYTHONPATH", str(TOOLS_DIR))
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            cwd=cwd or REPO_ROOT,
            timeout=120,
        )
        return result.returncode, result.stdout, result.stderr

    return _run
