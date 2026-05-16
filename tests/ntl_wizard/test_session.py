"""T011 — Session persistence tests (US3, data-model.md).

Verifies `ntl_wizard.session.WizardSession`:
- JSON persistence under $XDG_CACHE_HOME / NTL_WIZARD_CACHE_DIR
- atomic write (no partial file at the session path)
- host_fingerprint mismatch invalidates resume
- sessions > 7 days old are auto-discarded

RED until T019 (session.py) lands.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import pytest


def test_session_module_exists():
    try:
        from ntl_wizard import session  # noqa: F401
    except ImportError as exc:
        pytest.fail(f"ntl_wizard.session not importable: {exc}")


def test_create_session_writes_to_cache_dir(tmp_cache_dir):
    from ntl_wizard.session import WizardSession
    s = WizardSession.create(host_fingerprint="sha256:test-fingerprint", mode="batch")
    s.save()
    # The session file should be findable under the cache dir
    found = list(tmp_cache_dir.rglob("session-*.json"))
    assert found, f"No session file under {tmp_cache_dir}"
    raw = json.loads(found[0].read_text(encoding="utf-8"))
    assert raw["host_fingerprint"] == "sha256:test-fingerprint"
    assert raw["mode"] == "batch"
    assert raw["session_id"]


def test_resume_same_fingerprint_succeeds(tmp_cache_dir):
    from ntl_wizard.session import WizardSession
    s = WizardSession.create(host_fingerprint="sha256:abc", mode="tui")
    s.save()
    resumed = WizardSession.try_resume(host_fingerprint="sha256:abc")
    assert resumed is not None
    assert resumed.session_id == s.session_id


def test_resume_different_fingerprint_returns_none(tmp_cache_dir):
    from ntl_wizard.session import WizardSession
    s = WizardSession.create(host_fingerprint="sha256:abc", mode="tui")
    s.save()
    resumed = WizardSession.try_resume(host_fingerprint="sha256:xyz")
    assert resumed is None, (
        "Resume across host_fingerprint mismatch must return None / not match"
    )


def test_old_session_is_auto_discarded(tmp_cache_dir):
    """Sessions older than 7 days are not eligible for resume."""
    from ntl_wizard.session import WizardSession
    s = WizardSession.create(host_fingerprint="sha256:abc", mode="batch")
    # Hand-set started_at to 10 days ago
    ten_days_ago = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=10)
    s.started_at = ten_days_ago
    s.save()
    resumed = WizardSession.try_resume(host_fingerprint="sha256:abc")
    assert resumed is None, (
        "Session > 7 days old should not be resumable"
    )


def test_session_write_is_atomic(tmp_cache_dir):
    """Many concurrent saves should never leave a partial file at the
    target session path (atomic via tempfile + os.replace)."""
    from ntl_wizard.session import WizardSession
    import threading
    s = WizardSession.create(host_fingerprint="sha256:atomic-test", mode="batch")

    observations: list[bool] = []
    stop = threading.Event()

    session_file = tmp_cache_dir / f"session-{s.host_fingerprint}.json"

    def watcher() -> None:
        while not stop.is_set():
            if session_file.exists():
                try:
                    json.loads(session_file.read_text(encoding="utf-8"))
                    observations.append(True)
                except (json.JSONDecodeError, OSError):
                    observations.append(False)

    t = threading.Thread(target=watcher, daemon=True)
    t.start()
    try:
        for _ in range(20):
            s.save()
    finally:
        stop.set()
        t.join(timeout=2)
    assert all(observations), (
        "Observed a partial / non-JSON file at the session path; "
        "session save is not atomic."
    )
