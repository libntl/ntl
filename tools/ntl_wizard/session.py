"""Pause/resume session state for the Wizard.

Persisted as JSON under `${XDG_CACHE_HOME:-~/.cache}/ntl-wizard/`
(or `$NTL_WIZARD_CACHE_DIR` if set). One file per host_fingerprint.
Sessions older than 7 days are auto-discarded on resume attempts.
"""
from __future__ import annotations

import datetime as dt
import enum
import hashlib
import json
import os
import platform
import tempfile
import uuid
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any


class PhaseStatus(enum.Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    INTERRUPTED = "interrupted"


SESSION_TTL = dt.timedelta(days=7)
SESSION_FILE_PREFIX = "session-"


@dataclass
class PhaseState:
    """Per-phase tracking embedded inside a WizardSession."""
    status: PhaseStatus = PhaseStatus.PENDING
    started_at: dt.datetime | None = None
    completed_at: dt.datetime | None = None
    measurements: list[dict] = field(default_factory=list)
    derived_values: dict[str, Any] = field(default_factory=dict)
    error: str | None = None

    def to_dict(self) -> dict:
        return {
            "status": self.status.value,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
            "measurements": self.measurements,
            "derived_values": self.derived_values,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "PhaseState":
        return cls(
            status=PhaseStatus(raw["status"]),
            started_at=dt.datetime.fromisoformat(raw["started_at"]) if raw["started_at"] else None,
            completed_at=dt.datetime.fromisoformat(raw["completed_at"]) if raw["completed_at"] else None,
            measurements=list(raw.get("measurements", [])),
            derived_values=dict(raw.get("derived_values", {})),
            error=raw.get("error"),
        )


def cache_dir() -> Path:
    """Resolve the directory where session files live.

    Precedence: $NTL_WIZARD_CACHE_DIR > $XDG_CACHE_HOME/ntl-wizard > ~/.cache/ntl-wizard.
    """
    if (override := os.environ.get("NTL_WIZARD_CACHE_DIR")):
        return Path(override)
    if (xdg := os.environ.get("XDG_CACHE_HOME")):
        return Path(xdg) / "ntl-wizard"
    return Path.home() / ".cache" / "ntl-wizard"


def compute_host_fingerprint() -> str:
    """Stable hash of the build host's identity for matching session
    resumes. If anything in the (CPU + OS + Python version) tuple
    changes, the fingerprint changes and the session is invalidated.
    """
    parts = (
        platform.machine(),
        platform.system(),
        platform.release(),
        platform.python_version(),
    )
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


@dataclass
class WizardSession:
    """A single Wizard invocation's state. Persisted to disk for resume."""
    session_id: str
    host_fingerprint: str
    started_at: dt.datetime
    mode: str  # "tui" or "batch"
    output_artifact_path: str = "src/meson/tune-tables/host-tuned.ini"
    phases: dict[str, PhaseState] = field(default_factory=dict)

    @classmethod
    def create(
        cls,
        host_fingerprint: str,
        mode: str = "batch",
        output_artifact_path: str = "src/meson/tune-tables/host-tuned.ini",
    ) -> "WizardSession":
        """Fresh session with a new UUID."""
        return cls(
            session_id=str(uuid.uuid4()),
            host_fingerprint=host_fingerprint,
            started_at=dt.datetime.now(dt.timezone.utc),
            mode=mode,
            output_artifact_path=output_artifact_path,
            phases={},
        )

    def to_dict(self) -> dict:
        return {
            "session_id": self.session_id,
            "host_fingerprint": self.host_fingerprint,
            "started_at": self.started_at.isoformat(),
            "mode": self.mode,
            "output_artifact_path": self.output_artifact_path,
            "phases": {pid: ps.to_dict() for pid, ps in self.phases.items()},
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "WizardSession":
        return cls(
            session_id=raw["session_id"],
            host_fingerprint=raw["host_fingerprint"],
            started_at=dt.datetime.fromisoformat(raw["started_at"]),
            mode=raw["mode"],
            output_artifact_path=raw.get(
                "output_artifact_path", "src/meson/tune-tables/host-tuned.ini"
            ),
            phases={
                pid: PhaseState.from_dict(ps_raw)
                for pid, ps_raw in raw.get("phases", {}).items()
            },
        )

    def _path(self) -> Path:
        # Derive filename from host_fingerprint (NOT session_id) so
        # subsequent resume attempts can find this session without
        # passing the session_id around.
        return cache_dir() / f"{SESSION_FILE_PREFIX}{self.host_fingerprint}.json"

    def save(self) -> None:
        """Atomically persist to disk."""
        target = self._path()
        target.parent.mkdir(parents=True, exist_ok=True)
        rendered = json.dumps(self.to_dict(), indent=2, sort_keys=True)
        fd, tmp_name = tempfile.mkstemp(
            prefix=target.name + ".", suffix=".tmp", dir=str(target.parent),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
                f.write(rendered + "\n")
            os.replace(tmp_name, target)
        except Exception:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise

    @classmethod
    def try_resume(cls, host_fingerprint: str) -> "WizardSession | None":
        """Return a saved session if one exists for this fingerprint
        AND is younger than SESSION_TTL. Otherwise None.
        """
        candidate = cache_dir() / f"{SESSION_FILE_PREFIX}{host_fingerprint}.json"
        if not candidate.exists():
            return None
        try:
            raw = json.loads(candidate.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            return None
        if raw.get("host_fingerprint") != host_fingerprint:
            return None  # fingerprint mismatch (safety)
        try:
            session = cls.from_dict(raw)
        except (KeyError, ValueError):
            return None
        age = dt.datetime.now(dt.timezone.utc) - session.started_at
        if age > SESSION_TTL:
            return None
        return session

    def discard(self) -> None:
        """Remove this session's file from disk."""
        target = self._path()
        try:
            target.unlink()
        except FileNotFoundError:
            pass


def list_sessions() -> list[WizardSession]:
    """Return all currently-saved sessions (any host fingerprint).
    Used by `--status`. Sessions that fail to parse are silently
    skipped (they get rediscovered as missing on next try_resume)."""
    out: list[WizardSession] = []
    d = cache_dir()
    if not d.exists():
        return out
    for path in sorted(d.glob(f"{SESSION_FILE_PREFIX}*.json")):
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
            out.append(WizardSession.from_dict(raw))
        except (json.JSONDecodeError, KeyError, ValueError, OSError):
            continue
    return out
