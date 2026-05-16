"""Textual TUI for ntl-wizard.

A minimal but functional terminal interface that runs the same
measurement orchestration as `cli._run_batch` but with live progress
feedback. Designed to degrade to batch mode if Textual cannot create
a terminal app (e.g. unusual TTY setups).
"""
from __future__ import annotations

import argparse
import contextlib
import os
import platform
import shutil
import sys
from pathlib import Path

from . import __version__ as WIZARD_VERSION
from . import platform_check as plat
from .artifacts import write_artifact, now_utc_iso
from .cli import (
    EXIT_OK,
    EXIT_GENERIC,
    EXIT_CROSS_REFUSAL,
    EXIT_COMPILE_FAILURE,
    EXIT_RUNTIME_FAILURE,
    EXIT_MEASUREMENT_NOISE,
    _resolve_output_path,
    _resolve_phase_ids,
    _resolve_source_dir,
    _err,
)
from .measure import (
    CompileFailure,
    MeasureContext,
    MeasurementNoiseTooHigh,
    PHASES,
    PHASES_BY_ID,
    RuntimeFailure,
    make_ephemeral_build_dir,
    run_phase,
)
from .parameters import PARAMETERS
from .search import candidate_sets_for_phase, derive_values
from .session import (
    PhaseState,
    PhaseStatus,
    WizardSession,
    compute_host_fingerprint,
)


def run(args: argparse.Namespace) -> int:
    """Entry point invoked from `cli._run_tui`."""
    try:
        from textual.app import App, ComposeResult
        from textual.widgets import Header, Footer, Static, ProgressBar, Log
        from textual.containers import Vertical
        from textual.binding import Binding
    except ImportError as exc:
        _err(
            f"Textual not available ({exc}); fall back with --batch or "
            f"`pip install textual`."
        )
        return EXIT_GENERIC

    class WizardApp(App):
        CSS = """
        Screen {
            background: $surface;
        }
        #phases-list {
            border: solid $accent;
            padding: 1 2;
            height: auto;
            margin: 1 2;
        }
        #log {
            border: solid $primary;
            height: 1fr;
            margin: 0 2 1 2;
        }
        ProgressBar {
            margin: 1 2;
        }
        """
        BINDINGS = [
            Binding("ctrl+c", "abort", "Abort", priority=True),
            Binding("q", "quit", "Quit"),
        ]
        TITLE = f"ntl-wizard {WIZARD_VERSION}"

        def __init__(self, args: argparse.Namespace) -> None:
            super().__init__()
            self.args = args
            self.exit_code = EXIT_OK
            self.phase_ids: list[str] = []
            self.source_dir: Path | None = None
            self.output_path: Path | None = None
            self.session: WizardSession | None = None

        def compose(self) -> ComposeResult:
            yield Header()
            with Vertical():
                yield Static("Initializing…", id="status")
                yield Static("Phases:", id="phases-list")
                yield ProgressBar(id="progress", total=100, show_eta=False)
                yield Log(id="log", highlight=True)
            yield Footer()

        async def on_mount(self) -> None:
            self.run_worker(self._main(), exclusive=True)

        async def _main(self) -> None:
            log = self.query_one("#log", Log)
            status = self.query_one("#status", Static)
            progress = self.query_one("#progress", ProgressBar)

            # Platform check
            result = plat.check_native(target=self.args.target)
            if result.kind == plat.CheckResult.CROSS_REFUSAL:
                log.write_line(f"REFUSED: {result.message.splitlines()[0]}")
                _err(result.message)
                self.exit_code = EXIT_CROSS_REFUSAL
                self.exit()
                return
            log.write_line(f"platform: {result.build_host_arch} (native, OK)")

            # Source / output
            try:
                self.source_dir = _resolve_source_dir(self.args.ntl_source_dir)
                if not (self.source_dir / "src").exists():
                    raise FileNotFoundError(f"src/ not under {self.source_dir}")
                self.output_path = _resolve_output_path(self.args.output, self.source_dir)
                self.phase_ids = _resolve_phase_ids(self.args.phases)
            except (FileNotFoundError, ValueError) as exc:
                log.write_line(f"FATAL: {exc}")
                self.exit_code = EXIT_GENERIC
                self.exit()
                return

            log.write_line(f"ntl-source-dir: {self.source_dir}")
            log.write_line(f"output: {self.output_path}")
            log.write_line(f"phases: {self.phase_ids}")

            # Phases list display
            phases_widget = self.query_one("#phases-list", Static)
            phases_widget.update(
                "Phases:\n" + "\n".join(
                    f"  [{i+1}/{len(self.phase_ids)}] {pid}"
                    for i, pid in enumerate(self.phase_ids)
                )
            )

            if self.args.dry_run:
                log.write_line("dry-run: no measurements taken")
                self.exit_code = EXIT_OK
                self.exit()
                return

            # Session
            fingerprint = compute_host_fingerprint()
            if self.args.resume:
                self.session = WizardSession.try_resume(host_fingerprint=fingerprint)
            if self.session is None:
                self.session = WizardSession.create(
                    host_fingerprint=fingerprint, mode="tui",
                    output_artifact_path=str(self.output_path),
                )
            self.session.save()

            # Measure
            sub_build = make_ephemeral_build_dir()
            context = MeasureContext(
                ntl_source_dir=self.source_dir,
                sub_build_dir=sub_build,
                compiler=(
                    os.environ.get("CXX")
                    or shutil.which("c++")
                    or shutil.which("g++")
                    or "c++"
                ),
            )

            try:
                derived: dict[str, dict] = {}
                total_phases = len(self.phase_ids)
                for i, pid in enumerate(self.phase_ids):
                    phase = PHASES_BY_ID[pid]
                    status.update(f"Phase {i+1}/{total_phases}: {pid}")
                    progress.update(progress=int((i / total_phases) * 100))

                    if (pid in self.session.phases
                            and self.session.phases[pid].status == PhaseStatus.COMPLETED):
                        log.write_line(f"phase {pid}: already completed (resume)")
                        derived[pid] = self.session.phases[pid].derived_values
                        continue

                    log.write_line(f"phase {pid} starting")
                    ps = self.session.phases.setdefault(pid, PhaseState())
                    ps.status = PhaseStatus.RUNNING
                    self.session.save()

                    candidates = candidate_sets_for_phase(pid)
                    try:
                        measurements = run_phase(
                            context, phase, candidates,
                            repeats=max(1, self.args.iterations),
                        )
                    except (CompileFailure, RuntimeFailure, MeasurementNoiseTooHigh) as exc:
                        # Stay open so the user can read the error.
                        # Setting exit_code now means `q` (or the
                        # Footer's quit binding) returns the right
                        # status when the user dismisses the screen.
                        kind = type(exc).__name__
                        kind_to_code = {
                            "CompileFailure":          EXIT_COMPILE_FAILURE,
                            "RuntimeFailure":          EXIT_RUNTIME_FAILURE,
                            "MeasurementNoiseTooHigh": EXIT_MEASUREMENT_NOISE,
                        }
                        log.write_line(f"{kind.upper()}: {exc}")
                        ps.status = PhaseStatus.FAILED
                        ps.error = str(exc)
                        self.session.save()
                        self.exit_code = kind_to_code[kind]
                        status.update(
                            f"Phase {pid} failed ({kind}). Press Q to quit "
                            f"(exit code {self.exit_code})."
                        )
                        return  # leave the app running; user dismisses

                    chosen = derive_values(pid, measurements)
                    ps.derived_values = chosen
                    ps.status = PhaseStatus.COMPLETED
                    self.session.save()
                    derived[pid] = chosen
                    log.write_line(f"phase {pid}: {chosen}")

                # Write artifact
                progress.update(progress=100)
                status.update("Writing tune artifact…")
                all_values: dict[str, int | bool | str] = {p.name: p.default_value for p in PARAMETERS}
                for chosen in derived.values():
                    all_values.update(chosen)
                provenance = {
                    "ntl_version": (self.source_dir / "version.txt").read_text().strip(),
                    "wizard_version": WIZARD_VERSION,
                    "host_fingerprint": fingerprint,
                    "host_cpu": platform.processor() or platform.machine(),
                    "host_os": f"{platform.system()} {platform.release()}",
                    "compiler": context.compiler,
                    "generated_utc": now_utc_iso(),
                    "session_id": self.session.session_id,
                }
                write_artifact(self.output_path, all_values, provenance)
                log.write_line(f"OK: artifact={self.output_path}")
                status.update("Complete. Press Q to quit.")
            finally:
                with contextlib.suppress(OSError):
                    shutil.rmtree(sub_build, ignore_errors=True)

        def action_abort(self) -> None:
            self.exit_code = 130  # SIGINT exit code per contract
            self.exit()

    app = WizardApp(args)
    app.run()
    return app.exit_code
