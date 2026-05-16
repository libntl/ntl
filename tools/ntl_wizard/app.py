"""Textual TUI for ntl-wizard.

A minimal but functional terminal interface that runs the same
measurement orchestration as `cli._run_batch` but with live progress
feedback.

Debugging note: every line emitted by the TUI worker is mirrored to
${NTL_WIZARD_CACHE_DIR or ~/.cache/ntl-wizard}/last-tui.log so a
"black screen" symptom can be diagnosed by inspecting that file after
the app exits. Unhandled worker exceptions are also captured to the
same file with a traceback.
"""
from __future__ import annotations

import argparse
import asyncio
import contextlib
import datetime as dt
import os
import platform
import shutil
import sys
import traceback
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
)
from .measure import (
    CompileFailure,
    MeasureContext,
    MeasurementNoiseTooHigh,
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
    cache_dir,
    compute_host_fingerprint,
)


def _debug_log_path() -> Path:
    p = cache_dir() / "last-tui.log"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def run(args: argparse.Namespace) -> int:
    try:
        from textual.app import App, ComposeResult
        from textual.widgets import Header, Footer, Static, Log
        from textual.containers import Vertical
        from textual.binding import Binding
    except ImportError as exc:
        print(
            f"error: TUI mode requires Textual; install or run with --batch: {exc}",
            file=sys.stderr,
        )
        return EXIT_GENERIC

    # Wipe / start the debug log up front so users always have a fresh
    # record of the last TUI session, even if Textual exits with a
    # "black screen" before any compose() output.
    debug_log = _debug_log_path()
    debug_log.write_text(
        f"ntl-wizard {WIZARD_VERSION} TUI session started {now_utc_iso()}\n",
        encoding="utf-8",
    )

    def _trace(msg: str) -> None:
        with debug_log.open("a", encoding="utf-8") as f:
            f.write(f"[{now_utc_iso()}] {msg}\n")

    class WizardApp(App):
        CSS = """
        Screen { layout: vertical; }
        #status { padding: 1 2; background: $boost; color: $text; height: auto; }
        #log    { height: 1fr; margin: 0 1 1 1; border: round $primary; }
        """
        BINDINGS = [
            Binding("ctrl+c", "abort", "Abort", priority=True),
            Binding("q",      "quit",  "Quit"),
        ]
        TITLE = f"ntl-wizard {WIZARD_VERSION}"

        def __init__(self, args: argparse.Namespace) -> None:
            super().__init__()
            self.args = args
            self.exit_code = EXIT_OK

        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            with Vertical():
                yield Static(
                    f"Starting ntl-wizard {WIZARD_VERSION}…\n"
                    f"Logs mirror to: {debug_log}",
                    id="status",
                )
                yield Log(id="log", highlight=True)
            yield Footer()

        def _log(self, message: str) -> None:
            """Mirror a line to both the TUI Log widget and the debug
            file. Safe to call even if the widget isn't ready yet."""
            _trace(message)
            try:
                self.query_one("#log", Log).write_line(message)
            except Exception:  # widget not ready yet
                pass

        async def on_mount(self) -> None:
            _trace("on_mount fired; starting worker")
            self.run_worker(self._main(), exclusive=True, name="wizard-main")

        async def _main(self) -> None:
            try:
                await self._main_inner()
            except Exception:
                tb = traceback.format_exc()
                _trace("UNCAUGHT in _main_inner:\n" + tb)
                self._log(f"FATAL (uncaught): {tb.splitlines()[-1]}")
                self._log(f"see {debug_log} for full traceback")
                self._log("press Q to quit")
                self.exit_code = EXIT_GENERIC

        async def _main_inner(self) -> None:
            status = self.query_one("#status", Static)

            # Platform check
            self._log("checking platform…")
            result = plat.check_native(target=self.args.target)
            if result.kind == plat.CheckResult.CROSS_REFUSAL:
                self._log(f"REFUSED: {result.message.splitlines()[0]}")
                status.update("Cross-compile context — refused. Press Q to quit.")
                self.exit_code = EXIT_CROSS_REFUSAL
                return
            self._log(f"platform OK: {result.build_host_arch} (native)")

            # Source / output / phases
            try:
                source_dir = _resolve_source_dir(self.args.ntl_source_dir)
                if not (source_dir / "src").exists():
                    raise FileNotFoundError(f"src/ not under {source_dir}")
                output_path = _resolve_output_path(self.args.output, source_dir)
                phase_ids = _resolve_phase_ids(self.args.phases)
            except (FileNotFoundError, ValueError) as exc:
                self._log(f"FATAL setup: {exc}")
                status.update(f"Setup error: {exc}. Press Q to quit.")
                self.exit_code = EXIT_GENERIC
                return

            self._log(f"source dir: {source_dir}")
            self._log(f"output:     {output_path}")
            self._log(f"phases:     {', '.join(phase_ids)}")

            if self.args.dry_run:
                self._log("dry-run: no measurements taken")
                status.update("Dry-run complete. Press Q to quit.")
                return

            # Session
            fingerprint = compute_host_fingerprint()
            session: WizardSession | None = None
            if self.args.resume:
                session = WizardSession.try_resume(host_fingerprint=fingerprint)
            if session is None:
                session = WizardSession.create(
                    host_fingerprint=fingerprint, mode="tui",
                    output_artifact_path=str(output_path),
                )
            session.save()
            self._log(f"session: {session.session_id} ({session.mode})")

            # Compile + measurement loop
            sub_build = make_ephemeral_build_dir()
            context = MeasureContext(
                ntl_source_dir=source_dir,
                sub_build_dir=sub_build,
                compiler=(
                    os.environ.get("CXX")
                    or shutil.which("c++")
                    or shutil.which("g++")
                    or "c++"
                ),
            )
            self._log(f"compiler: {context.compiler}")
            self._log(f"sub-build dir: {sub_build}")

            try:
                derived: dict[str, dict] = {}
                total = len(phase_ids)
                for i, pid in enumerate(phase_ids):
                    phase = PHASES_BY_ID[pid]
                    status.update(f"Phase {i+1}/{total}: {pid}")

                    if (pid in session.phases
                            and session.phases[pid].status == PhaseStatus.COMPLETED):
                        self._log(f"phase {pid}: resume — already completed")
                        derived[pid] = session.phases[pid].derived_values
                        continue

                    self._log(f"phase {pid} starting")
                    ps = session.phases.setdefault(pid, PhaseState())
                    ps.status = PhaseStatus.RUNNING
                    session.save()

                    candidates = candidate_sets_for_phase(pid)
                    self._log(f"  {len(candidates)} candidate parameter set(s) to measure")

                    try:
                        # CRITICAL: run_phase is synchronous (subprocess.run
                        # to compile and run timing binaries — can take
                        # minutes per candidate). Running it directly in
                        # this coroutine would block the Textual event
                        # loop, freezing the UI and ignoring Ctrl-C / Q.
                        # asyncio.to_thread() moves it to a worker thread
                        # so the event loop keeps pumping.
                        measurements = await asyncio.to_thread(
                            run_phase,
                            context, phase, candidates,
                            repeats=max(1, self.args.iterations),
                        )
                    except (CompileFailure, RuntimeFailure, MeasurementNoiseTooHigh) as exc:
                        kind = type(exc).__name__
                        kind_to_code = {
                            "CompileFailure":          EXIT_COMPILE_FAILURE,
                            "RuntimeFailure":          EXIT_RUNTIME_FAILURE,
                            "MeasurementNoiseTooHigh": EXIT_MEASUREMENT_NOISE,
                        }
                        self._log(f"{kind.upper()}: {exc}")
                        ps.status = PhaseStatus.FAILED
                        ps.error = str(exc)
                        session.save()
                        self.exit_code = kind_to_code[kind]
                        status.update(
                            f"Phase {pid} failed ({kind}). Press Q to quit "
                            f"(exit {self.exit_code})."
                        )
                        return

                    chosen = derive_values(pid, measurements)
                    ps.derived_values = chosen
                    ps.status = PhaseStatus.COMPLETED
                    session.save()
                    derived[pid] = chosen
                    self._log(f"phase {pid} ok: {chosen}")

                # Write artifact
                status.update("Writing tune artifact…")
                all_values: dict[str, int | bool | str] = {p.name: p.default_value for p in PARAMETERS}
                for chosen in derived.values():
                    all_values.update(chosen)
                provenance = {
                    "ntl_version": (source_dir / "version.txt").read_text().strip(),
                    "wizard_version": WIZARD_VERSION,
                    "host_fingerprint": fingerprint,
                    "host_cpu": platform.processor() or platform.machine(),
                    "host_os": f"{platform.system()} {platform.release()}",
                    "compiler": context.compiler,
                    "generated_utc": now_utc_iso(),
                    "session_id": session.session_id,
                }
                write_artifact(output_path, all_values, provenance)
                self._log(f"OK: artifact={output_path}")
                status.update("Complete. Press Q to quit.")

            finally:
                with contextlib.suppress(OSError):
                    shutil.rmtree(sub_build, ignore_errors=True)

        def action_abort(self) -> None:
            self.exit_code = 130
            self.exit()

    app = WizardApp(args)
    app.run()
    return app.exit_code
