"""Textual TUI for ntl-wizard.

Multi-screen wizard:

  SetupScreen
      │
      ├── (mode=auto)   ── AutoMeasureScreen ──┐
      │                                          ├── ReviewScreen ── (write or discard)
      └── (mode=manual) ── ManualEditScreen ──┘

Every state-mutating callback writes a trace line to
${NTL_WIZARD_CACHE_DIR or ~/.cache/ntl-wizard}/last-tui.log so a
"black screen" / "frozen UI" symptom can be diagnosed by inspecting
that file after the app exits. Long-running synchronous work (compile
+ run timing binaries) is moved to a worker thread via
asyncio.to_thread() so the Textual event loop keeps pumping.
"""
from __future__ import annotations

import argparse
import asyncio
import contextlib
import os
import platform
import shutil
import sys
import traceback
from pathlib import Path
from typing import Optional

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
    PHASES,
    PHASES_BY_ID,
    RuntimeFailure,
    make_ephemeral_build_dir,
    run_phase,
)
from .parameters import PARAMETERS, PARAMETERS_BY_NAME, ValueType
from .search import candidate_sets_for_phase, derive_values
from .session import (
    PhaseState,
    PhaseStatus,
    WizardSession,
    cache_dir,
    compute_host_fingerprint,
)


# ---------------------------------------------------------------------------
# Debug log mirror
# ---------------------------------------------------------------------------


def _debug_log_path() -> Path:
    p = cache_dir() / "last-tui.log"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _trace_init() -> Path:
    p = _debug_log_path()
    p.write_text(
        f"ntl-wizard {WIZARD_VERSION} TUI session started {now_utc_iso()}\n",
        encoding="utf-8",
    )
    return p


def _trace(msg: str) -> None:
    try:
        with _debug_log_path().open("a", encoding="utf-8") as f:
            f.write(f"[{now_utc_iso()}] {msg}\n")
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Shared app state — the four screens read/write this on the App instance.
# ---------------------------------------------------------------------------


def run(args: argparse.Namespace) -> int:
    try:
        from textual.app import App, ComposeResult
        from textual.containers import Vertical, Horizontal
        from textual.screen import Screen
        from textual.widgets import (
            Header, Footer, Static, Button, Checkbox, Input, Log,
            DataTable, RadioSet, RadioButton,
        )
        from textual.binding import Binding
    except ImportError as exc:
        print(
            f"error: TUI mode requires Textual; install or run with --batch: {exc}",
            file=sys.stderr,
        )
        return EXIT_GENERIC

    debug_log = _trace_init()
    _trace(f"args: target={args.target!r} source={args.ntl_source_dir!r} "
           f"output={args.output!r} phases={args.phases!r} "
           f"iterations={args.iterations} dry_run={args.dry_run} "
           f"resume={args.resume}")

    # ---- Setup screen ----------------------------------------------------

    class SetupScreen(Screen):
        BINDINGS = [
            Binding("q", "quit_app", "Quit"),
            Binding("a", "start_auto", "Auto-tune"),
            Binding("m", "start_manual", "Manual edit"),
        ]

        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            with Vertical(id="setup-root"):
                yield Static(
                    f"[b]ntl-wizard[/b]  v{WIZARD_VERSION}\n"
                    f"Logs: {debug_log}",
                    id="banner",
                )
                yield Static("", id="preflight")
                yield Static("[b]Phases to run[/b]  (toggle with Space)", classes="hdr")
                with Vertical(id="phases-box"):
                    for ph in PHASES:
                        yield Checkbox(
                            f"{ph.id}  ({ph.timing_program}, ~{ph.expected_duration_seconds}s)",
                            value=True, id=f"phase-{ph.id}",
                        )
                yield Static("[b]Iterations per measurement[/b]", classes="hdr")
                yield Input(value=str(max(1, self.app.cli_args.iterations)),
                            id="iterations", restrict=r"\d*")
                yield Static("", id="setup-status")
                with Horizontal(id="setup-buttons"):
                    yield Button("Auto-tune (a)",  id="btn-auto",   variant="success")
                    yield Button("Manual edit (m)", id="btn-manual", variant="primary")
                    yield Button("Quit (q)",       id="btn-quit",   variant="error")
            yield Footer()

        async def on_mount(self) -> None:
            _trace("SetupScreen.on_mount")
            # Pre-flight diagnostics — surface state, not gate (gate is in
            # cli._run_tui before reaching here).
            try:
                source_dir = _resolve_source_dir(self.app.cli_args.ntl_source_dir)
                result = plat.check_native(target=self.app.cli_args.target)
                build_dir = source_dir / "build" / "src"
                libntl_present = any(
                    (build_dir / n).exists()
                    for n in ("libntl.so", "libntl.so.0", "libntl.dylib", "libntl.0.dylib")
                )
                lines = [
                    f"source: [b]{source_dir}[/b]",
                    f"host:   {result.build_host_arch} ({'native' if result.kind == plat.CheckResult.OK else 'CROSS'})",
                    f"libntl: {'present' if libntl_present else '[b red]MISSING[/b red] — run `meson compile -C build` first'}",
                ]
                self.query_one("#preflight", Static).update("\n".join(lines))
            except Exception as exc:
                _trace(f"SetupScreen preflight error: {exc}")

        def _selected_phases(self) -> list[str]:
            chosen: list[str] = []
            for ph in PHASES:
                if self.query_one(f"#phase-{ph.id}", Checkbox).value:
                    chosen.append(ph.id)
            return chosen

        def _iterations(self) -> int:
            try:
                return max(1, int(self.query_one("#iterations", Input).value or "1"))
            except ValueError:
                return 1

        def _commit_setup(self) -> bool:
            phases = self._selected_phases()
            if not phases:
                self.query_one("#setup-status", Static).update(
                    "[b red]Pick at least one phase before proceeding[/b red]"
                )
                return False
            self.app.selected_phases = phases
            self.app.iterations = self._iterations()
            _trace(f"setup committed: phases={phases} iterations={self.app.iterations}")
            return True

        async def action_start_auto(self) -> None:
            if not self._commit_setup():
                return
            await self.app.push_screen(AutoMeasureScreen())

        async def action_start_manual(self) -> None:
            if not self._commit_setup():
                return
            await self.app.push_screen(ManualEditScreen())

        async def action_quit_app(self) -> None:
            _trace("user quit from SetupScreen")
            self.app.exit_code = EXIT_OK
            self.app.exit()

        async def on_button_pressed(self, event: Button.Pressed) -> None:
            if event.button.id == "btn-auto":
                await self.action_start_auto()
            elif event.button.id == "btn-manual":
                await self.action_start_manual()
            elif event.button.id == "btn-quit":
                await self.action_quit_app()

    # ---- Auto-measure screen -----------------------------------------------

    class AutoMeasureScreen(Screen):
        BINDINGS = [
            Binding("ctrl+c", "abort", "Abort", priority=True),
            Binding("q", "quit_app", "Quit"),
        ]

        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            with Vertical():
                yield Static("Starting…", id="status")
                yield Log(id="log", highlight=True)
            yield Footer()

        def _log(self, msg: str) -> None:
            _trace(msg)
            try:
                self.query_one("#log", Log).write_line(msg)
            except Exception:
                pass

        async def on_mount(self) -> None:
            _trace("AutoMeasureScreen.on_mount; starting worker")
            self.run_worker(self._main(), exclusive=True, name="auto-measure")

        async def _main(self) -> None:
            try:
                await self._main_inner()
            except Exception:
                tb = traceback.format_exc()
                _trace("UNCAUGHT in AutoMeasureScreen._main:\n" + tb)
                self._log(f"FATAL: {tb.splitlines()[-1]}")
                self._log(f"see {debug_log} for full traceback")
                self.app.exit_code = EXIT_GENERIC

        async def _main_inner(self) -> None:
            status = self.query_one("#status", Static)
            args = self.app.cli_args

            source_dir = _resolve_source_dir(args.ntl_source_dir)
            output_path = _resolve_output_path(args.output, source_dir)
            phase_ids = self.app.selected_phases
            iterations = self.app.iterations
            self._log(f"phases: {phase_ids}; iterations: {iterations}")

            fingerprint = compute_host_fingerprint()
            session: Optional[WizardSession] = None
            if args.resume:
                session = WizardSession.try_resume(host_fingerprint=fingerprint)
            if session is None:
                session = WizardSession.create(
                    host_fingerprint=fingerprint, mode="tui",
                    output_artifact_path=str(output_path),
                )
            session.save()
            self._log(f"session: {session.session_id}")

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

            try:
                derived: dict[str, dict] = {}
                total = len(phase_ids)
                for i, pid in enumerate(phase_ids):
                    phase = PHASES_BY_ID[pid]
                    status.update(f"Phase {i+1}/{total}: {pid}")

                    if (pid in session.phases
                            and session.phases[pid].status == PhaseStatus.COMPLETED):
                        self._log(f"phase {pid}: already complete (resume)")
                        derived[pid] = session.phases[pid].derived_values
                        continue

                    ps = session.phases.setdefault(pid, PhaseState())
                    ps.status = PhaseStatus.RUNNING
                    session.save()

                    candidates = candidate_sets_for_phase(pid)
                    # Phase-budget header so the user has a rough
                    # sense of expected duration.
                    if len(candidates) > 0:
                        per_cand = phase.expected_duration_seconds // max(1, len(candidates))
                        self._log(
                            f"phase {pid}: {len(candidates)} candidate(s) to measure"
                            f"  (expected ~{phase.expected_duration_seconds}s total, "
                            f"~{per_cand}s per candidate)"
                        )
                    else:
                        self._log(f"phase {pid}: 0 candidates (skipped)")

                    # Per-candidate progress callback — bounces each
                    # event back to the Textual event loop via
                    # call_from_thread, so live updates appear in the
                    # Log widget while the measurement work runs in a
                    # worker thread.
                    def _on_progress(idx: int, total_c: int, stage: str, payload) -> None:
                        if stage == "compile":
                            line = f"  [{idx+1}/{total_c}] compiling…  params={payload}"
                        elif stage == "run":
                            line = f"  [{idx+1}/{total_c}] running…    params={payload}"
                        elif stage == "tick":
                            sub_stage, elapsed = payload
                            line = (
                                f"  [{idx+1}/{total_c}] {sub_stage} still running "
                                f"({elapsed:.0f}s elapsed)…"
                            )
                        elif stage == "line":
                            # Live streaming of subprocess stdout/stderr
                            # (gcc/clang output, our instrumented
                            # `[Poly1 pass X/5]` lines, etc.). Indent
                            # so it's visually distinct from the
                            # wrapper logs.
                            sub_stage, raw = payload
                            line = f"    │ {raw}"
                        elif stage == "done":
                            line = (
                                f"  [{idx+1}/{total_c}] done in "
                                f"{payload.wall_clock_seconds:.3f}s  "
                                f"(stddev {payload.noise_estimate:.3f}s)"
                            )
                        else:
                            line = f"  [{idx+1}/{total_c}] {stage}"
                        self.app.call_from_thread(self._log, line)

                    try:
                        # asyncio.to_thread keeps Textual's event loop
                        # responsive while the synchronous compile +
                        # run subprocess work runs in a worker thread.
                        measurements = await asyncio.to_thread(
                            run_phase,
                            context, phase, candidates,
                            repeats=iterations,
                            progress_callback=_on_progress,
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
                        self.app.exit_code = kind_to_code[kind]
                        status.update(
                            f"Phase {pid} failed ({kind}). Press Q to quit "
                            f"(exit {self.app.exit_code})."
                        )
                        return

                    chosen = derive_values(pid, measurements)
                    ps.derived_values = chosen
                    ps.status = PhaseStatus.COMPLETED
                    session.save()
                    derived[pid] = chosen
                    self._log(f"phase {pid} ok: {chosen}")

                # Hand off to ReviewScreen with the derived values.
                all_values: dict[str, int | bool | str] = {
                    p.name: p.default_value for p in PARAMETERS
                }
                for chosen in derived.values():
                    all_values.update(chosen)
                self.app.candidate_values = all_values
                self.app.session = session
                status.update("Measurement complete — switching to review…")
                self._log("measurement done; opening ReviewScreen")
                await self.app.push_screen(ReviewScreen())
            finally:
                with contextlib.suppress(OSError):
                    shutil.rmtree(sub_build, ignore_errors=True)

        async def action_abort(self) -> None:
            _trace("user aborted AutoMeasureScreen")
            self.app.exit_code = 130
            self.app.exit()

        async def action_quit_app(self) -> None:
            _trace("user quit from AutoMeasureScreen")
            self.app.exit()

    # ---- Manual edit screen -----------------------------------------------

    class ManualEditScreen(Screen):
        BINDINGS = [
            Binding("q",      "quit_app", "Quit"),
            Binding("ctrl+s", "save",     "Save"),
            Binding("escape", "back",     "Back"),
        ]

        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            with Vertical():
                yield Static(
                    "[b]Manual parameter editor[/b]  —  Tab into the table, "
                    "click a value to cycle through the domain. Ctrl-S to save.",
                    id="manual-banner",
                )
                yield DataTable(id="params-table", zebra_stripes=True)
                yield Static("", id="manual-status")
                with Horizontal():
                    yield Button("Save (Ctrl-S)", id="btn-save", variant="success")
                    yield Button("Back (Esc)",    id="btn-back", variant="default")
                    yield Button("Quit (q)",      id="btn-quit", variant="error")
            yield Footer()

        async def on_mount(self) -> None:
            _trace("ManualEditScreen.on_mount")
            table = self.query_one("#params-table", DataTable)
            table.add_columns("Parameter", "Family", "Value", "Domain")
            # Seed from current candidate_values if AutoMeasure ran
            # before, else from defaults.
            seed = dict(getattr(self.app, "candidate_values", None) or {})
            self.app.candidate_values = {
                p.name: seed.get(p.name, p.default_value) for p in PARAMETERS
            }
            for p in PARAMETERS:
                value = self.app.candidate_values[p.name]
                table.add_row(
                    p.name, p.family, str(value),
                    ",".join(str(v) for v in p.value_domain),
                    key=p.name,
                )

        async def on_data_table_cell_selected(self, event) -> None:  # type: ignore[no-untyped-def]
            row_key = event.cell_key.row_key.value if event.cell_key.row_key else None
            if not row_key:
                return
            param = PARAMETERS_BY_NAME[row_key]
            domain = list(param.value_domain)
            current = self.app.candidate_values[param.name]
            try:
                idx = domain.index(current)
            except ValueError:
                idx = -1
            new_value = domain[(idx + 1) % len(domain)]
            self.app.candidate_values[param.name] = new_value
            _trace(f"manual: {param.name} {current} -> {new_value}")
            table = self.query_one("#params-table", DataTable)
            table.update_cell(row_key, "Value", str(new_value))

        async def action_save(self) -> None:
            self.app.session = None  # manual mode → no measurements
            await self.app.push_screen(ReviewScreen())

        async def action_back(self) -> None:
            self.app.pop_screen()

        async def action_quit_app(self) -> None:
            self.app.exit()

        async def on_button_pressed(self, event: Button.Pressed) -> None:
            if event.button.id == "btn-save":
                await self.action_save()
            elif event.button.id == "btn-back":
                await self.action_back()
            elif event.button.id == "btn-quit":
                await self.action_quit_app()

    # ---- Review / Confirm screen ------------------------------------------

    class ReviewScreen(Screen):
        BINDINGS = [
            Binding("w",      "write",   "Write artifact"),
            Binding("e",      "edit",    "Edit manually"),
            Binding("d",      "discard", "Discard"),
            Binding("q",      "quit_app", "Quit"),
        ]

        def compose(self) -> ComposeResult:
            yield Header(show_clock=False)
            with Vertical():
                yield Static("[b]Review tune values[/b]  —  W to write, E to edit, D to discard.",
                             id="review-banner")
                yield DataTable(id="review-table", zebra_stripes=True)
                yield Static("", id="review-status")
                with Horizontal():
                    yield Button("Write (w)",   id="btn-write",   variant="success")
                    yield Button("Edit (e)",    id="btn-edit",    variant="primary")
                    yield Button("Discard (d)", id="btn-discard", variant="error")
            yield Footer()

        async def on_mount(self) -> None:
            _trace("ReviewScreen.on_mount")
            table = self.query_one("#review-table", DataTable)
            table.add_columns("Parameter", "Family", "Chosen value")
            values = self.app.candidate_values or {}
            for p in PARAMETERS:
                table.add_row(p.name, p.family, str(values.get(p.name, p.default_value)))

        async def action_write(self) -> None:
            args = self.app.cli_args
            source_dir = _resolve_source_dir(args.ntl_source_dir)
            output_path = _resolve_output_path(args.output, source_dir)
            try:
                session_id = self.app.session.session_id if self.app.session else "manual"
                provenance = {
                    "ntl_version": (source_dir / "version.txt").read_text().strip(),
                    "wizard_version": WIZARD_VERSION,
                    "host_fingerprint": compute_host_fingerprint(),
                    "host_cpu": platform.processor() or platform.machine(),
                    "host_os": f"{platform.system()} {platform.release()}",
                    "compiler": (os.environ.get("CXX") or shutil.which("c++")
                                 or shutil.which("g++") or "c++"),
                    "generated_utc": now_utc_iso(),
                    "session_id": session_id,
                }
                write_artifact(output_path, self.app.candidate_values, provenance)
                _trace(f"wrote artifact at {output_path}")
                # Hand the paths off to the CLI wrapper so it prints
                # the next-steps message on the REAL terminal after
                # the TUI exits — anything we put in the Log widget
                # vanishes with the alternate-screen when the user
                # presses Q.
                self.app.wrote_artifact_path = Path(output_path)
                self.app.wrote_source_dir = Path(source_dir)
                self.query_one("#review-status", Static).update(
                    f"[b green]✓ Wrote tune artifact:[/b green] {output_path}\n"
                    f"Press Q to quit — next-steps will print on the terminal."
                )
                self.app.exit_code = EXIT_OK
            except Exception as exc:
                _trace(f"write_artifact failed: {exc}")
                self.query_one("#review-status", Static).update(
                    f"[b red]Write failed: {exc}[/b red]\nPress E to edit or D to discard."
                )

        async def action_edit(self) -> None:
            await self.app.push_screen(ManualEditScreen())

        async def action_discard(self) -> None:
            _trace("user discarded review (no write)")
            self.query_one("#review-status", Static).update(
                "Discarded. Press Q to quit (nothing written)."
            )
            self.app.exit_code = EXIT_OK

        async def action_quit_app(self) -> None:
            self.app.exit()

        async def on_button_pressed(self, event: Button.Pressed) -> None:
            if event.button.id == "btn-write":
                await self.action_write()
            elif event.button.id == "btn-edit":
                await self.action_edit()
            elif event.button.id == "btn-discard":
                await self.action_discard()

    # ---- Main App ---------------------------------------------------------

    class WizardApp(App):
        CSS = """
        Screen { layout: vertical; }
        #banner, #manual-banner, #review-banner { padding: 1 2; background: $boost; }
        #preflight { padding: 1 2; }
        .hdr { padding: 0 2; color: $accent; }
        #phases-box { padding: 0 2; height: auto; }
        #setup-buttons, Horizontal { height: auto; padding: 1 2; }
        Button { margin: 0 1; }
        Input { width: 12; margin: 0 2; }
        #log { height: 1fr; margin: 0 1 1 1; border: round $primary; }
        DataTable { height: 1fr; margin: 0 1; }
        """
        BINDINGS = [Binding("ctrl+c", "force_quit", "Force quit", priority=True)]
        TITLE = f"ntl-wizard {WIZARD_VERSION}"

        def __init__(self, cli_args: argparse.Namespace) -> None:
            super().__init__()
            self.cli_args = cli_args
            self.exit_code = EXIT_OK
            # Shared mutable state passed between screens
            self.selected_phases: list[str] = [p.id for p in PHASES]
            self.iterations: int = max(1, cli_args.iterations)
            self.candidate_values: Optional[dict] = None
            self.session: Optional[WizardSession] = None
            # Set by ReviewScreen.action_write on success. Read by the
            # CLI wrapper AFTER app.run() returns so the next-steps
            # message is printed on the real terminal (and therefore
            # remains in the user's shell scrollback) rather than
            # vanishing with the TUI alternate-screen.
            self.wrote_artifact_path: Optional[Path] = None
            self.wrote_source_dir: Optional[Path] = None

        def on_mount(self) -> None:
            _trace("WizardApp.on_mount; pushing SetupScreen")
            self.push_screen(SetupScreen())

        def action_force_quit(self) -> None:
            _trace("Ctrl-C force quit")
            self.exit_code = 130
            self.exit()

    app = WizardApp(args)
    app.run()

    # Print next-steps on the real terminal AFTER Textual has torn
    # down the alternate screen. This way the message persists in the
    # user's shell scrollback. Skipped if the user didn't write an
    # artifact (e.g. quit from setup, discarded review).
    if app.wrote_artifact_path is not None and app.exit_code == EXIT_OK:
        rel_build = "build"
        # ANSI: bold green for headings, cyan for commands. Falls back
        # to plain text on terminals that don't honor it.
        BOLD = "\033[1m"
        GREEN = "\033[32m"
        CYAN  = "\033[36m"
        DIM   = "\033[2m"
        RESET = "\033[0m"
        print()
        print(f"{BOLD}{GREEN}✓ ntl-wizard wrote the tune artifact:{RESET}")
        print(f"  {app.wrote_artifact_path}")
        print()
        print(f"{BOLD}Next steps{RESET} {DIM}(from {app.wrote_source_dir}){RESET}:")
        print()
        print(f"  1. Rebuild NTL with the host-tuned table:")
        print(f"       {CYAN}meson setup --buildtype=release "
              f"-Dtune=host {rel_build}{RESET}")
        print(f"       {CYAN}meson compile -C {rel_build}{RESET}")
        print()
        print(f"  2. (optional) Verify with the test suite:")
        print(f"       {CYAN}meson test -C {rel_build}{RESET}")
        print()
        print(f"  3. Install where you want it:")
        print(f"       {CYAN}meson install -C {rel_build} "
              f"--destdir <target>{RESET}")
        print()
        print(f"{DIM}The host-tuned.ini lives under src/meson/tune-tables/ "
              f"and is gitignored by default.{RESET}")
        print(f"{DIM}If you want this tuning to be reproducible across a "
              f"team, commit it explicitly.{RESET}")
        print()
    return app.exit_code
