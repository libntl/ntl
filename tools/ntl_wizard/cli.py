"""ntl-wizard CLI surface (contracts/ntl-wizard-cli.md).

Typer-driven CLI that dispatches to either the TUI (default) or the
batch driver (`--batch`). Exit codes per contract. The contract is
CLI-library agnostic — tests in tests/ntl_wizard/test_cli_contract.py
invoke as subprocess and verify behavior only.
"""
from __future__ import annotations

import contextlib
import datetime as dt
import json
import os
import platform
import shutil
import signal
import sys
from pathlib import Path
from typing import NoReturn, Optional

import typer

from . import __version__ as WIZARD_VERSION
from . import platform_check as plat
from .artifacts import write_artifact, now_utc_iso
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
    list_sessions,
)


# ----- Exit codes (per contracts/ntl-wizard-cli.md) -----
EXIT_OK = 0
EXIT_GENERIC = 1
EXIT_CROSS_REFUSAL = 2
EXIT_COMPILE_FAILURE = 3
EXIT_RUNTIME_FAILURE = 4
EXIT_MEASUREMENT_NOISE = 5
EXIT_SESSION_CONFLICT = 6
EXIT_USER_INTERRUPT = 130


# ---------------------------------------------------------------------------
# Helpers (unchanged between argparse / typer revs)
# ---------------------------------------------------------------------------


def _resolve_source_dir(arg: Optional[Path]) -> Path:
    if arg is not None:
        return arg.resolve()
    if (env := os.environ.get("NTL_WIZARD_SOURCE_DIR")):
        return Path(env).resolve()
    return Path.cwd().resolve()


def _resolve_output_path(arg: Optional[Path], source_dir: Path) -> Path:
    if arg is not None:
        return arg.resolve()
    if (env := os.environ.get("NTL_WIZARD_OUTPUT")):
        return Path(env).resolve()
    return (source_dir / "src" / "meson" / "tune-tables" / "host-tuned.ini").resolve()


def _resolve_phase_ids(raw: Optional[str]) -> list[str]:
    if not raw:
        return [p.id for p in PHASES]
    requested = [s.strip() for s in raw.split(",") if s.strip()]
    unknown = [r for r in requested if r not in PHASES_BY_ID]
    if unknown:
        raise ValueError(f"Unknown phase id(s): {unknown}; known: {list(PHASES_BY_ID)}")
    return requested


def _log(level: str, message: str, *, stream=sys.stdout) -> None:
    print(f"[{now_utc_iso()}] {level} {message}", file=stream, flush=True)


def _err(message: str) -> None:
    """First-stderr-line `error: ...` convention per CLI contract."""
    print(f"error: {message}", file=sys.stderr, flush=True)


def _version_callback(value: bool) -> None:
    if value:
        print(f"ntl-wizard {WIZARD_VERSION}")
        raise typer.Exit()


def _check_tty_or_recommend_batch() -> Optional[int]:
    """Refuse if no TTY and --batch not requested."""
    if not sys.stdout.isatty() and not sys.stderr.isatty():
        _err(
            "No TTY detected and --batch was not requested. "
            "Re-run as `ntl-wizard --batch` for non-interactive use."
        )
        return EXIT_GENERIC
    return None


# ---------------------------------------------------------------------------
# Dispatchers
# ---------------------------------------------------------------------------


def _do_status() -> int:
    sessions = list_sessions()
    payload = {
        "sessions": [
            {
                "session_id": s.session_id,
                "host_fingerprint": s.host_fingerprint,
                "started_at": s.started_at.isoformat(),
                "mode": s.mode,
                "phases": {
                    pid: {"status": ps.status.value} for pid, ps in s.phases.items()
                },
            }
            for s in sessions
        ],
    }
    json.dump(payload, sys.stdout, indent=2)
    print()
    return EXIT_OK


def _run_batch(
    target: Optional[str],
    ntl_source_dir: Optional[Path],
    output: Optional[Path],
    phases_arg: Optional[str],
    iterations: int,
    dry_run: bool,
    resume: bool,
    quiet: bool,
) -> int:
    """The batch / non-interactive driver."""
    plat_result = plat.check_native(target=target)
    if plat_result.kind == plat.CheckResult.CROSS_REFUSAL:
        _err(plat_result.message)
        return EXIT_CROSS_REFUSAL
    if not quiet:
        _log("INFO", f"platform: {plat_result.build_host_arch} (native, OK)")

    source_dir = _resolve_source_dir(ntl_source_dir)
    if not (source_dir / "src").exists():
        _err(f"NTL source dir not found or missing src/: {source_dir}")
        return EXIT_GENERIC
    if not quiet:
        _log("INFO", f"ntl-source-dir: {source_dir}")

    try:
        phase_ids = _resolve_phase_ids(phases_arg)
    except ValueError as exc:
        _err(str(exc))
        return EXIT_GENERIC

    output_path = _resolve_output_path(output, source_dir)
    if not quiet:
        _log("INFO", f"artifact target: {output_path}")

    if dry_run:
        for pid in phase_ids:
            phase = PHASES_BY_ID[pid]
            src = source_dir / "src" / phase.timing_program
            exists = "OK" if src.exists() else "MISSING"
            _log("INFO", f"phase {pid}: {phase.timing_program} {exists}")
        _log("INFO", "dry-run: no measurements taken")
        return EXIT_OK

    fingerprint = compute_host_fingerprint()
    session: Optional[WizardSession] = None
    if resume:
        session = WizardSession.try_resume(host_fingerprint=fingerprint)
        if session is None:
            _err("No resumable session found for this host.")
            return EXIT_SESSION_CONFLICT
    if session is None:
        session = WizardSession.create(
            host_fingerprint=fingerprint, mode="batch",
            output_artifact_path=str(output_path),
        )
    session.save()

    sub_build = make_ephemeral_build_dir()
    context = MeasureContext(
        ntl_source_dir=source_dir,
        sub_build_dir=sub_build,
        compiler=os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++") or "c++",
    )

    try:
        derived: dict[str, dict] = {}
        for pid in phase_ids:
            phase = PHASES_BY_ID[pid]
            if pid in session.phases and session.phases[pid].status == PhaseStatus.COMPLETED:
                _log("INFO", f"phase {pid}: already completed (resume)")
                derived[pid] = session.phases[pid].derived_values
                continue
            _log("INFO", f"phase {pid} starting")
            ps = session.phases.setdefault(pid, PhaseState())
            ps.status = PhaseStatus.RUNNING
            ps.started_at = dt.datetime.now(dt.timezone.utc)
            session.save()
            candidates = candidate_sets_for_phase(pid)
            try:
                measurements = run_phase(
                    context, phase, candidates, repeats=max(1, iterations),
                )
            except CompileFailure as exc:
                _err(f"compile failure in phase {pid}: {exc}")
                ps.status = PhaseStatus.FAILED
                ps.error = str(exc)
                session.save()
                return EXIT_COMPILE_FAILURE
            except RuntimeFailure as exc:
                _err(f"runtime failure in phase {pid}: {exc}")
                ps.status = PhaseStatus.FAILED
                ps.error = str(exc)
                session.save()
                return EXIT_RUNTIME_FAILURE
            except MeasurementNoiseTooHigh as exc:
                _err(f"measurement noise too high in phase {pid}: {exc}")
                ps.status = PhaseStatus.FAILED
                ps.error = str(exc)
                session.save()
                return EXIT_MEASUREMENT_NOISE
            chosen = derive_values(pid, measurements)
            ps.derived_values = chosen
            ps.measurements = [
                {
                    "parameter_set": m.parameter_set,
                    "wall_clock_seconds": m.wall_clock_seconds,
                    "iteration_count": m.iteration_count,
                } for m in measurements
            ]
            ps.completed_at = dt.datetime.now(dt.timezone.utc)
            ps.status = PhaseStatus.COMPLETED
            session.save()
            derived[pid] = chosen
            _log("INFO", f"phase {pid} completed: {chosen}")

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
        print(f"OK: artifact={output_path} sessionId={session.session_id}", flush=True)
        return EXIT_OK

    finally:
        with contextlib.suppress(OSError):
            shutil.rmtree(sub_build, ignore_errors=True)


def _run_tui(
    target: Optional[str],
    ntl_source_dir: Optional[Path],
    output: Optional[Path],
    phases_arg: Optional[str],
    iterations: int,
    dry_run: bool,
    resume: bool,
) -> int:
    # Pre-flight checks happen BEFORE we enter the Textual alternate
    # screen, so any failure surfaces as a normal stderr message
    # instead of a TUI flash-and-exit that looks like "the TUI didn't
    # open."

    # 1. Platform / cross-vs-native.
    plat_result = plat.check_native(target=target)
    if plat_result.kind == plat.CheckResult.CROSS_REFUSAL:
        _err(plat_result.message)
        return EXIT_CROSS_REFUSAL

    # 2. NTL source directory must contain src/.
    source_dir = _resolve_source_dir(ntl_source_dir)
    if not (source_dir / "src").exists():
        _err(
            f"NTL source dir not found or missing src/: {source_dir}. "
            f"Run `ntl-wizard` from inside an NTL source tree, or pass "
            f"--ntl-source-dir=PATH explicitly."
        )
        return EXIT_GENERIC

    # 3. Phase list.
    try:
        _resolve_phase_ids(phases_arg)
    except ValueError as exc:
        _err(str(exc))
        return EXIT_GENERIC

    # 4. libntl must already be built — `measure.py` links each timing
    #    program against the existing libntl.{so,dylib} produced by
    #    `meson compile -C build`. Without it, the very first compile
    #    fails with undefined references and the TUI would flash open
    #    and close before the user can read the error. Skip the check
    #    in dry-run mode where no compile happens.
    if not dry_run:
        libntl_candidates = [
            source_dir / "build" / "src" / "libntl.so",
            source_dir / "build" / "src" / "libntl.so.0",
            source_dir / "build" / "src" / "libntl.dylib",
            source_dir / "build" / "src" / "libntl.0.dylib",
        ]
        if not any(p.exists() for p in libntl_candidates):
            _err(
                f"libntl shared library not found under {source_dir}/build/src/. "
                f"Build NTL first so the Wizard can link timing programs against it:\n"
                f"    meson setup --buildtype=release {source_dir}/build\n"
                f"    meson compile -C {source_dir}/build\n"
                f"Then re-run ntl-wizard.\n"
                f"IMPORTANT: use --buildtype=release. The Meson default is "
                f"debug (-O0), which makes the timing benchmarks 5-10x slower "
                f"than they should be and the Wizard will take an unreasonable "
                f"amount of time to complete."
            )
            return EXIT_GENERIC

        # Probe the buildtype Meson actually used. The default 'debug'
        # buildtype produces -O0 binaries; the Wizard's timing
        # measurements become 5-10x slower and the user thinks the
        # Wizard is hung. Issue a non-fatal warning so the user can
        # rebuild before committing to a 30-minute run.
        meson_info = source_dir / "build" / "meson-info" / "intro-buildoptions.json"
        if meson_info.exists():
            try:
                import json as _json
                opts = _json.loads(meson_info.read_text(encoding="utf-8"))
                bt = next(
                    (o["value"] for o in opts if o.get("name") == "buildtype"),
                    None,
                )
                if bt is not None and bt != "release" and bt != "debugoptimized":
                    print(
                        f"warning: libntl was built with buildtype={bt!r} "
                        f"(no optimization). The Wizard's timing measurements "
                        f"will be 5-10x slower than they should be. To get "
                        f"meaningful results in reasonable time:\n"
                        f"    /bin/rm -rf {source_dir}/build\n"
                        f"    meson setup --buildtype=release {source_dir}/build\n"
                        f"    meson compile -C {source_dir}/build\n"
                        f"Continue anyway? (Ctrl-C to abort, or wait 5s.)",
                        file=sys.stderr,
                    )
                    import time as _time
                    _time.sleep(5)
            except Exception:
                pass  # warning is best-effort

    # 5. Textual must be importable.
    try:
        from . import app as _app
    except ImportError as exc:
        _err(f"TUI mode requires Textual: {exc}")
        return EXIT_GENERIC

    # Pack the args into a simple namespace-like object that app.py
    # already understands (preserve the existing app.run signature).
    import types
    ns = types.SimpleNamespace(
        target=target,
        ntl_source_dir=ntl_source_dir,
        output=output,
        phases=phases_arg,
        iterations=iterations,
        dry_run=dry_run,
        resume=resume,
    )
    return _app.run(ns)


# ---------------------------------------------------------------------------
# Typer entrypoint
# ---------------------------------------------------------------------------


cli_app = typer.Typer(
    name="ntl-wizard",
    help=(
        "NTL auto-tuning Wizard. Measures host-CPU-specific tunable "
        "parameters and writes src/meson/tune-tables/host-tuned.ini, "
        "consumed by `meson setup -Dtune=host`."
    ),
    add_completion=True,  # adds `--install-completion` and `--show-completion`
    rich_markup_mode="rich",
    # Disable Rich's fancy traceback box on errors; we want stderr to
    # follow the CLI contract's `error: <one-line>` convention.
    pretty_exceptions_enable=False,
)


@cli_app.callback(invoke_without_command=True)
def _entry(
    ctx: typer.Context,
    batch: bool = typer.Option(False, "--batch", help="Non-interactive mode; suitable for CI / headless / containers. No TUI."),
    config: Optional[Path] = typer.Option(None, "--config", help="INI file with parameter overrides or constraints."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate platform + sources + Meson tree, then exit without measuring."),
    resume: bool = typer.Option(False, "--resume", help="Continue last incomplete session for this host."),
    status: bool = typer.Option(False, "--status", help="Print last session info to stdout (as JSON), exit."),
    output: Optional[Path] = typer.Option(None, "--output", help="Override artifact path. Default: src/meson/tune-tables/host-tuned.ini."),
    ntl_source_dir: Optional[Path] = typer.Option(None, "--ntl-source-dir", help="Override the NTL source tree root. Default: autodetect from cwd."),
    target: Optional[str] = typer.Option(None, "--target", help="Assert target triplet for the cross-vs-native check."),
    phases: Optional[str] = typer.Option(None, "--phases", help="Comma-separated subset of phase ids (poly1,poly2,poly3,gf2x)."),
    iterations: int = typer.Option(1, "--iterations", help="Per-measurement iteration count."),
    seed: Optional[int] = typer.Option(None, "--seed", help="RNG seed for repeatable runs."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output (batch mode only)."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress non-error output (batch mode only)."),
    version: Optional[bool] = typer.Option(None, "--version", callback=_version_callback, is_eager=True, help="Print version and exit."),
) -> None:
    """ntl-wizard entrypoint. Dispatches to TUI (default), batch mode (--batch),
    or status query (--status)."""
    # Signal handling: Ctrl-C exits 130 (per contract).
    def _on_sigint(_signum, _frame) -> NoReturn:
        _err("interrupted by user (session state saved)")
        sys.exit(EXIT_USER_INTERRUPT)
    signal.signal(signal.SIGINT, _on_sigint)

    if status:
        raise typer.Exit(_do_status())

    if not batch:
        refusal = _check_tty_or_recommend_batch()
        if refusal is not None:
            raise typer.Exit(refusal)
        raise typer.Exit(_run_tui(
            target=target, ntl_source_dir=ntl_source_dir, output=output,
            phases_arg=phases, iterations=iterations, dry_run=dry_run, resume=resume,
        ))

    raise typer.Exit(_run_batch(
        target=target, ntl_source_dir=ntl_source_dir, output=output,
        phases_arg=phases, iterations=iterations, dry_run=dry_run, resume=resume,
        quiet=quiet,
    ))


def main(argv: Optional[list[str]] = None) -> int:
    """Entry point referenced by `[project.scripts] ntl-wizard = ...` and
    by `python -m ntl_wizard`. In normal use Typer manages sys.exit, so
    this function does not return — the SystemExit is what carries the
    exit code.

    For test contexts that want a programmatic returncode (no SystemExit),
    pass `argv` explicitly; we catch the SystemExit and return its code.
    """
    try:
        cli_app(argv if argv is not None else None)
        return EXIT_OK  # unreachable under typer in standalone mode
    except SystemExit as exc:
        return int(exc.code) if isinstance(exc.code, int) else EXIT_OK
    except KeyboardInterrupt:
        _err("interrupted by user")
        return EXIT_USER_INTERRUPT


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
