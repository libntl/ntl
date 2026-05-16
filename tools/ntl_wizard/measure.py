"""Compile + run + parse one Wizard measurement phase.

Each `MeasurementPhase` (poly1, poly2, poly3, gf2x) corresponds to a
`*TimeTest.cpp` program kept from the legacy build. The Python Wizard
orchestrates a Meson sub-build per candidate parameter set: configure
with `-DNTL_<KEY>=<VALUE>` flags, compile, run the binary, parse the
wall-clock from stdout.

Heavy operations are isolated here so the rest of the codebase remains
testable without a real toolchain.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Optional

from .parameters import PARAMETERS_BY_NAME, ValueT


# Exit codes per contracts/ntl-wizard-cli.md
class MeasureError(Exception):
    """Base class for measurement-time errors."""

    exit_code: int = 1


class CompileFailure(MeasureError):
    exit_code = 3


class RuntimeFailure(MeasureError):
    exit_code = 4


class MeasurementNoiseTooHigh(MeasureError):
    exit_code = 5


@dataclass(frozen=True)
class MeasurementPhase:
    """One compile+run step. Matches data-model.md's entity."""
    id: str  # "poly1", "poly2", "poly3", "gf2x"
    timing_program: str  # source file basename, e.g. "Poly1TimeTest.cpp"
    expected_duration_seconds: int


PHASES: tuple[MeasurementPhase, ...] = (
    MeasurementPhase("poly1", "Poly1TimeTest.cpp", 300),
    MeasurementPhase("poly2", "Poly2TimeTest.cpp", 120),
    MeasurementPhase("poly3", "Poly3TimeTest.cpp", 120),
    MeasurementPhase("gf2x",  "GF2XTimeTest.cpp",  240),
)


PHASES_BY_ID: dict[str, MeasurementPhase] = {p.id: p for p in PHASES}


@dataclass(frozen=True)
class Measurement:
    """One raw timing data point."""
    parameter_set: dict[str, ValueT]
    wall_clock_seconds: float
    iteration_count: int = 0
    noise_estimate: float = 0.0


def _stringify_for_flag(value: ValueT) -> str:
    """Render a parameter value as it appears on the gcc command line."""
    if isinstance(value, bool):
        return "1" if value else "0"
    return str(value)


def _flags_for(parameter_set: Mapping[str, ValueT]) -> list[str]:
    """Convert a parameter dict into `-DNTL_KEY=VALUE` flags."""
    flags: list[str] = []
    for name, value in parameter_set.items():
        if name not in PARAMETERS_BY_NAME:
            raise ValueError(f"Unknown parameter: {name!r}")
        flags.append(f"-D{name}={_stringify_for_flag(value)}")
    return flags


_DURATION_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*$", re.MULTILINE)


def _parse_timing_stdout(stdout: str) -> float:
    """The legacy `*TimeTest.cpp` programs emit a single number on
    stdout (microseconds-of-elapsed-time, integer). Best-effort
    parser: take the first numeric line."""
    m = _DURATION_RE.search(stdout)
    if not m:
        raise RuntimeFailure(
            f"Could not parse a numeric duration from timing-test stdout: "
            f"{stdout[:200]!r}"
        )
    return float(m.group(1))


@dataclass
class MeasureContext:
    """All the locations measure.py needs to find / write to.

    Constructed once at Wizard start; passed into `run_phase` for each
    candidate parameter set.
    """
    ntl_source_dir: Path
    sub_build_dir: Path  # ephemeral, e.g. <tmp>/ntl-wizard-build/
    compiler: str = "c++"
    extra_cxxflags: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        self.ntl_source_dir = Path(self.ntl_source_dir)
        self.sub_build_dir = Path(self.sub_build_dir)


def _run_with_streaming(
    cmd: list[str],
    timeout_seconds: int,
    line_callback: Optional[Callable[[str], None]] = None,
    tick_callback: Optional[Callable[[float], None]] = None,
    tick_interval_seconds: float = 3.0,
) -> subprocess.CompletedProcess:
    """Run `cmd` similar to subprocess.run, but:
    - stream each stdout/stderr line to `line_callback` as it arrives
      (useful for showing live compile output in the TUI),
    - call `tick_callback(elapsed)` every ~tick_interval_seconds even
      if the child is silent (heartbeat for the "is it stuck?" case),
    - kill + raise subprocess.TimeoutExpired at `timeout_seconds`.

    stderr is merged into stdout so the caller sees the natural
    interleaving (gcc warnings + linker output + program output).
    Returns a CompletedProcess whose `stdout` holds the full merged
    output (joined from the streamed lines).
    """
    start = time.monotonic()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # merge for easier streaming
        text=True,
        bufsize=1,  # line-buffered
    )

    collected: list[str] = []

    def _reader() -> None:
        # Drain the child's stdout line by line. Runs in a daemon
        # thread; survives any callback exception.
        try:
            assert proc.stdout is not None
            for raw in proc.stdout:
                collected.append(raw)
                if line_callback:
                    try:
                        line_callback(raw.rstrip("\n"))
                    except Exception:
                        pass
        except Exception:
            pass

    reader_thread = threading.Thread(target=_reader, daemon=True)
    reader_thread.start()

    last_tick = start
    try:
        while True:
            rc = proc.poll()
            now = time.monotonic()
            if rc is not None:
                break
            elapsed = now - start
            if elapsed > timeout_seconds:
                proc.kill()
                proc.wait(timeout=5)
                reader_thread.join(timeout=2)
                raise subprocess.TimeoutExpired(cmd=cmd, timeout=timeout_seconds)
            if tick_callback and (now - last_tick) >= tick_interval_seconds:
                try:
                    tick_callback(elapsed)
                except Exception:
                    pass
                last_tick = now
            time.sleep(0.1)
    finally:
        if proc.poll() is None:
            proc.kill()
    reader_thread.join(timeout=2)
    return subprocess.CompletedProcess(
        cmd, proc.returncode, "".join(collected), ""
    )


def _build_one(
    context: MeasureContext,
    phase: MeasurementPhase,
    parameter_set: Mapping[str, ValueT],
    tick_callback: Optional[Callable[[float], None]] = None,
    line_callback: Optional[Callable[[str], None]] = None,
) -> Path:
    """Compile the phase's timing program with `parameter_set` flags.
    Returns the path to the built binary. Raises CompileFailure on
    error."""
    src = context.ntl_source_dir / "src" / phase.timing_program
    if not src.exists():
        raise CompileFailure(f"Missing timing source: {src}")

    binary = context.sub_build_dir / f"{phase.id}_{int(time.time()*1000)}"

    cmd = [
        context.compiler,
        "-O2", "-std=c++11",
        f"-I{context.ntl_source_dir / 'include'}",
        f"-I{context.sub_build_dir}",  # for any generated headers Meson placed here
        *_flags_for(parameter_set),
        *context.extra_cxxflags,
        str(src),
        "-o", str(binary),
    ]
    # Also include other NTL sources that the timing program depends on.
    # The legacy build copied a fixed list (see Wizard:30-58); we
    # approximate by linking against the surrounding libntl.so if
    # present, else falling back to the full source set.
    libntl = context.ntl_source_dir / "build" / "src" / "libntl.so"
    if libntl.exists():
        cmd.extend([f"-Wl,-rpath,{libntl.parent}", str(libntl)])
    # else: leave the user to pre-build libntl. CompileFailure will fire below if missing.

    # Echo the compile invocation up-front: gcc/clang is silent on a
    # clean successful build (only diagnostics appear). Without this
    # echo the user sees only "compiling…" + elapsed ticks and might
    # think nothing is happening.
    if line_callback:
        try:
            line_callback("$ " + " ".join(cmd))
        except Exception:
            pass

    compile_start = time.monotonic()
    try:
        result = _run_with_streaming(cmd, timeout_seconds=600,
                                     tick_callback=tick_callback,
                                     line_callback=line_callback)
    except subprocess.TimeoutExpired as exc:
        raise CompileFailure(
            f"Compile of {phase.timing_program} timed out after 600s"
        ) from exc
    compile_elapsed = time.monotonic() - compile_start
    if result.returncode != 0:
        raise CompileFailure(
            f"Compile failed for {phase.id} with params {dict(parameter_set)}:\n"
            f"{result.stdout[-2000:]}"
        )
    if not binary.exists():
        raise CompileFailure(
            f"Compiler reported success but binary {binary} not produced"
        )
    if line_callback:
        try:
            line_callback(f"compile OK ({compile_elapsed:.1f}s)")
        except Exception:
            pass
    return binary


def _run_one(
    binary: Path,
    repeats: int = 1,
    tick_callback: Optional[Callable[[float], None]] = None,
    line_callback: Optional[Callable[[str], None]] = None,
) -> tuple[float, float]:
    """Execute the timing binary `repeats` times. Return
    (mean_wall_clock_seconds, stddev_seconds)."""
    times: list[float] = []
    for _ in range(repeats):
        try:
            result = _run_with_streaming([str(binary)], timeout_seconds=900,
                                         tick_callback=tick_callback,
                                         line_callback=line_callback)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeFailure(
                f"Timing binary {binary.name} timed out after 900s"
            ) from exc
        if result.returncode != 0:
            raise RuntimeFailure(
                f"Timing binary {binary.name} exited {result.returncode}:\n"
                f"{result.stdout[-1000:]}"
            )
        # legacy programs emit microseconds as an integer
        microseconds = _parse_timing_stdout(result.stdout)
        times.append(microseconds / 1_000_000.0)
    mean = sum(times) / len(times)
    if len(times) > 1:
        variance = sum((t - mean) ** 2 for t in times) / (len(times) - 1)
        stddev = variance ** 0.5
    else:
        stddev = 0.0
    return mean, stddev


def run_phase(
    context: MeasureContext,
    phase: MeasurementPhase,
    parameter_sets: list[dict[str, ValueT]],
    *,
    repeats: int = 1,
    noise_threshold: float = 0.20,
    progress_callback=None,
) -> list[Measurement]:
    """Compile + run the given phase once per parameter set.

    Returns one Measurement per parameter set. Raises:
    - CompileFailure if any compile fails (caller decides whether to
      retry-on-different-params or abort).
    - RuntimeFailure if any binary crashes.
    - MeasurementNoiseTooHigh if the relative stddev across repeats
      exceeds `noise_threshold` (default 20% — generous; the legacy
      WizardAux had no explicit noise gate).

    `progress_callback`, if provided, is called as
        progress_callback(index, total, stage, payload)
    where `stage` ∈ {"compile", "run", "done"} and `payload` is the
    parameter set being measured (for "compile" / "run") or the
    Measurement instance (for "done"). It is invoked from the calling
    thread; if the caller is the Textual worker thread, the callback
    should use call_from_thread to bounce back to the UI loop.
    """
    context.sub_build_dir.mkdir(parents=True, exist_ok=True)
    measurements: list[Measurement] = []
    total = len(parameter_sets)
    def _make_tick(stage_name: str, idx_: int) -> Optional[Callable[[float], None]]:
        if not progress_callback:
            return None
        return lambda elapsed, _i=idx_, _s=stage_name: progress_callback(
            _i, total, "tick", (_s, elapsed)
        )

    def _make_line(stage_name: str, idx_: int) -> Optional[Callable[[str], None]]:
        if not progress_callback:
            return None
        return lambda line, _i=idx_, _s=stage_name: progress_callback(
            _i, total, "line", (_s, line)
        )

    for idx, params in enumerate(parameter_sets):
        if progress_callback:
            progress_callback(idx, total, "compile", params)
        binary = _build_one(context, phase, params,
                            tick_callback=_make_tick("compile", idx),
                            line_callback=_make_line("compile", idx))
        try:
            if progress_callback:
                progress_callback(idx, total, "run", params)
            mean, stddev = _run_one(binary, repeats=repeats,
                                    tick_callback=_make_tick("run", idx),
                                    line_callback=_make_line("run", idx))
        finally:
            try:
                binary.unlink()
            except OSError:
                pass
        if mean > 0 and (stddev / mean) > noise_threshold:
            raise MeasurementNoiseTooHigh(
                f"Phase {phase.id} param set {params}: "
                f"stddev/mean = {stddev/mean:.2f} exceeds {noise_threshold:.2f}"
            )
        m = Measurement(
            parameter_set=dict(params),
            wall_clock_seconds=mean,
            iteration_count=repeats,
            noise_estimate=stddev,
        )
        measurements.append(m)
        if progress_callback:
            progress_callback(idx, total, "done", m)
    return measurements


def make_ephemeral_build_dir() -> Path:
    """Create a fresh tempdir for sub-build artifacts. Caller is
    responsible for cleanup (`shutil.rmtree`)."""
    return Path(tempfile.mkdtemp(prefix="ntl-wizard-build-"))
