"""Native-vs-cross detection (FR-005d).

The Wizard runs measurement binaries on the build host. Cross-compile
contexts (build host arch ≠ target arch) MUST be refused with a clear
error pointing the user at the static tune tables.
"""
from __future__ import annotations

import enum
import platform
import subprocess
from dataclasses import dataclass


class CheckResult(enum.Enum):
    OK = "ok"
    CROSS_REFUSAL = "cross_refusal"


@dataclass(frozen=True)
class Result:
    kind: CheckResult
    build_host_arch: str
    target_arch: str
    message: str


# CPU-family normalization. Kept consistent with feature 001's
# src/meson/pick-abi.py::normalize_cpu_family. If those rules
# change, both copies need updating (single-source-of-truth refactor
# is a candidate follow-up).
_NORMALIZE_MAP = {
    # Linux uname variants
    "i386": "x86",
    "i486": "x86",
    "i586": "x86",
    "i686": "x86",
    "x86_64": "x86_64",
    "amd64": "x86_64",
    "armv7l": "arm",
    "armv6l": "arm",
    "arm": "arm",
    "aarch64": "aarch64",
    "arm64": "aarch64",
    "ppc64le": "ppc64",
    "ppc64": "ppc64",
    "riscv64": "riscv64",
    "s390x": "s390x",
}


def normalize_cpu_family(name: str) -> str:
    """Normalize a CPU name to a canonical family. Returns the input
    lowercased if no mapping is known."""
    return _NORMALIZE_MAP.get(name.lower(), name.lower())


def _extract_cpu_from_triplet(triplet: str) -> str:
    """`x86_64-linux-gnu` → `x86_64`. Best-effort."""
    return triplet.split("-", 1)[0]


def check_native(target: str | None = None) -> Result:
    """Determine whether the Wizard can run for `target`.

    Args:
        target: Either a full triplet (`aarch64-linux-gnu`), a bare
            CPU name (`aarch64`), or None. None means "no target
            specified — assume native."

    Returns:
        Result.kind == OK when target == build host arch (after
        normalization). Otherwise CROSS_REFUSAL with both arches
        named in the message.
    """
    build_host = normalize_cpu_family(platform.machine())

    if target is None:
        return Result(
            kind=CheckResult.OK,
            build_host_arch=build_host,
            target_arch=build_host,
            message=f"native build host: {build_host}",
        )

    target_cpu = _extract_cpu_from_triplet(target)
    target_norm = normalize_cpu_family(target_cpu)

    if target_norm == build_host:
        return Result(
            kind=CheckResult.OK,
            build_host_arch=build_host,
            target_arch=target_norm,
            message=f"native build host: {build_host} (target alias {target!r})",
        )

    return Result(
        kind=CheckResult.CROSS_REFUSAL,
        build_host_arch=build_host,
        target_arch=target_norm,
        message=(
            f"cross-compile context detected (build host {build_host}, "
            f"target {target_norm}). The Wizard runs measurement "
            f"binaries locally and cannot meaningfully tune for a "
            f"different architecture.\n"
            f"Use a static tune table instead, e.g.:\n"
            f"    meson setup --cross-file=... -Dtune=generic build\n"
            f"Or run the Wizard on {target_norm} hardware and commit "
            f"the resulting src/meson/tune-tables/host-tuned.ini in "
            f"your fork; then build with -Dtune=host."
        ),
    )


def cxx_compiler_dumpmachine() -> str | None:
    """Optional helper: ask the local C++ compiler for its target
    triplet via `-dumpmachine`. Used as a sanity fallback for cross
    detection. Returns None if no compiler is found.
    """
    for candidate in ("c++", "g++", "clang++"):
        try:
            result = subprocess.run(
                [candidate, "-dumpmachine"],
                capture_output=True, text=True, timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return None
