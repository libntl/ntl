#!/usr/bin/env python3
"""
T018: Load and validate an ABI table entry for a given target triplet.

Invoked from src/meson.build to resolve the per-triplet properties needed
to configure the build. Emits a series of `key=value` lines on stdout for
Meson to ingest via run_command(...).stdout(), or a clear error on stderr
with non-zero exit.

Schema is fixed (see specs/001-meson-cross-compile/contracts/abi-table.schema.md):
every key must be present; values must be drawn from the allowed sets.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Schema definition. Keep in lock-step with contracts/abi-table.schema.md.
# ---------------------------------------------------------------------------

REQUIRED_KEYS: list[tuple[str, set[str] | None]] = [
    # (key, allowed_values | None for free-form)
    ("bits_per_long", {"32", "64"}),
    ("arith_right_shift", {"0", "1"}),
    ("fma_policy", {"auto", "off"}),
    ("long_double", {"target_native", "disable"}),
    ("x86_specializations", {"true", "false"}),
    ("tune_table", {"generic", "x86", "linux-s390x"}),
    ("exec_mode", {"native", "qemu-user", "wine", "cross-only"}),
    ("exe_wrapper", None),
    ("shlib_style", {"elf", "dylib", "dll"}),
    ("threading", {"pthread", "winpthread", "none"}),
    ("tls_hack", {"true", "false"}),
]


def normalize_cpu_family(cpu: str) -> str:
    """Normalize architecture token to Meson's cpu_family vocabulary.

    Meson's host_machine.cpu_family() returns 'x86' for i386/i486/i586/i686,
    'arm' for armv6/armv7l/armv7, etc. The ABI table triplets use the
    longer form (e.g. i686-linux-gnu) but cross-key validation needs to
    compare against Meson's vocabulary.
    """
    if cpu in {"i386", "i486", "i586", "i686"}:
        return "x86"
    if cpu.startswith("armv") or cpu == "arm":
        return "arm"
    if cpu in {"powerpc64le", "ppc64le"}:
        return "ppc64"
    return cpu


def parse_triplet(triplet: str) -> tuple[str, str, str]:
    """Return (cpu_family, os, libc) given a Meson-style triplet.

    Best-effort: handles the FR-008 forms. For exotic triplets, caller is
    expected to override via the cross-file `[properties]` section.
    """

    parts = triplet.split("-")
    if len(parts) < 2:
        raise ValueError(f"Triplet {triplet!r} is not in canonical form")
    cpu = normalize_cpu_family(parts[0])
    if "apple" in parts:
        return cpu, "darwin", "darwin"
    if "w64" in parts and parts[-1].startswith("mingw"):
        return cpu, "windows", "mingw"
    if "freebsd" in triplet:
        return cpu, "freebsd", "freebsd"
    # linux-gnu, linux-musl, linux-gnueabihf-musl variants
    if "linux" in parts:
        libc = parts[-1] if parts[-1] in {"gnu", "musl"} else (
            "musl" if parts[-1].endswith("musl") else parts[-1]
        )
        return cpu, "linux", libc
    raise ValueError(f"Cannot parse triplet {triplet!r}")


def validate(
    abi: configparser.ConfigParser,
    triplet: str,
    cpu_family: str,
    os_name: str,
) -> dict[str, str]:
    if not abi.has_section("properties"):
        raise ValueError("[properties] section is missing")

    out: dict[str, str] = {}
    for key, allowed in REQUIRED_KEYS:
        if not abi.has_option("properties", key):
            raise ValueError(f"required key {key!r} is missing")
        value = abi.get("properties", key).strip()
        if allowed is not None and value not in allowed:
            raise ValueError(
                f"key {key!r} has value {value!r} which is not in "
                f"{sorted(allowed)}"
            )
        out[key] = value

    # Cross-key consistency
    if out["exec_mode"] in {"qemu-user", "wine"} and not out["exe_wrapper"]:
        raise ValueError(
            f"exec_mode={out['exec_mode']!r} requires exe_wrapper to be set"
        )
    if out["x86_specializations"] == "true" and cpu_family not in {"x86", "x86_64"}:
        raise ValueError(
            f"x86_specializations=true is incompatible with cpu_family="
            f"{cpu_family!r}"
        )
    style_for_os = {
        "linux": "elf",
        "freebsd": "elf",
        "darwin": "dylib",
        "windows": "dll",
    }
    expected_style = style_for_os.get(os_name)
    if expected_style and out["shlib_style"] != expected_style:
        raise ValueError(
            f"shlib_style={out['shlib_style']!r} is inconsistent with "
            f"os={os_name!r}; expected {expected_style!r}"
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--triplet", required=True, help="Canonical target triplet"
    )
    parser.add_argument(
        "--abi-file",
        type=Path,
        default=None,
        help="Explicit ABI table file (default: derived from --triplet and "
        "--abi-dir)",
    )
    parser.add_argument(
        "--abi-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "abi-tables",
        help="Directory containing per-triplet ABI INI files",
    )
    args = parser.parse_args()

    abi_path = args.abi_file or (args.abi_dir / f"{args.triplet}.ini")
    if not abi_path.is_file():
        print(
            f"ERROR: No ABI table entry for triplet {args.triplet!r}. "
            f"Expected file: {abi_path}. "
            f"See specs/001-meson-cross-compile/contracts/abi-table.schema.md "
            f"for the schema.",
            file=sys.stderr,
        )
        return 1

    try:
        cpu_family, os_name, _libc = parse_triplet(args.triplet)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    abi = configparser.ConfigParser()
    abi.read(abi_path, encoding="utf-8")

    try:
        properties = validate(abi, args.triplet, cpu_family, os_name)
    except ValueError as exc:
        print(
            f"ERROR: ABI table {abi_path} is invalid: {exc}",
            file=sys.stderr,
        )
        return 1

    # Emit as key=value lines, suitable for run_command().stdout() parsing.
    for key, _allowed in REQUIRED_KEYS:
        sys.stdout.write(f"{key}={properties[key]}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
