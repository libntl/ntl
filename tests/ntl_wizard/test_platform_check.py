"""T010 — Platform-check tests (FR-005d, US3).

Verifies `ntl_wizard.platform_check.check_native()` correctly accepts
native build contexts and rejects cross-compile attempts.

RED until T017 (platform_check.py) lands.
"""
from __future__ import annotations

import platform

import pytest


def test_platform_check_module_exists():
    try:
        from ntl_wizard import platform_check  # noqa: F401
    except ImportError as exc:
        pytest.fail(f"ntl_wizard.platform_check not importable: {exc}")


def test_native_x86_64_accepted():
    """When target == build host arch, check_native returns OK."""
    from ntl_wizard.platform_check import check_native, CheckResult
    host_arch = platform.machine()
    result = check_native(target=host_arch)
    assert result.kind == CheckResult.OK, (
        f"Native check ({host_arch}) should return OK; got {result!r}"
    )


def test_cross_aarch64_from_x86_64_refused():
    """Asserting target=aarch64-linux-gnu from any non-aarch64 host
    MUST return CROSS_REFUSAL with both arches named in the message."""
    from ntl_wizard.platform_check import check_native, CheckResult
    host_arch = platform.machine()
    if host_arch == "aarch64":
        pytest.skip("Can't test cross-from-x86_64 on an aarch64 host")
    result = check_native(target="aarch64-linux-gnu")
    assert result.kind == CheckResult.CROSS_REFUSAL
    assert "aarch64" in result.message
    assert host_arch in result.message or host_arch.replace("_", "-") in result.message


def test_no_target_defaults_to_build_host():
    """If --target is unset, check_native treats it as native by
    using the build host's machine."""
    from ntl_wizard.platform_check import check_native, CheckResult
    result = check_native(target=None)
    assert result.kind == CheckResult.OK


def test_normalize_armv7l_to_arm():
    """The cpu-family normalization MUST recognize armv7l <-> arm
    (consistent with feature 001's pick-abi.py normalization)."""
    from ntl_wizard.platform_check import check_native, CheckResult
    host_arch = platform.machine()
    if host_arch not in ("arm", "armv7l"):
        pytest.skip("Test runs only on arm/armv7l hosts")
    # Whichever spelling we got, the other should still be 'native'.
    other = "arm" if host_arch == "armv7l" else "armv7l"
    result = check_native(target=other)
    assert result.kind == CheckResult.OK


def test_message_includes_static_table_fallback_suggestion():
    """Cross-refusal message MUST point at the static-table fallback
    (-Dtune=...) per quickstart.md Scenario 4."""
    from ntl_wizard.platform_check import check_native, CheckResult
    result = check_native(target="riscv64-linux-gnu")
    if result.kind == CheckResult.OK:
        pytest.skip("Test runs only when build host != riscv64")
    assert "-Dtune=" in result.message or "tune=" in result.message
