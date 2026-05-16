"""T007 — Parameter parity test (FR-005a, US3).

The Python Wizard MUST cover EXACTLY the same set of tunable
parameters as the legacy `src/WizardAux`. Coverage is enforced by
comparing the Python `PARAMETERS` tuple (declared in
`tools/ntl_wizard/parameters.py`) against the frozen ground-truth
list captured in `specs/002-remove-legacy-build/captured-legacy-params.txt`
(see T002).

This test is RED until T016 (parameters.py) lands.
"""
from __future__ import annotations

import pytest


def test_parameters_module_exists():
    """The Wizard package must expose a PARAMETERS attribute."""
    try:
        from ntl_wizard import parameters
    except ImportError as exc:
        pytest.fail(f"ntl_wizard.parameters not importable: {exc}")
    assert hasattr(parameters, "PARAMETERS"), (
        "ntl_wizard.parameters must export a PARAMETERS tuple "
        "of TunableParameter instances."
    )


def test_parameter_name_set_matches_legacy(captured_legacy_params):
    """Every legacy-Wizard tunable parameter MUST appear in PARAMETERS,
    and PARAMETERS MUST NOT add parameters not present in legacy.

    This is the hard FR-005a invariant: "every tunable parameter the
    legacy `src/Wizard.cpp` exposes — none may be silently dropped".
    Both directions matter: dropping is the obvious harm, but adding
    parameters that the legacy did not have would also be a parity
    violation worth flagging (might be a typo or scope creep).
    """
    from ntl_wizard import parameters

    python_names = sorted(p.name for p in parameters.PARAMETERS)
    legacy_names = sorted(captured_legacy_params)

    missing_in_python = set(legacy_names) - set(python_names)
    extra_in_python = set(python_names) - set(legacy_names)

    if missing_in_python or extra_in_python:
        msg = ["Parameter parity violation:"]
        if missing_in_python:
            msg.append(
                f"  Missing in Python (Wizard regression): {sorted(missing_in_python)}"
            )
        if extra_in_python:
            msg.append(
                f"  Extra in Python (drift / typo / scope creep): {sorted(extra_in_python)}"
            )
        pytest.fail("\n".join(msg))


def test_parameter_names_are_valid_ntl_macros():
    """Every parameter name must look like an NTL_-prefixed macro."""
    import re
    from ntl_wizard import parameters

    pattern = re.compile(r"^NTL_[A-Z0-9_]+$")
    bad = [p.name for p in parameters.PARAMETERS if not pattern.match(p.name)]
    assert not bad, f"Parameter names not matching NTL_[A-Z0-9_]+: {bad}"


def test_parameter_value_domains_non_empty():
    """No parameter may have an empty value_domain — the Wizard's
    search loop would have nothing to try."""
    from ntl_wizard import parameters
    empty = [p.name for p in parameters.PARAMETERS if not p.value_domain]
    assert not empty, f"Parameters with empty value_domain: {empty}"


def test_parameter_defaults_are_in_domain():
    """Each parameter's default_value MUST be a member of its value_domain.
    Catches accidents like default=1 but domain=[0,2]."""
    from ntl_wizard import parameters
    mismatches = []
    for p in parameters.PARAMETERS:
        if p.default_value not in p.value_domain:
            mismatches.append((p.name, p.default_value, p.value_domain))
    assert not mismatches, (
        f"Parameters whose default_value is not in value_domain: {mismatches}"
    )


def test_parameter_families_cover_known_phases():
    """Every parameter's family must match a real MeasurementPhase id.
    The Wizard's data-model declares four phases: poly1, poly2, poly3, gf2x.
    Anything else is a parameter that no measurement phase will pick up.
    """
    from ntl_wizard import parameters
    valid_families = {"poly1", "poly2", "poly3", "gf2x"}
    bad = [
        (p.name, p.family) for p in parameters.PARAMETERS
        if p.family not in valid_families
    ]
    assert not bad, (
        f"Parameters with unknown family (not in {valid_families}): {bad}"
    )
