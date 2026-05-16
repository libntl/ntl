"""Frozen list of legacy-Wizard tunable parameters.

The canonical reference is `specs/002-remove-legacy-build/captured-legacy-params.txt`
(captured pre-deletion from `src/WizardAux`). This module re-declares
that list as Python data so the Wizard, the tune-table writer, and the
Meson reader all share one source of truth at runtime.

FR-005a forbids silent drops: the parity test
(`tests/ntl_wizard/test_parameters_parity.py::test_parameter_name_set_matches_legacy`)
keeps this module honest.
"""
from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Union


class ValueType(enum.Enum):
    """How the Wizard searches a parameter's value space."""
    INT_THRESHOLD = "int_threshold"
    BOOL_FLAG = "bool_flag"
    CHOICE = "choice"


ValueT = Union[int, bool, str]


@dataclass(frozen=True)
class TunableParameter:
    """One tunable Wizard parameter. Matches `data-model.md`'s entity."""
    name: str
    family: str  # phase id: "poly1", "poly2", "poly3", or "gf2x"
    value_type: ValueType
    value_domain: tuple[ValueT, ...]
    legacy_source_ref: str
    default_value: ValueT = 0

    def __post_init__(self) -> None:
        # Light internal sanity checks (the test suite is the strict
        # contract; these guard against typos at module load time).
        import re
        if not re.match(r"^NTL_[A-Z0-9_]+$", self.name):
            raise ValueError(
                f"TunableParameter name must match NTL_[A-Z0-9_]+: {self.name!r}"
            )
        if not self.value_domain:
            raise ValueError(
                f"TunableParameter {self.name}: value_domain must be non-empty"
            )
        if self.default_value not in self.value_domain:
            raise ValueError(
                f"TunableParameter {self.name}: default_value "
                f"{self.default_value!r} not in value_domain {self.value_domain!r}"
            )


# ---------------------------------------------------------------------------
# The frozen parity list. Order here defines the declaration order used by
# the artifact writer (artifacts.py) when emitting [parameters]. Reordering
# this list will produce diff churn on every Wizard run — only do it
# intentionally.
#
# Source citations point into `src/WizardAux` lines, valid as of the
# pre-deletion snapshot of NTL 11.6.0.
# ---------------------------------------------------------------------------

PARAMETERS: tuple[TunableParameter, ...] = (
    # ----- Poly1TimeTest phase (lines 80-188 of WizardAux) -----
    TunableParameter(
        name="NTL_FFT_LAZYMUL",
        family="poly1",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:117 (foreach $aflag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_SPMM_ULL",
        family="poly1",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:118 (foreach $bflag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_AVOID_BRANCHING",
        family="poly1",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:119 (foreach $cflag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_FFT_BIGTAB",
        family="poly1",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:180-188 (BIGTAB tradeoff)",
        default_value=1,
    ),

    # ----- GF2XTimeTest phase (lines 192-222 of WizardAux) -----
    TunableParameter(
        name="NTL_GF2X_NOINLINE",
        family="gf2x",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:201 (foreach $aflag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_GF2X_ALTCODE",
        family="gf2x",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:202 (foreach $bflag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_GF2X_ALTCODE1",
        family="gf2x",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:202 (foreach $bflag1)",
        default_value=0,
    ),

    # ----- Poly2TimeTest phase (lines 226-247 of WizardAux) -----
    TunableParameter(
        name="NTL_TBL_REM",
        family="poly2",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:233 (foreach $flag1)",
        default_value=0,
    ),

    # ----- Poly3TimeTest phase (lines 250-281 of WizardAux) -----
    TunableParameter(
        name="NTL_CRT_ALTCODE",
        family="poly3",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:257 (foreach $flag1)",
        default_value=0,
    ),
    TunableParameter(
        name="NTL_CRT_ALTCODE_SMALL",
        family="poly3",
        value_type=ValueType.BOOL_FLAG,
        value_domain=(0, 1),
        legacy_source_ref="src/WizardAux:275-281 (small-CRT tradeoff)",
        default_value=0,
    ),
)


# Convenience lookups.
PARAMETERS_BY_NAME: dict[str, TunableParameter] = {p.name: p for p in PARAMETERS}


def parameter_names() -> tuple[str, ...]:
    """The declared parameter names, in declaration order. Used by
    the artifact writer to determine [parameters] key ordering."""
    return tuple(p.name for p in PARAMETERS)


def parameters_for_phase(phase_id: str) -> tuple[TunableParameter, ...]:
    """All parameters whose `family` matches `phase_id`."""
    return tuple(p for p in PARAMETERS if p.family == phase_id)
