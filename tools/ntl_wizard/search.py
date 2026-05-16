"""Parameter-search algorithm for each measurement phase.

Replicates the legacy WizardAux logic (src/WizardAux:107-281): for each
phase, try a small Cartesian product of value sets and pick the one
that minimizes wall-clock. The legacy code uses a few heuristic skips
(e.g. NTL_FFT_LAZYMUL requires NTL_SPMM_ULL or NTL_LONGLONG_SP_MULMOD)
which we reproduce faithfully — the goal is parameter parity, not
algorithmic novelty.
"""
from __future__ import annotations

import itertools
from typing import Iterable, Mapping

from .measure import Measurement
from .parameters import (
    PARAMETERS_BY_NAME,
    TunableParameter,
    ValueT,
    parameters_for_phase,
)


def candidate_sets_for_phase(
    phase_id: str,
    host_features: Mapping[str, int] | None = None,
) -> list[dict[str, ValueT]]:
    """Generate every candidate parameter set the Wizard should try
    for a given phase.

    Args:
        phase_id: "poly1", "poly2", "poly3", or "gf2x".
        host_features: Detected host-feature flags from
            `InitSettings` (`NTL_HAVE_LL_TYPE`, `NTL_HAVE_PCLMUL`,
            `NTL_LONGLONG_SP_MULMOD`, `NTL_GMP_LIP`). Some legacy
            phases are skipped entirely when prerequisites aren't met.

    Returns:
        A list of parameter dicts. Each dict has every PARAMETERS
        entry for the phase, with one of its value_domain choices.
    """
    host_features = dict(host_features) if host_features else {}
    params = parameters_for_phase(phase_id)
    if not params:
        return []

    # Brute-force Cartesian product, then filter by phase-specific rules
    # that match the legacy WizardAux logic.
    domains = [(p, p.value_domain) for p in params]
    candidates: list[dict[str, ValueT]] = []
    for combo in itertools.product(*[d for _, d in domains]):
        candidate = {p.name: v for (p, _), v in zip(domains, combo)}
        if not _legacy_skip(phase_id, candidate, host_features):
            candidates.append(candidate)
    return candidates


def _legacy_skip(
    phase_id: str,
    candidate: Mapping[str, ValueT],
    host_features: Mapping[str, int],
) -> bool:
    """Replicate the `$skipit` heuristics from src/WizardAux.

    Returns True if this candidate should be skipped (i.e. not
    measured) — the legacy logic uses these to avoid known-bad
    or pointless-to-measure combinations.
    """
    if phase_id == "poly1":
        # WizardAux:128-143
        ll_type = host_features.get("NTL_HAVE_LL_TYPE", 0)
        long_sp = host_features.get("NTL_LONGLONG_SP_MULMOD", 0)
        spmm = candidate.get("NTL_SPMM_ULL", 0)
        lazy = candidate.get("NTL_FFT_LAZYMUL", 0)
        if spmm == 1 and ll_type == 0:
            return True  # skip1
        if lazy == 1 and not (long_sp == 1 or spmm == 1):
            return True  # skip2
        if long_sp == 1 and spmm == 1:
            return True  # skip3

    if phase_id == "gf2x":
        # WizardAux:195: "but not if we have PCLMUL"
        if host_features.get("NTL_HAVE_PCLMUL", 0) == 1:
            return True  # entire phase skipped

    if phase_id == "poly2":
        # WizardAux:228: only runs if NTL_HAVE_LL_TYPE == 1
        if host_features.get("NTL_HAVE_LL_TYPE", 0) != 1:
            return True

    if phase_id == "poly3":
        # WizardAux:252: only runs if NTL_HAVE_LL_TYPE == 1 AND NTL_GMP_LIP == 1
        if host_features.get("NTL_HAVE_LL_TYPE", 0) != 1:
            return True
        if host_features.get("NTL_GMP_LIP", 0) != 1:
            return True

    return False


def derive_values(
    phase_id: str,
    measurements: list[Measurement],
) -> dict[str, ValueT]:
    """Pick the parameter values for `phase_id` from a list of
    measurements. The winning parameter set is the one with the
    smallest wall-clock.

    Returns:
        dict mapping parameter name → chosen value, for every
        TunableParameter in the phase's family.
    """
    if not measurements:
        # Phase was entirely skipped (e.g. gf2x when PCLMUL is present).
        # Fall back to defaults.
        return {p.name: p.default_value for p in parameters_for_phase(phase_id)}

    best = min(measurements, key=lambda m: m.wall_clock_seconds)
    chosen: dict[str, ValueT] = dict(best.parameter_set)

    # Special case: WizardAux:275-281 — even if NTL_CRT_ALTCODE didn't
    # win outright, set NTL_CRT_ALTCODE_SMALL if the altcode result
    # wasn't *too* bad (within 15% of best). We approximate by
    # comparing the altcode vs no-altcode measurements.
    if phase_id == "poly3":
        no_alt = next((m for m in measurements if m.parameter_set.get("NTL_CRT_ALTCODE") == 0), None)
        alt = next((m for m in measurements if m.parameter_set.get("NTL_CRT_ALTCODE") == 1), None)
        if no_alt and alt and chosen.get("NTL_CRT_ALTCODE") == 0:
            if alt.wall_clock_seconds <= 1.15 * no_alt.wall_clock_seconds:
                chosen["NTL_CRT_ALTCODE_SMALL"] = 1

    return chosen
