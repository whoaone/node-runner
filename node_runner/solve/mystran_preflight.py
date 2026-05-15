"""MYSTRAN pre-flight compatibility scanner.

Walks the in-memory pyNastran BDF and produces a :class:`PreflightReport`
that the Run dialog inspects before launching MYSTRAN.

Designed for a < 1 s budget on a 2.4M-element deck:

  - Card-family checks (aero/nonlinear/contact/opt) go through the
    pyNastran BDF's dedicated attribute dicts: O(1) per family.
  - Element-type checks iterate ``model.elements.values()`` once but
    do nothing more than a set-membership lookup per element.
  - PARAM checks iterate ``model.params`` dict: O(family-count).
  - DOF estimate is ``6 * len(model.nodes)``; no SPC walk required.
"""

from __future__ import annotations

import time as _time
from typing import Iterable

from node_runner.solve import (
    PreflightIssue, PreflightReport, SUPPORTED_SOL_NUMBERS,
)
try:
    from node_runner.profiling import perf_event
except Exception:  # pragma: no cover
    def perf_event(*_a, **_kw):
        pass


# ---------------------------------------------------------------------------
# Lookup tables
# ---------------------------------------------------------------------------

# Element types MYSTRAN does not support. Match by `elem.type`.
_UNSUPPORTED_ELEMENT_TYPES: frozenset[str] = frozenset({
    "CQUADR", "CTRIAR", "CWELD", "CSEAM",
})

# (attribute_on_BDF, severity, code, message_fmt)
# message_fmt is .format(count=...) -- safe to call with the dict
# length even when the attribute is missing (we treat None as 0).
_CARD_FAMILY_CHECKS: tuple[tuple[str, str, str, str], ...] = (
    # Aerodynamics
    ("caeros",   "blocking", "AERO_CAERO",
     "Deck contains {count} CAERO* aero panel(s); MYSTRAN has no aero solver."),
    ("splines",  "blocking", "AERO_SPLINE",
     "Deck contains {count} SPLINE* aero spline(s); MYSTRAN has no aero solver."),
    ("aestats",  "blocking", "AERO_AESTAT",
     "Deck contains {count} AESTAT control-surface state(s); MYSTRAN has no aero solver."),
    ("trims",    "blocking", "AERO_TRIM",
     "Deck contains {count} TRIM card(s); MYSTRAN has no aero solver."),
    ("flutters", "blocking", "AERO_FLUTTER",
     "Deck contains {count} FLUTTER card(s); MYSTRAN has no aero solver."),
    ("mkaeros",  "blocking", "AERO_MKAERO",
     "Deck contains {count} MKAERO card(s); MYSTRAN has no aero solver."),
    ("gusts",    "blocking", "AERO_GUST",
     "Deck contains {count} GUST card(s); MYSTRAN has no aero solver."),
    # Nonlinear
    ("nlparms",  "blocking", "NL_NLPARM",
     "Deck contains {count} NLPARM card(s); MYSTRAN does not run nonlinear analysis."),
    ("nlpcis",   "blocking", "NL_NLPCI",
     "Deck contains {count} NLPCI card(s); MYSTRAN does not run nonlinear analysis."),
    # Contact
    ("bcsets",   "blocking", "CONTACT_BCSET",
     "Deck contains {count} BCSET / contact card(s); MYSTRAN does not run contact."),
    ("bctables", "blocking", "CONTACT_BCTABLE",
     "Deck contains {count} BCTABLE card(s); MYSTRAN does not run contact."),
    ("bcontacts", "blocking", "CONTACT_BCONTACT",
     "Deck contains {count} BCONTACT card(s); MYSTRAN does not run contact."),
    # Optimization
    ("desvars",  "blocking", "OPT_DESVAR",
     "Deck contains {count} DESVAR card(s); MYSTRAN does not run optimization."),
    ("dconstrs", "blocking", "OPT_DCONSTR",
     "Deck contains {count} DCONSTR card(s); MYSTRAN does not run optimization."),
    ("dresps",   "blocking", "OPT_DRESP",
     "Deck contains {count} DRESP card(s); MYSTRAN does not run optimization."),
)

# PARAMs that MYSTRAN doesn't accept; we surface a warning so the user
# knows the translation layer will drop them.
_WARN_PARAMS: frozenset[str] = frozenset({
    "POST", "COUPMASS", "K6ROT", "AUTOSPC", "BAILOUT",
    "MAXRATIO", "KROT", "PRGPST", "SRCOMPS", "NOCOMPS",
})

# Rough wall-time heuristic for SOL 1: ~50 ns per DOF for IntMKL.
_WALL_TIME_PER_DOF_SECONDS = 5.0e-8


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def scan_for_mystran(model, sol: int | None = None) -> PreflightReport:
    """Run all MYSTRAN-compatibility checks on ``model``.

    ``sol`` is the solution number the user has selected for this run
    (may be None; defaults to whatever model.sol reports).
    """
    t0 = _time.perf_counter()
    report = PreflightReport(solution_supported=False)

    # ---- 1. SOL ----
    actual_sol = sol if sol is not None else getattr(model, "sol", None)
    if actual_sol is None:
        report.issues.append(PreflightIssue(
            severity="warning",
            code="SOL_MISSING",
            message="No SOL number found in the deck or run options; "
                    "defaulting to SOL 1 (linear static).",
        ))
        actual_sol = 1
    if actual_sol not in SUPPORTED_SOL_NUMBERS:
        report.issues.append(PreflightIssue(
            severity="blocking",
            code="SOL_UNSUPPORTED",
            message=f"SOL {actual_sol} is not supported by MYSTRAN. "
                    f"Supported: SOL 1 / 3 / 5 (or rigid-format "
                    f"equivalents 101 / 103 / 105).",
        ))
    else:
        report.solution_supported = True

    # ---- 2. Card-family checks ----
    for attr, severity, code, fmt in _CARD_FAMILY_CHECKS:
        try:
            container = getattr(model, attr, None) or {}
            count = len(container)
        except Exception:
            count = 0
        if count > 0:
            report.issues.append(PreflightIssue(
                severity=severity,
                code=code,
                message=fmt.format(count=count),
                card_type=attr,
                card_count=count,
            ))

    # ---- 3. PARAMs (one issue per WARN param present) ----
    try:
        params = getattr(model, "params", {}) or {}
        for name in params.keys():
            upper = str(name).upper()
            if upper in _WARN_PARAMS:
                report.issues.append(PreflightIssue(
                    severity="warning",
                    code=f"PARAM_{upper}",
                    message=f"PARAM,{upper} is MSC/NX-specific and will "
                            f"be dropped from the MYSTRAN deck.",
                    card_type=f"PARAM,{upper}",
                ))
    except Exception:
        pass

    # ---- 4. Element-type single-pass ----
    unsupported_counts: dict[str, int] = {}
    try:
        elements = getattr(model, "elements", {}) or {}
        for elem in elements.values():
            t = getattr(elem, "type", None)
            if t in _UNSUPPORTED_ELEMENT_TYPES:
                unsupported_counts[t] = unsupported_counts.get(t, 0) + 1
    except Exception:
        pass
    for etype, count in sorted(unsupported_counts.items()):
        report.issues.append(PreflightIssue(
            severity="blocking",
            code=f"ELEM_{etype}",
            message=f"Deck contains {count:,} {etype} element(s); "
                    f"MYSTRAN does not implement this element type.",
            card_type=etype,
            card_count=count,
        ))

    # ---- 5. DOF / wall-time estimate ----
    try:
        n_nodes = len(getattr(model, "nodes", {}) or {})
    except Exception:
        n_nodes = 0
    report.estimated_dof = 6 * n_nodes
    if report.estimated_dof > 0 and actual_sol in (1, 101):
        report.estimated_wall_time_seconds = (
            report.estimated_dof * _WALL_TIME_PER_DOF_SECONDS)

    # ---- 6. Large-deck heads-up (info-level) ----
    if report.estimated_dof > 500_000:
        report.issues.append(PreflightIssue(
            severity="info",
            code="LARGE_DOF",
            message=f"Deck has ~{report.estimated_dof:,} DOF. "
                    f"MYSTRAN will run but expect a longer solve time.",
        ))

    report.scan_wall_s = round(_time.perf_counter() - t0, 4)
    perf_event('mystran.preflight', 'scan',
               wall_s=report.scan_wall_s,
               n_blocking=report.blocking_count,
               n_warnings=report.warning_count,
               n_info=report.info_count,
               estimated_dof=report.estimated_dof)
    return report
