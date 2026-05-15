"""MYSTRAN-targeted BDF export translator.

Most of the work is already done by the existing
:meth:`NastranModelGenerator._write_bdf`, which honors field formats and
re-emits skipped analysis-only cards. This module adds a thin
**translation layer** on top:

  - Drops MSC/NX-only PARAMs that MYSTRAN doesn't accept.
  - Injects MYSTRAN-specific PARAMs from user preferences (SOLLIB,
    QUAD4TYP, EPSIL, WTMASS, GRDPNT).
  - Refuses to export decks with cards MYSTRAN can't handle
    (aero, nonlinear, contact, optimization) - those should have been
    flagged by pre-flight; this is the belt-and-suspenders layer.

Output is a self-contained BDF file written to the MYSTRAN run
directory. INCLUDEs are flattened (the existing read pipeline already
does this on import; on export we write what's in-memory).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

try:
    from node_runner.profiling import perf_event
except Exception:  # pragma: no cover
    def perf_event(*_a, **_kw):
        pass


# ---------------------------------------------------------------------------
# Translation table (single source of truth -- also rendered into the
# README's MYSTRAN compatibility section).
# ---------------------------------------------------------------------------

# PARAMs we silently drop (MSC/NX-only, no MYSTRAN equivalent).
_DROP_PARAMS_SILENT: frozenset[str] = frozenset({
    "POST", "K6ROT",
})

# PARAMs we drop with a warning in the export report.
_DROP_PARAMS_WARN: frozenset[str] = frozenset({
    "COUPMASS", "AUTOSPC", "BAILOUT", "MAXRATIO", "KROT",
    "PRGPST", "SRCOMPS", "NOCOMPS",
})

# PARAMs we pass through unchanged.
_PASSTHROUGH_PARAMS: frozenset[str] = frozenset({
    "WTMASS", "GRDPNT", "EPSIL",
})

# Cards MYSTRAN cannot run. If present at export time it means
# pre-flight was skipped or buggy; we drop them defensively + raise a
# warning so the user knows their deck is being modified.
_BLOCKING_CARD_FAMILIES: tuple[str, ...] = (
    "caeros", "splines", "aestats", "aesurf", "aesurfs", "trims",
    "flutters", "mkaeros", "aero", "aeros", "gusts", "paeros",
    "nlparms", "nlpcis", "bcsets", "bctables", "bcontacts", "bsurfs",
    "desvars", "dconstrs", "dresps", "dvprels", "dvcrels", "dvmrels",
    "doptprm", "dlinks",
)

# Element types MYSTRAN does not support.
_BLOCKING_ELEMENT_TYPES: frozenset[str] = frozenset({
    "CQUADR", "CTRIAR", "CWELD", "CSEAM",
})


# ---------------------------------------------------------------------------
# Translation report
# ---------------------------------------------------------------------------

@dataclass
class ExportReport:
    """Result of one ``translate_for_mystran(model, ...)`` pass."""
    n_cards_written: int = 0
    n_params_dropped: int = 0
    n_blocking_cards_dropped: int = 0
    dropped_params: list[str] = field(default_factory=list)
    injected_params: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        parts = [f"{self.n_cards_written:,} cards written"]
        if self.n_params_dropped:
            parts.append(f"{self.n_params_dropped} PARAM(s) dropped")
        if self.n_blocking_cards_dropped:
            parts.append(
                f"{self.n_blocking_cards_dropped} MYSTRAN-incompatible "
                f"card(s) dropped")
        if self.injected_params:
            parts.append(
                f"{len(self.injected_params)} MYSTRAN PARAM(s) injected")
        return "; ".join(parts)


# ---------------------------------------------------------------------------
# Core entry point
# ---------------------------------------------------------------------------

def translate_for_mystran(
        model,
        output_path: Path,
        field_format: str = "short",
        sollib: str = "SPARSE",
        quad4typ: str = "MIN4T",
        wtmass: Optional[float] = None,
        grdpnt: Optional[int] = None,
        ) -> ExportReport:
    """Apply MYSTRAN-specific translations to ``model`` in-place, then
    delegate to ``model.write_bdf(...)``.

    The model object is mutated (PARAMs dropped/added). Callers that
    need the original deck back should either work on a deep copy or
    re-apply the inverse PARAM tweaks afterwards. In practice the
    Run dialog hands us a *throwaway* model that's been deep-copied
    from the live one (see :mod:`mystran_runner`).
    """
    import time as _time
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report = ExportReport()
    _t0 = _time.perf_counter()

    # ---- 1. drop MSC/NX-only PARAMs ----
    try:
        params = getattr(model, "params", {}) or {}
        for key in list(params.keys()):
            upper = str(key).upper()
            if upper in _DROP_PARAMS_SILENT:
                del params[key]
                report.dropped_params.append(upper)
                report.n_params_dropped += 1
            elif upper in _DROP_PARAMS_WARN:
                del params[key]
                report.dropped_params.append(upper)
                report.n_params_dropped += 1
                report.warnings.append(
                    f"Dropped MSC/NX-only PARAM,{upper} "
                    f"(no MYSTRAN equivalent).")
    except Exception as exc:
        report.warnings.append(f"PARAM scan failed: {exc}")

    # ---- 2. defensive drop of blocking card families ----
    for fam in _BLOCKING_CARD_FAMILIES:
        try:
            container = getattr(model, fam, None)
            if container and len(container) > 0:
                count = len(container)
                container.clear() if hasattr(container, 'clear') else None
                report.n_blocking_cards_dropped += count
                report.warnings.append(
                    f"Dropped {count} '{fam}' card(s) - "
                    f"MYSTRAN does not support them.")
        except Exception:
            continue

    # ---- 3. drop blocking element types in-place ----
    try:
        elements = getattr(model, "elements", {}) or {}
        to_remove = [eid for eid, elem in elements.items()
                     if getattr(elem, "type", None) in _BLOCKING_ELEMENT_TYPES]
        for eid in to_remove:
            del elements[eid]
        if to_remove:
            report.n_blocking_cards_dropped += len(to_remove)
            report.warnings.append(
                f"Dropped {len(to_remove)} element(s) of types "
                f"MYSTRAN does not support "
                f"({sorted(_BLOCKING_ELEMENT_TYPES)}).")
    except Exception:
        pass

    # ---- 4. inject MYSTRAN-specific PARAMs ----
    inj = _inject_mystran_params(
        model, sollib=sollib, quad4typ=quad4typ,
        wtmass=wtmass, grdpnt=grdpnt,
    )
    report.injected_params.extend(inj)

    # ---- 5. delegate to the existing field-format-aware writer ----
    try:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        gen._write_bdf(str(output_path), field_format=field_format)
    except Exception as exc:
        report.warnings.append(f"BDF write failed: {exc}")
        raise

    # Count what we wrote (rough -- includes header lines).
    try:
        report.n_cards_written = sum(
            1 for _ in output_path.read_text(encoding='utf-8',
                                             errors='replace').splitlines())
    except Exception:
        report.n_cards_written = -1

    perf_event('mystran.export', 'translated',
               wall_s=round(_time.perf_counter() - _t0, 3),
               n_cards_written=report.n_cards_written,
               n_params_dropped=report.n_params_dropped,
               n_blocking_dropped=report.n_blocking_cards_dropped,
               n_warnings=len(report.warnings))
    return report


def _inject_mystran_params(
        model,
        sollib: str,
        quad4typ: str,
        wtmass: Optional[float],
        grdpnt: Optional[int],
        ) -> list[str]:
    """Add the MYSTRAN-required PARAMs that aren't present already.

    Returns the names of PARAMs we injected (for the export report).
    """
    from pyNastran.bdf.cards.params import PARAM
    injected: list[str] = []
    params = getattr(model, "params", None)
    if params is None:
        return injected

    def _ensure(name: str, value):
        if name in params:
            return
        try:
            params[name] = PARAM(name, [value])
            injected.append(name)
        except Exception:
            pass

    # SOLLIB / QUAD4TYP are the standard MYSTRAN-tuning PARAMs.
    #
    # PARAM,EPSIL was previously injected here but the two-field MYSTRAN
    # syntax (PARAM,EPSIL,<eqn_set_id>,<tol>) doesn't survive pyNastran's
    # single-value PARAM writer -- the value ends up in field 3 where
    # MYSTRAN expects an integer, fails on the decimal point, and the
    # solver aborts before solving. MYSTRAN's default EPSIL is fine, so
    # we no longer inject it; users who need a custom epsilon can put
    # PARAM,EPSIL on their own bulk-data card.
    _ensure("SOLLIB", str(sollib))
    _ensure("QUAD4TYP", str(quad4typ))
    if wtmass is not None:
        _ensure("WTMASS", float(wtmass))
    if grdpnt is not None:
        _ensure("GRDPNT", int(grdpnt))
    return injected


# ---------------------------------------------------------------------------
# Compatibility table (also rendered into the README at doc-gen time).
# ---------------------------------------------------------------------------

def compatibility_table() -> list[dict]:
    """Returns the translation rules as a list of dict rows. Used by
    the README auto-generator + by tests that want to assert a card is
    handled.
    """
    rows: list[dict] = []
    for name in sorted(_DROP_PARAMS_SILENT):
        rows.append({"category": "PARAM",
                     "card": f"PARAM,{name}",
                     "action": "drop (silent)",
                     "rationale": "MSC/NX-only, no MYSTRAN equivalent"})
    for name in sorted(_DROP_PARAMS_WARN):
        rows.append({"category": "PARAM",
                     "card": f"PARAM,{name}",
                     "action": "drop (warn)",
                     "rationale": "MSC/NX-only, surfaced in export report"})
    for name in sorted(_PASSTHROUGH_PARAMS):
        rows.append({"category": "PARAM",
                     "card": f"PARAM,{name}",
                     "action": "pass-through",
                     "rationale": "supported"})
    for fam in _BLOCKING_CARD_FAMILIES:
        rows.append({"category": "Card family",
                     "card": fam,
                     "action": "drop (blocking)",
                     "rationale": "MYSTRAN does not support this analysis"})
    for et in sorted(_BLOCKING_ELEMENT_TYPES):
        rows.append({"category": "Element",
                     "card": et,
                     "action": "drop (blocking)",
                     "rationale": "element not implemented in MYSTRAN"})
    return rows
