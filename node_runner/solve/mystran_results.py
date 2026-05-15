"""MYSTRAN result-file adapter.

OP2-first, F06-fallback. The function :func:`load_mystran_results`
returns the same dict shape that :func:`node_runner.model.load_op2_results`
returns (``{'filepath', 'subcases': {id: {displacements, stresses,
eigenvectors, eigenvalues, frequencies, spc_forces}}}``) so the existing
Result Browser / Color-By-Results / Animation / Vector Overlay paths
work without modification.

F06 fallback is intentionally minimal in v5.0.0: only DISPLACEMENT
VECTOR and EIGENVALUE TABLE are parsed. Per-element stresses/strains
in F06 are deferred to a follow-up release. If the user's MYSTRAN
build doesn't emit OP2 and they need stresses, the Result Browser
will show a banner explaining the limitation.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from node_runner.solve import MystranRun

try:
    from node_runner.profiling import perf_event
except Exception:  # pragma: no cover
    def perf_event(*_a, **_kw):
        pass


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def load_mystran_results(run: MystranRun) -> Optional[dict]:
    """Load whatever results are available for this run.

    Tries OP2 first via pyNastran with MYSTRAN dialect enabled, falls
    back to a minimal F06 parser, returns ``None`` only if neither path
    can extract anything useful. Mutates ``run.results_source`` to
    ``"op2"`` or ``"f06"`` so callers can show a degraded-mode banner
    when relevant.
    """
    # ---- OP2 path ----
    if run.op2_path and Path(run.op2_path).exists():
        try:
            bundle = _load_op2(run.op2_path)
            run.results_source = "op2"
            perf_event('mystran.results', 'op2_loaded',
                       path=str(run.op2_path),
                       n_subcases=len(bundle.get('subcases', {})))
            return bundle
        except Exception as exc:
            perf_event('mystran.results', 'op2_failed',
                       err=str(exc)[:160])
            # fall through to F06

    # ---- F06 fallback ----
    if run.f06_path and Path(run.f06_path).exists():
        try:
            bundle = _load_f06_minimal(run.f06_path)
            run.results_source = "f06"
            perf_event('mystran.results', 'f06_fallback',
                       path=str(run.f06_path),
                       n_subcases=len(bundle.get('subcases', {})))
            return bundle
        except Exception as exc:
            perf_event('mystran.results', 'f06_failed',
                       err=str(exc)[:160])

    perf_event('mystran.results', 'no_results',
               f06=bool(run.f06_path), op2=bool(run.op2_path))
    return None


# ---------------------------------------------------------------------------
# OP2 path
# ---------------------------------------------------------------------------

def _load_op2(op2_path) -> dict:
    """Read MYSTRAN's OP2 via pyNastran and emit the shared dict shape."""
    from pyNastran.op2.op2 import OP2
    op2 = OP2(debug=False)
    # pyNastran's MYSTRAN dialect hook. Newer pyNastran versions expose
    # this as set_as_mystran(); older versions only set table-code
    # quirks via op2.set_subcases / etc. We try the modern call first,
    # then fall through.
    set_dialect = getattr(op2, 'set_as_mystran', None)
    if callable(set_dialect):
        try:
            set_dialect()
        except Exception:
            pass
    op2.read_op2(str(op2_path))
    return _op2_object_to_bundle(op2, str(op2_path))


def _op2_object_to_bundle(op2, filepath: str) -> dict:
    """Mirror the shape produced by :func:`node_runner.model.load_op2_results`."""
    results: dict = {'filepath': filepath, 'subcases': {}}

    all_sc_ids: set = set()
    all_sc_ids.update(getattr(op2, 'displacements', {}).keys())
    all_sc_ids.update(getattr(op2, 'eigenvectors', {}).keys())
    all_sc_ids.update(getattr(op2, 'spc_forces', {}).keys())
    for attr in ('cquad4_stress', 'ctria3_stress',
                 'cbar_stress', 'cbeam_stress', 'crod_stress',
                 'chexa_stress', 'ctetra_stress', 'cpenta_stress'):
        all_sc_ids.update(getattr(op2, attr, {}).keys())

    for sc_id in sorted(all_sc_ids):
        sc_data = {
            'displacements': {},
            'stresses': {},
            'eigenvectors': [],
            'spc_forces': {},
            'eigenvalues': [],
            'frequencies': [],
        }

        # Displacements (use the last time step for static SOLs).
        disp_dict = getattr(op2, 'displacements', {})
        if sc_id in disp_dict:
            disp_obj = disp_dict[sc_id]
            try:
                nids = disp_obj.node_gridtype[:, 0]
                last = disp_obj.data[-1]
                for i, nid in enumerate(nids):
                    sc_data['displacements'][int(nid)] = last[i].tolist()
            except Exception:
                pass

        # Eigenvectors / frequencies / eigenvalues (SOL 3 + SOL 5).
        eig_dict = getattr(op2, 'eigenvectors', {})
        if sc_id in eig_dict:
            eig_obj = eig_dict[sc_id]
            try:
                nids = eig_obj.node_gridtype[:, 0]
                data = eig_obj.data
                n_modes = data.shape[0]
                for mode_idx in range(n_modes):
                    mode_data = data[mode_idx]
                    mode_dict = {int(nid): mode_data[i].tolist()
                                 for i, nid in enumerate(nids)}
                    sc_data['eigenvectors'].append(mode_dict)
                eigrs = getattr(eig_obj, 'eigrs', None)
                if eigrs is not None:
                    sc_data['eigenvalues'] = [float(v) for v in eigrs]
                cycles = getattr(eig_obj, 'mode_cycles', None)
                if cycles is not None:
                    sc_data['frequencies'] = [float(v) for v in cycles]
            except Exception:
                pass

        # SPC forces (rarely needed but available).
        spc_dict = getattr(op2, 'spc_forces', {})
        if sc_id in spc_dict:
            spc_obj = spc_dict[sc_id]
            try:
                nids = spc_obj.node_gridtype[:, 0]
                last = spc_obj.data[-1]
                for i, nid in enumerate(nids):
                    sc_data['spc_forces'][int(nid)] = last[i].tolist()
            except Exception:
                pass

        # Stresses (CQUAD4 / CTRIA3 in v5.0.0; expand later).
        for attr, label in (
                ('cquad4_stress', 'CQUAD4'),
                ('ctria3_stress', 'CTRIA3')):
            stress_dict = getattr(op2, attr, {})
            if sc_id in stress_dict:
                s_obj = stress_dict[sc_id]
                try:
                    eids = s_obj.element_node[:, 0] if hasattr(
                        s_obj, 'element_node') else s_obj.element
                    data = s_obj.data[-1]
                    # MYSTRAN/MSC convention: data columns include
                    # oxx, oyy, txy, ...; we record von Mises in
                    # the last column if present.
                    for i, eid in enumerate(eids):
                        # Aggregate per-element (just the first
                        # fiber/centroid sample to keep the dict flat).
                        row = data[i].tolist()
                        sc_data['stresses'].setdefault(int(eid), {})
                        sc_data['stresses'][int(eid)]['raw'] = row
                except Exception:
                    pass

        results['subcases'][int(sc_id)] = sc_data

    return results


# ---------------------------------------------------------------------------
# F06 fallback (minimal: DISPLACEMENT VECTOR + EIGENVALUE TABLE)
# ---------------------------------------------------------------------------

# MYSTRAN's F06 uses "D I S P L A C E M E N T S" (trailing S).
# MSC/Femap uses "D I S P L A C E M E N T   V E C T O R". Accept either
# so the parser also handles hand-crafted MSC-style fixtures.
_F06_DISP_HEADER = re.compile(
    r"D\s*I\s*S\s*P\s*L\s*A\s*C\s*E\s*M\s*E\s*N\s*T\s*"
    r"(?:S|V\s*E\s*C\s*T\s*O\s*R)",
    re.IGNORECASE,
)
_F06_EIGS_HEADER = re.compile(
    r"R\s*E\s*A\s*L\s+E\s*I\s*G\s*E\s*N\s*V\s*A\s*L\s*U\s*E\s*S",
    re.IGNORECASE,
)
# Sections that share the disp row format but are NOT displacements
# -- detecting them as their own section ensures we don't accidentally
# overwrite the disp data when MYSTRAN emits SPC / MPC forces inline.
_F06_OTHER_NODE_TABLE = re.compile(
    r"(S\s*P\s*C\s+F\s*O\s*R\s*C\s*E\s*S|"
    r"M\s*P\s*C\s+F\s*O\s*R\s*C\s*E\s*S|"
    r"L\s*O\s*A\s*D\s+V\s*E\s*C\s*T\s*O\s*R|"
    r"A\s*P\s*P\s*L\s*I\s*E\s*D\s+L\s*O\s*A\s*D)",
    re.IGNORECASE,
)
_F06_SUBCASE = re.compile(
    r"SUBCASE\s+(\d+)\s*$", re.IGNORECASE)
# MYSTRAN's F06 also writes the subcase id in the OUTPUT FOR SUBCASE
# block right above each result table. Capture that form too so we
# don't fall back to subcase 1 when MYSTRAN dropped the case-control
# echo.
_F06_OUTPUT_FOR_SUBCASE = re.compile(
    r"OUTPUT\s+FOR\s+SUBCASE\s+(\d+)", re.IGNORECASE)

# Float pattern accepts plain (0.0), fixed (1.234), AND scientific
# (1.234E-05). MYSTRAN writes exact zeros as bare "0.0" and nonzero
# values in scientific notation, so a stricter "always-scientific"
# regex misses every row that has any zero component.
_F = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"

# Row formats supported:
#
# 1) MYSTRAN: "<GID>  <COORD_SYS>  <T1>  <T2>  <T3>  <R1>  <R2>  <R3>"
#       1        0  0.0           0.0          -6.181010E-04  0.0    ...
#    Field 2 is the integer coord-system id (0 = global).
#
# 2) MSC/Femap: "<GID>  G  <T1>  <T2>  <T3>  <R1>  <R2>  <R3>"
#       1       G   1.234567E-04  ...
#    Field 2 is the literal letter "G" marking grid type (vs "S" for
#    scalar point).
_DISP_ROW = re.compile(
    r"^\s*(\d+)\s+"          # GID
    r"(?:G|\d+)\s+"          # MSC type letter OR MYSTRAN coord-sys id
    r"(" + _F + r")\s+"      # T1
    r"(" + _F + r")\s+"      # T2
    r"(" + _F + r")\s+"      # T3
    r"(" + _F + r")\s+"      # R1
    r"(" + _F + r")\s+"      # R2
    r"(" + _F + r")"         # R3
    r"\s*$"
)

# Eigenvalue row: MODE  EXTRACT_ORDER  EIGENVALUE  RADIANS  CYCLES  ...
# Same generic float pattern -- MYSTRAN's eigenvalue table also mixes
# plain and scientific floats.
_EIG_ROW = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+"
    r"(" + _F + r")\s+"
    r"(" + _F + r")\s+"
    r"(" + _F + r")"
)


def _load_f06_minimal(f06_path) -> dict:
    """Parse displacements + eigenvalues from a MYSTRAN F06.

    This is intentionally minimal -- per-element stress/strain/force
    parsers are out of scope for v5.0.0. Callers should rely on OP2 for
    stresses; the Result Browser will note when only F06 fallback is
    available.
    """
    f06_path = Path(f06_path)
    text = f06_path.read_text(encoding='utf-8', errors='replace')

    subcases: dict[int, dict] = {}
    cur_sc = 1
    cur_section: Optional[str] = None  # 'disp' or 'eig'
    blank_streak = 0  # >=2 consecutive blank lines ends a section

    for line in text.splitlines():
        # Subcase tracking: case-control echo "SUBCASE N" form...
        m = _F06_SUBCASE.search(line)
        if m:
            cur_sc = int(m.group(1))
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        # ...and MYSTRAN's "OUTPUT FOR SUBCASE N" header above each
        # result block.
        m2 = _F06_OUTPUT_FOR_SUBCASE.search(line)
        if m2:
            cur_sc = int(m2.group(1))
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        # Section headers (also reset the blank streak).
        if _F06_DISP_HEADER.search(line):
            cur_section = 'disp'
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        if _F06_EIGS_HEADER.search(line):
            cur_section = 'eig'
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        if _F06_OTHER_NODE_TABLE.search(line):
            # SPC FORCES / MPC FORCES / LOAD VECTOR share the
            # disp-row syntax; treat them as "not displacements"
            # so we don't overwrite. v5.0.0 only consumes disp +
            # eigenvalues from F06; expand here when adding more.
            cur_section = 'other_node_table'
            blank_streak = 0
            continue
        # Multi-blank gap ends a section (page break / table boundary).
        if not line.strip():
            blank_streak += 1
            if blank_streak >= 2:
                cur_section = None
            continue
        blank_streak = 0

        if cur_section == 'disp':
            md = _DISP_ROW.match(line)
            if md:
                nid = int(md.group(1))
                vals = [float(md.group(i)) for i in range(2, 8)]
                subcases[cur_sc]['displacements'][nid] = vals
        elif cur_section == 'eig':
            me = _EIG_ROW.match(line)
            if me:
                eigenvalue = float(me.group(3))
                radians = float(me.group(4))
                cycles = float(me.group(5))
                subcases[cur_sc]['eigenvalues'].append(eigenvalue)
                subcases[cur_sc]['frequencies'].append(cycles)
                # No eigenvector data in this minimal path.
                subcases[cur_sc]['eigenvectors'].append({})

    # Drop subcases that contain no real data -- the parser would
    # otherwise return a "phantom" subcase 1 just because it saw a
    # SUBCASE header in the case-control echo at the top of the F06.
    # Callers use the emptiness check to distinguish "MYSTRAN ran +
    # produced output" from "MYSTRAN echoed the deck but bailed".
    real = {sid: sc for sid, sc in subcases.items()
            if sc['displacements'] or sc['eigenvalues']
            or sc['eigenvectors'] or sc['spc_forces']}
    return {'filepath': str(f06_path), 'subcases': real}


def _empty_subcase() -> dict:
    return {
        'displacements': {},
        'stresses': {},
        'eigenvectors': [],
        'spc_forces': {},
        'eigenvalues': [],
        'frequencies': [],
    }
