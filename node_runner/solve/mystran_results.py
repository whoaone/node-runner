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
# MSC Nastran uses "D I S P L A C E M E N T   V E C T O R". Accept either
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
# v5.3.0 patch 2 Q5: split SPC FORCES into its own section so the
# disp-style rows actually land in the bundle's `spc_forces` dict.
# Pre-patch the parser bucketed everything in this header into a
# generic "other_node_table" state that just skipped the rows.
_F06_SPCFORCE_HEADER = re.compile(
    r"S\s*P\s*C\s+F\s*O\s*R\s*C\s*E\s*S",
    re.IGNORECASE,
)
# Sections that share the disp row format but are NOT displacements
# and NOT SPC forces. Detecting them as their own section ensures we
# don't accidentally overwrite the disp data when MYSTRAN emits MPC
# forces or load vectors inline. SPC FORCES is handled separately
# above; MPC forces + load vector + applied load remain skipped until
# the bundle grows slots for them.
_F06_OTHER_NODE_TABLE = re.compile(
    r"(M\s*P\s*C\s+F\s*O\s*R\s*C\s*E\s*S|"
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

# v5.3.0 P1B: MYSTRAN F06 stress section banners. The element-type
# subheader determines which parser engages. MYSTRAN's banner text:
#   E L E M E N T   S T R E S S E S   I N   L O C A L   E L E M E N T
#         C O O R D I N A T E   S Y S T E M
#         F O R   E L E M E N T   T Y P E   Q U A D 4
# We accept the broad header as the gate and then look for the per-
# element-type subheader.
_F06_STRESS_HEADER = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+"
    r"S\s*T\s*R\s*E\s*S\s*S\s*E\s*S",
    re.IGNORECASE,
)
_F06_STRESS_QUAD4 = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+T\s*Y\s*P\s*E\s+"
    r"Q\s*U\s*A\s*D\s*4",
    re.IGNORECASE,
)
_F06_STRESS_TRIA3 = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+T\s*Y\s*P\s*E\s+"
    r"T\s*R\s*I\s*A\s*3",
    re.IGNORECASE,
)
_F06_STRESS_BAR = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+T\s*Y\s*P\s*E\s+"
    r"B\s*A\s*R",
    re.IGNORECASE,
)
# v5.3.0 patch 2 Q4: element engineering forces section.
_F06_FORCES_HEADER = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+"
    r"E\s*N\s*G\s*I\s*N\s*E\s*E\s*R\s*I\s*N\s*G\s+"
    r"F\s*O\s*R\s*C\s*E\s*S",
    re.IGNORECASE,
)
_F06_FORCES_BAR = re.compile(
    r"E\s*L\s*E\s*M\s*E\s*N\s*T\s+T\s*Y\s*P\s*E\s+"
    r"B\s*A\s*R",
    re.IGNORECASE,
)

# Float pattern accepts plain (0.0), fixed (1.234), AND scientific
# (1.234E-05). MYSTRAN writes exact zeros as bare "0.0" and nonzero
# values in scientific notation, so a stricter "always-scientific"
# regex misses every row that has any zero component.
_F = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"

# v5.3.0 P1B (row regexes that depend on _F).
# QUAD4 stress data rows have two flavors:
#   Top fiber:    "  EID  CENTER  -2.5E-02  Nx  Ny  Sxy  Angle  Max  Min  VonMises  Sxz  Syz"
#   Bottom fiber: "          2.5E-02  Nx  Ny  Sxy  Angle  Max  Min  VonMises"
# The top row starts with an integer EID; the bottom row starts with
# whitespace + the fiber distance float. Both have 6-9 stress columns.
_F06_QUAD4_TOP_ROW = re.compile(
    r"^\s*(\d+)\s+(?:CENTER|CORNER\s+\d+)\s+"
    r"(" + _F + r")\s+"          # fiber distance (Z)
    r"(" + _F + r")\s+"          # Normal-X
    r"(" + _F + r")\s+"          # Normal-Y
    r"(" + _F + r")\s+"          # Shear-XY
    r"(" + _F + r")\s+"          # Angle
    r"(" + _F + r")\s+"          # Major principal
    r"(" + _F + r")\s+"          # Minor principal
    r"(" + _F + r")"             # Von Mises
)
_F06_QUAD4_BOT_ROW = re.compile(
    r"^\s+"
    r"(" + _F + r")\s+"          # fiber distance (Z), no EID
    r"(" + _F + r")\s+"          # Normal-X
    r"(" + _F + r")\s+"          # Normal-Y
    r"(" + _F + r")\s+"          # Shear-XY
    r"(" + _F + r")\s+"          # Angle
    r"(" + _F + r")\s+"          # Major principal
    r"(" + _F + r")\s+"          # Minor principal
    r"(" + _F + r")"             # Von Mises
    r"\s*$"
)

# v5.3.0 patch 2 Q2: CTRIA3 stress rows. Top row uses "Anywhere"
# keyword (vs CQUAD4's "CENTER"); bottom row uses "in elem" continuation
# marker.
_F06_TRIA3_TOP_ROW = re.compile(
    r"^\s*(\d+)\s+Anywhere\s+"
    r"(" + _F + r")\s+"          # fiber distance (Z)
    r"(" + _F + r")\s+"          # Normal-X
    r"(" + _F + r")\s+"          # Normal-Y
    r"(" + _F + r")\s+"          # Shear-XY
    r"(" + _F + r")\s+"          # Angle
    r"(" + _F + r")\s+"          # Major principal
    r"(" + _F + r")\s+"          # Minor principal
    r"(" + _F + r")"             # Von Mises
)
_F06_TRIA3_BOT_ROW = re.compile(
    r"^\s+in\s+elem\s+"
    r"(" + _F + r")\s+"          # fiber distance
    r"(" + _F + r")\s+"          # Normal-X
    r"(" + _F + r")\s+"          # Normal-Y
    r"(" + _F + r")\s+"          # Shear-XY
    r"(" + _F + r")\s+"          # Angle
    r"(" + _F + r")\s+"          # Major principal
    r"(" + _F + r")\s+"          # Minor principal
    r"(" + _F + r")"             # Von Mises
    r"\s*$"
)

# v5.3.0 patch 2 Q3: CBAR stress rows. End-A row carries EID + 4 fiber
# stresses + axial + max/min + margin of safety; End-B row drops the
# axial column and uses M.S.-C instead of M.S.-T. Both can have a blank
# M.S. column when MYSTRAN couldn't compute it. We accept either a
# float OR the literal "----" placeholder for that column.
_MS_FIELD = r"(?:" + _F + r"|-+)"

_F06_BAR_END_A_ROW = re.compile(
    r"^\s*(\d+)\s+"
    r"(" + _F + r")\s+"          # SA1
    r"(" + _F + r")\s+"          # SA2
    r"(" + _F + r")\s+"          # SA3
    r"(" + _F + r")\s+"          # SA4
    r"(" + _F + r")\s+"          # Axial stress
    r"(" + _F + r")\s+"          # SA-Max
    r"(" + _F + r")"             # SA-Min
    r"(?:\s+" + _MS_FIELD + r")?"  # M.S.-T (optional)
    r"\s*$"
)
_F06_BAR_END_B_ROW = re.compile(
    r"^\s+"
    r"(" + _F + r")\s+"          # SB1
    r"(" + _F + r")\s+"          # SB2
    r"(" + _F + r")\s+"          # SB3
    r"(" + _F + r")\s+"          # SB4
    r"(" + _F + r")\s+"          # SB-Max
    r"(" + _F + r")"             # SB-Min
    r"(?:\s+" + _MS_FIELD + r")?"  # M.S.-C (optional)
    r"\s*$"
)

# v5.3.0 patch 2 Q4: CBAR element engineering forces. Single row per
# element with 8 columns after EID.
_F06_BAR_FORCES_ROW = re.compile(
    r"^\s*(\d+)\s+"
    r"(" + _F + r")\s+"          # M1-A (Plane 1, End A)
    r"(" + _F + r")\s+"          # M2-A (Plane 2, End A)
    r"(" + _F + r")\s+"          # M1-B (Plane 1, End B)
    r"(" + _F + r")\s+"          # M2-B (Plane 2, End B)
    r"(" + _F + r")\s+"          # V1 (Shear Plane 1)
    r"(" + _F + r")\s+"          # V2 (Shear Plane 2)
    r"(" + _F + r")\s+"          # Axial force
    r"(" + _F + r")"             # Torque
    r"\s*$"
)

# Row formats supported:
#
# 1) MYSTRAN: "<GID>  <COORD_SYS>  <T1>  <T2>  <T3>  <R1>  <R2>  <R3>"
#       1        0  0.0           0.0          -6.181010E-04  0.0    ...
#    Field 2 is the integer coord-system id (0 = global).
#
# 2) MSC Nastran: "<GID>  G  <T1>  <T2>  <T3>  <R1>  <R2>  <R3>"
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

    v5.3.0 P1B: also parses CQUAD4 stresses (top/bottom fiber, with
    Von Mises taken as the max of the two fibers).
    Other element types (CTRIA3 / CBEAM / CROD) and strain/force
    extraction remain out of scope; expansion planned for v5.3.1.
    """
    f06_path = Path(f06_path)
    text = f06_path.read_text(encoding='utf-8', errors='replace')

    subcases: dict[int, dict] = {}
    cur_sc = 1
    cur_section: Optional[str] = None
    # ^ Possible states: 'disp', 'eig', 'spcforce', 'other_node_table',
    #   'stress_pending', 'stress_quad4', 'stress_tria3', 'stress_bar',
    #   'forces_pending', 'forces_bar'.
    cur_stress_eid: Optional[int] = None  # the EID currently being filled
    blank_streak = 0  # >=2 consecutive blank lines ends a section
    n_quad4_stress_parsed = 0
    n_quad4_stress_unmatched = 0
    n_tria3_stress_parsed = 0
    n_bar_stress_parsed = 0
    n_bar_forces_parsed = 0

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
        if _F06_STRESS_HEADER.search(line):
            # Tentative stress section -- subtype set by the next
            # element-type subheader line.
            cur_section = 'stress_pending'
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        # v5.3.0 patch 2: stress / force section dispatch. Each
        # broad-banner block is followed by an element-type subheader
        # that locks in the parser state for the rows that follow.
        if cur_section in ('stress_pending', 'stress_quad4',
                            'stress_tria3', 'stress_bar'):
            if _F06_STRESS_QUAD4.search(line):
                cur_section = 'stress_quad4'
                cur_stress_eid = None
                blank_streak = 0
                continue
            if _F06_STRESS_TRIA3.search(line):
                cur_section = 'stress_tria3'
                cur_stress_eid = None
                blank_streak = 0
                continue
            if _F06_STRESS_BAR.search(line):
                cur_section = 'stress_bar'
                cur_stress_eid = None
                blank_streak = 0
                continue
        if _F06_FORCES_HEADER.search(line):
            cur_section = 'forces_pending'
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        if cur_section in ('forces_pending', 'forces_bar'):
            if _F06_FORCES_BAR.search(line):
                cur_section = 'forces_bar'
                blank_streak = 0
                continue
        # v5.3.0 patch 2 Q5: SPC FORCES section -- same row format as
        # the displacement vector. Populate `spc_forces` dict.
        if _F06_SPCFORCE_HEADER.search(line):
            cur_section = 'spcforce'
            subcases.setdefault(cur_sc, _empty_subcase())
            blank_streak = 0
            continue
        if _F06_OTHER_NODE_TABLE.search(line):
            # MPC FORCES / LOAD VECTOR / APPLIED LOAD share the disp
            # row format but lack bundle slots; skip until they're
            # added.
            cur_section = 'other_node_table'
            blank_streak = 0
            continue
        # Multi-blank gap ends a section (page break / table boundary).
        if not line.strip():
            blank_streak += 1
            if blank_streak >= 2:
                cur_section = None
                cur_stress_eid = None
            continue
        blank_streak = 0

        if cur_section == 'disp':
            md = _DISP_ROW.match(line)
            if md:
                nid = int(md.group(1))
                vals = [float(md.group(i)) for i in range(2, 8)]
                subcases[cur_sc]['displacements'][nid] = vals
        elif cur_section == 'spcforce':
            # v5.3.0 patch 2 Q5: SPC forces share the disp row format.
            md = _DISP_ROW.match(line)
            if md:
                nid = int(md.group(1))
                vals = [float(md.group(i)) for i in range(2, 8)]
                subcases[cur_sc]['spc_forces'][nid] = vals
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
        elif cur_section == 'stress_quad4':
            # MYSTRAN writes 2 rows per QUAD4: top fiber (with EID) +
            # bottom fiber (no EID). We collapse to one entry per EID
            # keyed on von-Mises = max(top, bottom).
            mt = _F06_QUAD4_TOP_ROW.match(line)
            if mt:
                eid = int(mt.group(1))
                vals = [float(mt.group(i)) for i in range(2, 10)]
                # vals = [Z, Nx, Ny, Sxy, Angle, MaxP, MinP, VonMises]
                sx, sy, txy = vals[1], vals[2], vals[3]
                max_p, min_p, vm = vals[5], vals[6], vals[7]
                # Store top-fiber values; bottom-fiber row below will
                # update to max-of-top-bottom magnitude.
                subcases[cur_sc]['stresses'][eid] = {
                    'sx': sx, 'sy': sy, 'txy': txy,
                    'max_principal': max_p,
                    'min_principal': min_p,
                    'von_mises': vm,
                    '_top': vm, '_bot': None,
                }
                cur_stress_eid = eid
                n_quad4_stress_parsed += 1
                continue
            if cur_stress_eid is not None:
                mb = _F06_QUAD4_BOT_ROW.match(line)
                if mb:
                    vals = [float(mb.group(i)) for i in range(1, 9)]
                    bot_vm = vals[7]
                    cur = subcases[cur_sc]['stresses'].get(cur_stress_eid)
                    if cur is not None:
                        cur['_bot'] = bot_vm
                        # Use the worse of top / bottom for contouring.
                        try:
                            cur['von_mises'] = max(
                                abs(cur.get('_top') or 0.0),
                                abs(bot_vm))
                        except (TypeError, ValueError):
                            pass
                        # Track abs-max signed values for sx/sy/txy too.
                        try:
                            if abs(vals[1]) > abs(cur.get('sx', 0.0)):
                                cur['sx'] = vals[1]
                            if abs(vals[2]) > abs(cur.get('sy', 0.0)):
                                cur['sy'] = vals[2]
                            if abs(vals[3]) > abs(cur.get('txy', 0.0)):
                                cur['txy'] = vals[3]
                            if abs(vals[5]) > abs(cur.get('max_principal',
                                                          0.0)):
                                cur['max_principal'] = vals[5]
                            if abs(vals[6]) > abs(cur.get('min_principal',
                                                          0.0)):
                                cur['min_principal'] = vals[6]
                        except (TypeError, ValueError):
                            pass
                    cur_stress_eid = None
                    continue
                # The line wasn't a recognized stress row; bump the
                # unmatched counter only for rows that look like
                # data (start with whitespace + a digit).
                if line.lstrip()[:1].isdigit() or line.startswith(' '):
                    if any(c.isdigit() for c in line):
                        n_quad4_stress_unmatched += 1
        elif cur_section == 'stress_tria3':
            # v5.3.0 patch 2 Q2: same shape as QUAD4 but with "Anywhere"
            # keyword on top row + "in elem" prefix on bottom row.
            mt = _F06_TRIA3_TOP_ROW.match(line)
            if mt:
                eid = int(mt.group(1))
                vals = [float(mt.group(i)) for i in range(2, 10)]
                # vals = [Z, Nx, Ny, Sxy, Angle, MaxP, MinP, VonMises]
                subcases[cur_sc]['stresses'][eid] = {
                    'sx': vals[1], 'sy': vals[2], 'txy': vals[3],
                    'max_principal': vals[5],
                    'min_principal': vals[6],
                    'von_mises': vals[7],
                    '_top': vals[7], '_bot': None,
                }
                cur_stress_eid = eid
                n_tria3_stress_parsed += 1
                continue
            if cur_stress_eid is not None:
                mb = _F06_TRIA3_BOT_ROW.match(line)
                if mb:
                    vals = [float(mb.group(i)) for i in range(1, 9)]
                    bot_vm = vals[7]
                    cur = subcases[cur_sc]['stresses'].get(cur_stress_eid)
                    if cur is not None:
                        cur['_bot'] = bot_vm
                        try:
                            cur['von_mises'] = max(
                                abs(cur.get('_top') or 0.0),
                                abs(bot_vm))
                            for col_idx, key in (
                                    (1, 'sx'), (2, 'sy'), (3, 'txy'),
                                    (5, 'max_principal'),
                                    (6, 'min_principal')):
                                if abs(vals[col_idx]) > abs(
                                        cur.get(key, 0.0)):
                                    cur[key] = vals[col_idx]
                        except (TypeError, ValueError):
                            pass
                    cur_stress_eid = None
                    continue
        elif cur_section == 'stress_bar':
            # v5.3.0 patch 2 Q3: CBAR stress -- End-A row carries the
            # EID + 4 fiber stresses + axial + max/min + (optional)
            # margin of safety. End-B row drops the axial column.
            ma = _F06_BAR_END_A_ROW.match(line)
            if ma:
                eid = int(ma.group(1))
                sa1, sa2, sa3, sa4 = (
                    float(ma.group(2)), float(ma.group(3)),
                    float(ma.group(4)), float(ma.group(5)))
                axial = float(ma.group(6))
                sa_max = float(ma.group(7))
                sa_min = float(ma.group(8))
                subcases[cur_sc]['stresses'][eid] = {
                    'sa1': sa1, 'sa2': sa2, 'sa3': sa3, 'sa4': sa4,
                    'axial': axial,
                    'sa_max': sa_max, 'sa_min': sa_min,
                    'sb_max': sa_max, 'sb_min': sa_min,  # filled by end B
                    'von_mises': max(abs(sa_max), abs(sa_min), abs(axial)),
                }
                cur_stress_eid = eid
                n_bar_stress_parsed += 1
                continue
            if cur_stress_eid is not None:
                mb = _F06_BAR_END_B_ROW.match(line)
                if mb:
                    sb1, sb2, sb3, sb4 = (
                        float(mb.group(1)), float(mb.group(2)),
                        float(mb.group(3)), float(mb.group(4)))
                    sb_max = float(mb.group(5))
                    sb_min = float(mb.group(6))
                    cur = subcases[cur_sc]['stresses'].get(cur_stress_eid)
                    if cur is not None:
                        cur.update({
                            'sb1': sb1, 'sb2': sb2, 'sb3': sb3,
                            'sb4': sb4,
                            'sb_max': sb_max, 'sb_min': sb_min,
                        })
                        cur['von_mises'] = max(
                            abs(cur.get('sa_max', 0.0)),
                            abs(cur.get('sa_min', 0.0)),
                            abs(sb_max), abs(sb_min),
                            abs(cur.get('axial', 0.0)))
                    cur_stress_eid = None
                    continue
        elif cur_section == 'forces_bar':
            # v5.3.0 patch 2 Q4: CBAR engineering forces (one row).
            mf = _F06_BAR_FORCES_ROW.match(line)
            if mf:
                eid = int(mf.group(1))
                subcases[cur_sc]['forces'][eid] = {
                    'm1_a':  float(mf.group(2)),
                    'm2_a':  float(mf.group(3)),
                    'm1_b':  float(mf.group(4)),
                    'm2_b':  float(mf.group(5)),
                    'v1':    float(mf.group(6)),
                    'v2':    float(mf.group(7)),
                    'axial': float(mf.group(8)),
                    'torque': float(mf.group(9)),
                }
                n_bar_forces_parsed += 1
                continue

    # Drop subcases that contain no real data -- the parser would
    # otherwise return a "phantom" subcase 1 just because it saw a
    # SUBCASE header in the case-control echo at the top of the F06.
    # Callers use the emptiness check to distinguish "MYSTRAN ran +
    # produced output" from "MYSTRAN echoed the deck but bailed".
    real = {sid: sc for sid, sc in subcases.items()
            if sc['displacements'] or sc['eigenvalues']
            or sc['eigenvectors'] or sc['spc_forces']
            or sc['stresses'] or sc['forces']}
    try:
        from node_runner.profiling import perf_event
        perf_event('mystran', 'f06.parsed',
                   n_quad4_stress=n_quad4_stress_parsed,
                   n_quad4_unmatched=n_quad4_stress_unmatched,
                   n_tria3_stress=n_tria3_stress_parsed,
                   n_bar_stress=n_bar_stress_parsed,
                   n_bar_forces=n_bar_forces_parsed)
    except Exception:
        pass
    return {'filepath': str(f06_path), 'subcases': real}


def _empty_subcase() -> dict:
    # v5.3.0 patch 2 Q1: add `strains` and `forces` slots so the F06
    # bundle shape matches what the OP2 path already produces. Both stay
    # empty when the F06 doesn't carry the data (MYSTRAN doesn't write
    # strain to F06 at all; element-engineering forces are only present
    # for CBAR in current MYSTRAN builds). The Output Vector dropdown
    # gates each category on the dict being non-empty so empty slots
    # are invisible to the user.
    return {
        'displacements': {},
        'stresses': {},
        'strains': {},
        'forces': {},
        'eigenvectors': [],
        'spc_forces': {},
        'eigenvalues': [],
        'frequencies': [],
    }
