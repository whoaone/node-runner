"""Scope an in-memory pyNastran BDF to an AnalysisSet selection.

v5.1.0 item 26: when the user runs/exports a specific AnalysisSet, we
want the resulting deck to contain only the entities that set targets:

- if ``group_target`` is set, only the nodes/elements that belong to
  that group (plus the properties / materials / coord systems they
  depend on);
- only the LOAD / SPC SIDs in the set's ``load_sids`` / ``spc_sids``
  white-lists (empty list = "include everything").

The scoping is implemented as a **shallow copy** of the BDF object:
we never mutate ``self.current_generator.model``. The copy is then
passed to the existing ``_write_bdf`` translator chain. Catchers in
:mod:`node_runner.profiling` log the before/after counts so the user
can verify what got included.

The function below is pure (no Qt deps) so it's straightforward to
test in isolation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from node_runner.dialogs.analysis import AnalysisSet

try:
    from node_runner.profiling import perf_event
except Exception:  # pragma: no cover
    def perf_event(*_a, **_kw):
        pass


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


@dataclass
class ScopeReport:
    """Summary of what :func:`scope_model_to_analysis_set` did."""
    n_nodes_kept: int = 0
    n_nodes_dropped: int = 0
    n_elements_kept: int = 0
    n_elements_dropped: int = 0
    n_properties_kept: int = 0
    n_properties_dropped: int = 0
    n_materials_kept: int = 0
    n_materials_dropped: int = 0
    n_coords_kept: int = 0
    n_coords_dropped: int = 0
    n_loads_kept: int = 0
    n_loads_dropped: int = 0
    n_spcs_kept: int = 0
    n_spcs_dropped: int = 0
    n_subcases_kept: int = 0
    n_subcases_dropped: int = 0
    group_target: Optional[str] = None
    notes: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def scope_model_to_analysis_set(
        model,
        analysis_set,
        groups: Optional[dict] = None,
        ):
    """Return ``(scoped_model, ScopeReport)``.

    ``scoped_model`` is a shallow copy of ``model`` with non-matching
    entities removed; ``model`` itself is untouched.

    When ``analysis_set.group_target`` is set, the matching entry in
    ``groups`` (dict {name: dict-with-nodes/elements/properties/
    materials/coords}) is used to define the geometric scope. If the
    group is missing or empty, an empty scope is returned and a note
    is recorded.

    LOAD / SPC SIDs are filtered by ``analysis_set.load_sids`` /
    ``spc_sids`` white-lists when those are non-empty (empty lists
    mean "include all").
    """
    report = ScopeReport(group_target=getattr(analysis_set, 'group_target',
                                              None))

    # Start with a shallow copy. pyNastran's BDF object isn't safely
    # deepcopy-friendly across all card classes, but a shallow copy
    # is enough for our purposes: we mutate the top-level dicts
    # (.nodes, .elements, .loads, .spcs, .properties, .materials,
    # .coords) on the copy, leaving the original untouched.
    scoped = copy.copy(model)
    for attr in ('nodes', 'elements', 'rigid_elements', 'properties',
                 'materials', 'thermal_materials', 'loads', 'spcs',
                 'coords', 'mpcs'):
        original = getattr(model, attr, None)
        if original is None:
            continue
        try:
            setattr(scoped, attr, dict(original))
        except Exception:
            # Some pyNastran versions present these as OrderedDicts /
            # CardDicts; copy through their public API where possible.
            setattr(scoped, attr, type(original)(original))

    group_target = getattr(analysis_set, 'group_target', None)
    load_sids = list(getattr(analysis_set, 'load_sids', []) or [])
    spc_sids = list(getattr(analysis_set, 'spc_sids', []) or [])

    # ---- 1. Geometric scope (group) ----
    if group_target:
        if not groups or group_target not in groups:
            report.notes.append(
                f"group_target '{group_target}' not found in groups dict; "
                f"falling back to full-model scope")
        else:
            _apply_group_scope(scoped, model, groups[group_target], report)

    # No group_target: keep everything geometric, just count for the report.
    if not group_target or group_target not in (groups or {}):
        report.n_nodes_kept = len(getattr(scoped, 'nodes', {}) or {})
        report.n_elements_kept = len(getattr(scoped, 'elements', {}) or {})
        report.n_properties_kept = len(
            getattr(scoped, 'properties', {}) or {})
        report.n_materials_kept = len(getattr(scoped, 'materials', {}) or {})
        report.n_coords_kept = len(getattr(scoped, 'coords', {}) or {})

    # ---- 2. LOAD SID white-list ----
    if load_sids:
        loads = getattr(scoped, 'loads', {}) or {}
        before = len(loads)
        keep = {sid: cards for sid, cards in loads.items()
                if int(sid) in set(int(x) for x in load_sids)}
        scoped.loads = keep
        report.n_loads_kept = len(keep)
        report.n_loads_dropped = before - len(keep)
    else:
        report.n_loads_kept = len(getattr(scoped, 'loads', {}) or {})

    # ---- 3. SPC SID white-list ----
    if spc_sids:
        spcs = getattr(scoped, 'spcs', {}) or {}
        before = len(spcs)
        keep = {sid: cards for sid, cards in spcs.items()
                if int(sid) in set(int(x) for x in spc_sids)}
        scoped.spcs = keep
        report.n_spcs_kept = len(keep)
        report.n_spcs_dropped = before - len(keep)
    else:
        report.n_spcs_kept = len(getattr(scoped, 'spcs', {}) or {})

    # ---- 4. Subcase filtering ----
    # Drop subcases that reference a SID we've removed.
    kept_loads = set(int(sid) for sid in (getattr(scoped, 'loads', {}) or {}))
    kept_spcs = set(int(sid) for sid in (getattr(scoped, 'spcs', {}) or {}))
    subcases = list(getattr(analysis_set, 'subcases', []) or [])
    kept_sc = []
    for sc in subcases:
        load_sid = sc.get('load_sid') if isinstance(sc, dict) else None
        spc_sid = sc.get('spc_sid') if isinstance(sc, dict) else None
        load_ok = (load_sid is None or int(load_sid) in kept_loads
                   or not getattr(scoped, 'loads', {}))
        spc_ok = (spc_sid is None or int(spc_sid) in kept_spcs
                  or not getattr(scoped, 'spcs', {}))
        if load_ok and spc_ok:
            kept_sc.append(sc)
    report.n_subcases_kept = len(kept_sc)
    report.n_subcases_dropped = len(subcases) - len(kept_sc)

    perf_event('export', 'scoped',
               group_target=group_target or '',
               n_nodes_kept=report.n_nodes_kept,
               n_elements_kept=report.n_elements_kept,
               n_loads_kept=report.n_loads_kept,
               n_loads_dropped=report.n_loads_dropped,
               n_spcs_kept=report.n_spcs_kept,
               n_spcs_dropped=report.n_spcs_dropped,
               n_subcases_kept=report.n_subcases_kept,
               n_subcases_dropped=report.n_subcases_dropped)
    return scoped, report


# ---------------------------------------------------------------------------
# Internal: geometric (group) scoping
# ---------------------------------------------------------------------------

def _apply_group_scope(scoped, original_model, group_data, report):
    """Filter scoped.nodes / .elements / .properties / .materials /
    .coords to entities in ``group_data`` plus auto-collected
    dependencies (so the resulting BDF is self-consistent).

    Mutates ``scoped`` in-place; reads ``original_model`` for the
    auto-collect lookups (element->PID, property->MID, node->CP, etc.).
    """
    group_nodes = set(int(n) for n in (group_data.get('nodes') or []))
    group_elements = set(int(e) for e in (group_data.get('elements') or []))
    group_properties = set(
        int(p) for p in (group_data.get('properties') or []))
    group_materials = set(int(m) for m in (group_data.get('materials') or []))
    group_coords = set(int(c) for c in (group_data.get('coords') or []))

    # ---- Elements ----
    elements = getattr(scoped, 'elements', {}) or {}
    rigid_elements = getattr(scoped, 'rigid_elements', {}) or {}
    if group_elements:
        before_e = len(elements)
        before_r = len(rigid_elements)
        scoped.elements = {eid: e for eid, e in elements.items()
                           if int(eid) in group_elements}
        scoped.rigid_elements = {eid: e for eid, e in rigid_elements.items()
                                 if int(eid) in group_elements}
        report.n_elements_kept = (len(scoped.elements)
                                  + len(scoped.rigid_elements))
        report.n_elements_dropped = (
            before_e - len(scoped.elements)
            + before_r - len(scoped.rigid_elements))
    else:
        # No element list -- drop all elements (group is node-only).
        report.n_elements_dropped = len(elements) + len(rigid_elements)
        scoped.elements = {}
        scoped.rigid_elements = {}
        report.n_elements_kept = 0

    # ---- Auto-collect dependencies from the kept elements ----
    auto_pids: set = set()
    auto_nodes: set = set()
    auto_coords: set = set()
    all_kept_elements = {**scoped.elements, **scoped.rigid_elements}
    for eid, elem in all_kept_elements.items():
        pid = getattr(elem, 'pid', None)
        if pid:
            auto_pids.add(int(pid))
        cid = getattr(elem, 'cid', None)
        if cid:
            auto_coords.add(int(cid))
        for n in (getattr(elem, 'nodes', None) or []):
            if n:
                auto_nodes.add(int(n))

    keep_pids = group_properties | auto_pids

    # ---- Properties ----
    properties = getattr(scoped, 'properties', {}) or {}
    before_p = len(properties)
    if keep_pids:
        scoped.properties = {pid: p for pid, p in properties.items()
                             if int(pid) in keep_pids}
    else:
        scoped.properties = {}
    report.n_properties_kept = len(scoped.properties)
    report.n_properties_dropped = before_p - len(scoped.properties)

    # ---- Materials (auto-collected from kept properties) ----
    auto_mids: set = set()
    for pid, prop in scoped.properties.items():
        for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
            m_id = getattr(prop, attr, None)
            if m_id:
                auto_mids.add(int(m_id))
        ply_mids = getattr(prop, 'mids', None)
        if ply_mids:
            for m_id in ply_mids:
                if m_id:
                    auto_mids.add(int(m_id))
    keep_mids = group_materials | auto_mids

    materials = getattr(scoped, 'materials', {}) or {}
    before_m = len(materials)
    if keep_mids:
        scoped.materials = {mid: m for mid, m in materials.items()
                            if int(mid) in keep_mids}
    else:
        scoped.materials = {}
    report.n_materials_kept = len(scoped.materials)
    report.n_materials_dropped = before_m - len(scoped.materials)

    # ---- Nodes ----
    nodes_src = getattr(scoped, 'nodes', {}) or {}
    keep_nids = group_nodes | auto_nodes
    before_n = len(nodes_src)
    if keep_nids:
        scoped.nodes = {nid: n for nid, n in nodes_src.items()
                        if int(nid) in keep_nids}
    else:
        scoped.nodes = {}
    report.n_nodes_kept = len(scoped.nodes)
    report.n_nodes_dropped = before_n - len(scoped.nodes)

    # ---- Coord systems (auto-collected from node CPs + element CIDs) ----
    for nid, node in scoped.nodes.items():
        cp = getattr(node, 'cp', None)
        if cp:
            auto_coords.add(int(cp))
    keep_cids = group_coords | auto_coords
    coords = getattr(scoped, 'coords', {}) or {}
    before_c = len(coords)
    if keep_cids:
        scoped.coords = {cid: c for cid, c in coords.items()
                         if int(cid) in keep_cids or int(cid) == 0}
    # Always keep the basic coord system (CID 0) if present.
    if 0 in coords and 0 not in scoped.coords:
        scoped.coords[0] = coords[0]
    report.n_coords_kept = len(scoped.coords)
    report.n_coords_dropped = before_c - len(scoped.coords)
