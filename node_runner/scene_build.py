"""Vectorized scene-build helpers + display-time LOD.

v3.2.0 replaces the per-element Python loop in MainWindow._update_viewer
with bulk numpy operations. On a 600k-element deck this drops the scene-
build wall clock from ~30 s (which made Windows mark the app "Not
Responding") to roughly 1 s. The bottleneck was always:

  1. node.get_position() called once per node - each call does a matrix
     multiply against the node's cp_ref coord system. For 1M nodes
     that's 1M Python method calls + 1M matmuls = many seconds.
  2. per-element type dispatch with [node_map[nid] for nid in elem.nodes]
     building per-type connectivity lists. For 1M elements that's
     many millions of dict lookups.

The vectorized path:
  * Builds node.xyz into one np.array (no method calls).
  * Groups nodes by cp; cp=0 nodes (almost all) get used directly,
    non-basic cp gets one batched matmul per coord.
  * Sorts elements by type once, then builds (n_elem, n_nodes_per_elem)
    numpy arrays per type with np.fromiter.
  * Uses np.searchsorted for node-ID -> grid-index lookup in O(N log N)
    instead of O(N) per element.

The display-LOD helper (apply_display_lod) is what we hand to
plotter.add_mesh; the underlying current_grid keeps every cell for
picking. Femap-style: above 500k elements we stride-sample shells and
extract the outer surface of solids.
"""

from __future__ import annotations

import time

import numpy as np
import pyvista as pv

try:  # v5.0.0 item 1: diagnostic timers are gated on NR_PROFILE=1.
    from node_runner.profiling import perf_event as _perf_event
except Exception:  # pragma: no cover
    def _perf_event(*_a, **_kw):
        pass

# Element-type metadata. Maps the elem.type string seen on the pyNastran
# card object to:
#   * 'kind'         : 'shell' | 'solid' | 'beam' | 'rod' | 'bush' | 'gap' | 'shear' | 'rigid'
#   * 'n_nodes'      : how many of elem.nodes are used
#   * 'vtk_cell'     : pv.CellType for the resulting VTK cell
ETYPE_INFO = {
    'CQUAD4':   {'kind': 'shell', 'n_nodes': 4, 'vtk_cell': pv.CellType.QUAD},
    'CMEMBRAN': {'kind': 'shell', 'n_nodes': 4, 'vtk_cell': pv.CellType.QUAD},
    'CTRIA3':   {'kind': 'shell', 'n_nodes': 3, 'vtk_cell': pv.CellType.TRIANGLE},
    'CHEXA':    {'kind': 'solid', 'n_nodes': 8, 'vtk_cell': pv.CellType.HEXAHEDRON},
    'CHEXA8':   {'kind': 'solid', 'n_nodes': 8, 'vtk_cell': pv.CellType.HEXAHEDRON},
    'CHEXA20':  {'kind': 'solid', 'n_nodes': 8, 'vtk_cell': pv.CellType.HEXAHEDRON},
    'CTETRA':   {'kind': 'solid', 'n_nodes': 4, 'vtk_cell': pv.CellType.TETRA},
    'CTETRA4':  {'kind': 'solid', 'n_nodes': 4, 'vtk_cell': pv.CellType.TETRA},
    'CTETRA10': {'kind': 'solid', 'n_nodes': 4, 'vtk_cell': pv.CellType.TETRA},
    'CPENTA':   {'kind': 'solid', 'n_nodes': 6, 'vtk_cell': pv.CellType.WEDGE},
    'CPENTA6':  {'kind': 'solid', 'n_nodes': 6, 'vtk_cell': pv.CellType.WEDGE},
    'CPENTA15': {'kind': 'solid', 'n_nodes': 6, 'vtk_cell': pv.CellType.WEDGE},
    'CSHEAR':   {'kind': 'shear', 'n_nodes': 4, 'vtk_cell': pv.CellType.QUAD},
    'CBEAM':    {'kind': 'beam',  'n_nodes': 2, 'vtk_cell': pv.CellType.LINE},
    'CBAR':     {'kind': 'beam',  'n_nodes': 2, 'vtk_cell': pv.CellType.LINE},
    'CROD':     {'kind': 'rod',   'n_nodes': 2, 'vtk_cell': pv.CellType.LINE},
    'CBUSH':    {'kind': 'bush',  'n_nodes': 2, 'vtk_cell': pv.CellType.LINE},
    'CGAP':     {'kind': 'gap',   'n_nodes': 2, 'vtk_cell': pv.CellType.LINE},
}

# kinds we render as shells (1 in the is_shell flag, used for surface
# normals / lighting in the plotter).
#
# v5.0.0 item 11c (shell-bucket split): the legacy single 'shell' kind
# mixed 4-node (CQUAD4 / CMEMBRAN) and 3-node (CTRIA3) cells in the
# same bucket. With max_n=4 and CTRIA3 having actual_n=3 the
# `(actual_n == max_n).all()` check failed and the entire bucket fell
# through to a Python row-by-row loop -- ~11 s on a 2.4M-cell deck.
# The split partitions shells by n_nodes so each sub-bucket is
# uniform-width and goes through the vectorized block-stack path.
# Both sub-kinds remain shells for is_shell=1 purposes (downstream
# lighting / extract_cells / slicing only reads the binary flag).
SHELL_KINDS = {'shell', 'shell_q4', 'shell_tri3', 'shear'}

# v3.3.0: tree-grouping maps. These used to be local dicts inside
# MainWindow._populate_tree (which iterated every element twice to
# build Counters). Hoisting them up means we can piggy-back the
# counting onto build_element_arrays_vectorized's existing single
# pass over all elements.
TREE_BY_TYPE_MAP = {
    'CBEAM': 'Beams', 'CBAR': 'Bars', 'CROD': 'Rods',
    'CBUSH': 'Bushes',
    'CQUAD4': 'Plates', 'CMEMBRAN': 'Plates', 'CTRIA3': 'Plates',
    'RBE2': 'Rigid', 'RBE3': 'Rigid',
    'CHEXA': 'Solids', 'CHEXA8': 'Solids', 'CHEXA20': 'Solids',
    'CTETRA': 'Solids', 'CTETRA4': 'Solids', 'CTETRA10': 'Solids',
    'CPENTA': 'Solids', 'CPENTA6': 'Solids', 'CPENTA15': 'Solids',
    'CSHEAR': 'Shear', 'CGAP': 'Gap',
}
TREE_SHAPE_MAP = {
    'Line': ('CBEAM', 'CBAR', 'CROD', 'CBUSH', 'CGAP'),
    'Quad': ('CQUAD4', 'CMEMBRAN', 'CSHEAR'),
    'Tria': ('CTRIA3',),
    'Rigid': ('RBE2', 'RBE3'),
    'Hex': ('CHEXA', 'CHEXA8', 'CHEXA20'),
    'Tet': ('CTETRA', 'CTETRA4', 'CTETRA10'),
    'Wedge': ('CPENTA', 'CPENTA6', 'CPENTA15'),
}
# Inverted for O(1) lookup in the element loop.
_TYPE_TO_SHAPE = {
    etype: shape for shape, etypes in TREE_SHAPE_MAP.items()
    for etype in etypes
}


def build_node_coords_vectorized(model, sorted_node_ids):
    """Return (N, 3) numpy array of global node positions.

    Equivalent to ``np.array([model.nodes[nid].get_position() for nid in
    sorted_node_ids])`` but ~50x faster because it avoids the per-node
    method call and re-uses each coord's transformation matrix.

    Requires ``model`` to have already been passed through
    ``NastranModelGenerator._finalize_for_viewer`` so every node has a
    ``cp_ref``. Falls back to ``node.xyz`` if cp_ref is missing or the
    coord system can't be resolved.
    """
    n = len(sorted_node_ids)
    if n == 0:
        return np.zeros((0, 3), dtype=np.float64)

    # Build a flat xyz array in node-id order.
    xyz = np.empty((n, 3), dtype=np.float64)
    cp_ids = np.empty(n, dtype=np.int64)
    for i, nid in enumerate(sorted_node_ids):
        node = model.nodes[nid]
        xyz[i, 0] = node.xyz[0]
        xyz[i, 1] = node.xyz[1]
        xyz[i, 2] = node.xyz[2]
        cp_val = getattr(node, 'cp', 0)
        try:
            cp_ids[i] = int(cp_val) if cp_val is not None else 0
        except (TypeError, ValueError):
            cp_ids[i] = 0

    # Fast path: every node uses the basic coord. xyz is already global.
    unique_cps = np.unique(cp_ids)
    if unique_cps.size == 1 and unique_cps[0] == 0:
        return xyz

    # General path: for each non-basic cp, batch-transform via
    #   global = origin + xyz_local @ axes_matrix.T
    # where axes_matrix has the coord system's i/j/k axes as columns.
    # We loop over unique cp values, not over nodes - typical decks have
    # 1-20 coord systems vs 1M nodes, so this scales beautifully.
    coords = getattr(model, 'coords', {}) or {}
    out = xyz.copy()
    for cp_int in unique_cps:
        if cp_int == 0:
            continue
        cs = coords.get(int(cp_int))
        if cs is None:
            continue
        try:
            origin = np.asarray(cs.origin, dtype=np.float64).reshape(3)
            i_hat = np.asarray(cs.i, dtype=np.float64).reshape(3)
            j_hat = np.asarray(cs.j, dtype=np.float64).reshape(3)
            k_hat = np.asarray(cs.k, dtype=np.float64).reshape(3)
        except Exception:
            continue
        # Axes matrix: columns i, j, k. Multiply local xyz by this
        # to rotate into the global frame, then add the origin.
        axes = np.column_stack([i_hat, j_hat, k_hat])
        mask = cp_ids == cp_int
        out[mask] = origin + xyz[mask] @ axes.T
    return out


def build_element_arrays_vectorized(
        model, sorted_node_ids, node_map, progress=None):
    """Build per-element-kind connectivity arrays for VTK construction.

    Returns a dict with:
      'cells'      : flat int64 ndarray ready for pv.UnstructuredGrid
      'cell_types' : ndarray of pv.CellType values, one per cell
      'eid'        : ndarray of element IDs (one per cell)
      'pid'        : ndarray of property IDs (one per cell)
      'etype'      : ndarray of dtype object with the elem.type string
      'is_shell'   : ndarray of 0/1 for is_shell flag
      'nodes_used' : set() of node IDs actually referenced (drives
                     free-node detection at the caller)

    ``progress``, if provided, is called as
    ``progress(done, total, kind)`` periodically so the dialog can
    advance during a long extract. The helper is purely numpy + a
    single dict-iterate, so it should be ~50x faster than the per-
    element Python loop in the old ``_update_viewer``.
    """
    nodes_used = set()
    elem_dicts = (
        getattr(model, 'elements', {}) or {},
        getattr(model, 'rigid_elements', {}) or {},
    )
    n_total = sum(len(d) for d in elem_dicts)
    # v3.3.0: counters for the tree's By-Type / By-Shape groups,
    # accumulated in the same pass as the cell-building so we don't
    # have to iterate all elements a second time later in
    # MainWindow._populate_tree.
    by_type_counts: dict = {}
    by_shape_counts: dict = {}

    if n_total == 0:
        return {
            'cells': np.array([], dtype=np.int64),
            'cell_types': np.array([], dtype=np.uint8),
            'eid': np.array([], dtype=np.int64),
            'pid': np.array([], dtype=np.int64),
            'etype': np.array([], dtype=object),
            'is_shell': np.array([], dtype=np.int8),
            'nodes_used': nodes_used,
            'by_type_counts': by_type_counts,
            'by_shape_counts': by_shape_counts,
        }

    # Pre-sort element IDs by kind. RBE2/RBE3 are rigid; we record
    # their node usage but don't render them as cells (MainWindow has
    # a dedicated _create_rbe_actors path).
    #
    # v5.0.0 item 11c: 'shell' is split into uniform-width sub-buckets
    # so the vectorized block-stack path engages instead of the Python
    # row-by-row fallback. shell_q4 holds 4-node (CQUAD4, CMEMBRAN);
    # shell_tri3 holds 3-node (CTRIA3). 'shell' remains in the dicts
    # as a back-compat placeholder (always empty, always 0); a
    # synthetic 'shell' aggregate is emitted in the perf event below.
    buckets_by_kind = {
        'shell': [], 'shell_q4': [], 'shell_tri3': [],
        'shear': [], 'solid': [], 'beam': [],
        'rod': [], 'bush': [], 'gap': [],
    }
    # Per-bucket parallel lists: eids, pids, etype, node-id rows (each
    # row is a tuple of length info.n_nodes).
    eid_by_kind = {k: [] for k in buckets_by_kind}
    pid_by_kind = {k: [] for k in buckets_by_kind}
    etype_by_kind = {k: [] for k in buckets_by_kind}
    nodes_by_kind = {k: [] for k in buckets_by_kind}
    nodes_per_kind = {k: 0 for k in buckets_by_kind}

    processed = 0
    progress_every = max(50_000, n_total // 50)

    # v5.0.0 item 1: per-kind wall-time + element-count attribution for the
    # element-dict iteration loop. The dominant kind on production decks is the lever
    # for the v5.1 real fix.
    _iter_t0 = time.perf_counter()
    kind_walltime: dict = {k: 0.0 for k in buckets_by_kind}
    kind_count: dict = {k: 0 for k in buckets_by_kind}
    kind_count['rigid_or_other'] = 0
    kind_walltime['rigid_or_other'] = 0.0

    for elem_dict in elem_dicts:
        for eid, elem in elem_dict.items():
            _t_elem = time.perf_counter()
            etype = getattr(elem, 'type', None)
            if etype is None:
                continue
            # v3.3.0: bump tree-group counters in the same pass.
            group = TREE_BY_TYPE_MAP.get(etype, 'Other')
            by_type_counts[group] = by_type_counts.get(group, 0) + 1
            shape = _TYPE_TO_SHAPE.get(etype, 'Other')
            by_shape_counts[shape] = by_shape_counts.get(shape, 0) + 1
            if etype in ('RBE2', 'RBE3'):
                # Record node usage for free-node detection.
                try:
                    nodes_used.update(elem.independent_nodes)
                    nodes_used.update(elem.dependent_nodes)
                except (AttributeError, TypeError):
                    pass
                processed += 1
                kind_count['rigid_or_other'] += 1
                kind_walltime['rigid_or_other'] += (
                    time.perf_counter() - _t_elem)
                continue
            info = ETYPE_INFO.get(etype)
            if info is None:
                kind_count['rigid_or_other'] += 1
                kind_walltime['rigid_or_other'] += (
                    time.perf_counter() - _t_elem)
                continue
            n_nodes = info['n_nodes']
            try:
                raw_nodes = elem.nodes[:n_nodes]
            except (AttributeError, TypeError):
                continue
            # Filter falsy entries (None / 0) which pyNastran uses for
            # collapsed CHEXA/CTETRA forms - those reduce n_nodes
            # effectively.
            row = tuple(n for n in raw_nodes if n)
            if len(row) < 2:
                continue
            kind = info['kind']
            # v5.0.0 item 11c: route shell elements into uniform-width
            # sub-buckets so the vectorized block-stack composes in
            # numpy rather than falling through to the Python row-by-row
            # loop. Solid is left alone -- collapsed CHEXA/CTETRA forms
            # do trigger the ragged path occasionally but the row count
            # is bounded (~6k on production decks) so the cost is
            # negligible. Only shells were a problem (2.4M rows).
            if kind == 'shell':
                if len(row) == 4:
                    kind = 'shell_q4'
                elif len(row) == 3:
                    kind = 'shell_tri3'
                # else: leave as 'shell' (unexpected, but defensive)
            eid_by_kind[kind].append(eid)
            pid_by_kind[kind].append(getattr(elem, 'pid', 0))
            etype_by_kind[kind].append(etype)
            nodes_by_kind[kind].append(row)
            nodes_per_kind[kind] = max(nodes_per_kind[kind], len(row))
            nodes_used.update(row)

            processed += 1
            kind_count[kind] += 1
            kind_walltime[kind] += time.perf_counter() - _t_elem
            if progress is not None and processed % progress_every == 0:
                try:
                    progress(processed, n_total, kind)
                except Exception:
                    pass

    _iter_total_s = time.perf_counter() - _iter_t0
    try:
        # v5.0.0 item 11c: synthesize a 'shell' aggregate (shell_q4 +
        # shell_tri3) in the kind_walltime / kind_count dicts so
        # downstream log-grep queries comparing old vs new builds
        # continue to find a 'shell' entry. The legacy 'shell' bucket
        # itself is always empty after the split.
        shell_total_wall = (kind_walltime.get('shell_q4', 0.0)
                            + kind_walltime.get('shell_tri3', 0.0))
        shell_total_count = (kind_count.get('shell_q4', 0)
                             + kind_count.get('shell_tri3', 0))
        _perf_event(
            'scene_build', 'kind_iteration_walltime',
            wall_s=round(_iter_total_s, 3),
            n_total=n_total,
            kind_walltime={k: round(v, 3) for k, v in kind_walltime.items()},
            kind_count=kind_count,
            shell_total_wall_s=round(shell_total_wall, 3),
            shell_total_count=shell_total_count,
        )
    except Exception:
        pass

    # Build VTK cells per kind. The cells array format is:
    #   [n0, idx0_0, idx0_1, ..., n1, idx1_0, idx1_1, ...]
    # where ni is the node count for cell i. We use np.searchsorted on
    # sorted_node_ids for fast id-to-index lookup. node_map is also
    # accepted (dict) for back-compat with existing tests.
    cells_all = []
    cell_types_all = []
    eids_all = []
    pids_all = []
    etype_all = []
    is_shell_all = []

    sorted_nid_arr = np.asarray(sorted_node_ids, dtype=np.int64)

    # v5.0.0 item 11c: 'shell' kept in the iteration tuple for safety
    # (always empty after the split routes elements into shell_q4 /
    # shell_tri3) so any defensive lookup of nodes_by_kind['shell'] in
    # downstream code returns an empty list rather than KeyError.
    for kind in ('shell', 'shell_q4', 'shell_tri3', 'shear',
                 'solid', 'beam', 'rod', 'bush', 'gap'):
        rows = nodes_by_kind[kind]
        if not rows:
            continue
        _t_kind = time.perf_counter()
        # All rows in a kind should have the same length (n_nodes)
        # because each ETYPE_INFO entry pins it. Sanity-pad if a
        # collapsed solid produced a shorter row.
        max_n = max(len(r) for r in rows)
        # Build (n_rows, max_n) array of raw node IDs.
        raw = np.zeros((len(rows), max_n), dtype=np.int64)
        actual_n = np.empty(len(rows), dtype=np.int64)
        for i, r in enumerate(rows):
            raw[i, :len(r)] = r
            actual_n[i] = len(r)
        # Vectorized lookup: searchsorted finds the position of each
        # node ID in the sorted id array. We trust the caller that all
        # node IDs are present.
        flat_raw = raw.ravel()
        flat_idx = np.searchsorted(sorted_nid_arr, flat_raw).reshape(
            len(rows), max_n)
        uniform = bool((actual_n == max_n).all())
        # Compose the VTK cells block. For uniform-width rows this is
        # one column of n_nodes + n_nodes columns of indices.
        if uniform:
            block = np.empty((len(rows), max_n + 1), dtype=np.int64)
            block[:, 0] = max_n
            block[:, 1:] = flat_idx
            cells_all.append(block.ravel())
        else:
            # Mixed-width within a kind (e.g. degenerate CHEXAs).
            # Fall back to a Python loop just for these rows; rare.
            for i, r in enumerate(rows):
                cells_all.append(np.concatenate(
                    [[len(r)], flat_idx[i, :len(r)]]).astype(np.int64))
        # cell_types: pick the dominant cell type for this kind.
        info_proto = ETYPE_INFO[etype_by_kind[kind][0]]
        cell_types_all.append(
            np.full(len(rows), info_proto['vtk_cell'], dtype=np.uint8))
        eids_all.append(np.asarray(eid_by_kind[kind], dtype=np.int64))
        pids_all.append(np.asarray(pid_by_kind[kind], dtype=np.int64))
        etype_all.append(np.asarray(etype_by_kind[kind], dtype=object))
        is_shell_all.append(np.full(
            len(rows), 1 if kind in SHELL_KINDS else 0, dtype=np.int8))
        # v5.0.0 item 1: per-kind cell-composition timer so the v5.1
        # fix can target the dominant kind precisely.
        try:
            _perf_event(
                'scene_build', f'compose_{kind}',
                wall_s=round(time.perf_counter() - _t_kind, 3),
                n_rows=len(rows), max_n=int(max_n), uniform=uniform,
            )
        except Exception:
            pass

    if cells_all:
        cells = np.concatenate(cells_all)
        cell_types = np.concatenate(cell_types_all)
        eids = np.concatenate(eids_all)
        pids = np.concatenate(pids_all)
        etype = np.concatenate(etype_all)
        is_shell = np.concatenate(is_shell_all)
    else:
        cells = np.array([], dtype=np.int64)
        cell_types = np.array([], dtype=np.uint8)
        eids = np.array([], dtype=np.int64)
        pids = np.array([], dtype=np.int64)
        etype = np.array([], dtype=object)
        is_shell = np.array([], dtype=np.int8)

    return {
        'cells': cells,
        'cell_types': cell_types,
        'eid': eids,
        'pid': pids,
        'etype': etype,
        'is_shell': is_shell,
        'nodes_used': nodes_used,
        'by_type_counts': by_type_counts,
        'by_shape_counts': by_shape_counts,
    }


def apply_display_lod(grid, threshold=500_000, target_shells=200_000):
    """Return a (possibly decimated) grid suitable for fast rendering.

    The returned grid is what we hand to plotter.add_mesh. The original
    ``grid`` stays untouched and is what the picker / queries use.

    Rules (Femap-style):
      * If grid.n_cells <= threshold: return grid unchanged (no LOD).
      * Otherwise:
          - For shell/quad/tria cells, stride-sample so the count is
            <= target_shells.
          - For solid cells (HEXA/TETRA/WEDGE), apply
            vtkDataSetSurfaceFilter to keep only the outer surface
            (typically 25-50% the cell count, plus they're 2D triangles
            so much cheaper to render).
          - Beam/rod/bush/gap (LINE cells) and vertices are kept as-is.

    Returns ``(display_grid, info_dict)`` where info_dict has:
      'lod_active'        : bool
      'displayed_cells'   : int
      'total_cells'       : int
    """
    info = {
        'lod_active': False,
        'displayed_cells': grid.n_cells,
        'total_cells': grid.n_cells,
    }
    if grid.n_cells <= threshold:
        return grid, info

    # Pull the cell types out so we can mask by kind.
    try:
        cell_types = np.asarray(grid.celltypes)
    except Exception:
        return grid, info

    QUAD = int(pv.CellType.QUAD)
    TRI = int(pv.CellType.TRIANGLE)
    HEX = int(pv.CellType.HEXAHEDRON)
    TET = int(pv.CellType.TETRA)
    WEDGE = int(pv.CellType.WEDGE)
    LINE = int(pv.CellType.LINE)
    VERTEX = int(pv.CellType.VERTEX)

    shell_mask = (cell_types == QUAD) | (cell_types == TRI)
    solid_mask = (cell_types == HEX) | (cell_types == TET) | (cell_types == WEDGE)
    other_mask = ~(shell_mask | solid_mask)

    keep_cells = []

    # v4.0.1 (Stage 1): instrument the two suspect VTK ops.
    from node_runner.profiling import perf_stage, perf_event

    # Stride-sample shells.
    shell_idxs = np.where(shell_mask)[0]
    if shell_idxs.size:
        with perf_stage('lod', 'shell_stride_sample',
                        n_shells=int(shell_idxs.size),
                        target=int(target_shells)):
            stride = max(1, int(np.ceil(shell_idxs.size / target_shells)))
            keep_cells.append(shell_idxs[::stride])

    # For solids: extract outer surface via a sub-grid.
    solid_idxs = np.where(solid_mask)[0]
    solid_surface = None
    if solid_idxs.size:
        try:
            with perf_stage('lod', 'solid_extract_cells',
                            n_solids=int(solid_idxs.size)):
                solid_grid = grid.extract_cells(solid_idxs)
            with perf_stage('lod', 'solid_extract_surface',
                            n_solids=int(solid_idxs.size)):
                solid_surface = solid_grid.extract_surface()
        except Exception as exc:
            perf_event('lod', 'solid_extract_failed', exc=str(exc)[:120])
            # Fallback: stride-sample the solids too.
            stride = max(1, int(np.ceil(solid_idxs.size / target_shells)))
            keep_cells.append(solid_idxs[::stride])

    # Keep all beams / lines / vertices (typically small count).
    other_idxs = np.where(other_mask)[0]
    if other_idxs.size:
        keep_cells.append(other_idxs)

    if keep_cells:
        keep_arr = np.concatenate(keep_cells)
        with perf_stage('lod', 'extract_keep_cells',
                        n_keep=int(keep_arr.size)):
            partial = grid.extract_cells(keep_arr)
    else:
        partial = None

    if solid_surface is not None and partial is not None:
        # Append the solid surface onto the cell-extracted partial.
        # Both should be UnstructuredGrids (or compatible).
        try:
            with perf_stage('lod', 'merge_partial_and_solid_surface'):
                display = partial.merge(solid_surface)
        except Exception:
            display = partial
    elif solid_surface is not None:
        display = solid_surface
    else:
        display = partial if partial is not None else grid

    info['lod_active'] = True
    info['displayed_cells'] = display.n_cells if display is not None else 0
    return display, info
