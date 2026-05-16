"""v5.2.0 item 43: contour data-conversion modes.

Common tools expose a "Data Conversion" menu on its contour controls so the
user can choose how element-centric output maps to vertex / cell
scalars (and vice versa). Node Runner v5.2.0 adopts the four most
useful modes:

* ``average``      -- nodal averaging (incident cells weighted equally).
                      The default; smoothest contour.
* ``no_avg``       -- per-element values, no smoothing. Cells keep their
                      own scalar; reveals discontinuities at element
                      boundaries.
* ``max_node``     -- when multiple elements share a node, take the max
                      of their values. Conservative-stress view.
* ``min_node``     -- like max_node but with min. Lowest-stress view.

This module is intentionally pure-Python (numpy only) so it is
testable without a live VTK plotter.

Two entry points::

    convert_nodal_data(nodal_values, node_ids, conversion,
                       cell_node_ids, cell_count) -> (point, cell)
    convert_element_data(elem_values, eids, conversion,
                         cell_eids, cell_node_ids, n_points) -> (point, cell)

``nodal_values`` and ``elem_values`` are 1-D numpy arrays already
projected to the chosen component (e.g. T3 magnitude or von_mises).
``cell_node_ids`` is a list of arrays giving the node-IDs for each
visible cell. ``cell_eids`` is the array of element-IDs for each cell.

Returns a ``(point_scalars, cell_scalars)`` tuple where exactly ONE
member is non-None depending on the conversion mode:

* ``average``  -> nodal source: ``point_scalars`` set.
                  element source: ``point_scalars`` set (computed by
                  averaging incident cells per node).
* ``no_avg``   -> always returns ``cell_scalars`` set.
* ``max_node`` -> nodal source: cell_scalars set (max of corner nodes).
                  element source: point_scalars set (max across
                  incident cells per node).
* ``min_node`` -> same as max_node but min.

The caller passes the result to PyVista as either
``add_mesh(scalars=point_scalars)`` (associated with grid.points) or
``add_mesh(scalars=cell_scalars)`` (associated with grid.cell_data).
"""

from __future__ import annotations

from typing import Optional

import numpy as np


CONVERSION_MODES = ('average', 'no_avg', 'max_node', 'min_node')


def convert_nodal_data(
    nodal_values: np.ndarray,
    node_ids_sorted: np.ndarray,
    conversion: str,
    cell_node_ids: list[np.ndarray] | None = None,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Convert a nodal scalar array to point or cell scalars per mode.

    Parameters
    ----------
    nodal_values
        1-D numpy array of length n_points carrying the scalar value
        at each grid point.
    node_ids_sorted
        Companion array of node-IDs (same length as ``nodal_values``)
        in grid-point order.
    conversion
        One of ``CONVERSION_MODES``.
    cell_node_ids
        List of length ``n_cells``; element ``i`` is a 1-D array of
        the node-IDs that make up cell ``i``. Only needed for
        ``no_avg`` / ``max_node`` / ``min_node`` modes.

    Returns
    -------
    (point_scalars, cell_scalars)
        Exactly one of the two is non-None per the mode.
    """
    if conversion not in CONVERSION_MODES:
        raise ValueError(f"Unknown conversion: {conversion}")
    nodal = np.asarray(nodal_values, dtype=float)
    if conversion == 'average':
        return nodal, None
    # The other three modes need cell connectivity to project to cells.
    if cell_node_ids is None:
        # Without connectivity we can't compute cell-level results;
        # fall back to the average path so the renderer at least
        # produces something visible.
        return nodal, None
    nid_to_idx = {int(n): i for i, n in enumerate(node_ids_sorted)}
    n_cells = len(cell_node_ids)
    if conversion == 'no_avg' or conversion in ('max_node', 'min_node'):
        # All three project to cell scalars first. max_node / min_node
        # for a NODAL source means "per cell, take the max/min of the
        # corner node values" -- which is then displayed flat-shaded.
        cell_scalars = np.zeros(n_cells, dtype=float)
        for ci, nids in enumerate(cell_node_ids):
            corner_vals = _gather(nodal, nid_to_idx, nids)
            if corner_vals.size == 0:
                continue
            if conversion == 'no_avg':
                cell_scalars[ci] = float(np.mean(corner_vals))
            elif conversion == 'max_node':
                cell_scalars[ci] = float(np.max(corner_vals))
            else:  # min_node
                cell_scalars[ci] = float(np.min(corner_vals))
        return None, cell_scalars
    return nodal, None


def convert_element_data(
    elem_values_by_eid: dict[int, float],
    conversion: str,
    cell_eids: np.ndarray,
    cell_node_ids: list[np.ndarray] | None,
    node_ids_sorted: np.ndarray | None,
    n_points: int,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Convert an element-centred scalar dict to point or cell scalars.

    Parameters
    ----------
    elem_values_by_eid
        Dict ``{eid: value}`` covering elements with results.
    conversion
        One of ``CONVERSION_MODES``.
    cell_eids
        Array of element-IDs in grid-cell order (one entry per visible
        cell). Length == number of cells.
    cell_node_ids
        See ``convert_nodal_data``. Required for ``average`` /
        ``max_node`` / ``min_node`` (the nodal-projection paths).
    node_ids_sorted
        Required for the nodal-projection paths.
    n_points
        Length of the point-scalars array to allocate when needed.

    Returns
    -------
    (point_scalars, cell_scalars)
    """
    if conversion not in CONVERSION_MODES:
        raise ValueError(f"Unknown conversion: {conversion}")
    cell_eids = np.asarray(cell_eids, dtype=np.int64)
    n_cells = cell_eids.size
    cell_scalars = np.zeros(n_cells, dtype=float)
    for i, eid in enumerate(cell_eids):
        v = elem_values_by_eid.get(int(eid))
        cell_scalars[i] = float(v) if v is not None else 0.0

    if conversion == 'no_avg':
        return None, cell_scalars

    # Other three modes project cell values to nodes.
    if cell_node_ids is None or node_ids_sorted is None:
        return None, cell_scalars   # fall back to centroid render

    nid_to_idx = {int(n): i for i, n in enumerate(node_ids_sorted)}
    sums = np.zeros(n_points, dtype=float)
    counts = np.zeros(n_points, dtype=np.int64)
    maxs = np.full(n_points, -np.inf, dtype=float)
    mins = np.full(n_points, np.inf, dtype=float)
    for ci, nids in enumerate(cell_node_ids):
        v = cell_scalars[ci]
        for nid in nids:
            pi = nid_to_idx.get(int(nid))
            if pi is None or pi < 0 or pi >= n_points:
                continue
            sums[pi] += v
            counts[pi] += 1
            if v > maxs[pi]:
                maxs[pi] = v
            if v < mins[pi]:
                mins[pi] = v
    if conversion == 'average':
        out = np.where(counts > 0, sums / np.maximum(counts, 1), 0.0)
    elif conversion == 'max_node':
        out = np.where(counts > 0, maxs, 0.0)
    else:  # min_node
        out = np.where(counts > 0, mins, 0.0)
    return out, None


def _gather(nodal: np.ndarray, nid_to_idx: dict, nids: np.ndarray) -> np.ndarray:
    """Pull the values at the given node IDs (silently dropping unknowns)."""
    out = []
    for nid in nids:
        pi = nid_to_idx.get(int(nid))
        if pi is not None and 0 <= pi < nodal.size:
            out.append(nodal[pi])
    return np.asarray(out, dtype=float)


__all__ = [
    'CONVERSION_MODES',
    'convert_nodal_data',
    'convert_element_data',
]
