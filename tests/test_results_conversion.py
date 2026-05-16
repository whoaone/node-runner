"""v5.2.0 item 43: unit tests for results_conversion.

A small 4-element shell strip (3 nodes wide) gives us enough
connectivity to exercise nodal averaging and Max/Min-at-Node paths
without dragging in a real OP2 fixture.

Strip layout (3 rows of nodes, 4 quads sharing edges):

    n4 -- n5 -- n6
    |  e3 |  e4 |
    n2 -- n3 -- ?         (wait -- this is harder. Use a simpler 4-node 3-elem strip)

Actually use a simple 6-node, 4-element strip on a 3x2 grid:

    n4 -- n5 -- n6
    |  e3 |  e4 |
    n1 -- n2 -- n3

Wait, n4-n5-n6 top, n1-n2-n3 bottom. e3 = (n1,n2,n5,n4) quad; e4 =
(n2,n3,n6,n5) quad. Only 2 quads. Good enough.

For tests we lean on hand-computed expected outputs.
"""

from __future__ import annotations

import numpy as np
import pytest

from node_runner.results_conversion import (
    CONVERSION_MODES,
    convert_nodal_data,
    convert_element_data,
)


# ---------- fixture data ----------

# 6 grid points: ids 1..6, sorted same as point order.
NODE_IDS = np.array([1, 2, 3, 4, 5, 6])

# Two quads: e10 uses nodes (1,2,5,4); e20 uses nodes (2,3,6,5).
CELL_EIDS = np.array([10, 20])
CELL_NODE_IDS = [
    np.array([1, 2, 5, 4]),
    np.array([2, 3, 6, 5]),
]


# ----------------------------------------------------------------------
# Nodal source tests
# ----------------------------------------------------------------------

class TestConvertNodalData:

    def test_average_passes_through(self):
        nodal = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        point, cell = convert_nodal_data(nodal, NODE_IDS, 'average',
                                         CELL_NODE_IDS)
        assert point is not None
        assert cell is None
        assert np.allclose(point, nodal)

    def test_no_avg_projects_to_cell_mean(self):
        nodal = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        point, cell = convert_nodal_data(nodal, NODE_IDS, 'no_avg',
                                         CELL_NODE_IDS)
        assert point is None
        assert cell is not None
        # e10 corner mean = (1+2+5+4)/4 = 3.0
        # e20 corner mean = (2+3+6+5)/4 = 4.0
        assert np.allclose(cell, [3.0, 4.0])

    def test_max_node_for_nodal_source(self):
        nodal = np.array([1.0, 10.0, 3.0, 4.0, 5.0, 6.0])
        point, cell = convert_nodal_data(nodal, NODE_IDS, 'max_node',
                                         CELL_NODE_IDS)
        assert point is None
        # e10 corner max = max(1,10,5,4) = 10
        # e20 corner max = max(10,3,6,5) = 10
        assert np.allclose(cell, [10.0, 10.0])

    def test_min_node_for_nodal_source(self):
        nodal = np.array([1.0, 10.0, 3.0, 4.0, 5.0, 6.0])
        point, cell = convert_nodal_data(nodal, NODE_IDS, 'min_node',
                                         CELL_NODE_IDS)
        assert point is None
        # e10 corner min = min(1,10,5,4) = 1
        # e20 corner min = min(10,3,6,5) = 3
        assert np.allclose(cell, [1.0, 3.0])

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError):
            convert_nodal_data(np.zeros(6), NODE_IDS, 'bogus', None)


# ----------------------------------------------------------------------
# Element source tests
# ----------------------------------------------------------------------

class TestConvertElementData:

    def test_no_avg_passes_cell_values(self):
        elem_vals = {10: 100.0, 20: 200.0}
        point, cell = convert_element_data(
            elem_vals, 'no_avg', CELL_EIDS,
            CELL_NODE_IDS, NODE_IDS, n_points=6)
        assert point is None
        assert np.allclose(cell, [100.0, 200.0])

    def test_average_nodes_get_mean_of_incident_cells(self):
        elem_vals = {10: 100.0, 20: 200.0}
        point, cell = convert_element_data(
            elem_vals, 'average', CELL_EIDS,
            CELL_NODE_IDS, NODE_IDS, n_points=6)
        assert cell is None
        assert point is not None
        # n1: only e10 -> 100
        # n2: shared (e10,e20) -> avg 150
        # n3: only e20 -> 200
        # n4: only e10 -> 100
        # n5: shared (e10,e20) -> 150
        # n6: only e20 -> 200
        assert np.allclose(point, [100, 150, 200, 100, 150, 200])

    def test_max_node_takes_max_incident_cell(self):
        elem_vals = {10: 100.0, 20: 200.0}
        point, cell = convert_element_data(
            elem_vals, 'max_node', CELL_EIDS,
            CELL_NODE_IDS, NODE_IDS, n_points=6)
        assert cell is None
        # n2 and n5 see both elements -> max = 200
        assert np.allclose(point, [100, 200, 200, 100, 200, 200])

    def test_min_node_takes_min_incident_cell(self):
        elem_vals = {10: 100.0, 20: 200.0}
        point, cell = convert_element_data(
            elem_vals, 'min_node', CELL_EIDS,
            CELL_NODE_IDS, NODE_IDS, n_points=6)
        assert cell is None
        # n2 and n5 see both elements -> min = 100
        assert np.allclose(point, [100, 100, 200, 100, 100, 200])

    def test_missing_eid_yields_zero(self):
        elem_vals = {10: 100.0}     # e20 missing
        point, cell = convert_element_data(
            elem_vals, 'no_avg', CELL_EIDS,
            CELL_NODE_IDS, NODE_IDS, n_points=6)
        assert cell is not None
        # e20 -> 0.0 because not in dict
        assert np.allclose(cell, [100.0, 0.0])

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError):
            convert_element_data({}, 'bogus', CELL_EIDS,
                                 None, None, n_points=0)


def test_all_modes_constant():
    """Sanity: CONVERSION_MODES exposes exactly the four documented modes."""
    assert set(CONVERSION_MODES) == {
        'average', 'no_avg', 'max_node', 'min_node'}
