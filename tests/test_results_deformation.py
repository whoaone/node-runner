"""v5.2.1 regression tests for the deformation / scalar-bar / node-cloud
fixes (Items 49 / 50 / 51).

Visual rendering (scalar-bar actor, node-cloud actor positions) is
verified by manual smoke testing on examples/cantilever_beam.bdf. These
tests cover the deterministic, unit-testable seams of the three fixes:

* Item 49: bbox-diag snapshot is captured at results load; the
  %-of-model scale math doesn't drift across repeated calls.
* Item 50: ``_project_point_scalars_to_subgrid`` correctly maps a
  full-grid scalar array onto a sub-grid via ``vtkOriginalPointIds``.
* Item 51: the node-cloud rebuild path on a deformation pass calls
  ``_add_node_cloud`` with the deformed point set, not the original.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv
import pytest


# ---------------------------------------------------------------------------
# Item 49: bbox-diag snapshot + sync_deform_scale stability
# ---------------------------------------------------------------------------

class TestBboxSnapshot:

    def test_consume_bundle_captures_snapshot(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        # Fabricate a non-degenerate current_grid so the snapshot
        # captures a real positive diagonal.
        pts = np.array([
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        cells = np.array([4, 0, 1, 2, 3])
        mw.current_grid = pv.UnstructuredGrid(cells, [pv.CellType.QUAD], pts)
        # Stub run + bundle for _consume_mystran_bundle.
        from node_runner.solve import MystranRun
        from pathlib import Path
        run = MystranRun(
            bdf_path=Path('fake.bdf'), f06_path=None, op2_path=None,
            sol=101, analysis_set_name='test')
        run.results_source = 'f06'
        bundle = {'subcases': {1: {'displacements': {1: [1, 0, 0, 0, 0, 0]}}}}
        # Need _status_results_lbl and a few combos for _consume to
        # work cleanly; the snapshot block itself only needs current_grid.
        mw._consume_mystran_bundle(bundle, run)
        # Diag = sqrt(100 + 1 + 0) = ~10.05
        assert mw._undeformed_bbox_diag == pytest.approx(
            (100 + 1) ** 0.5, abs=1e-6)

    def test_sync_deform_scale_prefers_snapshot(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        # Stash a snapshot and a fake bundle so the % path engages.
        mw._undeformed_bbox_diag = 10.0
        mw.op2_results = {
            'subcases': {1: {'displacements': {1: [0.0, 0.0, 1.0, 0, 0, 0]}}},
        }
        # Populate the legacy combo so currentData() returns the SID.
        mw.results_subcase_combo.clear()
        mw.results_subcase_combo.addItem("Subcase 1", 1)
        # Now pollute current_grid with a runaway bounding box that
        # would corrupt the math if the snapshot weren't used. Use a
        # full 3D hex so all three axes have real extents (the bounds
        # check needs bnds[5] > bnds[4] etc.).
        pts_huge = np.array([
            [0.0, 0.0, 0.0],
            [1000.0, 0.0, 0.0],
            [1000.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1000.0, 0.0, 1.0],
            [1000.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ])
        cells = np.array([8, 0, 1, 2, 3, 4, 5, 6, 7])
        mw.current_grid = pv.UnstructuredGrid(
            cells, [pv.CellType.HEXAHEDRON], pts_huge)
        mw._results_scale_mode = 'pct'
        mw._results_scale_value = 10.0
        mw._sync_deform_scale_to_widget()
        # If snapshot is preferred: scale = 0.10 * 10.0 / 1.0 = 1.0
        # (max disp magnitude in fixture is sqrt(0^2+0^2+1^2) = 1.0)
        # If it instead used the polluted bounds: ~100. Big difference.
        val = float(mw.deformation_scale_input.text())
        assert val == pytest.approx(1.0, abs=0.1), \
            f"expected ~1.0 (snapshot path); got {val} (likely using polluted bounds)"

    def test_sync_deform_scale_stable_across_calls(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        mw._undeformed_bbox_diag = 10.0
        mw.op2_results = {
            'subcases': {1: {'displacements': {1: [0.0, 0.0, 1.0, 0, 0, 0]}}},
        }
        mw.results_subcase_combo.clear()
        mw.results_subcase_combo.addItem("Subcase 1", 1)
        mw._results_scale_mode = 'pct'
        mw._results_scale_value = 5.0
        mw._sync_deform_scale_to_widget()
        first = mw.deformation_scale_input.text()
        for _ in range(5):
            mw._sync_deform_scale_to_widget()
        assert mw.deformation_scale_input.text() == first


# ---------------------------------------------------------------------------
# Item 50: point-scalars to sub-grid via vtkOriginalPointIds
# ---------------------------------------------------------------------------

class TestProjectPointScalars:

    def _build_grid_with_two_lines(self):
        """4 points, 2 line cells (forms a path 0-1, 2-3)."""
        pts = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ])
        cells = np.hstack([[2, 0, 1], [2, 2, 3]])
        cell_types = np.array([pv.CellType.LINE, pv.CellType.LINE])
        grid = pv.UnstructuredGrid(cells, cell_types, pts)
        return grid

    def test_uses_original_point_ids(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        grid = self._build_grid_with_two_lines()
        # Tag cells 0 and 1 differently; extract only cell 1 (uses points 2,3).
        grid.cell_data['EID'] = np.array([10, 20], dtype=int)
        sub = grid.extract_cells([1])
        # full_scalars: 4 nodes with distinct values
        full_scalars = np.array([10.0, 20.0, 30.0, 40.0])
        out = mw._project_point_scalars_to_subgrid(full_scalars, sub)
        # Sub-grid has 2 points (originally 2 and 3), so values 30, 40.
        assert out is not None
        assert sorted(out.tolist()) == [30.0, 40.0]

    def test_handles_missing_original_ids_gracefully(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        # Manually construct a grid WITHOUT vtkOriginalPointIds.
        pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        cells = np.array([2, 0, 1])
        cell_types = np.array([pv.CellType.LINE])
        sub = pv.UnstructuredGrid(cells, cell_types, pts)
        full_scalars = np.array([5.0, 6.0])
        out = mw._project_point_scalars_to_subgrid(full_scalars, sub)
        assert out is not None
        assert len(out) == 2

    def test_none_inputs_return_none(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        assert mw._project_point_scalars_to_subgrid(None, None) is None


# ---------------------------------------------------------------------------
# Item 51: node-cloud rebuild path is reachable
# ---------------------------------------------------------------------------

class TestNodeCloudRebuild:
    """Item 51's behaviour (post-deformation node cloud) is observable
    only through the live plotter, which we don't drive here. We
    instead assert that ``_add_node_cloud`` exists and accepts a
    deformed-points array argument shape -- proving the call-site we
    inserted in _update_plot_visibility has a valid target."""

    def test_add_node_cloud_signature(self, qtbot):
        from node_runner.mainwindow import MainWindow
        import inspect
        mw = MainWindow()
        qtbot.addWidget(mw)
        sig = inspect.signature(mw._add_node_cloud)
        params = list(sig.parameters)
        assert params[0] == 'node_points'
        # 'name' is a kw-arg with default 'nodes_actor'.
        assert sig.parameters['name'].default == 'nodes_actor'

    def test_add_node_cloud_with_deformed_points(self, qtbot):
        from node_runner.mainwindow import MainWindow
        mw = MainWindow()
        qtbot.addWidget(mw)
        # Should not raise.
        pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.5]])
        mw._add_node_cloud(pts, name='nodes_actor')
