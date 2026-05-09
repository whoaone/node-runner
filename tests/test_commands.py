"""Tests for every Theme A / Theme C / convert-units Command class.

Each command is exercised end-to-end (execute + undo) on a synthetic
model. These are the regression tests for mesh editing, load/BC creation,
and the unit converter.
"""
from __future__ import annotations

import pytest

from node_runner.commands import (
    SplitQuadCommand, RefineElementsCommand, CombineTriasCommand,
    SmoothNodesCommand, MirrorElementsCommand, CopyElementsCommand,
    InsertEdgeNodeCommand,
    AddPload1Command, AddPload2Command, AddSpcdCommand, AddMpcCommand,
    AddRbarCommand, AddRbe1Command, AddRsplineCommand, AddBoltCommand,
    ConvertUnitsCommand,
)


# ---------------------------------------------------------------------------
# Theme A - mesh editing
# ---------------------------------------------------------------------------

class TestSplitQuad:

    def test_split_one_quad_to_two_trias(self, tiny_gen):
        m = tiny_gen.model
        cmd = SplitQuadCommand([101])
        cmd.execute(m)
        assert 101 not in m.elements
        trias = [e for e in m.elements.values() if e.type == 'CTRIA3']
        assert len(trias) == 2

    def test_undo_restores_quad(self, tiny_gen):
        m = tiny_gen.model
        cmd = SplitQuadCommand([101])
        cmd.execute(m)
        cmd.undo(m)
        assert 101 in m.elements
        assert m.elements[101].type == 'CQUAD4'

    def test_skip_non_quads(self, tiny_gen):
        m = tiny_gen.model
        # 201 is a CBAR - should be silently skipped
        cmd = SplitQuadCommand([201])
        cmd.execute(m)
        assert 201 in m.elements


class TestRefineElements:

    def test_refine_quad_yields_4_quads_5_new_nodes(self, tiny_gen):
        m = tiny_gen.model
        cmd = RefineElementsCommand([101])
        cmd.execute(m)
        quads = [e for e in m.elements.values() if e.type == 'CQUAD4']
        assert len(quads) == 4
        # 4 edge midnodes + 1 centroid for a single quad
        assert len(m.nodes) == 6 + 5

    def test_refine_two_adjacent_quads_dedups_edge_midnode(self, two_quads_gen):
        m = two_quads_gen.model
        before = len(m.nodes)
        cmd = RefineElementsCommand([1, 2])
        cmd.execute(m)
        # 4 edges + 1 center per quad, but the shared edge midnode is shared.
        # 2*(4+1) - 1 shared = 9 new nodes
        assert len(m.nodes) - before == 9

    def test_undo_restores(self, tiny_gen):
        m = tiny_gen.model
        before_n = len(m.nodes)
        before_e = len(m.elements)
        cmd = RefineElementsCommand([101])
        cmd.execute(m)
        cmd.undo(m)
        assert len(m.nodes) == before_n
        assert len(m.elements) == before_e
        assert 101 in m.elements


class TestCombineTrias:

    def test_combine_round_trip(self, tiny_gen):
        m = tiny_gen.model
        # First split the quad to get two trias, then combine them.
        SplitQuadCommand([101]).execute(m)
        tria_eids = [eid for eid, e in m.elements.items() if e.type == 'CTRIA3']
        cmd = CombineTriasCommand(tria_eids, angle_tol_deg=5.0)
        cmd.execute(m)
        quads = [e for e in m.elements.values() if e.type == 'CQUAD4']
        assert len(quads) == 1


class TestSmoothNodes:

    def test_pinned_free_edge_node_does_not_move(self, tiny_gen):
        m = tiny_gen.model
        # All 4 corners of a single quad are on the free edge -> all pinned
        cmd = SmoothNodesCommand([1, 2, 3, 4], iterations=5,
                                 factor=0.5, pin_free_edges=True)
        cmd.execute(m)
        # Positions should be unchanged
        for nid in (1, 2, 3, 4):
            assert tuple(m.nodes[nid].get_position()) in (
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
            )

    def test_unpinned_smoothing_moves_node(self, two_quads_gen):
        m = two_quads_gen.model
        # Move node 2 (interior to the 2-quad strip) off-center
        m.nodes[2].set_position(m, [1.5, 0.5, 0.0])
        cmd = SmoothNodesCommand([2], iterations=10, factor=0.5,
                                 pin_free_edges=False)
        cmd.execute(m)
        new_pos = m.nodes[2].get_position()
        assert tuple(new_pos) != (1.5, 0.5, 0.0)


class TestMirrorElements:

    def test_mirror_creates_one_new_element(self, tiny_gen):
        m = tiny_gen.model
        before = len(m.elements)
        cmd = MirrorElementsCommand([101], plane='Y', plane_value=2.0)
        cmd.execute(m)
        assert len(m.elements) == before + 1

    def test_undo_removes_mirror(self, tiny_gen):
        m = tiny_gen.model
        before = len(m.elements)
        before_n = len(m.nodes)
        cmd = MirrorElementsCommand([101], plane='Y', plane_value=2.0)
        cmd.execute(m)
        cmd.undo(m)
        assert len(m.elements) == before
        assert len(m.nodes) == before_n


class TestCopyElements:

    def test_translate_only(self, tiny_gen):
        m = tiny_gen.model
        before = len(m.elements)
        cmd = CopyElementsCommand([101], translate=(5.0, 0.0, 0.0))
        cmd.execute(m)
        assert len(m.elements) == before + 1
        assert len(cmd._created_nids) == 4

    def test_rotation_about_z(self, tiny_gen):
        """90 deg rotation about Z + translate +5x: (1,0,0) -> (5,1,0)."""
        m = tiny_gen.model
        cmd = CopyElementsCommand(
            [101], translate=(5.0, 0.0, 0.0),
            axis='Z', angle_deg=90.0, center=(0.0, 0.0, 0.0),
        )
        cmd.execute(m)
        positions = [list(m.nodes[n].get_position()) for n in cmd._created_nids]
        # Original quad corners at (0,0,0), (1,0,0), (1,1,0), (0,1,0).
        # Rotate 90 deg about Z: (x,y) -> (-y, x). Then translate +5x.
        expected = {(5.0, 0.0, 0.0), (5.0, 1.0, 0.0),
                    (4.0, 1.0, 0.0), (4.0, 0.0, 0.0)}
        rounded = {tuple(round(v, 6) for v in p) for p in positions}
        assert rounded == expected


class TestInsertEdgeNode:

    def test_insert_on_shared_edge_splits_both_quads(self, two_quads_gen):
        m = two_quads_gen.model
        before_n = len(m.nodes)
        before_e = len(m.elements)
        cmd = InsertEdgeNodeCommand(2, 3)  # the shared edge
        cmd.execute(m)
        assert len(m.nodes) == before_n + 1
        # Each affected quad becomes 3 trias -> +4 elements net
        assert len(m.elements) == before_e + 4

    def test_undo_restores(self, two_quads_gen):
        m = two_quads_gen.model
        before_n = len(m.nodes)
        before_e = len(m.elements)
        cmd = InsertEdgeNodeCommand(2, 3)
        cmd.execute(m)
        cmd.undo(m)
        assert len(m.nodes) == before_n
        assert len(m.elements) == before_e
        assert 1 in m.elements and 2 in m.elements


# ---------------------------------------------------------------------------
# Theme C - loads & BCs
# ---------------------------------------------------------------------------

class TestLoadCommands:

    def test_pload1_creates_card(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddPload1Command(10, [201], {
            'load_type': 'FZ', 'scale': 'LE',
            'x1': 0.0, 'p1': 100.0, 'x2': 1.0, 'p2': 100.0,
        })
        cmd.execute(m)
        assert 10 in m.loads
        assert any(c.type == 'PLOAD1' for c in m.loads[10])

    def test_pload2_creates_card(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddPload2Command(20, [101], pressure=500.0)
        cmd.execute(m)
        assert 20 in m.loads
        assert any(c.type == 'PLOAD2' for c in m.loads[20])

    def test_spcd_creates_card(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddSpcdCommand(30, [1, 2], "13", 0.001)
        cmd.execute(m)
        assert 30 in m.loads

    def test_load_undo_clears_sid(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddPload2Command(20, [101], pressure=500.0)
        cmd.execute(m)
        cmd.undo(m)
        assert 20 not in m.loads


class TestConstraintCommands:

    def test_mpc_creates_card(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddMpcCommand(40, [(1, 1, 1.0), (2, 1, -1.0), (3, 1, 0.5)])
        cmd.execute(m)
        assert 40 in m.mpcs

    def test_rbar_creates_rigid(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddRbarCommand({
            'eid': 999, 'ga': 1, 'gb': 2,
            'cna': '123456', 'cnb': '', 'cma': '', 'cmb': '123456',
        })
        cmd.execute(m)
        assert 999 in m.rigid_elements

    def test_rbe1_creates_rigid(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddRbe1Command({
            'eid': 998,
            'indep_nodes': [1, 2], 'indep_dofs': '123456',
            'dep_nodes': [3, 4], 'dep_dofs': '123',
            'alpha': 0.0,
        })
        cmd.execute(m)
        assert 998 in m.rigid_elements

    def test_rspline_creates_rigid(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddRsplineCommand({
            'eid': 997, 'nodes': [1, 2, 3, 4], 'dofs': '123', 'diam': 0.1,
        })
        cmd.execute(m)
        assert 997 in m.rigid_elements

    def test_bolt_graceful_when_unavailable(self, tiny_gen):
        m = tiny_gen.model
        cmd = AddBoltCommand({'bid': 1, 'preload': 1000.0, 'eids': [201]})
        # Should not raise even if pyNastran lacks add_bolt
        cmd.execute(m)


# ---------------------------------------------------------------------------
# Convert Units
# ---------------------------------------------------------------------------

class TestConvertUnits:

    def test_m_to_mm_scales_coordinates(self, tiny_gen):
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=1.0)
        cmd.execute(m)
        assert tuple(m.nodes[2].get_position()) == (1000.0, 0.0, 0.0)

    def test_m_to_mm_scales_pshell_thickness(self, tiny_gen):
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=1.0)
        cmd.execute(m)
        assert m.properties[1].t == pytest.approx(1.0)

    def test_m_to_mm_scales_pbeam_area_by_L_squared(self, tiny_gen):
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=1.0)
        cmd.execute(m)
        # PBEAM stores per-station list
        assert m.properties[2].A[0] == pytest.approx(1e-4 * 1e6)

    def test_m_to_mm_scales_pbeam_J_by_L_to_4(self, tiny_gen):
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=1.0)
        cmd.execute(m)
        assert m.properties[2].j[0] == pytest.approx(1e-9 * 1e12, rel=1e-6)

    def test_full_unit_conversion_m_to_mm_t(self, tiny_gen):
        """Pa -> MPa (stress factor 1e-6); kg/m^3 -> t/mm^3 (density factor 1e-12)."""
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=0.001)
        cmd.execute(m)
        assert m.materials[1].e == pytest.approx(2e11 * 1e-6)  # 2e5 MPa
        assert m.materials[1].rho == pytest.approx(7800 * 1e-12)

    def test_undo_restores_everything(self, tiny_gen):
        m = tiny_gen.model
        cmd = ConvertUnitsCommand(length=1000.0, force=1.0, mass=0.001)
        cmd.execute(m)
        cmd.undo(m)
        assert tuple(m.nodes[2].get_position()) == (1.0, 0.0, 0.0)
        assert m.properties[1].t == pytest.approx(0.001)
        assert m.materials[1].e == pytest.approx(2e11)
