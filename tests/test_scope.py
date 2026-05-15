"""v5.1.0 item 26: tests for scope_model_to_analysis_set.

The scoping logic is built on top of an in-memory pyNastran BDF
object. We construct minimal BDFs in-place (no file I/O) and exercise
the four interesting cases: full model, group-only, group + load
filter, group + SPC filter, missing-group fallback.
"""

from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_model():
    """A tiny BDF with two materials, two properties, two groups of
    elements (Skin / Frames), two load SIDs, two SPC SIDs.

    Returns (model, groups) where ``groups`` is the MainWindow-style
    dict of name -> {gid, nodes, elements, properties, materials,
    coords}.
    """
    from pyNastran.bdf.bdf import BDF
    m = BDF(debug=False)

    # Nodes 1..6
    for nid in range(1, 7):
        m.add_grid(nid, [float(nid), 0.0, 0.0])

    # Two materials
    m.add_mat1(1, 1.0e7, None, 0.3, 0.1)
    m.add_mat1(2, 2.0e7, None, 0.3, 0.2)

    # Two shell properties (Skin uses MAT1; Frames uses MAT2)
    m.add_pshell(1, 1, 0.05)
    m.add_pshell(2, 2, 0.08)

    # Two CQUAD4s on PSHELL 1 (Skin) -- EIDs 101, 102
    m.add_cquad4(101, 1, [1, 2, 5, 4])
    m.add_cquad4(102, 1, [2, 3, 6, 5])

    # Two CQUAD4s on PSHELL 2 (Frames) -- EIDs 201, 202
    m.add_cquad4(201, 2, [1, 4, 5, 2])
    m.add_cquad4(202, 2, [2, 5, 6, 3])

    # Two FORCE LOAD SIDs (1, 2), two SPC SIDs (11, 22)
    m.add_force(1, 1, 100.0, [0., 0., -1.])
    m.add_force(2, 4, 50.0, [0., 0., -1.])
    m.add_spc1(11, '123456', [1])
    m.add_spc1(22, '123', [4])

    groups = {
        'Skin': {
            'gid': 1,
            'nodes': [1, 2, 3, 4, 5, 6],
            'elements': [101, 102],
            'properties': [1],
            'materials': [1],
            'coords': [],
        },
        'Frames': {
            'gid': 2,
            'nodes': [1, 2, 3, 4, 5, 6],
            'elements': [201, 202],
            'properties': [2],
            'materials': [2],
            'coords': [],
        },
    }
    return m, groups


# ---------------------------------------------------------------------------
# Scope cases
# ---------------------------------------------------------------------------

class TestScopeFullModel:

    def test_no_group_no_filter_keeps_everything(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=1, name='full', group_target=None,
                           load_sids=[], spc_sids=[])
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert len(scoped.nodes) == 6
        assert len(scoped.elements) == 4
        assert len(scoped.properties) == 2
        assert len(scoped.materials) == 2
        assert len(scoped.loads) == 2
        assert len(scoped.spcs) == 2
        # Source model untouched (shallow copy).
        assert len(model.elements) == 4


class TestScopeGroupOnly:

    def test_group_target_skin_keeps_only_skin_elements(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=2, name='skin_only', group_target='Skin')
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert set(scoped.elements.keys()) == {101, 102}
        # PSHELL 1 + MAT1 1 auto-included via property/material walk;
        # PSHELL 2 + MAT1 2 must have been dropped.
        assert set(scoped.properties.keys()) == {1}
        assert set(scoped.materials.keys()) == {1}
        assert report.n_elements_kept == 2
        assert report.n_elements_dropped == 2
        assert report.group_target == 'Skin'

    def test_group_target_frames_keeps_only_frames(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=3, name='frames_only', group_target='Frames')
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert set(scoped.elements.keys()) == {201, 202}
        assert set(scoped.properties.keys()) == {2}
        assert set(scoped.materials.keys()) == {2}


class TestScopeLoadSpcFilter:

    def test_load_sid_filter_keeps_only_listed(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=4, name='load1_only',
                           load_sids=[1], spc_sids=[])
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert set(scoped.loads.keys()) == {1}
        assert report.n_loads_kept == 1
        assert report.n_loads_dropped == 1

    def test_spc_sid_filter_keeps_only_listed(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=5, name='spc22_only',
                           load_sids=[], spc_sids=[22])
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert set(scoped.spcs.keys()) == {22}
        assert report.n_spcs_kept == 1
        assert report.n_spcs_dropped == 1

    def test_combined_group_and_filters(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=6, name='skin_load1_spc11',
                           group_target='Skin',
                           load_sids=[1], spc_sids=[11])
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert set(scoped.elements.keys()) == {101, 102}
        assert set(scoped.loads.keys()) == {1}
        assert set(scoped.spcs.keys()) == {11}


class TestScopeMissingGroup:

    def test_missing_group_falls_back_to_full_model(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(id=7, name='ghost', group_target='DoesNotExist')
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert len(scoped.elements) == 4   # unchanged
        assert any('not found' in n for n in report.notes)


class TestSubcaseFiltering:

    def test_subcase_referencing_dropped_load_is_filtered(self, simple_model):
        from node_runner.scope import scope_model_to_analysis_set
        from node_runner.dialogs.analysis import AnalysisSet
        model, groups = simple_model
        aset = AnalysisSet(
            id=8, name='partial',
            load_sids=[1], spc_sids=[11],
            subcases=[
                {'id': 1, 'load_sid': 1, 'spc_sid': 11},
                # This one references the dropped LOAD SID 2:
                {'id': 2, 'load_sid': 2, 'spc_sid': 11},
            ],
        )
        scoped, report = scope_model_to_analysis_set(model, aset, groups)
        assert report.n_subcases_kept == 1
        assert report.n_subcases_dropped == 1
