"""v5.3.0 regression tests.

Each test pins one of the seven items from the v5.3.0 plan so a future
refactor that breaks them surfaces immediately. Uses pytest-qt's ``qapp``
and ``qtbot`` fixtures via the ``main_window`` fixture from test_ui.py
(we duplicate it here so this file can run in isolation).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture
def main_window(qapp, qtbot):
    from node_runner.theme import dark_palette, DARK_STYLESHEET
    from node_runner.mainwindow import MainWindow
    qapp.setStyle("Fusion")
    qapp.setPalette(dark_palette)
    qapp.setStyleSheet(DARK_STYLESHEET)
    w = MainWindow()
    qtbot.addWidget(w)
    return w


def _find_action(menu, path):
    """Walk a menu down the named path. Returns the leaf QAction or None."""
    if hasattr(menu, "menuBar"):
        menu = menu.menuBar()
    cur = menu
    for label in path:
        nxt = None
        for act in cur.actions():
            if act.text() == label:
                nxt = act.menu() if act.menu() is not None else act
                break
        if nxt is None:
            return None
        cur = nxt
    return cur if not hasattr(cur, "actions") else cur


# ---------------------------------------------------------------------------
# Item A - Model menu regrouped
# ---------------------------------------------------------------------------

class TestModelMenuRegrouped:

    @pytest.mark.parametrize("path", [
        ("&Model", "Geometry", "Nodes..."),
        ("&Model", "Geometry", "Coordinate System..."),
        ("&Model", "Sections", "Property..."),
        ("&Model", "Sections", "Material..."),
        ("&Model", "Loads", "Load..."),
        ("&Model", "Loads", "Distributed Load (PLOAD1)..."),
        ("&Model", "Constraints", "Constraint..."),
        ("&Model", "Constraints", "Multi-Point Constraint (MPC)..."),
        ("&Model", "Modify", "Flip Element Normals..."),
        ("&Model", "Modify", "Orient Normals (Auto)..."),
        ("&Model", "Modify", "Split Plates..."),
    ])
    def test_path_exists(self, main_window, path):
        act = _find_action(main_window, list(path))
        assert act is not None, f"missing path: {' > '.join(path)}"

    def test_no_duplicate_orient_normals_in_tools(self, main_window):
        """v5.3.0 Item A: the second Orient Normals entry under
        Tools > Model Check should be gone (canonical home is now
        Model > Modify)."""
        bar = main_window.menuBar()
        # Walk Tools > Model Check and ensure no Orient Normals action.
        for top in bar.actions():
            if top.text() != "&Tools":
                continue
            for sub in (top.menu() or []).actions() if top.menu() else []:
                if sub.text() == "Model Check" and sub.menu():
                    labels = [a.text() for a in sub.menu().actions()]
                    assert "Orient Normals (Auto)..." not in labels


# ---------------------------------------------------------------------------
# Item B - Solution Type submenu removed
# ---------------------------------------------------------------------------

class TestAnalysisMenuCleanup:

    def test_no_solution_type_submenu(self, main_window):
        bar = main_window.menuBar()
        for top in bar.actions():
            if top.text() != "&Analysis":
                continue
            for act in (top.menu() or []).actions() if top.menu() else []:
                assert act.text() != "Solution Type", (
                    "v5.3.0 Item B: Solution Type submenu must be gone.")
            return
        pytest.fail("Analysis menu not found.")


# ---------------------------------------------------------------------------
# Item C - Settings menu consolidated into Preferences
# ---------------------------------------------------------------------------

class TestSettingsConsolidation:

    def test_settings_menu_absent(self, main_window):
        bar = main_window.menuBar()
        labels = [a.text() for a in bar.actions()]
        assert "&Settings" not in labels

    def test_preferences_dialog_has_new_tabs(self):
        from node_runner.dialogs.preferences import PreferencesDialog
        for fn in ("_build_general_tab", "_build_export_tab",
                   "_build_performance_tab"):
            assert hasattr(PreferencesDialog, fn)

    def test_load_preferences_carries_new_keys(self):
        from node_runner.dialogs.preferences import load_preferences
        prefs = load_preferences()
        for key in ("export_format", "units_label", "adaptive_lod",
                    "elem_lod_threshold", "parallel_parsing"):
            assert key in prefs


# ---------------------------------------------------------------------------
# Item D - Groups bug fix + new tabs
# ---------------------------------------------------------------------------

class TestGroupsExtended:

    def test_create_group_command_has_full_schema(self):
        from node_runner.commands import CreateGroupCommand
        groups = {}
        cmd = CreateGroupCommand(groups, "Test", gid=1)
        cmd.execute(model=None)
        g = groups["Test"]
        for key in ("gid", "nodes", "elements",
                    "properties", "materials", "coords"):
            assert key in g, f"CreateGroupCommand schema missing {key!r}"

    def test_modify_group_command_supports_mpc(self):
        """The new add/remove kwargs for properties / materials / coords
        must round-trip cleanly."""
        from node_runner.commands import (
            CreateGroupCommand, ModifyGroupCommand)
        groups = {}
        CreateGroupCommand(groups, "G", gid=2).execute(None)
        ModifyGroupCommand(
            groups, "G",
            add_properties=[10, 20],
            add_materials=[1, 2, 3],
            add_coords=[100],
        ).execute(None)
        assert sorted(groups["G"]["properties"]) == [10, 20]
        assert sorted(groups["G"]["materials"]) == [1, 2, 3]
        assert sorted(groups["G"]["coords"]) == [100]
        ModifyGroupCommand(
            groups, "G",
            remove_materials=[2],
        ).execute(None)
        assert sorted(groups["G"]["materials"]) == [1, 3]

    def test_group_editor_dialog_has_five_tabs(self, qapp):
        from node_runner.dialogs.info import GroupEditorDialog
        # Stringified IDs should be coerced to int.
        group = {
            "gid": 1, "nodes": ["1", "2", 3], "elements": [10, "11"],
            "properties": [100], "materials": [1], "coords": [],
        }
        dlg = GroupEditorDialog("Test", group, model=None)
        try:
            assert dlg._tabs.count() == 5
            tab_labels = [dlg._tabs.tabText(i) for i in range(5)]
            for prefix in ("Nodes", "Elements", "Materials",
                           "Properties", "Coords"):
                assert any(t.startswith(prefix) for t in tab_labels), (
                    f"missing tab starting with {prefix}: {tab_labels}")
            # Int coercion.
            assert dlg._nodes == [1, 2, 3]
            assert dlg._elements == [10, 11]
        finally:
            dlg.deleteLater()

    def test_get_changes_returns_dict(self, qapp):
        from node_runner.dialogs.info import GroupEditorDialog
        group = {"gid": 1, "nodes": [], "elements": [],
                 "properties": [], "materials": [], "coords": []}
        dlg = GroupEditorDialog("Test", group, model=None)
        try:
            ch = dlg.get_changes()
            assert isinstance(ch, dict)
            for k in ("added_nodes", "added_elements", "removed_nodes",
                      "removed_elements", "added_properties",
                      "added_materials", "added_coords",
                      "removed_properties", "removed_materials",
                      "removed_coords"):
                assert k in ch
        finally:
            dlg.deleteLater()


# ---------------------------------------------------------------------------
# Item E - Probe Mode tooltip + shortcut
# ---------------------------------------------------------------------------

class TestProbeMode:

    def test_probe_action_has_tooltip_and_shortcut(self, main_window):
        act = main_window.probe_mode_action
        assert act.toolTip(), "Probe Mode should carry a tooltip."
        assert act.statusTip(), "Probe Mode should carry a status tip."
        assert act.shortcut().toString() == "Ctrl+Shift+P"


# ---------------------------------------------------------------------------
# Item F - Edit -> Select Entities + entity editor dispatch
# ---------------------------------------------------------------------------

class TestSelectEntitiesWorkflow:

    def test_show_entity_details_method_exists(self, main_window):
        assert hasattr(main_window, "_show_entity_details")
        assert hasattr(main_window, "_open_entity_editor")

    def test_open_find_entity_tool_does_not_use_qinputdialog(self):
        """The v5.3.0 rewrite must not call QInputDialog.getItem any
        more (which was the old "what type?" prompt). Static check
        for the specific .getItem call signature so docstring mentions
        of the historical behavior don't trip the test."""
        import inspect
        from node_runner.mainwindow import MainWindow
        src = inspect.getsource(MainWindow._open_find_entity_tool)
        assert "QInputDialog.getItem" not in src, (
            "v5.3.0 Item F: _open_find_entity_tool should no longer "
            "use QInputDialog.getItem; use the entity-type picker "
            "dialog instead.")


# ---------------------------------------------------------------------------
# Item G - View cube duplicate-actor fix
# ---------------------------------------------------------------------------

class TestViewCubeFix:

    def test_toggle_theme_cleans_orientation_widget(self, main_window):
        """Toggling theme repeatedly must not accumulate axes actors.
        We can't reliably count VTK-internal actors in headless mode,
        but we can verify the _toggle_theme code path runs through the
        teardown branch without raising."""
        main_window._toggle_theme()
        main_window._toggle_theme()
        main_window._toggle_theme()
        # If we got here without a NoneType / VTK crash, the teardown
        # path is intact.
        assert main_window.is_dark_theme in (True, False)


# ---------------------------------------------------------------------------
# v5.3.0 patch (2026-05-17): post-release fixes for plate_with_rbe2_spider
# ---------------------------------------------------------------------------

class TestGroupEditorRowPopulation:
    """P4: ensure all rows in the Edit Group dialog get their X/Y/Z (or
    Type/PID) populated, not just one row. Root cause was QTableWidget's
    setSortingEnabled(True) firing during the populate loop."""

    def test_all_node_rows_have_xyz(self, qapp):
        from node_runner.dialogs.info import GroupEditorDialog

        # Fake model with 4 nodes at distinct xyz.
        class _FakeNode:
            def __init__(self, xyz):
                self.xyz = xyz

        class _FakeModel:
            def __init__(self):
                self.nodes = {
                    1: _FakeNode([1.0, 2.0, 3.0]),
                    2: _FakeNode([4.0, 5.0, 6.0]),
                    3: _FakeNode([7.0, 8.0, 9.0]),
                    4: _FakeNode([10.0, 11.0, 12.0]),
                }
                self.elements = {}
                self.rigid_elements = {}
                self.properties = {}
                self.materials = {}
                self.coords = {}

        group = {
            "gid": 1,
            "nodes": [1, 2, 3, 4],
            "elements": [],
            "properties": [],
            "materials": [],
            "coords": [],
        }
        dlg = GroupEditorDialog("T", group, _FakeModel())
        try:
            from PySide6.QtCore import Qt
            t = dlg._nodes_table
            assert t.rowCount() == 4
            # Pre-sort: rows are 0..3 in NID-ascending order.
            for row in range(4):
                nid_item = t.item(row, 0)
                assert nid_item is not None
                nid = nid_item.data(Qt.DisplayRole)
                # X column for the row should match the fake-model lookup.
                x_item = t.item(row, 1)
                assert x_item is not None, (
                    f"P4 regression: row {row} (NID {nid}) has no X cell")
                x_val = x_item.data(Qt.DisplayRole)
                expected_x = float(nid * 3 - 2)  # 1,4,7,10 for NIDs 1..4
                assert abs(x_val - expected_x) < 1e-6, (
                    f"row {row} NID {nid}: expected X={expected_x}, "
                    f"got {x_val}")
        finally:
            dlg.deleteLater()


class TestDeformScaleDefaults:
    """P3: the default deform scale should be ('actual', 1.0)."""

    def test_default_scale_mode_actual_value_one(self, main_window):
        assert main_window._results_scale_mode == 'actual'
        assert main_window._results_scale_value == 1.0


class TestRbeActorsPointsSource:
    """P2: _create_rbe_actors gains a points_source kwarg so the spider
    geometry can follow the deformed mesh. The model-load call-site
    keeps the original (no-kwarg) signature."""

    def test_signature_accepts_points_source(self):
        import inspect
        from node_runner.mainwindow import MainWindow
        sig = inspect.signature(MainWindow._create_rbe_actors)
        assert 'points_source' in sig.parameters, (
            "P2: _create_rbe_actors must accept a points_source kwarg.")
        # Default must be None so existing call-sites keep working.
        assert (sig.parameters['points_source'].default is None), (
            "P2: points_source default must be None to preserve "
            "back-compat with the model-load path.")


class TestF06QuadStressParser:
    """P1B: _load_f06_minimal now parses CQUAD4 stress blocks."""

    def test_quad4_stress_from_fixture(self, tmp_path):
        # Hand-crafted minimal F06 with one QUAD4 stress block.
        fixture = tmp_path / "snip.f06"
        fixture.write_text(
            "\n"
            " OUTPUT FOR SUBCASE        1\n"
            " TITLE\n"
            "\n"
            "        E L E M E N T   S T R E S S E S   I N   "
            "L O C A L   E L E M E N T   C O O R D I N A T E   "
            "S Y S T E M\n"
            "                       F O R   E L E M E N T   "
            "T Y P E   Q U A D 4\n"
            "    Elem  Location       Fibre       Stresses In "
            "Element Coord System     Principal Stresses\n"
            "     ID                 Distance    Normal-X     "
            "Normal-Y     Shear-XY     Angle     Major        "
            "Minor      von Mises\n"
            "\n"
            "        1  CENTER      -2.500E-02 -3.22504E+05 "
            "-6.98562E+04 -1.33103E+04  -86.99 -6.91569E+04 "
            "-3.23204E+05  2.94774E+05\n"
            "                        2.500E-02  3.22504E+05  "
            "6.98562E+04  1.33103E+04    3.01  3.23204E+05  "
            "6.91569E+04  2.94774E+05\n"
            "\n"
        )
        from node_runner.solve.mystran_results import _load_f06_minimal
        bundle = _load_f06_minimal(str(fixture))
        sc = bundle['subcases'].get(1)
        assert sc is not None, (
            "P1B: subcase 1 should land in the bundle.")
        assert 1 in sc['stresses'], (
            f"P1B: EID 1 should be in stresses; got {list(sc['stresses'].keys())}")
        s = sc['stresses'][1]
        assert abs(s['von_mises'] - 2.94774e5) < 1.0, (
            f"P1B: von Mises expected ~2.948e5, got {s['von_mises']}")


class TestCaseControlOutputRequestsResolver:
    """P1A: _apply_case_control_from_set must fall back to analysis-set
    Output Requests when a subcase doesn't define disp/stress/etc."""

    def test_output_requests_reach_case_control(self, main_window):
        from node_runner.dialogs.analysis import AnalysisSet
        # Build a minimal AnalysisSet with stress=ALL set at the
        # analysis-set level and one subcase that does NOT declare
        # stress in its own dict.
        aset = AnalysisSet(id=1, name="T")
        aset.sol_type = 101
        aset.output_requests['stress'] = 'ALL'
        aset.output_requests['displacement'] = 'ALL'
        aset.subcases = [{'id': 1, 'load_sid': 1, 'spc_sid': 1}]

        # Apply against a fresh BDF model.
        from pyNastran.bdf.bdf import BDF
        model = BDF(debug=False)
        main_window._apply_case_control_from_set(model, aset)
        dump = str(model.case_control_deck)
        assert "STRESS" in dump.upper(), (
            f"P1A: STRESS=ALL should land in CASE CONTROL; got:\n{dump}")
        assert "DISPLACEMENT" in dump.upper(), (
            f"P1A: DISPLACEMENT=ALL should land in CASE CONTROL.")


# ---------------------------------------------------------------------------
# v5.3.0 patch 2 (2026-05-17): results-reading expansion
# ---------------------------------------------------------------------------

class TestBundleSchemaExtension:
    """Q1: _empty_subcase() carries strains + forces slots."""

    def test_keys_present(self):
        from node_runner.solve.mystran_results import _empty_subcase
        sc = _empty_subcase()
        for key in ('displacements', 'stresses', 'strains', 'forces',
                    'eigenvectors', 'spc_forces', 'eigenvalues',
                    'frequencies'):
            assert key in sc, f"Q1: bundle slot {key!r} missing"


class TestF06SpcForcesPopulated:
    """Q5: SPC FORCES section now lands in the bundle's spc_forces dict.
    Uses the existing fixture which carries an SPC FORCES block."""

    def test_spc_forces_from_fixture(self):
        from node_runner.solve.mystran_results import _load_f06_minimal
        bundle = _load_f06_minimal(
            "tests/fixtures/mystran/sol1_cantilever.f06")
        sc = bundle['subcases'].get(1)
        assert sc is not None
        assert sc['spc_forces'], (
            "Q5: spc_forces dict should be non-empty")


class TestF06Tria3StressParser:
    """Q2: CTRIA3 stress block parsed (Anywhere keyword + in elem
    continuation)."""

    def test_tria3_stress_from_fixture(self, tmp_path):
        fixture = tmp_path / "tria.f06"
        fixture.write_text(
            "\n OUTPUT FOR SUBCASE        1\n\n"
            "    E L E M E N T   S T R E S S E S   I N   L O C A L   "
            "E L E M E N T   C O O R D I N A T E   S Y S T E M\n"
            "                       F O R   E L E M E N T   T Y P E   "
            "T R I A 3\n"
            "  Element Location Fibre Stresses In Element Coord System "
            "Principal Stresses\n"
            "     ID         Distance Normal-X Normal-Y Shear-XY Angle "
            "Major Minor von Mises\n"
            "\n"
            "        4    Anywhere   -2.50000E-02 -9.75787E+04 "
            "-3.22010E+04  1.04774E+04   81.114 -3.05629E+04 "
            "-9.92167E+04  8.80096E+04 -4.24292E+02  2.80128E+02\n"
            "             in elem     2.50000E-02  9.75787E+04  "
            "3.22010E+04 -1.04774E+04   -8.886  9.92167E+04  "
            "3.05629E+04  8.80096E+04\n"
            "\n"
        )
        bundle = _load_minimal_for_test(str(fixture))
        sc = bundle['subcases'].get(1)
        assert sc is not None
        assert 4 in sc['stresses'], (
            f"Q2: EID 4 should be in stresses; got "
            f"{list(sc['stresses'].keys())}")
        s = sc['stresses'][4]
        assert abs(s['von_mises'] - 8.80096e4) < 1.0


class TestF06BarStressParser:
    """Q3: CBAR stress block parsed (End-A row with axial + End-B row
    without)."""

    def test_bar_stress_from_fixture(self, tmp_path):
        fixture = tmp_path / "bar.f06"
        fixture.write_text(
            "\n OUTPUT FOR SUBCASE        1\n\n"
            "      E L E M E N T   S T R E S S E S   I N   L O C A L   "
            "E L E M E N T   C O O R D I N A T E   S Y S T E M\n"
            "                       F O R   E L E M E N T   T Y P E   "
            "B A R\n"
            "  Element   SA1   SA2   SA3   SA4   Axial   SA-Max   "
            "SA-Min  M.S.-T\n"
            "    ID      SB1   SB2   SB3   SB4   Stress  SB-Max   "
            "SB-Min  M.S.-C\n"
            "\n"
            "       7  6.000E+03 -6.000E+03  6.000E+03 -6.000E+03  "
            "0.000E+00  6.000E+03 -6.000E+03\n"
            "          5.400E+03 -5.400E+03  5.400E+03 -5.400E+03  "
            "          5.400E+03 -5.400E+03\n"
            "\n"
        )
        bundle = _load_minimal_for_test(str(fixture))
        sc = bundle['subcases'].get(1)
        assert sc is not None
        assert 7 in sc['stresses']
        s = sc['stresses'][7]
        assert abs(s['sa_max'] - 6000.0) < 1e-3
        assert abs(s['sb_max'] - 5400.0) < 1e-3


class TestF06BarForcesParser:
    """Q4: CBAR engineering forces parsed (one row per element)."""

    def test_bar_forces_from_fixture(self, tmp_path):
        fixture = tmp_path / "barf.f06"
        fixture.write_text(
            "\n OUTPUT FOR SUBCASE        1\n\n"
            "      E L E M E N T   E N G I N E E R I N G   "
            "F O R C E S\n"
            "                       F O R   E L E M E N T   T Y P E   "
            "B A R\n"
            "  Element  Bend-Moment End A  Bend-Moment End B  - Shear - "
            "Axial  Torque\n"
            "    ID     Plane 1 Plane 2 Plane 1 Plane 2 Plane 1 Plane 2 "
            "Force\n"
            "\n"
            "       3  100.0  0.0  50.0  0.0  10.0  0.0  -20.0  5.0\n"
            "\n"
        )
        bundle = _load_minimal_for_test(str(fixture))
        sc = bundle['subcases'].get(1)
        assert sc is not None
        assert 3 in sc['forces']
        f = sc['forces'][3]
        assert abs(f['m1_a'] - 100.0) < 1e-6
        assert abs(f['axial'] - (-20.0)) < 1e-6
        assert abs(f['torque'] - 5.0) < 1e-6


class TestDropdownStrainForceCategories:
    """Q6: Output Vector dropdown auto-adds Strain + Force categories
    once their bundle slots are populated."""

    def test_strain_category_appears(self, qapp):
        from node_runner.dialogs.results_tab import ResultsTab
        bundle = {
            'subcases': {1: {
                'displacements': {1: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]},
                'stresses': {}, 'strains': {1: {'exx': 1e-3}},
                'forces': {}, 'eigenvectors': [], 'spc_forces': {},
                'eigenvalues': [], 'frequencies': [],
            }}
        }
        tab = ResultsTab()
        tab.populate_from_results(bundle, "test")
        labels = [tab.output_vector_combo.itemText(i)
                  for i in range(tab.output_vector_combo.count())]
        assert any(l.startswith("Strain -") for l in labels), (
            f"Q6: Strain category missing from dropdown. Got: {labels}")
        tab.deleteLater()

    def test_force_category_appears(self, qapp):
        from node_runner.dialogs.results_tab import ResultsTab
        bundle = {
            'subcases': {1: {
                'displacements': {1: [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]},
                'stresses': {}, 'strains': {},
                'forces': {1: {'axial': 100.0}},
                'eigenvectors': [], 'spc_forces': {},
                'eigenvalues': [], 'frequencies': [],
            }}
        }
        tab = ResultsTab()
        tab.populate_from_results(bundle, "test")
        labels = [tab.output_vector_combo.itemText(i)
                  for i in range(tab.output_vector_combo.count())]
        assert any(l.startswith("Force -") for l in labels), (
            f"Q6: Force category missing. Got: {labels}")
        tab.deleteLater()


class TestVectorInfoDialog:
    """Q8: VectorInfoDialog constructs + refreshes from a dict."""

    def test_dialog_constructs_and_renders_empty(self, qapp):
        from node_runner.dialogs.result_browser import VectorInfoDialog
        dlg = VectorInfoDialog()
        dlg.refresh(None)
        assert "no data" in dlg._max_label.text().lower()
        dlg.refresh({
            'subcase': 1, 'subcase_title': 'Static analysis',
            'label': 'displacement/Magnitude', 'kind': 'Node',
            'max_value': 43.45, 'max_id': 100,
            'min_value': 0.0, 'min_id': 1,
        })
        assert "43.45" in dlg._max_label.text()
        assert "Node 100" in dlg._max_id_label.text()
        assert "Node 1" in dlg._min_id_label.text()
        dlg.deleteLater()


# Helper used by the F06 parser tests above. Sits at module level so
# tests can share it without re-importing on every call.
def _load_minimal_for_test(path):
    from node_runner.solve.mystran_results import _load_f06_minimal
    return _load_f06_minimal(path)
