"""UI smoke tests: dialogs construct, MainWindow boots, menu wiring is intact.

These don't render anything visually - they confirm that every Theme A/B/C
widget instantiates without error and that the menu items the user expects
to find are actually wired up.

Uses pytest-qt's `qapp` and `qtbot` fixtures.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

# The order matters: a QApplication needs to exist before importing
# anything that constructs widgets at module level.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture
def main_window(qapp, qtbot):
    """Construct a fresh MainWindow per test, properly cleaned up."""
    from node_runner.theme import dark_palette, DARK_STYLESHEET
    from node_runner.mainwindow import MainWindow
    qapp.setStyle("Fusion")
    qapp.setPalette(dark_palette)
    qapp.setStyleSheet(DARK_STYLESHEET)
    w = MainWindow()
    qtbot.addWidget(w)
    return w


# ---------------------------------------------------------------------------
# MainWindow basics
# ---------------------------------------------------------------------------

class TestMainWindow:

    def test_window_title_has_version(self, main_window):
        # Match any vN.X family (v3.x through v4.x and beyond); specific
        # minor/patch version is stamped in node_runner/__init__.py.
        import re
        assert re.search(r"v\d+\.", main_window.windowTitle())

    def test_status_widgets_exist(self, main_window):
        assert hasattr(main_window, "_status_model_lbl")
        assert hasattr(main_window, "_status_nodes_lbl")
        assert hasattr(main_window, "_status_format_lbl")
        assert hasattr(main_window, "_status_units_lbl")

    def test_units_label_blank_by_default(self, main_window):
        # Femap-style: units field hides itself when no hint set
        assert main_window._status_units_lbl.text() == ""

    def test_command_palette_action_with_shortcut(self, main_window):
        # Command palette is reached via Tools > Command Palette... with
        # Ctrl+P bound to the action (not a separate QShortcut).
        assert main_window.command_palette_action is not None
        assert main_window.command_palette_action.shortcut().toString() == "Ctrl+P"

    def test_action_registry_populated(self, main_window):
        # We harvested at least 100 actions from the menu
        assert len(main_window._action_registry) > 100


# ---------------------------------------------------------------------------
# Menu wiring - verify every Theme A / C / etc menu item exists
# ---------------------------------------------------------------------------

def _find_action_path(menubar, top_name, *path):
    """Walk the menubar to find an action by its path. Returns the QAction
    or None if not found. `top_name` is e.g. '&Model'."""
    for top in menubar.actions():
        if top.text() != top_name:
            continue
        current_actions = top.menu().actions()
        for label in path:
            found = None
            for a in current_actions:
                if a.text() == label:
                    found = a
                    break
            if found is None:
                return None
            sub = found.menu()
            if sub is None:
                return found
            current_actions = sub.actions()
        return None
    return None


class TestThemeAMenuWiring:

    @pytest.mark.parametrize("label", [
        "Split QUAD into TRIA...",
        "Refine Elements (1-into-4)...",
        "Combine TRIA into QUAD...",
        "Smooth Nodes...",
        "Mirror Elements...",
        "Copy Elements...",
        "Insert Node on Edge...",
    ])
    def test_modify_menu_has(self, main_window, label):
        action = _find_action_path(
            main_window.menuBar(), "&Model", "Modify", label,
        )
        assert action is not None, f"missing Model > Modify > {label}"


class TestThemeCMenuWiring:

    @pytest.mark.parametrize("label", [
        "Distributed Load (PLOAD1)...",
        "Element Pressure (PLOAD2)...",
        "Enforced Displacement (SPCD)...",
        "Multi-Point Constraint (MPC)...",
        "Bolt Preload (BOLT)...",
    ])
    def test_model_menu_has(self, main_window, label):
        action = _find_action_path(main_window.menuBar(), "&Model", label)
        assert action is not None, f"missing Model > {label}"

    @pytest.mark.parametrize("label", [
        "Rigid Bar (RBAR)...",
        "General Rigid (RBE1)...",
        "Rigid Spline (RSPLINE)...",
    ])
    def test_connections_submenu_has(self, main_window, label):
        action = _find_action_path(
            main_window.menuBar(), "&Model", "Connections", label,
        )
        assert action is not None


class TestSettingsAndToolsMenu:

    def test_settings_export_defaults(self, main_window):
        action = _find_action_path(
            main_window.menuBar(), "&Settings", "Export Defaults...",
        )
        assert action is not None

    def test_settings_unit_label(self, main_window):
        # Renamed from "Units..." to be unitless-aware
        action = _find_action_path(main_window.menuBar(), "&Settings", "Units...")
        assert action is not None

    def test_tools_convert_units(self, main_window):
        action = _find_action_path(
            main_window.menuBar(), "&Tools", "Convert Units...",
        )
        assert action is not None

    def test_tools_probe_mode(self, main_window):
        action = _find_action_path(main_window.menuBar(), "&Tools", "Probe Mode")
        assert action is not None
        assert action.isCheckable()

    def test_view_cross_section(self, main_window):
        action = _find_action_path(
            main_window.menuBar(), "&View", "Cross Section...",
        )
        assert action is not None


# ---------------------------------------------------------------------------
# Cross-section dock + plane definitions
# ---------------------------------------------------------------------------

class TestCrossSectionDock:

    def test_dock_lazy_built_then_toggles(self, main_window, qtbot):
        # Built lazily on first menu trigger.
        assert main_window._cross_section_dock is None
        main_window.show()
        qtbot.waitExposed(main_window)
        main_window._open_cross_section_dialog()
        dock = main_window._cross_section_dock
        assert dock is not None
        assert dock.isVisible()
        # Second call hides it (toggle behavior).
        main_window._open_cross_section_dialog()
        assert not dock.isVisible()

    def test_plane_definition_axis_resolves(self):
        from node_runner.cross_section import PlaneDefinition, METHOD_AXIS
        d = PlaneDefinition(method=METHOD_AXIS, axis_label="+Y")
        origin, normal = d.resolve(None)
        assert origin == (0.0, 0.0, 0.0)
        assert normal == (0.0, 1.0, 0.0)

    def test_plane_definition_coord_value_resolves(self):
        from node_runner.cross_section import PlaneDefinition, METHOD_COORD_VALUE
        d = PlaneDefinition(method=METHOD_COORD_VALUE, axis_label="+X", value=1500.0)
        origin, normal = d.resolve(None)
        assert origin == (1500.0, 0.0, 0.0)
        assert normal == (1.0, 0.0, 0.0)

    def test_plane_definition_normal_point_resolves(self):
        from node_runner.cross_section import PlaneDefinition, METHOD_NORMAL_POINT
        d = PlaneDefinition(method=METHOD_NORMAL_POINT,
                            normal=(2.0, 0.0, 0.0), point=(5.0, 6.0, 7.0))
        origin, normal = d.resolve(None)
        assert origin == (5.0, 6.0, 7.0)
        assert normal == (1.0, 0.0, 0.0)  # normalized

    def test_plane_definition_three_point_needs_model(self):
        from node_runner.cross_section import (
            PlaneDefinition, METHOD_THREE_POINT, PlaneDefinitionError,
        )
        d = PlaneDefinition(method=METHOD_THREE_POINT, node_ids=(1, 2, 3))
        with pytest.raises(PlaneDefinitionError):
            d.resolve(None)

    def test_define_plane_dialog_constructs(self, main_window, qtbot):
        from node_runner.dialogs import DefinePlaneDialog
        dlg = DefinePlaneDialog(main_window)
        qtbot.addWidget(dlg)
        # Default method combo populated with 5 entries.
        assert dlg._method_combo.count() == 5


# ---------------------------------------------------------------------------
# Result browser dock
# ---------------------------------------------------------------------------

class TestResultBrowserDock:

    def test_dock_built(self, main_window):
        assert main_window._result_browser_dock is not None
        assert main_window._anim_timeline_widget is not None
        assert main_window._vector_overlay_widget is not None

    def test_dock_hidden_until_op2(self, main_window):
        # No OP2 loaded -> dock should be hidden
        assert not main_window._result_browser_dock.isVisible()


# ---------------------------------------------------------------------------
# Dialogs construct without error
# ---------------------------------------------------------------------------

class TestDialogConstruction:

    def test_export_options_dialog(self, qtbot):
        from node_runner.dialogs import ExportOptionsDialog
        d = ExportOptionsDialog(default_format="long")
        qtbot.addWidget(d)
        assert d.selected_format == "long"

    def test_smooth_dialog(self, qtbot):
        from node_runner.dialogs import SmoothNodesDialog
        d = SmoothNodesDialog()
        qtbot.addWidget(d)
        assert d.iterations == 5
        assert d.factor == 0.5
        assert d.pin_free_edges is True

    def test_mirror_dialog(self, qtbot):
        from node_runner.dialogs import MirrorElementsDialog
        d = MirrorElementsDialog(suggested_plane='X', suggested_value=1.5)
        qtbot.addWidget(d)
        assert d.plane == 'X'
        assert d.value == 1.5

    def test_copy_dialog_has_rotation(self, qtbot):
        from node_runner.dialogs import CopyElementsDialog
        d = CopyElementsDialog()
        qtbot.addWidget(d)
        assert hasattr(d, 'rotation_axis')
        assert hasattr(d, 'rotation_angle_deg')
        assert hasattr(d, 'rotation_center')

    def test_insert_edge_node_dialog(self, qtbot):
        from node_runner.dialogs import InsertEdgeNodeDialog
        d = InsertEdgeNodeDialog()
        qtbot.addWidget(d)
        # Default values for both endpoint spinboxes
        assert d.n1 == 1
        assert d.n2 == 1

    def test_unit_conversion_dialog_default_factors(self, qtbot):
        from node_runner.dialogs import UnitConversionDialog
        d = UnitConversionDialog()
        qtbot.addWidget(d)
        f = d.factors
        assert f == {'length': 1.0, 'force': 1.0, 'mass': 1.0}

    @pytest.mark.parametrize("preset_label, expected", [
        ("m -> mm  (length only)", (1000.0, 1.0, 1.0)),
        ("in -> mm", (25.4, 1.0, 1.0)),
        ("m,kg,N  -> mm,t,N", (1000.0, 1.0, 0.001)),
    ])
    def test_unit_conversion_presets(self, qtbot, preset_label, expected):
        from node_runner.dialogs import UnitConversionDialog
        d = UnitConversionDialog()
        qtbot.addWidget(d)
        d._preset_combo.setCurrentText(preset_label)
        f = d.factors
        L, F, M = expected
        assert f['length'] == pytest.approx(L)
        assert f['force'] == pytest.approx(F)
        assert f['mass'] == pytest.approx(M)
