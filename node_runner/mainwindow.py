import sys
import os
import json
import math
import random
import re
import time
from pathlib import Path
import numpy as np
import pyvista as pv
import vtk
import copy
from pyvistaqt import QtInteractor

from PySide6 import QtCore, QtGui
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QMessageBox, QFormLayout, QGroupBox, QFrame,
    QLineEdit, QComboBox, QLabel, QCheckBox, QColorDialog, QDialog,
    QScrollArea, QGridLayout, QTreeWidget, QTreeWidgetItem, QSplitter,
    QDialogButtonBox, QTreeWidgetItemIterator, QTableWidget, QTableWidgetItem,
    QHeaderView, QListWidget, QListWidgetItem, QTextEdit, QTabWidget, QRadioButton,
    QButtonGroup, QStackedLayout, QInputDialog, QMenu, QSlider, QSpinBox,
    QDockWidget, QToolButton,
)
from PySide6.QtGui import (QPalette, QColor, QAction, QActionGroup, QDoubleValidator,
                           QPainter, QPen, QBrush, QPolygonF, QPixmap, QIcon)
from collections import Counter

from pyNastran.bdf.bdf import BDF
from node_runner.model import NastranModelGenerator
from node_runner.utils import get_entity_title_from_comment, PROJECT_ROOT
from node_runner.theme import dark_palette, light_palette, DARK_STYLESHEET, BACKGROUND_PRESETS
from node_runner.interactor import ClickAndDragInteractor
from node_runner.geometry import GeometryStore
from node_runner.commands import (
    CommandManager, CompoundCommand,
    AddNodesCommand, DeleteNodesCommand, TransformNodesCommand, MergeNodesCommand,
    AddElementCommand, DeleteElementsCommand, FlipNormalsCommand, EditElementCommand,
    AddMaterialCommand, EditMaterialCommand, DeleteMaterialCommand,
    AddPropertyCommand, EditPropertyCommand, DeletePropertyCommand,
    AddLoadCommand, DeleteLoadCommand, AddConstraintCommand, DeleteConstraintCommand,
    AddCoordCommand, DeleteCoordCommand, AddMassCommand, DeleteMassCommand,
    AddPlotelCommand, DeletePlotelCommand,
    CreateGroupCommand, DeleteGroupCommand, RenameGroupCommand, ModifyGroupCommand,
    ReassignPropertyCommand, ReassignMaterialCommand,
    AddLoadCombinationCommand, EditLoadCombinationCommand,
    CopyLoadCombinationCommand, RenumberCommand,
    AddGeometryPointCommand, AddGeometryLineCommand, AddGeometryArcCommand,
    AddGeometryCircleCommand, AddGeometrySurfaceCommand,
    MeshCurveCommand, MeshSurfaceCommand, NodesAtGeometryPointsCommand,
    DeleteGeometryPointsCommand, DeleteGeometryCurvesCommand,
    DeleteGeometrySurfacesCommand,
    AutoOrientNormalsCommand,
    ImportCADCommand,
    AddSpiderCommand, AddWeldCommand,
    EditGeometryPointCommand, EditGeometryCurveCommand,
    EditGeometrySurfaceBoundariesCommand,
    TransformGeometryPointsCommand, CopyGeometryCommand,
    SplitCurveCommand, OffsetCurveCommand, ProjectPointToCurveCommand,
    FilletCommand,
    ReplaceLoadSetCommand, ReplaceConstraintSetCommand,
)
from node_runner.dialogs import (
    MaterialEditorDialog, PropertyEditorDialog, ElementEditorDialog,
    CreateMaterialDialog, CreatePropertyDialog,
    CreateNodesDialog, CreateLineElementDialog, CreatePlateElementDialog,
    CreateSolidElementDialog, CreateShearElementDialog, CreateGapElementDialog,
    CreateBushElementDialog, CreateRbeDialog, CreateConm2Dialog, CreateCoordDialog,
    CreatePlotelDialog,
    EntitySelectionDialog,
    CoincidentNodeDialog, GeneratorDialog, ColorManagerDialog,
    NodeTransformDialog, ImportOptionsDialog, DuplicateElementDialog,
    FindReplaceDialog, RenumberDialog, ImportCADDialog,
    CreateSpiderDialog, CreateWeldDialog,
    InfoDialog, LenientImportReportDialog, FreeBodyDiagramDialog,
    MassPropertiesReportDialog, FreeEdgeReportDialog,
    QualitySummaryReportDialog, OrphanCheckDialog,
    CreateLoadDialog, CreateConstraintDialog, CreateLoadCombinationDialog,
    SubcaseEditorDialog,
    CreateGeometryPointDialog, CreateGeometryLineDialog,
    CreateGeometryArcDialog, CreateGeometryCircleDialog,
    CreateGeometrySurfaceDialog, MeshCurveDialog, MeshSurfaceDialog,
)


# v5.0.0 item 10: tiny clickable QLabel for status-bar hints. Emits
# a `clicked` signal on mouse press so we can wire it to an existing
# action callback (e.g. command palette). Kept inline here so it lives
# next to the only consumer in `_build_status_widgets`.
class _ClickableLabel(QLabel):
    clicked = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCursor(QtCore.Qt.PointingHandCursor)

    def mousePressEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(ev)


# v5.0.0 item 10: tiny clickable QLabel for status-bar hints. Emits
# a `clicked` signal on mouse press so we can wire it to an existing
# action callback (e.g. command palette). Kept inline here so it lives
# next to the only consumer in `_build_status_widgets`.
class _ClickableLabel(QLabel):
    clicked = QtCore.Signal()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCursor(QtCore.Qt.PointingHandCursor)

    def mousePressEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(ev)


# v4.0.11 (Fix C): module-level Nastran-type ↔ category lookup tables.
# Pre-built per-category actors at viewer-build time are named
# `cat_<NTYPE>` and toggled via these mappings.
#
# Source of truth: the inline `type_to_nastran_map` and
# `shape_to_nastran_map` in `_update_plot_visibility` (which we leave
# alone to keep the legacy path intact). The forward maps below
# duplicate that data; the inverse maps are pre-computed so the fast
# path can resolve compound visibility (type ∧ shape) in O(1) per
# ntype.

_TYPE_CATEGORY_TO_NTYPES: dict[str, list[str]] = {
    'Beams': ['CBEAM'],
    'Bars': ['CBAR'],
    'Rods': ['CROD'],
    'Bushes': ['CBUSH'],
    'Plates': ['CQUAD4', 'CMEMBRAN', 'CTRIA3'],
    'Rigid': ['RBE2', 'RBE3'],
    'Solids': ['CHEXA', 'CHEXA8', 'CHEXA20',
               'CTETRA', 'CTETRA4', 'CTETRA10',
               'CPENTA', 'CPENTA6', 'CPENTA15'],
    'Shear': ['CSHEAR'],
    'Gap': ['CGAP'],
}
_SHAPE_CATEGORY_TO_NTYPES: dict[str, list[str]] = {
    'Line': ['CBEAM', 'CBAR', 'CROD', 'CBUSH', 'CGAP'],
    'Quad': ['CQUAD4', 'CMEMBRAN', 'CSHEAR'],
    'Tria': ['CTRIA3'],
    'Rigid': ['RBE2', 'RBE3'],
    'Hex': ['CHEXA', 'CHEXA8', 'CHEXA20'],
    'Tet': ['CTETRA', 'CTETRA4', 'CTETRA10'],
    'Wedge': ['CPENTA', 'CPENTA6', 'CPENTA15'],
}
_NTYPE_TO_TYPE_CATEGORY: dict[str, str] = {
    ntype: cat
    for cat, ntypes in _TYPE_CATEGORY_TO_NTYPES.items()
    for ntype in ntypes
}
_NTYPE_TO_SHAPE_CATEGORY: dict[str, str] = {
    ntype: cat
    for cat, ntypes in _SHAPE_CATEGORY_TO_NTYPES.items()
    for ntype in ntypes
}

# v4.0.15: color modes the pre-built per-Nastran-type actors can
# render at build time. Adding a new mode? If the data lives on
# cell_data and a colormap is per-cell (like PID), add the name
# here AND a branch in `_build_pre_built_category_actors`'s color-
# computation block. Otherwise the legacy `_rebuild_plot` path
# handles it (e.g., 'quality' / 'results' which depend on external
# data flow).
_PRE_BUILT_SUPPORTED_COLOR_MODES: frozenset = frozenset({'type', 'property'})


# v5.2.0 item 43: human-friendly suffix for the active data conversion
# mode shown in the scalar-bar title. ``average`` is the default and
# stays unannotated.
_CONV_LABEL: dict = {
    'average':  'avg',
    'no_avg':   'no avg',
    'max_node': 'max@node',
    'min_node': 'min@node',
}


class SelectionOverlay:
    """VTK 2D actor overlay for selection feedback: box, circle, polygon.

    Uses vtkActor2D + vtkPolyDataMapper2D so drawing happens inside
    VTK's own render pipeline - no Qt widget compositing issues.
    Coordinates are in VTK display space (origin at bottom-left).
    """

    ACCENT = (0.537, 0.706, 0.980)  # #89b4fa

    def __init__(self, plotter):
        self._plotter = plotter
        self._line_actor = None
        self._vert_actor = None
        self._active = False

    # ---- internal helpers ------------------------------------------------

    def _ensure_actors(self):
        """Create the reusable 2D actors on first use."""
        if self._line_actor is None:
            mapper = vtk.vtkPolyDataMapper2D()
            mapper.SetInputData(vtk.vtkPolyData())
            self._line_actor = vtk.vtkActor2D()
            self._line_actor.SetMapper(mapper)
            self._line_actor.GetProperty().SetColor(*self.ACCENT)
            self._line_actor.GetProperty().SetLineWidth(2)
            self._line_actor.GetProperty().SetOpacity(0.85)
        if self._vert_actor is None:
            mapper = vtk.vtkPolyDataMapper2D()
            mapper.SetInputData(vtk.vtkPolyData())
            self._vert_actor = vtk.vtkActor2D()
            self._vert_actor.SetMapper(mapper)
            self._vert_actor.GetProperty().SetColor(*self.ACCENT)
            self._vert_actor.GetProperty().SetPointSize(8)

    def _show(self):
        self._ensure_actors()
        if not self._active:
            self._plotter.renderer.AddActor2D(self._line_actor)
            self._plotter.renderer.AddActor2D(self._vert_actor)
            self._active = True
        self._line_actor.SetVisibility(True)
        self._vert_actor.SetVisibility(True)

    @staticmethod
    def _make_rect_pd(x0, y0, x1, y1):
        pts = vtk.vtkPoints()
        for xy in [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]:
            pts.InsertNextPoint(xy[0], xy[1], 0)
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(5)
        for i in range(4):
            lines.InsertCellPoint(i)
        lines.InsertCellPoint(0)
        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetLines(lines)
        return pd

    @staticmethod
    def _make_circle_pd(cx, cy, r, n_seg=64):
        pts = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(n_seg + 1)
        for i in range(n_seg):
            a = 2 * np.pi * i / n_seg
            pts.InsertNextPoint(cx + r * np.cos(a), cy + r * np.sin(a), 0)
            lines.InsertCellPoint(i)
        lines.InsertCellPoint(0)
        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetLines(lines)
        return pd

    # ---- public API ------------------------------------------------------

    def update_box(self, start_vtk, end_vtk):
        self._show()
        pd = self._make_rect_pd(
            float(start_vtk[0]), float(start_vtk[1]),
            float(end_vtk[0]), float(end_vtk[1]))
        self._line_actor.GetMapper().SetInputData(pd)
        self._vert_actor.GetMapper().SetInputData(vtk.vtkPolyData())
        self._plotter.render()

    def update_circle(self, center_vtk, radius):
        self._show()
        pd = self._make_circle_pd(
            float(center_vtk[0]), float(center_vtk[1]), float(radius))
        self._line_actor.GetMapper().SetInputData(pd)
        self._vert_actor.GetMapper().SetInputData(vtk.vtkPolyData())
        self._plotter.render()

    def update_polygon(self, points_vtk, cursor_vtk=None):
        self._show()
        all_pts = list(points_vtk)
        if cursor_vtk:
            all_pts.append(cursor_vtk)

        # Lines between consecutive points (+ dashed line to cursor)
        line_pd = vtk.vtkPolyData()
        if len(all_pts) >= 2:
            vpts = vtk.vtkPoints()
            lines = vtk.vtkCellArray()
            for px, py in all_pts:
                vpts.InsertNextPoint(float(px), float(py), 0)
            for i in range(len(all_pts) - 1):
                lines.InsertNextCell(2)
                lines.InsertCellPoint(i)
                lines.InsertCellPoint(i + 1)
            line_pd.SetPoints(vpts)
            line_pd.SetLines(lines)
        self._line_actor.GetMapper().SetInputData(line_pd)

        # Vertex dots (only committed points, not cursor)
        vert_pd = vtk.vtkPolyData()
        if points_vtk:
            vpts2 = vtk.vtkPoints()
            verts = vtk.vtkCellArray()
            for i, (px, py) in enumerate(points_vtk):
                vpts2.InsertNextPoint(float(px), float(py), 0)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(i)
            vert_pd.SetPoints(vpts2)
            vert_pd.SetVerts(verts)
        self._vert_actor.GetMapper().SetInputData(vert_pd)
        self._plotter.render()

    def hide_overlay(self):
        if self._active:
            if self._line_actor:
                self._line_actor.SetVisibility(False)
            if self._vert_actor:
                self._vert_actor.SetVisibility(False)
            self._active = False
            self._plotter.render()


class MainWindow(QMainWindow):
    # v5.2.0 item 45: signal emitted whenever the viewport selection
    # changes. The Data Table dialog (Tools -> Data Table) subscribes
    # to live-sync its rows. Payload: (kind, ids) where kind is
    # 'Node' or 'Element'.
    selection_changed = Signal(str, list)

    def __init__(self):
        super().__init__()
        from node_runner import __version__ as _nr_version
        self.setWindowTitle(f"Node Runner v{_nr_version}")
        self.setGeometry(100, 100, 1200, 800)
        self.is_dark_theme, self.current_generator, self.current_grid = True, None, None
        self.shell_opacity, self.color_mode, self.render_style = 1.0, "property", "surface"
        # Phase 3: smaller default size and theme accent (Catppuccin blue),
        # paired with vtkPointGaussianMapper-rendered nodes for a crisp,
        # anti-aliased look instead of the old neon-green blob.
        from node_runner.theme import ACCENT as _THEME_ACCENT
        self.node_size, self.node_color = 3, _THEME_ACCENT
        self.edge_color, self.edge_width = '#000000', 1
        self.beam_width = 2
        self.elem_shrink = 0.0
        self.load_scaling_info = {}
        self.type_color_map, self.pid_color_map = {"Shells": "#0077be", "Beams": "#f85a40", "RBE2": "#ff3131", "RBE3": "#ffd700", "Masses": "#00cc66", "Solids": "#8b5cf6", "Shear": "#ff8c00", "Gap": "#20b2aa", "Plotel": "#ff69b4"}, {}

        # v3.5.0 item 3: load preferences from QSettings. Falls back
        # to defaults defined in dialogs/preferences.py if no settings
        # saved yet.
        try:
            from node_runner.dialogs.preferences import (
                load_preferences, DEFAULT_HIGHLIGHT_COLOR)
            prefs = load_preferences()
            # Merge persisted colors over the hardcoded defaults.
            for key, value in (prefs.get('colors') or {}).items():
                self.type_color_map[key] = value
            self.highlight_color = prefs.get(
                'highlight_color', DEFAULT_HIGHLIGHT_COLOR)
            sizes = prefs.get('sizes') or {}
            # Mass glyph scale: stored as percent, used as fraction.
            # v4.0.7: default lowered from 1.5% to 0.75% (cubes were
            # too prominent on dense models).
            self.mass_glyph_scale = float(
                sizes.get('mass_glyph_scale_pct', 0.75)) / 100.0
            self.node_size = int(sizes.get('node_size', self.node_size))
            self.beam_width = int(sizes.get('beam_width', self.beam_width))
            self.edge_width = int(sizes.get('edge_width', self.edge_width))
            self.rbe_line_width = int(sizes.get('rbe_line_width', 3))
            self.free_edge_width = int(sizes.get('free_edge_width', 4))
            self.highlight_outline_width = int(
                sizes.get('highlight_outline_width', 5))
        except Exception:
            from node_runner.theme import SELECTION_ACCENT
            self.highlight_color = SELECTION_ACCENT
            self.mass_glyph_scale = 0.0075  # v4.0.7: half of pre-v4.0.7 default
            self.rbe_line_width = 3
            self.free_edge_width = 4
            self.highlight_outline_width = 5
        
        # --- FIX: Separated chained assignments onto individual lines ---
        self.active_selection_dialog = None
        self.current_node_ids_sorted = []
        self.default_interactor_style = None

        self.picking_target_callback = None
        self.active_sub_dialog = None
        self.active_creation_dialog = None
        self.current_selection_type = None
        self.quality_results = None

        # v3.3.0: element group counts cached by _update_viewer from
        # the scene-build pass so _populate_tree doesn't re-iterate all
        # elements to compute By-Type/By-Shape Counters.
        self._cached_element_group_counts = {'by_type': {}, 'by_shape': {}}

        # v3.4.2: lazy load-actor caches. _create_all_load_actors does
        # scaling + arrays only; per-SID actors are built on demand by
        # _ensure_load_actor_for_sid when the visibility update path
        # actually needs them. Caps import time on decks with hundreds
        # or thousands of load SIDs (a production deck we exercised had
        # 1,962 SIDs).
        self._lazy_load_centers = None
        self._lazy_load_normals = None
        self._lazy_load_eid_to_idx = None
        self._lazy_load_built_sids: set = set()
        # v4.0.10: maintained sets of currently-visible load / SPC SIDs.
        # Updated by _fast_path_sid_toggle (O(1)) and refreshed from
        # tree state at every full visibility cycle entry
        # (defensive). The load_sid_loop / spc_sid_loop iterate
        # ONLY these sets instead of all model.loads.keys() / spcs.keys(),
        # collapsing the 1,962-iteration cost to N (whatever the user
        # has checked).
        self._visible_load_sids: set = set()
        self._visible_spc_sids: set = set()
        # v4.0.14 (Fix C v2): pre-built per-Nastran-type actors,
        # redesigned after v4.0.11 was reverted. Built from the FULL
        # `current_grid` (not the LOD'd `current_display_grid`.
        # that was the v4.0.11 visual regression: ~10% of plates
        # showed because the LOD shell-stride sampler reduces shells
        # to 200k and drops RBE2/RBE3 entirely). With current_grid:
        # full-density visual matches legacy; RBE2/RBE3 actors are
        # built. Category toggles become SetVisibility calls.
        # PID/MID/isolate/hidden_groups force fall-through to legacy
        # _rebuild_plot via `_using_legacy_mesh`.
        self._category_actors: dict = {}
        # v5.0.0 item 11a: per-cell visibility via vtkGhostType cell-data
        # on each `cat_<ntype>` actor's underlying grid. Sub-100 ms PID
        # and MID toggles flip these bits in-place; the legacy
        # _rebuild_plot path no longer engages for plain PID/MID filters.
        self._category_actor_grids: dict = {}
        self._hidden_pids: set = set()
        self._hidden_mids: set = set()
        # v5.0.0 item 11a: per-cell visibility via vtkGhostType cell-data
        # on each `cat_<ntype>` actor's underlying grid. Sub-100 ms PID
        # and MID toggles flip these bits in-place; the legacy
        # _rebuild_plot path no longer engages for plain PID/MID filters.
        self._category_actor_grids: dict = {}
        self._hidden_pids: set = set()
        self._hidden_mids: set = set()
        self._category_visible_ntypes: set = set()
        self._using_legacy_mesh: bool = False

        # --- Measurement tools state ---
        self._measurement_picks = []
        self._measurement_pick_count = 0
        self._measurement_callback = None

        # --- Phase 1: export format / units / loaded-file metadata for status bar ---
        from PySide6.QtCore import QSettings
        _settings = QSettings("NodeRunner", "NodeRunner")
        self._export_format = _settings.value("export/default_format", "short")
        if self._export_format not in ("short", "long", "free"):
            self._export_format = "short"
        # Node Runner is intentionally unitless (professional). The "units"
        # value is just a free-text label shown in the status bar so you
        # remember which unit system you've been working in. To actually
        # rescale model values, use Tools > Convert Units.
        self._units = _settings.value("display/units", "")
        if self._units in ("SI", "English"):
            # Migrate legacy values from earlier 3.0.0 builds.
            self._units = ""
        self._loaded_filepath: str | None = None

        # --- Phase 4: performance and LOD settings ---
        self._adaptive_lod_enabled = _settings.value(
            "render/adaptive_lod", True, type=bool,
        )
        self._ghost_mode_enabled = _settings.value(
            "render/ghost_mode", False, type=bool,
        )
        # v3.2.0: element-level LOD threshold. Above this cell count,
        # the displayed mesh stride-samples shells and extracts surface
        # for solids; current_grid still has every cell for picking.
        self._elem_lod_threshold = _settings.value(
            "render/elem_lod_threshold", 500_000, type=int,
        )
        # The display-decimated grid. Equal to current_grid when LOD is
        # inactive; a smaller pv.UnstructuredGrid when active.
        self.current_display_grid = None
        # NOTE: Parallel parsing was an experimental opt-in. It turned out
        # to be slower than single-threaded due to GIL contention on
        # pyNastran's pure-Python card parsing, and a stale "True" left
        # in the user's QSettings caused massive import hangs. We reset
        # the flag to False on every launch and ignore it in the worker.
        # The flag stays in QSettings so existing reads don't crash.
        self._parallel_parsing_enabled = False
        if _settings.value("import/parallel_parsing", False, type=bool):
            _settings.setValue("import/parallel_parsing", False)

        self.command_manager = CommandManager(max_history=20)
        self.subcases = []  # Case control deck subcase definitions (legacy)
        self.geometry_store = GeometryStore()
        self.sol_type = None  # None, 101, 103, 105 (legacy)
        self.eigrl_cards = []  # [{sid, v1, v2, nd}] (legacy)

        # --- Analysis Set Manager ---
        from node_runner.dialogs.analysis import AnalysisSet
        self.analysis_sets: dict[int, AnalysisSet] = {}
        self.active_analysis_set_id: int | None = None
        self.op2_results = None
        self.deformation_scale = 0.0
        self._anim_override_scale = None
        self._anim_phase = 0.0
        self._anim_timer = QtCore.QTimer(self)
        self._anim_timer.timeout.connect(self._on_animation_tick)

        # --- Theme B: result interaction state ---
        self._probe_enabled = False
        self._probe_tooltip = None  # lazy-created QLabel tooltip
        self._result_browser_dock = None  # built in initUI when OP2 loads
        self._vector_overlay_widget = None
        self._anim_timeline_widget = None
        # Active expression for custom contour. Empty string -> normal contour.
        self._contour_expression = ""

        # --- Phase 8: Group management ---
        self.groups: dict[str, dict] = {}  # { "name": {"nodes": [], "elements": [], "gid": int} }
        self._next_group_id: int = 1
        self._hidden_groups: set[str] = set()
        self._isolate_mode: str | None = None
        self._current_selection_eids: list[int] = []
        self._current_selection_nids: list[int] = []
        self._highlighted_group: str | None = None

        # --- Enhanced Entity Selection state ---
        self._previous_node_selection: set = set()
        self._previous_element_selection: set = set()
        self._pick_selection_mode: str = 'box'  # persists across picking sessions

        # --- Display settings (persisted via QSettings) ---
        self._display_settings = self._load_display_settings()

        self.initUI()
        self._create_new_model()
        self._auto_load_materials()
        self.tree_widget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tree_widget.customContextMenuRequested.connect(self._show_tree_context_menu)


    @staticmethod
    def _get_entity_title_from_comment(comment, entity_type, entity_id):
        return get_entity_title_from_comment(comment, entity_type, entity_id)


    def initUI(self):
        self._create_menu_bar()
        self._build_status_widgets()
        # Action registry is built AFTER the menu so we can index every
        # QAction the menu created. The Ctrl+P shortcut is owned by
        # self.command_palette_action under Tools (set via setShortcut),
        # which makes it both discoverable in the menu and reachable via
        # the keyboard. We do not also install a stand-alone QShortcut -
        # having two bindings on the same key produces an "ambiguous
        # shortcut" warning and either path can win unpredictably.
        self._action_registry = self._build_action_registry()
        self._command_palette_shortcut = None  # legacy attribute, kept for tests
        self._build_result_browser_dock()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        left_panel_container = QWidget()
        left_panel_layout = QVBoxLayout(left_panel_container)
        self.tree_widget = QTreeWidget()
        self.tree_widget.setHeaderLabels(["Model Tree"])

        # --- Phase 8: Groups tab ---
        self.groups_list = QTreeWidget()
        self.groups_list.setHeaderLabels(["Groups"])
        self.groups_list.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.groups_list.customContextMenuRequested.connect(self._show_groups_context_menu)
        self.groups_list.itemChanged.connect(self._on_group_item_changed)
        self.groups_list.itemDoubleClicked.connect(self._edit_group_contents)

        groups_widget = QWidget()
        groups_layout = QVBoxLayout(groups_widget)
        groups_layout.setContentsMargins(4, 4, 4, 4)
        groups_layout.setSpacing(4)

        # v3.4.0 item 7: professional top toolbar (two rows). All
        # buttons that used to live below the group list move up here.
        tb_row1 = QHBoxLayout(); tb_row1.setSpacing(3)
        for label, slot, tip in (
                ("New",    self._create_group,             "Create a new (empty) group"),
                ("Delete", self._delete_group,             "Delete the highlighted group"),
                ("Rename", self._rename_selected_group,    "Rename the highlighted group"),
        ):
            b = QPushButton(label); b.setToolTip(tip)
            b.clicked.connect(slot)
            tb_row1.addWidget(b)
        sep1 = QFrame(); sep1.setFrameShape(QFrame.VLine); sep1.setFrameShadow(QFrame.Sunken)
        tb_row1.addWidget(sep1)

        # Add dropdown: rule-based + From Selection
        self._add_to_group_btn = QToolButton()
        self._add_to_group_btn.setText("Add ▾")
        self._add_to_group_btn.setToolTip(
            "Add entities to the highlighted group")
        self._add_to_group_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._add_to_group_btn.setPopupMode(QToolButton.MenuButtonPopup)
        add_menu = QMenu(self)
        add_menu.addAction(
            "From current selection", self._add_selected_to_group)
        add_menu.addAction(
            "From rules...", lambda: self._open_group_add_dialog('add'))
        self._add_to_group_btn.setMenu(add_menu)
        # Default click action: From current selection (the most-used path).
        self._add_to_group_btn.clicked.connect(self._add_selected_to_group)
        tb_row1.addWidget(self._add_to_group_btn)

        # Remove dropdown: mirror of Add
        self._rem_from_group_btn = QToolButton()
        self._rem_from_group_btn.setText("Remove ▾")
        self._rem_from_group_btn.setToolTip(
            "Remove entities from the highlighted group")
        self._rem_from_group_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._rem_from_group_btn.setPopupMode(QToolButton.MenuButtonPopup)
        rem_menu = QMenu(self)
        rem_menu.addAction(
            "From current selection", self._remove_selected_from_group)
        rem_menu.addAction(
            "From rules...",
            lambda: self._open_group_add_dialog('remove'))
        self._rem_from_group_btn.setMenu(rem_menu)
        self._rem_from_group_btn.clicked.connect(self._remove_selected_from_group)
        tb_row1.addWidget(self._rem_from_group_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.setToolTip("Remove all items from the highlighted group")
        clear_btn.clicked.connect(self._clear_group)
        tb_row1.addWidget(clear_btn)
        tb_row1.addStretch(1)
        groups_layout.addLayout(tb_row1)

        # Row 2: visibility + auto + boolean
        tb_row2 = QHBoxLayout(); tb_row2.setSpacing(3)
        show_all_btn = QPushButton("Show All")
        show_all_btn.setToolTip("Show all groups (clear isolation)")
        show_all_btn.clicked.connect(self._show_all_groups)
        isolate_btn = QPushButton("Isolate")
        isolate_btn.setToolTip("Show only the highlighted group")
        isolate_btn.clicked.connect(self._isolate_group)
        highlight_btn = QPushButton("Highlight")
        highlight_btn.setToolTip("Flash-highlight entities in this group")
        highlight_btn.clicked.connect(self._highlight_group)
        for b in (show_all_btn, isolate_btn, highlight_btn):
            tb_row2.addWidget(b)
        sep2 = QFrame(); sep2.setFrameShape(QFrame.VLine); sep2.setFrameShadow(QFrame.Sunken)
        tb_row2.addWidget(sep2)

        # Auto: rule-based bulk group creation.
        self._auto_btn = QToolButton()
        self._auto_btn.setText("Auto ▾")
        self._auto_btn.setToolTip(
            "Auto-create groups by Property / Material / Element Type")
        self._auto_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._auto_btn.setPopupMode(QToolButton.InstantPopup)
        auto_menu = QMenu(self)
        auto_menu.addAction("By Property", self._group_by_property)
        auto_menu.addAction("By Material", self._group_by_material)
        auto_menu.addAction("By Element Type", self._group_by_element_type)
        self._auto_btn.setMenu(auto_menu)
        tb_row2.addWidget(self._auto_btn)

        # Boolean: set-arithmetic between two groups.
        self._boolean_btn = QToolButton()
        self._boolean_btn.setText("Boolean ▾")
        self._boolean_btn.setToolTip(
            "Union / Intersect / Difference between two groups")
        self._boolean_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._boolean_btn.setPopupMode(QToolButton.InstantPopup)
        bool_menu = QMenu(self)
        bool_menu.addAction("Union (A + B)",
                            lambda: self._group_boolean_op('union'))
        bool_menu.addAction("Intersect (A ∩ B)",
                            lambda: self._group_boolean_op('intersect'))
        bool_menu.addAction("Difference (A - B)",
                            lambda: self._group_boolean_op('difference'))
        self._boolean_btn.setMenu(bool_menu)
        tb_row2.addWidget(self._boolean_btn)
        tb_row2.addStretch(1)
        groups_layout.addLayout(tb_row2)

        # Group list lives BELOW the toolbar (v3.4.0 - was below buttons).
        groups_layout.addWidget(self.groups_list, 1)

        # auto_group_combo kept as a hidden attribute for back-compat
        # with any external caller; not displayed.
        self.auto_group_combo = QComboBox()
        self.auto_group_combo.addItems(
            ["By Property", "By Material", "By Element Type"])
        self.auto_group_combo.hide()

        # v3.2.4: search bar above the model tree. Search supports two
        # modes via the radio toggles:
        #   Filter: hide tree items that don't match (default)
        #   Find:   scroll-and-highlight the next matching item
        # Match logic: case-insensitive substring match against item
        # display text, which includes IDs and labels (e.g. "GRID 12345",
        # "CID 10: Rectangular - my_csys").
        model_tab_container = QWidget()
        model_tab_layout = QVBoxLayout(model_tab_container)
        model_tab_layout.setContentsMargins(4, 4, 4, 4)
        model_tab_layout.setSpacing(3)
        search_row = QHBoxLayout()
        search_row.setSpacing(4)
        self._tree_search = QLineEdit()
        self._tree_search.setPlaceholderText(
            "Search tree (ID, type, label)...")
        self._tree_search.setClearButtonEnabled(True)
        search_row.addWidget(self._tree_search, 1)
        self._tree_search_filter_radio = QRadioButton("Filter")
        self._tree_search_filter_radio.setChecked(True)
        self._tree_search_filter_radio.setToolTip(
            "Hide tree items that don't match")
        self._tree_search_find_radio = QRadioButton("Find")
        self._tree_search_find_radio.setToolTip(
            "Scroll the next matching item into view")
        _tree_search_grp = QButtonGroup(self)
        _tree_search_grp.addButton(self._tree_search_filter_radio)
        _tree_search_grp.addButton(self._tree_search_find_radio)
        search_row.addWidget(self._tree_search_filter_radio)
        search_row.addWidget(self._tree_search_find_radio)
        model_tab_layout.addLayout(search_row)
        model_tab_layout.addWidget(self.tree_widget, 1)
        self._tree_search.textChanged.connect(self._on_tree_search_changed)
        self._tree_search.returnPressed.connect(self._on_tree_search_next)
        self._tree_search_filter_radio.toggled.connect(
            lambda on: on and self._on_tree_search_changed(
                self._tree_search.text()))
        self._tree_search_find_radio.toggled.connect(
            lambda on: on and self._tree_search_clear_filter())
        # Tracks where 'Find' mode is in the cycle.
        self._tree_search_find_cursor = 0

        self.sidebar_tabs = QTabWidget()
        self.sidebar_tabs.addTab(model_tab_container, "Model")

        # --- Loads Tab (unified tree, mirrors Model tab layout) ---
        loads_tab = QWidget()
        loads_tab_layout = QVBoxLayout(loads_tab)
        loads_tab_layout.setContentsMargins(4, 4, 4, 4)
        loads_tab_layout.setSpacing(4)

        self.loads_tab_tree = QTreeWidget()
        self.loads_tab_tree.setHeaderLabels(["Loads / Constraints"])
        self.loads_tab_tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.loads_tab_tree.customContextMenuRequested.connect(self._show_loads_tab_context_menu)
        loads_tab_layout.addWidget(self.loads_tab_tree, 1)  # stretch factor 1 = fill space

        # Button row
        loads_btn_layout = QHBoxLayout()
        add_load_btn = QPushButton("Add Load")
        add_load_btn.clicked.connect(self._create_load)
        add_constr_btn = QPushButton("Add Constraint")
        add_constr_btn.clicked.connect(self._create_constraint)
        add_combo_btn = QPushButton("Add Combo")
        add_combo_btn.clicked.connect(self._create_load_combination)
        loads_btn_layout.addWidget(add_load_btn)
        loads_btn_layout.addWidget(add_constr_btn)
        loads_btn_layout.addWidget(add_combo_btn)
        loads_tab_layout.addLayout(loads_btn_layout)

        self.sidebar_tabs.addTab(loads_tab, "Loads")

        self.sidebar_tabs.addTab(groups_widget, "Groups")

        # --- Analysis Tab ---
        analysis_tab = QWidget()
        analysis_tab_layout = QVBoxLayout(analysis_tab)
        analysis_tab_layout.setContentsMargins(4, 4, 4, 4)
        analysis_tab_layout.setSpacing(4)

        self.analysis_tab_tree = QTreeWidget()
        self.analysis_tab_tree.setHeaderLabels(["Analysis Sets"])
        self.analysis_tab_tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.analysis_tab_tree.customContextMenuRequested.connect(
            self._show_analysis_tab_context_menu)
        analysis_tab_layout.addWidget(self.analysis_tab_tree, 1)

        analysis_btn_layout = QHBoxLayout()
        new_set_btn = QPushButton("New Set")
        new_set_btn.clicked.connect(self._new_analysis_set_from_tab)
        edit_set_btn = QPushButton("Edit Set")
        edit_set_btn.clicked.connect(self._open_analysis_set_manager)
        set_active_btn = QPushButton("Set Active")
        set_active_btn.clicked.connect(self._set_active_from_tab)
        analysis_btn_layout.addWidget(new_set_btn)
        analysis_btn_layout.addWidget(edit_set_btn)
        analysis_btn_layout.addWidget(set_active_btn)
        analysis_tab_layout.addLayout(analysis_btn_layout)

        self.sidebar_tabs.addTab(analysis_tab, "Analysis")

        left_panel_layout.addWidget(self.sidebar_tabs)
        # --- Display Options Tab (organized with sections, scrollable) ---
        display_scroll = QScrollArea()
        display_scroll.setWidgetResizable(True)
        display_scroll.setFrameShape(QFrame.NoFrame)
        display_inner = QWidget()
        display_tab_layout = QVBoxLayout(display_inner)
        display_tab_layout.setContentsMargins(6, 6, 6, 6)
        display_tab_layout.setSpacing(8)

        # Section: Render
        render_group = QGroupBox("Render")
        render_layout = QFormLayout(render_group)
        render_layout.setSpacing(6)
        self.render_style_combo = QComboBox()
        self.render_style_combo.addItems(["Surface with Edges", "Wireframe", "Surface Only", "Points"])
        self.render_style_combo.currentIndexChanged.connect(self._on_render_style_combo_changed)
        render_layout.addRow("Style:", self.render_style_combo)
        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems(["Property", "Type", "Quality", "Results"])
        self.color_mode_combo.currentIndexChanged.connect(self._on_color_mode_combo_changed)
        render_layout.addRow("Color By:", self.color_mode_combo)
        # v3.4.0 item 8: Shading checkbox in the Display tab. Mirrors
        # the View > Shading menu action so both stay in sync.
        self.shading_display_check = QCheckBox()
        # v4.0.1: default OFF to match the View > Shading action.
        self.shading_display_check.setChecked(False)
        self.shading_display_check.toggled.connect(
            self._on_shading_display_toggled)
        render_layout.addRow("Shading:", self.shading_display_check)
        self.bg_preset_combo = QComboBox()
        self.bg_preset_combo.addItems(list(BACKGROUND_PRESETS.keys()))
        self.bg_preset_combo.currentTextChanged.connect(self._on_bg_preset_changed)
        render_layout.addRow("Background:", self.bg_preset_combo)
        display_tab_layout.addWidget(render_group)

        # Section: Labels
        labels_group = QGroupBox("Labels")
        labels_layout = QFormLayout(labels_group)
        labels_layout.setSpacing(6)
        self.show_node_labels_check = QCheckBox()
        self.show_elem_labels_check = QCheckBox()
        labels_layout.addRow("Node IDs:", self.show_node_labels_check)
        labels_layout.addRow("Element IDs:", self.show_elem_labels_check)
        display_tab_layout.addWidget(labels_group)

        # Section: Node Display
        node_disp_group = QGroupBox("Node Display")
        node_disp_layout = QFormLayout(node_disp_group)
        node_disp_layout.setSpacing(6)
        self.node_size_spin = QSpinBox()
        self.node_size_spin.setRange(1, 20)
        self.node_size_spin.setValue(self.node_size)
        self.node_size_spin.valueChanged.connect(self._on_node_size_changed)
        node_disp_layout.addRow("Node Size:", self.node_size_spin)
        self.node_color_btn = QPushButton()
        self.node_color_btn.setFixedSize(60, 22)
        self.node_color_btn.setStyleSheet(f"background-color: {self.node_color}; border: 1px solid #555;")
        self.node_color_btn.clicked.connect(self._pick_node_color)
        node_disp_layout.addRow("Node Color:", self.node_color_btn)
        display_tab_layout.addWidget(node_disp_group)

        # Section: Element Display
        elem_disp_group = QGroupBox("Element Display")
        elem_disp_layout = QFormLayout(elem_disp_group)
        elem_disp_layout.setSpacing(6)
        self.shell_opacity_slider = QSlider(QtCore.Qt.Horizontal)
        self.shell_opacity_slider.setRange(0, 100)
        self.shell_opacity_slider.setValue(100)
        self.shell_opacity_slider.setToolTip("Shell element opacity (0-100%)")
        self.shell_opacity_slider.valueChanged.connect(self._on_shell_opacity_changed)
        elem_disp_layout.addRow("Shell Opacity:", self.shell_opacity_slider)
        self.elem_shrink_slider = QSlider(QtCore.Qt.Horizontal)
        self.elem_shrink_slider.setRange(0, 50)
        self.elem_shrink_slider.setValue(0)
        self.elem_shrink_slider.setToolTip("Shrink elements toward centroids (0-50%)")
        self.elem_shrink_slider.valueChanged.connect(self._on_elem_shrink_changed)
        elem_disp_layout.addRow("Element Shrink:", self.elem_shrink_slider)
        self.edge_color_btn = QPushButton()
        self.edge_color_btn.setFixedSize(60, 22)
        self.edge_color_btn.setStyleSheet(f"background-color: {self.edge_color}; border: 1px solid #555;")
        self.edge_color_btn.clicked.connect(self._pick_edge_color)
        elem_disp_layout.addRow("Edge Color:", self.edge_color_btn)
        self.edge_width_spin = QSpinBox()
        self.edge_width_spin.setRange(1, 5)
        self.edge_width_spin.setValue(self.edge_width)
        self.edge_width_spin.valueChanged.connect(self._on_edge_width_changed)
        elem_disp_layout.addRow("Edge Width:", self.edge_width_spin)
        self.beam_width_spin = QSpinBox()
        self.beam_width_spin.setRange(1, 10)
        self.beam_width_spin.setValue(self.beam_width)
        self.beam_width_spin.valueChanged.connect(self._on_beam_width_changed)
        elem_disp_layout.addRow("Beam Width:", self.beam_width_spin)
        display_tab_layout.addWidget(elem_disp_group)

        # Section: Element Visualization
        elem_viz_group = QGroupBox("Element Visualization")
        elem_viz_layout = QFormLayout(elem_viz_group)
        elem_viz_layout.setSpacing(6)
        self.show_normals_check = QCheckBox()
        self.normal_arrow_scale_input = QLineEdit("5")
        self.normal_arrow_scale_input.setValidator(QDoubleValidator(0.1, 100.0, 2, self))
        self.show_free_edges_check = QCheckBox()
        elem_viz_layout.addRow("Element Normals:", self.show_normals_check)
        elem_viz_layout.addRow("Normal Arrow Size (%):", self.normal_arrow_scale_input)
        elem_viz_layout.addRow("Free Edges:", self.show_free_edges_check)
        display_tab_layout.addWidget(elem_viz_group)

        # Section: Orientations
        orient_group = QGroupBox("Orientations")
        orient_layout = QFormLayout(orient_group)
        orient_layout.setSpacing(6)
        self.show_bush_orient_check = QCheckBox()
        self.show_beam_orient_check = QCheckBox()
        orient_layout.addRow("Bush Orientations:", self.show_bush_orient_check)
        orient_layout.addRow("Beam Orientations:", self.show_beam_orient_check)
        self.show_beam_sections_check = QCheckBox()
        orient_layout.addRow("Beam 3D Sections:", self.show_beam_sections_check)
        display_tab_layout.addWidget(orient_group)

        # Section: Coordinate Systems
        csys_group = QGroupBox("Coordinate Systems")
        csys_layout = QFormLayout(csys_group)
        csys_layout.setSpacing(6)
        self.show_csys_check = QCheckBox()
        self.show_csys_check.setChecked(True)
        self.show_global_csys_check = QCheckBox()
        self.show_global_csys_check.setToolTip("Show/hide CID 0 (Rectangular), 1 (Cylindrical), 2 (Spherical)")
        self.show_csys_labels_check = QCheckBox()
        self.show_csys_labels_check.setChecked(True)
        self.csys_scale_slider = QSlider(QtCore.Qt.Horizontal)
        self.csys_scale_slider.setRange(1, 30)
        self.csys_scale_slider.setValue(8)
        self.csys_scale_slider.setToolTip("Axis arrow size (‰ of model diagonal)")
        csys_layout.addRow("Show Coord Systems:", self.show_csys_check)
        csys_layout.addRow("Show Global (0,1,2):", self.show_global_csys_check)
        csys_layout.addRow("Show Labels:", self.show_csys_labels_check)
        csys_layout.addRow("Axis Size:", self.csys_scale_slider)
        display_tab_layout.addWidget(csys_group)

        # Section: Loads & Constraints
        loads_group = QGroupBox("Loads && Constraints")
        loads_layout = QFormLayout(loads_group)
        loads_layout.setSpacing(6)
        self.arrow_scale_input = QLineEdit("10.0")
        self.arrow_scale_input.setValidator(QDoubleValidator(0.1, 100.0, 2, self))
        self.relative_scaling_check = QCheckBox()
        self.relative_scaling_check.setChecked(True)
        loads_layout.addRow("Arrow Size (%):", self.arrow_scale_input)
        loads_layout.addRow("Relative Scaling:", self.relative_scaling_check)
        display_tab_layout.addWidget(loads_group)

        # Hidden quality metric widget (shown when quality coloring is active)
        self.quality_widget = QWidget()
        quality_layout = QFormLayout(self.quality_widget)
        quality_layout.setContentsMargins(0, 0, 0, 0)
        self.quality_metric_combo = QComboBox()
        quality_layout.addRow("Quality Metric:", self.quality_metric_combo)
        display_tab_layout.addWidget(self.quality_widget)
        self.quality_widget.setVisible(False)

        # Hidden results widget (shown when OP2 results loaded)
        self.results_widget = QGroupBox("Results")
        results_layout = QFormLayout(self.results_widget)
        results_layout.setContentsMargins(6, 6, 6, 6)
        results_layout.setSpacing(6)
        self.results_subcase_combo = QComboBox()
        results_layout.addRow("Subcase:", self.results_subcase_combo)
        self.results_type_combo = QComboBox()
        self.results_type_combo.addItems(["Displacement", "Stress", "Eigenvector", "SPC Forces"])
        results_layout.addRow("Result Type:", self.results_type_combo)
        self.results_component_combo = QComboBox()
        results_layout.addRow("Component:", self.results_component_combo)
        self.deformation_scale_input = QLineEdit("0.0")
        self.deformation_scale_input.setValidator(QDoubleValidator(0.0, 1e12, 4, self))
        results_layout.addRow("Deform Scale:", self.deformation_scale_input)
        self.results_mode_combo = QComboBox()
        results_layout.addRow("Mode:", self.results_mode_combo)
        self.results_mode_combo.setVisible(False)
        anim_widget = QWidget()
        anim_layout = QHBoxLayout(anim_widget)
        anim_layout.setContentsMargins(0, 0, 0, 0)
        self.anim_play_btn = QPushButton("Play")
        self.anim_play_btn.setCheckable(True)
        self.anim_speed_slider = QSlider(QtCore.Qt.Horizontal)
        self.anim_speed_slider.setRange(10, 200)
        self.anim_speed_slider.setValue(80)
        self.anim_speed_slider.setToolTip("Animation speed")
        self.anim_export_gif_btn = QPushButton("GIF")
        self.anim_export_gif_btn.setToolTip("Export animation as GIF")
        anim_layout.addWidget(self.anim_play_btn)
        anim_layout.addWidget(self.anim_speed_slider, 1)
        anim_layout.addWidget(self.anim_export_gif_btn)
        results_layout.addRow("Animate:", anim_widget)
        display_tab_layout.addWidget(self.results_widget)
        # v5.2.0 item 47: the legacy "Results" QGroupBox in the Display
        # tab is kept as a hidden state holder -- the combos still feed
        # _update_plot_visibility, but the new Post-Processing Toolbox
        # (Results sidebar tab) is the visible interface. Hide
        # permanently so the user never sees both UIs at once.
        self.results_widget.setVisible(False)
        self._legacy_results_widget_hidden = True
        self.results_subcase_combo.currentIndexChanged.connect(self._on_results_changed)
        self.results_type_combo.currentIndexChanged.connect(self._on_results_type_changed)
        self.results_component_combo.currentIndexChanged.connect(self._on_results_changed)
        self.results_mode_combo.currentIndexChanged.connect(self._on_results_changed)
        self.deformation_scale_input.editingFinished.connect(self._on_results_changed)
        self.anim_play_btn.toggled.connect(self._toggle_animation)
        self.anim_speed_slider.valueChanged.connect(self._on_anim_speed_changed)
        self.anim_export_gif_btn.clicked.connect(self._export_results_gif)

        # v4.0.0 (E3): "Show All Entities" button - one-click recovery
        # from accidentally hiding everything. Re-checks every tree
        # checkbox + re-evaluates plot visibility.
        self.show_all_btn = QPushButton("Show All Entities")
        self.show_all_btn.setToolTip(
            "Re-enable visibility for every entity in the Model tree.")
        self.show_all_btn.clicked.connect(self._show_all_entities)
        display_tab_layout.addWidget(self.show_all_btn)

        display_tab_layout.addStretch()
        display_scroll.setWidget(display_inner)
        self.sidebar_tabs.addTab(display_scroll, "Display")

        # ─── v5.1.0 item 25: Results sidebar tab ────────────────────
        from node_runner.dialogs.results_tab import ResultsTab
        self.results_tab = ResultsTab(self)
        # v5.2.0 item 42: new Post-Processing Toolbox signal API.
        # ----- Selection -----
        self.results_tab.output_set_changed.connect(
            self._on_results_output_set_changed)
        # v5.3.2 item 67: single Output Vector dropdown. Drives both
        # the contour render AND the deformation (deformation always
        # uses the displacement triplet from the same subcase).
        self.results_tab.output_vector_changed.connect(
            self._on_results_output_vector_changed)
        # ----- Deform -----
        self.results_tab.deform_style_changed.connect(
            self._on_results_deform_style_changed)
        self.results_tab.deform_scale_changed.connect(
            self._on_results_deform_scale_changed)
        self.results_tab.animation_settings_changed.connect(
            self._on_results_animation_settings_changed)
        self.results_tab.show_undeformed_ref_changed.connect(
            self._on_results_show_undeformed_ref_changed)
        # ----- Contour -----
        self.results_tab.contour_style_changed.connect(
            self._on_results_contour_style_changed)
        self.results_tab.data_conversion_changed.connect(
            self._on_results_data_conversion_changed)
        self.results_tab.levels_changed.connect(
            self._on_results_levels_changed)
        self.results_tab.palette_changed.connect(
            self._on_results_palette_changed)
        self.results_tab.color_range_changed.connect(
            self._on_results_color_range_changed)
        self.results_tab.markers_changed.connect(
            self._on_results_markers_changed)
        self.results_tab.element_labels_changed.connect(
            self._on_results_element_labels_changed)
        self.sidebar_tabs.addTab(self.results_tab, "Results")
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'built', version='v5.2.0',
                       n_sections=len(self.results_tab.sections))
        except Exception:
            pass

        self.plotter = QtInteractor(self)
        self.plotter.set_background('#1e1e2e')
        self.plotter.installEventFilter(self)  # Catch ESC for presentation mode

        # v4.0.1 (Stage 1): register the profiling status sink so
        # `with perf_stage(..., status_msg=...)` calls land in the
        # import dialog / status bar.
        try:
            from node_runner.profiling import set_status_callback
            set_status_callback(self._emit_render_progress)
        except Exception:
            pass

        # v4.0.0 (A2): install a 3-light rig so phong shading produces
        # a visible delta vs flat. PyVista's default is a single
        # headlight, which on edge-overlaid meshes gives a near-flat
        # appearance regardless of interpolation mode.
        self._install_light_rig()

        # v4.0.0 (A2): wrap plotter.add_mesh so every new actor inherits
        # the user's current shading state. Pre-v4.0.0, actors created
        # after a Shading toggle kept the default phong state and looked
        # out of sync with the rest of the scene.
        _orig_add_mesh = self.plotter.add_mesh

        def _add_mesh_with_shading(*args, **kw):
            actor = _orig_add_mesh(*args, **kw)
            try:
                self._apply_shading_to_actor(actor)
            except Exception:
                pass
            return actor

        self.plotter.add_mesh = _add_mesh_with_shading

        # v4.0.9: wrap plotter.render() so every frame-present cost is
        # observable. Separates VTK render time from our Python work so
        # we can tell whether per-click latency is GPU-bound or CPU-bound.
        # Gated by the perf_stage emit logic on NR_PROFILE=1.
        _orig_render = self.plotter.render

        def _render_with_timing(*args, **kw):
            try:
                from node_runner.profiling import perf_stage
                with perf_stage('render', 'present'):
                    return _orig_render(*args, **kw)
            except Exception:
                return _orig_render(*args, **kw)

        self.plotter.render = _render_with_timing

        # Phase 2.3: Render debouncer. _request_render() coalesces a burst
        # of render() calls into one frame on a 16ms timer. Existing direct
        # self.plotter.render() calls keep working unchanged; new code can
        # opt in by calling self._request_render() instead.
        self._render_timer = QtCore.QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(16)
        self._render_timer.timeout.connect(lambda: self.plotter.render())

        # Cross-section / clipping plane controller. Inactive by default;
        # View > Cross Section opens the controls.
        from node_runner.cross_section import CrossSectionController
        self._cross_section = CrossSectionController(self.plotter)
        # Replaced the old non-modal dialog with a dockable widget that
        # tabifies alongside Results / Animation / Vectors. Build it
        # lazily after the result browser dock so we can tabify.
        self._cross_section_dock = None
        # v5.0.0 item 3: floating Define-Plane panel replaces the dock.
        self._clipping_panel = None
        # v5.0.0 items 17/18: MYSTRAN solver worker + progress dialog
        # references kept here so they can be GC'd cleanly after a run.
        self._solver_worker = None
        self._solver_dialog = None
        # v5.1.0 item 25: Results sidebar tab state. Defaults match the
        # widget defaults so initial behaviour is "auto range, no
        # markers, no element labels".
        self._results_user_auto_range = True
        self._results_user_vmin = 0.0
        self._results_user_vmax = 1.0
        self._results_n_levels = 9
        self._results_show_min_marker = False
        self._results_show_max_marker = False
        self._results_show_element_labels = False
        self._results_element_labels_top_n = 0
        # v5.3.2 item 67: single Output Vector state. Drives both the
        # contour render and the deformation (deformation reads the
        # displacement triplet from the active subcase).
        self._output_vector_kind = 'displacement'
        self._output_vector_component = 'Magnitude'
        self._output_vector_mode_idx = -1
        # v5.3.2 item 69: scalar bar position lives in QSettings now,
        # set via File > Preferences. Read on each render via
        # ``_build_scalar_bar_kw``.
        from PySide6.QtCore import QSettings
        try:
            _s = QSettings("NodeRunner", "NodeRunner")
            self._results_bar_orientation = str(
                _s.value("results/bar_orientation", "right", type=str))
        except Exception:
            self._results_bar_orientation = "right"
        # v5.2.0 item 42: Post-Processing Toolbox state.
        # v5.3.3: defaults match the a professional FEA tool "None - Model Only" model.
        # Deform stays undeformed and contour stays off until the user
        # explicitly engages either. This stops a Deform-Style toggle
        # from also painting a contour the user never asked for.
        self._results_deform_style = 'undeformed'      # undeformed|deformed|animate
        self._results_scale_mode = 'pct'               # pct | actual
        self._results_scale_value = 10.0
        self._results_show_undeformed_ref = False
        self._results_contour_style = 'none'           # none|filled|filled_edges|bands
        self._results_data_conversion = 'average'      # average|no_avg|max_node|min_node
        self._results_palette = 'jet'
        # v5.2.0 item 44: animation state extensions.
        self._anim_mode = 'sine'                       # 'sine' or 'modes'
        self._anim_n_frames = 20
        self._anim_delay_ms = 40
        self._anim_through_modes_current = 0
        # v5.0.0 item 3: floating Define-Plane panel replaces the dock.
        self._clipping_panel = None
        # v5.0.0 items 17/18: MYSTRAN solver worker + progress dialog
        # references kept here so they can be GC'd cleanly after a run.
        self._solver_worker = None
        self._solver_dialog = None
        # v5.1.0 item 25: Results sidebar tab state. Defaults match the
        # widget defaults so initial behaviour is "auto range, no
        # markers, no element labels".
        self._results_user_auto_range = True
        self._results_user_vmin = 0.0
        self._results_user_vmax = 1.0
        self._results_n_levels = 9
        self._results_show_min_marker = False
        self._results_show_max_marker = False
        self._results_show_element_labels = False
        self._results_element_labels_top_n = 0

        # Selection overlay (VTK 2D actors for box/circle/polygon feedback)
        self._selection_overlay = SelectionOverlay(self.plotter)

        self.plotter.add_axes()
        self.axes_actor = self.plotter.renderer.axes_actor
        self._set_axes_label_color((0.804, 0.839, 0.957))

        # Apply persisted background settings
        self._apply_background()

        # Selection bar - floating dialog, created once and shown/hidden per workflow
        from node_runner.dialogs.selection import EntitySelectionBar
        self.selection_bar = EntitySelectionBar(self)

        splitter = QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_panel_container)
        splitter.addWidget(self.plotter.interactor)
        splitter.setSizes([300, 900])
        main_layout.addWidget(splitter)

        # --- Standard-view quick-access buttons (top-left of 3D viewer) ---
        self._view_bar = QWidget(self.plotter)
        self._view_bar.setObjectName("view_bar")
        self._view_bar.setStyleSheet("""
            QWidget#view_bar {
                background: rgba(30, 30, 46, 180);
                border-radius: 4px;
                border: 1px solid rgba(137, 180, 250, 60);
            }
            QPushButton {
                background: transparent;
                border: 1px solid transparent; border-radius: 3px;
                padding: 2px;
            }
            QPushButton:hover {
                background: rgba(137, 180, 250, 60);
                border-color: rgba(137, 180, 250, 120);
            }
            QPushButton:pressed { background: rgba(137, 180, 250, 100); }
        """)
        vb_layout = QHBoxLayout(self._view_bar)
        vb_layout.setContentsMargins(4, 4, 4, 4)
        vb_layout.setSpacing(2)
        icon_sz = 24
        view_defs = [
            ("top",    "Top (+Z)",    lambda: self._set_standard_view('top')),
            ("bottom", "Bottom (-Z)", lambda: self._set_standard_view('bottom')),
            ("front",  "Front (+X)",  lambda: self._set_standard_view('front')),
            ("back",   "Back (-X)",   lambda: self._set_standard_view('back')),
            ("right",  "Right (+Y)",  lambda: self._set_standard_view('right')),
            ("left",   "Left (-Y)",   lambda: self._set_standard_view('left')),
            ("iso",    "Isometric",   lambda: (self.plotter.view_isometric(), self._update_status("View set to: Isometric"))),
        ]
        for i, (vname, tooltip, cb) in enumerate(view_defs):
            # Thin separator between paired groups: T/Bo | F/Bk | R/L | Iso
            if i in (2, 4, 6):
                sep = QFrame()
                sep.setFrameShape(QFrame.VLine)
                sep.setFixedWidth(1)
                sep.setStyleSheet("color: rgba(137, 180, 250, 40);")
                vb_layout.addWidget(sep)
            btn = QPushButton()
            btn.setIcon(self._make_view_cube_icon(vname, icon_sz))
            btn.setIconSize(QtCore.QSize(icon_sz, icon_sz))
            btn.setFixedSize(icon_sz + 6, icon_sz + 6)
            btn.setToolTip(tooltip)
            btn.setFocusPolicy(QtCore.Qt.NoFocus)
            btn.clicked.connect(cb)
            vb_layout.addWidget(btn)
        self._view_bar.adjustSize()
        self._view_bar.move(8, 8)
        self._view_bar.show()

        # --- Debounced glyph refresh after camera changes ---
        self._last_visibility_mask = None
        self._glyph_refresh_timer = QtCore.QTimer(self)
        self._glyph_refresh_timer.setSingleShot(True)
        self._glyph_refresh_timer.setInterval(200)
        self._glyph_refresh_timer.timeout.connect(self._refresh_orientation_glyphs)
        iren = self.plotter.interactor
        iren.AddObserver(vtk.vtkCommand.MouseWheelForwardEvent, lambda *_: self._schedule_glyph_refresh())
        iren.AddObserver(vtk.vtkCommand.MouseWheelBackwardEvent, lambda *_: self._schedule_glyph_refresh())
        iren.AddObserver(vtk.vtkCommand.EndInteractionEvent, lambda *_: self._schedule_glyph_refresh())
        # Theme B1: probe hover. We throttle to ~10 Hz via a lightweight
        # timer so we don't pick on every mouse-move event (60+ Hz).
        self._probe_throttle = QtCore.QTimer(self)
        self._probe_throttle.setSingleShot(True)
        self._probe_throttle.setInterval(80)
        self._probe_throttle.timeout.connect(self._do_probe_pick)
        iren.AddObserver(vtk.vtkCommand.MouseMoveEvent,
                         lambda *_: self._on_probe_mouse_move())

        self.statusBar().showMessage("Ready.")

        # --- Selection bar signal connections (connected once, persist for lifetime) ---
        self._selection_accept_callback = None
        self.selection_bar.request_show_selection.connect(self._highlight_entities)
        self.selection_bar.request_picking_mode.connect(self._set_picking_mode)
        self.selection_bar.request_advanced_selection.connect(self._on_advanced_selection)
        self.selection_bar.request_selection_mode.connect(self._on_selection_mode_change)
        self.selection_bar.request_previous_selection.connect(self._on_previous_selection)
        self.selection_bar.accepted.connect(self._on_selection_bar_accepted)
        self.selection_bar.rejected.connect(self._on_selection_bar_rejected)

        self.theme_button.clicked.connect(self._toggle_theme)
        # v4.0.4 (Stage R): use Qt.UniqueConnection so a duplicate
        # connect attempt RAISES instead of silently creating a
        # second connection (which would mute the first slot in
        # some Qt versions). Wrap in try/except so we don't crash if
        # something already connected this slot earlier.
        try:
            self.tree_widget.itemChanged.connect(
                self._handle_tree_item_changed, QtCore.Qt.UniqueConnection)
        except Exception as _exc:
            try:
                from node_runner.profiling import perf_event
                perf_event('tree', 'unique_connect_failed',
                           signal='itemChanged',
                           exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            except Exception:
                pass
            # Fall back to a normal connect (may silently double-bind).
            self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)
        # v4.0.10 hotfix: wire the loads-tab tree's itemChanged signal.
        # Pre-v4.0.10 this connection didn't exist, so clicking a
        # load_set / constraint_set checkbox visually toggled the
        # checkbox but fired NO Python handler. no visibility update.
        # Confirmed via v4.0.9 profile log
        # (profile_20260513_135348.log): zero fast_path.load_set_toggle
        # events despite the user clicking load_sets. v4.0.5's
        # default-OFF flip surfaced the bug; v4.0.9's fast path made
        # it loud.
        try:
            self.loads_tab_tree.itemChanged.connect(
                self._handle_tree_item_changed, QtCore.Qt.UniqueConnection)
        except Exception as _exc:
            try:
                from node_runner.profiling import perf_event
                perf_event('tree', 'unique_connect_failed',
                           signal='loads_tab.itemChanged',
                           exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            except Exception:
                pass
            self.loads_tab_tree.itemChanged.connect(self._handle_tree_item_changed)
        try:
            from node_runner.profiling import perf_event
            perf_event('tree', 'signal_connected',
                       model_tree_itemChanged=True,
                       loads_tab_itemChanged=True)
        except Exception:
            pass
        # v4.0.2 (Stage E): also wire itemClicked. This fires for EVERY
        # click on the tree regardless of whether the checkbox state
        # changed. Gives us:
        #   1. A diagnostic event in the profile log (we can tell
        #      whether a "I clicked but nothing happened" is a tree-
        #      handler issue vs. a click-on-row-text-not-checkbox issue).
        #   2. The chance to auto-toggle the checkbox when the user
        #      clicks anywhere on a category row.
        try:
            self.tree_widget.itemClicked.connect(
                self._handle_tree_item_clicked, QtCore.Qt.UniqueConnection)
        except Exception as _exc:
            try:
                from node_runner.profiling import perf_event
                perf_event('tree', 'unique_connect_failed',
                           signal='itemClicked',
                           exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            except Exception:
                pass
            self.tree_widget.itemClicked.connect(self._handle_tree_item_clicked)

        # v4.0.4 (Stage R): alternative signal connections. If the
        # signal pipeline is partly broken, ONE of these might still
        # fire and tell us where the break is. Each handler just emits
        # a perf_event.
        try:
            self.tree_widget.itemPressed.connect(
                self._handle_tree_item_pressed)
            self.tree_widget.itemSelectionChanged.connect(
                self._handle_tree_selection_changed)
        except Exception:
            pass

        # v4.0.3 (Stage K) / v4.0.5 (Stage V): event filter on
        # tree_widget only. Pre-v4.0.5 we ALSO installed it on the
        # viewport. The v4.0.4 log showed Qt's `itemPressed` /
        # `itemClicked` / `itemChanged` signals all dead while
        # `itemSelectionChanged` worked. a known PySide6 pattern
        # for "event filter on QAbstractItemView's viewport breaks
        # the view's mousePressEvent dispatch". v4.0.5 bisects by
        # removing the viewport filter.
        # Gated on NR_PROFILE=1; zero overhead in normal runs.
        try:
            from node_runner.profiling import enabled as _profile_on, perf_event
            if _profile_on():
                self.tree_widget.installEventFilter(self)
                # v4.0.5 (Stage V): viewport filter intentionally
                # NOT installed. If tree item signals fire after this
                # change, the viewport filter was the cause.
                perf_event('tree', 'viewport_filter_disabled',
                           reason='v4.0.5_bisect_test')
                perf_event('tree', 'signal_connected',
                           itemChanged=True,
                           itemClicked=True,
                           event_filter_installed=True,
                           viewport_filtered=False)
        except Exception:
            pass
        self.show_node_labels_check.toggled.connect(self._toggle_node_labels)
        self.show_elem_labels_check.toggled.connect(self._toggle_element_labels)

        # --- NEW: Connect scaling controls to the redraw function ---
        self.arrow_scale_input.textChanged.connect(self._update_plot_visibility)
        self.relative_scaling_check.stateChanged.connect(self._update_plot_visibility)

        #find_button.clicked.connect(self._find_and_zoom_to_entity)
        #clear_button.clicked.connect(self._clear_entity_highlight)


        self.show_normals_check.stateChanged.connect(self._toggle_element_normals_visibility)
        self.normal_arrow_scale_input.textChanged.connect(self._toggle_element_normals_visibility)
        self.show_free_edges_check.stateChanged.connect(self._toggle_free_edges_visibility)

        self.quality_metric_combo.currentIndexChanged.connect(self._update_plot_visibility)
        self.show_csys_check.stateChanged.connect(self._refresh_coord_actors)
        self.show_global_csys_check.stateChanged.connect(self._refresh_coord_actors)
        self.show_csys_labels_check.stateChanged.connect(self._refresh_coord_actors)
        self.csys_scale_slider.valueChanged.connect(self._refresh_coord_actors)
        self.show_bush_orient_check.stateChanged.connect(self._update_plot_visibility)
        self.show_beam_orient_check.stateChanged.connect(self._update_plot_visibility)
        self.show_beam_sections_check.stateChanged.connect(self._update_plot_visibility)
    
# START: Replacement for _create_new_model in main.py
    def _undo(self):
        model = self.current_generator.model if self.current_generator else None
        if not model or not self.command_manager.can_undo:
            self._update_status("Nothing to undo.")
            return
        desc = self.command_manager.undo_description
        self.command_manager.undo(model)
        self._update_viewer(self.current_generator, reset_camera=False)
        self._update_status(f"Undo: {desc}")

    def _redo(self):
        model = self.current_generator.model if self.current_generator else None
        if not model or not self.command_manager.can_redo:
            self._update_status("Nothing to redo.")
            return
        desc = self.command_manager.redo_description
        self.command_manager.redo(model)
        self._update_viewer(self.current_generator, reset_camera=False)
        self._update_status(f"Redo: {desc}")

    def _create_new_model(self):
        """Creates a new, empty model with default global coordinate systems."""
        self.current_generator = NastranModelGenerator()
        # The helper function now ensures the default coordinate systems are created.
        self._ensure_default_coords(self.current_generator.model)
        self.command_manager.clear()
        self.subcases = []
        self.geometry_store.clear()
        self.sol_type = None
        self.eigrl_cards = []
        self.analysis_sets.clear()
        self.active_analysis_set_id = None
        self.op2_results = None
        self.results_widget.setVisible(False)
        # Phase 8: Reset groups on new model
        self.groups.clear()
        self._hidden_groups.clear()
        self._isolate_mode = None
        self._populate_groups_list()
        self._update_viewer(self.current_generator)
        self._update_status("New model created.")
# END: Replacement for _create_new_model in main.py


    # START: New helper method in main.py
    def _ensure_default_coords(self, model):
        """
        Ensures the three default coordinate systems exist in the given model.
        This will overwrite any user-defined cards at CIDs 1 and 2.
        """
        # pyNastran's BDF() object creates CID 0 (Global Rectangular) by default.
        # We only need to ensure CIDs 1 and 2 are correct.
        model.add_cord2c(1, rid=0, origin=[0.,0.,0.], zaxis=[0.,0.,1.], xzplane=[1.,0.,0.], comment='$ Global Cylindrical')
        model.add_cord2s(2, rid=0, origin=[0.,0.,0.], zaxis=[0.,0.,1.], xzplane=[1.,0.,0.], comment='$ Global Spherical')
# END: New helper method in main.py



    def _set_axes_label_color(self, color_tuple):
        if not hasattr(self, 'axes_actor') or not self.axes_actor:
            return
        self.axes_actor.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)
        self.axes_actor.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)
        self.axes_actor.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)

    def _create_menu_bar(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("&File")
        
        new_action = QAction("&New Model", self)
        new_action.triggered.connect(self._create_new_model)
        file_menu.addAction(new_action)
        file_menu.addSeparator()

        open_action = QAction("&Open BDF/DAT...", self)
        open_action.triggered.connect(self._open_file_dialog)
        file_menu.addAction(open_action)
        
        # --- NEW: Add the Import menu action ---
        import_action = QAction("&Import BDF/DAT...", self)
        import_action.triggered.connect(self._import_file_dialog)
        file_menu.addAction(import_action)
        import_cad_action = QAction("Import CAD (STEP/IGES/STL)...", self)
        import_cad_action.triggered.connect(self._import_cad_file)
        file_menu.addAction(import_cad_action)

        save_action = QAction("&Save BDF File...", self)
        save_action.triggered.connect(self._save_file_dialog)
        file_menu.addAction(save_action)

        file_menu.addSeparator()
        load_results_action = QAction("Load Results (&OP2)...", self)
        load_results_action.triggered.connect(self._load_results_dialog)
        file_menu.addAction(load_results_action)
        # v5.2.2 issue 5: detach the current results bundle without
        # unloading the model.
        detach_results_action = QAction("Detach Results", self)
        detach_results_action.setStatusTip(
            "Clear the loaded results bundle but keep the model open.")
        detach_results_action.triggered.connect(self._detach_results)
        file_menu.addAction(detach_results_action)

        file_menu.addSeparator()
        save_conf_action = QAction("Save Fuselage Config...", self)
        save_conf_action.triggered.connect(self._save_configuration)
        file_menu.addAction(save_conf_action)
        load_conf_action = QAction("Load Fuselage Config...", self)
        load_conf_action.triggered.connect(self._load_configuration)
        file_menu.addAction(load_conf_action)
        file_menu.addSeparator()
        quit_action = QAction("&Quit", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        edit_menu = menu_bar.addMenu("&Edit")
        undo_action = QAction("&Undo", self)
        undo_action.setShortcut("Ctrl+Z")
        undo_action.triggered.connect(self._undo)
        edit_menu.addAction(undo_action)
        redo_action = QAction("&Redo", self)
        redo_action.setShortcut("Ctrl+Y")
        redo_action.triggered.connect(self._redo)
        edit_menu.addAction(redo_action)
        edit_menu.addSeparator()
        select_entities_action = QAction("Select &Entities...", self)
        select_entities_action.setShortcut("Ctrl+F")
        select_entities_action.triggered.connect(self._open_find_entity_tool)
        edit_menu.addAction(select_entities_action)

        # v3.5.0 item 3: Preferences dialog (colors / sizes / highlight)
        edit_menu.addSeparator()
        prefs_action = QAction("&Preferences...", self)
        prefs_action.setShortcut("Ctrl+,")
        prefs_action.triggered.connect(self._open_preferences_dialog)
        edit_menu.addAction(prefs_action)

        # --- Model menu (merged: old Model + Mesh + Edit items) ---
        model_menu = menu_bar.addMenu("&Model")

        # Section 1: Create - Nodes, Elements, Properties, Materials, Coords
        create_nodes_action = QAction("Nodes...", self)
        create_nodes_action.triggered.connect(self._create_nodes)
        model_menu.addAction(create_nodes_action)
        elements_menu = model_menu.addMenu("Elements")
        create_lines_action = QAction("Line...", self)
        create_lines_action.triggered.connect(self._create_line_elements)
        elements_menu.addAction(create_lines_action)
        create_plates_action = QAction("Plate...", self)
        create_plates_action.triggered.connect(self._create_plate_elements)
        elements_menu.addAction(create_plates_action)
        create_solid_action = QAction("Solid...", self)
        create_solid_action.triggered.connect(self._create_solid_elements)
        elements_menu.addAction(create_solid_action)
        create_bush_action = QAction("Bush...", self)
        create_bush_action.triggered.connect(self._create_bush_elements)
        elements_menu.addAction(create_bush_action)
        create_shear_action = QAction("Shear Panel...", self)
        create_shear_action.triggered.connect(self._create_shear_elements)
        elements_menu.addAction(create_shear_action)
        create_gap_action = QAction("Gap...", self)
        create_gap_action.triggered.connect(self._create_gap_elements)
        elements_menu.addAction(create_gap_action)
        create_rbes_action = QAction("Rigid...", self)
        create_rbes_action.triggered.connect(self._create_rbes)
        elements_menu.addAction(create_rbes_action)
        create_conm2_action = QAction("Concentrated Mass (CONM2)...", self)
        create_conm2_action.triggered.connect(self._create_conm2)
        elements_menu.addAction(create_conm2_action)
        create_plotel_action = QAction("Plot Element (PLOTEL)...", self)
        create_plotel_action.triggered.connect(self._create_plotels)
        elements_menu.addAction(create_plotel_action)
        create_prop_action = QAction("Property...", self)
        create_prop_action.triggered.connect(self._create_property)
        model_menu.addAction(create_prop_action)
        create_mat_action = QAction("Material...", self)
        create_mat_action.triggered.connect(self._create_material)
        model_menu.addAction(create_mat_action)
        create_coord_action = QAction("Coordinate System...", self)
        create_coord_action.triggered.connect(self._create_coord_system)
        model_menu.addAction(create_coord_action)

        # Section 2: Loads & Constraints
        model_menu.addSeparator()
        create_load_action = QAction("Load...", self)
        create_load_action.triggered.connect(self._create_load)
        model_menu.addAction(create_load_action)
        create_constraint_action = QAction("Constraint...", self)
        create_constraint_action.triggered.connect(self._create_constraint)
        model_menu.addAction(create_constraint_action)
        create_load_combo_action = QAction("Load Combination (LOAD)...", self)
        create_load_combo_action.triggered.connect(self._create_load_combination)
        model_menu.addAction(create_load_combo_action)
        # --- Theme C: additional load and BC creators ---
        create_pload1_action = QAction("Distributed Load (PLOAD1)...", self)
        create_pload1_action.triggered.connect(self._create_pload1)
        model_menu.addAction(create_pload1_action)
        create_pload2_action = QAction("Element Pressure (PLOAD2)...", self)
        create_pload2_action.triggered.connect(self._create_pload2)
        model_menu.addAction(create_pload2_action)
        create_spcd_action = QAction("Enforced Displacement (SPCD)...", self)
        create_spcd_action.triggered.connect(self._create_spcd)
        model_menu.addAction(create_spcd_action)
        create_mpc_action = QAction("Multi-Point Constraint (MPC)...", self)
        create_mpc_action.triggered.connect(self._create_mpc)
        model_menu.addAction(create_mpc_action)
        create_bolt_action = QAction("Bolt Preload (BOLT)...", self)
        create_bolt_action.triggered.connect(self._create_bolt)
        model_menu.addAction(create_bolt_action)

        # Section 3: Modify (moved from Edit menu)
        model_menu.addSeparator()
        modify_menu = model_menu.addMenu("Modify")
        move_nodes_action = QAction("Transform Nodes...", self)
        move_nodes_action.triggered.connect(self._open_node_move_tool)
        modify_menu.addAction(move_nodes_action)
        edit_element_action = QAction("Element...", self)
        edit_element_action.triggered.connect(self._open_element_editor)
        modify_menu.addAction(edit_element_action)
        edit_props_action = QAction("Properties...", self)
        edit_props_action.triggered.connect(self._open_property_editor)
        modify_menu.addAction(edit_props_action)
        edit_mats_action = QAction("Materials...", self)
        edit_mats_action.triggered.connect(self._open_material_editor)
        modify_menu.addAction(edit_mats_action)
        # --- Theme A: mesh-editing operations ---
        modify_menu.addSeparator()
        split_quad_action = QAction("Split QUAD into TRIA...", self)
        split_quad_action.triggered.connect(self._open_split_quad_tool)
        modify_menu.addAction(split_quad_action)
        refine_action = QAction("Refine Elements (1-into-4)...", self)
        refine_action.triggered.connect(self._open_refine_tool)
        modify_menu.addAction(refine_action)
        combine_action = QAction("Combine TRIA into QUAD...", self)
        combine_action.triggered.connect(self._open_combine_tool)
        modify_menu.addAction(combine_action)
        smooth_action = QAction("Smooth Nodes...", self)
        smooth_action.triggered.connect(self._open_smooth_tool)
        modify_menu.addAction(smooth_action)
        mirror_elem_action = QAction("Mirror Elements...", self)
        mirror_elem_action.triggered.connect(self._open_mirror_elements_tool)
        modify_menu.addAction(mirror_elem_action)
        copy_elem_action = QAction("Copy Elements...", self)
        copy_elem_action.triggered.connect(self._open_copy_elements_tool)
        modify_menu.addAction(copy_elem_action)
        insert_edge_node_action = QAction("Insert Node on Edge...", self)
        insert_edge_node_action.triggered.connect(self._open_insert_edge_node_tool)
        modify_menu.addAction(insert_edge_node_action)

        # Section 4: Mesh Operations (moved from Mesh menu)
        model_menu.addSeparator()
        flip_normals_action = QAction("Flip Element Normals...", self)
        flip_normals_action.triggered.connect(self._open_element_flip_tool)
        model_menu.addAction(flip_normals_action)
        orient_normals_action = QAction("Orient Normals (Auto)...", self)
        orient_normals_action.triggered.connect(self._auto_orient_normals)
        model_menu.addAction(orient_normals_action)
        split_plates_action = QAction("Split Plates...", self)
        split_plates_action.triggered.connect(self._open_element_split_tool)
        model_menu.addAction(split_plates_action)

        # Section 5: Connections (moved from Mesh menu)
        model_menu.addSeparator()
        connections_menu = model_menu.addMenu("Connections")
        spider_action = QAction("Spider Connection (RBE)...", self)
        spider_action.triggered.connect(self._create_spider_connection)
        connections_menu.addAction(spider_action)
        weld_action = QAction("Weld/Fastener (CWELD)...", self)
        weld_action.triggered.connect(self._create_weld_fastener)
        connections_menu.addAction(weld_action)
        # --- Theme C: rigid element creators ---
        connections_menu.addSeparator()
        rbar_action = QAction("Rigid Bar (RBAR)...", self)
        rbar_action.triggered.connect(self._create_rbar)
        connections_menu.addAction(rbar_action)
        rbe1_action = QAction("General Rigid (RBE1)...", self)
        rbe1_action.triggered.connect(self._create_rbe1)
        connections_menu.addAction(rbe1_action)
        rspline_action = QAction("Rigid Spline (RSPLINE)...", self)
        rspline_action.triggered.connect(self._create_rspline)
        connections_menu.addAction(rspline_action)

        # Section 6: Geometry Mesh (moved from Mesh menu)
        model_menu.addSeparator()
        geom_mesh_menu = model_menu.addMenu("Geometry Mesh")
        mesh_curve_action = QAction("Mesh Curve...", self)
        mesh_curve_action.triggered.connect(self._mesh_curve)
        geom_mesh_menu.addAction(mesh_curve_action)
        nodes_at_points_action = QAction("Nodes at Points...", self)
        nodes_at_points_action.triggered.connect(self._nodes_at_geometry_points)
        geom_mesh_menu.addAction(nodes_at_points_action)
        mesh_surface_action = QAction("Mesh Surface...", self)
        mesh_surface_action.triggered.connect(self._mesh_surface)
        geom_mesh_menu.addAction(mesh_surface_action)

        geom_menu = menu_bar.addMenu("&Geometry")
        geom_create_menu = geom_menu.addMenu("Create")
        geom_point_action = QAction("Point...", self)
        geom_point_action.triggered.connect(self._create_geometry_point)
        geom_create_menu.addAction(geom_point_action)
        geom_line_action = QAction("Line...", self)
        geom_line_action.triggered.connect(self._create_geometry_line)
        geom_create_menu.addAction(geom_line_action)
        geom_arc_action = QAction("Arc...", self)
        geom_arc_action.triggered.connect(self._create_geometry_arc)
        geom_create_menu.addAction(geom_arc_action)
        geom_circle_action = QAction("Circle...", self)
        geom_circle_action.triggered.connect(self._create_geometry_circle)
        geom_create_menu.addAction(geom_circle_action)
        geom_surface_action = QAction("Surface...", self)
        geom_surface_action.triggered.connect(self._create_geometry_surface)
        geom_create_menu.addAction(geom_surface_action)

        geom_modify_menu = geom_menu.addMenu("Modify")
        edit_point_action = QAction("Edit Point...", self)
        edit_point_action.triggered.connect(self._edit_geometry_point)
        geom_modify_menu.addAction(edit_point_action)
        edit_curve_action = QAction("Edit Curve...", self)
        edit_curve_action.triggered.connect(self._edit_geometry_curve)
        geom_modify_menu.addAction(edit_curve_action)
        edit_surface_action = QAction("Edit Surface Boundaries...", self)
        edit_surface_action.triggered.connect(self._edit_geometry_surface)
        geom_modify_menu.addAction(edit_surface_action)

        geom_transform_menu = geom_menu.addMenu("Transform")
        translate_geom_action = QAction("Translate...", self)
        translate_geom_action.triggered.connect(lambda: self._transform_geometry('translate'))
        geom_transform_menu.addAction(translate_geom_action)
        rotate_geom_action = QAction("Rotate...", self)
        rotate_geom_action.triggered.connect(lambda: self._transform_geometry('rotate'))
        geom_transform_menu.addAction(rotate_geom_action)
        mirror_geom_action = QAction("Mirror...", self)
        mirror_geom_action.triggered.connect(lambda: self._transform_geometry('mirror'))
        geom_transform_menu.addAction(mirror_geom_action)
        scale_geom_action = QAction("Scale...", self)
        scale_geom_action.triggered.connect(lambda: self._transform_geometry('scale'))
        geom_transform_menu.addAction(scale_geom_action)
        geom_transform_menu.addSeparator()
        copy_geom_action = QAction("Copy...", self)
        copy_geom_action.triggered.connect(self._copy_geometry)
        geom_transform_menu.addAction(copy_geom_action)

        geom_ops_menu = geom_menu.addMenu("Operations")
        split_curve_action = QAction("Split Curve...", self)
        split_curve_action.triggered.connect(self._split_curve)
        geom_ops_menu.addAction(split_curve_action)
        offset_curve_action = QAction("Offset Curve...", self)
        offset_curve_action.triggered.connect(self._offset_curve)
        geom_ops_menu.addAction(offset_curve_action)
        project_pt_action = QAction("Project Point to Curve...", self)
        project_pt_action.triggered.connect(self._project_point_to_curve)
        geom_ops_menu.addAction(project_pt_action)
        fillet_action = QAction("Fillet...", self)
        fillet_action.triggered.connect(self._fillet_curves)
        geom_ops_menu.addAction(fillet_action)

        delete_menu = menu_bar.addMenu("&Delete")
        delete_nodes_action = QAction("Nodes...", self)
        delete_nodes_action.triggered.connect(self._delete_nodes)
        delete_menu.addAction(delete_nodes_action)
        delete_elements_action = QAction("Elements...", self)
        delete_elements_action.triggered.connect(self._delete_elements)
        delete_menu.addAction(delete_elements_action)
        delete_coords_action = QAction("Coordinate Systems...", self)
        delete_coords_action.triggered.connect(self._delete_coord_system_menu)
        delete_menu.addAction(delete_coords_action)
        delete_menu.addSeparator()
        delete_geom_pts_action = QAction("Geometry Points...", self)
        delete_geom_pts_action.triggered.connect(self._delete_geometry_points)
        delete_menu.addAction(delete_geom_pts_action)
        delete_geom_curves_action = QAction("Geometry Curves...", self)
        delete_geom_curves_action.triggered.connect(self._delete_geometry_curves)
        delete_menu.addAction(delete_geom_curves_action)
        delete_geom_surfs_action = QAction("Geometry Surfaces...", self)
        delete_geom_surfs_action.triggered.connect(self._delete_geometry_surfaces)
        delete_menu.addAction(delete_geom_surfs_action)

        list_menu = menu_bar.addMenu("&List")
        list_node_action = QAction("Node Information...", self)
        list_node_action.triggered.connect(self._list_node_info)
        list_menu.addAction(list_node_action)
        list_element_action = QAction("Element Information...", self)
        list_element_action.triggered.connect(self._list_element_info)
        list_menu.addAction(list_element_action)
        list_menu.addSeparator()
        model_summary_action = QAction("Model Summary...", self)
        model_summary_action.triggered.connect(self._show_model_summary)
        list_menu.addAction(model_summary_action)
        mass_summary_action = QAction("Mass && CG Summary...", self)
        mass_summary_action.triggered.connect(self._show_mass_summary)
        list_menu.addAction(mass_summary_action)

        analysis_menu = menu_bar.addMenu("&Analysis")
        analysis_set_mgr_action = QAction("Analysis Set Manager...", self)
        analysis_set_mgr_action.triggered.connect(self._open_analysis_set_manager)
        analysis_menu.addAction(analysis_set_mgr_action)
        analysis_menu.addSeparator()
        sol_menu = analysis_menu.addMenu("Solution Type")
        try:
            sol_menu.setToolTipsVisible(True)
        except Exception:
            pass
        from node_runner.solve import mystran_tips as _stips
        sol_none_action = QAction("None (Punch)", self, checkable=True, checked=True)
        sol_none_action.setToolTip(
            "<b>None (Punch)</b><br>No SOL number written -- the deck "
            "exports as a punch fragment. Use only when handing the "
            "deck to a solver other than MSC / MYSTRAN.")
        sol_101_action = QAction("SOL 101 - Linear Statics", self, checkable=True)
        sol_101_action.setToolTip(_stips.SOL_101)
        sol_103_action = QAction("SOL 103 - Normal Modes", self, checkable=True)
        sol_103_action.setToolTip(_stips.SOL_103)
        sol_105_action = QAction("SOL 105 - Buckling", self, checkable=True)
        sol_105_action.setToolTip(_stips.SOL_105)
        sol_106_action = QAction("SOL 106 - Nonlinear Static", self, checkable=True)
        sol_106_action.setToolTip(
            "<b>SOL 106. Nonlinear Static</b><br>"
            "Iterative nonlinear (large-deformation, material "
            "nonlinearity, contact). <b>NOT supported by MYSTRAN.</b> "
            "Selecting this will cause MYSTRAN pre-flight to block "
            "the run. Use for MSC / NX export only.")
        sol_group = QActionGroup(self)
        for action in [sol_none_action, sol_101_action, sol_103_action, sol_105_action, sol_106_action]:
            sol_group.addAction(action)
            sol_menu.addAction(action)
        sol_none_action.triggered.connect(lambda: self._set_sol_type(None))
        sol_101_action.triggered.connect(lambda: self._set_sol_type(101))
        sol_103_action.triggered.connect(lambda: self._set_sol_type(103))
        sol_105_action.triggered.connect(lambda: self._set_sol_type(105))
        sol_106_action.triggered.connect(lambda: self._set_sol_type(106))
        analysis_menu.addSeparator()
        subcase_action = QAction("Subcases (Case Control)...", self)
        subcase_action.triggered.connect(self._open_subcase_editor)
        analysis_menu.addAction(subcase_action)
        eigrl_action = QAction("EIGRL (Eigenvalue Extraction)...", self)
        eigrl_action.triggered.connect(self._create_eigrl)
        analysis_menu.addAction(eigrl_action)
        analysis_menu.addSeparator()
        output_req_action = QAction("Output Requests...", self)
        output_req_action.triggered.connect(self._open_output_requests_dialog)
        analysis_menu.addAction(output_req_action)

        # ─── v5.0.0 items 17 + 18: MYSTRAN solver integration ───
        analysis_menu.addSeparator()
        from node_runner.solve import mystran_tips as _mtips
        run_mystran_action = QAction("Run Analysis (MYSTRAN)…", self)
        run_mystran_action.setShortcut("Ctrl+R")
        run_mystran_action.setStatusTip(
            "Pre-flight + solve via the open-source MYSTRAN binary.")
        run_mystran_action.setToolTip(_mtips.ANALYSIS_RUN_MENU)
        run_mystran_action.triggered.connect(self._run_mystran)
        analysis_menu.addAction(run_mystran_action)
        configure_mystran_action = QAction("Configure MYSTRAN…", self)
        configure_mystran_action.setStatusTip(
            "Open Preferences -> MYSTRAN tab.")
        configure_mystran_action.setToolTip(_mtips.ANALYSIS_CONFIGURE_MENU)
        configure_mystran_action.triggered.connect(
            self._configure_mystran)
        analysis_menu.addAction(configure_mystran_action)
        analysis_history_action = QAction("Analysis History…", self)
        analysis_history_action.setStatusTip(
            "Browse and reload prior MYSTRAN runs.")
        analysis_history_action.setToolTip(_mtips.ANALYSIS_HISTORY_MENU)
        analysis_history_action.triggered.connect(
            self._open_analysis_history)
        analysis_menu.addAction(analysis_history_action)
        # Qt menus suppress tooltips by default; enable them.
        try:
            analysis_menu.setToolTipsVisible(True)
        except Exception:
            pass

        tools_menu = menu_bar.addMenu("&Tools")
        fuselage_gen_action = QAction("Fuselage Generator...", self)
        fuselage_gen_action.triggered.connect(self._open_fuselage_generator)
        tools_menu.addAction(fuselage_gen_action)
        tools_menu.addSeparator()
        coincident_node_action = QAction("Coincident Node Checker...", self)
        coincident_node_action.triggered.connect(self._open_coincident_node_checker)
        tools_menu.addAction(coincident_node_action)
        duplicate_elem_action = QAction("Duplicate Element Checker...", self)
        duplicate_elem_action.triggered.connect(self._open_duplicate_element_checker)
        tools_menu.addAction(duplicate_elem_action)
        tools_menu.addSeparator()
        measure_dist_action = QAction("Measure Distance...", self)
        measure_dist_action.triggered.connect(self._measure_distance)
        tools_menu.addAction(measure_dist_action)
        measure_angle_action = QAction("Measure Angle...", self)
        measure_angle_action.triggered.connect(self._measure_angle)
        tools_menu.addAction(measure_angle_action)
        tools_menu.addSeparator()
        # v5.2.0 item 45: Data Table dialog (professional post-processing
        # data table that live-syncs with the viewport selection).
        data_table_action = QAction("Data Table…", self)
        data_table_action.setStatusTip(
            "Open the post-processing data table. Live-syncs with "
            "viewport selection.")
        data_table_action.triggered.connect(self._open_data_table_dialog)
        tools_menu.addAction(data_table_action)
        tools_menu.addSeparator()
        find_replace_action = QAction("Find/Replace Property/Material...", self)
        find_replace_action.triggered.connect(self._open_find_replace_dialog)
        tools_menu.addAction(find_replace_action)
        renumber_action = QAction("Renumber Nodes/Elements...", self)
        renumber_action.triggered.connect(self._open_renumber_dialog)
        tools_menu.addAction(renumber_action)
        tools_menu.addSeparator()
        fbd_action = QAction("Free Body Diagram...", self)
        fbd_action.triggered.connect(self._show_free_body_diagram)
        tools_menu.addAction(fbd_action)
        # Theme B: Probe tool. Hovering over nodes/elements while this is
        # checked surfaces a tooltip with their results.
        self.probe_mode_action = QAction("Probe Mode", self, checkable=True)
        self.probe_mode_action.toggled.connect(self._toggle_probe_mode)
        tools_menu.addAction(self.probe_mode_action)
        # Command palette - same shortcut as install_command_palette_shortcut,
        # but having a menu entry makes it discoverable and provides a
        # second way to fire it if the global shortcut gets eaten by a
        # focused widget.
        tools_menu.addSeparator()
        self.command_palette_action = QAction("Command Palette...", self)
        self.command_palette_action.setShortcut("Ctrl+P")
        self.command_palette_action.triggered.connect(self._open_command_palette)
        tools_menu.addAction(self.command_palette_action)
        tools_menu.addSeparator()
        # professional unit conversion - the tool itself is unitless, but
        # this rescales the whole model by user-supplied factors.
        convert_units_action = QAction("Convert Units...", self)
        convert_units_action.triggered.connect(self._open_convert_units_dialog)
        tools_menu.addAction(convert_units_action)
        tools_menu.addSeparator()
        model_check_menu = tools_menu.addMenu("Model Check")
        mass_props_action = QAction("Mass Properties...", self)
        mass_props_action.triggered.connect(self._show_mass_properties)
        model_check_menu.addAction(mass_props_action)
        free_edge_action = QAction("Free Edge Report...", self)
        free_edge_action.triggered.connect(self._show_free_edge_report)
        model_check_menu.addAction(free_edge_action)
        quality_summary_action = QAction("Quality Summary...", self)
        quality_summary_action.triggered.connect(self._show_quality_summary)
        model_check_menu.addAction(quality_summary_action)
        orphan_check_action = QAction("Orphan Check...", self)
        orphan_check_action.triggered.connect(self._show_orphan_check)
        model_check_menu.addAction(orphan_check_action)
        model_check_menu.addSeparator()
        orient_normals_action2 = QAction("Orient Normals (Auto)...", self)
        orient_normals_action2.triggered.connect(self._auto_orient_normals)
        model_check_menu.addAction(orient_normals_action2)

        groups_menu = menu_bar.addMenu("&Groups")
        groups_select_action = QAction("Select Entities...", self)
        groups_select_action.triggered.connect(self._open_find_entity_tool)
        groups_menu.addAction(groups_select_action)
        groups_menu.addSeparator()
        create_group_action = QAction("Create Group...", self)
        create_group_action.triggered.connect(self._create_group)
        groups_menu.addAction(create_group_action)
        delete_group_action = QAction("Delete Group", self)
        delete_group_action.triggered.connect(self._delete_group)
        groups_menu.addAction(delete_group_action)
        add_to_group_action = QAction("Add Selected to Group", self)
        add_to_group_action.triggered.connect(self._add_selected_to_group)
        groups_menu.addAction(add_to_group_action)
        isolate_group_action = QAction("Isolate Group", self)
        isolate_group_action.triggered.connect(self._isolate_group)
        groups_menu.addAction(isolate_group_action)
        groups_menu.addSeparator()
        group_by_prop_action = QAction("Group by Property", self)
        group_by_prop_action.triggered.connect(self._group_by_property)
        groups_menu.addAction(group_by_prop_action)

        # ─── v5.0.0 item 2: professional View menu reorganization ───
        view_menu = menu_bar.addMenu("&View")

        find_entities_action = QAction("Find Entities...", self)
        find_entities_action.setShortcut("Ctrl+F")
        find_entities_action.triggered.connect(self._open_find_entity_tool)
        view_menu.addAction(find_entities_action)
        view_menu.addSeparator()

        # --- Orientation submenu (was "Zoom" + standalone "Standard Views") ---
        orientation_menu = view_menu.addMenu("Orientation")
        views_menu = orientation_menu.addMenu("Standard Views")
        for label, key in [("Top (+Z)", 'top'), ("Bottom (-Z)", 'bottom'),
                           ("Front (+X)", 'front'), ("Back (-X)", 'back'),
                           ("Right (+Y)", 'right'), ("Left (-Y)", 'left')]:
            act = QAction(label, self)
            act.triggered.connect(
                lambda _checked=False, k=key: self._set_standard_view(k))
            views_menu.addAction(act)
        orientation_menu.addSeparator()
        zoom_model_action = QAction("Fit to Model", self)
        zoom_model_action.setShortcut("Ctrl+0")
        zoom_model_action.triggered.connect(self._zoom_to_model)
        orientation_menu.addAction(zoom_model_action)
        zoom_all_action = QAction("Fit to All", self)
        zoom_all_action.triggered.connect(self._zoom_to_all)
        orientation_menu.addAction(zoom_all_action)

        # --- Camera submenu (Perspective / Parallel radio pair) ---
        camera_menu = view_menu.addMenu("Camera")
        projection_group = QActionGroup(self)
        projection_group.setExclusive(True)
        self.perspective_action = QAction(
            "Perspective", self, checkable=True, checked=True)
        self.perspective_action.toggled.connect(self._toggle_perspective_view)
        projection_group.addAction(self.perspective_action)
        camera_menu.addAction(self.perspective_action)
        self.parallel_action = QAction(
            "Parallel (Orthographic)", self, checkable=True)
        # QActionGroup(exclusive=True) auto-unchecks perspective when
        # parallel becomes checked, which fires perspective_action.toggled
        # (-> _toggle_perspective_view) with False - that flips the
        # camera to parallel projection. Nothing more needed here.
        projection_group.addAction(self.parallel_action)
        camera_menu.addAction(self.parallel_action)

        view_menu.addSeparator()

        # --- Display Style submenu ---
        style_menu = view_menu.addMenu("Display Style")
        style_group = QActionGroup(self)
        self.style_wire_action = QAction("Wireframe", self, checkable=True)
        self.style_wire_action.triggered.connect(
            lambda: self._set_render_style("wireframe"))
        style_group.addAction(self.style_wire_action)
        self.style_surf_action = QAction(
            "Surface with Edges", self, checkable=True, checked=True)
        self.style_surf_action.triggered.connect(
            lambda: self._set_render_style("surface"))
        style_group.addAction(self.style_surf_action)
        style_menu.addAction(self.style_wire_action)
        style_menu.addAction(self.style_surf_action)
        style_menu.addSeparator()
        # v4.0.1: default OFF (user preference; matches typical FEA preprocessors).
        self.shading_action = QAction(
            "Shading (Phong)", self, checkable=True, checked=False)
        self.shading_action.toggled.connect(self._toggle_shading)
        style_menu.addAction(self.shading_action)
        self.transparent_shells_action = QAction(
            "Transparent Shells", self, checkable=True)
        self.transparent_shells_action.toggled.connect(
            self._toggle_shell_transparency)
        style_menu.addAction(self.transparent_shells_action)
        self.ghost_mode_action = QAction(
            "Ghost Hidden Groups", self, checkable=True)
        self.ghost_mode_action.setChecked(self._ghost_mode_enabled)
        self.ghost_mode_action.toggled.connect(self._toggle_ghost_mode)
        style_menu.addAction(self.ghost_mode_action)

        # --- Color By submenu ---
        color_menu = view_menu.addMenu("Color By")
        color_group = QActionGroup(self)
        color_type = QAction("Element Type", self, checkable=True)
        color_type.triggered.connect(lambda: self._set_coloring_mode("type"))
        color_group.addAction(color_type)
        color_prop = QAction("Property ID", self, checkable=True, checked=True)
        color_prop.triggered.connect(
            lambda: self._set_coloring_mode("property"))
        color_group.addAction(color_prop)
        color_qual = QAction("Quality", self, checkable=True)
        color_qual.triggered.connect(
            lambda: self._set_coloring_mode("quality"))
        color_group.addAction(color_qual)
        color_results = QAction("Results", self, checkable=True)
        color_results.triggered.connect(
            lambda: self._set_coloring_mode("results"))
        color_group.addAction(color_results)
        for a in (color_type, color_prop, color_qual, color_results):
            color_menu.addAction(a)

        view_menu.addSeparator()

        # --- Annotations submenu (promoted Show Origin, Show Axes) ---
        annotations_menu = view_menu.addMenu("Annotations")
        self.show_origin_action = QAction("Show Origin", self, checkable=True)
        self.show_origin_action.toggled.connect(self._toggle_origin)
        annotations_menu.addAction(self.show_origin_action)
        # v5.0.0 item 2: axes triad toggle promoted from Display Settings.
        self.show_axes_action = QAction(
            "Show Axes Triad", self, checkable=True,
            checked=bool(self._display_settings.get('show_axes', True)))
        self.show_axes_action.toggled.connect(self._toggle_axes_triad)
        annotations_menu.addAction(self.show_axes_action)
        self.show_normals_action = QAction(
            "Show Element Normals", self, checkable=True)
        self.show_normals_action.toggled.connect(
            self._toggle_element_normals_visibility)
        annotations_menu.addAction(self.show_normals_action)
        self.show_free_edges_action = QAction(
            "Show Free Edges", self, checkable=True)
        self.show_free_edges_action.toggled.connect(
            self._toggle_free_edges_visibility)
        annotations_menu.addAction(self.show_free_edges_action)
        # v5.2.0 item 46: displacement vector arrows toggle. Replaces
        # the removed Vectors sub-tab in the Results sidebar.
        self.show_disp_vectors_action = QAction(
            "Show Displacement Vectors", self, checkable=True, checked=False)
        self.show_disp_vectors_action.toggled.connect(
            lambda on: self._on_vector_overlay_toggled('disp', bool(on)))
        annotations_menu.addAction(self.show_disp_vectors_action)
        # SPC reaction force arrows (orange) -- complementary to disp.
        self.show_reac_vectors_action = QAction(
            "Show Reaction Vectors", self, checkable=True, checked=False)
        self.show_reac_vectors_action.toggled.connect(
            lambda on: self._on_vector_overlay_toggled('reac', bool(on)))
        annotations_menu.addAction(self.show_reac_vectors_action)
        # In-plane principal stress arrows (red).
        self.show_principal_vectors_action = QAction(
            "Show Principal Stress Vectors", self, checkable=True, checked=False)
        self.show_principal_vectors_action.toggled.connect(
            lambda on: self._on_vector_overlay_toggled('prin', bool(on)))
        annotations_menu.addAction(self.show_principal_vectors_action)

        view_menu.addSeparator()

        # --- Clipping Plane submenu (v5.0.0 item 3) ---
        clipping_menu = view_menu.addMenu("Clipping Plane")
        self.clipping_toggle_action = QAction(
            "Toggle Active", self, checkable=True, checked=False)
        self.clipping_toggle_action.setShortcut("Ctrl+Shift+C")
        self.clipping_toggle_action.triggered.connect(
            self._toggle_clipping_plane)
        clipping_menu.addAction(self.clipping_toggle_action)
        clipping_define_action = QAction("Define…", self)
        clipping_define_action.setShortcut("Ctrl+Shift+X")
        clipping_define_action.triggered.connect(
            self._show_clipping_define_panel)
        clipping_menu.addAction(clipping_define_action)
        from PySide6.QtCore import QSettings as _QS
        _show_outline_default = _QS("NodeRunner", "NodeRunner").value(
            "cross_section/show_outline", False, type=bool)
        self.clipping_outline_action = QAction(
            "Show Plane Outline", self, checkable=True,
            checked=bool(_show_outline_default))
        self.clipping_outline_action.toggled.connect(
            self._toggle_clipping_outline)
        clipping_menu.addAction(self.clipping_outline_action)

        view_menu.addSeparator()

        # --- Tools submenu ---
        tools_submenu = view_menu.addMenu("Tools")
        capture_screenshot_action = QAction("Capture Screenshot...", self)
        capture_screenshot_action.triggered.connect(self._capture_screenshot)
        tools_submenu.addAction(capture_screenshot_action)
        self.presentation_action = QAction(
            "Presentation Mode", self, checkable=True)
        self.presentation_action.setShortcut("F11")
        self.presentation_action.toggled.connect(self._toggle_presentation_mode)
        tools_submenu.addAction(self.presentation_action)
        tools_submenu.addSeparator()
        edit_colors_action = QAction("Edit Colors...", self)
        edit_colors_action.triggered.connect(self._open_color_manager)
        tools_submenu.addAction(edit_colors_action)
        display_settings_action = QAction("Display Settings...", self)
        display_settings_action.triggered.connect(self._open_display_settings)
        tools_submenu.addAction(display_settings_action)

        settings_menu = menu_bar.addMenu("&Settings")
        export_defaults_action = QAction("Export Defaults...", self)
        export_defaults_action.triggered.connect(self._open_export_defaults_dialog)
        settings_menu.addAction(export_defaults_action)
        units_action = QAction("Units...", self)
        units_action.triggered.connect(self._open_units_dialog)
        settings_menu.addAction(units_action)
        settings_menu.addSeparator()
        # Phase 4.1: adaptive node decimation toggle.
        self.adaptive_lod_action = QAction(
            "Adaptive Node LOD (decimate huge clouds)", self, checkable=True,
        )
        self.adaptive_lod_action.setChecked(self._adaptive_lod_enabled)
        self.adaptive_lod_action.toggled.connect(self._toggle_adaptive_lod)
        settings_menu.addAction(self.adaptive_lod_action)
        # v3.2.0: element-level LOD threshold (professional).
        elem_lod_action = QAction("Element LOD Threshold...", self)
        elem_lod_action.triggered.connect(self._open_elem_lod_dialog)
        settings_menu.addAction(elem_lod_action)
        # Parallel BDF parsing toggle has been retired (was slower than
        # single-threaded due to GIL contention; left stale-True values
        # in some users' QSettings causing huge import hangs). The
        # action is kept as a hidden no-op so any external code or
        # registry references don't NoneType-crash.
        self.parallel_parsing_action = QAction(
            "Parallel BDF Parsing (retired)", self, checkable=True,
        )
        self.parallel_parsing_action.setVisible(False)
        self.parallel_parsing_action.toggled.connect(self._toggle_parallel_parsing)

        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("&About...", self)
        about_action.triggered.connect(self._show_about_dialog)
        help_menu.addAction(about_action)
        license_action = QAction("&License...", self)
        license_action.triggered.connect(self._show_license_dialog)
        help_menu.addAction(license_action)

        self.theme_button = QPushButton("Switch to Light Mode")
        self.theme_button.setFlat(True)
        menu_bar.setCornerWidget(self.theme_button, QtCore.Qt.TopRightCorner)




# --- NEW: Handler for the "Save BDF File..." menu action ---
    def _save_file_dialog(self):
        """v5.1.0 item 28: a single richer dialog covering path, solver
        target, field format, scope, and NR-META options. Replaces the
        two-step v5.0.x flow (QFileDialog -> ExportOptionsDialog)."""
        if not self.current_generator:
            self._update_status("No model to save.", is_error=True)
            return

        from node_runner.dialogs import SaveBdfDialog
        from PySide6.QtCore import QSettings
        settings = QSettings("NodeRunner", "NodeRunner")
        default_format = self._export_format or settings.value(
            "export/default_format", "short")
        default_target = settings.value("export/default_target", "generic")
        default_path = self._loaded_filepath or ""

        # Look up the active analysis set name (if any) for the Scope combo.
        active_set_name = None
        try:
            active_set = (self.analysis_sets or {}).get(
                self.active_analysis_set_id)
            if active_set is not None:
                active_set_name = getattr(active_set, 'name', None)
        except Exception:
            active_set_name = None

        group_names = sorted((self.groups or {}).keys())

        dlg = SaveBdfDialog(
            default_path=default_path,
            default_format=default_format,
            default_target=default_target,
            active_analysis_set_name=active_set_name,
            group_names=group_names,
            parent=self,
        )
        if not dlg.exec():
            self._update_status("Save cancelled.")
            return
        opts = dlg.result_payload
        filepath = opts.get('path')
        if not filepath:
            self._update_status("Save cancelled (no path).")
            return

        chosen_format = opts.get('field_format', 'short')
        chosen_target = opts.get('target', 'generic')
        if opts.get('remember'):
            settings.setValue("export/default_format", chosen_format)
            settings.setValue("export/default_target", chosen_target)

        # v5.1.0 item 28: build NR-META block if requested.
        meta_block = ""
        if opts.get('include_nr_meta'):
            try:
                from node_runner import __version__ as _nr_ver
                from node_runner.nr_meta import dump_meta
                meta_block = dump_meta(
                    self.groups,
                    hidden_groups=self._hidden_groups,
                    isolate_mode=getattr(self, '_isolate_mode', None),
                    nr_version=_nr_ver,
                )
            except Exception as exc:
                self._update_status(
                    f"NR-META build failed (continuing without it): {exc}",
                    is_error=True)
                meta_block = ""

        # Scope handling: 'analysis_set' and 'group' require Item 26's
        # scope helper. v5.1.0 ships full-model scoping first; the
        # other paths will route to scope.py once Item 26 lands.
        if opts.get('scope_kind') in ('analysis_set', 'group'):
            self._update_status(
                "Scoped export (AnalysisSet / Group) lands in a later "
                "v5.1.0 pass -- writing full model for now.",
                is_error=False)

        try:
            self._apply_case_control()
            self.current_generator._write_bdf(
                filepath,
                field_format=chosen_format,
                target=chosen_target,
                nr_meta=meta_block,
                sidecar_nrmeta=opts.get('sidecar_nrmeta', False),
            )
            self._export_format = chosen_format
            self._refresh_status_widgets()
            self._update_status(
                f"Model saved to {os.path.basename(filepath)} "
                f"({chosen_target} target, {chosen_format} field"
                f"{', NR-META embedded' if meta_block and not opts.get('sidecar_nrmeta') else ''}"
                f"{', NR-META sidecar' if meta_block and opts.get('sidecar_nrmeta') else ''}).")
            try:
                from node_runner.profiling import perf_event
                perf_event('save_bdf', 'completed',
                           target=chosen_target,
                           field_format=chosen_format,
                           nr_meta=bool(meta_block),
                           sidecar=bool(opts.get('sidecar_nrmeta')),
                           scope_kind=opts.get('scope_kind', 'full'))
            except Exception:
                pass
        except Exception as e:
            self._update_status(f"File save failed: {e}", is_error=True)
            QMessageBox.critical(self, "Error", f"Could not save file: {e}")

    # --- Phase 3: Gaussian-mapped node rendering ---
    # --- Phase 4.1: adaptive decimation ---
    _LOD_DECIMATE_THRESHOLD = 50_000   # cloud size that triggers stride
    _LOD_DISPLAY_TARGET = 10_000       # target display size after stride

    def _add_node_cloud(self, node_points, name='nodes_actor'):
        """Render a node cloud as small fixed-pixel-size markers (professional).

        Each node is drawn as a flat anti-aliased screen-space point of size
        ``self.node_size`` PIXELS - the size never grows or shrinks with the
        camera distance. This matches how professional FEA preprocessors
        (Patran, HyperMesh) render nodes: they stay as crisp markers
        you can read past, instead of expanding into a fuzzy "cloud" when
        you zoom out on a large model.

        We deliberately do NOT use PyVista's ``style='points_gaussian'`` (the
        Gaussian splat is rendered in WORLD units, so points grow into
        blobs when zoomed out) or ``render_points_as_spheres=True`` (3D
        billboards with the same world-space scaling problem).

        Above ``_LOD_DECIMATE_THRESHOLD`` points, we stride-sample down to
        ~``_LOD_DISPLAY_TARGET`` to keep the per-point GPU bookkeeping
        cheap on huge models. The full set of nodes is still owned by
        ``self.current_grid`` for picking; this only thins the standalone
        display cloud.

        Returns the vtkActor.
        """
        n_total = len(node_points)
        decimated = False
        if (
            self._adaptive_lod_enabled
            and n_total > self._LOD_DECIMATE_THRESHOLD
        ):
            stride = max(1, n_total // self._LOD_DISPLAY_TARGET)
            node_points = node_points[::stride]
            decimated = True

        poly = pv.PolyData(node_points)
        actor = self.plotter.add_points(
            poly,
            color=self.node_color,
            point_size=self.node_size,
            render_points_as_spheres=False,
            name=name,
            reset_camera=False,  # v4.0.13: was missing
        )

        if decimated:
            self._update_status(
                f"Adaptive LOD: showing ~{len(node_points):,} of "
                f"{n_total:,} nodes (toggle in Settings).",
            )

        return actor

    # --- Phase 2.3: render request (debounced) ---
    def _request_render(self):
        """Schedule a plotter render on the next 16ms tick, coalescing bursts.

        Multiple calls within the timer window collapse into a single
        render() call. Direct ``self.plotter.render()`` calls elsewhere
        in the codebase still execute immediately and bypass this debouncer
        (intentional: code that needs an immediate render keeps working).
        """
        if hasattr(self, "_render_timer") and self._render_timer is not None:
            self._render_timer.start()
        else:
            # Defensive: render synchronously if the timer was not yet built.
            try:
                self.plotter.render()
            except Exception:
                pass

    # --- Phase 2: flat action registry harvested from the menu bar ---
    @staticmethod
    def _clean_menu_label(text):
        """Strip '&' Qt mnemonics, keeping a literal '&&' as one '&'."""
        return text.replace("&&", "\x00").replace("&", "").replace("\x00", "&")

    def _build_action_registry(self):
        """Walk the menu bar and return a flat list of [(path, kw, QAction)].

        Only leaf actions (those that do nothing but trigger a slot) are
        added. Separators and sub-menu containers are skipped. The path is
        a `Parent > Child > Leaf` label with mnemonics stripped.
        """
        registry = []

        def visit(menu, prefix):
            for act in menu.actions():
                if act.isSeparator():
                    continue
                label = self._clean_menu_label(act.text()).strip()
                if not label:
                    continue
                full = f"{prefix} > {label}" if prefix else label
                sub = act.menu()
                if sub is not None:
                    visit(sub, full)
                else:
                    keywords = " ".join((label, prefix or "")).lower()
                    registry.append((full, keywords, act))

        bar = self.menuBar()
        for top_act in bar.actions():
            sub = top_act.menu()
            if sub is None:
                continue
            top_label = self._clean_menu_label(top_act.text()).strip()
            visit(sub, top_label)
        return registry

    # --- Phase 1: persistent status-bar widgets (model / nodes / elements / units / palette hint) ---
    def _build_status_widgets(self):
        """Add permanent labels to the status bar.

        Permanent widgets coexist with transient `showMessage()` calls; Qt
        renders transient messages on the left and permanent widgets on the
        right. v5.0.0: Export-format label moved into Export Defaults dialog;
        Ctrl+P command-palette hint added as a clickable label on the far
        right.
        """
        from PySide6.QtWidgets import QLabel, QFrame
        bar = self.statusBar()
        self._status_model_lbl = QLabel("Model: Untitled")
        self._status_nodes_lbl = QLabel("Nodes: 0")
        self._status_elements_lbl = QLabel("Elements: 0")
        self._status_units_lbl = QLabel(
            f"Units: {self._units}" if self._units else ""
        )
        # v5.0.0 items 16/17: results state segment, populated by
        # _on_mystran_results_loaded / _update_results_status. Empty
        # string when no results are loaded.
        self._status_results_lbl = QLabel("")
        try:
            from node_runner.solve import mystran_tips as _mtips
            self._status_results_lbl.setToolTip(_mtips.STATUS_RESULTS)
        except Exception:
            pass
        try:
            from node_runner.solve import mystran_tips as _mtips
            self._status_results_lbl.setToolTip(_mtips.STATUS_RESULTS)
        except Exception:
            pass
        for lbl in (self._status_model_lbl, self._status_nodes_lbl,
                    self._status_elements_lbl, self._status_units_lbl,
                    self._status_results_lbl):
            lbl.setStyleSheet("padding: 0 8px;")
        sep1 = QFrame(); sep1.setFrameShape(QFrame.VLine); sep1.setFrameShadow(QFrame.Sunken)
        sep_el = QFrame(); sep_el.setFrameShape(QFrame.VLine); sep_el.setFrameShadow(QFrame.Sunken)
        sep_units = QFrame(); sep_units.setFrameShape(QFrame.VLine); sep_units.setFrameShadow(QFrame.Sunken)
        sep_results = QFrame(); sep_results.setFrameShape(QFrame.VLine); sep_results.setFrameShadow(QFrame.Sunken)
        sep_hint = QFrame(); sep_hint.setFrameShape(QFrame.VLine); sep_hint.setFrameShadow(QFrame.Sunken)
        bar.addPermanentWidget(self._status_model_lbl)
        bar.addPermanentWidget(sep1)
        bar.addPermanentWidget(self._status_nodes_lbl)
        bar.addPermanentWidget(sep_el)
        bar.addPermanentWidget(self._status_elements_lbl)
        bar.addPermanentWidget(sep_units)
        bar.addPermanentWidget(self._status_units_lbl)
        bar.addPermanentWidget(sep_results)
        bar.addPermanentWidget(self._status_results_lbl)

        # v5.0.0 item 10: Ctrl+P command-palette hint, clickable.
        self._status_palette_hint = _ClickableLabel("Ctrl+P: command palette")
        self._status_palette_hint.setStyleSheet(
            "color: #6c7086; padding: 0 8px; font-size: 11px;"
        )
        self._status_palette_hint.setToolTip(
            "Click to open the command palette (or press Ctrl+P)."
        )
        self._status_palette_hint.clicked.connect(self._open_command_palette)
        bar.addPermanentWidget(sep_hint)
        bar.addPermanentWidget(self._status_palette_hint)

    def _refresh_status_widgets(self):
        """Recompute the four status labels from current state."""
        if not hasattr(self, '_status_model_lbl'):
            return  # status widgets not built yet (early-init path)
        if self._loaded_filepath:
            name = os.path.basename(self._loaded_filepath)
        else:
            name = "Untitled"
        self._status_model_lbl.setText(f"Model: {name}")
        node_count = 0
        element_count = 0
        if self.current_generator is not None and self.current_generator.model is not None:
            try:
                node_count = len(self.current_generator.model.nodes)
            except Exception:
                node_count = 0
            # v4.1.0: element count matches the tree's "Elements (N)" total
            # (elements + rigid_elements).
            try:
                m = self.current_generator.model
                element_count = len(m.elements) + len(m.rigid_elements)
            except Exception:
                element_count = 0
        self._status_nodes_lbl.setText(f"Nodes: {node_count:,}")
        self._status_elements_lbl.setText(f"Elements: {element_count:,}")
        self._status_units_lbl.setText(
            f"Units: {self._units}" if self._units else ""
        )

    def _on_model_loaded(self, filepath, detected_format):
        """Called after a successful BDF open/import.

        Stores the filepath and applies format detection per the PRD rule:
        the user's QSettings default (if set) wins on save, but for the
        initial export-format display we honor what the file actually was
        so round-tripping is intuitive.

        v5.1.0 item 28: also restores groups + tree state from any
        ``$ NR-META v1`` block embedded in the deck.
        """
        self._loaded_filepath = filepath
        if detected_format in ("short", "long", "free"):
            self._export_format = detected_format
        self._refresh_status_widgets()
        # v4.0.1: removed v4.0.0's cross-session hidden-groups restore.
        # `_hidden_groups` stays empty on a fresh open so the user
        # always sees the full model on load. In-session toggles are
        # still remembered via the plain `self._hidden_groups` set.

        # v5.1.0 item 28: restore Node Runner metadata if present.
        try:
            self._restore_nr_meta_from_model()
        except Exception:
            pass
        # Also try a sidecar .nrmeta file if it exists.
        try:
            self._restore_nr_meta_from_sidecar(filepath)
        except Exception:
            pass

    def _restore_nr_meta_from_model(self):
        """v5.1.0 item 28: parse any ``$ NR-META`` block the streaming
        reader stripped off the deck and restore groups + tree state.
        Idempotent; safe to call multiple times."""
        if not self.current_generator or not self.current_generator.model:
            return
        block = getattr(self.current_generator.model, '_nr_meta_block', '')
        if not block:
            return
        from node_runner.nr_meta import parse_meta
        meta = parse_meta(block)
        self._apply_parsed_nr_meta(meta, source='embedded')

    def _restore_nr_meta_from_sidecar(self, filepath):
        """v5.1.0 item 28: look for a ``<deck>.nrmeta`` sidecar file
        and merge its contents. The sidecar wins over embedded when
        both are present -- the assumption being that a user who
        chose sidecar specifically wanted the metadata externalized.
        """
        if not filepath:
            return
        import os as _os
        sidecar = str(filepath) + '.nrmeta'
        if not _os.path.exists(sidecar):
            return
        try:
            with open(sidecar, 'r', encoding='utf-8',
                      errors='replace') as fh:
                text = fh.read()
        except OSError:
            return
        from node_runner.nr_meta import parse_meta
        meta = parse_meta(text)
        self._apply_parsed_nr_meta(meta, source='sidecar')

    def _apply_parsed_nr_meta(self, meta, source='embedded'):
        """Common application path for ParsedMeta from either the
        embedded block or the sidecar file."""
        if not meta.has_content:
            return
        # Merge groups -- existing groups (built automatically by
        # _group_by_property etc.) are preserved unless the metadata
        # explicitly overrides them by name.
        for name, data in meta.groups.items():
            self.groups[name] = self._make_group_data(
                gid=data.get('gid', self._next_gid()),
                nodes=data.get('nodes'),
                elements=data.get('elements'),
                properties=data.get('properties'),
                materials=data.get('materials'),
                coords=data.get('coords'),
            )
        if meta.hidden_groups:
            self._hidden_groups = set(meta.hidden_groups)
        if meta.isolate_mode is not None:
            self._isolate_mode = meta.isolate_mode
        try:
            self._populate_groups_list()
        except Exception:
            pass
        try:
            from node_runner.profiling import perf_event
            perf_event('nr_meta', 'restored',
                       source=source,
                       n_groups=len(meta.groups),
                       n_hidden=len(meta.hidden_groups),
                       version=meta.version or '')
        except Exception:
            pass
        if meta.groups:
            self._update_status(
                f"Restored {len(meta.groups)} group(s) from "
                f"$ NR-META ({source}).")

    # --- Phase 1: Settings menu handlers ---
    def _open_export_defaults_dialog(self):
        """Settings -> Export Defaults..."""
        from node_runner.dialogs import ExportDefaultsDialog
        from PySide6.QtCore import QSettings
        settings = QSettings("NodeRunner", "NodeRunner")
        current = settings.value("export/default_format", self._export_format)
        if current not in ("short", "long", "free"):
            current = "short"
        dlg = ExportDefaultsDialog(current_format=current, parent=self)
        if dlg.exec():
            settings.setValue("export/default_format", dlg.selected_format)
            # If no file is loaded, also update the live status widget so the
            # user sees their choice take effect immediately.
            if self._loaded_filepath is None:
                self._export_format = dlg.selected_format
                self._refresh_status_widgets()
            self._update_status(f"Default export format set to {dlg.selected_format}.")

    def _open_units_dialog(self):
        """Settings -> Units..."""
        from node_runner.dialogs import UnitsDialog
        from PySide6.QtCore import QSettings
        dlg = UnitsDialog(current_units=self._units, parent=self)
        if dlg.exec():
            self._units = dlg.selected_units
            QSettings("NodeRunner", "NodeRunner").setValue("display/units", self._units)
            self._refresh_status_widgets()
            self._update_status(f"Units set to {self._units}.")

    def _open_command_palette(self):
        """Open the command palette imperatively. Same path as Ctrl+P."""
        from node_runner.dialogs import CommandPaletteDialog
        registry = getattr(self, '_action_registry', None) or []
        if not registry:
            return
        dlg = CommandPaletteDialog(registry, parent=self)
        dlg.exec()

    # --- Theme B: result browser + probe + expression contour ---
    def _build_result_browser_dock(self):
        """v5.2.0 item 47: no-op stub. Results UI now lives in the
        Post-Processing Toolbox (results_tab.py); tabular workflow
        moved to Tools -> Data Table. Method kept so initUI's call
        site continues to compile.
        """
        self._result_browser_dock = None
        self._anim_dock = None
        self._vector_dock = None
        self._anim_timeline_widget = None
        self._vector_overlay_widget = None
        try:
            from node_runner.profiling import perf_event
            perf_event('results', 'dock.eliminated',
                       reason='v5.2.0_toolbox')
        except Exception:
            pass

    def _show_result_panels(self, on: bool):
        """v5.2.0: no-op for back-compat -- the v3.x docks no longer exist."""
        return

    def _populate_result_browser(self):
        """Push current OP2 + model data into the result browser table model."""
        if not self._result_browser_dock or not self.op2_results:
            return
        sc_id = self._current_subcase_id_for_results()
        if sc_id is None:
            return
        sc = self.op2_results['subcases'].get(sc_id, {})
        # Nodes table: NID, |U|, Ux, Uy, Uz
        disps = sc.get('displacements', {}) or {}
        if disps:
            nids = sorted(disps.keys())
            ux = np.array([disps[n][0] for n in nids])
            uy = np.array([disps[n][1] for n in nids])
            uz = np.array([disps[n][2] for n in nids])
            mag = np.sqrt(ux * ux + uy * uy + uz * uz)
            self._result_browser_dock.update_nodal_results(
                ['NID', '|U|', 'Ux', 'Uy', 'Uz'],
                {'NID': np.array(nids, dtype=int), '|U|': mag,
                 'Ux': ux, 'Uy': uy, 'Uz': uz},
            )
        # Elements table: EID, vonMises, MaxPrin, MinPrin, Oxx, Oyy, Txy
        stresses = sc.get('stresses', {}) or {}
        if stresses:
            eids = sorted(stresses.keys())
            vm = np.array([stresses[e].get('von_mises', 0.0) for e in eids])
            maxp = np.array([stresses[e].get('max_principal', 0.0) for e in eids])
            minp = np.array([stresses[e].get('min_principal', 0.0) for e in eids])
            oxx = np.array([stresses[e].get('oxx', 0.0) for e in eids])
            oyy = np.array([stresses[e].get('oyy', 0.0) for e in eids])
            txy = np.array([stresses[e].get('txy', 0.0) for e in eids])
            self._result_browser_dock.update_element_results(
                ['EID', 'vonMises', 'MaxPrin', 'MinPrin', 'Oxx', 'Oyy', 'Txy'],
                {'EID': np.array(eids, dtype=int), 'vonMises': vm,
                 'MaxPrin': maxp, 'MinPrin': minp,
                 'Oxx': oxx, 'Oyy': oyy, 'Txy': txy},
            )

    def _current_subcase_id_for_results(self):
        """Best-effort current subcase ID for the result browser."""
        if not self.op2_results:
            return None
        # The existing UI tracks this via self.results_subcase_combo
        try:
            sc_id = self.results_subcase_combo.currentData()
        except Exception:
            sc_id = None
        if sc_id is None:
            keys = sorted(self.op2_results['subcases'].keys())
            sc_id = keys[0] if keys else None
        return sc_id

    def _on_result_entity_picked(self, entity_type: str, entity_id: int):
        """Double-click on a row -> highlight and zoom to the entity in 3D."""
        if not self.current_generator:
            return
        try:
            self._highlight_entities(entity_type, [entity_id])
        except Exception:
            pass

    def _on_expression_changed(self, expr: str):
        """Apply / clear a custom contour expression. Empty string -> revert."""
        self._contour_expression = (expr or "").strip()
        if self._contour_expression:
            self._update_status(f"Contour expression: {self._contour_expression}")
        else:
            self._update_status("Contour expression cleared.")
        # Trigger a re-render in results-color mode so the contour updates.
        if self.color_mode == "results":
            self._update_plot_visibility()

    def _evaluate_contour_expression(self, eids):
        """Evaluate the active expression against the current OP2 stresses
        for the given element IDs. Returns a numpy array aligned to eids,
        or None if no active expression / no OP2 loaded."""
        if not self._contour_expression or not self.op2_results:
            return None
        sc_id = self._current_subcase_id_for_results()
        if sc_id is None:
            return None
        sc = self.op2_results['subcases'].get(sc_id, {})
        stresses = sc.get('stresses', {}) or {}
        strains = sc.get('strains', {}) or {}
        n = len(eids)
        # Build named arrays. Missing values default to 0 so users get a
        # sensible result even if some elements weren't in the OP2.
        def col(d, key):
            return np.array([d.get(int(e), {}).get(key, 0.0) for e in eids])
        variables = {
            'stress_vm':         col(stresses, 'von_mises'),
            'stress_max_prin':   col(stresses, 'max_principal'),
            'stress_min_prin':   col(stresses, 'min_principal'),
            'oxx':               col(stresses, 'oxx'),
            'oyy':               col(stresses, 'oyy'),
            'txy':               col(stresses, 'txy'),
            'strain_vm':         col(strains, 'von_mises'),
            'strain_max_prin':   col(strains, 'max_principal'),
            'strain_min_prin':   col(strains, 'min_principal'),
            'exx':               col(strains, 'exx'),
            'eyy':               col(strains, 'eyy'),
            'gxy':               col(strains, 'gxy'),
            'eid':               np.array(eids, dtype=float),
        }
        try:
            from node_runner.expression import evaluate, ExpressionError
            result = evaluate(self._contour_expression, variables)
            if result.ndim == 0:
                result = np.full(n, float(result))
            return result.astype(float)
        except Exception as exc:
            self._update_status(f"Expression error: {exc}", is_error=True)
            return None

    # --- Probe tool (hover) ---

    def _toggle_probe_mode(self, on: bool):
        self._probe_enabled = bool(on)
        if not on and self._probe_tooltip is not None:
            self._probe_tooltip.hide()
        self._update_status(f"Probe mode {'enabled' if on else 'disabled'}.")

    def _ensure_probe_tooltip(self):
        if self._probe_tooltip is None:
            from PySide6.QtWidgets import QLabel
            tip = QLabel(self)
            tip.setWindowFlag(QtCore.Qt.ToolTip, True)
            tip.setStyleSheet(
                "QLabel { background: #1e1e2e; color: #cdd6f4; "
                "padding: 6px 8px; border: 1px solid #45475a; "
                "font-family: monospace; font-size: 11px; }"
            )
            tip.hide()
            self._probe_tooltip = tip
        return self._probe_tooltip

    def _show_probe_tip(self, screen_pos, lines: list[str]):
        tip = self._ensure_probe_tooltip()
        tip.setText("\n".join(lines))
        tip.adjustSize()
        tip.move(screen_pos.x() + 16, screen_pos.y() + 16)
        tip.show()

    def _on_probe_mouse_move(self):
        """Re-arm the probe throttle. The actual pick happens in _do_probe_pick
        after the user has stopped moving for ~80 ms - keeps things smooth."""
        if not self._probe_enabled:
            return
        self._probe_throttle.start()

    def _do_probe_pick(self):
        """Run a point + cell pick at the current cursor position; show a
        tooltip if anything was hit."""
        if not self._probe_enabled or not self.current_generator:
            return
        try:
            iren = self.plotter.interactor
            x, y = iren.GetEventPosition()
            renderer = self.plotter.renderer

            # Try element pick first (cells visible on top); fall back to node
            cell_picker = vtk.vtkCellPicker()
            cell_picker.SetTolerance(0.005)
            cell_picker.Pick(x, y, 0, renderer)
            cell_id = cell_picker.GetCellId()
            picked_actor = cell_picker.GetActor()
            if picked_actor and cell_id != -1:
                dataset = picked_actor.GetMapper().GetInput()
                if dataset and 'EID' in dataset.cell_data:
                    eid = int(dataset.cell_data['EID'][cell_id])
                    self._probe_element(eid)
                    return

            # Fallback to node pick
            point_picker = vtk.vtkPointPicker()
            point_picker.Pick(x, y, 0, renderer)
            point_id = point_picker.GetPointId()
            picked_actor = point_picker.GetActor()
            if picked_actor and point_id != -1:
                dataset = picked_actor.GetMapper().GetInput()
                if dataset and 'vtkOriginalPointIds' in dataset.point_data:
                    orig = int(dataset.point_data['vtkOriginalPointIds'][point_id])
                    if 0 <= orig < len(self.current_node_ids_sorted):
                        self._probe_node(self.current_node_ids_sorted[orig])
                        return

            # Nothing hit - hide tooltip if it was up
            if self._probe_tooltip is not None:
                self._probe_tooltip.hide()
        except Exception:
            # Probe is best-effort; never let a hover failure break interaction.
            pass

    def _probe_node(self, nid: int):
        if not self._probe_enabled or not self.current_generator:
            return
        if nid not in self.current_generator.model.nodes:
            return
        xyz = self.current_generator.model.nodes[nid].get_position()
        lines = [f"Node {nid}",
                 f"  xyz = ({xyz[0]:.4g}, {xyz[1]:.4g}, {xyz[2]:.4g})"]
        if self.op2_results:
            sc_id = self._current_subcase_id_for_results()
            sc = self.op2_results['subcases'].get(sc_id, {}) if sc_id is not None else {}
            disp = sc.get('displacements', {}).get(nid)
            if disp:
                mag = (disp[0] ** 2 + disp[1] ** 2 + disp[2] ** 2) ** 0.5
                lines.append(f"  |U| = {mag:.4g}")
                lines.append(f"  U   = ({disp[0]:.4g}, {disp[1]:.4g}, {disp[2]:.4g})")
        from PySide6.QtGui import QCursor
        self._show_probe_tip(QCursor.pos(), lines)

    def _probe_element(self, eid: int):
        if not self._probe_enabled or not self.current_generator:
            return
        elem = self.current_generator.model.elements.get(eid)
        if elem is None:
            return
        lines = [f"Element {eid} ({elem.type})", f"  PID = {elem.pid}"]
        if self.op2_results:
            sc_id = self._current_subcase_id_for_results()
            sc = self.op2_results['subcases'].get(sc_id, {}) if sc_id is not None else {}
            s = sc.get('stresses', {}).get(eid)
            if s:
                lines.append(f"  vonMises = {s.get('von_mises', 0.0):.4g}")
                lines.append(f"  MaxPrin  = {s.get('max_principal', 0.0):.4g}")
                lines.append(f"  MinPrin  = {s.get('min_principal', 0.0):.4g}")
        from PySide6.QtGui import QCursor
        self._show_probe_tip(QCursor.pos(), lines)

    # --- Vector overlay handlers ---

    def _on_vector_overlay_toggled(self, key: str, on: bool):
        if not self.current_generator or not self.op2_results:
            return
        actor_name = f"vec_overlay_{key}"
        self.plotter.remove_actor(actor_name, render=False)
        if not on:
            self._request_render()
            return
        sc_id = self._current_subcase_id_for_results()
        sc = self.op2_results['subcases'].get(sc_id, {}) if sc_id is not None else {}
        # v5.2.0 item 46: the v3.x VectorOverlayWidget is gone -- use a
        # sensible auto-scale (5% of bbox) instead of reading a widget.
        try:
            scale = float(self._vector_overlay_widget.scale_spin.value())  # type: ignore[union-attr]
        except Exception:
            scale = 0.05
        diag = max(self.current_grid.length, 1e-6) if self.current_grid is not None else 1.0
        length = scale * diag

        if key == 'disp':
            disps = sc.get('displacements', {}) or {}
            if not disps:
                return
            nids = sorted(disps.keys())
            pts, vecs = [], []
            for nid in nids:
                if nid in self.current_generator.model.nodes:
                    pts.append(self.current_generator.model.nodes[nid].get_position())
                    d = disps[nid]
                    vecs.append([d[0], d[1], d[2]])
            if not pts:
                return
            arr_pts = np.array(pts)
            arr_vecs = np.array(vecs)
            mags = np.linalg.norm(arr_vecs, axis=1)
            mmax = mags.max() if mags.size else 1.0
            if mmax < 1e-30:
                return
            arr_vecs = arr_vecs * (length / mmax)
            self.plotter.add_arrows(arr_pts, arr_vecs, color='#89b4fa', name=actor_name, reset_camera=False)
        elif key == 'reac':
            spcf = sc.get('spc_forces', {}) or {}
            if not spcf:
                return
            nids = sorted(spcf.keys())
            pts, vecs = [], []
            for nid in nids:
                if nid in self.current_generator.model.nodes:
                    pts.append(self.current_generator.model.nodes[nid].get_position())
                    f = spcf[nid]
                    vecs.append([f[0], f[1], f[2]])
            if not pts:
                return
            arr_pts = np.array(pts); arr_vecs = np.array(vecs)
            mags = np.linalg.norm(arr_vecs, axis=1)
            mmax = mags.max() if mags.size else 1.0
            if mmax < 1e-30:
                return
            arr_vecs = arr_vecs * (length / mmax)
            self.plotter.add_arrows(arr_pts, arr_vecs, color='#fab387', name=actor_name, reset_camera=False)
        elif key == 'prin':
            stresses = sc.get('stresses', {}) or {}
            if not stresses or self.current_grid is None:
                return
            cell_data = self.current_grid.cell_data.get('EID') if hasattr(self.current_grid, 'cell_data') else None
            if cell_data is None:
                return
            centers = self.current_grid.cell_centers().points

            # Compute the in-plane principal direction at each element from
            # the (oxx, oyy, txy) stress tensor and re-express it in world
            # coordinates using the element's local plane normal. The shell
            # element's local x is taken as (n2 - n1); local z is the face
            # normal; local y closes the right-hand triad.
            pts, vecs = [], []
            for i, eid in enumerate(cell_data):
                s = stresses.get(int(eid))
                if not s:
                    continue
                elem = self.current_generator.model.elements.get(int(eid))
                if elem is None or elem.type not in ('CQUAD4', 'CMEMBRAN', 'CTRIA3'):
                    continue
                try:
                    nlist = elem.nodes
                    p0 = np.asarray(self.current_generator.model.nodes[nlist[0]].get_position())
                    p1 = np.asarray(self.current_generator.model.nodes[nlist[1]].get_position())
                    p2 = np.asarray(self.current_generator.model.nodes[nlist[2]].get_position())
                except (KeyError, IndexError):
                    continue
                ex = p1 - p0
                ex_n = np.linalg.norm(ex)
                if ex_n < 1e-12:
                    continue
                ex = ex / ex_n
                ez = np.cross(ex, p2 - p0)
                ez_n = np.linalg.norm(ez)
                if ez_n < 1e-12:
                    continue
                ez = ez / ez_n
                ey = np.cross(ez, ex)
                # In-plane stress tensor (oxx, oyy, txy) -> principal angle
                # measured from ex toward ey: theta = 0.5 * atan2(2*txy, oxx-oyy)
                oxx = float(s.get('oxx', 0.0))
                oyy = float(s.get('oyy', 0.0))
                txy = float(s.get('txy', 0.0))
                denom = oxx - oyy
                theta = 0.5 * math.atan2(2.0 * txy, denom) if abs(denom) + abs(txy) > 1e-30 else 0.0
                d_local = math.cos(theta) * ex + math.sin(theta) * ey
                mag = float(s.get('max_principal', 0.0))
                pts.append(centers[i])
                vecs.append([d_local[0] * mag, d_local[1] * mag, d_local[2] * mag])
            if not pts:
                return
            arr_pts = np.array(pts); arr_vecs = np.array(vecs)
            mags = np.linalg.norm(arr_vecs, axis=1)
            mmax = mags.max() if mags.size else 1.0
            if mmax < 1e-30:
                return
            arr_vecs = arr_vecs * (length / mmax)
            self.plotter.add_arrows(arr_pts, arr_vecs, color='#f38ba8', name=actor_name, reset_camera=False)
        self._request_render()

    def _on_vector_overlay_scale(self, _new_scale):
        # v5.2.0 item 46: the v3.x VectorOverlayWidget is gone; this
        # hook stays for back-compat but no longer reads widget state.
        # The View menu toggle drives the arrows directly.
        return

    def _on_timeline_phase_changed(self, phase: float):
        """User scrubbed/played the animation timeline. Update mode-shape phase."""
        # v5.1.2 item 39: any animation interaction should engage
        # Results mode so the user sees motion instead of a static
        # mesh. Cheap no-op if already in results mode.
        try:
            self._ensure_results_mode_active(trigger='animate')
        except Exception:
            pass
        self._anim_phase = float(phase) * 2 * math.pi
        # Reuse existing animation tick logic
        self._on_animation_tick()

    # --- v5.0.0 item 3: Clipping plane (floating Define panel) ---
    def _ensure_clipping_panel(self):
        """Lazy-build the floating Define-Plane panel."""
        if self._clipping_panel is None:
            from node_runner.dialogs.cross_section import CrossSectionPanel
            self._clipping_panel = CrossSectionPanel(
                self._cross_section, self, parent=self)
            self._clipping_panel.plane_committed.connect(
                self._sync_clipping_toggle_action)
        return self._clipping_panel

    def _show_clipping_define_panel(self):
        """View -> Clipping Plane -> Define...: floating panel."""
        panel = self._ensure_clipping_panel()
        panel.snapshot()
        try:
            from node_runner.profiling import perf_event
            perf_event('clipping', 'define_open',
                       had_definition=self._cross_section.definition is not None,
                       was_enabled=self._cross_section.is_enabled)
        except Exception:
            pass
        # Position at bottom-center of the main window so the user knows
        # where it'll show up (mirrors EntitySelectionBar's positioning).
        panel.show()
        panel.adjustSize()
        try:
            pg = self.frameGeometry()
            panel.move(
                pg.x() + max(0, (pg.width() - panel.sizeHint().width()) // 2),
                pg.y() + pg.height() - panel.sizeHint().height() - 80,
            )
        except Exception:
            pass
        panel.raise_()
        panel.activateWindow()

    def _toggle_clipping_plane(self, checked):
        """View -> Clipping Plane -> Toggle Active."""
        controller = self._cross_section
        if controller.definition is None:
            # No definition yet -> open the panel so the user can define
            # one first. Restore the menu checkbox to OFF in the
            # meantime; it flips on after the user hits OK.
            self.clipping_toggle_action.blockSignals(True)
            self.clipping_toggle_action.setChecked(False)
            self.clipping_toggle_action.blockSignals(False)
            self._show_clipping_define_panel()
            return
        try:
            if checked:
                controller.enable()
            else:
                controller.disable()
            try:
                from node_runner.profiling import perf_event
                perf_event('clipping', 'toggle',
                           on=bool(checked),
                           had_definition=True)
            except Exception:
                pass
        except Exception as exc:
            QMessageBox.warning(self, "Clipping Plane", str(exc))

    def _toggle_clipping_outline(self, checked):
        """View -> Clipping Plane -> Show Plane Outline."""
        from PySide6.QtCore import QSettings
        QSettings("NodeRunner", "NodeRunner").setValue(
            "cross_section/show_outline", bool(checked))
        try:
            self._cross_section.set_outline_visible(bool(checked))
        except Exception:
            pass
        try:
            from node_runner.profiling import perf_event
            perf_event('clipping', 'outline_toggle', on=bool(checked))
        except Exception:
            pass

    def _sync_clipping_toggle_action(self):
        """Refresh the menu's Toggle Active checkbox after a Define
        panel OK so its state reflects whether clipping is active.
        """
        try:
            self.clipping_toggle_action.blockSignals(True)
            self.clipping_toggle_action.setChecked(
                bool(self._cross_section.is_enabled))
        finally:
            self.clipping_toggle_action.blockSignals(False)

    def _toggle_axes_triad(self, checked):
        """View -> Annotations -> Show Axes Triad.

        v5.0.0 item 2: this toggle was promoted from the Display Settings
        dialog. We update the live actor immediately AND persist the
        choice via QSettings so the next viewer rebuild honours it.
        """
        from PySide6.QtCore import QSettings
        self._display_settings['show_axes'] = bool(checked)
        QSettings("NodeRunner", "NodeRunner").setValue(
            "display/show_axes", bool(checked))
        if hasattr(self, 'axes_actor') and self.axes_actor is not None:
            try:
                self.axes_actor.SetVisibility(bool(checked))
                self.plotter.render()
            except Exception:
                pass

    # --- legacy: kept so old shortcuts / commands still work if any
    # external caller references _open_cross_section_dialog by name. ---
    def _open_cross_section_dialog(self):
        self._show_clipping_define_panel()

    # --- Phase 4 toggles ---
    def _open_convert_units_dialog(self):
        """Tools > Convert Units... - run a professional multi-factor scaling."""
        if not self.current_generator:
            QMessageBox.warning(self, "No Model",
                                "Open or create a model before converting units.")
            return
        from node_runner.dialogs import UnitConversionDialog
        from node_runner.commands import ConvertUnitsCommand
        dlg = UnitConversionDialog(self)
        if not dlg.exec():
            return
        f = dlg.factors
        if any(v <= 0 for v in f.values()):
            QMessageBox.warning(self, "Invalid factors",
                                "All factors must be positive non-zero numbers.")
            return
        cmd = ConvertUnitsCommand(
            length=f['length'], force=f['force'], mass=f['mass'],
        )
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Unit conversion applied: L*{f['length']:g}, "
            f"F*{f['force']:g}, M*{f['mass']:g}. Use Edit > Undo to revert."
        )
        # If the user opted in, also refresh the unit hint label after
        # conversion so the status bar reflects the new system.
        if dlg.update_label_hint:
            self._open_units_dialog()
        # Force a viewer rebuild so the rescaled coordinates take effect.
        self._update_viewer(self.current_generator, reset_camera=True)

    def _toggle_ghost_mode(self, on):
        """View -> Ghost Hidden Groups."""
        from PySide6.QtCore import QSettings
        self._ghost_mode_enabled = bool(on)
        QSettings("NodeRunner", "NodeRunner").setValue(
            "render/ghost_mode", self._ghost_mode_enabled,
        )
        # Force a visibility refresh so the ghost actor appears or vanishes.
        if self.current_generator and self.current_grid:
            self._update_plot_visibility()
        self._update_status(
            f"Ghost mode {'enabled' if on else 'disabled'}.",
        )

    def _toggle_adaptive_lod(self, on):
        """Settings -> Adaptive Node LOD."""
        from PySide6.QtCore import QSettings
        self._adaptive_lod_enabled = bool(on)
        QSettings("NodeRunner", "NodeRunner").setValue(
            "render/adaptive_lod", self._adaptive_lod_enabled,
        )
        if self.current_generator and self.current_grid:
            self._update_plot_visibility()
        self._update_status(
            f"Adaptive LOD {'enabled' if on else 'disabled'}.",
        )

    def _toggle_parallel_parsing(self, on):
        """Settings -> Parallel BDF Parsing (experimental)."""
        from PySide6.QtCore import QSettings
        self._parallel_parsing_enabled = bool(on)
        QSettings("NodeRunner", "NodeRunner").setValue(
            "import/parallel_parsing", self._parallel_parsing_enabled,
        )
        self._update_status(
            f"Parallel parsing {'enabled' if on else 'disabled'}.",
        )

    def _open_elem_lod_dialog(self):
        """Settings -> Element LOD Threshold...

        Above this many cells the displayed mesh stride-samples shells
        and extracts the outer surface of solids. Picking and queries
        still see every cell (current_grid is unchanged).
        """
        from PySide6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox,
            QDialogButtonBox,
        )
        from PySide6.QtCore import QSettings
        dlg = QDialog(self)
        dlg.setWindowTitle("Element LOD Threshold")
        dlg.setMinimumWidth(420)
        v = QVBoxLayout(dlg)
        v.addWidget(QLabel(
            "Above this cell count the displayed mesh is decimated for "
            "speed. Picking and queries still see every cell.\n\n"
            "Default: 500,000. Range: 50,000 - 5,000,000."
        ))
        row = QHBoxLayout()
        row.addWidget(QLabel("Threshold:"))
        spin = QSpinBox()
        spin.setRange(50_000, 5_000_000)
        spin.setSingleStep(50_000)
        spin.setValue(int(self._elem_lod_threshold))
        spin.setSuffix(" cells")
        row.addWidget(spin)
        row.addStretch(1)
        v.addLayout(row)
        if getattr(self, '_lod_info', None) and self._lod_info.get('lod_active'):
            info = QLabel(
                f"Currently active: showing "
                f"{self._lod_info['displayed_cells']:,} of "
                f"{self._lod_info['total_cells']:,} cells."
            )
            info.setStyleSheet("color: #888;")
            v.addWidget(info)
        bb = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        v.addWidget(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        new_threshold = int(spin.value())
        if new_threshold == self._elem_lod_threshold:
            return
        self._elem_lod_threshold = new_threshold
        QSettings("NodeRunner", "NodeRunner").setValue(
            "render/elem_lod_threshold", new_threshold,
        )
        self._update_status(
            f"Element LOD threshold set to {new_threshold:,} cells. "
            f"Re-open or reload the model to apply.",
        )
        # If there's a current model, rebuild the display grid in place.
        if self.current_grid is not None and self.current_grid.n_cells > 0:
            try:
                from node_runner.scene_build import apply_display_lod
                self.current_display_grid, self._lod_info = apply_display_lod(
                    self.current_grid, threshold=new_threshold,
                )
                self._update_plot_visibility()
            except Exception:
                pass

    # --- NEW: Handler for the "License..." menu action ---
    def _show_license_dialog(self):
        license_text = """
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent

      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS
        """
        dialog = QDialog(self)
        dialog.setWindowTitle("License Information")
        dialog.setMinimumSize(600, 500)
        layout = QVBoxLayout(dialog)
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setText(license_text)
        text_edit.setFontFamily("monospace")
        layout.addWidget(text_edit)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(dialog.accept)
        layout.addWidget(button_box)
        dialog.exec()

    # --- BUG FIX: Implement the missing method to add material cards ---
    def _add_material_card_from_params(self, params):
        """Adds a MAT1, MAT8, or MAT9 card to the model from a parameters dict."""
        params_copy = params.copy()
        mat_type = params_copy.pop('type')
        mid = params_copy.pop('mid')
        comment = params_copy.pop('comment')
        model = self.current_generator.model

        try:
            if mat_type == 'MAT1':
                model.add_mat1(mid, params_copy['E'], params_copy['G'], params_copy['nu'],
                               rho=params_copy.get('rho', 0.0), a=params_copy.get('a', 0.0),
                               tref=params_copy.get('tref', 0.0), ge=params_copy.get('ge', 0.0), comment=comment)
            elif mat_type == 'MAT8':
                model.add_mat8(mid, params_copy['E1'], params_copy['E2'], params_copy['nu12'],
                               g12=params_copy.get('G12', 0.0), g1z=params_copy.get('G1z', 0.0), g2z=params_copy.get('G2z', 0.0),
                               rho=params_copy.get('rho', 0.0), a1=params_copy.get('a1', 0.0), a2=params_copy.get('a2', 0.0),
                               tref=params_copy.get('tref', 0.0), comment=comment)
            elif mat_type == 'MAT9':
                g_values = {key.lower(): val for key, val in params_copy.items() if key.startswith('G')}
                model.add_mat9(mid,
                               rho=params_copy.get('rho', 0.0), 
                               tref=params_copy.get('tref', 0.0),
                               comment=comment, **g_values)
        except KeyError as e:
            QMessageBox.warning(self, "Creation Error", f"Could not create material. Missing parameter: {e}")

    def _edit_material(self, mid, new_params):
        """Edit a material via undo-able command. Called by MaterialEditorDialog."""
        cmd = EditMaterialCommand(mid, new_params)
        self.command_manager.execute(cmd, self.current_generator.model)

    def _edit_property(self, pid, new_params):
        """Edit a property via undo-able command. Called by PropertyEditorDialog."""
        cmd = EditPropertyCommand(pid, new_params)
        self.command_manager.execute(cmd, self.current_generator.model)

    def _open_coincident_node_checker(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            self._update_status("No model loaded to check.", is_error=True)
            QMessageBox.warning(self, "No Model", "Please load a model before checking for coincident nodes.")
            return

        dialog = CoincidentNodeDialog(main_window=self, parent=self)
        
        # The dialog's exec() method is blocking. Code will only continue here after the
        # user has clicked "Merge Selected" (which calls accept()) or "Close".
        if dialog.exec() == QDialog.Accepted:
            groups_to_merge = dialog.get_groups_to_merge()

            if not groups_to_merge:
                self._update_status("Merge cancelled: No groups were selected.", is_error=True)
                return

            # Wrap the merge in a command for undo/redo
            cmd = MergeNodesCommand(groups_to_merge)
            self.command_manager.execute(cmd, self.current_generator.model)

            merged_count = len(cmd._deleted_nodes)
            group_count = len(groups_to_merge)
            QMessageBox.information(self, "Merge Complete",
                                    f"Successfully merged {merged_count} node(s) from {group_count} group(s).")
            
            self._update_status(f"Merged {merged_count} coincident nodes.")
            
            # Refresh the entire viewer and model tree to reflect the changes
            self._update_viewer(self.current_generator, reset_camera=False)

    def _open_duplicate_element_checker(self):
        if not self.current_generator or not self.current_generator.model.elements:
            self._update_status("No model loaded to check.", is_error=True)
            QMessageBox.warning(self, "No Model", "Please load a model before checking for duplicate elements.")
            return

        dialog = DuplicateElementDialog(main_window=self, parent=self)

        if dialog.exec() == QDialog.Accepted:
            eids_to_delete = dialog.get_eids_to_delete()

            if not eids_to_delete:
                self._update_status("Delete cancelled: No duplicates were selected.", is_error=True)
                return

            cmd = DeleteElementsCommand(eids_to_delete)
            self.command_manager.execute(cmd, self.current_generator.model)

            QMessageBox.information(self, "Delete Complete",
                                    f"Successfully deleted {len(eids_to_delete)} duplicate element(s).")

            self._update_status(f"Deleted {len(eids_to_delete)} duplicate elements.")
            self._update_viewer(self.current_generator, reset_camera=False)

    # ─── Measurement Tools ───────────────────────────────────────────────

    def _activate_multi_node_picker(self, count, callback):
        """Activates a sequential multi-node picker.

        Args:
            count: Number of nodes to pick.
            callback: Function to call with list of picked node IDs.
        """
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "Please load a model first.")
            return
        self._measurement_picks = []
        self._measurement_pick_count = count
        self._measurement_callback = callback
        self.picking_target_callback = self._on_measurement_node_picked
        self._set_picking_mode("Node", True)
        self._update_status(f"Pick node 1 of {count}...")

    def _on_measurement_node_picked(self, nid_str):
        """Called by the interactor when a node is picked during measurement."""
        nid = int(nid_str)
        self._measurement_picks.append(nid)
        picked_so_far = len(self._measurement_picks)

        # Highlight all picked nodes so far
        self._highlight_entities('Node', self._measurement_picks)

        if picked_so_far < self._measurement_pick_count:
            # Need more picks - re-arm on next event loop tick
            # (interactor.py lines 89-90 auto-disable picking after each pick)
            remaining = self._measurement_pick_count - picked_so_far
            self._update_status(f"Pick node {picked_so_far + 1} of {self._measurement_pick_count}...")
            QtCore.QTimer.singleShot(0, self._rearm_measurement_picker)
        else:
            # All picks collected - invoke the callback
            callback = self._measurement_callback
            self._measurement_callback = None
            self._measurement_pick_count = 0
            callback(list(self._measurement_picks))
            self._measurement_picks = []

    def _rearm_measurement_picker(self):
        """Re-enters picking mode for the next measurement node pick."""
        self.picking_target_callback = self._on_measurement_node_picked
        self._set_picking_mode("Node", True)

    def _measure_distance(self):
        """Activates 2-node distance measurement."""
        self._activate_multi_node_picker(2, self._compute_distance)

    def _compute_distance(self, nids):
        """Computes and displays distance between two picked nodes."""
        model = self.current_generator.model
        p1 = np.array(model.nodes[nids[0]].get_position())
        p2 = np.array(model.nodes[nids[1]].get_position())
        dx, dy, dz = p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]
        dist = np.linalg.norm(p2 - p1)

        # Draw a yellow line between the two points
        line = pv.Line(p1, p2)
        self.plotter.add_mesh(line, color='yellow', line_width=3,
                              name='measurement_line', reset_camera=False)

        sep = "-" * 30
        text = "\n".join([
            sep,
            f"  Node A:    {nids[0]}",
            f"  Node B:    {nids[1]}",
            sep,
            f"  Distance:  {dist:.6f}",
            sep,
            f"  dX:        {dx:+.6f}",
            f"  dY:        {dy:+.6f}",
            f"  dZ:        {dz:+.6f}",
            sep,
        ])
        self._update_status(f"Distance: {dist:.6f}  (Node {nids[0]} to Node {nids[1]})")
        InfoDialog("Measure Distance", text, self).exec()

    def _measure_angle(self):
        """Activates 3-node angle measurement (angle at middle node)."""
        self._activate_multi_node_picker(3, self._compute_angle)

    def _compute_angle(self, nids):
        """Computes and displays the angle at the middle node (node 2)."""
        model = self.current_generator.model
        p1 = np.array(model.nodes[nids[0]].get_position())
        p2 = np.array(model.nodes[nids[1]].get_position())  # vertex
        p3 = np.array(model.nodes[nids[2]].get_position())

        v1 = p1 - p2
        v2 = p3 - p2
        len1 = np.linalg.norm(v1)
        len2 = np.linalg.norm(v2)

        if len1 < 1e-12 or len2 < 1e-12:
            self._update_status("Cannot compute angle: two nodes are coincident.", is_error=True)
            return

        cos_angle = np.clip(np.dot(v1, v2) / (len1 * len2), -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)

        # Draw lines from vertex to endpoints
        line1 = pv.Line(p2, p1)
        line2 = pv.Line(p2, p3)
        merged = line1 + line2
        self.plotter.add_mesh(merged, color='yellow', line_width=3,
                              name='measurement_line', reset_camera=False)

        sep = "-" * 30
        text = "\n".join([
            sep,
            f"  Node A:     {nids[0]}",
            f"  Vertex B:   {nids[1]}",
            f"  Node C:     {nids[2]}",
            sep,
            f"  Angle:      {angle_deg:.4f} deg",
            sep,
        ])
        self._update_status(f"Angle at Node {nids[1]}: {angle_deg:.4f} deg")
        InfoDialog("Measure Angle", text, self).exec()


     # --- NEW: Handler to open the Element Normals Flip tool (Phase 4) ---

    def _open_element_flip_tool(self):
        """Launches the entity selection dialog to choose shell elements to flip."""
        if self.selection_bar.isVisible():
            return

        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model contains no elements to flip.")
            return

        all_shell_eids = [
            eid for eid, elem in self.current_generator.model.elements.items()
            if elem.type in ['CQUAD4', 'CTRIA3', 'CMEMBRAN']
        ]
        if not all_shell_eids:
            QMessageBox.information(self, "No Shells", "No shell elements found in the model.")
            return

        self._start_selection('Element', all_shell_eids,
                              self._on_flip_elements_accept,
                              self._build_select_by_data())




    def _on_flip_elements_accept(self):
        selected_eids = self.selection_bar.get_selected_ids()
        normals_were_visible = self.show_normals_check.isChecked()
        self._end_selection_mode()

        # --- The filtering logic that was here is now removed ---

        if not selected_eids:
            self._update_status("Flip operation cancelled: No elements selected.")
            return

        # Wrap the flip in a command for undo/redo
        cmd = FlipNormalsCommand(selected_eids)
        self.command_manager.execute(cmd, self.current_generator.model)
        # Count how many were actually flipped (filter to shells)
        flipped_count = sum(1 for eid in selected_eids
                           if eid in self.current_generator.model.elements
                           and self.current_generator.model.elements[eid].type in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'))

        # The confirmation message will now be correct because the flipped_count is correct.
        if flipped_count > 0:
            QMessageBox.information(self, "Operation Complete", f"Successfully flipped normals for {flipped_count} of {len(selected_eids)} selected element(s).")
            self._update_status(f"Flipped normals for {flipped_count} elements.")
        else:
            QMessageBox.warning(self, "No Valid Elements", f"You selected {len(selected_eids)} element(s), but none were flippable plates.")
            self._update_status("Flip operation cancelled: No valid shell elements were selected.")


        # Viewer update logic remains the same
        self._update_viewer(self.current_generator, reset_camera=False)
        if normals_were_visible:
            self._toggle_element_normals_visibility()

    # --- Theme A: mesh-edit handlers (split / refine / combine / smooth /
    #     mirror / copy). Each follows the existing flip-normals pattern:
    #     gather a selection via _start_selection, then run the matching
    #     Command via the command_manager so undo/redo works. ---

    def _eligible_eids_for_quads(self):
        if not self.current_generator:
            return []
        return [
            eid for eid, e in self.current_generator.model.elements.items()
            if e.type in ('CQUAD4', 'CMEMBRAN')
        ]

    def _eligible_eids_for_trias(self):
        if not self.current_generator:
            return []
        return [
            eid for eid, e in self.current_generator.model.elements.items()
            if e.type == 'CTRIA3'
        ]

    def _eligible_eids_for_shells(self):
        if not self.current_generator:
            return []
        return [
            eid for eid, e in self.current_generator.model.elements.items()
            if e.type in ('CQUAD4', 'CMEMBRAN', 'CTRIA3')
        ]

    def _eligible_eids_all(self):
        if not self.current_generator:
            return []
        return list(self.current_generator.model.elements.keys())

    def _open_split_quad_tool(self):
        if self.selection_bar.isVisible():
            return
        eids = self._eligible_eids_for_quads()
        if not eids:
            QMessageBox.warning(self, "No Quads", "No CQUAD4 elements in the model.")
            return
        self._start_selection(
            'Element', eids, self._on_split_quad_accept,
            self._build_select_by_data(),
        )

    def _on_split_quad_accept(self):
        from node_runner.commands import SplitQuadCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Split cancelled: no elements selected.")
            return
        cmd = SplitQuadCommand(selected)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Split {len(selected)} quad(s) into trias.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_refine_tool(self):
        if self.selection_bar.isVisible():
            return
        eids = self._eligible_eids_for_shells()
        if not eids:
            QMessageBox.warning(self, "No Shells", "No shell elements in the model.")
            return
        self._start_selection(
            'Element', eids, self._on_refine_accept,
            self._build_select_by_data(),
        )

    def _on_refine_accept(self):
        from node_runner.commands import RefineElementsCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Refine cancelled: no elements selected.")
            return
        cmd = RefineElementsCommand(selected)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Refined {len(selected)} element(s) (1 -> 4 each).")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_combine_tool(self):
        if self.selection_bar.isVisible():
            return
        eids = self._eligible_eids_for_trias()
        if not eids:
            QMessageBox.warning(self, "No Trias", "No CTRIA3 elements in the model.")
            return
        from node_runner.dialogs import CombineTriasDialog
        dlg = CombineTriasDialog(self)
        if not dlg.exec():
            return
        self._combine_angle_tol = dlg.angle_tol_deg
        self._start_selection(
            'Element', eids, self._on_combine_accept,
            self._build_select_by_data(),
        )

    def _on_combine_accept(self):
        from node_runner.commands import CombineTriasCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Combine cancelled: no elements selected.")
            return
        cmd = CombineTriasCommand(selected, angle_tol_deg=self._combine_angle_tol)
        self.command_manager.execute(cmd, self.current_generator.model)
        n_quads = len(cmd._new_quads)
        self._update_status(
            f"Combined {len(selected)} triangle(s) into {n_quads} quad(s)."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_smooth_tool(self):
        if self.selection_bar.isVisible():
            return
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Nodes", "The model has no nodes to smooth.")
            return
        from node_runner.dialogs import SmoothNodesDialog
        dlg = SmoothNodesDialog(self)
        if not dlg.exec():
            return
        self._smooth_iter = dlg.iterations
        self._smooth_factor = dlg.factor
        self._smooth_pin = dlg.pin_free_edges
        all_nids = list(self.current_generator.model.nodes.keys())
        self._start_selection(
            'Node', all_nids, self._on_smooth_accept,
            self._build_select_by_data(),
        )

    def _on_smooth_accept(self):
        from node_runner.commands import SmoothNodesCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Smooth cancelled: no nodes selected.")
            return
        cmd = SmoothNodesCommand(
            selected,
            iterations=self._smooth_iter,
            factor=self._smooth_factor,
            pin_free_edges=self._smooth_pin,
        )
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Smoothed {len(selected)} node(s) ({self._smooth_iter} iterations).",
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_mirror_elements_tool(self):
        if self.selection_bar.isVisible():
            return
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model has no elements to mirror.")
            return
        from node_runner.dialogs import MirrorElementsDialog
        dlg = MirrorElementsDialog(parent=self)
        if not dlg.exec():
            return
        self._mirror_plane = dlg.plane
        self._mirror_value = dlg.value
        self._mirror_weld_tol = dlg.weld_tol
        self._start_selection(
            'Element', self._eligible_eids_all(),
            self._on_mirror_accept, self._build_select_by_data(),
        )

    def _on_mirror_accept(self):
        from node_runner.commands import MirrorElementsCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Mirror cancelled: no elements selected.")
            return
        cmd = MirrorElementsCommand(
            selected,
            plane=self._mirror_plane,
            plane_value=self._mirror_value,
            weld_tol=self._mirror_weld_tol,
        )
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Mirrored {len(cmd._created_eids)} element(s) across "
            f"{self._mirror_plane}={self._mirror_value}."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_copy_elements_tool(self):
        if self.selection_bar.isVisible():
            return
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model has no elements to copy.")
            return
        from node_runner.dialogs import CopyElementsDialog
        dlg = CopyElementsDialog(self)
        if not dlg.exec():
            return
        self._copy_translate = dlg.translation
        self._copy_axis = dlg.rotation_axis
        self._copy_angle = dlg.rotation_angle_deg
        self._copy_center = dlg.rotation_center
        self._start_selection(
            'Element', self._eligible_eids_all(),
            self._on_copy_accept, self._build_select_by_data(),
        )

    def _on_copy_accept(self):
        from node_runner.commands import CopyElementsCommand
        selected = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not selected:
            self._update_status("Copy cancelled: no elements selected.")
            return
        cmd = CopyElementsCommand(
            selected,
            translate=self._copy_translate,
            axis=self._copy_axis,
            angle_deg=self._copy_angle,
            center=self._copy_center,
        )
        self.command_manager.execute(cmd, self.current_generator.model)
        rot_note = (f", rotate {self._copy_angle:.2f} deg about {self._copy_axis}"
                    if abs(self._copy_angle) > 1e-9 else "")
        self._update_status(
            f"Copied {len(cmd._created_eids)} element(s) by "
            f"{self._copy_translate}{rot_note}."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _open_insert_edge_node_tool(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Nodes", "The model has no nodes.")
            return
        from node_runner.dialogs import InsertEdgeNodeDialog
        from node_runner.commands import InsertEdgeNodeCommand
        dlg = InsertEdgeNodeDialog(self)
        if not dlg.exec():
            return
        n1, n2 = dlg.n1, dlg.n2
        if n1 == n2:
            QMessageBox.warning(self, "Same node",
                                "Edge endpoints must be two different node IDs.")
            return
        if n1 not in self.current_generator.model.nodes or \
                n2 not in self.current_generator.model.nodes:
            QMessageBox.warning(self, "Unknown nodes",
                                f"Node {n1} or {n2} not found in the model.")
            return
        cmd = InsertEdgeNodeCommand(n1, n2)
        self.command_manager.execute(cmd, self.current_generator.model)
        if cmd._created_nid is None:
            self._update_status("Edge insert: no elements share that edge.",
                                is_error=True)
            return
        self._update_status(
            f"Inserted node {cmd._created_nid} on edge ({n1}-{n2}); "
            f"split {len(cmd._old_elements)} element(s)."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _auto_orient_normals(self):
        """Run BFS auto-orient on all shell elements for consistent normals."""
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model has no elements.")
            return

        gen = self.current_generator
        shell_count = sum(1 for e in gen.model.elements.values()
                         if e.type in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'))
        if shell_count == 0:
            QMessageBox.information(self, "No Shells",
                                    "No shell elements found in the model.")
            return

        to_flip = gen.compute_auto_orient_flips()
        if not to_flip:
            QMessageBox.information(self, "Already Consistent",
                                    f"All {shell_count} shell element normals "
                                    "are already consistently oriented.")
            return

        normals_were_visible = self.show_normals_check.isChecked()

        cmd = AutoOrientNormalsCommand(list(to_flip))
        self.command_manager.execute(cmd, gen.model)

        QMessageBox.information(
            self, "Normals Oriented",
            f"Flipped {len(to_flip)} of {shell_count} shell elements "
            "for consistent normal orientation.")
        self._update_status(
            f"Auto-oriented normals: {len(to_flip)} elements flipped.")

        self._update_viewer(gen, reset_camera=False)
        if normals_were_visible:
            self._toggle_element_normals_visibility()

    def _on_advanced_selection(self, mode, angle):
        """Handle Select Adjacent / Connected / Grow / Shrink from selection dialog."""
        dialog = self.active_selection_dialog
        if not dialog or dialog.entity_type != 'Element':
            return
        if not self.current_generator:
            return
        seed_eids = list(dialog.selected_ids)
        if not seed_eids:
            self._update_status("Select at least one element first.")
            return
        gen = self.current_generator
        if mode == 'adjacent':
            result = gen.select_adjacent_elements(seed_eids, angle)
            dialog.add_selection(result)
            self._update_status(
                f"Adjacent selection: {len(result)} elements (angle ≤ {angle:.0f}°)")
        elif mode == 'connected':
            result = gen.select_connected_elements(seed_eids)
            dialog.add_selection(result)
            self._update_status(
                f"Connected selection: {len(result)} elements")
        elif mode == 'grow':
            result = gen.grow_element_selection(seed_eids, angle)
            dialog.replace_selection(result)
            added = len(result) - len(seed_eids)
            self._update_status(
                f"Grow selection: {len(result)} elements (+{max(0, added)}, angle ≤ {angle:.0f}°)")
        elif mode == 'shrink':
            result = gen.shrink_element_selection(seed_eids)
            dialog.replace_selection(result)
            removed = len(seed_eids) - len(result)
            self._update_status(
                f"Shrink selection: {len(result)} elements (-{max(0, removed)})")

    def _on_selection_mode_change(self, mode):
        """Handle selection mode change (box/polygon/circle) from dialog."""
        self._pick_selection_mode = mode  # Always store for next picking session
        style = self.plotter.interactor.GetInteractorStyle()
        if isinstance(style, ClickAndDragInteractor):
            style.set_selection_mode(mode)
        else:
            self._update_status(f"Pick mode: {mode.capitalize()}")

    def _on_previous_selection(self, entity_type):
        """v4.0.1: supply the previous selection to whichever surface
        is active. Pre-v4.0.1 this only handled the modal selection
        dialog (``self.active_selection_dialog``) and silently no-op'd
        when the user clicked Previous on the persistent selection bar
        (``self.selection_bar``). Now both surfaces are covered.
        """
        if entity_type == 'Node':
            prev = self._previous_node_selection
        elif entity_type == 'Element':
            prev = self._previous_element_selection
        else:
            prev = set()

        if not prev:
            self._update_status("No previous selection available.")
            return

        # Prefer the modal dialog when it's open; fall back to the bar.
        sink = self.active_selection_dialog
        if sink is None:
            bar = getattr(self, 'selection_bar', None)
            if bar is not None and bar.isVisible():
                sink = bar
        if sink is None:
            self._update_status("No selection surface is active.")
            return
        sink.add_selection(prev)
        self._update_status(f"Restored {len(prev)} previous {entity_type.lower()}(s)")

    # --- Model Check handlers ---

    def _show_mass_properties(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        mass_data = self.current_generator.calculate_mass_summary()
        dlg = MassPropertiesReportDialog(mass_data, self)
        dlg.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dlg.show()

    def _show_free_edge_report(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        free_edges = self.current_generator.find_free_edges()
        node_coords = {nid: node.xyz for nid, node
                       in self.current_generator.model.nodes.items()}
        dlg = FreeEdgeReportDialog(free_edges, node_coords, self)
        dlg.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        def zoom_to(n1, n2):
            try:
                p1 = np.array(node_coords[n1])
                p2 = np.array(node_coords[n2])
                mid = (p1 + p2) / 2.0
                length = float(np.linalg.norm(p2 - p1))
                radius = max(length * 3, 1.0)
                self.plotter.camera.focal_point = mid.tolist()
                self.plotter.camera.position = (mid + np.array([0, 0, radius])).tolist()
                self.plotter.render()
            except (KeyError, TypeError):
                pass

        dlg.zoom_to_edge = zoom_to
        dlg.show()

    def _show_quality_summary(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        quality_data = self.current_generator.calculate_element_quality()
        if not quality_data:
            QMessageBox.information(self, "No Data",
                                    "No shell elements with quality metrics found.")
            return
        dlg = QualitySummaryReportDialog(quality_data, self)
        dlg.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        def select_failing(eids):
            if eids:
                self._highlight_entities('Element', eids)
                self._update_status(
                    f"Highlighted {len(eids)} failing elements.")

        dlg.select_failing_callback = select_failing
        dlg.show()

    def _show_orphan_check(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        orphan_data = self.current_generator.find_orphans()
        dlg = OrphanCheckDialog(orphan_data, self)
        dlg.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dlg.show()

    # --- Connection Modeling handlers ---

    def _create_spider_connection(self):
        """Create a spider connection (center node + RBE2/RBE3)."""
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return

        dlg = CreateSpiderDialog(self)
        if not dlg.exec():
            return

        settings = dlg.get_settings()
        ring_nids = settings['ring_nids']
        if len(ring_nids) < 2:
            QMessageBox.warning(self, "Not Enough Nodes",
                                "Select at least 2 ring nodes.")
            return

        # Validate all ring nodes exist
        model = self.current_generator.model
        missing = [n for n in ring_nids if n not in model.nodes]
        if missing:
            QMessageBox.warning(self, "Invalid Nodes",
                                f"Nodes not found: {missing[:10]}")
            return

        try:
            cmd = AddSpiderCommand(ring_nids, settings['rbe_type'],
                                   settings['dof'])
            self.command_manager.execute(cmd, model)
            QMessageBox.information(
                self, "Spider Created",
                f"Created {settings['rbe_type']} spider with center node "
                f"{cmd._center_nid} and {len(ring_nids)} ring nodes.")
            self._update_status(
                f"Spider {settings['rbe_type']}: center N{cmd._center_nid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def _create_weld_fastener(self):
        """Create a weld/fastener connection (CWELD-like) between two nodes."""
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return

        available_pids = list(self.current_generator.model.properties.keys())
        dlg = CreateWeldDialog(available_pids, self)
        if not dlg.exec():
            return

        settings = dlg.get_settings()
        model = self.current_generator.model

        if settings['nid_a'] not in model.nodes:
            QMessageBox.warning(self, "Invalid Node",
                                f"Node A ({settings['nid_a']}) not found.")
            return
        if settings['nid_b'] not in model.nodes:
            QMessageBox.warning(self, "Invalid Node",
                                f"Node B ({settings['nid_b']}) not found.")
            return

        try:
            cmd = AddWeldCommand(settings['nid_a'], settings['nid_b'],
                                 settings['diameter'], settings['pid'])
            self.command_manager.execute(cmd, model)
            QMessageBox.information(
                self, "Weld Created",
                f"Created weld connector between nodes "
                f"{settings['nid_a']} and {settings['nid_b']}.")
            self._update_status(
                f"Weld: N{settings['nid_a']}-N{settings['nid_b']}")
            self._update_viewer(self.current_generator, reset_camera=False)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def _open_element_split_tool(self):
        """Launches the entity selection dialog to choose shell elements to split."""
        if self.selection_bar.isVisible():
            return

        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Plates", "The model contains no plates to split.")
            return

        all_shell_eids = [
            eid for eid, elem in self.current_generator.model.elements.items()
            if elem.type in ['CQUAD4', 'CTRIA3']
        ]
        if not all_shell_eids:
            QMessageBox.information(self, "No Plates",
                                    "There are no CQUAD4 or CTRIA3 plates to split.")
            return

        self._start_selection('Element', all_shell_eids,
                              self._on_split_elements_accept,
                              self._build_select_by_data())

    def _on_split_elements_accept(self):
        """Handle the accept signal from the split elements selection dialog."""
        selected_eids = self.selection_bar.get_selected_ids()
        self._end_selection_mode()

        if not selected_eids:
            self._update_status("Split cancelled: No elements selected.")
            return

        # Pure calculation - no model mutation yet
        result = self.current_generator.split_shell_elements(selected_eids)

        if not result['valid_eids']:
            QMessageBox.warning(self, "No Valid Elements",
                                "None of the selected elements could be split.")
            return

        # Build compound command: delete originals → add nodes → add elements
        sub_cmds = [DeleteElementsCommand(result['valid_eids'])]
        sub_cmds.append(AddNodesCommand(result['new_node_data']))
        for params in result['new_element_params']:
            sub_cmds.append(AddElementCommand('plate', params))

        n_orig = len(result['valid_eids'])
        n_new = len(result['new_element_params'])
        compound = CompoundCommand(sub_cmds, f"Split {n_orig} element(s) into {n_new}")
        self.command_manager.execute(compound, self.current_generator.model)

        self._update_status(f"Split {n_orig} element(s) into {n_new}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _toggle_element_normals_visibility(self):
        """
        Shows, hides, or redraws the element normals actor based on UI controls.
        This single handler manages the checkbox, menu item, and size input.
        """
        # --- 1. Synchronize the UI controls (Checkbox and Menu Item) ---
        # We block signals to prevent a recursive loop where a change to one
        # control triggers a signal that then tries to change the first one back.
        self.show_normals_action.blockSignals(True)
        self.show_normals_check.blockSignals(True)

        # Determine the correct checked state. The widget that sent the signal
        # (the sender) dictates the new state.
        is_checked = self.show_normals_check.isChecked()
        if self.sender() in [self.show_normals_action, self.show_normals_check]:
            is_checked = self.sender().isChecked()

        # Set both controls to the same state
        self.show_normals_check.setChecked(is_checked)
        self.show_normals_action.setChecked(is_checked)
        
        # Unblock signals now that they are in sync
        self.show_normals_action.blockSignals(False)
        self.show_normals_check.blockSignals(False)
        
        # --- 2. Update the 3D Plot ---
        # The simplest and most robust way to handle both visibility and size changes
        # is to always remove the old actor and recreate it if needed.
        self.plotter.remove_actor('element_normals_actor', render=False)

        # If the controls are checked on, call the method to create the new actor.
        if is_checked:
            self._create_element_normals_actor()

        # Render the final scene
        self.plotter.render()

    def _toggle_free_edges_visibility(self):
        """Shows or hides the free edge overlay. Syncs checkbox and menu action."""
        self.show_free_edges_action.blockSignals(True)
        self.show_free_edges_check.blockSignals(True)

        is_checked = self.show_free_edges_check.isChecked()
        if self.sender() in [self.show_free_edges_action, self.show_free_edges_check]:
            is_checked = self.sender().isChecked()

        self.show_free_edges_check.setChecked(is_checked)
        self.show_free_edges_action.setChecked(is_checked)

        self.show_free_edges_action.blockSignals(False)
        self.show_free_edges_check.blockSignals(False)

        self.plotter.remove_actor('free_edges_actor', render=False)
        if is_checked:
            self._create_free_edges_actor()
        self.plotter.render()

    def _show_about_dialog(self):
        # v5.0.0 item 21: pull the version tag from __version__ instead
        # of a hardcoded literal so the dialog never drifts from the
        # package version again.
        from node_runner import __version__ as _nr_version
        title = f"About Node Runner v{_nr_version}"
        text = """
        <style>
            body { font-family: Segoe UI, sans-serif; }
            p, li {
                font-size: 13px;
                line-height: 1.3;
            }
            td {
                font-size: 13px;
                line-height: 1.3;
                vertical-align: top;
                padding: 2px 6px;
            }
            h2 {
                color: #f38ba8;
                margin-bottom: 4px;
                margin-top: 18px;
                font-size: 16px;
            }
            h3 {
                color: #89b4fa;
                margin-bottom: 4px;
                margin-top: 14px;
                border-bottom: 1px solid #45475a;
                padding-bottom: 3px;
                font-size: 14px;
            }
            h4 {
                color: #a6e3a1;
                margin-bottom: 2px;
                margin-top: 10px;
                font-size: 13px;
            }
            .section { margin-left: 12px; margin-bottom: 6px; }
            .subtle { font-size: 12px; color: #a6adc8; }
            .tag {
                background-color: #313244;
                border: 1px solid #45475a;
                border-radius: 3px;
                padding: 1px 5px;
                font-size: 12px;
                color: #cdd6f4;
            }
            .warn { color: #fab387; }
            .card {
                background-color: #1e1e2e;
                border: 1px solid #45475a;
                border-radius: 4px;
                padding: 6px 10px;
                margin: 3px 0;
            }
            ul { margin-top: 2px; margin-bottom: 4px; }
            li { margin-bottom: 1px; }
        </style>

        <pre style="color: #cdd6f4; font-size: 11px; line-height: 1.1; margin: 0; font-family: 'Cascadia Code', 'Consolas', monospace;">
    _   __          __        ____
   / | / /___  ____/ /__     / __ \__  ______  ____  ___  _____
  /  |/ / __ \/ __  / _ \   / /_/ / / / / __ \/ __ \/ _ \/ ___/
 / /|  / /_/ / /_/ /  __/  / _, _/ /_/ / / / / / / /  __/ /
/_/ |_/\____/\__,_/\___/  /_/ |_|\__,_/_/ /_/_/ /_/\___/_/</pre>
        <p style="margin-top: 4px;"><span class="tag">v__NR_VERSION__</span></p>
        <p class="subtle">Created by Angel Linares<br>Escape Velocity Ventures, LLC</p>
        <hr>

        <!-- ============================================================ -->
        <h2>What Is Node Runner?</h2>
        <!-- ============================================================ -->

        <div class="section">
        <p>
            <b>Node Runner</b> is a standalone, open-source finite element
            model pre- and post-processor purpose-built for the
            <b>Nastran</b> solver family. It is designed from the ground up
            to give stress engineers a fast, visual,
            keyboard-and-mouse-driven environment for building, inspecting,
            editing, and exporting Nastran Bulk Data Files
            (<b>.bdf</b> / <b>.dat</b>), as well as viewing solver results
            from OP2 files.
        </p>
        <p>
            It is <i>not</i> a solver. It handles pre-processing and
            post-processing &mdash; and aims to do both well.
        </p>
        <p>
            Node Runner sits in the gap between heavyweight commercial
            packages (Patran, HyperMesh) and manual text editing of
            BDF files. It offers a graphical, interactive experience without
            the licensing cost, installation complexity, or learning curve
            of the full-featured tools &mdash; while still being far more
            productive than hand-editing bulk data cards in a text editor.
        </p>
        </div>

        <!-- ============================================================ -->
        <h2>Who Is It For?</h2>
        <!-- ============================================================ -->

        <div class="section">
        <p>Node Runner is aimed at:</p>
        <ul>
            <li><b>Stress engineers</b> in aerospace MRO, OEM, and STC
                organizations who need to build or modify small-to-medium
                Nastran models (panels, doublers, repairs, fuselage
                sections, brackets, fittings).</li>
            <li><b>Students and researchers</b> learning Nastran who want
                to visualize what their BDF cards actually produce, without
                needing a commercial pre-processor license.</li>
            <li><b>Anyone</b> who regularly edits <tt>.bdf</tt> files by
                hand and wants a faster, less error-prone workflow with
                real-time 3D feedback.</li>
        </ul>
        <p>
            If your workflow involves opening a BDF in a text editor,
            changing a node coordinate, saving, re-running, and hoping
            the mesh looks right &mdash; Node Runner is for you.
        </p>
        </div>

        <!-- ============================================================ -->
        <h2>Detailed Capabilities</h2>
        <!-- ============================================================ -->

        <!-- ---- File I/O ---- -->
        <h3>File Handling &amp; Solver Compatibility</h3>
        <div class="section">
            <table>
            <tr><td><b>Open / Import:</b></td>
                <td>Reads standard Nastran Bulk Data Files
                    (<tt>.bdf</tt>, <tt>.dat</tt>, <tt>.nas</tt>)
                    in <b>small field</b> (8-char),
                    <b>large field</b> (16-char), and <b>free field</b>
                    (comma-delimited) formats. Import can either replace
                    the current model or append cards to it.</td></tr>
            <tr><td><b>Save / Export:</b></td>
                <td>Writes to clean, standard <b>small-field BDF</b>,
                    optionally including a Case Control Deck with
                    subcases, output requests, and load/SPC
                    references.</td></tr>
            <tr><td><b>Solver Support:</b></td>
                <td>Output is compatible with <b>MSC Nastran</b>,
                    <b>Siemens NX Nastran</b>, and the open-source
                    <b>MYSTRAN</b> solver (new in v5.0.0. see the
                    integration block below). No solver-specific
                    extensions are used in the generic export
                    target.</td></tr>
            </table>
        </div>

        <!-- ============================================================ -->
        <h3>Post-Processing + Scoped Analysis <span class="tag">new in v5.1.0</span></h3>
        <!-- ============================================================ -->

        <div class="section">
        <p>
            v5.1.0 turns the existing pre/solve/post baseline into a
            full full FEA workflow:
        </p>
        <table>
        <tr><td><b>Results tab:</b></td>
            <td>A fifth sidebar tab with a tree of subcases &rarr;
                output vectors &rarr; components. Right-click for
                <i>Contour / Vector / Apply as deformed shape /
                Animate / Copy values to CSV</i>. Controls panel below
                the tree: Auto/Manual color range, level count
                (2&ndash;32), Show Min/Max Markers, Show Element Value
                Labels (with top-N filter), deformed-shape toggle and
                scale.</td></tr>
        <tr><td><b>AnalysisSet scope:</b></td>
            <td>AnalysisSets now carry a <i>solver_target</i> (MYSTRAN /
                MSC / NX / Generic), a <i>group_target</i>, explicit
                LOAD / SPC SID white-lists, and a per-set field
                format. <i>Analysis &rarr; Run Analysis (MYSTRAN)</i>
                requires an active set so the scope is always
                explicit. Right-click an AnalysisSet on the sidebar
                for <i>Run / Export / Set Active / Edit / Duplicate /
                Delete</i>.</td></tr>
        <tr><td><b>Scoped exports:</b></td>
            <td>Exports and MYSTRAN runs build a shallow-copy of the
                BDF scoped to the active AnalysisSet (group + load +
                SPC filtering + auto-collected dependencies) so the
                resulting deck is self-consistent. Source model is
                never mutated.</td></tr>
        <tr><td><b>Save BDF redesign:</b></td>
            <td><i>File &rarr; Save BDF&hellip;</i> is now a single
                richer dialog: output path, solver target, field
                format, scope picker, and an embedded
                <code>$ NR-META v1</code> metadata block (groups + tree
                state) that survives a MSC/NX/MYSTRAN round-trip as
                ordinary Nastran comments. Sidecar <code>.nrmeta</code>
                file available as an opt-in.</td></tr>
        <tr><td><b>Groups carry dependencies:</b></td>
            <td>Right-click a group &rarr; <i>Collect Dependencies</i>
                auto-fills the group's <i>properties / materials /
                coords</i> from the elements it contains. Required for
                the new group-scoped exports.</td></tr>
        </table>
        </div>

        <!-- ============================================================ -->
        <h3>Post-Processing + Scoped Analysis <span class="tag">new in v5.1.0</span></h3>
        <!-- ============================================================ -->

        <div class="section">
        <p>
            v5.1.0 turns the existing pre/solve/post baseline into a
            full full FEA workflow:
        </p>
        <table>
        <tr><td><b>Results tab:</b></td>
            <td>A fifth sidebar tab with a tree of subcases &rarr;
                output vectors &rarr; components. Right-click for
                <i>Contour / Vector / Apply as deformed shape /
                Animate / Copy values to CSV</i>. Controls panel below
                the tree: Auto/Manual color range, level count
                (2&ndash;32), Show Min/Max Markers, Show Element Value
                Labels (with top-N filter), deformed-shape toggle and
                scale.</td></tr>
        <tr><td><b>AnalysisSet scope:</b></td>
            <td>AnalysisSets now carry a <i>solver_target</i> (MYSTRAN /
                MSC / NX / Generic), a <i>group_target</i>, explicit
                LOAD / SPC SID white-lists, and a per-set field
                format. <i>Analysis &rarr; Run Analysis (MYSTRAN)</i>
                requires an active set so the scope is always
                explicit. Right-click an AnalysisSet on the sidebar
                for <i>Run / Export / Set Active / Edit / Duplicate /
                Delete</i>.</td></tr>
        <tr><td><b>Scoped exports:</b></td>
            <td>Exports and MYSTRAN runs build a shallow-copy of the
                BDF scoped to the active AnalysisSet (group + load +
                SPC filtering + auto-collected dependencies) so the
                resulting deck is self-consistent. Source model is
                never mutated.</td></tr>
        <tr><td><b>Save BDF redesign:</b></td>
            <td><i>File &rarr; Save BDF&hellip;</i> is now a single
                richer dialog: output path, solver target, field
                format, scope picker, and an embedded
                <code>$ NR-META v1</code> metadata block (groups + tree
                state) that survives a MSC/NX/MYSTRAN round-trip as
                ordinary Nastran comments. Sidecar <code>.nrmeta</code>
                file available as an opt-in.</td></tr>
        <tr><td><b>Groups carry dependencies:</b></td>
            <td>Right-click a group &rarr; <i>Collect Dependencies</i>
                auto-fills the group's <i>properties / materials /
                coords</i> from the elements it contains. Required for
                the new group-scoped exports.</td></tr>
        </table>
        </div>

        <!-- ============================================================ -->
        <h3>MYSTRAN Solver Integration <span class="tag">new in v5.0.0</span></h3>
        <!-- ============================================================ -->

        <div class="section">
        <p>
            Node Runner v5.0.0 adds an end-to-end
            <b>pre &rarr; solve &rarr; post</b> workflow that drives the
            free, open-source <a href="https://www.mystran.com/">MYSTRAN</a>
            binary as a subprocess. Configure the executable path in
            <i>Preferences &rarr; MYSTRAN</i>, then click
            <i>Analysis &rarr; Run Analysis (MYSTRAN)</i>&hellip;
            (<tt>Ctrl+R</tt>).
        </p>
        <table>
        <tr><td><b>Supported solutions:</b></td>
            <td><span class="tag">SOL 101</span> Linear Static,
                <span class="tag">SOL 103</span> Normal Modes,
                <span class="tag">SOL 105</span> Linear Buckling.
                Same numbering as MSC Nastran. SOL 106 (nonlinear),
                aero, transient, frequency response, and optimization
                solutions are <span class="warn">not supported</span>
                by MYSTRAN and are caught by the pre-flight scanner.
            </td></tr>
        <tr><td><b>Pre-flight scanner:</b></td>
            <td>Walks the in-memory model in &lt; 1 s and flags any
                cards MYSTRAN can&rsquo;t run (aero, nonlinear,
                contact, optimization, unsupported elements). Blocking
                issues stop the run; warnings let the export translator
                drop the offending cards.</td></tr>
        <tr><td><b>Export translator:</b></td>
            <td>Re-writes the deck for MYSTRAN: strips MSC/NX-only
                <tt>PARAM</tt>s (POST, COUPMASS, K6ROT, AUTOSPC, &hellip;),
                injects MYSTRAN-required PARAMs (SOLLIB, QUAD4TYP,
                optionally WTMASS / GRDPNT), and drops blocking card
                families. Output lands in a timestamped scratch
                directory; your source decks are never modified.</td></tr>
        <tr><td><b>Run UX:</b></td>
            <td>Non-blocking progress dialog with stage tracking from
                MYSTRAN&rsquo;s stdout (Reading bulk data &rarr;
                Assembling K &rarr; Factorizing &rarr; Solving &rarr;
                Recovering element results). Cancel button cleanly
                terminates the child process. Post-run summary dialog
                shows the run folder path, the <tt>.ERR</tt> file
                contents inline, and Open Folder buttons.</td></tr>
        <tr><td><b>Results adapter:</b></td>
            <td>OP2 binary first (when MYSTRAN emits one), F06 text
                fallback for displacements + eigenvalues. Results land
                in the existing Result Browser / Vector Overlay /
                Animation Timeline docks &mdash; same UX as loading a
                vendor OP2.</td></tr>
        <tr><td><b>Tooltips throughout:</b></td>
            <td>Every MYSTRAN-specific control (SOLLIB, QUAD4TYP,
                WTMASS, output requests, SOL radio, etc.) carries a
                hover tip that cross-references the MSC Nastran or
                common-tool equivalent.</td></tr>
        <tr><td><b>Analysis history:</b></td>
            <td><i>Analysis &rarr; Analysis History</i>&hellip; browses
                prior runs in the scratch directory; double-click
                reloads a run&rsquo;s OP2 / F06 results into the
                current model.</td></tr>
        </table>
        <p class="subtle">
            Node Runner does <b>not</b> bundle the MYSTRAN binary &mdash;
            download it from the upstream project and point Node Runner
            at it via <i>Preferences &rarr; MYSTRAN &rarr; MYSTRAN
            executable</i>.
        </p>
        </div>

        <!-- ---- Elements ---- -->
        <h3>Supported Elements</h3>
        <div class="section">
            <h4>Shell Elements</h4>
            <table>
            <tr><td><span class="tag">CQUAD4</span></td>
                <td>4-node isoparametric quadrilateral shell</td></tr>
            <tr><td><span class="tag">CTRIA3</span></td>
                <td>3-node triangular shell</td></tr>
            <tr><td><span class="tag">CMEMBRAN</span></td>
                <td>Membrane element (in-plane only, no bending)</td></tr>
            </table>

            <h4>Line Elements</h4>
            <table>
            <tr><td><span class="tag">CBEAM</span></td>
                <td>General beam with 6 DOF per node, orientation vector,
                    and WA/WB end offsets (OFFT)</td></tr>
            <tr><td><span class="tag">CBAR</span></td>
                <td>Simple beam/bar element with 6 DOF per node</td></tr>
            <tr><td><span class="tag">CROD</span></td>
                <td>Axial rod (tension/compression + torsion only)</td></tr>
            </table>

            <h4>Solid Elements</h4>
            <table>
            <tr><td><span class="tag">CHEXA</span></td>
                <td>8-node hexahedral solid</td></tr>
            <tr><td><span class="tag">CTETRA</span></td>
                <td>4-node tetrahedral solid</td></tr>
            <tr><td><span class="tag">CPENTA</span></td>
                <td>6-node pentahedral (wedge) solid</td></tr>
            </table>

            <h4>Connector &amp; Special Elements</h4>
            <table>
            <tr><td><span class="tag">RBE2</span></td>
                <td>Rigid body element &mdash; single independent node,
                    multiple dependent nodes</td></tr>
            <tr><td><span class="tag">RBE3</span></td>
                <td>Interpolation element &mdash; distributes loads from
                    dependent to independent nodes by weighted
                    averaging</td></tr>
            <tr><td><span class="tag">CBUSH</span></td>
                <td>Generalized spring/damper element with 6 independent
                    stiffness components</td></tr>
            <tr><td><span class="tag">CSHEAR</span></td>
                <td>4-node shear panel</td></tr>
            <tr><td><span class="tag">CGAP</span></td>
                <td>Contact/gap element with axial stiffness and
                    preload</td></tr>
            <tr><td><span class="tag">CONM2</span></td>
                <td>Concentrated mass with optional 6&times;6 inertia
                    matrix</td></tr>
            <tr><td><span class="tag">PLOTEL</span></td>
                <td>Visualization-only element (ignored by solver)</td></tr>
            </table>
        </div>

        <!-- ---- Properties ---- -->
        <h3>Supported Properties</h3>
        <div class="section">
            <table>
            <tr><td><span class="tag">PSHELL</span></td>
                <td>Homogeneous shell &mdash; thickness, material, NSM</td></tr>
            <tr><td><span class="tag">PCOMP</span></td>
                <td>Composite layup &mdash; per-ply MID, thickness, fiber
                    angle, stress output flag. Includes a dedicated
                    <b>composite layup editor</b> with add/remove/reorder
                    ply controls and failure theory selection
                    (Hill, Hoffman, Tsai-Wu, Max Strain).</td></tr>
            <tr><td><span class="tag">PBEAM</span></td>
                <td>Beam cross-section &mdash; A, I1, I2, J, MID</td></tr>
            <tr><td><span class="tag">PBAR</span></td>
                <td>Bar cross-section &mdash; A, I1, I2, J, MID</td></tr>
            <tr><td><span class="tag">PROD</span></td>
                <td>Rod cross-section &mdash; A, J, stress recovery
                    coefficient</td></tr>
            <tr><td><span class="tag">PBUSH</span></td>
                <td>Bush/spring stiffness &mdash; K1&ndash;K6 (3
                    translational + 3 rotational)</td></tr>
            <tr><td><span class="tag">PSOLID</span></td>
                <td>Solid property &mdash; material reference</td></tr>
            <tr><td><span class="tag">PSHEAR</span></td>
                <td>Shear panel property &mdash; thickness, MID</td></tr>
            <tr><td><span class="tag">PGAP</span></td>
                <td>Gap property &mdash; initial gap, preload, axial/contact
                    stiffness</td></tr>
            </table>
        </div>

        <!-- ---- Materials ---- -->
        <h3>Supported Materials</h3>
        <div class="section">
            <table>
            <tr><td><span class="tag">MAT1</span></td>
                <td><b>Isotropic</b> &mdash; E, G, &nu;, &rho;, &alpha;,
                    T<sub>ref</sub>, g<sub>e</sub></td></tr>
            <tr><td><span class="tag">MAT8</span></td>
                <td><b>2D Orthotropic</b> (composites) &mdash; E1, E2,
                    &nu;12, G12, G1z, G2z, &rho;, &alpha;1,
                    &alpha;2</td></tr>
            <tr><td><span class="tag">MAT9</span></td>
                <td><b>3D Anisotropic</b> &mdash; full 6&times;6 symmetric
                    stiffness matrix (21 independent G<sub>ij</sub>
                    terms)</td></tr>
            </table>
        </div>

        <!-- ---- Loads ---- -->
        <h3>Loads &amp; Constraints</h3>
        <div class="section">
            <h4>Loads</h4>
            <table>
            <tr><td><span class="tag">FORCE</span></td>
                <td>Concentrated force at nodes (Fx, Fy, Fz)</td></tr>
            <tr><td><span class="tag">MOMENT</span></td>
                <td>Concentrated moment at nodes (Mx, My, Mz)</td></tr>
            <tr><td><span class="tag">PLOAD4</span></td>
                <td>Surface pressure on shell elements</td></tr>
            <tr><td><span class="tag">GRAV</span></td>
                <td>Gravity / body acceleration &mdash; direction vector
                    &times; scale, with CID support</td></tr>
            <tr><td><span class="tag">TEMPD</span></td>
                <td>Default uniform temperature field</td></tr>
            <tr><td><span class="tag">LOAD</span></td>
                <td>Load combination card &mdash; linearly combine
                    multiple load sets with individual scale factors and
                    an overall multiplier</td></tr>
            </table>

            <h4>Constraints</h4>
            <table>
            <tr><td><span class="tag">SPC</span></td>
                <td>Single-point constraint &mdash; fix any combination of
                    the 6 DOFs (TX, TY, TZ, RX, RY, RZ) on selected
                    nodes</td></tr>
            </table>
        </div>

        <!-- ---- Coordinate Systems ---- -->
        <h3>Coordinate Systems</h3>
        <div class="section">
            <table>
            <tr><td><span class="tag">CORD2R</span></td>
                <td>Rectangular (Cartesian) coordinate system</td></tr>
            <tr><td><span class="tag">CORD2C</span></td>
                <td>Cylindrical coordinate system</td></tr>
            <tr><td><span class="tag">CORD2S</span></td>
                <td>Spherical coordinate system</td></tr>
            </table>
            <p class="subtle">Define by three points (origin, Z-axis, XZ-plane)
                or by translate-and-rotate from a reference system. Global
                Rectangular (CID 0), Cylindrical (CID 1), and Spherical
                (CID 2) are created automatically.</p>
        </div>

        <!-- ---- Case Control ---- -->
        <h3>Case Control &amp; Analysis Setup</h3>
        <div class="section">
        <p>
            Define one or more <b>subcases</b> through the Subcase Editor.
            Each subcase specifies:
        </p>
        <ul>
            <li>A <b>LOAD</b> set reference (any load SID or LOAD
                combination)</li>
            <li>An <b>SPC</b> set reference</li>
            <li>Output requests: <b>DISPLACEMENT</b>, <b>STRESS</b>,
                <b>FORCE</b>, <b>STRAIN</b> (each set to NONE or
                ALL)</li>
            <li>A <b>METHOD</b> reference (for eigenvalue extraction)</li>
            <li>A <b>STATSUB</b> reference (for buckling pre-load
                subcase)</li>
        </ul>

        <h4>Solution Types</h4>
        <p>
            Set the analysis type via the <b>Analysis</b> menu:
        </p>
        <ul>
            <li><span class="tag">SOL 101</span> &mdash; Linear Static</li>
            <li><span class="tag">SOL 103</span> &mdash; Normal Modes (real eigenvalue)</li>
            <li><span class="tag">SOL 105</span> &mdash; Linear Buckling
                (requires a static pre-load subcase referenced via STATSUB
                and an eigenvalue extraction subcase with METHOD pointing to
                an EIGRL card)</li>
        </ul>

        <h4>Eigenvalue Extraction</h4>
        <table>
        <tr><td><span class="tag">EIGRL</span></td>
            <td>Lanczos method eigenvalue extraction card &mdash; define
                SID, frequency range (V1/V2), and number of desired
                roots (ND). Create and manage via
                <b>Analysis &gt; EIGRL</b>.</td></tr>
        </table>

        <p>
            The Case Control Deck is automatically written at the top of
            the BDF file on save. If no subcases are defined, no case
            control section is written (bulk-data-only output).
        </p>
        </div>

        <!-- ---- Editing ---- -->
        <h3>Editing &amp; Interaction</h3>
        <div class="section">
            <b>Full CRUD Lifecycle:</b> Create, view, edit, and delete all
            major entity types (nodes, elements, properties, materials,
            loads, constraints, coordinate systems, mass elements, plot
            elements).<br><br>

            <b>Undo / Redo:</b> All model mutations are tracked by a
            command-pattern undo system (Ctrl+Z / Ctrl+Y) with a 20-step
            history. Every operation &mdash; from adding a single node to
            renumbering the entire model &mdash; can be reversed.<br><br>

            <b>Selection Tools:</b> Pick entities directly from the 3D
            viewer, or use the Entity Selection Dialog to select by ID,
            ID range, or a pasted list. Advanced selection modes include
            <b>By Property</b>, <b>By Material</b>, <b>By Element Type</b>,
            and <b>By Shape</b>. The Find Entities tool can search and
            zoom to any entity by ID.<br><br>

            <b>Node Transform Tool:</b> Translate (delta), Move-to-Point,
            Scale, and Rotate selected nodes. Supports global origin,
            selection centroid, or a custom point as the transformation
            center.<br><br>

            <b>Element Editor:</b> Load any element by EID to edit its
            property assignment, node connectivity, beam orientation
            (vector / G0 node / CID), and CBEAM end offsets
            (WA, WB, OFFT).<br><br>

            <b>Find/Replace Property &amp; Material:</b> Bulk-reassign the
            property ID across all matching elements, or the material ID
            across all matching properties, with a live preview of affected
            counts.<br><br>

            <b>Renumber Nodes &amp; Elements:</b> Sequentially renumber all
            node and/or element IDs starting from a user-defined value.
            All internal references (element connectivity, loads, SPCs) are
            automatically remapped. Fully undoable.
        </div>

        <!-- ---- Visualization ---- -->
        <h3>3D Visualization</h3>
        <div class="section">
            <b>Rendering Engine:</b> Real-time 3D powered by
            <b>PyVista</b> (VTK). Handles tens of thousands of elements
            interactively.<br><br>

            <b>Rendering Styles:</b> Wireframe or Surface-with-Edges, with
            an optional transparent shell mode.<br><br>

            <b>Coloring Modes:</b><br>
            <ul>
                <li><b>By Element Type</b> &mdash; shells, beams, RBEs,
                    masses, etc. each get a distinct color</li>
                <li><b>By Property ID</b> &mdash; each PID mapped to a
                    unique color (default mode)</li>
                <li><b>By Element Quality</b> &mdash; heat-map visualization
                    of mesh quality metrics (see below)</li>
                <li><b>Results</b> &mdash; contour plot of loaded OP2
                    results (displacement, stress, eigenvectors, SPC
                    forces) with a jet colormap and scalar bar</li>
            </ul>

            <b>Overlay Actors:</b> Load vectors (with relative or uniform
            scaling), SPC constraint markers, RBE spider plots, mass
            element glyphs, beam/bush orientation vectors, element normal
            arrows, free-edge highlights, and PLOTEL lines &mdash; each
            independently togglable from the model tree.<br><br>

            <b>Standard Views:</b> Top, Bottom, Front, Back, Left, Right
            with one-click access. Perspective / Orthographic toggle.<br><br>

            <b>Node &amp; Element Labels:</b> Toggle numeric ID labels
            on all nodes and/or elements.<br><br>

            <b>Color Manager:</b> Customize or randomize the color
            palette for any coloring mode.<br><br>

            <b>Themes:</b> Dark mode (Catppuccin Mocha) and Light mode.
        </div>

        <!-- ---- Quality Metrics ---- -->
        <h3>Element Quality Metrics</h3>
        <div class="section">
        <p>
            Real-time element quality analysis with per-element heat-map
            visualization. Select a metric from the dropdown to color
            the mesh by quality:
        </p>
            <h4>Quad Metrics (CQUAD4 / CMEMBRAN)</h4>
            <table>
            <tr><td><b>Warping</b></td>
                <td>Out-of-plane angle between triangle normals of a
                    split quad. 0&deg; = perfectly planar.</td></tr>
            <tr><td><b>Aspect Ratio</b></td>
                <td>Longest edge / shortest edge. 1.0 = square.</td></tr>
            <tr><td><b>Skew</b></td>
                <td>Deviation of the diagonal crossing angle from
                    90&deg;. 0&deg; = perfectly orthogonal
                    diagonals.</td></tr>
            <tr><td><b>Jacobian Ratio</b></td>
                <td>min(J)/max(J) evaluated at the four isoparametric
                    corners. 1.0 = no distortion.</td></tr>
            <tr><td><b>Taper</b></td>
                <td>Ratio of triangle areas when quad is split along
                    each diagonal. 1.0 = symmetric.</td></tr>
            <tr><td><b>Min / Max Angle</b></td>
                <td>Interior corner angles. Ideal quad = 90&deg; at
                    every corner.</td></tr>
            </table>

            <h4>Triangle Metrics (CTRIA3)</h4>
            <table>
            <tr><td><b>Aspect Ratio</b></td>
                <td>Longest edge / shortest edge. 1.0 = equilateral.</td></tr>
            <tr><td><b>Skew</b></td>
                <td>Max deviation of any interior angle from the ideal
                    60&deg;.</td></tr>
            <tr><td><b>Min / Max Angle</b></td>
                <td>Interior corner angles. Ideal triangle = 60&deg;
                    at every corner.</td></tr>
            </table>
        </div>

        <!-- ---- Geometry ---- -->
        <h3>Geometry Engine</h3>
        <div class="section">
        <p>
            Create and manage parametric geometry primitives that exist
            independently of the FE mesh. Geometry entities are displayed
            as colored overlays (orange points, magenta lines, cyan
            arcs/circles) and can be meshed into nodes and elements.
        </p>
            <h4>Primitives</h4>
            <table>
            <tr><td><b>Points</b></td>
                <td>Define by X, Y, Z coordinates. Supports single-point
                    creation or batch import from Excel
                    (<b>From Excel</b> button &mdash; columns: ID, CSys,
                    X, Y, Z). Batch import is a single undoable
                    operation.</td></tr>
            <tr><td><b>Lines</b></td>
                <td>Straight line between two geometry points.</td></tr>
            <tr><td><b>Arcs</b></td>
                <td>Circular arc through three geometry points
                    (start, mid, end).</td></tr>
            <tr><td><b>Circles</b></td>
                <td>Full circle defined by a center point and
                    radius.</td></tr>
            <tr><td><b>Surfaces</b></td>
                <td>Bounded region defined by a closed loop of boundary
                    curves.</td></tr>
            </table>

            <h4>Geometry Meshing</h4>
            <table>
            <tr><td><b>Mesh Curve</b></td>
                <td>Discretize a line or arc into evenly spaced nodes
                    connected by line elements (CBAR, CBEAM, or CROD)
                    with a user-specified number of segments and
                    property ID.</td></tr>
            <tr><td><b>Nodes at Points</b></td>
                <td>Create FE grid nodes at the locations of selected
                    geometry points.</td></tr>
            </table>

            <h4>Deletion</h4>
            <p>
                Delete geometry points, curves, or surfaces via the
                <b>Delete</b> menu. Deletion cascades automatically:
                deleting a point also removes dependent curves and their
                dependent surfaces. All deletions are fully undoable.
            </p>
        </div>

        <!-- ---- Mesh Tools ---- -->
        <h3>Mesh &amp; Model Tools</h3>
        <div class="section">
            <table>
            <tr><td><b>Coincident Node Checker</b></td>
                <td>Detect duplicate nodes within a user-defined
                    tolerance. Review results, remove false positives,
                    then merge. All connectivity is remapped
                    automatically.</td></tr>
            <tr><td><b>Duplicate Element Checker</b></td>
                <td>Find elements that share the exact same node set.
                    Review and batch-delete duplicates.</td></tr>
            <tr><td><b>Flip Element Normals</b></td>
                <td>Reverse the winding order of selected shell elements
                    to flip their normal direction.</td></tr>
            <tr><td><b>Split Plates</b></td>
                <td>1-to-4 subdivision: CQUAD4&rarr;4 CQUAD4s,
                    CTRIA3&rarr;4 CTRIA3s. Shared edges reuse midside
                    nodes.</td></tr>
            <tr><td><b>Measure Distance</b></td>
                <td>Pick two nodes to compute and display the
                    straight-line distance between them, with a visual
                    line in the 3D viewer.</td></tr>
            <tr><td><b>Measure Angle</b></td>
                <td>Pick three nodes to compute and display the angle at
                    the central node.</td></tr>
            <tr><td><b>Mass &amp; CG Summary</b></td>
                <td>Computes total mass, structural vs. CONM2 mass
                    breakdown, center of gravity, and per-property mass
                    distribution. Accounts for shells
                    (area&times;t&times;&rho;), beams/bars
                    (A&times;L&times;&rho;), rods, and CONM2
                    elements.</td></tr>
            <tr><td><b>Model Summary</b></td>
                <td>Entity count overview: nodes, elements by type,
                    properties, materials, loads, SPCs, coordinate
                    systems.</td></tr>
            <tr><td><b>Free Edge Display</b></td>
                <td>Highlight shell element edges that belong to only one
                    element &mdash; useful for finding mesh gaps and
                    boundaries.</td></tr>
            </table>
        </div>

        <!-- ---- Groups ---- -->
        <h3>Group Management</h3>
        <div class="section">
        <p>
            Organize nodes and elements into named groups for targeted
            visibility control:
        </p>
        <ul>
            <li><b>Auto-Group by Property</b> &mdash; one-click grouping
                based on PID assignment</li>
            <li><b>Manual Groups</b> &mdash; create, rename, delete, and
                populate groups from the current selection</li>
            <li><b>Isolate Mode</b> &mdash; hide everything except the
                selected group, then restore with one click</li>
            <li><b>Visibility Checkboxes</b> &mdash; toggle individual
                groups on/off in the model tree</li>
        </ul>
        </div>

        <!-- ---- Post-Processing ---- -->
        <h3>Post-Processing (OP2 Results)</h3>
        <div class="section">
        <p>
            Load Nastran OP2 result files via <b>File &gt; Load Results
            (OP2)</b> to visualize solver output directly on the model.
        </p>
            <h4>Supported Result Types</h4>
            <table>
            <tr><td><b>Displacement</b></td>
                <td>Translational (T1, T2, T3) and rotational (R1, R2, R3)
                    displacements per node. Components: X, Y, Z, magnitude,
                    RX, RY, RZ.</td></tr>
            <tr><td><b>Stress</b></td>
                <td>Element-centroid shell stresses for CQUAD4 and CTRIA3.
                    Components: Von Mises, max/min principal,
                    &sigma;<sub>xx</sub>, &sigma;<sub>yy</sub>,
                    &tau;<sub>xy</sub>.</td></tr>
            <tr><td><b>Eigenvector</b></td>
                <td>Modal displacement shapes from SOL 103 or SOL 105,
                    with eigenvalue and frequency metadata.</td></tr>
            <tr><td><b>SPC Forces</b></td>
                <td>Reaction forces at constrained nodes (6 DOF).</td></tr>
            </table>

            <h4>Visualization Features</h4>
            <ul>
                <li><b>Contour plots</b> with jet colormap and scalar bar
                    showing min/max values</li>
                <li><b>Deformed shape display</b> with adjustable
                    deformation scale factor</li>
                <li><b>Subcase selector</b> to switch between load cases
                    or modes</li>
                <li><b>Component selector</b> to choose which result
                    component to display</li>
            </ul>
        </div>

        <!-- ---- Fuselage Generator ---- -->
        <h3>Fuselage Generator</h3>
        <div class="section">
        <p>
            A parametric tool that generates a complete cylindrical fuselage
            finite element model from a small set of inputs:
        </p>
        <ul>
            <li>Fuselage radius, number of bays, frame spacing</li>
            <li>Skin element type (CQUAD4 or CMEMBRAN), thickness, and
                material</li>
            <li>Stringer count, element type (CBEAM/CBAR/CROD), section
                shape (L/T/Z/J) with dimension inputs, and material</li>
            <li>Frame element type, section shape, and material</li>
            <li>Optional <b>floor structure</b> with floor beams and
                stanchions, using either snap-to-node or mesh-modify
                connection methods</li>
        </ul>
        <p>
            Section properties (A, I1, I2, J) are computed automatically
            from the cross-section dimensions. Generator configurations can
            be saved to and loaded from JSON files.
        </p>
        </div>

        <!-- ============================================================ -->
        <h2>Current Limitations</h2>
        <!-- ============================================================ -->

        <div class="section">
        <p>The following are current limitations or explicitly
            <b>out of scope</b>:</p>
        <ul>
            <li class="warn"><b>No solver.</b> Node Runner does not run
                Nastran or any other FEA solver. It produces input files
                only.</li>
            <li class="warn"><b>Post-processing is OP2-only.</b> F06, H5,
                and PCH result files are not supported. Only OP2
                is read.</li>
            <li class="warn"><b>Limited surface meshing.</b> Geometry
                curves can be meshed into line elements and nodes can be
                created at geometry points, but there is no automatic
                surface mesher, tet mesher, or mapped meshing algorithm.
                Shell elements are created individually or via the
                fuselage generator.</li>
            <li class="warn"><b>No nonlinear cards.</b> SOL 400 / SOL 600
                specific cards (MATS1, MATEP, BCTABLE, etc.) are not
                supported.</li>
            <li class="warn"><b>No thermal analysis cards.</b> CHBDYE,
                PHBDY, CONV, QBDY, and related thermal cards are not
                supported.</li>
            <li class="warn"><b>Limited dynamic analysis.</b> EIGRL is
                supported for SOL 103 and SOL 105 (buckling). EIGB,
                FREQ, TSTEP, DLOAD, RLOAD, TLOAD, and related
                transient/frequency-response cards are not
                supported.</li>
            <li class="warn"><b>Solid element mass</b> is not included in
                the Mass &amp; CG summary (volume integration not
                implemented).</li>
            <li class="warn"><b>Fuselage generator</b> produces only
                cylindrical (constant-radius) sections. Tapered,
                double-bubble, or complex fuselage shapes are not
                supported.</li>
            <li class="warn"><b>Large models</b> (100k+ elements) may
                experience slower rendering. Node Runner is optimized
                for small-to-medium models typical of detail stress
                analysis.</li>
        </ul>
        </div>

        <!-- ============================================================ -->
        <h2>Technology Stack</h2>
        <!-- ============================================================ -->

        <div class="section">
        <table>
            <tr><td><b>Language:</b></td><td>Python 3.11+</td></tr>
            <tr><td><b>GUI Framework:</b></td><td>PySide6 (Qt 6)</td></tr>
            <tr><td><b>3D Rendering:</b></td><td>PyVista / VTK</td></tr>
            <tr><td><b>Nastran I/O:</b></td>
                <td>pyNastran (BDF read/write + OP2 results reader)</td></tr>
            <tr><td><b>Numerics:</b></td><td>NumPy</td></tr>
            <tr><td><b>Excel Import:</b></td><td>openpyxl (optional, for geometry point import)</td></tr>
        </table>
        </div>

        <hr>
        <p class="subtle" style="text-align: center;">
            Licensed under the Apache License 2.0<br>
            &copy; 2025&ndash;2026 Angel Linares / Escape Velocity Ventures, LLC<br>
            <a href="https://escvelocity.space" style="color: #5dade2;">escvelocity.space</a>
        </p>
        """
        # v5.0.0 item 21: substitute the version tag from __version__
        # so the About dialog never drifts from the package version.
        text = text.replace("__NR_VERSION__", _nr_version)

        # v5.0.0 item 21: substitute the version tag from __version__
        # so the About dialog never drifts from the package version.
        text = text.replace("__NR_VERSION__", _nr_version)

        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        layout = QVBoxLayout(dialog)

        text_edit = QTextEdit(text)
        text_edit.setReadOnly(True)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(dialog.accept)

        layout.addWidget(text_edit)
        layout.addWidget(button_box)

        dialog.setMinimumSize(700, 800)
        dialog.exec()
    
    
    
    
# START: Replacement for _show_tree_context_menu in main.py
    def _show_tree_context_menu(self, position):
        item = self.tree_widget.itemAt(position)
        if not item:
            return

        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple) or len(item_data) < 2:
            return

        entity_type, entity_id = item_data
        menu = QMenu()

        if entity_type in ['property', 'material', 'load_set', 'constraint_set', 'coord']:
            # --- MODIFIED LOGIC: Prevent editing of default coordinate systems ---
            is_default_coord = (entity_type == 'coord' and entity_id in [0, 1, 2])
            
            if not is_default_coord:
                action_text = f"Edit {entity_type.replace('_', ' ').capitalize()} {entity_id}..."
                edit_action = QAction(action_text, self)
                edit_action.triggered.connect(lambda: self._handle_edit_from_tree(entity_type, entity_id))
                menu.addAction(edit_action)

        if entity_type == 'load_set':
            add_to_set_action = QAction(f"Add to Load Set {entity_id}...", self)
            add_to_set_action.triggered.connect(
                lambda: self._add_to_load_set(entity_id))
            menu.addAction(add_to_set_action)

        if entity_type in ['load_set', 'constraint_set']:
            if menu.actions():
                menu.addSeparator()
            delete_action = QAction(f"Delete Set {entity_id}...", self)
            delete_action.triggered.connect(lambda: self._handle_delete_from_tree(entity_type, entity_id))
            menu.addAction(delete_action)

        if entity_type in ['property', 'material']:
            if menu.actions():
                menu.addSeparator()
            delete_action = QAction(f"Delete {entity_type.capitalize()} {entity_id}...", self)
            delete_action.triggered.connect(lambda: self._handle_delete_from_tree(entity_type, entity_id))
            menu.addAction(delete_action)

        if entity_type == 'coord' and entity_id not in [0, 1, 2]:
            if menu.actions():
                menu.addSeparator()
            delete_action = QAction(f"Delete Coord System {entity_id}...", self)
            delete_action.triggered.connect(lambda: self._handle_delete_from_tree(entity_type, entity_id))
            menu.addAction(delete_action)


        if menu.actions():
            menu.exec(self.tree_widget.viewport().mapToGlobal(position))
# END: Replacement for _show_tree_context_menu in main.py

    def _handle_delete_from_tree(self, entity_type, entity_id):
        if not self.current_generator: return
        model = self.current_generator.model

        if entity_type == 'load_set':
            reply = QMessageBox.question(self, "Confirm Deletion", f"Are you sure you want to delete Load Set {entity_id}?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteLoadCommand(entity_id)
                self.command_manager.execute(cmd, model)
                self._create_all_load_actors()
                self._update_status(f"Deleted Load Set {entity_id}.")
                self._populate_tree(); self._update_plot_visibility()

        elif entity_type == 'constraint_set':
            reply = QMessageBox.question(self, "Confirm Deletion", f"Are you sure you want to delete Constraint Set {entity_id}?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteConstraintCommand(entity_id)
                self.command_manager.execute(cmd, model)
                self._create_all_constraint_actors()
                self._update_status(f"Deleted Constraint Set {entity_id}.")
                self._populate_tree(); self._update_plot_visibility()

        elif entity_type == 'property':
            using_eids = [eid for eid, elem in model.elements.items() if hasattr(elem, 'pid') and elem.pid == entity_id]
            msg = f"Are you sure you want to delete Property {entity_id}?"
            if using_eids:
                msg += f"\n\nWarning: {len(using_eids)} element(s) reference this property."
            reply = QMessageBox.question(self, "Confirm Deletion", msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeletePropertyCommand(entity_id)
                self.command_manager.execute(cmd, model)
                self._update_status(f"Deleted Property {entity_id}.")
                self._populate_tree()
                self._update_viewer(self.current_generator, reset_camera=False)

        elif entity_type == 'material':
            using_pids = [pid for pid, prop in model.properties.items() if hasattr(prop, 'mid') and prop.mid == entity_id]
            msg = f"Are you sure you want to delete Material {entity_id}?"
            if using_pids:
                msg += f"\n\nWarning: {len(using_pids)} property/properties reference this material."
            reply = QMessageBox.question(self, "Confirm Deletion", msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteMaterialCommand(entity_id)
                self.command_manager.execute(cmd, model)
                self._update_status(f"Deleted Material {entity_id}.")
                self._populate_tree()

        elif entity_type == 'coord':
            reply = QMessageBox.question(self, "Confirm Deletion",
                f"Are you sure you want to delete Coordinate System {entity_id}?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteCoordCommand(entity_id)
                self.command_manager.execute(cmd, model)
                self._update_status(f"Deleted Coordinate System {entity_id}.")
                self._update_viewer(self.current_generator, reset_camera=False)

    # ------------------------------------------------------------------
    # Phase 8 – Group management
    # ------------------------------------------------------------------

    def _next_gid(self):
        """Return the next available group ID number."""
        existing = {g.get("gid", 0) for g in self.groups.values()}
        gid = self._next_group_id
        while gid in existing:
            gid += 1
        self._next_group_id = gid + 1
        return gid

    def _populate_groups_list(self):
        """Rebuild the groups tree widget from self.groups.

        v3.4.0: summary line now shows props/mats counts too. Legacy
        groups (only nodes/elements) read 0 for the new fields."""
        self.groups_list.blockSignals(True)
        self.groups_list.clear()
        for name, data in sorted(self.groups.items(),
                                 key=lambda kv: kv[1].get("gid", 0)):
            gid = data.get("gid", 0)
            n_e = len(data.get("elements", []))
            n_n = len(data.get("nodes", []))
            n_p = len(data.get("properties", []) or [])
            n_m = len(data.get("materials", []) or [])
            n_c = len(data.get("coords", []) or [])
            parts = [f"{n_e} elems", f"{n_n} nodes"]
            if n_p:
                parts.append(f"{n_p} props")
            if n_m:
                parts.append(f"{n_m} mats")
            if n_c:
                parts.append(f"{n_c} coords")
            item = QTreeWidgetItem(
                self.groups_list,
                [f"{gid}: {name}  ({', '.join(parts)})"])
            item.setData(0, QtCore.Qt.UserRole, ("group_entry", name))
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(
                0, QtCore.Qt.Unchecked if name in self._hidden_groups else QtCore.Qt.Checked
            )
        self.groups_list.blockSignals(False)

    def _open_group_add_dialog(self, mode='add'):
        """v3.4.0 item 7: open the tabbed professional Add/Remove dialog
        for the currently-highlighted group."""
        from node_runner.dialogs.group_add import GroupAddDialog

        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(
                self, "No Selection",
                "Highlight a group in the list before opening Add/Remove.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        gname = item_data[1]
        if not self.current_generator:
            QMessageBox.information(self, "No Model", "Open a model first.")
            return
        model = self.current_generator.model
        all_elements = {**model.elements, **model.rigid_elements}

        # Build the context dicts the dialog uses for live previews.
        eids_by_pid: dict[int, set[int]] = {}
        eids_by_type: dict[str, set[int]] = {}
        for eid, elem in all_elements.items():
            etype = getattr(elem, 'type', None)
            if etype is not None:
                eids_by_type.setdefault(etype, set()).add(int(eid))
            pid = getattr(elem, 'pid', None)
            if pid:
                eids_by_pid.setdefault(int(pid), set()).add(int(eid))

        pids_by_mid: dict[int, set[int]] = {}
        for pid, prop in (model.properties or {}).items():
            for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                m_id = getattr(prop, attr, None)
                if m_id:
                    pids_by_mid.setdefault(int(m_id), set()).add(int(pid))
                    break

        # eids_by_mid: walk PIDs in pids_by_mid and union their elements.
        eids_by_mid: dict[int, set[int]] = {}
        for mid, pids in pids_by_mid.items():
            s: set[int] = set()
            for pid in pids:
                s |= eids_by_pid.get(pid, set())
            eids_by_mid[mid] = s

        # nids_for_eids: for the "By connected to selected element"
        # Node-tab selector.
        nids_for_eids: dict[int, set[int]] = {}
        for eid, elem in all_elements.items():
            ns: set[int] = set()
            t = elem.type
            if t == 'RBE2':
                if getattr(elem, 'gn', None):
                    ns.add(elem.gn)
                for n in (elem.Gmi or []):
                    if n:
                        ns.add(n)
            elif t == 'RBE3':
                if getattr(elem, 'refgrid', None):
                    ns.add(elem.refgrid)
                for wt in (elem.wt_cg_groups or []):
                    for n in (wt[2] or []):
                        if n:
                            ns.add(n)
            elif hasattr(elem, 'nodes'):
                ns.update(n for n in elem.nodes if n)
            nids_for_eids[int(eid)] = ns

        ctx = {
            'group_name': gname,
            'mode': mode,
            'node_ids': set(model.nodes.keys()),
            'element_ids': set(all_elements.keys()),
            'property_ids': set(model.properties.keys()),
            'material_ids': set(model.materials.keys()),
            'coord_ids': set(model.coords.keys()),
            'eids_by_pid': eids_by_pid,
            'eids_by_mid': eids_by_mid,
            'eids_by_type': eids_by_type,
            'pids_by_mid': pids_by_mid,
            'nids_for_eids': nids_for_eids,
            'groups': self.groups,
            'current_selection_nids': set(
                getattr(self, '_current_selection_nids', []) or []),
            'current_selection_eids': set(
                getattr(self, '_current_selection_eids', []) or []),
            'current_selection_pids': set(),
            'current_selection_mids': set(),
            'current_selection_cids': set(),
        }
        dlg = GroupAddDialog(ctx, self)
        if dlg.exec() != QDialog.Accepted:
            return
        result = dlg.ids_by_entity()
        grp = self.groups[gname]
        # Ensure modern shape (legacy groups missing new fields).
        for key in ('nodes', 'elements', 'properties', 'materials', 'coords'):
            grp.setdefault(key, [])
        verb = 'Added' if mode == 'add' else 'Removed'
        total = 0
        for key in ('nodes', 'elements', 'properties', 'materials', 'coords'):
            picked = result.get(key) or set()
            if not picked:
                continue
            cur = set(grp.get(key, []) or [])
            if mode == 'add':
                cur |= picked
            else:
                cur -= picked
            grp[key] = sorted(cur)
            total += len(picked)
        self._populate_groups_list()
        self._update_plot_visibility()
        self._update_status(
            f"{verb} {total} entit{'ies' if total != 1 else 'y'} "
            f"{'to' if mode == 'add' else 'from'} group '{gname}'.")

    def _group_boolean_op(self, op: str):
        """v3.4.0 item 7: Boolean ops (union / intersect / difference)
        between two existing groups. Opens a small dialog to pick A
        and B, writes the result into a new group A_op_B."""
        if len(self.groups) < 2:
            QMessageBox.information(
                self, "Need 2 groups",
                "At least two groups required for a Boolean operation.")
            return
        names = sorted(self.groups.keys(),
                       key=lambda n: self.groups[n].get('gid', 0))
        # Two-step QInputDialog: pick A then B.
        a, ok = QInputDialog.getItem(
            self, f"Boolean: {op}", "Group A:", names, 0, False)
        if not ok:
            return
        b, ok = QInputDialog.getItem(
            self, f"Boolean: {op}", "Group B:", names, 0, False)
        if not ok:
            return
        if a == b:
            QMessageBox.warning(self, "Same group",
                                "Choose two different groups.")
            return
        ga, gb = self.groups[a], self.groups[b]
        new_data = self._make_group_data(self._next_gid())
        for key in ('nodes', 'elements', 'properties', 'materials', 'coords'):
            sa = set(ga.get(key, []) or [])
            sb = set(gb.get(key, []) or [])
            if op == 'union':
                new_data[key] = sorted(sa | sb)
            elif op == 'intersect':
                new_data[key] = sorted(sa & sb)
            else:  # difference
                new_data[key] = sorted(sa - sb)
        op_symbol = {'union': '+', 'intersect': '&', 'difference': '-'}[op]
        new_name = f"{a} {op_symbol} {b}"
        # Avoid name collisions.
        base = new_name; suffix = 1
        while new_name in self.groups:
            suffix += 1
            new_name = f"{base} ({suffix})"
        self.groups[new_name] = new_data
        self._populate_groups_list()
        self._update_status(f"Created group '{new_name}'.")

    def _create_group(self):
        name, ok = QInputDialog.getText(self, "Create Group", "Group name:")
        if ok and name:
            if name in self.groups:
                QMessageBox.warning(self, "Duplicate", f"Group '{name}' already exists.")
                return
            gid = self._next_gid()
            cmd = CreateGroupCommand(self.groups, name, gid=gid)
            model = self.current_generator.model if self.current_generator else None
            self.command_manager.execute(cmd, model)
            self._populate_groups_list()
            self._update_status(f"Created group {gid}: '{name}'.")

    def _delete_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group to delete.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        cmd = DeleteGroupCommand(self.groups, name)
        model = self.current_generator.model if self.current_generator else None
        self.command_manager.execute(cmd, model)
        self._hidden_groups.discard(name)
        self._populate_groups_list()
        self._update_plot_visibility()
        self._update_status(f"Deleted group '{name}'.")

    def _add_selected_to_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group first.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        selected_eids = list(self._get_selected_element_ids())
        selected_nids = list(self._get_selected_node_ids())
        if not selected_eids and not selected_nids:
            QMessageBox.information(
                self, "Nothing Selected",
                "Select elements or nodes first (use Edit > Select Entities, or Ctrl+F)."
            )
            return
        cmd = ModifyGroupCommand(
            self.groups, name,
            add_nodes=selected_nids, add_elements=selected_eids,
        )
        model = self.current_generator.model if self.current_generator else None
        self.command_manager.execute(cmd, model)
        self._populate_groups_list()
        self._update_status(
            f"Added {len(selected_eids)} elements, {len(selected_nids)} nodes to '{name}'."
        )

    def _isolate_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group to isolate.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        self._isolate_group_by_name(item_data[1])

    def _isolate_group_by_name(self, name):
        """Isolate: show only the given group's elements."""
        self.groups_list.blockSignals(True)
        for i in range(self.groups_list.topLevelItemCount()):
            grp_item = self.groups_list.topLevelItem(i)
            grp_data = grp_item.data(0, QtCore.Qt.UserRole)
            if not isinstance(grp_data, tuple):
                continue
            grp_name = grp_data[1]
            if grp_name == name:
                grp_item.setCheckState(0, QtCore.Qt.Checked)
                self._hidden_groups.discard(grp_name)
            else:
                grp_item.setCheckState(0, QtCore.Qt.Unchecked)
                self._hidden_groups.add(grp_name)
        self.groups_list.blockSignals(False)
        self._isolate_mode = name
        self._update_plot_visibility()
        self._isolate_mode = None
        self._update_status(f"Isolated group '{name}'.")

    def _group_by_property(self):
        """v3.4.0 item 4: auto-create one group per PID. Walks both
        model.elements AND model.rigid_elements (rigids typically have
        no PID and won't end up in any PID group, but the iteration is
        symmetric for safety).

        v5.1.0 item 27: also auto-fills `materials` for each group from
        the property's MID(s). Previously only `properties` was set, so
        downstream scoping (AnalysisSet group_target -> material
        carry-along) had to recompute the lookup."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        # Remove existing auto-generated PID groups
        for name in [n for n in self.groups if n.startswith("PID ")]:
            self.groups.pop(name)
        all_elements = {**model.elements, **model.rigid_elements}
        for pid in sorted(model.properties.keys()):
            eids = [
                eid for eid, elem in all_elements.items()
                if hasattr(elem, "pid") and elem.pid == pid
            ]
            nids: set[int] = set()
            for eid in eids:
                elem = all_elements[eid]
                if hasattr(elem, "nodes"):
                    nids.update(n for n in elem.nodes if n)
            # v5.1.0 item 27: collect referenced MIDs from this property
            # (MAT1/MAT8 properties use 'mid' / 'mid1'; PCOMP uses
            # 'mids' list; rigid-element pseudo-properties use 'mid').
            mids: set[int] = set()
            prop = model.properties.get(pid)
            if prop is not None:
                for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                    m_id = getattr(prop, attr, None)
                    if m_id:
                        mids.add(int(m_id))
                ply_mids = getattr(prop, 'mids', None)
                if ply_mids:
                    for m_id in ply_mids:
                        if m_id:
                            mids.add(int(m_id))
            gid = self._next_gid()
            self.groups[f"PID {pid}"] = self._make_group_data(
                gid, nodes=list(nids), elements=eids,
                properties=[pid], materials=sorted(mids))
        self._populate_groups_list()
        self._update_status(f"Created {len(model.properties)} property groups.")

    def _group_collect_dependencies(self, name):
        """v5.1.0 item 27: walk a group's elements, auto-fill its
        ``properties`` from the elements' PIDs, ``materials`` from those
        properties' MIDs, and ``coords`` from element CIDs + node CPs.

        Used by the Groups context menu so the user can build a group
        of elements and then one-click attach everything needed to
        export it as a self-contained sub-model.
        """
        if not self.current_generator:
            return
        if name not in self.groups:
            return
        model = self.current_generator.model
        data = self.groups[name]
        eids = list(data.get("elements", []))
        nids = set(data.get("nodes", []))

        all_elements = {**model.elements, **model.rigid_elements}
        pids: set[int] = set()
        coord_ids: set[int] = set()
        for eid in eids:
            elem = all_elements.get(eid)
            if elem is None:
                continue
            pid = getattr(elem, 'pid', None)
            if pid:
                pids.add(int(pid))
            # CBAR/CBEAM/CTUBE/CROD/etc. may carry a CID via the offt /
            # cid attribute; CHEXA/CTETRA don't.
            cid = getattr(elem, 'cid', None)
            if cid:
                coord_ids.add(int(cid))
            # Auto-include the nodes the element references so the
            # group is self-consistent.
            for n in (getattr(elem, 'nodes', None) or []):
                if n:
                    nids.add(int(n))

        mids: set[int] = set()
        for pid in pids:
            prop = model.properties.get(pid)
            if prop is None:
                continue
            for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                m_id = getattr(prop, attr, None)
                if m_id:
                    mids.add(int(m_id))
            ply_mids = getattr(prop, 'mids', None)
            if ply_mids:
                for m_id in ply_mids:
                    if m_id:
                        mids.add(int(m_id))

        # CP from nodes (the reference coord for each node).
        for nid in nids:
            node = (model.nodes or {}).get(nid)
            cp = getattr(node, 'cp', None) if node is not None else None
            if cp:
                coord_ids.add(int(cp))

        before_p = len(data.get("properties", []) or [])
        before_m = len(data.get("materials", []) or [])
        before_c = len(data.get("coords", []) or [])
        data["nodes"] = sorted(nids)
        data["properties"] = sorted(pids)
        data["materials"] = sorted(mids)
        data["coords"] = sorted(coord_ids)
        added_p = len(data["properties"]) - before_p
        added_m = len(data["materials"]) - before_m
        added_c = len(data["coords"]) - before_c
        self._populate_groups_list()
        try:
            from node_runner.profiling import perf_event
            perf_event('groups', 'collect_dependencies',
                       group=name,
                       n_props=len(data["properties"]),
                       n_mats=len(data["materials"]),
                       n_coords=len(data["coords"]),
                       added_props=added_p,
                       added_mats=added_m,
                       added_coords=added_c)
        except Exception:
            pass
        self._update_status(
            f"Group '{name}': +{added_p} props, +{added_m} mats, "
            f"+{added_c} coords (now {len(data['properties'])} props, "
            f"{len(data['materials'])} mats, {len(data['coords'])} coords).")

    def _get_selected_element_ids(self):
        """Return set of currently highlighted/selected element IDs."""
        return set(self._current_selection_eids)

    def _get_selected_node_ids(self):
        """Return set of currently highlighted/selected node IDs."""
        return set(self._current_selection_nids)

    def _show_groups_context_menu(self, position):
        item = self.groups_list.itemAt(position)
        if not item:
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        menu = QMenu()
        rename_action = QAction(f"Rename '{name}'...", self)
        rename_action.triggered.connect(lambda checked=False, n=name: self._rename_group(n))
        menu.addAction(rename_action)
        edit_action = QAction(f"Edit Contents...", self)
        edit_action.triggered.connect(lambda checked=False, i=item: self._edit_group_contents(i))
        menu.addAction(edit_action)
        # v5.1.0 item 27: one-click "make this group self-contained"
        # action. Walks the group's elements and auto-fills properties
        # / materials / coords from their dependencies so the group is
        # ready to be the target of a scoped AnalysisSet export.
        collect_deps_action = QAction("Collect Dependencies", self)
        collect_deps_action.setToolTip(
            "Walk this group's elements and auto-attach the properties, "
            "materials, and coord systems they reference. Required for a "
            "self-contained scoped export.")
        collect_deps_action.triggered.connect(
            lambda checked=False, n=name: self._group_collect_dependencies(n))
        menu.addAction(collect_deps_action)
        list_action = QAction("List Contents", self)
        list_action.triggered.connect(lambda checked=False, n=name: self._list_group_contents(n))
        menu.addAction(list_action)
        export_bdf_action = QAction("Export as BDF...", self)
        export_bdf_action.triggered.connect(lambda checked=False, n=name: self._export_group_as_bdf(n))
        menu.addAction(export_bdf_action)
        export_set_action = QAction("Export as SET...", self)
        export_set_action.triggered.connect(lambda checked=False, n=name: self._export_group_as_set(n))
        menu.addAction(export_set_action)
        menu.addSeparator()
        highlight_action = QAction("Highlight", self)
        highlight_action.triggered.connect(self._highlight_group)
        menu.addAction(highlight_action)
        isolate_action = QAction(f"Isolate '{name}'", self)
        isolate_action.triggered.connect(lambda checked=False, n=name: self._isolate_group_by_name(n))
        menu.addAction(isolate_action)
        select_action = QAction("Select in Viewport", self)
        select_action.triggered.connect(lambda checked=False, n=name: self._select_group_in_viewport(n))
        menu.addAction(select_action)
        menu.addSeparator()
        delete_action = QAction(f"Delete '{name}'...", self)
        delete_action.triggered.connect(lambda checked=False, n=name: self._delete_group_by_name(n))
        menu.addAction(delete_action)
        menu.exec(self.groups_list.viewport().mapToGlobal(position))

    def _rename_group(self, old_name):
        new_name, ok = QInputDialog.getText(self, "Rename Group", "New name:", text=old_name)
        if ok and new_name and new_name != old_name:
            if new_name in self.groups:
                QMessageBox.warning(self, "Duplicate", f"Group '{new_name}' already exists.")
                return
            cmd = RenameGroupCommand(self.groups, old_name, new_name)
            model = self.current_generator.model if self.current_generator else None
            self.command_manager.execute(cmd, model)
            if old_name in self._hidden_groups:
                self._hidden_groups.discard(old_name)
                self._hidden_groups.add(new_name)
            self._populate_groups_list()
            self._update_status(f"Renamed group '{old_name}' to '{new_name}'.")

    def _delete_group_by_name(self, name):
        cmd = DeleteGroupCommand(self.groups, name)
        model = self.current_generator.model if self.current_generator else None
        self.command_manager.execute(cmd, model)
        self._hidden_groups.discard(name)
        self._populate_groups_list()
        self._update_plot_visibility()
        self._update_status(f"Deleted group '{name}'.")

    def _on_group_item_changed(self, item, column):
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        if item.checkState(0) == QtCore.Qt.Checked:
            self._hidden_groups.discard(name)
        else:
            self._hidden_groups.add(name)
        # v4.0.1: removed cross-session QSettings persistence. The user
        # asked for in-session memory only; reopening the deck should
        # show all groups. ``self._hidden_groups`` is a plain set that
        # gets cleared when the model is reloaded.
        self._update_plot_visibility()

    def _rename_selected_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group to rename.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if isinstance(item_data, tuple):
            self._rename_group(item_data[1])

    def _remove_selected_from_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group first.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        selected_eids = list(self._get_selected_element_ids())
        selected_nids = list(self._get_selected_node_ids())
        if not selected_eids and not selected_nids:
            QMessageBox.information(self, "Nothing Selected",
                                    "Select elements or nodes first (use Edit > Select Entities, or Ctrl+F).")
            return
        cmd = ModifyGroupCommand(
            self.groups, name,
            remove_nodes=selected_nids, remove_elements=selected_eids,
        )
        model = self.current_generator.model if self.current_generator else None
        self.command_manager.execute(cmd, model)
        self._populate_groups_list()
        self._update_plot_visibility()
        self._update_status(f"Removed {len(selected_eids)} elements, {len(selected_nids)} nodes from '{name}'.")

    def _clear_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group to clear.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        group_data = self.groups.get(name, {})
        all_nodes = list(group_data.get("nodes", []))
        all_elements = list(group_data.get("elements", []))
        if not all_nodes and not all_elements:
            QMessageBox.information(self, "Empty", f"Group '{name}' is already empty.")
            return
        cmd = ModifyGroupCommand(
            self.groups, name,
            remove_nodes=all_nodes, remove_elements=all_elements,
        )
        model = self.current_generator.model if self.current_generator else None
        self.command_manager.execute(cmd, model)
        self._populate_groups_list()
        self._update_plot_visibility()
        self._update_status(f"Cleared group '{name}'.")

    def _show_all_groups(self):
        self._isolate_mode = None
        self._hidden_groups.clear()
        self.groups_list.blockSignals(True)
        for i in range(self.groups_list.topLevelItemCount()):
            self.groups_list.topLevelItem(i).setCheckState(0, QtCore.Qt.Checked)
        self.groups_list.blockSignals(False)
        self._update_plot_visibility()
        self._update_status("Showing all groups.")

    def _highlight_group(self):
        item = self.groups_list.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection", "Select a group to highlight.")
            return
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        # Toggle: if already highlighting this group, clear it
        if getattr(self, '_highlighted_group', None) == name:
            self._clear_highlight()
            return
        group_data = self.groups.get(name, {})
        nids = group_data.get("nodes", [])
        eids = group_data.get("elements", [])
        if eids:
            self._highlight_entities('Element', eids)
        elif nids:
            self._highlight_entities('Node', nids)
        self._highlighted_group = name
        self._update_status(f"Highlighted group '{name}' ({len(eids)} elements, {len(nids)} nodes).")

    def _clear_highlight(self):
        """Remove all selection highlight actors and reset state."""
        self.plotter.remove_actor('selection_highlight', render=False)
        self._highlighted_group = None
        self.plotter.render()
        self._update_status("Highlight cleared.")

    def _edit_group_contents(self, item, column=0):
        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple):
            return
        name = item_data[1]
        group_data = self.groups.get(name, {})
        model = self.current_generator.model if self.current_generator else None
        from node_runner.dialogs.info import GroupEditorDialog
        dialog = GroupEditorDialog(name, group_data, model, self)
        if dialog.exec():
            added_n, added_e, removed_n, removed_e = dialog.get_changes()
            if added_n or added_e or removed_n or removed_e:
                cmd = ModifyGroupCommand(
                    self.groups, name,
                    add_nodes=added_n, add_elements=added_e,
                    remove_nodes=removed_n, remove_elements=removed_e,
                )
                self.command_manager.execute(cmd, model)
                self._populate_groups_list()
                self._update_plot_visibility()
                self._update_status(f"Updated group '{name}'.")

    def _auto_group(self):
        if not self.current_generator:
            return
        mode = self.auto_group_combo.currentText()
        if mode == "By Property":
            self._group_by_property()
        elif mode == "By Material":
            self._group_by_material()
        elif mode == "By Element Type":
            self._group_by_element_type()

    def _group_by_material(self):
        """v3.4.0 (item 4): auto-create one group per material, with
        every element whose property's material matches. Walks both
        model.elements AND model.rigid_elements so rigid bodies are
        considered (though they typically have no PID and won't end up
        in a MAT group)."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        for name in [n for n in self.groups if n.startswith("MAT ")]:
            self.groups.pop(name)
        all_elements = {**model.elements, **model.rigid_elements}
        for mid in sorted(model.materials.keys()):
            pids_for_mid = set()
            for pid, prop in model.properties.items():
                for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                    if getattr(prop, attr, None) == mid:
                        pids_for_mid.add(pid)
                        break
            eids = [eid for eid, elem in all_elements.items()
                    if hasattr(elem, "pid") and elem.pid in pids_for_mid]
            nids: set[int] = set()
            for eid in eids:
                elem = all_elements[eid]
                if hasattr(elem, "nodes"):
                    nids.update(n for n in elem.nodes if n)
            gid = self._next_gid()
            self.groups[f"MAT {mid}"] = self._make_group_data(
                gid, nodes=list(nids), elements=eids, materials=[mid])
        self._populate_groups_list()
        self._update_status(f"Created {len(model.materials)} material groups.")

    def _group_by_element_type(self):
        """v3.4.0 (item 4): auto-create one group per element type.
        Walks both model.elements AND model.rigid_elements - previously
        RBE2/RBE3 were ignored because the loop only saw
        model.elements, so 'TYPE RBE2' / 'TYPE RBE3' groups never
        appeared."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        for name in [n for n in self.groups if n.startswith("TYPE ")]:
            self.groups.pop(name)
        type_groups: dict[str, dict] = {}
        all_elements = {**model.elements, **model.rigid_elements}
        for eid, elem in all_elements.items():
            etype = elem.type
            if etype not in type_groups:
                type_groups[etype] = {"nodes": set(), "elements": []}
            type_groups[etype]["elements"].append(eid)
            # RBE2: nodes live on .gn + .Gmi. RBE3: .refgrid + flatten
            # of wt_cg_groups[*][2]. Everything else exposes .nodes.
            if etype == 'RBE2':
                if getattr(elem, 'gn', None):
                    type_groups[etype]["nodes"].add(elem.gn)
                for g in (getattr(elem, 'Gmi', None) or []):
                    if g:
                        type_groups[etype]["nodes"].add(g)
            elif etype == 'RBE3':
                if getattr(elem, 'refgrid', None):
                    type_groups[etype]["nodes"].add(elem.refgrid)
                for wt in (getattr(elem, 'wt_cg_groups', None) or []):
                    for g in (wt[2] or []):
                        if g:
                            type_groups[etype]["nodes"].add(g)
            elif hasattr(elem, "nodes"):
                type_groups[etype]["nodes"].update(
                    n for n in elem.nodes if n)
        for etype, data in sorted(type_groups.items()):
            gid = self._next_gid()
            self.groups[f"TYPE {etype}"] = self._make_group_data(
                gid, nodes=list(data["nodes"]),
                elements=data["elements"])
        self._populate_groups_list()
        self._update_status(f"Created {len(type_groups)} element type groups.")

    @staticmethod
    def _make_group_data(gid, nodes=None, elements=None,
                         properties=None, materials=None, coords=None):
        """v3.4.0 item 7: factory for the extended group dict shape.

        Pre-v3.4.0 groups were ``{gid, nodes, elements}``. We now also
        track properties / materials / coords so a group can express
        professional multi-entity membership. Older code paths that read
        only nodes/elements still work because the new fields default
        to empty lists."""
        return {
            "gid": int(gid),
            "nodes": list(nodes or []),
            "elements": list(elements or []),
            "properties": list(properties or []),
            "materials": list(materials or []),
            "coords": list(coords or []),
        }

    def _select_group_in_viewport(self, name):
        group_data = self.groups.get(name, {})
        eids = group_data.get("elements", [])
        nids = group_data.get("nodes", [])
        self._current_selection_eids = list(eids)
        self._current_selection_nids = list(nids)
        # v5.2.0 item 45: keep Data Table live-sync working.
        try:
            self.selection_changed.emit('Element', list(eids))
            self.selection_changed.emit('Node', list(nids))
        except Exception:
            pass
        if eids:
            self._highlight_entities('Element', eids)
        elif nids:
            self._highlight_entities('Node', nids)
        self._highlighted_group = name
        self._update_status(f"Selected group '{name}' ({len(eids)} elements, {len(nids)} nodes).")

    def _list_group_contents(self, name):
        """Show group contents in a tabular info dialog."""
        group_data = self.groups.get(name, {})
        model = self.current_generator.model if self.current_generator else None
        if not model:
            return
        from node_runner.dialogs.info import NodeElementInfoDialog

        nids = group_data.get("nodes", [])
        eids = group_data.get("elements", [])
        all_elements = {**model.elements, **model.rigid_elements}

        # Build node rows
        node_rows = []
        node_columns = ['NID', 'X', 'Y', 'Z', 'CP', 'CD']
        node_tooltips = {
            'NID': 'Node ID', 'X': 'X coordinate', 'Y': 'Y coordinate', 'Z': 'Z coordinate',
            'CP': 'Position Coordinate System (0 = basic rectangular)',
            'CD': 'Displacement Coordinate System (0 = basic rectangular)',
        }
        for nid in sorted(nids):
            if nid in model.nodes:
                node = model.nodes[nid]
                node_rows.append({
                    'NID': node.nid,
                    'X': round(node.xyz[0], 6), 'Y': round(node.xyz[1], 6), 'Z': round(node.xyz[2], 6),
                    'CP': node.cp, 'CD': node.cd,
                })

        # Build element rows
        elem_rows = []
        elem_columns = ['EID', 'Type', 'PID', 'MID', 'Nodes']
        elem_tooltips = {
            'EID': 'Element ID', 'Type': 'Element type',
            'PID': 'Property ID', 'MID': 'Material ID',
            'Nodes': 'Grid point connectivity',
        }
        for eid in sorted(eids):
            if eid in all_elements:
                elem = all_elements[eid]
                pid = getattr(elem, 'pid', '')
                mid = ''
                if pid and pid in model.properties:
                    prop = model.properties[pid]
                    mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None) or ''
                nodes_str = ', '.join(str(n) for n in elem.nodes) if hasattr(elem, 'nodes') else ''
                elem_rows.append({
                    'EID': elem.eid, 'Type': elem.type,
                    'PID': pid, 'MID': mid, 'Nodes': nodes_str,
                })

        # Show nodes if available, otherwise elements
        if node_rows:
            title = f"Group '{name}' - Nodes ({len(node_rows)})"
            dlg = NodeElementInfoDialog(title, node_rows, node_columns, self, column_tooltips=node_tooltips)
            dlg.exec()
        if elem_rows:
            title = f"Group '{name}' - Elements ({len(elem_rows)})"
            dlg = NodeElementInfoDialog(title, elem_rows, elem_columns, self, column_tooltips=elem_tooltips)
            dlg.exec()
        if not node_rows and not elem_rows:
            QMessageBox.information(self, "Empty Group", f"Group '{name}' has no contents.")

    def _export_group_as_bdf(self, name):
        """Export the group's nodes, elements, properties, and materials as a standalone BDF."""
        group_data = self.groups.get(name, {})
        model = self.current_generator.model if self.current_generator else None
        if not model:
            return
        nids = group_data.get("nodes", [])
        eids = group_data.get("elements", [])
        if not nids and not eids:
            QMessageBox.information(self, "Empty Group", f"Group '{name}' has no contents to export.")
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self, f"Export Group '{name}' as BDF", f"{name}.bdf",
            "Nastran Files (*.bdf *.dat);;All Files (*)")
        if not filepath:
            return

        try:
            new_model = BDF(debug=False)
            all_elements = {**model.elements, **model.rigid_elements}

            # Copy nodes
            for nid in nids:
                if nid in model.nodes:
                    new_model.nodes[nid] = model.nodes[nid]

            # Copy elements and collect PIDs
            pids_needed = set()
            for eid in eids:
                if eid in all_elements:
                    elem = all_elements[eid]
                    if eid in model.elements:
                        new_model.elements[eid] = elem
                    else:
                        new_model.rigid_elements[eid] = elem
                    pid = getattr(elem, 'pid', None)
                    if pid and isinstance(pid, int):
                        pids_needed.add(pid)

            # Copy properties and collect MIDs
            mids_needed = set()
            for pid in pids_needed:
                if pid in model.properties:
                    prop = model.properties[pid]
                    new_model.properties[pid] = prop
                    for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                        val = getattr(prop, attr, None)
                        if val and isinstance(val, int):
                            mids_needed.add(val)

            # Copy materials
            for mid in mids_needed:
                if mid in model.materials:
                    new_model.materials[mid] = model.materials[mid]

            # Copy coordinate systems referenced by nodes
            cids_needed = set()
            for nid in nids:
                if nid in model.nodes:
                    node = model.nodes[nid]
                    if node.cp and node.cp > 0:
                        cids_needed.add(node.cp)
                    if node.cd and node.cd > 0:
                        cids_needed.add(node.cd)
            for cid in cids_needed:
                if cid in model.coords:
                    new_model.coords[cid] = model.coords[cid]

            # Respect the user's export-format preference (Phase 1).
            fmt = self._export_format if self._export_format in ("short", "long", "free") else "short"
            if fmt == "long":
                new_model.write_bdf(filepath, size=16, is_double=False)
            else:
                new_model.write_bdf(filepath, size=8, is_double=False)
                if fmt == "free":
                    from node_runner.model import _convert_bdf_to_free
                    _convert_bdf_to_free(filepath)
            self._update_status(f"Exported group '{name}' as BDF ({fmt}): {os.path.basename(filepath)}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export group as BDF:\n{e}")

    def _export_group_as_set(self, name):
        """Export the group's node/element IDs as Nastran SET1 cards."""
        group_data = self.groups.get(name, {})
        nids = sorted(group_data.get("nodes", []))
        eids = sorted(group_data.get("elements", []))
        if not nids and not eids:
            QMessageBox.information(self, "Empty Group", f"Group '{name}' has no contents to export.")
            return

        filepath, _ = QFileDialog.getSaveFileName(
            self, f"Export Group '{name}' as SET", f"{name}_sets.bdf",
            "Nastran Files (*.bdf *.dat);;Text Files (*.txt);;All Files (*)")
        if not filepath:
            return

        def compress_with_thru(ids):
            """Compress consecutive ID ranges using THRU notation."""
            if not ids:
                return []
            result = []
            start = ids[0]
            prev = ids[0]
            for val in ids[1:]:
                if val == prev + 1:
                    prev = val
                else:
                    if prev - start >= 2:
                        result.extend([start, 'THRU', prev])
                    elif prev - start == 1:
                        result.extend([start, prev])
                    else:
                        result.append(start)
                    start = val
                    prev = val
            # Flush last range
            if prev - start >= 2:
                result.extend([start, 'THRU', prev])
            elif prev - start == 1:
                result.extend([start, prev])
            else:
                result.append(start)
            return result

        def write_set1_card(f, sid, ids, comment):
            """Write a SET1 card in small-field format (8-char fields)."""
            compressed = compress_with_thru(ids)
            f.write(f"$ {comment}\n")
            # First line: SET1 + SID + up to 6 values
            fields_per_first_line = 6
            fields_per_cont_line = 7
            line = f"{'SET1':8s}{sid:>8d}"
            count = 0
            for item in compressed:
                if count < fields_per_first_line and count == 0 or count < fields_per_first_line:
                    line += f"{str(item):>8s}"
                    count += 1
                else:
                    f.write(line + "\n")
                    line = f"{'':8s}{str(item):>8s}"
                    count = 1
            if line.strip():
                f.write(line + "\n")

        try:
            with open(filepath, 'w') as f:
                f.write(f"$ SET cards exported from group: {name}\n")
                f.write(f"$ Generated by Node Runner\n")
                sid = 1
                if nids:
                    write_set1_card(f, sid, nids, f"Nodes for group: {name} ({len(nids)} nodes)")
                    sid += 1
                if eids:
                    write_set1_card(f, sid, eids, f"Elements for group: {name} ({len(eids)} elements)")
            self._update_status(f"Exported group '{name}' as SET: {os.path.basename(filepath)}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export group as SET:\n{e}")

    def _handle_edit_from_tree(self, entity_type, entity_id):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator: return
        
        dialog = None
        if entity_type == 'load_set':
            # Open the Load Set Manager for editing existing sets
            self._open_load_set_manager(entity_id)
            return
        elif entity_type == 'constraint_set':
            dialog = CreateConstraintDialog(self.current_generator.model, self, existing_sid=entity_id)
            dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
            dialog.accepted.connect(self._on_constraint_creation_accept)
            dialog.rejected.connect(self._on_creation_reject)
        elif entity_type == 'coord':
            coord_to_edit = self.current_generator.model.coords.get(entity_id)
            if coord_to_edit:
                dialog = CreateCoordDialog(entity_id, self.current_generator.model, self, existing_coord=coord_to_edit)
                dialog.accepted.connect(self._on_coord_creation_accept)
                dialog.rejected.connect(self._on_creation_reject)
        elif entity_type == 'property':
            prop_to_edit = self.current_generator.model.properties.get(entity_id)
            if prop_to_edit:
                edit_dialog = CreatePropertyDialog(entity_id, self.current_generator.model, self, prop_to_edit)
                if edit_dialog.exec():
                    if params := edit_dialog.get_parameters():
                        self._edit_property(entity_id, params)
                        self._populate_tree()
                        self._update_viewer(self.current_generator, reset_camera=False)
            return
        elif entity_type == 'material':
            mat_to_edit = self.current_generator.model.materials.get(entity_id)
            if mat_to_edit:
                edit_dialog = CreateMaterialDialog(entity_id, self, mat_to_edit)
                if edit_dialog.exec():
                    if params := edit_dialog.get_parameters():
                        self._edit_material(entity_id, params)
                        self._populate_tree()
            return

        if dialog:
            dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            dialog.show()
            self.active_creation_dialog = dialog

    def _create_material(self):
        if not self.current_generator:
            self.current_generator = NastranModelGenerator()
            self._ensure_default_coords(self.current_generator.model)

        next_mid = max(self.current_generator.model.materials.keys()) + 1 if self.current_generator.model.materials else 1
        
        dialog = CreateMaterialDialog(next_mid, self)
        dialog.import_requested.connect(self._import_materials_from_bdf)
        if dialog.exec():
            if params := dialog.get_parameters():
                mat_type = params['type']
                cmd = AddMaterialCommand(params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created Material MID {params['mid']} ({mat_type}).")
                self._populate_tree()

    @staticmethod
    def _mat_object_to_params(mat):
        """Convert a pyNastran material object to the params dict used by AddMaterialCommand."""
        comment = mat.comment.strip() if mat.comment else ''
        params = {'mid': mat.mid, 'comment': comment}
        if mat.type == 'MAT1':
            params.update({
                'type': 'MAT1',
                'E': mat.e, 'G': mat.g, 'nu': mat.nu,
                'rho': mat.rho, 'a': mat.a, 'tref': mat.tref, 'ge': mat.ge,
            })
        elif mat.type == 'MAT8':
            params.update({
                'type': 'MAT8',
                'E1': mat.e1, 'E2': mat.e2, 'nu12': mat.nu12,
                'G12': mat.g12, 'G1z': mat.g1z, 'G2z': mat.g2z,
                'rho': mat.rho, 'a1': mat.a1, 'a2': mat.a2, 'tref': mat.tref,
            })
        elif mat.type == 'MAT9':
            params['type'] = 'MAT9'
            for label in [
                'G11', 'G12', 'G13', 'G14', 'G15', 'G16',
                'G22', 'G23', 'G24', 'G25', 'G26',
                'G33', 'G34', 'G35', 'G36',
                'G44', 'G45', 'G46',
                'G55', 'G56', 'G66',
            ]:
                params[label] = getattr(mat, label.lower(), 0.0)
            params.update({'rho': mat.rho, 'tref': mat.tref})
        else:
            return None  # Unsupported material type
        return params

    def _import_materials_from_bdf(self):
        """Import materials from an external BDF file into the current model."""
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Import Materials from BDF", "",
            "Nastran Files (*.bdf *.dat *.nas);;All Files (*)")
        if not filepath:
            return

        if not self.current_generator:
            self.current_generator = NastranModelGenerator()
            self._ensure_default_coords(self.current_generator.model)

        # Parse the BDF file using the existing lenient pipeline
        try:
            temp_model, lenient_result = NastranModelGenerator._read_bdf_robust(filepath)
        except RuntimeError as e:
            QMessageBox.critical(self, "Import Failed", str(e))
            return

        if not temp_model.materials:
            QMessageBox.information(self, "No Materials Found",
                                   "The selected file does not contain any material cards.")
            return

        existing_mids = set(self.current_generator.model.materials.keys())
        imported = []       # list of params dicts
        skipped_dup = []    # list of (mid, mat_type) for MID conflicts
        skipped_type = []   # list of (mid, mat_type_str) for unsupported types

        for mid, mat in temp_model.materials.items():
            if mid in existing_mids:
                skipped_dup.append((mid, mat.type))
                continue
            params = self._mat_object_to_params(mat)
            if params is None:
                skipped_type.append((mid, mat.type))
                continue
            imported.append(params)

        if not imported:
            QMessageBox.information(self, "Nothing to Import",
                                   "No new materials to import (all MIDs already exist or unsupported types).")
            return

        # Build compound command for single undo step
        import os
        fname = os.path.basename(filepath)
        cmds = [AddMaterialCommand(p) for p in imported]
        compound = CompoundCommand(cmds, f"Import {len(cmds)} materials from {fname}")
        self.command_manager.execute(compound, self.current_generator.model)

        self._populate_tree()
        self._rebuild_plot()

        # Build and show report
        from node_runner.dialogs.info import MaterialImportReportDialog
        parse_failures = []
        if lenient_result and lenient_result.skipped:
            for sc in lenient_result.skipped:
                if sc.card.startswith('MAT'):
                    parse_failures.append(sc)

        report = MaterialImportReportDialog(
            fname, imported, skipped_dup, skipped_type,
            parse_failures, parent=self)
        report.exec()

        total = len(imported)
        self._update_status(f"Imported {total} material{'s' if total != 1 else ''} from {fname}.")

    def _create_property(self):
        if not self.current_generator:
            self.current_generator = NastranModelGenerator()
            self._ensure_default_coords(self.current_generator.model)

        if not self.current_generator.model.materials:
            QMessageBox.warning(self, "No Materials", "You must create a material before creating a property.")
            return

        next_pid = max(self.current_generator.model.properties.keys()) + 1 if self.current_generator.model.properties else 1
        
        dialog = CreatePropertyDialog(next_pid, self.current_generator.model, self)
        if dialog.exec():
            if params := dialog.get_parameters():
                cmd = AddPropertyCommand(params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created Property PID {params['pid']} ({params['type']}).")
                self._populate_tree()
                self._update_viewer(self.current_generator, reset_camera=False)

    def _create_coord(self):
        QMessageBox.information(self, "Not Implemented", "Coordinate system creation will be added in a future update.")

    def _add_property_card_from_params(self, params):
        params_copy = params.copy()
        prop_type = params_copy.pop('type')
        model = self.current_generator.model
        pid = params_copy.pop('pid')
        comment = params_copy.pop('comment')

        if prop_type == 'PSHELL':
            model.add_pshell(pid, mid1=params_copy['mid1'], t=params_copy['t'], nsm=params_copy['nsm'], comment=comment)
        elif prop_type == 'PCOMP':
            plies = params_copy['plies']
            if plies:
                mids, thicknesses, thetas, souts = zip(*plies)
                model.add_pcomp(pid, list(mids), list(thicknesses), list(thetas), souts=list(souts),
                                nsm=params_copy['nsm'], ft=params_copy['ft'], comment=comment)
            else:
                model.add_pcomp(pid, [], [], [], nsm=params_copy['nsm'], ft=params_copy['ft'], comment=comment)

        elif prop_type == 'PBAR':
            model.add_pbar(pid, params_copy['mid'], params_copy['A'], params_copy['i1'], params_copy['i2'], params_copy['j'], comment=comment)
        elif prop_type == 'PBEAM':
            model.add_pbeam(pid, params_copy['mid'], [0.0], ['C'], [params_copy['A']], [params_copy['i1']],
                            [params_copy['i2']], [0.0], [params_copy['j']], comment=comment)
        elif prop_type == 'PBUSH':
            model.add_pbush(pid, k=params_copy['k'], b=[], ge=[], comment=comment)

    def _set_standard_view(self, view):
        if view == 'top': self.plotter.view_xy()
        elif view == 'bottom': self.plotter.view_xy(negative=True)
        elif view == 'front': self.plotter.view_yz()
        elif view == 'back': self.plotter.view_yz(negative=True)
        elif view == 'right': self.plotter.view_xz()
        elif view == 'left': self.plotter.view_xz(negative=True)
        self._update_status(f"View set to: {view.capitalize()}")

    # ------------------------------------------------------------------
    # View-cube icon generation (CATIA / SolidWorks style)
    # ------------------------------------------------------------------

    @staticmethod
    def _make_view_cube_icon(view_name, size=24):
        """Generate a small isometric 3D-cube icon with one face highlighted.

        Blue highlight = positive direction (top/front/right).
        Orange highlight = negative direction (bottom/back/left).
        Isometric = all faces tinted equally.
        """
        from PySide6.QtCore import QPointF

        pixmap = QPixmap(size, size)
        pixmap.fill(QColor(0, 0, 0, 0))

        p = QPainter(pixmap)
        p.setRenderHint(QPainter.Antialiasing)

        cx, cy = size * 0.5, size * 0.5
        w = size * 0.34       # horizontal half-extent of the hexagon
        h = size * 0.19       # vertical step per axis

        # Six vertices of the hexagonal isometric-cube projection
        vtop = QPointF(cx, cy - 2 * h)
        vul  = QPointF(cx - w, cy - h)
        vur  = QPointF(cx + w, cy - h)
        vctr = QPointF(cx, cy)
        vll  = QPointF(cx - w, cy + h)
        vlr  = QPointF(cx + w, cy + h)
        vbot = QPointF(cx, cy + 2 * h)

        # 3 visible cube faces (back-top-right corner is at vctr).
        top_face   = QPolygonF([vtop, vur, vctr, vul])
        left_face  = QPolygonF([vul, vctr, vbot, vll])    # Y+ face → "right" view
        right_face = QPolygonF([vur, vlr, vbot, vctr])    # X+ face → "front" view

        # v4.0.1: 3 hidden cube faces. opposite each visible face, all
        # meeting at the front-bottom-left cube corner (which projects
        # to the same screen point as the back-top-right: vctr). These
        # are rendered with a dashed outline + lower alpha so the user
        # sees the negative-direction face "through" the cube.
        bottom_face_h = QPolygonF([vctr, vlr, vbot, vll])  # Z- → "bottom"
        back_face_h   = QPolygonF([vctr, vll, vul, vtop])  # X- → "back"
        leftv_face_h  = QPolygonF([vctr, vlr, vur, vtop])  # Y- → "left" view

        is_neg = view_name in ('bottom', 'back', 'left')
        hidden_for_neg = {
            'bottom': bottom_face_h,
            'back':   back_face_h,
            'left':   leftv_face_h,
        }

        # Map positive-direction views to which visible face highlights.
        face_key_pos = {
            'top': 'top',
            'front': 'right',
            'right': 'left',
            'iso': 'all',
        }
        active = face_key_pos.get(view_name, None)
        # Negative views: no visible face is highlighted; the hidden
        # face does the work.
        if is_neg:
            active = None

        # Highlight colour
        hi_blue   = QColor(100, 165, 255, 215)
        hi_orange = QColor(248, 163, 76, 230)

        # Dim face shades (visible faces)
        dim_top   = QColor(78, 83, 112, 175)
        dim_right = QColor(62, 67, 97, 165)
        dim_left  = QColor(50, 55, 82, 155)

        # Isometric: pleasant blue tints
        iso_top   = QColor(92, 152, 215, 195)
        iso_right = QColor(72, 122, 190, 180)
        iso_left  = QColor(58, 102, 170, 170)

        if active == 'all':
            colors = {'top': iso_top, 'right': iso_right, 'left': iso_left}
        elif active is None:
            # Negative view: visible faces are dim + reduced alpha so the
            # hidden orange face shows through.
            colors = {
                'top':   QColor(dim_top.red(),   dim_top.green(),   dim_top.blue(),   90),
                'right': QColor(dim_right.red(), dim_right.green(), dim_right.blue(), 85),
                'left':  QColor(dim_left.red(),  dim_left.green(),  dim_left.blue(),  80),
            }
        else:
            colors = {
                'top':   hi_blue if active == 'top'   else dim_top,
                'right': hi_blue if active == 'right' else dim_right,
                'left':  hi_blue if active == 'left'  else dim_left,
            }

        solid_pen = QPen(QColor(170, 190, 230, 195), 1.0)
        # v4.0.1: dashed/hidden-line pen for the negative-direction face.
        dashed_pen = QPen(QColor(255, 200, 130, 240), 1.2, QtCore.Qt.DashLine)

        # v4.0.1: for negative views, paint the hidden face FIRST so
        # the visible faces (drawn next with reduced alpha) sit on top.
        # The hidden face's edges then poke through as a dashed outline.
        if is_neg:
            hidden_poly = hidden_for_neg[view_name]
            p.setPen(QtCore.Qt.NoPen)
            p.setBrush(QBrush(hi_orange))
            p.drawPolygon(hidden_poly)

        # Draw the 3 visible cube faces back-to-front: left, right, top
        for fname, poly in [('left', left_face),
                            ('right', right_face),
                            ('top', top_face)]:
            p.setPen(solid_pen)
            p.setBrush(QBrush(colors[fname]))
            p.drawPolygon(poly)

        # v4.0.1: overlay the hidden-face outline as dashed lines so
        # the orange face reads as a "hidden line" peek-through.
        if is_neg:
            hidden_poly = hidden_for_neg[view_name]
            p.setPen(dashed_pen)
            p.setBrush(QtCore.Qt.NoBrush)
            p.drawPolygon(hidden_poly)

        p.end()
        return QIcon(pixmap)

    def _toggle_origin(self, state):
        if state: self.plotter.add_mesh(pv.Sphere(center=(0,0,0), radius=0.01 * self.current_grid.length if self.current_grid else 0.05), color='cyan', name='origin_actor', reset_camera=False)
        else: self.plotter.remove_actor('origin_actor')

    def _toggle_shading(self, state):
        """v3.4.0 item 8: flip between phong (lit) and a flat unlit
        look that's obviously different to the eye.

        Before v3.4.0 the toggle changed ambient/diffuse but kept the
        same interpolation, so on a coloured-by-property deck with
        edges drawn the visible difference was tiny. Now:
          - ON  : lit phong, low ambient + high diffuse (3D shaded)
          - OFF : unlit flat, full ambient, no diffuse (flat colours,
                  no gradients across faces)

        v4.0.0 (A2): the per-actor state-set is split into
        ``_apply_shading_to_actor`` so that actors created later by
        ``add_mesh`` can inherit the current state without waiting for
        the next user toggle.
        """
        from node_runner.profiling import perf_event, perf_stage, enabled as _profile_on
        import os
        debug = os.environ.get('NR_DEBUG_SHADING') == '1' or _profile_on()

        # v4.0.1 (Stage 4 diagnostic): inspect BOTH actor containers.
        # plotter.renderer.actors is the legacy ctor-time iteration;
        # plotter.actors is what the rest of the codebase reads (see
        # `plotter.actors.get('rbe_actors')` for the working precedent).
        renderer_actors = list(self.plotter.renderer.actors.values())
        plotter_actors_dict = dict(self.plotter.actors)

        if debug:
            perf_event('shading', 'container_counts',
                       toggle_state=state,
                       renderer_actors=len(renderer_actors),
                       plotter_actors=len(plotter_actors_dict),
                       plotter_keys=','.join(sorted(plotter_actors_dict.keys())))
            # v4.0.8: light_count catcher. plotter.clear() in
            # _update_viewer removes lights as well as actors; if the
            # rig is missing, lighting=True/False render identically
            # and the shading toggle has no visible effect. n_lights=0
            # here is the smoking gun for that case.
            try:
                _lights = list(self.plotter.renderer.lights)
                perf_event('shading', 'light_count',
                           n_lights=len(_lights),
                           n_lights_switched_on=sum(
                               1 for l in _lights
                               if hasattr(l, 'GetSwitch') and l.GetSwitch()))
            except Exception as _exc:
                perf_event('shading', 'light_count_failed',
                           exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            # v4.0.3 (Stage L): per-actor `_log_shading_actor` calls
            # are now wrapped in try/except. v4.0.2 had this loop
            # raising silently on vtkActor2D HUD entries which killed
            # the rest of `_toggle_shading` (including the status bar
            # update). Now we keep going on any single-actor failure
            # and emit a perf_event so we know which actor blew up.
            for name, actor in plotter_actors_dict.items():
                # Skip HUD / 2D actors that don't participate in
                # shading (orientation marker, scalar bar, axes text).
                if isinstance(name, str) and (
                        name.startswith('vtkActor2D')
                        or name.startswith('vtkScalarBarActor')
                        or name.startswith('vtkAxisActor2D')
                        or name.startswith('vtkOrientationMarker')):
                    continue
                try:
                    if hasattr(actor, 'prop') and actor.prop is not None:
                        self._log_shading_actor("BEFORE", actor, actor_name=name)
                except Exception as _exc:
                    perf_event('shading', 'log_actor_failed',
                               name=name,
                               phase='BEFORE',
                               exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            print(f"[shading] toggle -> {state!r}, "
                  f"renderer_actors={len(renderer_actors)}, "
                  f"plotter_actors={len(plotter_actors_dict)}", flush=True)

        # v4.0.1: iterate the plotter.actors dict (which contains the
        # named mesh actors), not just renderer.actors. The union
        # captures both code paths so we don't miss any actor regardless
        # of which container PyVista actually populated.
        seen = set()
        applied = 0
        failed = 0
        with perf_stage('shading', 'apply_to_all_actors',
                        n_renderer=len(renderer_actors),
                        n_plotter=len(plotter_actors_dict)):
            for actor in list(plotter_actors_dict.values()) + renderer_actors:
                key = id(actor)
                if key in seen:
                    continue
                seen.add(key)
                # v4.0.3 (Stage L): never let a single actor's prop
                # accessor kill the whole toggle. vtkActor2D and other
                # HUD actors don't have .prop as a PyVista wrapper;
                # the access raises and used to propagate out before
                # the status-bar update.
                try:
                    self._apply_shading_to_actor(actor, enabled=state)
                    applied += 1
                except Exception as _exc:
                    failed += 1
                    perf_event('shading', 'apply_failed',
                               actor_class=type(actor).__name__,
                               exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
        perf_event('shading', 'apply_summary',
                   applied=applied, failed=failed, total_seen=len(seen))

        if debug and plotter_actors_dict:
            for name, actor in plotter_actors_dict.items():
                if isinstance(name, str) and (
                        name.startswith('vtkActor2D')
                        or name.startswith('vtkScalarBarActor')
                        or name.startswith('vtkAxisActor2D')
                        or name.startswith('vtkOrientationMarker')):
                    continue
                try:
                    if hasattr(actor, 'prop') and actor.prop is not None:
                        self._log_shading_actor("AFTER ", actor, actor_name=name)
                except Exception as _exc:
                    perf_event('shading', 'log_actor_failed',
                               name=name,
                               phase='AFTER',
                               exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")

        # Keep both UI surfaces in sync (View menu + Display tab).
        try:
            self.shading_action.blockSignals(True)
            self.shading_action.setChecked(state)
        finally:
            self.shading_action.blockSignals(False)
        try:
            self.shading_display_check.blockSignals(True)
            self.shading_display_check.setChecked(state)
        finally:
            self.shading_display_check.blockSignals(False)

        self._update_status(f"Shading: {'phong (lit)' if state else 'flat (unlit)'}.")
        self.plotter.render()

    def _apply_shading_to_actor(self, actor, enabled=None):
        """Apply current shading state to a single actor.

        Called by ``_toggle_shading`` (in a loop) and by ``add_mesh``
        sites after they create a new actor so new actors inherit the
        user's current shading preference. ``enabled=None`` reads the
        live state from the View-menu action.
        """
        if actor is None or not hasattr(actor, 'prop') or actor.prop is None:
            return
        if enabled is None:
            try:
                enabled = self.shading_action.isChecked()
            except AttributeError:
                # v4.0.1: default to flat/unlit. This is the user's
                # preferred default; matches the View > Shading action.
                enabled = False
        if enabled:
            # v4.0.0 (A2): bump ambient from 0.15 -> 0.35. v5.0.0 item 9
            # bumps it further to 0.40 (in concert with bumped key/fill/
            # headlight intensities) to address the "mesh appears dark"
            # user feedback. Multiplied by brightness mode (Subdued /
            # Normal / Bright) so the user can tune it.
            mult = self._brightness_multiplier()
            actor.prop.lighting = True
            actor.prop.interpolation = 'phong'
            actor.prop.ambient = 0.40 * mult
            actor.prop.diffuse = 0.65
            actor.prop.specular = 0.0
        else:
            actor.prop.lighting = False
            actor.prop.interpolation = 'flat'
            actor.prop.ambient = 1.0
            actor.prop.diffuse = 0.0
            actor.prop.specular = 0.0

    def _brightness_multiplier(self):
        """v5.0.0 item 9: brightness multiplier read from QSettings.

        Lazy-cached on ``self`` so the per-actor shading hook isn't
        hitting QSettings on every paint. Cache is invalidated on
        Preferences Apply via :meth:`_invalidate_brightness_cache`.
        """
        mult = getattr(self, '_brightness_mult_cache', None)
        if mult is not None:
            return mult
        try:
            from PySide6.QtCore import QSettings
            mode = QSettings("NodeRunner", "NodeRunner").value(
                "render/brightness_mode", "Normal", type=str)
        except Exception:
            mode = "Normal"
        mult = {"Subdued": 0.85, "Normal": 1.0, "Bright": 1.15}.get(mode, 1.0)
        self._brightness_mult_cache = mult
        return mult

    def _invalidate_brightness_cache(self):
        self._brightness_mult_cache = None

    def _install_light_rig(self):
        """v4.0.0 (A2): replace PyVista's default single headlight with
        a warm key + cool fill + dim back rim. Makes the phong vs flat
        toggle clearly visible and gives the viewport a premium look.

        v5.0.0 item 9: bumped fill 0.45 -> 0.55 and headlight 0.35 ->
        0.55 to address "mesh appears dark" feedback. Intensities are
        further scaled by the user's brightness preference.
        Failures are swallowed: rendering with the default light setup
        is still acceptable.
        """
        mult = self._brightness_multiplier()
        key_i  = 0.85 * mult
        fill_i = 0.55 * mult
        rim_i  = 0.25 * mult
        head_i = 0.55 * mult
        try:
            import pyvista as pv
            renderer = self.plotter.renderer
            renderer.remove_all_lights()

            key = pv.Light(
                position=(1.0, 1.0, 1.0),
                light_type='scene light',
                color=(1.0, 0.96, 0.90),  # warm
                intensity=key_i,
            )
            fill = pv.Light(
                position=(-1.0, -0.5, 1.0),
                light_type='scene light',
                color=(0.85, 0.90, 1.0),  # cool
                intensity=fill_i,
            )
            rim = pv.Light(
                position=(0.0, -1.0, -0.5),
                light_type='scene light',
                color=(1.0, 1.0, 1.0),
                intensity=rim_i,
            )
            head = pv.Light(light_type='headlight', intensity=head_i)

            for lt in (key, fill, rim, head):
                renderer.add_light(lt)
        except Exception:
            pass
        # v4.0.8: emit count so profile log proves the rig is present
        # at every install. plotter.clear() in _update_viewer wipes
        # lights along with actors; without this catcher the
        # "no visible shading delta" symptom looked identical to the
        # "lights installed but shader doesn't honour them" case.
        try:
            from node_runner.profiling import perf_event
            perf_event('lights', 'rig_installed',
                       n_lights=len(self.plotter.renderer.lights))
            perf_event('lights', 'intensity_summary',
                       mult=mult, key=key_i, fill=fill_i,
                       rim=rim_i, headlight=head_i,
                       total=key_i + fill_i + rim_i + head_i)
        except Exception:
            pass

    def _emit_mem_event(self, where):
        """v4.0.9: emit a memory snapshot to the profile log. Zero-dep
        on Windows via ctypes + GetProcessMemoryInfo. Failure is
        swallowed so this never blocks the calling path. The
        ``where`` field tags the snapshot location (e.g.
        ``viewer.total``, ``visibility.cycle``).

        Establishes the "Lightweight" pillar baseline before
        pre-built per-category actors permanently retain mesh
        polydatas in memory. v4.1.0 acts on the data this emits.
        """
        try:
            from node_runner.profiling import perf_event
            rss_mb = None
            try:
                import psutil  # pragma: no cover - optional dep
                rss_mb = psutil.Process().memory_info().rss / (1024 ** 2)
            except Exception:
                pass
            if rss_mb is None:
                try:
                    import ctypes
                    from ctypes import wintypes

                    class _PMC(ctypes.Structure):
                        _fields_ = [
                            ('cb', wintypes.DWORD),
                            ('PageFaultCount', wintypes.DWORD),
                            ('PeakWorkingSetSize', ctypes.c_size_t),
                            ('WorkingSetSize', ctypes.c_size_t),
                            ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                            ('PagefileUsage', ctypes.c_size_t),
                            ('PeakPagefileUsage', ctypes.c_size_t),
                        ]

                    pmc = _PMC()
                    pmc.cb = ctypes.sizeof(_PMC)
                    handle = ctypes.windll.kernel32.GetCurrentProcess()
                    ok = ctypes.windll.psapi.GetProcessMemoryInfo(
                        handle, ctypes.byref(pmc), pmc.cb)
                    if ok:
                        rss_mb = pmc.WorkingSetSize / (1024 ** 2)
                except Exception:
                    rss_mb = None
            if rss_mb is not None:
                perf_event('mem', 'rss_mb',
                           where=where, mb=round(rss_mb, 1))
        except Exception:
            pass

    def _project_point_scalars_to_subgrid(self, full_scalars, sub_grid):
        """v5.2.1 item 50: map a full-grid per-node scalar array onto a
        sub-grid produced by ``extract_cells``.

        ``extract_cells`` attaches a ``vtkOriginalPointIds`` array to
        the sub-grid that lists the original-grid point indices for
        each retained point. We index ``full_scalars`` through that
        mapping so the sub-grid receives a scalar of the right length
        and ordering -- avoiding the off-by-one trap of naive
        ``full_scalars[:sub_grid.n_points]``.
        """
        import numpy as _np
        if full_scalars is None or sub_grid is None or sub_grid.n_points == 0:
            return None
        try:
            orig_ids = sub_grid.point_data.get('vtkOriginalPointIds')
        except Exception:
            orig_ids = None
        if orig_ids is None:
            n = min(sub_grid.n_points, len(full_scalars))
            out = _np.zeros(sub_grid.n_points, dtype=float)
            out[:n] = _np.asarray(full_scalars, dtype=float)[:n]
            return out
        idx = _np.asarray(orig_ids, dtype=_np.int64)
        idx = _np.clip(idx, 0, len(full_scalars) - 1)
        return _np.asarray(full_scalars, dtype=float)[idx]

    def _attach_per_cell_rgb_for_property_mode(self, grid):
        """v4.0.17: vectorized per-cell RGB lookup for
        ``color_mode='property'``. Writes a ``'cell_rgb'`` uint8
        ``(N, 3)`` array onto ``grid.cell_data`` so
        ``add_mesh(scalars='cell_rgb', rgb=True)`` renders each
        cell with its own ``pid_color_map[PID]`` color, bypassing
        the LUT entirely.

        Replaces the v4.0.5-era
        ``scalars='PID' + cmap + categories=True`` pattern. That
        pattern rendered the same PID with DIFFERENT colors
        across actors with different unique-PID sets. PyVista's
        categorical LUT appears to do range-normalized
        positional mapping rather than value-keyed lookup, so PID
        1234 in pids=[100..1234] mapped to the LAST cmap slot
        while PID 1234 in pids=[1234, 5000] mapped to the FIRST
        slot.

        Per-cell RGB sidesteps the LUT entirely. Same PID, same
        color, every render, every actor.

        Cost: O(n_unique_pids) Python work + numpy fan-out
        (~10 ms on a deck with 2.4M cells). Memory: 3 bytes/cell extra
        (uint8 × 3) vs the int32 PID array → +7 MB per actor on
        production-deck. Acceptable.
        """
        if 'PID' not in grid.cell_data:
            return None
        from node_runner.profiling import perf_event
        pids_arr = np.asarray(grid.cell_data['PID'])
        unique_pids, inverse = np.unique(pids_arr, return_inverse=True)
        # Defensive populate-on-demand. After v4.0.16's reorder
        # pid_color_map is always fully populated before any
        # build runs, but this guarantees correctness regardless
        # of future call-order changes.
        n_populated = 0
        for p in unique_pids:
            pkey = int(p)
            if pkey not in self.pid_color_map:
                self.pid_color_map[pkey] = QColor(
                    random.randint(50, 220),
                    random.randint(50, 220),
                    random.randint(50, 220)).name()
                n_populated += 1
        if n_populated:
            perf_event('rgb_lookup', 'pid_populated_on_demand',
                       n_populated=n_populated)

        def _hex_to_rgb(h):
            h = h.lstrip('#')
            return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))

        unique_rgb = np.array(
            [_hex_to_rgb(self.pid_color_map.get(int(p), "#FFFFFF"))
             for p in unique_pids],
            dtype=np.uint8)
        cell_rgb = unique_rgb[inverse]
        grid.cell_data['cell_rgb'] = cell_rgb
        return cell_rgb

    def _get_moment_glyph_geom(self):
        """v4.0.13: physics-textbook double-tipped moment glyph.

        A moment vector is drawn as a single arrow shaft with a
        DOUBLE arrowhead at one end (▶▶) to distinguish it from a
        force vector (single arrowhead ▶). The shaft direction is
        the right-hand-rule thumb direction; the moment rotates
        with fingers curling around the shaft.

        Geometry: pv.Arrow (shaft + tip cone) merged with a second
        pv.Cone behind the first tip. Both are pure-triangle
        PolyDatas. no tube() filter (which caused the v4.0.11/12
        SEGFAULTs in vtkGlyph3D). The merge of two pure-triangle
        PolyDatas reliably produces a clean PolyData.

        History:
        - v4.0.5: `fwd.merge(back)` symmetric ←→ double-arrow. Right-
          hand rule was ambiguous (axis encoded but not direction).
        - v4.0.11/12: composite (arrow + 270° tube-arc curl + tangent
          cone). Crashed in vtkGlyph3D. `tube()` produces cells
          the glyph filter can't reliably handle.
        - v4.0.13 (this): arrow + second cone at +X end. Simplest
          possible composite that's directional and visually distinct
          from force arrows.
        """
        from node_runner.profiling import perf_event
        cached = getattr(self, '_moment_glyph_geom_cache', None)
        if cached is not None:
            self._moment_glyph_cache_hits = getattr(
                self, '_moment_glyph_cache_hits', 0) + 1
            if self._moment_glyph_cache_hits <= 5:
                perf_event('glyph', 'moment_double_tip_cache_hit',
                           hits=self._moment_glyph_cache_hits)
            return cached
        try:
            import pyvista as pv
            # Main arrow: shaft along +X with single tip at x=1.0.
            # Tip cone occupies x ∈ [0.8, 1.0].
            arrow = pv.Arrow(
                start=(0.0, 0.0, 0.0),
                direction=(1.0, 0.0, 0.0),
                tip_length=0.2, tip_radius=0.08,
                shaft_radius=0.03,
            )
            # Second tip cone behind the first. pv.Cone places its
            # base at `center - direction*height/2` and tip at
            # `center + direction*height/2`. With center=(0.65,0,0)
            # and direction=+X, height=0.2: base at x=0.55, tip at
            # x=0.75. This sits just behind the arrow's tip cone
            # (base at x=0.8, tip at x=1.0). Visually: "▶▶" stack
            # at the +X end with a small gap between bases.
            second_cone = pv.Cone(
                center=(0.65, 0.0, 0.0),
                direction=(1.0, 0.0, 0.0),
                height=0.2, radius=0.08,
            )
            composite_raw = arrow.merge(second_cone)
            if isinstance(composite_raw, pv.PolyData):
                composite = composite_raw
            else:
                composite = composite_raw.extract_surface()
                perf_event('glyph', 'moment_double_tip_extract_surface',
                           raw_type=type(composite_raw).__name__)
            self._moment_glyph_geom_cache = composite
            perf_event('glyph', 'moment_double_tip_built',
                       n_points=int(composite.n_points),
                       n_cells=int(composite.n_cells))
            return composite
        except Exception as _exc:
            perf_event('glyph', 'moment_double_tip_fallback',
                       exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            try:
                import pyvista as pv
                fallback = pv.Arrow(
                    start=(0.0, 0.0, 0.0),
                    direction=(1.0, 0.0, 0.0),
                    tip_length=0.2, tip_radius=0.08,
                    shaft_radius=0.03,
                )
                self._moment_glyph_geom_cache = fallback
                return fallback
            except Exception:
                return None

    @staticmethod
    def _log_shading_actor(tag, actor, actor_name=''):
        """v4.0.1 / v4.0.2: emit actor prop + mapper state to BOTH
        stdout and the profile log so we can diagnose why the
        toggle doesn't change visual appearance.

        Fields logged:
          - lighting (bool: VTK lighting on/off)
          - interpolation (str: 'phong'/'flat'/'gouraud')
          - interp_int (int: VTK enum 0=flat, 1=gouraud, 2=phong)
          - ambient / diffuse / specular (float: material coefficients)
          - normals (str: which arrays carry 'Normals')
          - scalar_visibility (bool: True means the mapper is
            coloring by scalars - lighting/material color is then
            largely cosmetic)
          - color_mode (int: VTK ScalarMode)
        """
        from node_runner.profiling import perf_event
        p = actor.prop
        normals = "n/a"
        scalar_vis = None
        scalar_mode = None
        try:
            mapper = actor.mapper if hasattr(actor, 'mapper') else None
            if mapper is None and hasattr(actor, 'GetMapper'):
                mapper = actor.GetMapper()
            if mapper is not None:
                try:
                    scalar_vis = bool(mapper.GetScalarVisibility())
                except Exception:
                    pass
                try:
                    scalar_mode = int(mapper.GetScalarMode())
                except Exception:
                    pass
                ds = getattr(mapper, 'dataset', None)
                if ds is None and hasattr(mapper, 'GetInput'):
                    try:
                        ds = mapper.GetInput()
                    except Exception:
                        ds = None
                if ds is not None:
                    try:
                        pt_keys = list(ds.point_data.keys()) if hasattr(ds, 'point_data') else []
                        cell_keys = list(ds.cell_data.keys()) if hasattr(ds, 'cell_data') else []
                    except Exception:
                        pt_keys, cell_keys = [], []
                    normals = f"pt={'Normals' in pt_keys} cell={'Normals' in cell_keys}"
        except Exception:
            pass
        try:
            interp_int = int(p.GetInterpolation()) if hasattr(p, 'GetInterpolation') else None
        except Exception:
            interp_int = None
        line = (
            f"[shading] {tag} name={actor_name!r} lighting={p.lighting} "
            f"interp={p.interpolation} interp_int={interp_int} "
            f"ambient={p.ambient:.2f} diffuse={p.diffuse:.2f} "
            f"specular={p.specular:.2f} normals=[{normals}] "
            f"scalar_vis={scalar_vis} scalar_mode={scalar_mode}"
        )
        print(line, flush=True)
        perf_event('shading', tag.strip().lower(),
                   name=actor_name,
                   lighting=p.lighting,
                   interp=p.interpolation,
                   interp_int=interp_int,
                   ambient=f"{p.ambient:.2f}",
                   diffuse=f"{p.diffuse:.2f}",
                   specular=f"{p.specular:.2f}",
                   normals=normals,
                   scalar_vis=scalar_vis,
                   scalar_mode=scalar_mode)

    def _on_shading_display_toggled(self, state):
        """Bridge: Display-tab checkbox -> _toggle_shading (which also
        keeps the View menu action in sync)."""
        self._toggle_shading(state)
        
    def _zoom_to_model(self):
        if self.current_grid: self.plotter.reset_camera(bounds=self.current_grid.bounds)
        
    def _zoom_to_all(self):
        if self.current_grid:
            # --- FIX: Manually calculate the expanded bounding box ---
            # Get the original bounds (xmin, xmax, ymin, ymax, zmin, zmax)
            bounds = np.array(self.current_grid.bounds).reshape(3, 2)
            
            # Calculate the center and size (extent) of the box
            center = np.mean(bounds, axis=1)
            extents = bounds[:, 1] - bounds[:, 0]
            
            # Increase the extents by a 1.1x factor
            new_extents = extents * 1.1
            
            # Calculate the new min/max points from the new center and extents
            new_bounds_arr = np.vstack([
                center - new_extents / 2,
                center + new_extents / 2
            ]).T.ravel()
            
            # Apply the new, expanded bounds to the camera
            self.plotter.reset_camera(bounds=new_bounds_arr)
        
    def _create_nodes(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator:
            self.current_generator = NastranModelGenerator()
            self._ensure_default_coords(self.current_generator.model)
        next_id = self.current_generator.get_next_available_id('node'); dialog = CreateNodesDialog(next_id, self)
        dialog.accepted.connect(self._on_node_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        
    def _on_node_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_creation_parameters()
        if params:
            model = self.current_generator.model
            if params['type'] == 'single':
                nid = params['id'] or self.current_generator.get_next_available_id('node')
                cmd = AddNodesCommand([(nid, *params['coords'])])
                self.command_manager.execute(cmd, model)
                self._update_status(f"Created Node: {nid}")
            elif params['type'] == 'between':
                # Pre-compute the interpolated positions so we can wrap them
                start = model.nodes[params['start_nid']].get_position()
                end = model.nodes[params['end_nid']].get_position()
                num = params['num_nodes']
                node_data = []
                for i in range(1, num + 1):
                    frac = i / (num + 1.0)
                    coord = start + frac * (end - start)
                    nid = self.current_generator.get_next_available_id('node')
                    node_data.append((nid, float(coord[0]), float(coord[1]), float(coord[2])))
                cmd = AddNodesCommand(node_data)
                self.command_manager.execute(cmd, model)
                self._update_status(f"Created {len(node_data)} nodes.")
            self._update_viewer(self.current_generator, reset_camera=False)
        self.active_creation_dialog = None
        
    def _create_line_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateLineElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_line_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        
    def _on_line_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_parameters()
        if params:
            if not params.get('eid'):
                params['eid'] = self.current_generator.get_next_available_id('element')
            cmd = AddElementCommand('line', params)
            self.command_manager.execute(cmd, self.current_generator.model)
            self._update_status(f"Created {params['type']}: {cmd._created_eid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        self._highlight_nodes([])
        self.active_creation_dialog = None
        
    def _create_plate_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreatePlateElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_plate_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        

    def _on_plate_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddElementCommand('plate', params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created {params['type']}: {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None
        
    

# START: New handler methods in main.py
    def _create_bush_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator: return
        
        # Check if any PBUSH properties exist
        has_pbush = any(prop.type == 'PBUSH' for prop in self.current_generator.model.properties.values())
        if not has_pbush:
            QMessageBox.warning(self, "No Bush Properties", "You must create a PBUSH property before creating a CBUSH element.")
            return

        next_eid = self.current_generator.get_next_available_id('element')
        dialog = CreateBushElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_bush_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _on_bush_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddElementCommand('bush', params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created CBUSH: {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None
# END: New handler methods in main.py






    def _create_rbes(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.nodes: QMessageBox.warning(self, "No Nodes", "Create nodes first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateRbeDialog(next_eid, self.current_generator.model.nodes.keys(), self)
        dialog.accepted.connect(self._on_rbe_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _delete_nodes(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "No nodes to delete.")
            return
        if self.selection_bar.isVisible():
            return
        all_node_ids = self.current_generator.model.nodes.keys()
        self._start_selection('Node', all_node_ids, self._on_delete_nodes_accept)

    def _on_delete_nodes_accept(self):
        node_ids = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not node_ids: return

        reply = QMessageBox.question(self, "Confirm Deletion",
                                   f"Are you sure you want to delete {len(node_ids)} selected node(s)? "
                                   "All connected elements will also be deleted.",
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            cmd = DeleteNodesCommand(node_ids)
            self.command_manager.execute(cmd, self.current_generator.model)
            n_nodes = len(cmd._deleted_nodes)
            n_elems = len(cmd._deleted_elements) + len(cmd._deleted_rigid_elements)
            self._update_status(f"Deleted {n_nodes} nodes and {n_elems} elements.")
            self._update_viewer(self.current_generator, reset_camera=False)

    def _delete_elements(self):
        if not self.current_generator or not self.current_grid or self.current_grid.n_cells == 0:
            QMessageBox.warning(self, "No Model", "No elements to delete.")
            return
        if self.selection_bar.isVisible():
            return
        all_eids = self.current_grid.cell_data['EID']
        self._start_selection('Element', all_eids,
                              self._on_delete_elements_accept,
                              self._build_select_by_data())

    def _on_delete_elements_accept(self):
        eids_to_delete = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if not eids_to_delete: return
        
        reply = QMessageBox.question(self, "Confirm Deletion",
                                f"Are you sure you want to delete {len(eids_to_delete)} element(s)?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            cmd = DeleteElementsCommand(eids_to_delete)
            self.command_manager.execute(cmd, self.current_generator.model)
            deleted_count = len(cmd._deleted_elements) + len(cmd._deleted_rigid_elements)
            self._update_status(f"Deleted {deleted_count} elements.")
            self._update_viewer(self.current_generator, reset_camera=False)

    def _delete_mesh(self):
        if not self.current_generator: return
        reply = QMessageBox.question(self, "Confirm Delete Mesh",
                                   "Are you sure you want to delete all nodes and elements?",
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.current_generator.model = BDF(debug=False)
            self._ensure_default_coords(self.current_generator.model)
            self._update_viewer(self.current_generator)
            self._update_status("Mesh deleted.")

    def _delete_coord_system_menu(self):
        """Opens a dialog to select and delete a user-created coordinate system."""
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "No model loaded.")
            return
        model = self.current_generator.model
        user_cids = [cid for cid in sorted(model.coords.keys()) if cid > 2]
        if not user_cids:
            QMessageBox.information(self, "No Coordinate Systems", "There are no user-created coordinate systems to delete.")
            return
        from PySide6.QtWidgets import QInputDialog
        items = [f"CID {cid}" for cid in user_cids]
        item, ok = QInputDialog.getItem(self, "Delete Coordinate System", "Select coordinate system to delete:", items, 0, False)
        if ok and item:
            cid = int(item.split()[1])
            reply = QMessageBox.question(self, "Confirm Deletion",
                f"Are you sure you want to delete Coordinate System {cid}?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteCoordCommand(cid)
                self.command_manager.execute(cmd, model)
                self._update_status(f"Deleted Coordinate System {cid}.")
                self._update_viewer(self.current_generator, reset_camera=False)

    def _show_model_summary(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please load or create a model first.")
            return
        model = self.current_generator.model

        # --- Nodes ---
        num_nodes = len(model.nodes)

        # --- Elements by type ---
        all_elems = {**model.elements, **model.rigid_elements}
        elem_counts = Counter(e.type for e in all_elems.values())

        # --- Properties, Materials, Coords ---
        num_props = len(model.properties)
        num_mats = len(model.materials)
        num_coords = len(model.coords) - 1  # exclude default CID 0

        # --- Loads & Constraints ---
        load_sids = sorted(model.loads.keys()) if model.loads else []
        spc_sids = sorted(model.spcs.keys()) if model.spcs else []

        # --- Masses (CONM2) ---
        total_mass = sum(m.mass for m in model.masses.values()) if model.masses else 0.0

        # --- Build text ---
        sep = "-" * 34
        lines = []
        lines.append(sep)
        lines.append(f"  Nodes          {num_nodes:>10d}")
        lines.append(f"  Elements       {len(all_elems):>10d}")
        lines.append(sep)

        if elem_counts:
            for etype, count in sorted(elem_counts.items()):
                lines.append(f"    {etype:<16s}{count:>8d}")
            lines.append(sep)

        lines.append(f"  Properties     {num_props:>10d}")
        lines.append(f"  Materials      {num_mats:>10d}")
        lines.append(f"  Coord Systems  {max(num_coords, 0):>10d}")

        if model.masses:
            lines.append(sep)
            lines.append(f"  Mass Elements  {len(model.masses):>10d}")
            lines.append(f"  Total Mass     {total_mass:>14.4f}")

        if load_sids or spc_sids:
            lines.append(sep)
            if load_sids:
                lines.append(f"  Load SIDs:  {', '.join(str(s) for s in load_sids)}")
            if spc_sids:
                lines.append(f"  SPC SIDs:   {', '.join(str(s) for s in spc_sids)}")

        lines.append(sep)

        InfoDialog("Model Summary", "\n".join(lines), self).exec()

    def _show_mass_summary(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please load or create a model first.")
            return
        model = self.current_generator.model
        if not model.elements and not model.masses:
            QMessageBox.warning(self, "No Data", "Model has no elements or mass elements.")
            return

        result = self.current_generator.calculate_mass_summary()
        sep = "-" * 38
        lines = []
        lines.append(sep)
        lines.append(f"  Total Mass       {result['total_mass']:>14.4f}")
        lines.append(sep)
        lines.append(f"  Structural Mass  {result['structural_mass']:>14.4f}")
        lines.append(f"  CONM2 Mass       {result['conm2_mass']:>14.4f}")
        lines.append(sep)
        cg = result['cg']
        lines.append(f"  CG  X            {cg[0]:>14.4f}")
        lines.append(f"  CG  Y            {cg[1]:>14.4f}")
        lines.append(f"  CG  Z            {cg[2]:>14.4f}")

        if result['by_pid']:
            lines.append(sep)
            lines.append("  Mass by Property:")
            for pid in sorted(result['by_pid'].keys()):
                lines.append(f"    PID {pid:<8d}   {result['by_pid'][pid]:>14.4f}")

        if result['notes']:
            lines.append(sep)
            for note in result['notes']:
                lines.append(f"  * {note}")

        lines.append(sep)
        InfoDialog("Mass & CG Summary", "\n".join(lines), self).exec()

    # ─── F16: Find/Replace Property/Material ─────────────────────────────

    def _open_data_table_dialog(self):
        """v5.2.0 item 45: open (or focus) the professional Data Table."""
        try:
            from node_runner.dialogs.data_table import DataTableDialog
            if getattr(self, '_data_table_dialog', None) is None:
                self._data_table_dialog = DataTableDialog(self)
            self._data_table_dialog.show()
            self._data_table_dialog.raise_()
            self._data_table_dialog.activateWindow()
            try:
                from node_runner.profiling import perf_event
                perf_event(
                    'data_table', 'open',
                    n_selected_nids=len(
                        getattr(self, '_current_selection_nids', []) or []),
                    n_selected_eids=len(
                        getattr(self, '_current_selection_eids', []) or []))
            except Exception:
                pass
        except Exception as e:
            QMessageBox.critical(self, "Data Table failed", str(e))

    def _open_find_replace_dialog(self):
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Model", "Please load a model with elements first.")
            return

        model = self.current_generator.model
        dialog = FindReplaceDialog(model, self)
        if dialog.exec() == QDialog.Accepted:
            params = dialog.get_parameters()
            if params is None:
                return

            if params['mode'] == 'property':
                # Find all elements with matching PID
                eids = [eid for eid, e in model.elements.items()
                        if hasattr(e, 'pid') and e.pid == params['from_id']]
                if not eids:
                    self._update_status("No elements found with that property.", is_error=True)
                    return
                cmd = ReassignPropertyCommand(eids, params['to_id'])
                self.command_manager.execute(cmd, model)
                self._update_status(f"Reassigned {len(eids)} elements from PID {params['from_id']} to PID {params['to_id']}.")
            else:
                # Find all properties with matching MID
                pids = [pid for pid, p in model.properties.items()
                        if (hasattr(p, 'mid') and p.mid == params['from_id']) or
                           (hasattr(p, 'mid1') and p.mid1 == params['from_id'])]
                if not pids:
                    self._update_status("No properties found with that material.", is_error=True)
                    return
                cmd = ReassignMaterialCommand(pids, params['to_id'])
                self.command_manager.execute(cmd, model)
                self._update_status(f"Reassigned {len(pids)} properties from MID {params['from_id']} to MID {params['to_id']}.")

            self._update_viewer(self.current_generator, reset_camera=False)
            self._populate_tree()

    # ─── v3.5.0 item 3: Preferences (colors / sizes / highlight) ────────

    def _open_preferences_dialog(self):
        """Open the Preferences dialog. On accept, write to QSettings,
        apply the new values to live attributes, and rebuild the viewer
        so colors/sizes/highlight update immediately."""
        from node_runner.dialogs.preferences import (
            PreferencesDialog, save_preferences, load_preferences)
        current = load_preferences()
        # Overlay any in-memory overrides from the running session so
        # the dialog reflects the live state.
        for key, value in (self.type_color_map or {}).items():
            current.setdefault('colors', {})[key] = value
        if getattr(self, 'highlight_color', None):
            current['highlight_color'] = self.highlight_color
        dlg = PreferencesDialog(current, self)
        if dlg.exec() != QDialog.Accepted:
            return
        payload = dlg.result_payload()
        save_preferences(payload)
        # Apply: live attributes -> viewer rebuild.
        for key, value in (payload.get('colors') or {}).items():
            self.type_color_map[key] = value
        self.highlight_color = payload.get(
            'highlight_color', getattr(self, 'highlight_color', '#fab387'))
        sizes = payload.get('sizes') or {}
        # v4.0.7: default percent lowered from 1.5% to 0.75%.
        self.mass_glyph_scale = float(
            sizes.get('mass_glyph_scale_pct', 0.75)) / 100.0
        self.node_size = int(sizes.get('node_size', self.node_size))
        self.beam_width = int(sizes.get('beam_width', self.beam_width))
        self.edge_width = int(sizes.get('edge_width', self.edge_width))
        self.rbe_line_width = int(
            sizes.get('rbe_line_width', getattr(self, 'rbe_line_width', 3)))
        self.free_edge_width = int(
            sizes.get('free_edge_width', getattr(self, 'free_edge_width', 4)))
        self.highlight_outline_width = int(
            sizes.get('highlight_outline_width',
                      getattr(self, 'highlight_outline_width', 5)))
        # v5.0.0 item 9: brightness mode invalidates the cached
        # multiplier so the next light-rig install reads the new value.
        self._invalidate_brightness_cache()
        try:
            from node_runner.profiling import perf_event
            perf_event('lights', 'brightness_changed',
                       to=payload.get('brightness_mode', 'Normal'))
        except Exception:
            pass
        # v5.3.2 item 69: scalar bar position preference.
        new_orient = payload.get('results_bar_orientation', 'right')
        if new_orient not in ('right', 'bottom'):
            new_orient = 'right'
        self._results_bar_orientation = new_orient
        # Re-render the results pipeline so the bar moves immediately.
        if self.color_mode == 'results':
            try:
                self._clear_scalar_bars()
                self._on_results_changed()
            except Exception:
                pass
        if self.current_generator:
            self._update_viewer(self.current_generator, reset_camera=False)
        self._update_status("Preferences applied - viewer updated.")

    # ─── F8: LOAD Combination ────────────────────────────────────────────

    def _create_load_combination(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please load a model first.")
            return
        model = self.current_generator.model
        if not model.loads:
            QMessageBox.warning(self, "No Loads", "No load sets exist to combine.")
            return

        dialog = CreateLoadCombinationDialog(model, self)
        if dialog.exec() == QDialog.Accepted:
            params = dialog.get_parameters()
            if params is None:
                return
            cmd = AddLoadCombinationCommand(params)
            self.command_manager.execute(cmd, model)
            self._populate_tree()
            self._populate_loads_tab()
            self._update_status(f"Created LOAD combination SID {params['sid']}.")

    def _edit_load_combination(self, combo_sid):
        """v3.5.0 item 2: open the professional edit dialog on an
        existing load combination. Tab + Model tree refresh on accept."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        if combo_sid not in (model.load_combinations or {}):
            return
        from node_runner.dialogs.load_combination import LoadCombinationDialog
        dlg = LoadCombinationDialog(
            model, mode=LoadCombinationDialog.MODE_EDIT,
            combo_sid=combo_sid, parent=self)
        if dlg.exec() != QDialog.Accepted:
            return
        new_sid, payload = dlg.result_payload()
        cmd = EditLoadCombinationCommand(
            combo_sid, new_sid, payload)
        self.command_manager.execute(cmd, model)
        self._populate_tree()
        self._populate_loads_tab()
        if new_sid == combo_sid:
            self._update_status(f"Edited LOAD combination SID {new_sid}.")
        else:
            self._update_status(
                f"Renamed LOAD combination SID {combo_sid} -> {new_sid}.")

    def _copy_load_combination(self, combo_sid):
        """v3.5.0 item 2: open the dialog pre-populated from an existing
        combo but with a fresh SID. Source untouched on accept."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        if combo_sid not in (model.load_combinations or {}):
            return
        from node_runner.dialogs.load_combination import LoadCombinationDialog
        dlg = LoadCombinationDialog(
            model, mode=LoadCombinationDialog.MODE_COPY,
            combo_sid=combo_sid, parent=self)
        if dlg.exec() != QDialog.Accepted:
            return
        new_sid, payload = dlg.result_payload()
        cmd = CopyLoadCombinationCommand(new_sid, payload)
        self.command_manager.execute(cmd, model)
        self._populate_tree()
        self._populate_loads_tab()
        self._update_status(
            f"Copied LOAD combination SID {combo_sid} -> {new_sid}.")

    def _rename_load_combination(self, combo_sid):
        """v3.5.0 item 2: quick SID-only rename via QInputDialog. For
        a full edit use _edit_load_combination."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        if combo_sid not in (model.load_combinations or {}):
            return
        existing = set((model.load_combinations or {}).keys()) \
                   | set((model.loads or {}).keys())
        new_sid, ok = QInputDialog.getInt(
            self, f"Rename LOAD SID {combo_sid}",
            "New SID:",
            value=combo_sid, min=1, max=99999999)
        if not ok or new_sid == combo_sid:
            return
        if new_sid in existing:
            QMessageBox.warning(
                self, "SID conflict",
                f"SID {new_sid} is already used.")
            return
        payload = dict(model.load_combinations[combo_sid])
        cmd = EditLoadCombinationCommand(combo_sid, new_sid, payload)
        self.command_manager.execute(cmd, model)
        self._populate_tree()
        self._populate_loads_tab()
        self._update_status(
            f"Renamed LOAD combination SID {combo_sid} -> {new_sid}.")

    # ─── F5: Subcase / Case Control ──────────────────────────────────────

    def _open_subcase_editor(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please load a model first.")
            return
        model = self.current_generator.model
        eigrl_sids = [e['sid'] for e in self.eigrl_cards]
        dialog = SubcaseEditorDialog(model, self.subcases, eigrl_sids=eigrl_sids, parent=self)
        if dialog.exec() == QDialog.Accepted:
            result = dialog.get_subcases()
            if result is not None:
                self.subcases = result
                self._update_status(f"Defined {len(self.subcases)} subcase(s) for case control.")

    def _apply_case_control(self):
        """Applies the subcase definitions to the model's case control deck before saving."""
        if not self.current_generator:
            return
        model = self.current_generator.model

        # Check if we have an active analysis set - use it if available
        active_set = self.analysis_sets.get(self.active_analysis_set_id)
        if active_set:
            return self._apply_case_control_from_set(model, active_set)

        # --- Legacy behavior (no analysis set) ---

        # Set SOL type
        if self.sol_type:
            model.sol = self.sol_type

        # Add EIGRL cards to bulk data
        for eigrl in self.eigrl_cards:
            try:
                model.add_eigrl(eigrl['sid'], v1=eigrl.get('v1'),
                                v2=eigrl.get('v2'), nd=eigrl.get('nd'))
            except Exception:
                pass  # May already exist from previous save

        if not self.subcases:
            model.case_control_deck = None
            return

        from pyNastran.bdf.case_control_deck import CaseControlDeck
        lines = ['CEND']
        model.case_control_deck = CaseControlDeck(lines, log=None)
        for sc in self.subcases:
            subcase = model.case_control_deck.create_new_subcase(sc['id'])
            if sc.get('load_sid'):
                subcase.add('LOAD', sc['load_sid'], [], 'STRESS-type')
            if sc.get('spc_sid'):
                subcase.add('SPC', sc['spc_sid'], [], 'STRESS-type')
            if sc.get('method_sid'):
                subcase.add('METHOD', sc['method_sid'], [], 'STRESS-type')
            if sc.get('statsub_id'):
                subcase.add('STATSUB', sc['statsub_id'], [], 'STRESS-type')
            # Output requests
            for req_key, req_name in [('disp', 'DISPLACEMENT'), ('stress', 'STRESS'),
                                       ('force', 'FORCE'), ('strain', 'STRAIN')]:
                val = sc.get(req_key, 'NONE')
                if val == 'ALL':
                    subcase.add(req_name, 'ALL', ['SORT1'], 'STRESS-type')

    def _apply_case_control_from_set(self, model, active_set):
        """Apply case control from an AnalysisSet object."""
        # 1. SOL type
        if active_set.sol_type:
            model.sol = active_set.sol_type

        # 2. EIGRL card
        if active_set.eigrl:
            try:
                model.add_eigrl(
                    active_set.eigrl['sid'],
                    v1=active_set.eigrl.get('v1'),
                    v2=active_set.eigrl.get('v2'),
                    nd=active_set.eigrl.get('nd'))
            except Exception:
                pass

        # 3. PARAM cards
        for param_name, param_value in active_set.params.items():
            try:
                model.add_card(['PARAM', param_name, param_value], 'PARAM')
            except Exception:
                pass

        # 4. Case control deck
        subcases = active_set.subcases
        if not subcases:
            # No subcases - write global output requests only
            from pyNastran.bdf.case_control_deck import CaseControlDeck
            lines = ['CEND']
            model.case_control_deck = CaseControlDeck(lines, log=None)
            # Global output requests
            output_map = {
                'displacement': 'DISPLACEMENT', 'velocity': 'VELOCITY',
                'acceleration': 'ACCELERATION', 'stress': 'STRESS',
                'strain': 'STRAIN', 'force': 'FORCE',
                'spcforce': 'SPCFORCES', 'mpcforce': 'MPCFORCES',
                'gpforce': 'GPFORCE', 'oload': 'OLOAD',
            }
            for key, cc_name in output_map.items():
                val = active_set.output_requests.get(key, 'NONE')
                if val == 'ALL':
                    model.case_control_deck.add_parameter_to_global_subcase(
                        f'{cc_name}(PLOT) = ALL')
            return

        from pyNastran.bdf.case_control_deck import CaseControlDeck
        lines = ['CEND']
        model.case_control_deck = CaseControlDeck(lines, log=None)

        for sc in subcases:
            subcase = model.case_control_deck.create_new_subcase(sc['id'])
            if sc.get('load_sid'):
                subcase.add('LOAD', sc['load_sid'], [], 'STRESS-type')
            if sc.get('spc_sid'):
                subcase.add('SPC', sc['spc_sid'], [], 'STRESS-type')
            if sc.get('method_sid'):
                subcase.add('METHOD', sc['method_sid'], [], 'STRESS-type')
            if sc.get('statsub_id'):
                subcase.add('STATSUB', sc['statsub_id'], [], 'STRESS-type')

            # Per-subcase output requests (from subcase dict if present)
            for req_key, req_name in [('disp', 'DISPLACEMENT'), ('stress', 'STRESS'),
                                       ('force', 'FORCE'), ('strain', 'STRAIN')]:
                val = sc.get(req_key, 'NONE')
                if val == 'ALL':
                    subcase.add(req_name, 'ALL', ['SORT1'], 'STRESS-type')

            # Additional global output requests from the analysis set
            extra_map = {
                'velocity': 'VELOCITY', 'acceleration': 'ACCELERATION',
                'spcforce': 'SPCFORCES', 'mpcforce': 'MPCFORCES',
                'gpforce': 'GPFORCE', 'oload': 'OLOAD',
            }
            for key, cc_name in extra_map.items():
                val = active_set.output_requests.get(key, 'NONE')
                if val == 'ALL':
                    subcase.add(cc_name, 'ALL', ['SORT1'], 'STRESS-type')

    # ─── F6: Renumber Nodes/Elements ─────────────────────────────────────

    def _open_renumber_dialog(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "Please load a model first.")
            return

        model = self.current_generator.model
        dialog = RenumberDialog(model, self)
        if dialog.exec() == QDialog.Accepted:
            params = dialog.get_parameters()
            if params is None:
                return
            if not params['do_nodes'] and not params['do_elements']:
                self._update_status("Nothing selected to renumber.", is_error=True)
                return

            cmd = RenumberCommand(
                params['start_nid'], params['start_eid'],
                params['do_nodes'], params['do_elements']
            )
            self.command_manager.execute(cmd, model)
            # Clear groups since all IDs changed
            self.groups.clear()
            self._hidden_groups.clear()
            self._isolate_mode = None
            self._populate_groups_list()
            self._update_viewer(self.current_generator, reset_camera=False)
            self._populate_tree()
            self._update_status("IDs renumbered successfully.")

    def _list_node_info(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "No nodes to list.")
            return
        if self.selection_bar.isVisible():
            return
        all_node_ids = self.current_generator.model.nodes.keys()
        self._start_selection('Node', all_node_ids, self._on_list_nodes_accept)

    def _list_element_info(self):
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Model", "No elements to list.")
            return
        if self.selection_bar.isVisible():
            return
        all_elements = {**self.current_generator.model.elements, **self.current_generator.model.rigid_elements}
        if not all_elements:
            QMessageBox.warning(self, "No Model", "No elements to list.")
            return
        self._start_selection('Element', all_elements.keys(),
                              self._on_list_elements_accept,
                              self._build_select_by_data())

    def _on_rbe_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_parameters()
        if params:
            if not params.get('eid'):
                params['eid'] = self.current_generator.get_next_available_id('element')
            cmd = AddElementCommand('rbe', params)
            self.command_manager.execute(cmd, self.current_generator.model)
            self._update_status(f"Created {params.get('elem_type', 'RBE')}: {cmd._created_eid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        self.active_creation_dialog = None

    def _create_conm2(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Nodes", "Create nodes first."); return
        next_eid = self.current_generator.get_next_available_id('element')
        dialog = CreateConm2Dialog(next_eid, self.current_generator.model.nodes.keys(), self)
        dialog.accepted.connect(self._on_conm2_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _on_conm2_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_parameters()
        if params:
            if not params.get('eid'):
                params['eid'] = self.current_generator.get_next_available_id('element')
            cmd = AddMassCommand(params)
            self.command_manager.execute(cmd, self.current_generator.model)
            self._update_status(f"Created CONM2 (mass={params['mass']}) at node {params['nid']}: EID {cmd._created_eid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        self.active_creation_dialog = None

    def _create_solid_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateSolidElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_solid_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _on_solid_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddElementCommand('solid', params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created {params['type']}: EID {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None

    def _create_shear_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateShearElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_shear_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _on_shear_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddElementCommand('shear', params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created CSHEAR: EID {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None

    def _create_gap_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateGapElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_gap_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _on_gap_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddElementCommand('gap', params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created CGAP: EID {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None

    def _create_plotels(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Nodes", "Create nodes first."); return
        next_eid = self.current_generator.get_next_available_id('element')
        dialog = CreatePlotelDialog(next_eid, self.current_generator.model.nodes.keys(), self)
        dialog.accepted.connect(self._on_plotel_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _on_plotel_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                if not params.get('eid'):
                    params['eid'] = self.current_generator.get_next_available_id('element')
                cmd = AddPlotelCommand(params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_status(f"Created PLOTEL: Nodes {params['nodes'][0]}-{params['nodes'][1]}, EID {cmd._created_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None

    def _on_creation_reject(self):
        self._update_status("Creation cancelled.")
        self._highlight_entities('Node', [])
        self.active_creation_dialog = None

    def _open_node_move_tool(self):
        if self.selection_bar.isVisible():
            return
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        all_node_ids = self.current_generator.model.nodes.keys()
        self._start_selection('Node', all_node_ids, self._on_move_nodes_accept)

    # ------------------------------------------------------------------
    # Selection bar helpers
    # ------------------------------------------------------------------

    def _start_selection(self, entity_type, all_ids, accept_callback,
                         select_by_data=None):
        """Start a selection workflow using the floating selection window.

        v3.5.0 item 1: by default the bar's pool is intersected with
        the tree-visible entity set so users can't accidentally pick
        something hidden. The full pool is passed as
        ``all_entity_ids_unfiltered`` so the bar's
        "Include hidden entities" checkbox can restore it on demand.
        """
        if self.selection_bar.isVisible():
            return  # already in a selection workflow
        self.current_selection_type = entity_type
        full_pool = set(all_ids) if all_ids is not None else set()
        try:
            visible_pool = full_pool & set(self.visible_entity_ids(entity_type))
        except Exception:
            visible_pool = full_pool
        # If nothing is visible (e.g. user turned everything off and
        # then triggered a list action), fall back to the full pool
        # rather than handing the bar an empty selection.
        default_pool = visible_pool if visible_pool else full_pool
        self.selection_bar.configure(
            entity_type, default_pool,
            select_by_data=select_by_data,
            groups=self.groups,
            all_entity_ids_unfiltered=full_pool,
        )
        self.selection_bar.show()
        self.active_selection_dialog = self.selection_bar
        self._selection_accept_callback = accept_callback

    def _on_selection_bar_accepted(self):
        """Route the bar's OK signal to the stored callback."""
        if self._selection_accept_callback:
            self._selection_accept_callback()

    def _on_selection_bar_rejected(self):
        """Handle the bar's Cancel signal."""
        # If we were selecting for a creation dialog, re-show it
        if getattr(self, '_selection_reject_restore_creation', False):
            self._selection_reject_restore_creation = False
            if self.active_creation_dialog:
                self.active_creation_dialog.show()
        self._update_status("Selection cancelled.")
        self._end_selection_mode()

    def _end_selection_mode(self):
        """Hide the selection bar and clean up selection state."""
        bar = self.selection_bar
        if bar.isVisible():
            entity_type = bar.entity_type
            if entity_type == 'Node' and bar.selected_ids:
                self._previous_node_selection = set(bar.selected_ids)
            elif entity_type == 'Element' and bar.selected_ids:
                self._previous_element_selection = set(bar.selected_ids)
            bar.hide()

        self.active_selection_dialog = None
        self._selection_accept_callback = None

        if self.current_selection_type:
            self._highlight_entities(self.current_selection_type, [])
            self._set_picking_mode(self.current_selection_type, False)

        self.current_selection_type = None

    def _on_list_nodes_accept(self):
        selected_ids = self.selection_bar.get_selected_ids()
        if selected_ids:
            from node_runner.dialogs.info import NodeElementInfoDialog
            columns = ['NID', 'X', 'Y', 'Z', 'CP', 'CD', 'PS', 'SEID']
            tooltips = {
                'NID': 'Node ID - unique identifier for this grid point',
                'X': 'X coordinate in the position coordinate system (CP)',
                'Y': 'Y coordinate in the position coordinate system (CP)',
                'Z': 'Z coordinate in the position coordinate system (CP)',
                'CP': 'Position Coordinate System - the coordinate system in which X, Y, Z are defined. 0 = basic rectangular.',
                'CD': 'Displacement Coordinate System - the coordinate system for displacement, force, and constraint output. 0 = basic rectangular.',
                'PS': 'Permanent Single-point Constraints - DOFs permanently constrained at this node (digits 1-6 map to T1, T2, T3, R1, R2, R3).',
                'SEID': 'Superelement ID - 0 means the node belongs to the residual structure.',
            }
            rows = []
            for nid in selected_ids:
                if nid in self.current_generator.model.nodes:
                    node = self.current_generator.model.nodes[nid]
                    rows.append({
                        'NID': node.nid,
                        'X': round(node.xyz[0], 8), 'Y': round(node.xyz[1], 8), 'Z': round(node.xyz[2], 8),
                        'CP': node.cp, 'CD': node.cd,
                        'PS': str(node.ps) if node.ps else '',
                        'SEID': node.seid,
                    })
            title = f"Node Information ({len(rows)} node{'s' if len(rows) != 1 else ''})"
            info_dialog = NodeElementInfoDialog(title, rows, columns, self, column_tooltips=tooltips)
            info_dialog.exec()

        self._end_selection_mode()

    def _on_list_elements_accept(self):
        selected_ids = self.selection_bar.get_selected_ids()
        if selected_ids:
            from node_runner.dialogs.info import NodeElementInfoDialog
            model = self.current_generator.model
            all_elements = {**model.elements, **model.rigid_elements,
                            **getattr(model, 'masses', {})}
            columns = ['EID', 'Type', 'PID', 'MID', 'Nodes']
            tooltips = {
                'EID': 'Element ID - unique identifier for this element',
                'Type': 'Element type (e.g. CQUAD4, CTRIA3, CBEAM, CBUSH, RBE2)',
                'PID': 'Property ID - references the property card that defines element thickness, material, etc.',
                'MID': 'Material ID - the material referenced by this element\'s property',
                'Nodes': 'Grid point IDs defining this element\'s connectivity',
            }
            rows = []
            for eid in selected_ids:
                if eid not in all_elements:
                    continue
                elem = all_elements[eid]
                t = elem.type
                pid_disp = '-'
                mid_disp = '-'
                # v3.4.0 item 5: type-dispatched PID + MID derivation.
                pid_raw = getattr(elem, 'pid', None)
                if t not in ('RBE2', 'RBE3', 'CONM2') and pid_raw:
                    pid_disp = pid_raw
                    if pid_raw in model.properties:
                        prop = model.properties[pid_raw]
                        m_id = (getattr(prop, 'mid', None)
                                or getattr(prop, 'mid1', None))
                        if m_id:
                            mid_disp = m_id
                # v3.4.0 item 5: type-dispatched node-string. RBE2/RBE3
                # don't expose elem.nodes; CONM2 has a single nid.
                nodes_str = self._format_element_nodes(elem)
                rows.append({
                    'EID': elem.eid, 'Type': t,
                    'PID': pid_disp, 'MID': mid_disp, 'Nodes': nodes_str,
                })
            title = f"Element Information ({len(rows)} element{'s' if len(rows) != 1 else ''})"
            info_dialog = NodeElementInfoDialog(title, rows, columns, self, column_tooltips=tooltips)
            info_dialog.exec()

        self._end_selection_mode()

    @staticmethod
    def _format_element_nodes(elem):
        """v3.4.0 item 5: type-dispatched node-string for the Element
        Information dialog. RBE2 / RBE3 / CONM2 don't have a uniform
        `.nodes` attribute, so reading it blindly (the v3.3.x bug)
        produced empty Nodes columns for the most common rigid-body
        constraints in any aerospace deck."""
        t = getattr(elem, 'type', None)
        try:
            if t == 'RBE2':
                legs = ', '.join(str(n) for n in (elem.Gmi or []) if n)
                return f"{elem.gn} -> [{legs}]"
            if t == 'RBE3':
                legs = ', '.join(
                    str(n) for wt in (elem.wt_cg_groups or [])
                    for n in (wt[2] or []) if n)
                return f"[{legs}] -> {elem.refgrid}"
            if t == 'CONM2':
                return str(getattr(elem, 'nid', ''))
            if hasattr(elem, 'nodes'):
                return ', '.join(str(n) for n in elem.nodes if n)
        except Exception:
            return ''
        return ''

    def _on_move_nodes_accept(self):
        """Handle node move: open transform dialog after selection."""
        selected_ids = self.selection_bar.get_selected_ids()
        self._end_selection_mode()
        if selected_ids:
            transform_dialog = NodeTransformDialog(self)
            if transform_dialog.exec() and (params := transform_dialog.get_transform_parameters()):
                cmd = TransformNodesCommand(selected_ids, params)
                self.command_manager.execute(cmd, self.current_generator.model)
                self._update_viewer(self.current_generator)
                self._update_status(f"Moved {len(selected_ids)} nodes.")
            else:
                self._update_status("Node transform cancelled.")

    def _open_element_editor(self):
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model contains no elements to edit.")
            return

        dialog = ElementEditorDialog(self.current_generator.model, self)
        # Snapshot the element before the dialog mutates it
        pre_edit_snapshot = None
        original_accept = dialog.accept
        def accept_with_snapshot():
            nonlocal pre_edit_snapshot
            if dialog.current_eid and dialog.current_eid in self.current_generator.model.elements:
                pre_edit_snapshot = copy.deepcopy(self.current_generator.model.elements[dialog.current_eid])
            original_accept()
        dialog.accept = accept_with_snapshot

        if dialog.exec():
            eid = dialog.current_eid
            if eid and pre_edit_snapshot is not None:
                # The dialog already mutated the element. Create a command
                # that stores the old snapshot and can redo/undo.
                elem = self.current_generator.model.elements[eid]
                cmd = EditElementCommand(eid, elem.pid, list(elem.nodes),
                                         dialog.orientation_widget.get_orientation() if dialog.orientation_widget else None)
                # Manually set old state instead of executing (already done by dialog)
                cmd._old_element = pre_edit_snapshot
                self.command_manager.undo_stack.append(cmd)
                if len(self.command_manager.undo_stack) > self.command_manager.max_history:
                    self.command_manager.undo_stack.pop(0)
                self.command_manager.redo_stack.clear()

            self._update_status(f"Element {eid} modified.")
            self._update_viewer(self.current_generator)

# Location: main.py, around line 2862
# --- MODIFY THE METHOD SIGNATURE AND ADD ONE LINE ---
    def _activate_single_node_picker(self, target_callback, calling_dialog=None):
        self.picking_target_callback = target_callback
        self.active_sub_dialog = calling_dialog # Store a reference to the dialog
        self._set_picking_mode("Node", True)
        self._update_status("PICKING MODE: Click a single node.")

    def _highlight_entities(self, entity_type, entity_ids):
        camera_before = self.plotter.camera.copy()

        self.plotter.remove_actor('selection_highlight', render=False)

        # Phase 8: Track current selection for group management
        if entity_type == 'Node':
            self._current_selection_nids = list(entity_ids) if entity_ids else []
        elif entity_type == 'Element':
            self._current_selection_eids = list(entity_ids) if entity_ids else []
        # v5.2.0 item 45: notify subscribers (Data Table dialog) that
        # the live selection changed.
        try:
            self.selection_changed.emit(
                entity_type, list(entity_ids) if entity_ids else [])
        except Exception:
            pass

        # Theme B3: sync to result browser. When a single entity is
        # highlighted, scroll the matching row into view + select it.
        # v5.1.2: the tables now live inside the Results sidebar tab;
        # no need to gate by .isVisible() since it's always available.
        if (self._result_browser_dock is not None
                and entity_ids and len(entity_ids) == 1):
            try:
                if entity_type == 'Node':
                    self._result_browser_dock.select_node(int(entity_ids[0]))
                elif entity_type == 'Element':
                    self._result_browser_dock.select_element(int(entity_ids[0]))
            except Exception:
                pass

        if not entity_ids or not self.current_generator or not self.current_grid:
            self.plotter.render()
            return

        if entity_type == 'Node':
            coords = []
            for nid in entity_ids:
                if nid in self.current_generator.model.nodes:
                    coords.append(self.current_generator.model.nodes[nid].get_position())

            if coords:
                # v3.5.0 item 3: highlight color is now user-configurable
                # via Preferences > Highlight.
                hl_color = getattr(self, 'highlight_color', '#fab387')
                self.plotter.add_points(
                    np.array(coords),
                    color=hl_color,
                    point_size=12,
                    render_points_as_spheres=True,
                    name='selection_highlight',
                    reset_camera=False,
                )

        elif entity_type == 'Element':
            # v3.2.4: extract_cells deep-copies the selected cells +
            # their points. On a 700k-element 'select all' it's many
            # seconds and looks like a hang. Performance guard:
            #   * If the selection size approaches the full element
            #     count, skip the deep copy and overlay the WHOLE
            #     current_grid as a wireframe with the accent color.
            #     Visually equivalent for the user.
            #   * Otherwise extract the subset as before.
            # v3.5.0 item 3: highlight color + outline width are user prefs.
            hl_color = getattr(self, 'highlight_color', '#fab387')
            hl_width = int(getattr(self, 'highlight_outline_width', 5))
            n_cells = self.current_grid.n_cells
            n_sel = len(entity_ids)
            if n_cells > 0 and n_sel >= n_cells * 0.5:
                # Big selection - just overlay the full grid.
                self.plotter.add_mesh(
                    self.current_grid, style='wireframe',
                    color=hl_color, line_width=max(3, hl_width - 2),
                    name='selection_highlight', pickable=False,
                    reset_camera=False,
                )
            else:
                indices = np.isin(self.current_grid.cell_data['EID'], entity_ids)
                if np.any(indices):
                    highlight_grid = self.current_grid.extract_cells(indices)
                    self.plotter.add_mesh(
                        highlight_grid, style='wireframe',
                        color=hl_color, line_width=hl_width,
                        name='selection_highlight', pickable=False,
                        reset_camera=False,
                    )

        self.plotter.camera = camera_before
        self.plotter.render()

    def _update_rubber_band(self, start_vtk, end_vtk):
        """Show/update the box selection overlay (VTK display coords)."""
        self._selection_overlay.update_box(start_vtk, end_vtk)

    def _update_circle_overlay(self, center_vtk, radius):
        """Show/update the circle selection overlay (VTK display coords)."""
        self._selection_overlay.update_circle(center_vtk, radius)

    def _update_polygon_overlay(self, polygon_points_vtk, cursor_vtk=None):
        """Show/update the polygon selection overlay (VTK display coords)."""
        self._selection_overlay.update_polygon(polygon_points_vtk, cursor_vtk)

    def _hide_rubber_band(self):
        """Hide the selection overlay."""
        if hasattr(self, '_selection_overlay'):
            self._selection_overlay.hide_overlay()

    def _set_picking_mode(self, entity_type, enabled):
        if enabled:
            if not self.default_interactor_style:
                self.default_interactor_style = self.plotter.interactor_style
            mode_label = self._pick_selection_mode.capitalize()
            self._update_status(f"PICKING MODE ({mode_label}): Left=Select, Middle=Rotate, Ctrl+Middle=Pan, Scroll=Zoom. ENTER=accept, ESC=cancel.")
            self.plotter.setCursor(QtCore.Qt.CrossCursor)
            self.plotter.add_text(f"PICKING MODE ({mode_label}): Left=Select, Middle=Rotate, Ctrl+Mid=Pan. ENTER=accept, ESC=cancel", position='upper_right', color='yellow', font_size=8, name='picking_mode_text')
            interactor = ClickAndDragInteractor(main_window=self)
            interactor._selection_mode = self._pick_selection_mode  # Apply stored mode
            self.plotter.interactor.SetInteractorStyle(interactor)
        else:
            self.plotter.setCursor(QtCore.Qt.ArrowCursor)
            self.plotter.remove_actor('picking_mode_text')
            if self.default_interactor_style:
                self.plotter.interactor.SetInteractorStyle(self.default_interactor_style)
            # Selection bar is embedded - no need to show/activate it after picking
            if self.active_sub_dialog and not self.active_sub_dialog.isVisible():
                self.active_sub_dialog.show()
                self.active_sub_dialog.activateWindow()
                self.active_sub_dialog = None # Clear reference after showing
            if self.active_creation_dialog and not self.active_creation_dialog.isVisible():
                self.active_creation_dialog.show()
                self.active_creation_dialog.activateWindow()
            self.picking_target_callback = None
            if self.current_selection_type is None:
                self._update_status("Picking mode disabled.")

    def _request_disable_picking_mode(self):
        QtCore.QTimer.singleShot(0, self._disable_picking_mode)

    def _request_accept_picking_mode(self):
        QtCore.QTimer.singleShot(0, self._accept_picking_mode)

    def _accept_picking_mode(self):
        if self.current_selection_type:
            self._set_picking_mode(self.current_selection_type, False)
        if self.selection_bar.isVisible():
            self._update_status("Selection complete. Click OK to proceed.")

    def _disable_picking_mode(self):
        if self.current_selection_type:
            self._set_picking_mode(self.current_selection_type, False)
        else:
            self._set_picking_mode('Node', False)


# START: Corrected replacement for _update_viewer method in main.py
    def _update_viewer(self, generator, reset_camera=True):
        # v4.0.1 (Stage 1): perf instrumentation. Imports the lightweight
        # profiling helper; gated on NR_PROFILE=1.
        from node_runner.profiling import perf_stage, perf_event, log_deck_header
        _t_total = time.perf_counter()

        with perf_stage('viewer', 'clear_and_axes'):
            self.quality_results = None
            self.plotter.clear()
            self.plotter.add_axes()
            # v4.0.8: plotter.clear() above removes all renderer lights
            # along with actors. Without this re-install the scene
            # has zero lights after every model load and the shading
            # toggle (lighting=True/False) renders identically. The
            # rig in __init__ is no longer sufficient on its own.
            self._install_light_rig()
            # v4.0.14 (Fix C v2): plotter.clear() also dropped the
            # pre-built per-category actor handles. Invalidate the
            # dict so the build site below starts fresh.
            self._category_actors = {}
            self._category_actor_grids = {}
            self._category_visible_ntypes = set()
            self._using_legacy_mesh = False
            # v5.0.0 item 11a: PID/MID hidden sets are PER-deck; clear
            # on viewer rebuild so a freshly opened model starts with
            # every property/material visible.
            self._hidden_pids = set()
            self._hidden_mids = set()
            self.axes_actor = self.plotter.renderer.axes_actor
            self._set_axes_label_color((0.804, 0.839, 0.957) if self.is_dark_theme else (0, 0, 0))

        self.current_generator, self.current_grid = generator, None
        model = generator.model
        # Phase 1: refresh persistent status-bar widgets whenever the viewer
        # is rebuilt - covers undo/redo, generator runs, and CAD imports.
        self._refresh_status_widgets()

        if not model.nodes:
            self.tree_widget.blockSignals(True)
            self._populate_tree()
            self.tree_widget.blockSignals(False)
            self.current_grid = pv.UnstructuredGrid()
            # Coord actors are handled centrally by _refresh_coord_actors
            self._update_plot_visibility()
            return

        self.current_node_ids_sorted = sorted(model.nodes.keys())
        node_map = {nid: i for i, nid in enumerate(self.current_node_ids_sorted)}
        # v3.2.0: bulk-build node coords via the vectorized helper. The
        # old code called node.get_position() per node which is ~50x
        # slower because of pyNastran's per-call matmul. The vectorized
        # path groups nodes by cp and does one batched matmul per
        # unique coord system.
        from PySide6.QtWidgets import QApplication
        from node_runner.scene_build import (
            build_node_coords_vectorized,
            build_element_arrays_vectorized,
        )
        n_nodes = len(self.current_node_ids_sorted)
        n_elems_total = len(model.elements) + len(model.rigid_elements)
        log_deck_header(
            getattr(self, '_loaded_filepath', '') or '',
            n_nodes=n_nodes,
            n_elements=n_elems_total,
        )
        if n_nodes > 50_000:
            self._update_status(f"Loading {n_nodes:,} nodes (vectorized)...")
            QApplication.processEvents()
        with perf_stage('viewer', 'build_node_coords', n_nodes=n_nodes):
            node_coords = build_node_coords_vectorized(
                model, self.current_node_ids_sorted)
        if n_nodes > 50_000:
            QApplication.processEvents()

        # v3.2.0: bulk-build per-element-kind connectivity arrays via
        # the vectorized helper. The old per-element Python loop here
        # was the dominant time sink on big decks (~30 s on 600k
        # elements). The new helper sorts elements by kind, builds
        # uniform (n_elem, n_nodes) numpy arrays, and does the node-ID
        # -> grid-index lookup via np.searchsorted in a single
        # vectorized op.
        n_elems = len(model.elements) + len(model.rigid_elements)
        if n_elems > 50_000:
            self._update_status(
                f"Building scene from {n_elems:,} elements (vectorized)...")
            QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            QApplication.processEvents()

        def _scene_build_progress(done, total, kind):
            self._update_status(
                f"Building scene: {done:,} / {total:,} elements ({kind})...")
            QApplication.processEvents()

        try:
            with perf_stage('viewer', 'build_element_arrays', n_elems=n_elems):
                elem_arrays = build_element_arrays_vectorized(
                    model,
                    self.current_node_ids_sorted,
                    node_map,
                    progress=(_scene_build_progress if n_elems > 50_000 else None),
                )
        finally:
            if n_elems > 50_000:
                QApplication.restoreOverrideCursor()
                self._update_status(
                    f"Scene built: {n_elems:,} elements. Rendering...")
                QApplication.processEvents()

        nodes_used_by_elements = elem_arrays['nodes_used']

        # v3.3.0: stash the tree-group counts the helper accumulated
        # for us during its single pass over all elements. The tree
        # population code reads these via self._cached_element_group_counts
        # instead of running its own 1.2M-iteration Counter loop.
        self._cached_element_group_counts = {
            'by_type': dict(elem_arrays.get('by_type_counts') or {}),
            'by_shape': dict(elem_arrays.get('by_shape_counts') or {}),
        }

        # PLOTELs are rendered separately but their nodes should not appear as free nodes
        for eid, elem in model.plotels.items():
            try: nodes_used_by_elements.update(elem.nodes)
            except (AttributeError, KeyError): pass

        all_node_ids = set(model.nodes.keys())
        free_node_ids = all_node_ids - nodes_used_by_elements

        # Free-node VERTEX cells. Build via numpy directly.
        free_node_idx = np.array(
            [node_map[nid] for nid in free_node_ids if nid in node_map],
            dtype=np.int64,
        )
        n_free = free_node_idx.size

        # Compose the master cells / cell_types arrays.
        main_cells = elem_arrays['cells']
        main_cell_types = elem_arrays['cell_types']
        if n_free:
            vertex_cells = np.empty(n_free * 2, dtype=np.int64)
            vertex_cells[0::2] = 1  # one node per VERTEX cell
            vertex_cells[1::2] = free_node_idx
            vertex_types = np.full(n_free, pv.CellType.VERTEX, dtype=np.uint8)
            cells = np.concatenate([main_cells, vertex_cells])
            cell_types = np.concatenate([main_cell_types, vertex_types])
        else:
            cells = main_cells
            cell_types = main_cell_types

        with perf_stage('viewer', 'assemble_unstructured_grid',
                        n_cells=int(cell_types.size if cell_types is not None else 0)):
            if cells.size == 0 and node_coords.size > 0:
                # No elements, no free-node vertices yet - render all nodes
                # as vertex cells so something appears.
                num_nodes = node_coords.shape[0]
                cells = np.empty(num_nodes * 2, dtype=np.int64)
                cells[0::2] = 1
                cells[1::2] = np.arange(num_nodes, dtype=np.int64)
                cell_types = np.full(num_nodes, pv.CellType.VERTEX, dtype=np.uint8)
                self.current_grid = pv.UnstructuredGrid(cells, cell_types, node_coords)
            elif cells.size > 0:
                self.current_grid = pv.UnstructuredGrid(cells, cell_types, node_coords)
            else:
                self.current_grid = pv.UnstructuredGrid()

            self.current_grid.point_data['vtkOriginalPointIds'] = np.arange(self.current_grid.n_points)

        # Attach per-cell metadata. We had main element cells + free-node
        # vertex cells; pad vertex cells with -1 EID / -1 PID / VTK_VERTEX.
        if main_cells.size > 0 or n_free > 0:
            n_main = elem_arrays['eid'].size
            if n_free:
                eid_full = np.concatenate([
                    elem_arrays['eid'],
                    np.full(n_free, -1, dtype=np.int64)])
                pid_full = np.concatenate([
                    elem_arrays['pid'],
                    np.full(n_free, -1, dtype=np.int64)])
                etype_full = np.concatenate([
                    elem_arrays['etype'],
                    np.array(['VTK_VERTEX'] * n_free, dtype=object)])
                is_shell_full = np.concatenate([
                    elem_arrays['is_shell'],
                    np.zeros(n_free, dtype=np.int8)])
            else:
                eid_full = elem_arrays['eid']
                pid_full = elem_arrays['pid']
                etype_full = elem_arrays['etype']
                is_shell_full = elem_arrays['is_shell']
            self.current_grid.cell_data['EID'] = eid_full
            self.current_grid.cell_data['PID'] = pid_full
            self.current_grid.cell_data['type'] = etype_full
            self.current_grid.cell_data['is_shell'] = is_shell_full

        # v3.2.0: build the display-LOD'd grid. This is what
        # _rebuild_plot passes to plotter.add_mesh; current_grid keeps
        # every cell for picking and queries.
        try:
            from node_runner.scene_build import apply_display_lod
            with perf_stage('viewer', 'apply_display_lod',
                            status_msg='Building display LOD...',
                            status=True,
                            n_input_cells=int(self.current_grid.n_cells)):
                display_grid, lod_info = apply_display_lod(
                    self.current_grid,
                    threshold=self._elem_lod_threshold,
                    target_shells=200_000,
                )
            self.current_display_grid = display_grid
            self._lod_info = lod_info
            perf_event('viewer', 'lod_result',
                       active=lod_info.get('lod_active'),
                       displayed=lod_info.get('displayed_cells'),
                       total=lod_info.get('total_cells'))
            if lod_info.get('lod_active'):
                self._update_status(
                    f"LOD active: showing {lod_info['displayed_cells']:,} "
                    f"of {lod_info['total_cells']:,} cells (full mesh kept "
                    f"for picking)."
                )
        except Exception as _lod_exc:
            # Fail-safe: render the full grid.
            self.current_display_grid = self.current_grid
            self._lod_info = {'lod_active': False}
            perf_event('viewer', 'lod_exception', exc=str(_lod_exc)[:120])

        # v4.0.16 (reorder): populate pid_color_map BEFORE building
        # pre-built per-category actors. The build site reads
        # self.pid_color_map to compute its categorical cmap; if a
        # PID is missing, the cmap entry falls back to "#FFFFFF"
        # (white). Pre-v4.0.16 the build ran first → fresh-load
        # PIDs had no entries → user saw a mostly-white mesh on
        # first open. v4.0.15's _set_coloring_mode rebuild
        # accidentally hid this by rebuilding the actors AFTER the
        # map was populated; re-clicking "Color by Property ID"
        # appeared to "fix" the colors when really it was just
        # rebuilding with the now-populated map.
        for pid in model.properties.keys():
            if pid not in self.pid_color_map:
                self.pid_color_map[pid] = QColor(random.randint(50, 220), random.randint(50, 220), random.randint(50, 220)).name()

        # v4.0.14 (Fix C v2): build pre-built per-Nastran-type actors
        # from the FULL current_grid. Category toggles (Plates,
        # Beams, Solids, etc.) then become SetVisibility calls in
        # the fast-path router. PID/MID/isolate filtering still
        # routes through legacy _rebuild_plot via _has_active_cell_filter.
        try:
            self._build_pre_built_category_actors()
        except Exception as _exc:
            perf_event('actors', 'build_per_category_failed',
                       exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")

        self._emit_render_progress("Populating model tree...")
        with perf_stage('viewer', 'populate_tree',
                        status_msg='Populating model tree...', status=True):
            self.tree_widget.blockSignals(True)
            self._populate_tree()
            self.tree_widget.blockSignals(False)

        # Coord actors are created centrally via _refresh_coord_actors (called by _update_plot_visibility)
        # Each of these can iterate the full bulk (e.g. 268k PLOAD4s);
        # show a status hint between each so the user sees motion.
        self._emit_render_progress("Building load actors (FORCE/MOMENT/PLOAD4)...")
        with perf_stage('viewer', 'create_load_actors',
                        status_msg='Building load actors (FORCE/MOMENT/PLOAD4)...',
                        status=True,
                        n_load_sids=len(model.loads)):
            self._create_all_load_actors()
        self._emit_render_progress("Building constraint actors (SPC/SPC1)...")
        with perf_stage('viewer', 'create_constraint_actors',
                        status_msg='Building constraint actors (SPC/SPC1)...',
                        status=True,
                        n_spc_sids=len(model.spcs)):
            self._create_all_constraint_actors()
        self._emit_render_progress("Building rigid-element actors (RBE/RBAR)...")
        with perf_stage('viewer', 'create_rbe_actors',
                        status_msg='Building rigid-element actors (RBE/RBAR)...',
                        status=True,
                        n_rigid=len(model.rigid_elements)):
            self._create_rbe_actors()
        self._emit_render_progress("Building mass actors (CONM2)...")
        with perf_stage('viewer', 'create_mass_actors',
                        status_msg='Building mass actors (CONM2)...',
                        status=True,
                        n_mass=len(model.masses)):
            self._create_mass_actors()
        self._emit_render_progress("Building plot-element actors (PLOTEL)...")
        with perf_stage('viewer', 'create_plotel_actors',
                        status_msg='Building plot-element actors (PLOTEL)...',
                        status=True,
                        n_plotel=len(model.plotels)):
            self._create_plotel_actors()
        self._emit_render_progress("Final render...")
        with perf_stage('viewer', 'first_update_plot_visibility',
                        status_msg='Initial visibility update...', status=True):
            self._update_plot_visibility()

        with perf_stage('viewer', 'reset_camera_and_render'):
            if reset_camera and self.current_grid and self.current_grid.n_points > 0:
                self.plotter.reset_camera()
                self.plotter.view_isometric()
            else:
                self.plotter.render()
        perf_event('viewer', 'total', wall_s=f"{time.perf_counter()-_t_total:.3f}")
        self._emit_mem_event('viewer.total')
# END: Corrected replacement for _update_viewer method in main.py


    # v4.0.4 (Stage Q): UserRole-tuple key prefixes that belong to the
    # Loads tab tree (loads_tab_tree), NOT the Model tree (tree_widget).
    # `_find_tree_items` routes lookups for these prefixes to the
    # loads_tab index. Pre-v4.0.4 they all queried tree_widget and
    # always returned [] -> load actors were never visible.
    _LOADS_TAB_KEY_PREFIXES = frozenset({
        'load_set', 'constraint_set',
        'load_set_placeholder', 'constraint_set_placeholder',
        'load_entry', 'constraint_entry',
        'tempd_entry',
        'loads_root', 'constraints_root', 'combos_root',
        'load_combo',
    })

    def _find_tree_items(self, data_tuple):
        # v4.0.3 (Stage J) / v4.0.4 (Stage Q): O(1) lookup via the
        # pre-built tree index. Now ROUTES based on key prefix:
        # load_set / constraint_set / etc. → loads_tab index; the
        # rest → model_tree index. The pre-v4.0.4 path always
        # queried tree_widget regardless of key, so all load_set /
        # constraint_set lookups returned [] (the wrong-tree bug
        # the tree_index.miss logger exposed).
        self._find_tree_calls_total = getattr(self, '_find_tree_calls_total', 0) + 1
        # Normalize the lookup key the same way the index was built:
        # list/tuple inputs get hashed as tuples.
        try:
            key = tuple(data_tuple) if isinstance(data_tuple, (list, tuple)) else data_tuple
        except TypeError:
            key = data_tuple

        # v4.0.4 (Stage Q): route by prefix.
        prefix = key[0] if isinstance(key, tuple) and key else None
        if prefix in self._LOADS_TAB_KEY_PREFIXES:
            index = getattr(self, '_loads_tab_item_index', None)
            index_name = 'loads_tab_index'
            tree = getattr(self, 'loads_tab_tree', None)
        else:
            index = getattr(self, '_tree_item_index', None)
            index_name = 'tree_index'
            tree = getattr(self, 'tree_widget', None)

        if index is None:
            from node_runner.profiling import perf_event
            perf_event(index_name, 'fallback_walk',
                       reason='no_index_yet', key=str(key)[:80])
            if tree is None:
                return []
            return [it.value() for it in QTreeWidgetItemIterator(tree)
                    if it.value().data(0, QtCore.Qt.UserRole) == data_tuple]
        hits = index.get(key)
        if hits is None:
            # Real miss. key isn't in the index. Could indicate a
            # caller asking for an item that doesn't exist (normal)
            # OR a stale index (bug). Log a sample so we can spot
            # patterns.
            counter_attr = (
                '_loads_tab_index_misses' if prefix in self._LOADS_TAB_KEY_PREFIXES
                else '_tree_index_misses')
            cur = getattr(self, counter_attr, 0) + 1
            setattr(self, counter_attr, cur)
            # v4.0.9: gate elem_by_type_group / elem_shape_group misses
            # to types/shapes that ARE present in this deck. Callers
            # iterate a fixed list of known categories on every
            # visibility cycle, so misses for absent categories are
            # expected noise (e.g. Shear/Gap/Tet/Wedge on a deck that
            # only has plates+solids). Real surprises (stale index,
            # genuine bug) still emit through the `unexpected_miss`
            # path below.
            expected_absent = False
            counts = getattr(self, '_cached_element_group_counts', None) or {}
            if (isinstance(key, tuple) and len(key) == 2
                    and key[0] == 'elem_by_type_group'):
                expected_absent = key[1] not in (counts.get('by_type') or {})
            elif (isinstance(key, tuple) and len(key) == 2
                    and key[0] == 'elem_shape_group'):
                expected_absent = key[1] not in (counts.get('by_shape') or {})
            if cur <= 25:  # cap log spam
                from node_runner.profiling import perf_event
                if expected_absent:
                    # Quieter signal kept so we still see it if we
                    # need to debug "why didn't this match?".
                    perf_event(index_name, 'expected_miss',
                               key=str(key)[:80], total_misses=cur)
                else:
                    perf_event(index_name, 'miss', key=str(key)[:80],
                               total_misses=cur)
            return []
        # Filter out items that have been deleted (Qt can leave dangling
        # pointers in the index after partial tree mutations). Skip on
        # any access exception.
        live = []
        for it in hits:
            try:
                _ = it.text(0)  # access fails on a deleted item
                live.append(it)
            except RuntimeError:
                continue
        return live

    def _build_tree_item_index(self):
        """v4.0.3 (Stage J): walk the model tree once and bucket every
        item by its UserRole tuple. Called from ``_populate_tree``
        after the rebuild completes.

        Index shape: ``{key: [QTreeWidgetItem, ...]}``. Keys are
        normalized via ``tuple(...)`` so list/tuple equality works.
        Items with no UserRole are skipped (no way to look them up).
        """
        from node_runner.profiling import perf_event
        import time as _time
        _t0 = _time.perf_counter()
        index = {}
        n_entries = 0
        it = QTreeWidgetItemIterator(self.tree_widget)
        while it.value():
            item = it.value()
            try:
                data = item.data(0, QtCore.Qt.UserRole)
            except Exception:
                data = None
            if data is not None:
                try:
                    key = tuple(data) if isinstance(data, (list, tuple)) else data
                    index.setdefault(key, []).append(item)
                    n_entries += 1
                except TypeError:
                    pass  # unhashable; skip
            it += 1
        self._tree_item_index = index
        self._tree_index_misses = 0
        perf_event('tree_index', 'build',
                   n_entries=n_entries,
                   n_keys=len(index),
                   wall_s=f"{_time.perf_counter()-_t0:.3f}")

    def _build_loads_tab_item_index(self):
        """v4.0.4 (Stage Q): mirror of ``_build_tree_item_index`` but
        for the LOADS tab tree (``loads_tab_tree``). Required because
        all ``load_set`` / ``constraint_set`` UserRole keys live in
        that tree, not the model tree. pre-v4.0.4 we were querying
        the wrong tree and getting [] for every lookup.
        """
        from node_runner.profiling import perf_event
        import time as _time
        _t0 = _time.perf_counter()
        index = {}
        n_entries = 0
        tree = getattr(self, 'loads_tab_tree', None)
        if tree is not None:
            it = QTreeWidgetItemIterator(tree)
            while it.value():
                item = it.value()
                try:
                    data = item.data(0, QtCore.Qt.UserRole)
                except Exception:
                    data = None
                if data is not None:
                    try:
                        key = tuple(data) if isinstance(data, (list, tuple)) else data
                        index.setdefault(key, []).append(item)
                        n_entries += 1
                    except TypeError:
                        pass
                it += 1
        self._loads_tab_item_index = index
        self._loads_tab_index_misses = 0
        perf_event('loads_tab_index', 'build',
                   n_entries=n_entries,
                   n_keys=len(index),
                   wall_s=f"{_time.perf_counter()-_t0:.3f}")

    def _build_select_by_data(self):
        """Build property/material/type lookup dicts for element selection.

        v3.2.4: 'By Quality' is computed lazily on demand instead of
        up front. Element-quality metrics (aspect, skew, warp,
        Jacobian, taper) for a 700k-element deck are many minutes of
        pure-Python work; computing them every time the user opens
        any selection dialog made `List > Element Information` hang
        for 5+ minutes. The Method menu now shows 'By Quality
        (compute...)' which triggers the calc only if the user
        actually picks it.
        """
        model = self.current_generator.model
        all_elements = {**model.elements, **model.rigid_elements}
        by_property, by_material, by_type = {}, {}, {}
        for eid, elem in all_elements.items():
            by_type.setdefault(elem.type, []).append(eid)
            pid = getattr(elem, 'pid', None)
            if pid is not None:
                by_property.setdefault(pid, []).append(eid)
                prop = model.properties.get(pid)
                if prop:
                    mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None)
                    if mid: by_material.setdefault(mid, []).append(eid)

        result = {'By Property': by_property, 'By Material': by_material,
                  'By Type': by_type}
        # 'By Quality' is a deferred-compute placeholder. The selection
        # bar shows it as a method but only triggers the actual quality
        # calc when the user picks it. See _compute_quality_select_data.
        result['By Quality'] = '__deferred_quality__'
        return result

    def _compute_quality_select_data(self):
        """Lazy producer of the 'By Quality' selection categories.

        Called from EntitySelectionBar when the user picks 'By Quality'
        method. Wraps the slow calculate_element_quality() call in a
        wait cursor + status message so the user can see something is
        happening. Returns a dict {label: [eid, ...]}.
        """
        from PySide6.QtWidgets import QApplication
        try:
            self._update_status(
                "Computing element-quality metrics... (one-time, may "
                "take 30-60s on big decks)")
            QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            QApplication.processEvents()
            quality = self.current_generator.calculate_element_quality()
        finally:
            QApplication.restoreOverrideCursor()
        thresholds = {'Aspect > 5': ('aspect', 5.0, 'gt'),
                      'Skew > 45 deg': ('skew', 45.0, 'gt'),
                      'Skew > 60 deg': ('skew', 60.0, 'gt'),
                      'Warp > 5 deg': ('warp', 5.0, 'gt'),
                      'Warp > 10 deg': ('warp', 10.0, 'gt'),
                      'Jacobian < 0.6': ('jacobian', 0.6, 'lt'),
                      'Taper > 0.5': ('taper', 0.5, 'gt')}
        out = {}
        for label, (metric, thresh, op) in thresholds.items():
            eids = []
            for eid, metrics in quality.items():
                val = metrics.get(metric)
                if val is None:
                    continue
                if (op == 'gt' and val > thresh) or (op == 'lt' and val < thresh):
                    eids.append(eid)
            if eids:
                out[label] = eids
        self._update_status(
            f"Quality computed for {len(quality):,} elements.")
        return out

    def _open_property_editor(self):
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "No properties to edit."); return
        dialog = PropertyEditorDialog(self.current_generator.model, self)
        if dialog.exec(): self._update_plot_visibility(); self._populate_tree(); self._update_status("Properties updated.")
        
    def _open_material_editor(self):
        if not self.current_generator or not self.current_generator.model.materials: QMessageBox.warning(self, "No Materials", "No materials to edit."); return
        dialog = MaterialEditorDialog(self.current_generator.model, self)
        if dialog.exec(): self._update_status("Materials updated.")


    def _show_all_entities(self):
        """v4.0.0 (E3): re-check every checkable item in the Model tree
        and the Loads tab tree, then refresh visibility. The user's
        recovery path when they've hidden things and want to start
        clean without rebuilding the model.
        """
        from PySide6 import QtCore as _QtCore

        def _walk_check(parent):
            for i in range(parent.childCount()):
                c = parent.child(i)
                if c.flags() & _QtCore.Qt.ItemIsUserCheckable:
                    c.setCheckState(0, _QtCore.Qt.Checked)
                _walk_check(c)

        for tw_name in ('tree_widget', 'loads_tab_tree'):
            tw = getattr(self, tw_name, None)
            if tw is None:
                continue
            tw.blockSignals(True)
            try:
                _walk_check(tw.invisibleRootItem())
            finally:
                tw.blockSignals(False)

        # Also clear any hidden_groups flags so groups come back on.
        try:
            self._hidden_groups.clear()
            self._populate_groups_list()
        except Exception:
            pass

        try:
            self._update_plot_visibility()
        except Exception:
            pass
        self._update_status("Show All: every entity is visible.")

    def _dispatch_refresh(self, hint):
        """v4.0.0 (B2): map a ``RefreshHint`` to the matching narrow
        refresh helper. Callers can opt in:

            self.command_manager.execute(cmd, model)
            self._dispatch_refresh(cmd.refresh_hint)

        Falls back to ``_update_viewer`` for FULL (the pre-v4.0.0
        default) so unmigrated commands keep working unchanged.
        """
        from node_runner.commands import RefreshHint
        if hint == RefreshHint.NONE:
            return
        if hint == RefreshHint.CONSTRAINTS:
            self._refresh_constraints()
            return
        if hint == RefreshHint.LOADS:
            self._refresh_loads()
            return
        if hint == RefreshHint.GRID_COLOR:
            self._refresh_grid_color_data()
            return
        if hint == RefreshHint.COORDS:
            self._refresh_coord_actors()
            return
        if hint == RefreshHint.GROUPS:
            self._refresh_groups()
            return
        # FULL or any unknown hint → blanket refresh (pre-v4.0.0 default).
        self._update_viewer()

    def _refresh_constraints(self):
        """v4.0.0 (B2): narrow refresh after a constraint add/remove.
        Rebuilds only the constraint actors + the Loads-tab tree
        (preserving check state via B1 snapshot/restore). The main
        Model tree's "Boundary Conditions" count may show as stale by
        one until the next full _update_viewer; that's acceptable for
        the user-visible interactive flow.
        """
        try:
            self._populate_loads_tab()
        except Exception:
            pass
        try:
            self._create_all_constraint_actors()
        except Exception:
            pass
        try:
            self._update_plot_visibility()
        except Exception:
            pass

    def _refresh_loads(self):
        """v4.0.0 (B2): narrow refresh after a load add/remove."""
        try:
            self._populate_loads_tab()
        except Exception:
            pass
        try:
            self._create_all_load_actors()
        except Exception:
            pass
        try:
            self._update_plot_visibility()
        except Exception:
            pass

    def _refresh_grid_color_data(self):
        """v4.0.0 (B2): refresh grid coloring after material/property
        edits. Currently a thin wrapper around _update_plot_visibility
        which rebuilds the colored mesh actors; the grid geometry is
        unchanged so no full _update_viewer is needed.
        """
        try:
            self._update_plot_visibility()
        except Exception:
            pass

    def _refresh_coord_actors(self):
        """v4.0.0 (B2): refresh coord-system axis arrows after a
        coord add/edit/delete. Falls back to _update_viewer for now
        (axis arrows are created inline with _update_viewer). A future
        pass can extract a dedicated _create_all_coord_actors helper.
        """
        self._update_viewer()

    def _refresh_groups(self):
        """v4.0.0 (B2): refresh groups list + visibility after a
        group create/modify/rename. Skips the full model rebuild.
        """
        try:
            self._populate_groups_list()
        except Exception:
            pass
        try:
            self._update_plot_visibility()
        except Exception:
            pass

    @staticmethod
    def _snapshot_tree_state(tree, tree_id='?'):
        """v4.0.0 (B1) / v4.0.3 (Stage O): walk a QTreeWidget and
        capture per-item state keyed by the item's UserRole data tuple.

        v4.0.3 adds a diagnostic perf_event with the snapshot size
        and a few sample keys so we can verify the restore side picks
        the right matches (debug for the constraint-set uncheck bug).

        Returns ``{key: {'checked': CheckState, 'expanded': bool}}``.
        Keys come from ``item.data(0, UserRole)`` set when the item was
        created. Items without a UserRole key are skipped (no way to
        match them on rebuild).
        """
        from PySide6 import QtCore as _QtCore
        out = {}
        # v4.0.3 (Stage O): also capture the per-item state for any
        # constraint_set keys so we can correlate snapshot↔restore.
        constraint_snapshot = {}

        def walk(parent):
            n = parent.childCount()
            for i in range(n):
                child = parent.child(i)
                key = child.data(0, _QtCore.Qt.UserRole)
                if key is not None:
                    try:
                        hashable = tuple(key) if isinstance(key, (list, tuple)) else key
                        out[hashable] = {
                            'checked': child.checkState(0),
                            'expanded': child.isExpanded(),
                        }
                        if (isinstance(hashable, tuple)
                                and len(hashable) >= 1
                                and hashable[0] == 'constraint_set'):
                            # v4.0.5b: int(CheckState) fails on
                            # PySide6 strict enums. Use safe fallback
                            # so the snapshot diagnostic doesn't lose
                            # data via the outer TypeError except.
                            _cs = child.checkState(0)
                            try:
                                _cs_int = int(_cs)
                            except Exception:
                                try:
                                    _cs_int = int(_cs.value)
                                except Exception:
                                    _cs_int = repr(_cs)
                            constraint_snapshot[hashable] = _cs_int
                    except TypeError:
                        pass  # unhashable key; skip
                walk(child)

        root = tree.invisibleRootItem()
        walk(root)

        try:
            from node_runner.profiling import perf_event
            sample_keys = list(out.keys())[:8]
            perf_event('tree_snapshot', 'capture',
                       tree_id=tree_id,
                       n_items=len(out),
                       sample_keys=','.join(str(k)[:30] for k in sample_keys),
                       constraint_sets=str(constraint_snapshot)[:200])
        except Exception:
            pass
        return out

    @staticmethod
    def _restore_tree_state(tree, snapshot, tree_id='?'):
        """v4.0.0 (B1) / v4.0.3 (Stage O): apply a snapshot from
        ``_snapshot_tree_state`` back onto the rebuilt tree.

        Items added since snapshot use whatever default the build set;
        items removed since snapshot are silently dropped. Signals must
        be blocked by the caller during clear/rebuild. restore itself
        sets check state via ``setCheckState`` which can fire
        ``itemChanged``; callers that care should block.

        v4.0.3 adds match/miss counters so we can see whether the
        constraint_set checkbox uncheck symptom is a snapshot-side
        problem (missing key) or a restore-side problem (key present
        but ignored).
        """
        from PySide6 import QtCore as _QtCore
        if not snapshot:
            try:
                from node_runner.profiling import perf_event
                perf_event('tree_snapshot', 'restore_skip_empty',
                           tree_id=tree_id)
            except Exception:
                pass
            return

        # v4.0.3 (Stage O) counters.
        matched = 0
        skipped_not_in_snap = 0
        constraint_restore = {}

        # v4.0.5 (Stage W): keys we care about specifically for the
        # top-level-row checkbox regression diagnostic.
        TOP_ROW_KEYS = {
            ('loads_root',),
            ('constraints_root',),
            ('combos_root',),
        }

        def walk(parent):
            nonlocal matched, skipped_not_in_snap
            n = parent.childCount()
            for i in range(n):
                child = parent.child(i)
                key = child.data(0, _QtCore.Qt.UserRole)
                if key is not None:
                    try:
                        hashable = tuple(key) if isinstance(key, (list, tuple)) else key
                        saved = snapshot.get(hashable)
                        # v4.0.5 (Stage W): log top-level row state
                        # before/after restore touches them.
                        if hashable in TOP_ROW_KEYS:
                            try:
                                from node_runner.profiling import perf_event as _w2_perf_event
                                # v4.0.5b: int(ItemFlag) raises on
                                # PySide6 6.x strict enums. Fall back
                                # to .value / repr safely.
                                try:
                                    _flags_int = int(child.flags())
                                except Exception:
                                    try:
                                        _flags_int = int(child.flags().value)
                                    except Exception:
                                        _flags_int = repr(child.flags())
                                try:
                                    _checkable_before = bool(
                                        child.flags() & _QtCore.Qt.ItemIsUserCheckable)
                                except Exception:
                                    _checkable_before = None
                                try:
                                    _saved_checked = (int(saved['checked'])
                                                      if saved else None)
                                except Exception:
                                    _saved_checked = None
                                _w2_perf_event(
                                    'top_row_flags', 'restore_before',
                                    key=str(hashable),
                                    flags=_flags_int,
                                    checkable_before=_checkable_before,
                                    saved_present=(saved is not None),
                                    saved_checked=_saved_checked,
                                )
                            except Exception:
                                pass
                        if saved is not None:
                            if child.flags() & _QtCore.Qt.ItemIsUserCheckable:
                                child.setCheckState(0, saved['checked'])
                            child.setExpanded(saved['expanded'])
                            matched += 1
                            if (isinstance(hashable, tuple)
                                    and len(hashable) >= 1
                                    and hashable[0] == 'constraint_set'):
                                # v4.0.5b: safe int conversion for
                                # CheckState (strict-enum compat).
                                def _ck_int(v):
                                    try:
                                        return int(v)
                                    except Exception:
                                        try:
                                            return int(v.value)
                                        except Exception:
                                            return repr(v)
                                constraint_restore[hashable] = (
                                    _ck_int(saved['checked']),
                                    _ck_int(child.checkState(0)),
                                )
                        else:
                            skipped_not_in_snap += 1
                        if hashable in TOP_ROW_KEYS:
                            try:
                                from node_runner.profiling import perf_event as _w3_perf_event
                                # v4.0.5b: safe int conversion (see
                                # restore_before block above).
                                try:
                                    _flags_int = int(child.flags())
                                except Exception:
                                    try:
                                        _flags_int = int(child.flags().value)
                                    except Exception:
                                        _flags_int = repr(child.flags())
                                try:
                                    _checkable_after = bool(
                                        child.flags() & _QtCore.Qt.ItemIsUserCheckable)
                                except Exception:
                                    _checkable_after = None
                                _w3_perf_event(
                                    'top_row_flags', 'restore_after',
                                    key=str(hashable),
                                    flags=_flags_int,
                                    checkable_after=_checkable_after,
                                )
                            except Exception:
                                pass
                    except TypeError:
                        pass
                walk(child)

        root = tree.invisibleRootItem()
        walk(root)

        try:
            from node_runner.profiling import perf_event
            perf_event('tree_snapshot', 'restore_match',
                       tree_id=tree_id,
                       n_snapshot=len(snapshot),
                       matched=matched,
                       skipped_not_in_snap=skipped_not_in_snap,
                       missed=max(0, len(snapshot) - matched),
                       constraint_restore=str(constraint_restore)[:200])
        except Exception:
            pass

    def _populate_tree(self):
        """Populate the model tree.

        v3.3.0: speed up by (a) reading element By-Type / By-Shape
        counts from ``self._cached_element_group_counts`` (filled by
        ``_update_viewer`` from the scene-build pass) instead of running
        two more 600k-element Counter loops here, (b) wrapping the
        whole build with ``setUpdatesEnabled(False)`` so Qt does one
        layout pass at the end instead of one per insert, (c) building
        QTreeWidgetItem lists then calling ``addChildren()`` once per
        high-cardinality parent, (d) replacing ``expandAll()`` + collapse
        with explicit per-group expand so the coord group never gets
        the full layout treatment.

        v4.0.0 (B1): snapshot tree state before clear, restore after
        build. Without this, small edits (add SPC, edit material) that
        end with a tree refresh would wipe the user's checkbox
        visibility selections.
        """
        from node_runner.profiling import perf_stage
        # v4.0.3 (Stage J): invalidate the lookup index before the
        # rebuild so any in-flight `_find_tree_items` call falls back
        # to a legacy walk rather than reading stale pointers.
        self._tree_item_index = None
        with perf_stage('tree', 'snapshot_state'):
            snapshot = self._snapshot_tree_state(self.tree_widget, tree_id='model_tree')
        try:
            self.tree_widget.itemChanged.disconnect(self._handle_tree_item_changed)
        except (RuntimeError, TypeError): pass
        with perf_stage('tree', 'clear_widget'):
            self.tree_widget.clear()

        if not self.current_generator:
            self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)
            # v4.0.3 (Stage J): empty tree gets an empty index so
            # callers don't fall back to the legacy walk forever.
            self._tree_item_index = {}
            return

        model = self.current_generator.model

        # v3.3.0: batch all the work between setUpdatesEnabled toggles.
        self.tree_widget.setUpdatesEnabled(False)
        try:
            with perf_stage('tree', 'build_contents',
                            n_props=len(model.properties),
                            n_mats=len(model.materials)):
                self._build_tree_contents(model)
            with perf_stage('tree', 'restore_state'):
                self._restore_tree_state(self.tree_widget, snapshot, tree_id='model_tree')
            # v4.0.3 (Stage J): rebuild the lookup index now that the
            # tree contents are final.
            with perf_stage('tree', 'build_item_index'):
                self._build_tree_item_index()
        finally:
            self.tree_widget.setUpdatesEnabled(True)

        self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)
        with perf_stage('tree', 'populate_loads_tab'):
            self._populate_loads_tab()
        with perf_stage('tree', 'populate_analysis_tab'):
            self._populate_analysis_tab()

    def _build_tree_contents(self, model):
        """Actual tree-population body. Called with updates disabled."""
        def _make_item(text, data=None, checked=True):
            item = QTreeWidgetItem([text])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(
                0, QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            item.setData(0, QtCore.Qt.UserRole, data)
            return item

        def add_top(text, data=None, checked=True):
            item = _make_item(text, data, checked)
            self.tree_widget.addTopLevelItem(item)
            return item

        # ---- Coordinate Systems (always shown; v3.4.0 item 6) ----
        # Group defaults to UNCHECKED + COLLAPSED. Always created so
        # the user sees N=0 explicitly when there are no coord systems
        # (rather than the branch silently disappearing).
        n_coords = len(model.coords or {})
        coords_top = add_top(
            f"Coordinate Systems ({n_coords})",
            data=('group', 'coords'), checked=False)
        if model.coords:
            coord_children = []
            for cid in sorted(model.coords.keys()):
                if cid == 0:
                    label = "Global Rectangular"
                elif cid == 1:
                    label = "Global Cylindrical"
                elif cid == 2:
                    label = "Global Spherical"
                else:
                    coord = model.coords[cid]
                    type_map = {'CORD2R': 'Rectangular',
                                'CORD2C': 'Cylindrical',
                                'CORD2S': 'Spherical'}
                    type_name = type_map.get(coord.type, coord.type)
                    title = get_entity_title_from_comment(
                        getattr(coord, 'comment', ''), "Coord", cid)
                    label = (f"{type_name} - {title}"
                             if title != f"Coord {cid}" else type_name)
                coord_children.append(_make_item(
                    f"CID {cid}: {label}",
                    data=('coord', cid), checked=False))
            coords_top.addChildren(coord_children)

        # ---- Materials (always shown; v3.4.0 item 6) ----
        n_materials = len(model.materials or {})
        mats_item = add_top(
            f"Materials ({n_materials})", data=('group', 'materials'))
        if model.materials:
            mat_children = [
                _make_item(
                    f"{mid}: {self._get_entity_title_from_comment(m.comment, 'Material', mid)}",
                    data=('material', mid))
                for mid, m in sorted(model.materials.items())
            ]
            mats_item.addChildren(mat_children)
            mats_item.setExpanded(True)

        # ---- Nodes header ----
        nodes_item = add_top(
            f"Nodes ({len(model.nodes)})", data=('group', 'nodes'))
        nodes_item.setExpanded(True)

        # ---- Elements (use cached counts) ----
        all_elements_n = len(model.elements) + len(model.rigid_elements)
        if all_elements_n:
            elems_item = add_top(
                f"Elements ({all_elements_n})", data=('group', 'elements'))
            type_item = _make_item("By Type", data=('group', 'elem_by_type'))
            shape_item = _make_item("By Shape", data=('group', 'elem_by_shape'))
            elems_item.addChildren([type_item, shape_item])

            counts = getattr(self, '_cached_element_group_counts', None) or {}
            by_type_counts = counts.get('by_type') or {}
            by_shape_counts = counts.get('by_shape') or {}
            type_children = [
                _make_item(f"{name} ({count})",
                           data=('elem_by_type_group', name))
                for name, count in sorted(by_type_counts.items())
            ]
            if type_children:
                type_item.addChildren(type_children)
            shape_children = [
                _make_item(f"{shape} ({count})",
                           data=('elem_shape_group', shape))
                for shape, count in sorted(by_shape_counts.items())
            ]
            if shape_children:
                shape_item.addChildren(shape_children)
            elems_item.setExpanded(True)
            type_item.setExpanded(True)
            shape_item.setExpanded(True)

        # ---- Masses (batch insert) ----
        if model.masses:
            masses_item = add_top(
                f"Masses ({len(model.masses)})", data=('group', 'masses'))
            mass_children = [
                _make_item(
                    f"CONM2 {eid}: Node {mass_elem.nid}, M={mass_elem.mass}",
                    data=('mass', eid))
                for eid, mass_elem in sorted(model.masses.items())
            ]
            masses_item.addChildren(mass_children)
            masses_item.setExpanded(True)

        # ---- Plot Elements (batch insert) ----
        if model.plotels:
            plotels_item = add_top(
                f"Plot Elements ({len(model.plotels)})",
                data=('group', 'plotels'))
            plotel_children = [
                _make_item(
                    f"PLOTEL {eid}: Nodes {elem.nodes[0]}-{elem.nodes[1]}",
                    data=('plotel', eid))
                for eid, elem in sorted(model.plotels.items())
            ]
            plotels_item.addChildren(plotel_children)
            plotels_item.setExpanded(True)

        # ---- Properties (always shown; v3.4.0 item 6) ----
        n_props = len(model.properties or {})
        props_item = add_top(
            f"Properties ({n_props})", data=('group', 'properties'))
        if model.properties:
            prop_children = [
                _make_item(
                    f"{pid}: {self._get_entity_title_from_comment(p.comment, p.type, pid)}",
                    data=('property', pid))
                for pid, p in sorted(model.properties.items())
            ]
            props_item.addChildren(prop_children)
            props_item.setExpanded(True)

        # ---- Geometry tree section ----
        if not self.geometry_store.is_empty:
            geom_item = add_top("Geometry", data=('group', 'geometry'))
            geom_subitems = []
            if self.geometry_store.points:
                pts_item = _make_item(
                    f"Points ({len(self.geometry_store.points)})",
                    data=('geom_group', 'points'))
                pt_children = [
                    _make_item(
                        f"P{pid}: ({pt.xyz[0]:.4g}, "
                        f"{pt.xyz[1]:.4g}, {pt.xyz[2]:.4g})",
                        data=('geom_point', pid))
                    for pid, pt in sorted(self.geometry_store.points.items())
                ]
                pts_item.addChildren(pt_children)
                geom_subitems.append(pts_item)
            all_curves = self.geometry_store.all_curve_ids()
            if all_curves:
                curves_item = _make_item(
                    f"Curves ({len(all_curves)})",
                    data=('geom_group', 'curves'))
                curve_children = []
                for cid in all_curves:
                    curve = self.geometry_store.get_curve(cid)
                    if (hasattr(curve, 'start_point_id')
                            and hasattr(curve, 'end_point_id')
                            and not hasattr(curve, 'mid_point_id')):
                        label = (f"Line {cid}: P{curve.start_point_id} -> "
                                 f"P{curve.end_point_id}")
                    elif hasattr(curve, 'mid_point_id'):
                        label = (f"Arc {cid}: P{curve.start_point_id} -> "
                                 f"P{curve.mid_point_id} -> "
                                 f"P{curve.end_point_id}")
                    elif hasattr(curve, 'radius'):
                        label = (f"Circle {cid}: Center P{curve.center_point_id}, "
                                 f"R={curve.radius:.4g}")
                    else:
                        label = f"Curve {cid}"
                    curve_children.append(_make_item(
                        label, data=('geom_curve', cid)))
                curves_item.addChildren(curve_children)
                geom_subitems.append(curves_item)
            if self.geometry_store.surfaces:
                surfs_item = _make_item(
                    f"Surfaces ({len(self.geometry_store.surfaces)})",
                    data=('geom_group', 'surfaces'))
                surf_children = [
                    _make_item(
                        f"Surface {sid}: {len(surf.boundary_curve_ids)} curves",
                        data=('geom_surface', sid))
                    for sid, surf in sorted(self.geometry_store.surfaces.items())
                ]
                surfs_item.addChildren(surf_children)
                geom_subitems.append(surfs_item)
            if geom_subitems:
                geom_item.addChildren(geom_subitems)
                geom_item.setExpanded(True)

        # v3.3.0: explicit per-group expand rather than expandAll()
        # then collapse-just-coords (which forced Qt to compute layout
        # for every nested branch only to undo it). Coord group stays
        # collapsed - the only group we never want auto-expanded.
        if coords_top is not None:
            coords_top.setExpanded(False)


    def _emit_render_progress(self, message: str) -> None:
        """v3.2.4: push a status update into the import dialog (if one
        is still open) during the post-parse scene-build steps.

        The import dialog stays visible through ``on_success`` so the
        user sees progress while we build PyVista actors, populate
        the tree, and rebuild the plot. Without this, the dialog
        would have closed at 'Imported N nodes' and the user would
        see a several-second freeze while the scene assembled.
        """
        dialog = getattr(self, '_active_import_dialog', None)
        if dialog is None:
            self._update_status(message)
            return
        try:
            from node_runner.model import ImportProgress
            dialog.apply_progress(ImportProgress(
                stage='render',
                label='Stage 6/6: Building 3D scene',
                detail=message,
                fraction=1.0,
            ))
            from PySide6.QtWidgets import QApplication
            QApplication.processEvents()
        except Exception:
            pass

    # ----- v3.2.4 tree search bar -----

    def _on_tree_search_changed(self, text: str):
        """Search box content changed.

        Filter mode: walk the tree, hide every item whose display
        text doesn't contain the query (case-insensitive), and unhide
        any item on the path from a matching child up to root.

        Find mode: just resets the find cursor; actual scroll happens
        on Enter via _on_tree_search_next.
        """
        if self._tree_search_filter_radio.isChecked():
            self._apply_tree_filter(text.strip().lower())
        else:
            # Find mode - reset the cursor so the next Enter starts
            # from the top.
            self._tree_search_find_cursor = 0

    def _on_tree_search_next(self):
        """Enter pressed in the search box.

        Filter mode: re-applies the filter (no-op if text unchanged).
        Find mode: scroll the next matching item into view + select.
        Wraps to the top after the last match.
        """
        text = self._tree_search.text().strip().lower()
        if not text:
            return
        if self._tree_search_filter_radio.isChecked():
            self._apply_tree_filter(text)
            return
        # Find mode: walk the tree, find the next match after the cursor.
        matches = []
        it = QTreeWidgetItemIterator(self.tree_widget)
        while it.value():
            item = it.value()
            label = item.text(0).lower()
            if text in label:
                matches.append(item)
            it += 1
        if not matches:
            self._update_status(f"Tree search: no match for '{text}'.")
            return
        idx = self._tree_search_find_cursor % len(matches)
        target = matches[idx]
        # Expand ancestors so the item is actually visible.
        parent = target.parent()
        while parent is not None:
            parent.setExpanded(True)
            parent = parent.parent()
        self.tree_widget.setCurrentItem(target)
        self.tree_widget.scrollToItem(target)
        self._tree_search_find_cursor = idx + 1
        self._update_status(
            f"Tree search: match {idx + 1} / {len(matches)} for '{text}'.")

    def _tree_search_clear_filter(self):
        """Switching back to Find mode - unhide everything."""
        it = QTreeWidgetItemIterator(self.tree_widget)
        while it.value():
            it.value().setHidden(False)
            it += 1

    def _apply_tree_filter(self, query: str):
        """Hide tree items that don't contain ``query`` (lowercase).

        Two-pass:
          1) Walk every item; mark "matches" if its own text contains
             query OR any descendant's text does.
          2) Hide items that are NOT marked.
        Empty query unhides everything.
        """
        if not query:
            self._tree_search_clear_filter()
            return

        # Collect all items so we can do a child-then-parent pass.
        all_items = []
        it = QTreeWidgetItemIterator(self.tree_widget)
        while it.value():
            all_items.append(it.value())
            it += 1

        # Pass 1: each item is "matched" if its own text matches OR
        # any descendant did. Process bottom-up so descendants are
        # resolved first.
        matched = {}
        for item in reversed(all_items):
            self_match = query in item.text(0).lower()
            any_child = False
            for k in range(item.childCount()):
                if matched.get(id(item.child(k)), False):
                    any_child = True
                    break
            matched[id(item)] = self_match or any_child

        # Pass 2: hide unmatched items; expand parents of matches.
        for item in all_items:
            is_match = matched[id(item)]
            item.setHidden(not is_match)
            if is_match and item.childCount():
                # Auto-expand so matches are visible.
                item.setExpanded(True)

    def _handle_tree_item_clicked(self, item, column):
        """v4.0.2 (Stage E): fires on every tree click regardless of
        check-state delta. Two roles:

        1. Diagnostic: emit a perf_event so the profile log records
           that the click reached the tree, even if itemChanged
           never fires (e.g. user clicked the row text not the
           checkbox).
        2. UX: if the user clicked on a category row with a
           checkbox (anything with UserRole data starting with
           'elem_by_type_group' / 'elem_shape_group' / 'group'),
           auto-toggle the checkbox. Cuts the "miss-the-checkbox"
           frustration.
        """
        from node_runner.profiling import perf_event
        try:
            data = item.data(0, QtCore.Qt.UserRole)
            label = item.text(0)
            cur_state = item.checkState(0)
        except Exception:
            return
        # v4.0.5b: safe int conversion (PySide6 strict-enum compat).
        try:
            _sb = int(cur_state)
        except Exception:
            try:
                _sb = int(cur_state.value)
            except Exception:
                _sb = repr(cur_state)
        perf_event('tree', 'item_clicked',
                   label=label,
                   data=str(data)[:80],
                   state_before=_sb,
                   checkable=bool(item.flags() & QtCore.Qt.ItemIsUserCheckable))
        # v4.0.6 (auto-toggle disabled): the v4.0.2 Stage E auto-toggle
        # caused a double-fire when the user clicked the checkbox
        # itself. Qt emits itemChanged (toggling once), THEN emits
        # itemClicked which entered this branch and toggled again,
        # leaving the visual state unchanged but firing TWO visibility
        # updates. The fix is to leave the click alone: itemChanged
        # fires natively on checkbox clicks (now that v4.0.5 removed
        # the viewport eventFilter that was blocking it). The
        # "click-on-label-doesn't-toggle" UX issue is acceptable for
        # now. the user can click the checkbox directly. Diagnostic
        # logging stays so we can still see clicks in the log.

    def _handle_tree_item_pressed(self, item, column):
        """v4.0.4 (Stage R): diagnostic-only handler for `itemPressed`.
        Fires before `itemClicked` for any mouse press on an item.
        If this fires but `itemClicked`/`itemChanged` don't, the bug
        is later in the Qt signal pipeline.
        """
        try:
            from node_runner.profiling import perf_event
            data = item.data(0, QtCore.Qt.UserRole) if item else None
            label = item.text(0) if item else ''
            perf_event('tree', 'item_pressed',
                       label=label,
                       data=str(data)[:80],
                       column=column)
        except Exception:
            pass

    def _handle_tree_selection_changed(self):
        """v4.0.4 (Stage R): diagnostic-only handler for
        `itemSelectionChanged`. Fires every time the tree selection
        changes (independent of itemClicked/itemChanged). Gives us a
        third probe point in the signal pipeline.
        """
        try:
            from node_runner.profiling import perf_event
            selected = self.tree_widget.selectedItems()
            sample = selected[0] if selected else None
            label = sample.text(0) if sample else ''
            data = sample.data(0, QtCore.Qt.UserRole) if sample else None
            perf_event('tree', 'selection_changed',
                       n_selected=len(selected),
                       sample_label=label,
                       sample_data=str(data)[:80])
        except Exception:
            pass

    def _handle_tree_item_changed(self, item, column):
        # v4.0.4 (Stage R): stdout sentinel. fires unconditionally so
        # we know whether Qt is dispatching the slot at all, even when
        # NR_PROFILE is unset. Visible in dev runs via `python -u run.py`
        # or in a console-attached build.
        print('[slot-fire] item_changed entered', flush=True)
        from node_runner.profiling import perf_stage, perf_event
        # Identify what the user clicked so the perf log is readable.
        try:
            data = item.data(0, QtCore.Qt.UserRole)
            label = item.text(0)
        except Exception:
            data, label = None, ''
        # v4.0.5b: int(CheckState) raises on PySide6 strict enums.
        try:
            _state = int(item.checkState(0))
        except Exception:
            try:
                _state = int(item.checkState(0).value)
            except Exception:
                _state = repr(item.checkState(0))
        perf_event('tree', 'item_changed', label=label,
                   data=str(data)[:60],
                   state=_state)
        with perf_stage('tree', 'handle_item_changed',
                        status_msg=f'Updating visibility for: {label}…',
                        status=True, label=label):
            # v4.0.10: source-tree-aware blockSignals. Pre-v4.0.10 this
            # only ever blocked the model tree, but with the loads_tab
            # connection now wired the same handler runs for both trees.
            # Block the tree that actually fired so the cascade
            # setCheckState calls below don't re-enter this handler.
            try:
                src_tree = item.treeWidget()
            except Exception:
                src_tree = None
            if src_tree is None:
                src_tree = self.tree_widget
            src_tree.blockSignals(True)
            try:
                def set_child_states(parent_item):
                    if parent_item.checkState(0) == QtCore.Qt.PartiallyChecked: return
                    for i in range(parent_item.childCount()):
                        child = parent_item.child(i); child.setCheckState(0, parent_item.checkState(0)); set_child_states(child)
                set_child_states(item)
                parent = item.parent()
                while parent:
                    child_states = [parent.child(i).checkState(0) for i in range(parent.childCount())]
                    if all(s == QtCore.Qt.Checked for s in child_states): parent.setCheckState(0, QtCore.Qt.Checked)
                    elif all(s == QtCore.Qt.Unchecked for s in child_states): parent.setCheckState(0, QtCore.Qt.Unchecked)
                    else: parent.setCheckState(0, QtCore.Qt.PartiallyChecked)
                    parent = parent.parent()
            finally:
                src_tree.blockSignals(False)
            # v4.0.9: fast path for toggles that affect a single named
            # actor. The full _update_plot_visibility cycle costs ~4-5 s
            # on a large production deck; if the click only toggles one independent actor
            # we skip the rebuild_plot / shells_compute_normals / SID
            # loops entirely. The router returns False on any
            # unexpected condition so we fall through safely.
            try:
                if self._route_fast_path_visibility(data, item.checkState(0)):
                    return
            except Exception as _exc:
                perf_event('fast_path', 'router_failed',
                           data=str(data)[:60],
                           exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            self._update_plot_visibility()

    def _route_fast_path_visibility(self, data, check_state):
        """v4.0.9: dispatch a tree-checkbox change to the right
        targeted handler when one exists. Returns True if the handler
        completed; False to let the caller fall through to
        _update_plot_visibility.

        Covered cases (each affects exactly one named actor):
          - ('load_set', sid)           → fast_path_sid_toggle
          - ('constraint_set', sid)     → fast_path_sid_toggle
          - ('elem_by_type_group', 'Rigid')  → rbe_actors.SetVisibility
          - ('group', 'masses')         → mass_actors.SetVisibility
          - ('group', 'plotels')        → plotel_actors.SetVisibility
          - ('group', 'coords')         → _refresh_coord_actors()
          - ('coord', cid)              → _refresh_coord_actors()
        """
        from node_runner.profiling import perf_event
        if not isinstance(data, tuple) or not data:
            return False
        is_checked = (check_state == QtCore.Qt.Checked)
        # Per-SID fast paths (load + constraint glyphs).
        if len(data) == 2 and data[0] in ('load_set', 'constraint_set'):
            return self._fast_path_sid_toggle(data[0], data[1], check_state)
        # Single-actor visibility toggles.
        # v4.0.14: the Rigid special case is again subsumed by
        # _fast_path_category_toggle below. v4.0.11 had this bug
        # because the LOD path dropped RBE2/RBE3; v4.0.14 builds
        # pre-built actors from current_grid (full mesh) so RBE2/
        # RBE3 ARE included. Plus the fast path now derives
        # rbe_actors visibility from the Rigid checkbox state
        # directly (not via target_states.get), so the v4.0.11
        # failure mode is defensively guarded against.
        if data == ('group', 'masses'):
            actor = self.plotter.actors.get('mass_actors')
            if actor is None:
                return False
            actor.SetVisibility(is_checked)
            perf_event('fast_path', 'masses_toggle',
                       checked=is_checked)
            self.plotter.render()
            return True
        if data == ('group', 'plotels'):
            actor = self.plotter.actors.get('plotel_actors')
            if actor is None:
                return False
            actor.SetVisibility(is_checked)
            perf_event('fast_path', 'plotels_toggle',
                       checked=is_checked)
            self.plotter.render()
            return True
        # Coord-system toggles only affect the coord arrows + labels,
        # nothing else in the scene. Refresh the coord actors directly.
        if data == ('group', 'coords') or (
                len(data) == 2 and data[0] == 'coord'):
            self._refresh_coord_actors(render=True)
            perf_event('fast_path', 'coords_toggle',
                       data=str(data)[:40])
            return True
        # v4.0.14 (Fix C v2): category toggles map to pre-built
        # per-Nastran-type actors built from full current_grid. The
        # handler returns False when isolate/hidden filter is active,
        # falling through to legacy _rebuild_plot.
        if (len(data) == 2
                and data[0] in ('elem_by_type_group', 'elem_shape_group')):
            return self._fast_path_category_toggle(data[0], data[1])
        # v5.0.0 item 11a: PID and MID toggles flip per-cell vtkGhostType
        # in-place on the pre-built actors. Sub-100 ms on a large production deck.
        if len(data) == 2 and data[0] == 'property':
            return self._fast_path_pid_toggle(data[1], is_checked)
        if len(data) == 2 and data[0] == 'material':
            return self._fast_path_mid_toggle(data[1], is_checked)
        return False

    # --- Loads Tab Methods ---

    def _format_load_entry(self, card) -> str:
        """Format a load card for display in the loads tree."""
        try:
            card_type = card.type
            if card_type in ('FORCE', 'FORCE1', 'FORCE2'):
                return f"FORCE  Node {card.node_id}  F={card.mag:.4g}  [{card.xyz[0]:.3g}, {card.xyz[1]:.3g}, {card.xyz[2]:.3g}]"
            elif card_type in ('MOMENT', 'MOMENT1', 'MOMENT2'):
                return f"MOMENT  Node {card.node_id}  M={card.mag:.4g}  [{card.xyz[0]:.3g}, {card.xyz[1]:.3g}, {card.xyz[2]:.3g}]"
            elif card_type == 'PLOAD4':
                return f"PLOAD4  Elem {card.eids[0]}  P={card.pressures[0]:.4g}"
            elif card_type == 'GRAV':
                return f"GRAV  Scale={card.scale:.4g}  [{card.N[0]:.3g}, {card.N[1]:.3g}, {card.N[2]:.3g}]"
            elif card_type == 'TEMPD':
                return f"TEMPD  T={card.temperature:.4g}"
            else:
                return f"{card_type}"
        except Exception:
            return str(card.type) if hasattr(card, 'type') else "Unknown"

    def _format_spc_entry(self, spc_card, node_id) -> str:
        """Format an SPC entry for a specific node.

        v4.0.6: SPC1 cards store a SINGLE components string that
        applies to ALL nodes in the card (Nastran semantics). The
        v3.5.0 code indexed ``spc_card.components`` by node index,
        which only happened to work when components had as many chars
        as nodes (e.g. a 3-node card with ``components='123'`` would
        show TX/TY/TZ for nodes 0/1/2). For a card with 7 nodes
        and ``components='123'``, nodes 3-6 raised IndexError and
        showed without a DOF label.

        Correct: use the whole components string for every node.
        """
        try:
            dof_names = {1: 'TX', 2: 'TY', 3: 'TZ', 4: 'RX', 5: 'RY', 6: 'RZ'}
            comp_raw = getattr(spc_card, 'components', None)
            if comp_raw is None:
                return f"Node {node_id}"
            # SPC1: components is typically a single str shared by all
            # nodes. SPC: components may be per-node (list/tuple). Try
            # per-node indexing only if components looks indexable AND
            # has at least as many entries as the node count.
            if isinstance(comp_raw, (list, tuple)):
                try:
                    idx = list(spc_card.nodes).index(node_id)
                    comp_for_node = comp_raw[idx]
                except (ValueError, IndexError):
                    comp_for_node = comp_raw[0] if comp_raw else ''
            else:
                # Single string applies to all nodes.
                comp_for_node = comp_raw
            comp_str = str(comp_for_node)
            dofs = " ".join(dof_names.get(int(c), c) for c in comp_str if c.isdigit())
            return f"Node {node_id}: {dofs}" if dofs else f"Node {node_id}"
        except Exception:
            return f"Node {node_id}"

    def _populate_loads_tab(self):
        """Populate the Loads sidebar tab tree.

        v3.5.0 item 4 (perf): lazy child population. Each Load Set /
        Constraint Set header gets a placeholder child '(N entries -
        expand to load)' and the real entries are built on first
        expand via _on_load_tab_item_expanded. On a large production
        deck (512k load entries) this drops Loads-tab populate from
        5-30s to <100ms.

        v4.0.0 (B1): preserve per-item check + expansion state across
        the clear/rebuild so unrelated edits don't wipe the user's
        visibility selections.
        """
        tree = self.loads_tab_tree
        # v4.0.4 (Stage Q): invalidate the loads-tab index before
        # the rebuild so any in-flight `_find_tree_items` lookup
        # falls back to the legacy walk (correct, just slower) rather
        # than reading stale pointers.
        self._loads_tab_item_index = None
        snapshot = self._snapshot_tree_state(tree, tree_id='loads_tab')
        # Disconnect itemExpanded to avoid firing during the clear/rebuild
        try:
            tree.itemExpanded.disconnect(self._on_load_tab_item_expanded)
        except (RuntimeError, TypeError):
            pass
        # v4.0.10: also block itemChanged. Since v4.0.10 wired this
        # signal at __init__ time, setCheckState calls inside the
        # rebuild (especially _restore_tree_state) would fire a
        # _handle_tree_item_changed cascade for every restored
        # checkbox state. blockSignals on the tree suppresses these
        # without disturbing the connection itself.
        tree.blockSignals(True)
        tree.setUpdatesEnabled(False)
        try:
            tree.clear()
            self._build_loads_tab_contents()
            self._restore_tree_state(tree, snapshot, tree_id='loads_tab')
            # v4.0.4 (Stage Q): rebuild the loads-tab lookup index
            # now that the tree contents are final.
            self._build_loads_tab_item_index()
        finally:
            tree.setUpdatesEnabled(True)
            tree.blockSignals(False)

        # Reconnect lazy-expand hook after the tree is rebuilt.
        tree.itemExpanded.connect(self._on_load_tab_item_expanded)
        # v4.0.10: snapshot restore re-set checkstates while signals
        # were blocked, so the fast-path-maintained sets are stale.
        # Refresh from tree truth now.
        try:
            self._refresh_visible_sid_sets()
        except Exception:
            pass

    def _build_loads_tab_contents(self):
        """Actual tree-population body. Called with updates disabled."""
        from node_runner.dialogs.load_combination import read_combo_payload
        tree = self.loads_tab_tree
        gen = getattr(self, 'current_generator', None)
        model = gen.model if gen else None
        if not model:
            return

        # --- Load Sets (combine model.loads and model.tempds) ---
        all_load_sids = set(model.loads.keys()) | set(model.tempds.keys())
        # v4.0.5 (Stage W): log flags of top-level rows so we can see
        # which step (if any) adds ItemIsUserCheckable to them.
        from node_runner.profiling import perf_event as _w_perf_event

        # v4.0.5b (hotfix): PySide6 6.x with strict enums raises
        # "int() argument must be ... not 'ItemFlag'" when you call
        # int(item.flags()) directly. Use this safe helper so the
        # diagnostic doesn't crash the post-parse handler.
        def _flag_int(flags):
            try:
                return int(flags)
            except Exception:
                try:
                    return int(flags.value)
                except Exception:
                    try:
                        return repr(flags)
                    except Exception:
                        return '?'

        def _is_checkable(flags):
            try:
                return bool(flags & QtCore.Qt.ItemIsUserCheckable)
            except Exception:
                # As a fallback, render the bitwise-AND string.
                try:
                    return ('ItemIsUserCheckable' in repr(flags))
                except Exception:
                    return None

        if all_load_sids:
            n_loads = sum(len(v) for v in model.loads.values()) + len(model.tempds)
            loads_root = QTreeWidgetItem(tree, [f"Load Sets ({n_loads} entries)"])
            loads_root.setData(0, QtCore.Qt.UserRole, ('loads_root',))
            # v4.0.9: PySide6 6.x default flags include ItemIsUserCheckable
            # (flags=61). The header rows shouldn't be checkable. clear the
            # bit so they look like proper section headers.
            loads_root.setFlags(loads_root.flags() & ~QtCore.Qt.ItemIsUserCheckable)
            _w_perf_event('top_row_flags', 'after_loads_root_create',
                          flags=_flag_int(loads_root.flags()),
                          checkable=_is_checkable(loads_root.flags()))
            # Build all Load Set headers without children first (batched).
            set_items = []
            for sid in sorted(all_load_sids):
                load_list = model.loads.get(sid, [])
                tempd_card = model.tempds.get(sid, None)
                counts = {}
                for load in load_list:
                    counts[load.type] = counts.get(load.type, 0) + 1
                if tempd_card:
                    counts['TEMPD'] = counts.get('TEMPD', 0) + 1
                summary = ", ".join(f"{c} {t}" for t, c in counts.items())
                n_children = len(load_list) + (1 if tempd_card else 0)
                set_item = QTreeWidgetItem(
                    [f"SID {sid}: {summary}"])
                set_item.setFlags(
                    set_item.flags() | QtCore.Qt.ItemIsUserCheckable)
                # v4.0.5b: default OFF. Pre-Stage Q (v4.0.4), the
                # visibility lookup queried the wrong tree and returned
                # [] for every load_set, so every SID was effectively
                # invisible. Stage Q routed the lookup to the correct
                # tree, exposing the fact that the build-time default
                # was Checked. which, on a large production deck (1,962 SIDs / 512k
                # entries), made the post-parse pass build glyph
                # geometry for every load entry. Default Unchecked
                # matches the user's expected UX ("loads/constraints
                # off unless explicitly enabled") AND keeps post-parse
                # render fast.
                set_item.setCheckState(0, QtCore.Qt.Unchecked)
                set_item.setData(0, QtCore.Qt.UserRole, ('load_set', sid))
                # v3.5.0: placeholder child so the disclosure triangle
                # appears. Real entries are built on first expand.
                if n_children > 0:
                    placeholder = QTreeWidgetItem(
                        [f"({n_children} entries - expand to load)"])
                    placeholder.setData(
                        0, QtCore.Qt.UserRole, ('load_set_placeholder', sid))
                    set_item.addChild(placeholder)
                set_items.append(set_item)
            loads_root.addChildren(set_items)
            loads_root.setExpanded(True)
            _w_perf_event('top_row_flags', 'after_loads_root_addchildren',
                          flags=_flag_int(loads_root.flags()),
                          checkable=_is_checkable(loads_root.flags()))

        # --- Constraint Sets (same lazy pattern) ---
        if model.spcs:
            n_constr = sum(
                len(s.nodes) for sl in model.spcs.values() for s in sl)
            constr_root = QTreeWidgetItem(
                tree, [f"Constraint Sets ({n_constr} nodes)"])
            constr_root.setData(0, QtCore.Qt.UserRole, ('constraints_root',))
            # v4.0.9: clear ItemIsUserCheckable from header row (see loads_root).
            constr_root.setFlags(constr_root.flags() & ~QtCore.Qt.ItemIsUserCheckable)
            _w_perf_event('top_row_flags', 'after_constr_root_create',
                          flags=_flag_int(constr_root.flags()),
                          checkable=_is_checkable(constr_root.flags()))
            constr_items = []
            for sid, spc_list in sorted(model.spcs.items()):
                all_nodes = set()
                for spc in spc_list:
                    all_nodes.update(spc.nodes)
                set_item = QTreeWidgetItem(
                    [f"SID {sid}: SPC ({len(all_nodes)} Nodes)"])
                set_item.setFlags(
                    set_item.flags() | QtCore.Qt.ItemIsUserCheckable)
                # v4.0.5b: default OFF. same rationale as load_set
                # above. Constraint glyph geometry only builds when
                # the user explicitly turns on a SID.
                set_item.setCheckState(0, QtCore.Qt.Unchecked)
                set_item.setData(0, QtCore.Qt.UserRole, ('constraint_set', sid))
                if all_nodes:
                    placeholder = QTreeWidgetItem(
                        [f"({len(all_nodes)} nodes - expand to load)"])
                    placeholder.setData(
                        0, QtCore.Qt.UserRole,
                        ('constraint_set_placeholder', sid))
                    set_item.addChild(placeholder)
                constr_items.append(set_item)
            constr_root.addChildren(constr_items)
            constr_root.setExpanded(True)
            _w_perf_event('top_row_flags', 'after_constr_root_addchildren',
                          flags=_flag_int(constr_root.flags()),
                          checkable=_is_checkable(constr_root.flags()))

        # --- Load Combinations (batched; combos themselves are flat) ---
        if hasattr(model, 'load_combinations') and model.load_combinations:
            combos_root = QTreeWidgetItem(
                tree,
                [f"Load Combinations ({len(model.load_combinations)})"])
            combos_root.setData(0, QtCore.Qt.UserRole, ('combos_root',))
            # v4.0.9: clear ItemIsUserCheckable from header row (see loads_root).
            combos_root.setFlags(combos_root.flags() & ~QtCore.Qt.ItemIsUserCheckable)
            _w_perf_event('top_row_flags', 'after_combos_root_create',
                          flags=_flag_int(combos_root.flags()),
                          checkable=_is_checkable(combos_root.flags()))
            combo_items = []
            for combo_sid, combo_data in sorted(
                    model.load_combinations.items()):
                try:
                    payload = read_combo_payload(combo_data)
                    scale = payload['scale']
                    components = payload['scale_factors']
                    load_ids = payload['load_ids']
                    n_members = len(load_ids)
                    parts = [f"SID {lid} x{sf:.3g}"
                             for sf, lid in zip(components, load_ids)]
                    detail = ", ".join(parts) if parts else "empty"
                    label = (f"LOAD SID {combo_sid}: S={scale:.3g}  "
                             f"({n_members} member{'s' if n_members != 1 else ''}) "
                             f"[{detail}]")
                except Exception:
                    label = f"LOAD SID {combo_sid}"
                combo_item = QTreeWidgetItem([label])
                combo_item.setData(
                    0, QtCore.Qt.UserRole, ('load_combo', combo_sid))
                combo_items.append(combo_item)
            combos_root.addChildren(combo_items)
            combos_root.setExpanded(False)  # collapsed by default

    def _on_load_tab_item_expanded(self, item):
        """v3.5.0 item 4: lazy children for Load / Constraint sets.

        Replaces a single ('..._placeholder', sid) child with the real
        entry items on first expand. Idempotent - subsequent expands
        are no-ops because the placeholder is gone."""
        if item is None or item.childCount() != 1:
            return
        first = item.child(0)
        data = first.data(0, QtCore.Qt.UserRole)
        if not isinstance(data, tuple) or len(data) < 2:
            return
        kind, sid = data[0], data[1]
        gen = getattr(self, 'current_generator', None)
        model = gen.model if gen else None
        if not model:
            return
        tree = self.loads_tab_tree
        tree.setUpdatesEnabled(False)
        try:
            if kind == 'load_set_placeholder':
                load_list = model.loads.get(sid, [])
                tempd_card = model.tempds.get(sid, None)
                kids = []
                for idx, card in enumerate(load_list):
                    entry_text = self._format_load_entry(card)
                    e = QTreeWidgetItem([entry_text])
                    e.setData(0, QtCore.Qt.UserRole,
                              ('load_entry', sid, idx))
                    kids.append(e)
                if tempd_card:
                    e = QTreeWidgetItem(
                        [self._format_load_entry(tempd_card)])
                    e.setData(0, QtCore.Qt.UserRole, ('tempd_entry', sid))
                    kids.append(e)
                # Remove placeholder, batch-add real children.
                item.removeChild(first)
                if kids:
                    item.addChildren(kids)
            elif kind == 'constraint_set_placeholder':
                spc_list = model.spcs.get(sid, [])
                kids = []
                for spc_card in spc_list:
                    for nid in spc_card.nodes:
                        e = QTreeWidgetItem(
                            [self._format_spc_entry(spc_card, nid)])
                        e.setData(0, QtCore.Qt.UserRole,
                                  ('constraint_entry', sid, nid))
                        kids.append(e)
                item.removeChild(first)
                if kids:
                    item.addChildren(kids)
        finally:
            tree.setUpdatesEnabled(True)

    def _populate_analysis_tab(self):
        """Populate the Analysis sidebar tab tree."""
        tree = self.analysis_tab_tree
        tree.clear()

        if not self.analysis_sets:
            return

        for sid, aset in sorted(self.analysis_sets.items()):
            prefix = "\u2605 " if sid == self.active_analysis_set_id else ""
            sol_str = f"SOL {aset.sol_type}" if aset.sol_type else "Punch"
            n_sc = len(aset.subcases)
            label = (f"{prefix}Set {sid}: {aset.name} ({sol_str}, "
                     f"{n_sc} subcase{'s' if n_sc != 1 else ''})")

            set_item = QTreeWidgetItem(tree, [label])
            set_item.setData(0, QtCore.Qt.UserRole, ('analysis_set', sid))

            # Show subcases as children
            for sc in aset.subcases:
                sc_parts = [f"SC {sc['id']}"]
                if sc.get('load_sid'):
                    sc_parts.append(f"LOAD={sc['load_sid']}")
                if sc.get('spc_sid'):
                    sc_parts.append(f"SPC={sc['spc_sid']}")
                sc_label = ", ".join(sc_parts)
                sc_item = QTreeWidgetItem(set_item, [sc_label])
                sc_item.setData(0, QtCore.Qt.UserRole,
                                ('analysis_subcase', sid, sc['id']))

            # Show params summary
            if aset.params:
                params_item = QTreeWidgetItem(
                    set_item, [f"Parameters ({len(aset.params)})"])
                params_item.setData(0, QtCore.Qt.UserRole,
                                    ('analysis_params', sid))

        tree.expandAll()

    def _show_loads_tab_context_menu(self, pos):
        """Unified context menu for the loads tab tree."""
        item = self.loads_tab_tree.itemAt(pos)
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if not data:
            return

        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        kind = data[0]

        if kind == 'loads_root':
            add_action = menu.addAction("Add Load...")
            add_action.triggered.connect(self._create_load)
        elif kind == 'load_set':
            sid = data[1]
            add_action = menu.addAction(f"Add to Load Set {sid}...")
            add_action.triggered.connect(lambda: self._add_to_load_set(sid))
            edit_action = menu.addAction(f"Edit Load Set {sid}...")
            edit_action.triggered.connect(lambda: self._open_load_set_manager(sid))
            menu.addSeparator()
            delete_action = menu.addAction(f"Delete Set {sid}")
            delete_action.triggered.connect(lambda: self._delete_load_set_from_tab(sid))
        elif kind == 'load_entry':
            sid, idx = data[1], data[2]
            remove_action = menu.addAction("Remove Entry")
            remove_action.triggered.connect(
                lambda: self._remove_specific_load_entry(sid, idx))
        elif kind == 'tempd_entry':
            sid = data[1]
            remove_action = menu.addAction("Remove TEMPD Entry")
            remove_action.triggered.connect(
                lambda: self._remove_tempd_entry(sid))
        elif kind == 'constraints_root':
            add_action = menu.addAction("Add Constraint...")
            add_action.triggered.connect(self._create_constraint)
        elif kind == 'constraint_set':
            sid = data[1]
            add_action = menu.addAction(f"Add to Constraint Set {sid}...")
            add_action.triggered.connect(self._create_constraint)
            menu.addSeparator()
            delete_action = menu.addAction(f"Delete Set {sid}")
            delete_action.triggered.connect(
                lambda: self._delete_constraint_set_from_tab(sid))
        elif kind == 'constraint_entry':
            sid, nid = data[1], data[2]
            remove_action = menu.addAction(f"Remove Node {nid} from Set")
            remove_action.triggered.connect(
                lambda: self._remove_specific_constraint_entry(sid, nid))
        elif kind == 'combos_root':
            add_action = menu.addAction("Add Combination...")
            add_action.triggered.connect(self._create_load_combination)
        elif kind == 'load_combo':
            # v3.5.0 item 2: full professional context menu
            # (Edit / Copy / Rename / Delete).
            combo_sid = data[1]
            edit_action = menu.addAction(f"Edit Combination SID {combo_sid}...")
            edit_action.triggered.connect(
                lambda: self._edit_load_combination(combo_sid))
            copy_action = menu.addAction(f"Copy Combination SID {combo_sid}...")
            copy_action.triggered.connect(
                lambda: self._copy_load_combination(combo_sid))
            rename_action = menu.addAction(f"Rename SID {combo_sid}...")
            rename_action.triggered.connect(
                lambda: self._rename_load_combination(combo_sid))
            menu.addSeparator()
            delete_action = menu.addAction(f"Delete Combination SID {combo_sid}")
            delete_action.triggered.connect(
                lambda: self._delete_load_combo_from_tab(combo_sid))

        if menu.actions():
            menu.exec(self.loads_tab_tree.viewport().mapToGlobal(pos))

    def _remove_load_entry_from_tab(self):
        """Remove the selected load entry from the Loads tab tree."""
        item = self.loads_tab_tree.currentItem()
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if data and data[0] == 'load_entry':
            self._remove_specific_load_entry(data[1], data[2])
        elif data and data[0] == 'tempd_entry':
            self._remove_tempd_entry(data[1])
        elif data and data[0] == 'load_set':
            QMessageBox.information(self, "Info",
                                    "Select a specific entry to remove, or use 'Remove Set'.")

    def _remove_specific_load_entry(self, sid, entry_index):
        """Remove a specific load entry by SID and index from model.loads."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        model = gen.model
        if sid not in model.loads:
            return
        load_list = list(model.loads[sid])
        if entry_index < 0 or entry_index >= len(load_list):
            return

        import copy
        new_list = [copy.deepcopy(c) for i, c in enumerate(load_list) if i != entry_index]
        if not new_list and sid not in model.tempds:
            # Removing last entry and no TEMPD for this SID - delete entire set
            from node_runner.commands import DeleteLoadCommand
            cmd = DeleteLoadCommand(sid)
            self.command_manager.execute(cmd, model)
        else:
            # Replace with remaining entries (or empty list if TEMPD still present)
            from node_runner.commands import ReplaceLoadSetCommand
            cmd = ReplaceLoadSetCommand(sid, new_list)
            self.command_manager.execute(cmd, model)
        self._populate_loads_tab()
        self._update_viewer()

    def _remove_tempd_entry(self, sid):
        """Remove a TEMPD entry by SID from model.tempds."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        model = gen.model
        if sid not in model.tempds:
            return
        # If SID has no other loads (only TEMPD), delete entire set (cleans both dicts)
        has_other_loads = sid in model.loads and model.loads[sid]
        if not has_other_loads:
            from node_runner.commands import DeleteLoadCommand
            cmd = DeleteLoadCommand(sid)
        else:
            # Only remove the TEMPD; keep model.loads[sid] intact
            from node_runner.commands import DeleteTempdCommand
            cmd = DeleteTempdCommand(sid)
        self.command_manager.execute(cmd, model)
        self._populate_loads_tab()
        self._populate_tree()
        self._update_viewer()

    def _remove_load_set_from_tab(self):
        """Remove the selected load set from the Loads tab."""
        item = self.loads_tab_tree.currentItem()
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if not data:
            return
        # Navigate to parent if an entry is selected
        sid = None
        if data[0] == 'load_set':
            sid = data[1]
        elif data[0] == 'load_entry':
            sid = data[1]
        if sid is None:
            return

        reply = QMessageBox.question(
            self, "Delete Load Set",
            f"Delete entire load set SID {sid}?",
            QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            self._delete_load_set_from_tab(sid)

    def _delete_load_set_from_tab(self, sid):
        """Delete a load set by SID."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        from node_runner.commands import DeleteLoadCommand
        cmd = DeleteLoadCommand(sid)
        self.command_manager.execute(cmd, gen.model)
        self._populate_loads_tab()
        self._populate_tree()
        self._update_viewer()

    def _remove_constraint_entry_from_tab(self):
        """Remove the selected constraint entry from the Loads tab."""
        item = self.loads_tab_tree.currentItem()
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if data and data[0] == 'constraint_entry':
            self._remove_specific_constraint_entry(data[1], data[2])
        elif data and data[0] == 'constraint_set':
            QMessageBox.information(self, "Info",
                                    "Select a specific entry to remove, or use 'Remove Set'.")

    def _remove_specific_constraint_entry(self, sid, node_id):
        """Remove a specific constraint entry (by node) from an SPC set."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        model = gen.model
        if sid not in model.spcs:
            return

        import copy
        # Build new SPC list filtering out the specified node
        new_list = []
        for spc_card in model.spcs[sid]:
            new_nodes = [n for n in spc_card.nodes if n != node_id]
            if new_nodes:
                new_spc = copy.deepcopy(spc_card)
                # Rebuild with filtered nodes
                new_list.append(new_spc)

        if not new_list:
            from node_runner.commands import DeleteConstraintCommand
            cmd = DeleteConstraintCommand(sid)
            self.command_manager.execute(cmd, model)
        else:
            from node_runner.commands import ReplaceConstraintSetCommand
            cmd = ReplaceConstraintSetCommand(sid, new_list)
            self.command_manager.execute(cmd, model)
        self._populate_loads_tab()
        self._update_viewer()

    def _remove_constraint_set_from_tab(self):
        """Remove the selected constraint set from the Loads tab."""
        item = self.loads_tab_tree.currentItem()
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        sid = None
        if data and data[0] == 'constraint_set':
            sid = data[1]
        elif data and data[0] == 'constraint_entry':
            sid = data[1]
        if sid is None:
            return

        reply = QMessageBox.question(
            self, "Delete Constraint Set",
            f"Delete entire constraint set SID {sid}?",
            QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            self._delete_constraint_set_from_tab(sid)

    def _delete_constraint_set_from_tab(self, sid):
        """Delete a constraint set by SID."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        from node_runner.commands import DeleteConstraintCommand
        cmd = DeleteConstraintCommand(sid)
        self.command_manager.execute(cmd, gen.model)
        self._populate_loads_tab()
        self._populate_tree()
        self._update_viewer()

    def _remove_load_combo_from_tab(self):
        """Remove the selected load combination from the Loads tab."""
        item = self.loads_tab_tree.currentItem()
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if not data or data[0] != 'load_combo':
            return
        self._delete_load_combo_from_tab(data[1])

    def _delete_load_combo_from_tab(self, combo_sid):
        """Delete a load combination by SID."""
        gen = getattr(self, 'current_generator', None)
        if not gen:
            return
        model = gen.model
        if hasattr(model, 'load_combinations') and combo_sid in model.load_combinations:
            reply = QMessageBox.question(
                self, "Delete Load Combination",
                f"Delete load combination SID {combo_sid}?",
                QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                del model.load_combinations[combo_sid]
                self._populate_loads_tab()
                self._populate_tree()
                self._update_status(f"Deleted load combination SID {combo_sid}.")

    def _open_find_entity_tool(self):
        """Launches a tool to find and zoom to selected entities."""
        if not self.current_generator or not self.current_generator.model.nodes:
            self._update_status("No model loaded to find entities in.", is_error=True)
            return

        entity_types = ["Node", "Element"]
        entity_type, ok = QInputDialog.getItem(self, "Select Entity Type", 
                                               "What type of entity do you want to find?", entity_types, 0, False)

        if ok and entity_type:
            if entity_type == 'Node':
                all_ids = self.current_generator.model.nodes.keys()
            else:
                all_ids = self.current_generator.model.elements.keys()
            select_by = self._build_select_by_data() if entity_type == 'Element' else None
            self._start_selection(entity_type, all_ids,
                                  self._on_find_entities_accept, select_by)

    def _on_find_entities_accept(self):
        """Callback for when the user confirms their selection in the Find tool."""
        selected_ids = self.selection_bar.get_selected_ids()
        entity_type = self.selection_bar.entity_type
        self._end_selection_mode()

        if not selected_ids:
            # --- FIX: Call the correct highlighting function with an empty list ---
            # This correctly removes any existing highlight actor from the scene.
            self._highlight_entities(entity_type, [])
            self._update_status("Highlight cleared.")
            return
        
        # --- Calculate bounding box of selection and zoom to it ---
        all_points = []
        if entity_type == 'Node':
            for nid in selected_ids:
                if nid in self.current_generator.model.nodes:
                    all_points.append(self.current_generator.model.nodes[nid].xyz)
        else: # Element
            nodes_to_get = set()
            for eid in selected_ids:
                if eid in self.current_generator.model.elements:
                    nodes_to_get.update(self.current_generator.model.elements[eid].nodes)
            for nid in nodes_to_get:
                 if nid in self.current_generator.model.nodes:
                    all_points.append(self.current_generator.model.nodes[nid].xyz)
        
        if all_points:
            bounds = pv.PolyData(np.array(all_points)).bounds
            self.plotter.reset_camera(bounds=bounds, render=False)
            self._highlight_entities(entity_type, selected_ids) # Re-highlight after zoom
        
        self._update_status(f"Found and zoomed to {len(selected_ids)} {entity_type.lower()}(s).")
    



    def _open_fuselage_generator(self):
        mat_model = None
        if self.current_generator:
            mat_model = self.current_generator.model
        else:
            mat_model = self._auto_load_materials()

        if not mat_model or not mat_model.materials:
            QMessageBox.warning(
                self, "No Materials",
                "Load a model with materials, or configure a materials "
                "library path in Preferences -> Library.")
            return

        dialog = GeneratorDialog(self); dialog.set_materials(mat_model)
        if dialog.exec():
            try:
                params = dialog.get_parameters()
                generator = NastranModelGenerator(params=params)
                for mid, mat in mat_model.materials.items():
                    generator.model.add_mat1(mid, mat.e, mat.g, mat.nu, comment=mat.comment)
                
                generator._generate_nodes()
                generator._generate_elements()
                generator.generate_floor_structure()
                self._ensure_default_coords(generator.model)

                self.command_manager.clear()
                self._update_viewer(generator); self._update_status("Fuselage model generated.")
            except Exception as e: self._update_status(f"Generation failed: {e}", is_error=True); QMessageBox.critical(self, "Error", f"Failed to generate model: {e}")



    def _parse_bdf_to_generator(self, filepath):
        """
        Robustly parses a BDF file into a new NastranModelGenerator object.
        Uses multiple fallback strategies including lenient card-by-card
        parsing that skips unparseable cards.

        Args:
            filepath (str): The path to the BDF file.

        Returns:
            tuple: (NastranModelGenerator, LenientResult or None, detected_format)
                Strict success → (generator, None, 'short'|'long'|'free')
                Lenient fallback → (generator, LenientResult, 'short'|'long'|'free')

        Raises:
            RuntimeError: If the file cannot be parsed by any strategy.
        """
        from node_runner.model import detect_bdf_field_format
        detected_format = detect_bdf_field_format(filepath)
        generator = NastranModelGenerator()
        model, lenient_result = NastranModelGenerator._read_bdf_robust(filepath)
        generator.model = model
        return generator, lenient_result, detected_format









    def _open_file_dialog(self):
        # Accept the full Nastran family: bulk data (.bdf), executive
        # data (.dat), free-form decks (.nas), punch (.pch), and
        # assembly (.asm). The latter two show up via INCLUDE statements
        # in real aerospace decks; users sometimes want to open them
        # standalone too.
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Open Nastran File", "",
            "Nastran Files (*.bdf *.dat *.nas *.pch *.asm);;All Files (*)")
        if not filepath:
            return

        from node_runner.workers import run_bdf_import_threaded

        def on_success(generator, lenient_result, detected_format):
            try:
                self._ensure_default_coords(generator.model)
                self.command_manager.clear()
                self.subcases = []
                self.geometry_store.clear()
                self.sol_type = None
                self.eigrl_cards = []
                self.analysis_sets.clear()
                self.active_analysis_set_id = None
                self.op2_results = None
                self.results_widget.setVisible(False)
                self.groups.clear()
                self._hidden_groups.clear()
                self._isolate_mode = None
                self._populate_groups_list()
                self._update_viewer(generator)
                self._on_model_loaded(filepath, detected_format)

                fname = os.path.basename(filepath)
                # Build a short summary of advanced-entity counts so
                # the user can see at a glance whether superelement
                # matrices and rigid/multi-point connections came
                # through cleanly.
                m = generator.model
                advanced = []
                dmig_total = sum(len(getattr(m, d, {}) or {})
                                 for d in ('dmig', 'dmij', 'dmiji', 'dmik', 'dmi'))
                if dmig_total:
                    advanced.append(f"{dmig_total} DMIG matrix(es)")
                rigid_n = len(getattr(m, 'rigid_elements', {}) or {})
                if rigid_n:
                    advanced.append(f"{rigid_n} rigid elements (RBE/RBAR)")
                mpc_n = sum(len(v) for v in (getattr(m, 'mpcs', {}) or {}).values())
                if mpc_n:
                    advanced.append(f"{mpc_n} MPC equation(s)")
                # Superelement bulk: count SEBULK/SECONCT entries if
                # pyNastran has them on the model.
                se_n = 0
                for attr in ('superelement_models', 'seloc', 'sebulk',
                             'seconct', 'setrans'):
                    v = getattr(m, attr, None)
                    if isinstance(v, dict):
                        se_n += len(v)
                if se_n:
                    advanced.append(f"{se_n} superelement card(s)")
                # v3.2.2: analysis-only cards stashed by the import path.
                # They never went through pyNastran so they're not in
                # any model dict; the count lives on the model itself.
                stash_counts = getattr(m, '_skipped_card_counts', {}) or {}
                if stash_counts:
                    total = sum(stash_counts.values())
                    top = ', '.join(
                        f'{n}: {c:,}' for n, c in
                        sorted(stash_counts.items(),
                               key=lambda kv: -kv[1])[:3])
                    advanced.append(
                        f"{total:,} analysis-only cards stashed "
                        f"({top}) - written back on export")
                # v3.3.0: surface coord systems whose rid pointed at a
                # coord not in the deck (typically defined in a missing
                # INCLUDE). _finalize_for_viewer falls them back to basic
                # so the import doesn't fail.
                coord_warnings = getattr(
                    m, '_coord_resolution_warnings', None) or []
                if coord_warnings:
                    advanced.append(
                        f"{len(coord_warnings)} coord system(s) with "
                        f"dangling parent - treated as basic")
                # v3.4.0 item 6: always show the model summary so the
                # user can spot when an entity class came back empty
                # (e.g. materials in a missing INCLUDE).
                n_nodes = len(getattr(m, 'nodes', {}) or {})
                n_elems = (len(getattr(m, 'elements', {}) or {})
                           + len(getattr(m, 'rigid_elements', {}) or {}))
                n_props = len(getattr(m, 'properties', {}) or {})
                n_mats = len(getattr(m, 'materials', {}) or {})
                n_coords = len(getattr(m, 'coords', {}) or {})
                advanced.insert(0,
                    f"{n_nodes:,} nodes, {n_elems:,} elements, "
                    f"{n_props:,} properties, {n_mats:,} materials, "
                    f"{n_coords:,} coord systems")
                # v3.4.0 item 5: RBE integrity warnings (populated by
                # _create_rbe_actors during scene-build).
                rbe_warnings = getattr(
                    m, '_rbe_integrity_warnings', None) or {}
                if rbe_warnings:
                    parts = []
                    if rbe_warnings.get('center_missing'):
                        parts.append(
                            f"{rbe_warnings['center_missing']} missing center")
                    if rbe_warnings.get('partial_legs'):
                        parts.append(
                            f"{rbe_warnings['partial_legs']} with partial legs")
                    if rbe_warnings.get('empty'):
                        parts.append(
                            f"{rbe_warnings['empty']} with no valid legs")
                    if parts:
                        advanced.append(
                            "rigid integrity: " + ", ".join(parts))
                advanced_suffix = ""
                if advanced:
                    advanced_suffix = " - includes " + ", ".join(advanced)

                if lenient_result and lenient_result.skipped:
                    n = len(lenient_result.skipped)
                    self._update_status(
                        f"Opened {fname} "
                        f"({n} card{'s' if n != 1 else ''} skipped)"
                        f"{advanced_suffix}")
                    LenientImportReportDialog(fname, lenient_result, self).exec()
                elif lenient_result:
                    self._update_status(
                        f"Opened {fname} "
                        f"(lenient mode, all cards OK){advanced_suffix}")
                else:
                    self._update_status(
                        f"Displayed {fname}{advanced_suffix}")
            except Exception as e:
                self._update_status("File open failed.", is_error=True)
                QMessageBox.critical(self, "Error", f"Post-parse handling failed: {e}")

        def on_failure(msg):
            self._update_status("File open failed.", is_error=True)
            QMessageBox.critical(self, "Error", f"Could not open/parse file: {msg}")

        thread, worker, dialog = run_bdf_import_threaded(
            self, filepath, on_success, on_failure,
        )
        # Hold references so QThread/Qt doesn't garbage-collect them.
        self._active_import_thread = thread
        self._active_import_worker = worker
        self._active_import_dialog = dialog


    def _import_cad_file(self):
        """Import a STEP/IGES/STL file via gmsh meshing."""
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Import CAD File", "",
            "CAD Files (*.step *.stp *.iges *.igs *.stl);;All Files (*)")
        if not filepath:
            return

        # Ensure we have a model
        if not self.current_generator:
            from node_runner.model import NastranModelGenerator
            gen = NastranModelGenerator()
            gen.create_new_model()
            self._ensure_default_coords(gen.model)
            self.command_manager.clear()
            self._update_viewer(gen)

        available_pids = list(self.current_generator.model.properties.keys()) or [1]
        dlg = ImportCADDialog(filepath, available_pids, self)
        if not dlg.exec():
            return

        settings = dlg.get_settings()

        # Ensure the PID exists
        pid = settings['pid']
        model = self.current_generator.model
        if pid not in model.properties:
            # Create a default PSHELL with thickness 1.0
            mid = min(model.materials.keys()) if model.materials else 1
            if mid not in model.materials:
                model.add_mat1(mid, 70000.0, 0.33, rho=2.7e-9,
                               comment='Auto-created for CAD import')
            model.add_pshell(pid, mid1=mid, t=1.0,
                             comment='Auto-created for CAD import')

        try:
            cmd = ImportCADCommand(
                settings['filepath'], settings['mesh_size'],
                settings['elem_preference'], pid)
            self.command_manager.execute(cmd, model)
            n_nodes = len(cmd._created_nids)
            n_elems = len(cmd._created_eids)
            QMessageBox.information(
                self, "CAD Import Complete",
                f"Created {n_nodes} nodes and {n_elems} elements "
                f"from {os.path.basename(filepath)}.")
            self._update_status(
                f"Imported CAD: {n_nodes} nodes, {n_elems} elements.")
            self._update_viewer(self.current_generator, reset_camera=True)
        except Exception as e:
            QMessageBox.critical(self, "Import Error", str(e))
            self._update_status(f"CAD import failed: {e}", is_error=True)

    def _import_file_dialog(self):
        """Handles the 'Import' file action, showing an options dialog."""
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Import Nastran File", "",
            "Nastran Files (*.bdf *.dat *.nas *.pch *.asm);;All Files (*)")
        if not filepath:
            self._update_status("Import cancelled.")
            return

        # If there's no model loaded, the only option is to create a new one.
        if not self.current_generator or not self.current_generator.model.nodes:
            option = 'new'
        else:
            dialog = ImportOptionsDialog(self)
            if dialog.exec():
                option = dialog.get_selected_option()
            else:
                self._update_status("Import cancelled.")
                return

        filename = os.path.basename(filepath)

        if option == 'new':
            from node_runner.workers import run_bdf_import_threaded

            def on_success(generator, lenient_result, detected_format):
                try:
                    self._ensure_default_coords(generator.model)
                    self.command_manager.clear()
                    self.groups.clear()
                    self._hidden_groups.clear()
                    self._isolate_mode = None
                    self._populate_groups_list()
                    self._update_viewer(generator)
                    self._on_model_loaded(filepath, detected_format)

                    if lenient_result and lenient_result.skipped:
                        n = len(lenient_result.skipped)
                        self._update_status(
                            f"Opened {filename} ({n} card{'s' if n != 1 else ''} skipped)")
                        LenientImportReportDialog(filename, lenient_result, self).exec()
                    elif lenient_result:
                        self._update_status(
                            f"Opened {filename} (lenient mode, all cards OK)")
                    else:
                        self._update_status(f"Opened new model from {filename}")
                except Exception as e:
                    self._update_status("File import failed.", is_error=True)
                    QMessageBox.critical(self, "Error", f"Post-parse handling failed: {e}")

            def on_failure(msg):
                self._update_status("File import failed.", is_error=True)
                QMessageBox.critical(self, "Error", f"Could not import file: {msg}")

            thread, worker, dlg = run_bdf_import_threaded(
                self, filepath, on_success, on_failure,
            )
            self._active_import_thread = thread
            self._active_import_worker = worker
            self._active_import_dialog = dlg
            return

        # Append path stays synchronous: it operates on the current model in
        # place, which is harder to thread safely without rework. It is also
        # typically fast (small/medium decks).
        try:
            if option == 'append':
                summary = self.current_generator.import_and_append_bdf(filepath)
                self._ensure_default_coords(self.current_generator.model)

                summary_parts = [f"{count} {name}" for name, count in summary.items() if count > 0]
                if summary_parts:
                    message = f"Appended {', '.join(summary_parts)} from {filename}."
                else:
                    message = f"No new cards were appended from {filename}."

                self._update_status(message)
                self._update_viewer(self.current_generator, reset_camera=False)
        except Exception as e:
            self._update_status("File import failed.", is_error=True)
            QMessageBox.critical(self, "Error", f"Could not import file: {e}")



    def _auto_load_materials(self):
        """v5.0.0 item 7: load materials from the user-configured library
        path (Preferences -> Library) instead of a hardcoded
        ``materials.bdf`` in the install directory. The legacy
        ``materials.bdf`` is auto-migrated to the new setting on first
        run for users upgrading from v4.x.
        """
        from PySide6.QtCore import QSettings
        from node_runner.profiling import perf_event
        settings = QSettings("NodeRunner", "NodeRunner")
        materials_path = settings.value("materials/library_path", "", type=str)

        # One-time migration from the legacy auto-search to a saved path.
        legacy = os.path.join(PROJECT_ROOT, "materials.bdf")
        if not materials_path and os.path.exists(legacy):
            materials_path = legacy
            settings.setValue("materials/library_path", materials_path)
            try:
                perf_event('materials', 'library_migrated',
                           from_path=legacy)
            except Exception:
                pass
            self._update_status(
                "Adopted legacy materials.bdf as the materials library "
                "(change via Preferences -> Library).")

        if not materials_path:
            try:
                perf_event('materials', 'library_not_configured')
            except Exception:
                pass
            return None
        if not os.path.exists(materials_path):
            try:
                perf_event('materials', 'library_missing',
                           path=materials_path)
            except Exception:
                pass
            self._update_status(
                f"Material library not found: {materials_path} "
                f"(configure via Preferences -> Library).",
                is_error=True)
            return None
        try:
            mat_model = BDF(debug=False)
            mat_model.read_bdf(materials_path, punch=True)
            try:
                perf_event('materials', 'library_loaded',
                           path=materials_path,
                           n_materials=len(mat_model.materials))
            except Exception:
                pass
            self._update_status(
                f"Loaded {len(mat_model.materials)} materials from library.")
            return mat_model
        except Exception as e:
            self._update_status(
                f"Could not parse material library: {e}", is_error=True)
        return None

    def _save_configuration(self):
        if not (self.current_generator and self.current_generator.params): QMessageBox.warning(self, "No Data", "Only generator parameters can be saved."); return
        if not (save_path := QFileDialog.getSaveFileName(self, "Save Configuration", "", "JSON Files (*.json)")[0]): return
        try:
            with open(save_path, 'w') as f: json.dump(self.current_generator.params, f, indent=4)
            self._update_status(f"Config saved to {os.path.basename(save_path)}")
        except Exception as e: self._update_status(f"Error saving config: {e}", is_error=True)

    def _load_configuration(self):
        if not (load_path := QFileDialog.getOpenFileName(self, "Load Configuration", "", "JSON Files (*.json)")[0]): return
        try:
            with open(load_path, 'r') as f: config = json.load(f)
            dialog = GeneratorDialog(self);
            
            mat_model = self.current_generator.model if self.current_generator else self._auto_load_materials()
            if not mat_model:
                 QMessageBox.warning(self, "No Materials", "Could not find materials to use for generator."); return
            dialog.set_materials(mat_model)

            for key, w in [('fuselage_radius',dialog.radius_input),('num_bays',dialog.bays_input),('frame_spacing',dialog.spacing_input),('skin_thickness',dialog.skin_thick_input),('skin_material_id',dialog.skin_mat_combo),('num_stringers',dialog.num_stringers_input),('stringer_section_type',dialog.stringer_type_combo),('stringer_material_id',dialog.stringer_mat_combo),('frame_section_type',dialog.frame_type_combo),('frame_material_id',dialog.frame_mat_combo)]:
                if key in config: getattr(w, 'setText' if isinstance(w, QLineEdit) else 'setCurrentText')(str(config[key]))
            self._update_status(f"Config loaded from {os.path.basename(load_path)}.")
            if dialog.exec():
                params = dialog.get_parameters()
                generator = NastranModelGenerator(params=params)
                for mid, mat in mat_model.materials.items():
                    generator.model.add_mat1(mid, mat.e, mat.g, mat.nu, comment=mat.comment)
                
                generator._generate_nodes(); generator._generate_elements()
                self._ensure_default_coords(generator.model)
                self.command_manager.clear()
                self._update_viewer(generator); self._update_status("Model generated from loaded configuration.")
        except Exception as e: self._update_status(f"Error loading config: {e}", is_error=True); QMessageBox.critical(self, "Error", f"Could not load config: {e}")
        
    def _open_color_manager(self):
        if not self.current_grid: self._update_status("No model loaded.", is_error=True); return
        current_map = self.type_color_map if self.color_mode == "type" else self.pid_color_map; dialog = ColorManagerDialog(current_map, self)
        if dialog.exec():
            (self.type_color_map if self.color_mode == "type" else self.pid_color_map).update(dialog.color_map)
            self._update_plot_visibility(); self._update_status("Colors updated.")
            
    # --- Display Settings, Screenshot, Presentation Mode ---

    def _load_display_settings(self) -> dict:
        """Load display settings from QSettings."""
        from PySide6.QtCore import QSettings
        s = QSettings("NodeRunner", "NodeRunner")
        return {
            'bg_mode': s.value("display/bg_mode", "solid"),
            'bg_color': s.value("display/bg_color", "#1e1e2e"),
            'bg_top_color': s.value("display/bg_top", "#1a2a4a"),
            'bg_bottom_color': s.value("display/bg_bottom", "#0a0a1e"),
            'bg_image_path': s.value("display/bg_image", ""),
            'show_axes': s.value("display/show_axes", True, type=bool),
            'axes_color': s.value("display/axes_color", "#cdd6f4"),
            'edge_color': s.value("display/edge_color", "#000000"),
            'node_size': int(s.value("display/node_size", 8)),
            'line_width': int(s.value("display/line_width", 2)),
            'selection_color': s.value("display/selection_color", "#ffff00"),
            'watermark_text': s.value("display/watermark_text", ""),
            'watermark_color': s.value("display/watermark_color", "#ffffff"),
            'watermark_position': s.value("display/watermark_position", "Bottom-Right"),
            'screenshot_resolution': s.value("display/screenshot_res", "Current Viewport"),
            'custom_width': int(s.value("display/custom_width", 1920)),
            'custom_height': int(s.value("display/custom_height", 1080)),
            'transparent_bg': s.value("display/transparent_bg", False, type=bool),
            'antialiasing': s.value("display/antialiasing", True, type=bool),
            'file_format': s.value("display/file_format", "PNG"),
        }

    def _save_display_settings(self):
        """Save display settings to QSettings."""
        from PySide6.QtCore import QSettings
        s = QSettings("NodeRunner", "NodeRunner")
        ds = self._display_settings
        s.setValue("display/bg_mode", ds.get('bg_mode', 'solid'))
        s.setValue("display/bg_color", ds.get('bg_color', '#1e1e2e'))
        s.setValue("display/bg_top", ds.get('bg_top_color', '#1a2a4a'))
        s.setValue("display/bg_bottom", ds.get('bg_bottom_color', '#0a0a1e'))
        s.setValue("display/bg_image", ds.get('bg_image_path', ''))
        s.setValue("display/show_axes", ds.get('show_axes', True))
        s.setValue("display/axes_color", ds.get('axes_color', '#cdd6f4'))
        s.setValue("display/edge_color", ds.get('edge_color', '#000000'))
        s.setValue("display/node_size", ds.get('node_size', 8))
        s.setValue("display/line_width", ds.get('line_width', 2))
        s.setValue("display/selection_color", ds.get('selection_color', '#ffff00'))
        s.setValue("display/watermark_text", ds.get('watermark_text', ''))
        s.setValue("display/watermark_color", ds.get('watermark_color', '#ffffff'))
        s.setValue("display/watermark_position", ds.get('watermark_position', 'Bottom-Right'))
        s.setValue("display/screenshot_res", ds.get('screenshot_resolution', 'Current Viewport'))
        s.setValue("display/custom_width", ds.get('custom_width', 1920))
        s.setValue("display/custom_height", ds.get('custom_height', 1080))
        s.setValue("display/transparent_bg", ds.get('transparent_bg', False))
        s.setValue("display/antialiasing", ds.get('antialiasing', True))
        s.setValue("display/file_format", ds.get('file_format', 'PNG'))

    def _apply_background(self):
        """Apply current background settings to the plotter."""
        import os
        ds = self._display_settings
        bg_mode = ds.get('bg_mode', 'solid')

        # Remove any existing background image
        try:
            self.plotter.renderer.RemoveAllViewProps()
        except Exception:
            pass

        if bg_mode == 'gradient':
            top = ds.get('bg_top_color', '#1a2a4a')
            bottom = ds.get('bg_bottom_color', '#0a0a1e')
            self.plotter.set_background(bottom, top=top)
        elif bg_mode == 'image':
            fallback = ds.get('bg_color', '#1e1e2e')
            self.plotter.set_background(fallback)
            img_path = ds.get('bg_image_path', '')
            if img_path and os.path.exists(img_path):
                try:
                    self.plotter.add_background_image(img_path, scale=1.0, as_global=True)
                except Exception as e:
                    self._update_status(f"Failed to load background image: {e}", is_error=True)
        else:
            color = ds.get('bg_color', '#1e1e2e')
            self.plotter.set_background(color)
        self.plotter.render()

    def _open_display_settings(self):
        """Open the Display Settings dialog."""
        from node_runner.dialogs.tools import DisplaySettingsDialog
        dialog = DisplaySettingsDialog(settings=self._display_settings.copy(), parent=self)
        if dialog.exec():
            new_settings = dialog.get_settings()
            capture = new_settings.pop('capture_requested', False)

            # Update display settings
            self._display_settings.update(new_settings)
            self._save_display_settings()
            self._apply_background()

            # Apply axes visibility
            if hasattr(self, 'axes_actor') and self.axes_actor:
                self.axes_actor.SetVisibility(new_settings.get('show_axes', True))

            # Apply axes label color
            axes_color_hex = new_settings.get('axes_color', '#cdd6f4')
            color = QtGui.QColor(axes_color_hex)
            self._set_axes_label_color(
                (color.redF(), color.greenF(), color.blueF()))

            self._update_plot_visibility()
            self._update_status("Display settings updated.")

            # If capture was requested, do it now
            if capture:
                self._capture_screenshot()

    def _capture_screenshot(self):
        """Take a screenshot with current display settings."""
        from PySide6.QtWidgets import QFileDialog
        ds = self._display_settings
        fmt = ds.get('file_format', 'PNG').lower()
        filters = {
            'png': "PNG Files (*.png)",
            'jpg': "JPEG Files (*.jpg *.jpeg)",
            'bmp': "BMP Files (*.bmp)",
            'tiff': "TIFF Files (*.tiff *.tif)",
        }
        file_filter = filters.get(fmt, "PNG Files (*.png)")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Screenshot", f"screenshot.{fmt}", file_filter)
        if not path:
            return

        # Determine resolution
        res_text = ds.get('screenshot_resolution', 'Current Viewport')
        if res_text == '1920 x 1080 (Full HD)':
            w, h = 1920, 1080
        elif res_text == '2560 x 1440 (2K)':
            w, h = 2560, 1440
        elif res_text == '3840 x 2160 (4K)':
            w, h = 3840, 2160
        elif res_text == 'Custom':
            w = ds.get('custom_width', 1920)
            h = ds.get('custom_height', 1080)
        else:
            # Current Viewport
            w, h = None, None

        try:
            transparent = ds.get('transparent_bg', False)
            kwargs = {'transparent_background': transparent}
            if w and h:
                kwargs['window_size'] = (w, h)

            self.plotter.screenshot(path, **kwargs)
            self._update_status(f"Screenshot saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Screenshot Error",
                                f"Failed to save screenshot:\n{str(e)}")

    def _toggle_presentation_mode(self, checked):
        """Toggle presentation mode - fullscreen with UI hidden."""
        self.menuBar().setVisible(not checked)
        self.statusBar().setVisible(not checked)
        # Hide the sidebar tabs container
        sidebar_parent = self.sidebar_tabs.parent()
        if sidebar_parent:
            sidebar_parent.setVisible(not checked)
        if checked:
            self.showFullScreen()
        else:
            self.showNormal()

    def eventFilter(self, obj, event):
        """Catch ESC on the plotter to exit presentation mode.

        v4.0.3 (Stage K): also catch every event reaching tree_widget
        (and its viewport) when profiling is enabled. Pure diagnostic
       . we never consume the event. Lets us prove or disprove
        whether user clicks even reach Qt before reaching our slots.
        """
        # v4.0.3 (Stage K) / v4.0.4 (Stage R): tree-widget click trace.
        # v4.0.4 enriches MouseButtonPress with `itemAt(pos)` + tree
        # state so we can prove or disprove whether the click lands
        # on a real item and whether Qt is suppressing signals.
        try:
            tree = getattr(self, 'tree_widget', None)
            if tree is not None and (obj is tree or obj is tree.viewport()):
                etype = int(event.type())
                # Only log the few event types we actually care about
                # so we don't drown the log in MouseMove etc.
                CARED_ABOUT = (
                    QtCore.QEvent.MouseButtonPress,
                    QtCore.QEvent.MouseButtonRelease,
                    QtCore.QEvent.MouseButtonDblClick,
                    QtCore.QEvent.KeyPress,
                    QtCore.QEvent.FocusIn,
                    QtCore.QEvent.FocusOut,
                    QtCore.QEvent.Enter,
                    QtCore.QEvent.Leave,
                )
                if event.type() in CARED_ABOUT:
                    from node_runner.profiling import perf_event
                    pos = None
                    px = py = None
                    btn = None
                    try:
                        if hasattr(event, 'position'):
                            p = event.position()
                            px, py = int(p.x()), int(p.y())
                            pos = f"({px},{py})"
                        if hasattr(event, 'button'):
                            btn = int(event.button())
                    except Exception:
                        pass

                    # v4.0.4 (Stage R). for MouseButtonPress on the
                    # viewport, log full tree state to nail down why
                    # itemChanged/itemClicked don't fire.
                    extra = {}
                    if (event.type() == QtCore.QEvent.MouseButtonPress
                            and obj is tree.viewport()
                            and px is not None):
                        try:
                            # itemAt expects viewport coords.
                            hit_item = tree.itemAt(px, py)
                        except Exception:
                            hit_item = None
                        if hit_item is not None:
                            try:
                                hit_label = hit_item.text(0)
                            except Exception:
                                hit_label = '?'
                            try:
                                hit_data = hit_item.data(0, QtCore.Qt.UserRole)
                            except Exception:
                                hit_data = None
                            try:
                                hit_state = int(hit_item.checkState(0))
                            except Exception:
                                hit_state = -1
                            try:
                                hit_checkable = bool(
                                    hit_item.flags() & QtCore.Qt.ItemIsUserCheckable)
                            except Exception:
                                hit_checkable = None
                            extra.update(
                                item_at_hit=True,
                                hit_label=hit_label,
                                hit_data=str(hit_data)[:80],
                                hit_state=hit_state,
                                hit_checkable=hit_checkable,
                            )
                        else:
                            extra['item_at_hit'] = False
                        try:
                            extra['tree_signals_blocked'] = bool(tree.signalsBlocked())
                        except Exception:
                            pass
                        try:
                            extra['viewport_signals_blocked'] = bool(
                                tree.viewport().signalsBlocked())
                        except Exception:
                            pass
                        try:
                            extra['tree_enabled'] = bool(tree.isEnabled())
                            extra['tree_visible'] = bool(tree.isVisible())
                            extra['viewport_enabled'] = bool(
                                tree.viewport().isEnabled())
                            extra['updates_enabled'] = bool(
                                tree.updatesEnabled())
                        except Exception:
                            pass

                    perf_event(
                        'tree_event_filter', 'event',
                        target=('viewport' if obj is tree.viewport() else 'widget'),
                        type=etype,
                        pos=pos,
                        button=btn,
                        **extra)
        except Exception:
            pass

        if (event.type() == QtCore.QEvent.KeyPress
                and event.key() == QtCore.Qt.Key_Escape
                and self.presentation_action.isChecked()):
            self.presentation_action.setChecked(False)
            return True  # consumed
        return super().eventFilter(obj, event)

    def keyPressEvent(self, event):
        """Handle key press events for presentation mode exit."""
        if event.key() == QtCore.Qt.Key_Escape and self.presentation_action.isChecked():
            self.presentation_action.setChecked(False)
        else:
            super().keyPressEvent(event)

    def _set_render_style(self, style):
        self.render_style = style
        self.transparent_shells_action.setEnabled(style == "surface")
        # Sync sidebar combo
        style_map = {"surface": 0, "wireframe": 1, "surface_only": 2, "points": 3}
        idx = style_map.get(style, 0)
        self.render_style_combo.blockSignals(True)
        self.render_style_combo.setCurrentIndex(idx)
        self.render_style_combo.blockSignals(False)
        self._update_plot_visibility(); self._refresh_labels()
        self._update_status(f"Render style: {style}.")

    def _on_render_style_combo_changed(self, index):
        style_map = {0: "surface", 1: "wireframe", 2: "surface_only", 3: "points"}
        style = style_map.get(index, "surface")
        self.render_style = style
        self.transparent_shells_action.setEnabled(style == "surface")
        # Sync menu actions
        if style == "wireframe":
            self.style_wire_action.setChecked(True)
        else:
            self.style_surf_action.setChecked(True)
        self._update_plot_visibility(); self._refresh_labels()
        self._update_status(f"Render style: {style}.")

    def _on_color_mode_combo_changed(self, index):
        mode_map = {0: "property", 1: "type", 2: "quality", 3: "results"}
        mode = mode_map.get(index, "property")
        self._set_coloring_mode(mode)

    def _toggle_shell_transparency(self, state):
        self.shell_opacity = 0.4 if state else 1.0
        self.shell_opacity_slider.blockSignals(True)
        self.shell_opacity_slider.setValue(int(self.shell_opacity * 100))
        self.shell_opacity_slider.blockSignals(False)
        self._update_plot_visibility(); self._refresh_labels()
        self._update_status(f"Shell transparency {'ON' if state else 'OFF'}.")

    def _toggle_perspective_view(self, state):
        self.plotter.camera.SetParallelProjection(not state); self.plotter.render(); self._update_status(f"Perspective view {'ON' if state else 'OFF'}.")

    # --- Display tab handlers ---
    def _on_node_size_changed(self, value):
        self.node_size = value
        self._update_plot_visibility()

    def _pick_node_color(self):
        color = QColorDialog.getColor(QColor(self.node_color), self, "Node Color")
        if color.isValid():
            self.node_color = color.name()
            self.node_color_btn.setStyleSheet(f"background-color: {self.node_color}; border: 1px solid #555;")
            self._update_plot_visibility()

    def _on_shell_opacity_changed(self, value):
        self.shell_opacity = value / 100.0
        self._update_plot_visibility()

    def _on_elem_shrink_changed(self, value):
        self.elem_shrink = value / 100.0
        self._update_plot_visibility()

    def _pick_edge_color(self):
        color = QColorDialog.getColor(QColor(self.edge_color), self, "Edge Color")
        if color.isValid():
            self.edge_color = color.name()
            self.edge_color_btn.setStyleSheet(f"background-color: {self.edge_color}; border: 1px solid #555;")
            self._update_plot_visibility()

    def _on_edge_width_changed(self, value):
        self.edge_width = value
        self._update_plot_visibility()

    def _on_beam_width_changed(self, value):
        self.beam_width = value
        self._update_plot_visibility()

    def _on_bg_preset_changed(self, preset_name):
        preset = BACKGROUND_PRESETS.get(preset_name)
        if not preset:
            return
        if preset["mode"] == "solid":
            self.plotter.set_background(preset["color"])
        elif preset["mode"] == "gradient":
            self.plotter.set_background(preset["bottom"], top=preset["top"])
        # Debounced render: a burst of preset changes coalesces to one frame.
        self._request_render()
        
    def _set_coloring_mode(self, mode):
        from node_runner.profiling import perf_event as _cm_perf_event
        self.color_mode = mode
        # v4.0.15: pre-built per-category actors are color-mode-
        # dependent. Tear them down + rebuild so the next fast-path
        # category toggle renders the right colors. Skipped when
        # no pre-built actors exist (before first viewer build or
        # when the prior mode was unsupported so they were never
        # built). Future-proof: the build site reads `self.color_mode`
        # and consults `_PRE_BUILT_SUPPORTED_COLOR_MODES` so adding
        # a new color mode "just works" here.
        if getattr(self, '_category_actors', None):
            try:
                for ntype in list(self._category_actors.keys()):
                    self.plotter.remove_actor(f'cat_{ntype}', render=False)
                self._category_actors = {}
                self._category_visible_ntypes = set()
                self._using_legacy_mesh = False
                if (self.current_grid is not None
                        and 'type' in self.current_grid.cell_data):
                    self._build_pre_built_category_actors()
                _cm_perf_event('actors', 'rebuild_on_color_mode',
                               mode=mode,
                               n_built=len(self._category_actors))
            except Exception as _exc:
                _cm_perf_event('actors', 'rebuild_on_color_mode_failed',
                               mode=mode,
                               exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
        elif self.current_grid is not None:
            # No pre-built actors yet (prior mode was unsupported, or
            # they were never built). If the new mode IS supported,
            # build them now so the next category toggle hits the
            # fast path with the right colors.
            if mode in _PRE_BUILT_SUPPORTED_COLOR_MODES:
                try:
                    self._build_pre_built_category_actors()
                    _cm_perf_event('actors', 'rebuild_on_color_mode',
                                   mode=mode,
                                   n_built=len(self._category_actors))
                except Exception as _exc:
                    _cm_perf_event('actors', 'rebuild_on_color_mode_failed',
                                   mode=mode,
                                   exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
        # Sync sidebar combo
        mode_map = {"property": 0, "type": 1, "quality": 2, "results": 3}
        idx = mode_map.get(mode, 0)
        self.color_mode_combo.blockSignals(True)
        self.color_mode_combo.setCurrentIndex(idx)
        self.color_mode_combo.blockSignals(False)
        is_quality_mode = (mode == "quality")
        is_results_mode = (mode == "results")
        self.quality_widget.setVisible(is_quality_mode)
        # v5.2.0 item 47: the legacy results_widget stays hidden -- the
        # Post-Processing Toolbox (Results sidebar tab) is the canonical
        # UI now. The widget itself still exists as a state holder.
        if not getattr(self, '_legacy_results_widget_hidden', False):
            self.results_widget.setVisible(is_results_mode and self.op2_results is not None)

        if is_quality_mode and self.current_generator:
            # Populate the dropdown with available metrics the first time
            if self.quality_metric_combo.count() == 0:
                 has_quads = any(e.type in ['CQUAD4', 'CMEMBRAN'] for e in self.current_generator.model.elements.values())
                 metrics = ["Aspect Ratio", "Skew", "Min Angle", "Max Angle"]
                 if has_quads:
                     metrics = ["Warping", "Jacobian", "Taper"] + metrics
                 self.quality_metric_combo.addItems(metrics)

        if is_results_mode and not self.op2_results:
            QMessageBox.information(self, "No Results",
                                    "No results loaded. Use File > Load Results (OP2) first.")

        self._update_plot_visibility()
        self._update_status(f"Color mode: {mode}.")
        
    def _rebuild_plot(self, visibility_mask=None):
        from node_runner.profiling import perf_event as _rb_perf_event
        # v4.0.14 (Fix C v2): coexistence with pre-built per-category
        # actors. When pre-built actors exist AND no cell-level
        # filter is active (visibility_mask is all-True), skip the
        # legacy shells/beams extract+split+add. the pre-built
        # actors already render the mesh at full density. When a
        # filter IS active, hide pre-built actors and add legacy
        # filtered shells/beams. The user's category toggle on the
        # next click (after filter clears) restores pre-built
        # visibility via the fast path.
        _has_pre_built = bool(getattr(self, '_category_actors', None))
        _filter_active = False
        if visibility_mask is not None:
            try:
                _filter_active = not bool(np.all(visibility_mask))
            except Exception:
                _filter_active = True
        _skip_legacy_mesh = _has_pre_built and not _filter_active
        if _skip_legacy_mesh:
            # Restore pre-built actor visibility per category state.
            for ntype, actor in self._category_actors.items():
                vis = self._compute_category_actor_visibility(ntype)
                actor.SetVisibility(bool(vis) if vis is not None else True)
            self._using_legacy_mesh = False
            _rb_perf_event('mesh', 'rebuild_plot_skipped',
                           reason='pre_built_no_filter')
        else:
            # Hide pre-built so they don't double-render with the
            # filtered shells/beams added below.
            for actor in (getattr(self, '_category_actors', None) or {}).values():
                actor.SetVisibility(False)
            self._using_legacy_mesh = True
            if _has_pre_built:
                _rb_perf_event('mesh', 'rebuild_plot_entered',
                               reason='filter_active' if _filter_active else 'no_pre_built')
        self.plotter.remove_actor(['nodes_actor', 'shells', 'beams', 'beam_sections_3d',
                                   'selection_highlight'])
        # v4.0.14: when skipping the legacy mesh, still let the
        # nodes_actor / selection_highlight code below run, then
        # early-return before the extract_cells + add_mesh.

        nodes_visible = False
        if (nodes_item_list := self._find_tree_items(('group', 'nodes'))):
            if nodes_item_list[0].checkState(0) == QtCore.Qt.Checked:
                nodes_visible = True

        # v5.2.1 Round-2 fix: predict whether the deformation block
        # later in this method will re-issue the node cloud at deformed
        # positions. If yes, skip the early draw -- otherwise the user
        # sees the cloud snap from undeformed-to-deformed within a
        # single render frame on each animation tick (perceived as a
        # flicker against the bars).
        _will_deform_nodes = False
        try:
            if (self.color_mode == 'results'
                    and getattr(self, 'op2_results', None) is not None
                    and hasattr(self, 'results_type_combo')
                    and self.results_type_combo.currentText() in (
                        'Displacement', 'Eigenvector')):
                if getattr(self, '_anim_override_scale', None) is not None:
                    _will_deform_nodes = abs(self._anim_override_scale) > 1e-30
                else:
                    try:
                        _ds = float(self.deformation_scale_input.text() or 0)
                    except (ValueError, AttributeError):
                        _ds = 0.0
                    _will_deform_nodes = _ds != 0.0
        except Exception:
            _will_deform_nodes = False

        if (nodes_visible
                and self.current_grid
                and self.current_grid.n_points > 0
                and not _will_deform_nodes):
            node_points = self.current_grid.points
            if self.current_node_ids_sorted and (
                (hasattr(self, '_isolate_mode') and self._isolate_mode) or
                (hasattr(self, '_hidden_groups') and self._hidden_groups)
            ):
                all_nids = np.array(self.current_node_ids_sorted)
                if self._isolate_mode:
                    group_data = self.groups.get(self._isolate_mode, {})
                    show_nids = set(group_data.get("nodes", []))
                    if show_nids:
                        point_mask = np.isin(all_nids, list(show_nids))
                        if point_mask.any():
                            node_points = self.current_grid.points[point_mask]
                        else:
                            node_points = None
                    else:
                        node_points = None
                elif self._hidden_groups:
                    hidden_nids = set()
                    for grp_name in self._hidden_groups:
                        grp = self.groups.get(grp_name, {})
                        hidden_nids.update(grp.get("nodes", []))
                    if hidden_nids:
                        point_mask = ~np.isin(all_nids, list(hidden_nids))
                        if point_mask.any():
                            node_points = self.current_grid.points[point_mask]
                        else:
                            node_points = None
            if node_points is not None and len(node_points) > 0:
                self._add_node_cloud(node_points, name='nodes_actor')

        # v4.0.14: when pre-built actors handle the mesh and no
        # filter is active, the legacy extract+split+add chain
        # below is wasted work. Pre-built actors are already
        # visible from the gating above. Render and return.
        if _skip_legacy_mesh:
            self.plotter.render()
            return

        # v3.2.0: prefer the display-LOD'd grid when no visibility mask
        # is active. With a visibility mask we have to extract from the
        # full grid because cell indices are full-grid indices. Cell
        # picking always targets self.current_grid so EID resolution
        # is unaffected by display LOD.
        from node_runner.profiling import perf_stage as _perf_stage
        if visibility_mask is not None:
            visible_indices = np.where(visibility_mask)[0]
            if visible_indices.size > 0:
                with _perf_stage('rebuild_plot', 'extract_cells_visible',
                                 n_visible=int(visible_indices.size),
                                 n_total=int(visibility_mask.size)):
                    grid_to_render = self.current_grid.extract_cells(visible_indices)
            else:
                self.plotter.remove_actor(['shells', 'beams'])
                self.plotter.render()
                return
        else:
            grid_to_render = (
                self.current_display_grid
                if self.current_display_grid is not None
                else self.current_grid
            )

        # --- FIX: Check if cell data (like 'is_shell') exists before trying to use it ---
        if grid_to_render and grid_to_render.n_cells > 0 and 'is_shell' in grid_to_render.cell_data:
            with _perf_stage('rebuild_plot', 'split_shells_beams',
                             n_cells=int(grid_to_render.n_cells)):
                shells = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 1)
                beams = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 0)

            # Apply element shrink if set
            if self.elem_shrink > 0 and shells.n_cells > 0:
                shells = shells.shrink(1.0 - self.elem_shrink)

            # v4.0.7: compute point normals on the shells surface so
            # phong shading has data to interpolate. Without normals,
            # VTK falls back to flat-equivalent shading regardless of
            # the actor's `interpolation` property, which is why the
            # v4.0.5 / v4.0.6 shading toggle visibly did nothing. We
            # convert to PolyData first (UnstructuredGrid doesn't
            # auto-compute normals) via extract_surface, then run the
            # vtkPolyDataNormals filter via .compute_normals(). The
            # cost is one-shot per visibility update on the LOD'd
            # mesh (~248k cells on a large production deck, fast in C++).
            if shells.n_cells > 0:
                # v4.0.9: cache the computed-normals shells PolyData
                # across visibility cycles. The compute_normals call on
                # 2.36M LOD'd cells costs ~1.1 s per click on a large production deck, but
                # the inputs (visibility-filtered shells subset + shrink
                # factor) are identical across many consecutive clicks
                # when the user is toggling unrelated things like loads
                # or coord systems. Cache key is (id(current_grid),
                # n_cells, EID-array hash, shrink). invalidates
                # automatically when _update_viewer runs (new
                # current_grid) or the user filters by PID/Type/Material
                # (different EID set). Hashing the EID array is ~2 ms;
                # saves ~1.1 s on every cache hit.
                from node_runner.profiling import perf_event as _norm_perf_event
                cache_key = None
                try:
                    eid_hash = None
                    if 'EID' in shells.cell_data:
                        eid_hash = hash(bytes(shells.cell_data['EID']))
                    cache_key = (
                        id(self.current_grid),
                        int(shells.n_cells),
                        eid_hash,
                        float(self.elem_shrink),
                    )
                except Exception:
                    cache_key = None
                cached = getattr(self, '_shells_normals_cache', None)
                if cache_key is not None and cached is not None and cached[0] == cache_key:
                    shells = cached[1]
                    _norm_perf_event('shells_normals', 'cache_hit',
                                     n_cells=int(shells.n_cells))
                else:
                    try:
                        with _perf_stage('rebuild_plot', 'shells_compute_normals',
                                         n_cells=int(shells.n_cells)):
                            shells = shells.extract_surface().compute_normals(
                                cell_normals=False,
                                point_normals=True,
                                consistent_normals=True,
                                auto_orient_normals=False,
                                split_vertices=False,
                                non_manifold_traversal=True,
                                inplace=False,
                            )
                        if cache_key is not None:
                            self._shells_normals_cache = (cache_key, shells)
                            _norm_perf_event('shells_normals', 'cache_miss_recompute',
                                             n_cells=int(shells.n_cells))
                    except Exception as _exc:
                        _norm_perf_event('rebuild_plot', 'shells_normals_failed',
                                         exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")

            # Determine render style for pyvista
            pv_style = self.render_style
            if pv_style == "surface_only":
                pv_style = "surface"
            show_edges = self.render_style == "surface"

            # Common shell kwargs
            # v4.0.13 (Fix A): reset_camera=False on every legacy-mesh
            # add_mesh path. Without it, every Plates/Beams/Solids/etc.
            # toggle was snapping the camera back to fit the full
            # model. v4.0.11 fixed the GLYPH paths but missed the
            # MESH paths (shell_kw / beam_kw feed all ~10 add_mesh
            # sites in _rebuild_plot's color-mode branches).
            shell_kw = dict(style=pv_style, show_edges=show_edges,
                            opacity=self.shell_opacity, name="shells",
                            reset_camera=False)
            if show_edges:
                shell_kw['edge_color'] = self.edge_color
                shell_kw['line_width'] = self.edge_width
            beam_kw = dict(render_lines_as_tubes=True,
                           line_width=self.beam_width, name="beams",
                           reset_camera=False)

            if self.color_mode == "type":
                if shells.n_cells > 0:
                    self.plotter.add_mesh(shells, color=self.type_color_map["Shells"], **shell_kw)
                if beams.n_cells > 0:
                    self.plotter.add_mesh(beams, color=self.type_color_map["Beams"], **beam_kw)
            elif self.color_mode == "property":
                # v4.0.17: per-cell RGB instead of categorical LUT.
                # Identical pid_color_map[PID] color for every cell,
                # consistent across legacy and pre-built code paths,
                # consistent across rebuilds. See
                # `_attach_per_cell_rgb_for_property_mode` for details.
                if shells.n_cells > 0:
                    self._attach_per_cell_rgb_for_property_mode(shells)
                    self.plotter.add_mesh(
                        shells, scalars='cell_rgb', rgb=True,
                        show_scalar_bar=False, **shell_kw)
                if beams.n_cells > 0:
                    self._attach_per_cell_rgb_for_property_mode(beams)
                    self.plotter.add_mesh(
                        beams, scalars='cell_rgb', rgb=True,
                        show_scalar_bar=False, **beam_kw)
            elif self.color_mode == "quality":
                if self.quality_results is None:
                    self.quality_results = self.current_generator.calculate_element_quality()

                metric_text = self.quality_metric_combo.currentText()
                metric_key_map = {
                    "Warping": "warp", "Aspect Ratio": "aspect", "Skew": "skew",
                    "Jacobian": "jacobian", "Taper": "taper",
                    "Min Angle": "min_angle", "Max Angle": "max_angle",
                }
                metric_key = metric_key_map.get(metric_text, "aspect")

                if shells.n_cells > 0 and self.quality_results and metric_key:
                    scalars = np.full(shells.n_cells, 0.0, dtype=float)

                    for i, eid in enumerate(shells.cell_data['EID']):
                        scalars[i] = self.quality_results.get(eid, {}).get(metric_key, 0.0)

                    self.plotter.add_mesh(shells, scalars=scalars, cmap='viridis',
                                          scalar_bar_args={'title': metric_text}, **shell_kw)
                if beams.n_cells > 0:
                    self.plotter.add_mesh(beams, color=self.type_color_map["Beams"], **beam_kw)
            elif self.color_mode == "results" and self.op2_results:
                # v5.3.0 item 60: capture animation-tick context so the
                # catchers downstream can be cross-referenced against
                # the actor_snapshot events.
                _anim_tick = getattr(self, '_anim_tick_counter', -1)
                _anim_active = bool(getattr(self, '_anim_timer', None)
                                    and self._anim_timer.isActive())
                # v5.2.1 Round-3 fix: the legacy contour code at
                # lines 12084/12113/12146/12213 references
                # ``model.elements`` but ``model`` was never defined
                # in this scope -- a latent bug that only fires on
                # beam-only decks (the shells branch happens to be a
                # dead path for those decks). Assign it once here.
                model = (self.current_generator.model
                         if self.current_generator is not None else None)
                sc_id = self.results_subcase_combo.currentData()
                sc_data = self.op2_results['subcases'].get(sc_id, {}) if sc_id is not None else {}
                result_type = self.results_type_combo.currentText()
                component = self.results_component_combo.currentText()
                # v5.2.2 issue 3: honor the Contour Style dropdown.
                # The four styles map to add_mesh kwargs:
                #   filled        -> style='surface',   show_edges=False
                #   filled_edges  -> style='surface',   show_edges=True
                #   bands         -> n_colors=<levels>  (discrete LUT)
                #   none          -> short-circuit to solid type-color
                # v5.3.2 item 67: deformation always uses the
                # displacement triplet from the current subcase --
                # regardless of which Output Vector the user picked.
                # Independent of contour entirely. If the active
                # vector is an eigenvector mode, use that mode's
                # displacement; otherwise use 'displacements'.
                if self._anim_override_scale is not None:
                    _deform_scale = self._anim_override_scale
                else:
                    try:
                        _deform_scale = float(self.deformation_scale_input.text() or 0)
                    except ValueError:
                        _deform_scale = 0.0
                _dvkind = getattr(self, '_output_vector_kind', 'displacement')
                _dvmode = int(getattr(self, '_output_vector_mode_idx', -1))
                _deform_data = None
                if _deform_scale != 0:
                    if _dvkind == 'eigenvector':
                        eigs = sc_data.get('eigenvectors') or []
                        if 0 <= _dvmode < len(eigs):
                            _deform_data = eigs[_dvmode] or {}
                    else:
                        # Default deformation source = displacements.
                        # Applies even when the user picked a stress
                        # vector for contour -- the deformation is the
                        # geometric displacement.
                        _deform_data = sc_data.get('displacements') or {}
                if _deform_data:
                    # v5.2.1 item 49: deep-copy before mutating points.
                    import time as _time
                    _t0 = _time.perf_counter()
                    grid_to_render = grid_to_render.copy(deep=True)
                    try:
                        from node_runner.profiling import perf_event
                        perf_event('results', 'deform_grid_copied',
                                   wall_s=round(_time.perf_counter() - _t0, 4),
                                   n_points=int(grid_to_render.n_points))
                    except Exception:
                        pass
                    _deformed_points = grid_to_render.points.copy()
                    for _i, _nid in enumerate(self.current_node_ids_sorted):
                        _vals = _deform_data.get(_nid)
                        if _vals and _i < len(_deformed_points):
                            _deformed_points[_i] += np.array(_vals[:3]) * _deform_scale
                    grid_to_render.points = _deformed_points
                    shells = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 1)
                    beams = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 0)
                    try:
                        if nodes_visible and len(_deformed_points) > 0:
                            self._add_node_cloud(_deformed_points, name='nodes_actor')
                            try:
                                from node_runner.profiling import perf_event
                                perf_event('results', 'node_cloud',
                                           points_source='deformed',
                                           n=int(len(_deformed_points)))
                            except Exception:
                                pass
                    except Exception:
                        pass

                _cstyle = getattr(self, '_results_contour_style', 'none')
                if _cstyle == 'none':
                    # v5.3.1 item 62: render property RGB on the
                    # (possibly deformed) sub-grids. Deformation was
                    # already applied above; here we just paint colour.
                    self._clear_scalar_bars()
                    if shells.n_cells > 0:
                        self._attach_per_cell_rgb_for_property_mode(shells)
                        self.plotter.add_mesh(
                            shells, scalars='cell_rgb', rgb=True,
                            show_scalar_bar=False, **shell_kw)
                    if beams.n_cells > 0:
                        self._attach_per_cell_rgb_for_property_mode(beams)
                        self.plotter.add_mesh(
                            beams, scalars='cell_rgb', rgb=True,
                            show_scalar_bar=False, **beam_kw)
                    if self.show_beam_sections_check.isChecked() and self.current_generator:
                        self._render_beam_sections()
                    self.plotter.render()
                    return
                # Apply the contour style by overriding shell_kw /
                # beam_kw before any add_mesh call below uses them.
                if _cstyle == 'filled_edges':
                    shell_kw['show_edges'] = True
                elif _cstyle == 'filled':
                    shell_kw['show_edges'] = False
                # 'bands' is handled later via n_colors on the mapper.
                _bands_n_colors = (int(getattr(self, '_results_n_levels', 9))
                                   if _cstyle == 'bands' else None)
                # v5.3.1 item 65: sweep ghost scalar bars before adding
                # the fresh one. Without this, switching Contour Vector
                # accumulates one bar per title we've ever rendered.
                self._clear_scalar_bars()
                # v5.3.1 item 64: scalar-bar kwargs built from the
                # user's selected orientation preference (Right vertical
                # or Bottom horizontal). Title carries a newline so the
                # top tick doesn't overlap (item 63).
                _scalar_bar_kw = self._build_scalar_bar_kw()

                if result_type in ("Displacement", "Eigenvector", "SPC Forces"):
                    key_map = {"Displacement": "displacements",
                               "Eigenvector": "eigenvectors",
                               "SPC Forces": "spc_forces"}
                    raw = sc_data.get(key_map.get(result_type, ""), {})
                    if result_type == "Eigenvector" and isinstance(raw, list):
                        mode_idx = max(0, self.results_mode_combo.currentIndex())
                        data_dict = raw[mode_idx] if mode_idx < len(raw) else {}
                    else:
                        data_dict = raw if isinstance(raw, dict) else {}
                    comp_idx = {"T1 (X)": 0, "T2 (Y)": 1, "T3 (Z)": 2,
                                "R1 (RX)": 3, "R2 (RY)": 4, "R3 (RZ)": 5}.get(component)

                    # Build point scalars on the full grid
                    point_scalars = np.zeros(self.current_grid.n_points)
                    for i, nid in enumerate(self.current_node_ids_sorted):
                        vals = data_dict.get(nid)
                        if vals:
                            if component == "Magnitude":
                                point_scalars[i] = np.linalg.norm(vals[:3])
                            elif comp_idx is not None:
                                point_scalars[i] = vals[comp_idx]

                    # Apply deformation if scale != 0
                    if self._anim_override_scale is not None:
                        deform_scale = self._anim_override_scale
                    else:
                        try:
                            deform_scale = float(self.deformation_scale_input.text() or 0)
                        except ValueError:
                            deform_scale = 0.0
                    # v5.3.0 item 59: deformation is now applied
                    # earlier in this method using the DEFORM VECTOR
                    # (independent of result_type). The block below is
                    # kept as a no-op fallback for cases where the
                    # deform vector wasn't applicable.
                    if False and deform_scale != 0 and result_type in ("Displacement", "Eigenvector"):
                        # v5.2.1 item 49: deep-copy before mutating points.
                        # Without this, ``grid_to_render.points = deformed_points``
                        # mutates ``self.current_grid.points`` in place (the
                        # else branch at line ~11788 assigned grid_to_render
                        # = current_grid by reference). That made the model
                        # compound on every toggle.
                        import time as _time
                        _t0 = _time.perf_counter()
                        grid_to_render = grid_to_render.copy(deep=True)
                        try:
                            from node_runner.profiling import perf_event
                            perf_event('results', 'deform_grid_copied',
                                       wall_s=round(_time.perf_counter() - _t0, 4),
                                       n_points=int(grid_to_render.n_points))
                        except Exception:
                            pass
                        disp_data = data_dict
                        deformed_points = grid_to_render.points.copy()
                        for i, nid in enumerate(self.current_node_ids_sorted):
                            vals = disp_data.get(nid)
                            if vals and i < len(deformed_points):
                                deformed_points[i] += np.array(vals[:3]) * deform_scale
                        grid_to_render.points = deformed_points
                        shells = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 1)
                        beams = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 0)
                        # v5.2.1 item 51: the early node-cloud draw
                        # (line ~11760) used undeformed current_grid
                        # points and is now stale. Re-issue with the
                        # deformed points so the dots ride the bars.
                        try:
                            if nodes_visible and len(deformed_points) > 0:
                                # Honor any isolate/hide filtering by
                                # re-using the same mask logic. For the
                                # common all-visible case just pass the
                                # full deformed point set.
                                self._add_node_cloud(
                                    deformed_points, name='nodes_actor')
                                try:
                                    from node_runner.profiling import perf_event
                                    perf_event(
                                        'results', 'node_cloud',
                                        points_source='deformed',
                                        n=int(len(deformed_points)))
                                except Exception:
                                    pass
                        except Exception:
                            pass

                    self.current_grid.point_data['result_scalars'] = point_scalars
                    # v5.2.0 item 43: data conversion + palette.
                    cmap_name = getattr(self, '_results_palette', 'jet')
                    conv_mode = getattr(self, '_results_data_conversion', 'average')
                    if shells.n_cells > 0:
                        from node_runner.results_conversion import (
                            convert_nodal_data)
                        # Build cell_node_ids list for the visible shells.
                        cell_node_ids = []
                        for eid in shells.cell_data['EID']:
                            elem = (model.elements.get(int(eid))
                                    or model.rigid_elements.get(int(eid)))
                            cell_node_ids.append(
                                np.asarray(elem.nodes if elem else [],
                                           dtype=np.int64))
                        # Need the point_scalars restricted to the shells'
                        # extracted point set, but extract_cells preserves
                        # original points so passing the full nodal array
                        # paired with the original sorted node IDs works.
                        node_ids_sorted = np.asarray(
                            self.current_node_ids_sorted, dtype=np.int64)
                        point_scalars_shells, cell_scalars_shells = (
                            convert_nodal_data(
                                point_scalars, node_ids_sorted, conv_mode,
                                cell_node_ids))
                        # v5.3.2 item 70: title now has a trailing blank
                        # line so VTK reserves extra vertical space
                        # between the title block and the top tick.
                        title = f"{result_type}\n{component}\n"
                        if conv_mode != 'average':
                            title = (f"{result_type}\n{component}\n"
                                     f"({_CONV_LABEL.get(conv_mode, conv_mode)})\n")
                        # v5.3.1 item 65: explicit clim so changing
                        # vector / data conversion actually refreshes
                        # the scalar bar range. Honors the user's
                        # Auto-Range checkbox.
                        _user_auto = getattr(self, '_results_user_auto_range', True)
                        if _user_auto:
                            _clim_disp = None  # PyVista auto-computes
                        else:
                            _clim_disp = (float(self._results_user_vmin),
                                          float(self._results_user_vmax))
                        _add_kw = {
                            'cmap': cmap_name,
                            'n_colors': (_bands_n_colors or 256),
                            'scalar_bar_args': {
                                'title': title,
                                'n_labels': (_bands_n_colors + 1
                                             if _bands_n_colors else 5),
                                **_scalar_bar_kw},
                        }
                        if _clim_disp is not None:
                            _add_kw['clim'] = _clim_disp
                        if cell_scalars_shells is not None:
                            self.plotter.add_mesh(
                                shells, scalars=cell_scalars_shells,
                                **_add_kw, **shell_kw)
                        else:
                            # Average path: project to per-cell mean so the
                            # legacy shell mapping stays cell-scalar based.
                            shell_scalars = np.zeros(shells.n_cells)
                            for i, eid in enumerate(shells.cell_data['EID']):
                                elem = (model.elements.get(eid)
                                        or model.rigid_elements.get(eid))
                                if elem:
                                    node_ids_elem = elem.nodes
                                    vals_list = [data_dict.get(nid)
                                                 for nid in node_ids_elem
                                                 if data_dict.get(nid)]
                                    if vals_list:
                                        if component == "Magnitude":
                                            shell_scalars[i] = np.mean(
                                                [np.linalg.norm(v[:3])
                                                 for v in vals_list])
                                        elif comp_idx is not None:
                                            shell_scalars[i] = np.mean(
                                                [v[comp_idx] for v in vals_list])
                            self.plotter.add_mesh(
                                shells, scalars=shell_scalars,
                                **_add_kw, **shell_kw)
                    if beams.n_cells > 0:
                        # v5.2.1 Round-2 fix: use cell_data scalars
                        # (mean of 2 corner nodes per bar) instead of
                        # point_data. PyVista's render_lines_as_tubes
                        # mode -- on by default in beam_kw -- extrudes
                        # line cells into 3D tube surfaces, and the
                        # tube filter doesn't pass point arrays
                        # through cleanly. Cell scalars survive the
                        # extrusion. Wrap in try/except + telemetry so
                        # any future regression here is loud, not silent.
                        try:
                            cell_scalars_beams = np.zeros(beams.n_cells, dtype=float)
                            for i, eid in enumerate(beams.cell_data['EID']):
                                elem = (model.elements.get(int(eid))
                                        or model.rigid_elements.get(int(eid)))
                                if elem is None:
                                    continue
                                vals_list = []
                                for nid in elem.nodes:
                                    v = data_dict.get(nid)
                                    if not v:
                                        continue
                                    if component == "Magnitude" and len(v) >= 3:
                                        vals_list.append(
                                            (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5)
                                    elif comp_idx is not None and comp_idx < len(v):
                                        vals_list.append(float(v[comp_idx]))
                                if vals_list:
                                    cell_scalars_beams[i] = (
                                        sum(vals_list) / len(vals_list))
                            self.plotter.add_mesh(
                                beams, scalars=cell_scalars_beams,
                                **_add_kw, **beam_kw)
                            try:
                                from node_runner.profiling import perf_event
                                perf_event('results', 'beam_contour',
                                           applied=True,
                                           n_cells=int(beams.n_cells))
                            except Exception:
                                pass
                        except Exception as _beam_exc:
                            try:
                                from node_runner.profiling import perf_event
                                perf_event('results', 'beam_contour_failed',
                                           err=str(_beam_exc)[:200])
                            except Exception:
                                pass
                            # Fall back to solid color so bars still show.
                            self.plotter.add_mesh(
                                beams,
                                color=self.type_color_map["Beams"],
                                **beam_kw)

                elif result_type == "Stress":
                    stress_data = sc_data.get("stresses", {})
                    comp_key_map = {
                        "von Mises": "von_mises", "Max Principal": "max_principal",
                        "Min Principal": "min_principal",
                        "XX": "oxx", "YY": "oyy", "XY": "txy",
                    }
                    comp_key = comp_key_map.get(component, "von_mises")
                    cmap_name = getattr(self, '_results_palette', 'jet')
                    conv_mode = getattr(self, '_results_data_conversion', 'average')

                    if shells.n_cells > 0:
                        from node_runner.results_conversion import (
                            convert_element_data)
                        # Build the eid -> value dict from stress data
                        # (filtered to our visible cell EIDs).
                        cell_eids = np.asarray(
                            shells.cell_data['EID'], dtype=np.int64)
                        elem_vals = {int(eid): float(
                            stress_data.get(int(eid), {}).get(comp_key, 0.0))
                            for eid in cell_eids}
                        # Cell -> nodes connectivity, needed for averaging
                        # paths.
                        cell_node_ids = []
                        for eid in cell_eids:
                            elem = (model.elements.get(int(eid))
                                    or model.rigid_elements.get(int(eid)))
                            cell_node_ids.append(
                                np.asarray(elem.nodes if elem else [],
                                           dtype=np.int64))
                        node_ids_sorted = np.asarray(
                            self.current_node_ids_sorted, dtype=np.int64)
                        point_scalars_s, cell_scalars_s = convert_element_data(
                            elem_vals, conv_mode, cell_eids, cell_node_ids,
                            node_ids_sorted, n_points=self.current_grid.n_points)
                        # v5.3.2 item 70: trailing blank line for
                        # title breathing room.
                        title = f"Stress\n{component}\n"
                        if conv_mode != 'no_avg':
                            title = (f"Stress\n{component}\n"
                                     f"({_CONV_LABEL.get(conv_mode, conv_mode)})\n")
                        # v5.3.1 item 65: explicit clim + scalar bar
                        # refresh on every stress render too.
                        _user_auto = getattr(self, '_results_user_auto_range', True)
                        if _user_auto:
                            _clim_stress = None
                        else:
                            _clim_stress = (float(self._results_user_vmin),
                                            float(self._results_user_vmax))
                        _add_kw = {
                            'cmap': cmap_name,
                            'n_colors': (_bands_n_colors or 256),
                            'scalar_bar_args': {
                                'title': title,
                                'n_labels': (_bands_n_colors + 1
                                             if _bands_n_colors else 5),
                                **_scalar_bar_kw},
                        }
                        if _clim_stress is not None:
                            _add_kw['clim'] = _clim_stress
                        if point_scalars_s is not None:
                            # v5.3.3 fix 2: attach the scalar array to
                            # the SHELLS sub-grid directly. The prior
                            # code attached to ``current_grid`` but
                            # ``extract_cells`` returns a NEW grid that
                            # doesn't share its point_data dict -- so
                            # PyVista's ``scalars='stress_point_scalars'``
                            # lookup on shells silently failed and the
                            # whole mesh rendered with no scalars
                            # (invisible). Map full-grid scalars onto
                            # the shells sub-grid via the helper.
                            shell_pts_scalars = self._project_point_scalars_to_subgrid(
                                point_scalars_s, shells)
                            if shell_pts_scalars is not None:
                                shells.point_data['stress_point_scalars'] = shell_pts_scalars
                                self.plotter.add_mesh(
                                    shells, scalars='stress_point_scalars',
                                    **_add_kw, **shell_kw)
                            else:
                                # Fallback to a per-cell mean if the
                                # projection failed.
                                self.plotter.add_mesh(
                                    shells, scalars=cell_scalars_s,
                                    **_add_kw, **shell_kw)
                        else:
                            self.plotter.add_mesh(
                                shells, scalars=cell_scalars_s,
                                **_add_kw, **shell_kw)
                    if beams.n_cells > 0:
                        # v5.2.1 Round-2 fix: cell_data scalars on
                        # bars so render_lines_as_tubes works. Bar
                        # stress per-cell comes from looking up the
                        # element's stress dict directly (CBEAM/CBAR
                        # stresses are element-centered, same as shells).
                        try:
                            cell_scalars_beams_s = np.zeros(
                                beams.n_cells, dtype=float)
                            for i, eid in enumerate(beams.cell_data['EID']):
                                s = stress_data.get(int(eid), {}) or {}
                                cell_scalars_beams_s[i] = float(
                                    s.get(comp_key, 0.0))
                            self.plotter.add_mesh(
                                beams, scalars=cell_scalars_beams_s,
                                **_add_kw, **beam_kw)
                            try:
                                from node_runner.profiling import perf_event
                                perf_event('results', 'beam_contour',
                                           applied=True,
                                           n_cells=int(beams.n_cells),
                                           kind='stress')
                            except Exception:
                                pass
                        except Exception as _beam_exc:
                            try:
                                from node_runner.profiling import perf_event
                                perf_event('results', 'beam_contour_failed',
                                           err=str(_beam_exc)[:200],
                                           kind='stress')
                            except Exception:
                                pass
                            self.plotter.add_mesh(
                                beams,
                                color=self.type_color_map["Beams"],
                                **beam_kw)
                else:
                    # Fallback: property coloring (v4.0.17 per-cell RGB).
                    if shells.n_cells > 0:
                        self._attach_per_cell_rgb_for_property_mode(shells)
                        self.plotter.add_mesh(
                            shells, scalars='cell_rgb', rgb=True,
                            show_scalar_bar=False, **shell_kw)
                    if beams.n_cells > 0:
                        self._attach_per_cell_rgb_for_property_mode(beams)
                        self.plotter.add_mesh(
                            beams, scalars='cell_rgb', rgb=True,
                            show_scalar_bar=False, **beam_kw)
        else:
             # If there's no 'is_shell' data, it means there are no elements to draw.
             # Ensure any old shell/beam actors are removed.
             self.plotter.remove_actor(['shells', 'beams'])

        # Beam 3D section rendering
        if self.show_beam_sections_check.isChecked() and self.current_generator:
            self._render_beam_sections()

        self.plotter.render()

    def _render_beam_sections(self):
        """Render 3D extruded beam cross-sections for all CBAR/CBEAM elements."""
        import pyvista as pv

        gen = self.current_generator
        model = gen.model
        all_node_coords = {nid: node.xyz for nid, node in model.nodes.items()}

        all_verts = []
        all_faces = []
        vertex_offset = 0

        beam_types = {'CBEAM', 'CBAR'}
        for eid, elem in model.elements.items():
            if elem.type not in beam_types:
                continue

            pid = getattr(elem, 'pid', None)
            if pid is None:
                continue

            polygon_2d = gen.get_beam_section_polygon(pid)
            if polygon_2d is None:
                continue

            # Filter out hole markers
            outer_pts = [(y, z) for y, z in polygon_2d
                         if y is not None and z is not None]
            if len(outer_pts) < 3:
                continue

            # Get beam axis
            try:
                p1 = np.array(all_node_coords[elem.nodes[0]], dtype=float)
                p2 = np.array(all_node_coords[elem.nodes[1]], dtype=float)
            except (KeyError, IndexError):
                continue

            axis = p2 - p1
            length = np.linalg.norm(axis)
            if length < 1e-12:
                continue
            axis_dir = axis / length

            # Orientation vector: g0 node or x vector
            try:
                if hasattr(elem, 'g0') and elem.g0 is not None and elem.g0 > 0:
                    g0_pos = np.array(all_node_coords[elem.g0], dtype=float)
                    v_orient = g0_pos - p1
                elif hasattr(elem, 'x') and elem.x is not None:
                    v_orient = np.array(elem.x, dtype=float)
                else:
                    # Default: pick arbitrary perpendicular
                    ref = np.array([0.0, 0.0, 1.0])
                    if abs(np.dot(axis_dir, ref)) > 0.9:
                        ref = np.array([0.0, 1.0, 0.0])
                    v_orient = ref
            except (KeyError, TypeError):
                ref = np.array([0.0, 0.0, 1.0])
                if abs(np.dot(axis_dir, ref)) > 0.9:
                    ref = np.array([0.0, 1.0, 0.0])
                v_orient = ref

            # Build local coordinate frame
            v_y = v_orient - np.dot(v_orient, axis_dir) * axis_dir
            v_y_len = np.linalg.norm(v_y)
            if v_y_len < 1e-12:
                ref = np.array([1.0, 0.0, 0.0])
                v_y = ref - np.dot(ref, axis_dir) * axis_dir
                v_y_len = np.linalg.norm(v_y)
                if v_y_len < 1e-12:
                    continue
            v_y = v_y / v_y_len
            v_z = np.cross(axis_dir, v_y)

            # Extrude polygon: create quads connecting p1-section to p2-section
            n_pts = len(outer_pts)
            for i, (y, z) in enumerate(outer_pts):
                pt1 = p1 + y * v_y + z * v_z
                pt2 = p2 + y * v_y + z * v_z
                all_verts.append(pt1)
                all_verts.append(pt2)

            # Create quad faces along the extrusion
            for i in range(n_pts - 1):
                i0 = vertex_offset + i * 2
                i1 = vertex_offset + i * 2 + 1
                i2 = vertex_offset + (i + 1) * 2 + 1
                i3 = vertex_offset + (i + 1) * 2
                all_faces.extend([4, i0, i1, i2, i3])

            # End caps (triangulate as fan)
            cap1_center = vertex_offset + n_pts * 2
            cap2_center = vertex_offset + n_pts * 2 + 1
            c1 = np.mean([p1 + y * v_y + z * v_z
                          for y, z in outer_pts], axis=0)
            c2 = np.mean([p2 + y * v_y + z * v_z
                          for y, z in outer_pts], axis=0)
            all_verts.append(c1)
            all_verts.append(c2)

            for i in range(n_pts - 1):
                # Cap at p1
                all_faces.extend([3, cap1_center,
                                  vertex_offset + (i + 1) * 2,
                                  vertex_offset + i * 2])
                # Cap at p2
                all_faces.extend([3, cap2_center,
                                  vertex_offset + i * 2 + 1,
                                  vertex_offset + (i + 1) * 2 + 1])

            vertex_offset = len(all_verts)

        if all_verts:
            mesh = pv.PolyData(np.array(all_verts), faces=all_faces)
            color = self.type_color_map.get("Beams", "#89b4fa")
            self.plotter.add_mesh(mesh, color=color, opacity=0.7,
                                  name='beam_sections_3d',
                                  reset_camera=False)

    def _refresh_labels(self):
        if self.show_node_labels_check.isChecked(): self._toggle_node_labels(False); self._toggle_node_labels(True)
        if self.show_elem_labels_check.isChecked(): self._toggle_element_labels(False); self._toggle_element_labels(True)
        
    def _toggle_node_labels(self, state):
        self.plotter.remove_actor('node_labels', render=False)
        if state and self.current_grid and self.current_grid.n_points > 0:
            color = 'white' if self.is_dark_theme else 'black'
            self.plotter.add_point_labels(
                self.current_grid.points, self.current_node_ids_sorted,
                name='node_labels', font_size=14, text_color=color,
                shape=None, bold=False, always_visible=True,
                reset_camera=False, render=True
            )
        self._update_status(f"Node labels {'ON' if state else 'OFF'}")

    def _toggle_element_labels(self, state):
        self.plotter.remove_actor('elem_labels', render=False)
        if state and self.current_grid and self.current_grid.n_cells > 0 and 'EID' in self.current_grid.cell_data:
            color = 'white' if self.is_dark_theme else 'black'
            self.plotter.add_point_labels(
                self.current_grid.cell_centers().points, self.current_grid.cell_data["EID"],
                name='elem_labels', font_size=14, text_color=color,
                shape=None, bold=False, always_visible=True,
                reset_camera=False, render=True
            )
        self._update_status(f"Element labels {'ON' if state else 'OFF'}")
    
    def _find_and_zoom_to_entity(self):
        if not self.current_generator or not self.current_grid: self._update_status("No model.", is_error=True); return
        try: entity_id = int(self.entity_id_input.text()); etype = self.entity_type_combo.currentText()
        except ValueError: self._update_status("Invalid ID.", is_error=True); return
        
        coords, found = None, False
        model = self.current_generator.model

        if etype == "Node":
            if entity_id in model.nodes:
                coords = model.nodes[entity_id].get_position()
                found = True
        elif etype == "Element":
            if entity_id in self.current_grid.cell_data["EID"]:
                idx = np.where(self.current_grid.cell_data["EID"] == entity_id)[0][0]
                coords = self.current_grid.get_cell(idx).center
                found = True

        self.plotter.remove_actor('highlight_actor', render=False)
        if found:
            # v3.5.0 item 3: highlight color via user preference.
            hl_color = getattr(self, 'highlight_color', '#fab387')
            self.plotter.add_points(np.array(coords), color=hl_color, point_size=15, render_points_as_spheres=True, name='highlight_actor'); self.plotter.fly_to(coords)
            self._update_status(f"Found {etype} {entity_id}.")
        else: self._update_status(f"{etype} {entity_id} not found.", is_error=True)

    def _clear_entity_highlight(self): 
        self.plotter.remove_actor('highlight_actor', render=True); self.entity_id_input.clear(); self._update_status("Highlight cleared.")
        
    def _update_status(self, msg, is_error=False):
        self.statusBar().showMessage(msg); self.statusBar().setStyleSheet(f"color: {'#f38ba8' if is_error else '#a6adc8'};")
        
    def _toggle_theme(self):
        app = QApplication.instance()
        self.is_dark_theme = not self.is_dark_theme
        app.setPalette(dark_palette if self.is_dark_theme else light_palette)
        app.setStyleSheet(DARK_STYLESHEET if self.is_dark_theme else "")

        bg_color = '#1e1e2e' if self.is_dark_theme else 'white'
        label_color_tuple = (0.804, 0.839, 0.957) if self.is_dark_theme else (0, 0, 0)

        self.plotter.set_background(bg_color)

        self.plotter.remove_axes()
        self.plotter.add_axes()
        self.axes_actor = self.plotter.renderer.axes_actor
        self._set_axes_label_color(label_color_tuple)

        self.theme_button.setText(f"Switch to {'Light' if self.is_dark_theme else 'Dark'} Mode")
        self._update_plot_visibility()


    def _create_load(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please open or generate a model first.")
            return

        dialog = CreateLoadDialog(self.current_generator.model, self)
        dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
        dialog.accepted.connect(self._on_load_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _create_constraint(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please open or generate a model first.")
            return

        dialog = CreateConstraintDialog(self.current_generator.model, self)
        # --- FIX: Connect directly to the handler, since the signal now carries the entity type ---
        dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
        dialog.accepted.connect(self._on_constraint_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    # --- Theme C: load and BC creators (PLOAD1/2, SPCD, MPC, RBAR, RBE1, RSPLINE, BOLT) ---

    def _require_model(self):
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please open or generate a model first.")
            return False
        return True

    def _create_pload1(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreatePload1Dialog
        from node_runner.commands import AddPload1Command
        dlg = CreatePload1Dialog(self)
        if not dlg.exec():
            return
        if not dlg.eids:
            QMessageBox.warning(self, "No EIDs", "Please specify at least one element ID.")
            return
        cmd = AddPload1Command(dlg.sid, dlg.eids, dlg.values)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Added PLOAD1 to {len(dlg.eids)} element(s) in SID {dlg.sid}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_pload2(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreatePload2Dialog
        from node_runner.commands import AddPload2Command
        dlg = CreatePload2Dialog(self)
        if not dlg.exec():
            return
        if not dlg.eids:
            QMessageBox.warning(self, "No EIDs", "Please specify at least one shell element ID.")
            return
        cmd = AddPload2Command(dlg.sid, dlg.eids, dlg.pressure_value)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Added PLOAD2 to {len(dlg.eids)} shell(s) in SID {dlg.sid}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_spcd(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateSpcdDialog
        from node_runner.commands import AddSpcdCommand
        dlg = CreateSpcdDialog(self)
        if not dlg.exec():
            return
        if not dlg.nids:
            QMessageBox.warning(self, "No nodes", "Please specify at least one node ID.")
            return
        cmd = AddSpcdCommand(dlg.sid, dlg.nids, dlg.dof_string, dlg.enforced_value)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Added SPCD on {len(dlg.nids)} node(s), DOFs {dlg.dof_string}, "
            f"value={dlg.enforced_value} in SID {dlg.sid}."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_mpc(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateMpcDialog
        from node_runner.commands import AddMpcCommand
        dlg = CreateMpcDialog(self)
        if not dlg.exec():
            return
        cmd = AddMpcCommand(dlg.sid, dlg.terms)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Added MPC SID {dlg.sid} with {len(dlg.terms)} terms."
        )
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_rbar(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateRbarDialog
        from node_runner.commands import AddRbarCommand
        dlg = CreateRbarDialog(self)
        if not dlg.exec():
            return
        cmd = AddRbarCommand(dlg.values)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Added RBAR {dlg.values['eid']}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_rbe1(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateRbe1Dialog
        from node_runner.commands import AddRbe1Command
        dlg = CreateRbe1Dialog(self)
        if not dlg.exec():
            return
        v = dlg.values
        if not v['indep_nodes'] or not v['dep_nodes']:
            QMessageBox.warning(self, "Missing nodes",
                                "RBE1 needs at least one independent and one dependent node.")
            return
        cmd = AddRbe1Command(v)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Added RBE1 {v['eid']}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_rspline(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateRsplineDialog
        from node_runner.commands import AddRsplineCommand
        dlg = CreateRsplineDialog(self)
        if not dlg.exec():
            return
        v = dlg.values
        if len(v['nodes']) < 2:
            QMessageBox.warning(self, "Need >= 2 nodes",
                                "RSPLINE requires at least two spline nodes.")
            return
        cmd = AddRsplineCommand(v)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(f"Added RSPLINE {v['eid']}.")
        self._update_viewer(self.current_generator, reset_camera=False)

    def _create_bolt(self):
        if not self._require_model():
            return
        from node_runner.dialogs import CreateBoltDialog
        from node_runner.commands import AddBoltCommand
        dlg = CreateBoltDialog(self)
        if not dlg.exec():
            return
        v = dlg.values
        cmd = AddBoltCommand(v)
        self.command_manager.execute(cmd, self.current_generator.model)
        self._update_status(
            f"Added BOLT {v['bid']} (preload={v['preload']}, EIDs={v['eids']})."
        )

    def _on_load_creation_accept(self):
        dialog = self.active_creation_dialog
        if not dialog: return

        params = dialog.get_parameters()
        if params:
            model = self.current_generator.model
            sid = params['sid']
            load_type = params.pop('type')

            # Always append to the load set - never delete-and-replace.
            # Use the Load Set Manager to edit/remove individual entries.
            cmd = AddLoadCommand(load_type, params)
            self.command_manager.execute(cmd, model)

            self._update_status(f"Added load entries to SID {sid}.")
            # v4.0.0 (B2): narrow LOADS refresh via the command's hint.
            self._dispatch_refresh(cmd.refresh_hint)

        self.active_creation_dialog = None

    def _on_constraint_creation_accept(self):
        dialog = self.active_creation_dialog
        if not dialog: return

        params = dialog.get_parameters()
        if params:
            model = self.current_generator.model
            sid = params['sid']

            # Always append to the constraint set - never delete-and-replace.
            cmd = AddConstraintCommand(params)
            self.command_manager.execute(cmd, model)

            self._update_status(f"Added constraint entries to SID {sid}.")
            # v4.0.0 (B2): use the command's narrow refresh hint
            # (CONSTRAINTS) instead of the pre-v4.0.0 _populate_tree
            # blanket rebuild. On large production decks this drops the
            # post-edit pause from multi-second tree rebuilds to a
            # near-instant constraint-actors-only refresh, and B1's
            # tree snapshot/restore preserves the user's checkbox state.
            self._dispatch_refresh(cmd.refresh_hint)

        self.active_creation_dialog = None

    def _open_load_set_manager(self, sid):
        """Open the Load Set Manager dialog for viewing/editing entries in a set."""
        if not self.current_generator:
            return
        model = self.current_generator.model
        if sid not in model.loads and sid not in model.tempds:
            QMessageBox.information(self, "Empty Set",
                                    f"Load set SID {sid} has no entries.")
            return
        from node_runner.dialogs.loads import LoadSetManagerDialog
        dialog = LoadSetManagerDialog(model, sid, self)
        if dialog.exec():
            new_entries = dialog.get_modified_entries()
            tempd_removed = dialog.tempd_was_removed()
            if new_entries is not None:
                cmd = ReplaceLoadSetCommand(sid, new_entries)
                self.command_manager.execute(cmd, model)
            if tempd_removed:
                from node_runner.commands import DeleteTempdCommand
                cmd2 = DeleteTempdCommand(sid)
                self.command_manager.execute(cmd2, model)
            if new_entries is not None or tempd_removed:
                self._populate_tree()
                self._populate_loads_tab()
                self._create_all_load_actors()
                self._update_plot_visibility()
                n = len(new_entries) if new_entries is not None else 0
                self._update_status(
                    f"Updated load set SID {sid} ({n} entries)")

    def _add_to_load_set(self, sid):
        """Open CreateLoadDialog pre-filled with the given SID to add entries."""
        if not self.current_generator:
            return
        if self.active_creation_dialog:
            self.active_creation_dialog.activateWindow()
            return
        dialog = CreateLoadDialog(self.current_generator.model, self,
                                  existing_sid=sid)
        dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
        dialog.accepted.connect(self._on_load_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _handle_sub_dialog_selection_request(self, entity_type):
        if self.active_creation_dialog:
            self.active_creation_dialog.hide()

        if entity_type == 'Node':
            all_ids = self.current_generator.model.nodes.keys()
        else:
            all_ids = [eid for eid, elem in self.current_generator.model.elements.items() if elem.type in ['CQUAD4', 'CTRIA3']]

        select_by = self._build_select_by_data() if entity_type == 'Element' else None

        # Store entity_type for the closure
        _et = entity_type

        def on_accept():
            selected_ids = self.selection_bar.get_selected_ids()
            if self.active_creation_dialog:
                self.active_creation_dialog.update_selection_list(_et, selected_ids)
                self.active_creation_dialog.show()
            self._end_selection_mode()

        self._start_selection(entity_type, all_ids, on_accept, select_by)
        # Override reject to also re-show the creation dialog
        self._selection_reject_restore_creation = True

    def _compute_grid_centers_and_normals(self, grid):
        """v3.4.1: precompute cell centers + shell normals once per grid
        in fully vectorized numpy.

        Before this, _create_all_load_actors did per-PLOAD4 cell access
        (``grid.get_cell(idx).points`` + ``np.cross`` + normalize), which
        wraps a VTK cell in a pyvista.Cell object on every iteration.
        On a large pressure-load production deck (~600k+ PLOAD4 face refs)
        that loop blocked the main thread long enough for Windows to
        mark the app "Not Responding".

        Returns ``(centers, normals)`` as ``(n_cells, 3)`` float64
        numpy arrays. Shell-face normals (QUAD + TRIANGLE) are computed
        via a single batched ``np.cross`` over the first 3 vertices.
        Non-shell cells get the +Z fallback so PLOAD4-on-solid arrows
        still get *some* direction (the magnitude/coloring is correct
        even if the direction is approximate).
        """
        if grid is None or grid.n_cells == 0:
            return (np.zeros((0, 3)), np.zeros((0, 3)))
        n = int(grid.n_cells)
        # Bulk VTK call - O(n_cells) but inside C++, not Python.
        centers = np.asarray(grid.cell_centers().points)

        normals = np.zeros((n, 3), dtype=np.float64)
        normals[:, 2] = 1.0  # default
        try:
            celltypes = np.asarray(grid.celltypes)
            conn = np.asarray(grid.cell_connectivity)
            offsets = np.asarray(grid.offset)
            points = np.asarray(grid.points)
        except Exception:
            return centers, normals
        QUAD = int(pv.CellType.QUAD)
        TRI = int(pv.CellType.TRIANGLE)
        shell_mask = (celltypes == QUAD) | (celltypes == TRI)
        shell_idx = np.where(shell_mask)[0]
        if shell_idx.size:
            # First 3 vertex indices for each shell cell.
            p0_pi = conn[offsets[shell_idx]]
            p1_pi = conn[offsets[shell_idx] + 1]
            p2_pi = conn[offsets[shell_idx] + 2]
            p0 = points[p0_pi]
            p1 = points[p1_pi]
            p2 = points[p2_pi]
            cross = np.cross(p1 - p0, p2 - p0)
            norm = np.linalg.norm(cross, axis=1)
            mask = norm > 1e-12
            cross[mask] /= norm[mask, None]
            cross[~mask] = np.array([0.0, 0.0, 1.0])
            normals[shell_idx] = cross
        return centers, normals

    def _build_load_glyphs_for_sid(self, sid):
        """v4.0.9: build force/moment/pressure glyph actors for ONE
        load SID. Extracted from the load_sid_loop inside
        ``_update_plot_visibility`` so the per-SID fast path in
        ``_handle_tree_item_changed`` can build glyphs for a single
        SID without iterating all 1,962. The lazily-built per-SID
        data actors (``force_actors_{sid}`` etc.) must already exist.
        callers are responsible for calling
        ``_ensure_load_actor_for_sid(sid)`` first.
        """
        if not self.current_grid:
            return
        try:
            arrow_size_percent = float(self.arrow_scale_input.text())
            use_relative_scaling = self.relative_scaling_check.isChecked()
        except (ValueError, AttributeError):
            arrow_size_percent = 10.0
            use_relative_scaling = True
        base_scale = self.current_grid.length * (arrow_size_percent / 100.0)
        max_force_mag = self.load_scaling_info.get('force_max_mag', 0.0)
        max_moment_mag = self.load_scaling_info.get('moment_max_mag', 0.0)
        max_pressure_mag = self.load_scaling_info.get('pressure_max_mag', 0.0)

        from node_runner.profiling import perf_event
        for kind, max_mag in [('force', max_force_mag),
                              ('moment', max_moment_mag),
                              ('pressure', max_pressure_mag)]:
            actor_name = f"{kind}_actors_{sid}"
            glyph_name = f"{kind}_glyphs_{sid}"
            if glyph_name in self.plotter.actors:
                self.plotter.remove_actor(glyph_name, render=False)
            actor = self.plotter.actors.get(actor_name)
            if not actor:
                # v4.0.10: distinguish "no lazy actor for this SID" from
                # "lazy actor exists but dataset empty". Both produce
                # no glyphs, but the diagnostic story is different.
                perf_event('glyph', f'{kind}_no_actor', sid=sid)
                continue
            dataset = actor.mapper.dataset
            if not dataset.n_points:
                # v4.0.10: this SID has no cards of this kind (e.g.
                # SID with only forces → moment dataset empty). Quiet
                # signal so users can tell "moments not built because
                # SID has no moment cards" from "build failed".
                perf_event('glyph', f'{kind}_empty_dataset', sid=sid)
                continue
            if kind == 'force':
                vectors = dataset['vectors']
                mags = np.linalg.norm(vectors, axis=1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                if use_relative_scaling and max_mag > 1e-9:
                    scaled_vectors = unit_vectors * (mags / max_mag * base_scale)[:, np.newaxis]
                else:
                    scaled_vectors = unit_vectors * base_scale
                self.plotter.add_arrows(dataset.points, scaled_vectors,
                                        color='red', name=glyph_name,
                                        reset_camera=False)
                # v4.0.10: settle "did N arrows actually land on screen?".
                perf_event('glyph', 'force_added',
                           sid=sid, n_arrows=int(dataset.n_points))
            elif kind == 'moment':
                # v4.0.11: right-hand-rule-respecting composite moment
                # glyph (single arrow + curl + curl arrowhead). Cached
                # via _get_moment_glyph_geom.
                # v4.0.12 (Bug 4 fix): two-tier fallback. The v4.0.11
                # a production-deck crash was a VTK SEGFAULT in vtkGlyph3D when the
                # composite resolved to an UnstructuredGrid (mixed
                # cell types). We attempt the composite first; on any
                # Python-side exception OR when the composite returns
                # None, retry with a plain single-tipped pv.Arrow. If
                # THAT fails too, log and skip this SID's moment
                # glyph rather than crashing the whole app.
                def _try_build_moment_glyph(geom_polydata):
                    """Inner: do the glyph instance + add_mesh. May raise."""
                    if use_relative_scaling and max_mag > 1e-9:
                        mags = np.linalg.norm(dataset['vectors'], axis=1)
                        dataset['relative_mag'] = mags / max_mag
                        out_glyphs = dataset.glyph(orient='vectors',
                                                   scale='relative_mag',
                                                   factor=base_scale,
                                                   geom=geom_polydata)
                    else:
                        out_glyphs = dataset.glyph(orient='vectors',
                                                   scale=False,
                                                   factor=base_scale,
                                                   geom=geom_polydata)
                    self.plotter.add_mesh(out_glyphs, color='cyan',
                                          name=glyph_name,
                                          reset_camera=False)
                    return out_glyphs
                glyphs = None
                composite_attempted = False
                try:
                    composite_geom = self._get_moment_glyph_geom()
                    if composite_geom is not None:
                        composite_attempted = True
                        glyphs = _try_build_moment_glyph(composite_geom)
                except Exception as _exc:
                    perf_event('glyph', 'moment_glyph_failed', sid=sid,
                               reason='composite_attempt',
                               exc=f"{_exc.__class__.__name__}: {str(_exc)[:120]}")
                    glyphs = None
                if glyphs is None:
                    # Fallback: single-tipped pv.Arrow. Still right-hand
                    # rule by convention (vector direction = thumb),
                    # just without the curl decoration.
                    try:
                        fallback_arrow = pv.Arrow(
                            start=(0.0, 0.0, 0.0),
                            direction=(1.0, 0.0, 0.0),
                            tip_length=0.2, tip_radius=0.08,
                            shaft_radius=0.03)
                        glyphs = _try_build_moment_glyph(fallback_arrow)
                        perf_event('glyph', 'moment_fallback_arrow',
                                   sid=sid,
                                   composite_attempted=composite_attempted)
                    except Exception as _exc:
                        perf_event('glyph', 'moment_build_failed', sid=sid,
                                   reason='fallback_arrow_too',
                                   exc=f"{_exc.__class__.__name__}: {str(_exc)[:120]}")
                        glyphs = None
                if glyphs is not None:
                    try:
                        n_glyph_cells = int(glyphs.n_cells)
                    except Exception:
                        n_glyph_cells = -1
                    perf_event('glyph', 'moment_added',
                               sid=sid,
                               n_arrows=int(dataset.n_points),
                               n_glyph_cells=n_glyph_cells)
            elif kind == 'pressure':
                vectors = dataset['vectors']
                mags = np.linalg.norm(vectors, axis=1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                uniform_size = base_scale * 0.75
                if use_relative_scaling and max_mag > 1e-9:
                    scaled_vectors = unit_vectors * (mags / max_mag * uniform_size)[:, np.newaxis]
                else:
                    scaled_vectors = unit_vectors * uniform_size
                self.plotter.add_arrows(dataset.points, scaled_vectors,
                                        color='yellow', name=glyph_name,
                                        reset_camera=False)
                perf_event('glyph', 'pressure_added',
                           sid=sid, n_arrows=int(dataset.n_points))

    def _remove_load_glyphs_for_sid(self, sid):
        """v4.0.9: remove the 3 glyph actors for a single load SID."""
        for prefix in ('force_glyphs_', 'moment_glyphs_', 'pressure_glyphs_'):
            name = f"{prefix}{sid}"
            if name in self.plotter.actors:
                self.plotter.remove_actor(name, render=False)

    def _build_constraint_glyphs_for_sid(self, sid):
        """v4.0.9: build the constraint glyph actor for ONE SPC SID.
        Extracted from the spc_sid_loop inside
        ``_update_plot_visibility``.
        """
        from node_runner.profiling import perf_event
        if not self.current_grid:
            return
        glyph_name = f"constraint_glyphs_{sid}"
        if glyph_name in self.plotter.actors:
            self.plotter.remove_actor(glyph_name, render=False)
        actor = self.plotter.actors.get(f"constraint_actors_{sid}")
        if not actor:
            perf_event('glyph', 'constraint_no_actor', sid=sid)
            return
        dataset = actor.mapper.dataset
        if not dataset.n_points:
            perf_event('glyph', 'constraint_empty_dataset', sid=sid)
            return
        glyph_geom = pv.Cone(direction=(0, 0, 1), height=1.0, radius=0.3)
        scale_factor = self.current_grid.length * 0.015
        glyphs = dataset.glyph(scale=False, factor=scale_factor,
                               orient=False, geom=glyph_geom)
        self.plotter.add_mesh(glyphs, color='cyan', name=glyph_name,
                              reset_camera=False)
        try:
            n_glyph_cells = int(glyphs.n_cells)
        except Exception:
            n_glyph_cells = -1
        perf_event('glyph', 'constraint_added',
                   sid=sid,
                   n_arrows=int(dataset.n_points),
                   n_glyph_cells=n_glyph_cells)

    def _remove_constraint_glyphs_for_sid(self, sid):
        """v4.0.9: remove the constraint glyph actor for one SPC SID."""
        name = f"constraint_glyphs_{sid}"
        if name in self.plotter.actors:
            self.plotter.remove_actor(name, render=False)

    def _fast_path_sid_toggle(self, kind, sid, check_state):
        """v4.0.9: O(1) per-SID visibility update. Returns True on
        success so the caller can skip the full
        ``_update_plot_visibility`` cycle (~4-5 s on a large production deck). Returns
        False on any failure. caller MUST fall through to the full
        cycle so the user never sees a stale scene.

        Safe to skip the full cycle because load_set / constraint_set
        checkbox toggles don't affect any other scene element (main
        mesh, coord systems, RBE spiders, etc.). Only the per-SID
        glyphs change.
        """
        from node_runner.profiling import perf_event
        if not self.current_generator or not self.current_grid:
            return False
        is_checked = (check_state == QtCore.Qt.Checked)
        try:
            if kind == 'load_set':
                if is_checked:
                    self._ensure_load_actor_for_sid(sid)
                    self._build_load_glyphs_for_sid(sid)
                    self._visible_load_sids.add(sid)
                else:
                    self._remove_load_glyphs_for_sid(sid)
                    self._visible_load_sids.discard(sid)
            elif kind == 'constraint_set':
                if is_checked:
                    self._build_constraint_glyphs_for_sid(sid)
                    self._visible_spc_sids.add(sid)
                else:
                    self._remove_constraint_glyphs_for_sid(sid)
                    self._visible_spc_sids.discard(sid)
            else:
                return False
        except Exception as _exc:
            perf_event('fast_path', f'{kind}_toggle_failed', sid=sid,
                       exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            return False
        perf_event('fast_path', f'{kind}_toggle',
                   sid=sid, checked=is_checked,
                   n_visible_load_sids=len(self._visible_load_sids),
                   n_visible_spc_sids=len(self._visible_spc_sids))
        self.plotter.render()
        return True

    # v4.0.14 (Fix C v2): pre-built per-Nastran-type actors.
    # Redesigned after v4.0.11 was reverted in v4.0.12. Key
    # differences from v4.0.11:
    #   1. Build from current_grid (full mesh). not LOD'd display_grid.
    #      Restores full visual density on Plates and brings RBE2/RBE3
    #      into the dict so the Rigid spider can toggle correctly.
    #   2. RBE spider visibility derived from the Rigid checkbox state
    #      directly (not target_states.get). defensive against
    #      `_category_actors` not containing every ntype.
    #   3. Same filter-mode coexistence pattern as v4.0.11 (legacy
    #      _rebuild_plot path triggered when PID/MID/isolate active;
    #      pre-built actors hidden in legacy mode).

    def _build_pre_built_category_actors(self):
        """v4.0.14: build one named PyVista actor per unique Nastran
        element type present in ``self.current_grid`` (the full,
        non-LOD'd mesh). Category checkbox toggles in the model tree
        (Plates, Beams, Solids, etc.) become ``actor.SetVisibility``
        calls. sub-100 ms on a large production deck. bypassing the legacy
        ``_rebuild_plot`` extract_cells + split + compute_normals
        chain.

        PID/MID/isolate/hidden filtering still routes through legacy
        ``_rebuild_plot``; when it runs, pre-built actors are hidden
        to avoid double-rendering. When the filter clears, the
        category fast-path restores pre-built visibility.
        """
        from node_runner.profiling import perf_stage, perf_event
        if self.current_grid is None:
            return
        if 'type' not in self.current_grid.cell_data:
            perf_event('actors', 'build_per_category_skipped',
                       reason='no_type_cell_data')
            return
        # v4.0.15: gate by color_mode. Modes that need external data
        # (quality / results) route through legacy `_rebuild_plot`
        # so its existing scalars-from-OP2 / scalars-from-quality
        # pipelines stay intact. The set is module-level + checked
        # here so adding a new mode is a single-line addition there
        # plus a branch below.
        if self.color_mode not in _PRE_BUILT_SUPPORTED_COLOR_MODES:
            perf_event('actors', 'build_per_category_skipped',
                       reason=f'color_mode_{self.color_mode}_unsupported')
            return
        # v4.0.15: capture color_mode at build time. `_set_coloring_mode`
        # tears these actors down and rebuilds on mode change so the
        # cached coloring stays in sync with self.color_mode.
        build_color_mode = self.color_mode
        with perf_stage('actors', 'build_per_category'):
            types_arr = np.asarray(self.current_grid.cell_data['type'])
            unique_types = np.unique(types_arr)
            shell_color = self.type_color_map.get('Shells', '#cccccc')
            beam_color = self.type_color_map.get('Beams', '#88c0d0')
            show_edges = (self.render_style == 'surface')
            for ntype in unique_types:
                if not isinstance(ntype, str) or ntype in ('VTK_VERTEX',):
                    continue
                cell_idx = np.where(types_arr == ntype)[0]
                if cell_idx.size == 0:
                    continue
                try:
                    with perf_stage('actors', 'add_one_category',
                                    ntype=ntype, n_cells=int(cell_idx.size)):
                        slice_grid = self.current_grid.extract_cells(cell_idx)
                        is_shell_arr = slice_grid.cell_data.get('is_shell')
                        is_shell = (
                            is_shell_arr is not None
                            and len(is_shell_arr)
                            and bool(np.all(np.asarray(is_shell_arr) == 1)))
                        if is_shell:
                            try:
                                slice_render = (slice_grid.extract_surface()
                                                .compute_normals(
                                    cell_normals=False, point_normals=True,
                                    consistent_normals=True, inplace=False))
                            except Exception:
                                slice_render = slice_grid
                            shape_kw = dict(
                                show_edges=show_edges,
                                opacity=self.shell_opacity)
                        else:
                            slice_render = slice_grid
                            shape_kw = dict(render_lines_as_tubes=True,
                                            line_width=self.beam_width)
                        # v4.0.15: per-mode color kwargs.
                        # v4.0.17: 'property' mode now uses per-cell
                        # RGB via _attach_per_cell_rgb_for_property_mode
                        # instead of the categorical LUT (which
                        # rendered the same PID with different colors
                        # across actors because PyVista's
                        # `categories=True` is range-positional, not
                        # value-keyed). Same PID, same color, every
                        # render. across all actors and all rebuilds.
                        if (build_color_mode == 'property'
                                and 'PID' in slice_render.cell_data):
                            self._attach_per_cell_rgb_for_property_mode(
                                slice_render)
                            color_kw = dict(
                                scalars='cell_rgb', rgb=True,
                                show_scalar_bar=False)
                        else:
                            # 'type' mode (or any safe-uniform fallback):
                            # one color per shell/beam family.
                            color_kw = dict(
                                color=shell_color if is_shell else beam_color)
                        actor = self.plotter.add_mesh(
                            slice_render,
                            name=f'cat_{ntype}',
                            pickable=False, reset_camera=False,
                            **color_kw, **shape_kw)
                        self._category_actors[ntype] = actor
                        # v5.0.0 item 11a: attach per-cell vtkGhostType
                        # for sub-100 ms PID/MID toggles. Initialized
                        # all-visible; mutated in-place by the new
                        # fast-path PID/MID toggles.
                        self._category_actor_grids[ntype] = slice_render
                        try:
                            n_cells = int(slice_render.n_cells)
                            if n_cells:
                                slice_render.cell_data['vtkGhostType'] = (
                                    np.zeros(n_cells, dtype=np.uint8))
                        except Exception:
                            pass
                        perf_event('actors', 'category_added',
                                   ntype=ntype,
                                   n_cells=int(cell_idx.size),
                                   is_shell=is_shell,
                                   mode=build_color_mode)
                except Exception as _exc:
                    perf_event('actors', 'add_one_category_failed',
                               ntype=ntype,
                               exc=f"{type(_exc).__name__}: {str(_exc)[:120]}")
            self._category_visible_ntypes = set(self._category_actors.keys())
            self._using_legacy_mesh = False
            # v5.0.0 item 11a: replay any pre-existing PID/MID hidden
            # state (e.g. after a color-mode change that tore down and
            # rebuilt the actors). Applies in-place; no-op when both
            # sets are empty.
            try:
                self._apply_pid_mid_visibility_to_actors()
            except Exception:
                pass
            try:
                from node_runner.profiling import perf_event
                perf_event('pid_visibility', 'cache_init',
                           n_pids_hidden=len(self._hidden_pids),
                           n_mids_hidden=len(self._hidden_mids),
                           n_actors=len(self._category_actor_grids))
            except Exception:
                pass

    def _has_active_cell_filter(self):
        """v4.0.14: return True if any tree-state setting forces
        cell-level filtering routed through legacy ``_rebuild_plot``.

        v5.0.0 item 11a: PID and MID checkboxes are NO LONGER counted
        as "active cell filter" - they're handled in-place via the per-
        cell ``vtkGhostType`` array on each ``cat_<ntype>`` actor's
        grid, which keeps the pre-built path engaged. The remaining
        cell-level filters (isolate mode, hidden groups) still force
        legacy.
        """
        if getattr(self, '_isolate_mode', None):
            return True
        if getattr(self, '_hidden_groups', None):
            return True
        return False

    def _compute_category_actor_visibility(self, ntype):
        """v4.0.14: compound (by_type ∧ by_shape) visibility for the
        per-Nastran-type actor with name ``cat_<ntype>``. Returns
        True/False if resolvable from category state alone, or None
        if cell-level filtering is active.
        """
        type_cat = _NTYPE_TO_TYPE_CATEGORY.get(ntype)
        shape_cat = _NTYPE_TO_SHAPE_CATEGORY.get(ntype)
        if type_cat:
            items = self._find_tree_items(('elem_by_type_group', type_cat))
            if items and items[0].checkState(0) != QtCore.Qt.Checked:
                return False
        if shape_cat:
            items = self._find_tree_items(('elem_shape_group', shape_cat))
            if items and items[0].checkState(0) != QtCore.Qt.Checked:
                return False
        if self._has_active_cell_filter():
            return None
        return True

    def _fast_path_category_toggle(self, group_kind, group_name):
        """v4.0.14: O(N_categories) visibility update for a category
        (by_type or by_shape) checkbox click. Returns True if all
        per-type actors resolve from category state alone; False if
        any falls into cell-level filter territory. caller falls
        through to legacy ``_rebuild_plot``.

        Differs from v4.0.11: rbe_actors visibility is derived from
        the Rigid checkbox state DIRECTLY (not from target_states),
        defensively guarding against `_category_actors` missing
        RBE2/RBE3 entries (e.g., a deck with no rigids).
        """
        from node_runner.profiling import perf_event
        if not self._category_actors:
            return False
        target_states = {}
        new_visible = set()
        for ntype in self._category_actors:
            vis = self._compute_category_actor_visibility(ntype)
            if vis is None:
                perf_event('fast_path', 'category_toggle_filter_active',
                           group_kind=group_kind, group_name=group_name,
                           ntype=ntype)
                return False
            target_states[ntype] = bool(vis)
            if vis:
                new_visible.add(ntype)
        n_changed = 0
        for ntype, actor in self._category_actors.items():
            new_vis = target_states[ntype]
            if bool(actor.GetVisibility()) != new_vis:
                actor.SetVisibility(new_vis)
                n_changed += 1
        # Hide any legacy shells/beams actors that may remain from a
        # prior filter cycle. v4.0.15 fix: do NOT hide nodes_actor.
        # it's the user-controlled Nodes display, not a filter-cycle
        # leftover. The v4.0.14 version was hiding it on every
        # category toggle so the first click made nodes disappear
        # and subsequent clicks left them hidden.
        for legacy_name in ('shells', 'beams'):
            if legacy_name in self.plotter.actors:
                la = self.plotter.actors[legacy_name]
                if bool(la.GetVisibility()):
                    la.SetVisibility(False)
                    n_changed += 1
        # v4.0.14: rbe_actors (RBE spider visualization) visibility is
        # derived from the Rigid type+shape checkbox state DIRECTLY,
        # not from target_states. The v4.0.11 bug was that the LOD
        # path dropped RBE2/RBE3 from _category_actors entirely, so
        # `target_states.get('RBE2', True)` returned True and the
        # spider stayed visible. Defensive query of the tree avoids
        # this regardless of whether cat_RBE2/RBE3 exist.
        rigid_items = self._find_tree_items(('elem_by_type_group', 'Rigid'))
        rigid_checked = bool(
            rigid_items and rigid_items[0].checkState(0) == QtCore.Qt.Checked)
        shape_rigid_items = self._find_tree_items(('elem_shape_group', 'Rigid'))
        shape_rigid_checked = bool(
            not shape_rigid_items
            or shape_rigid_items[0].checkState(0) == QtCore.Qt.Checked)
        rbe_should_show = rigid_checked and shape_rigid_checked
        rbe_actor = self.plotter.actors.get('rbe_actors')
        if rbe_actor is not None:
            if bool(rbe_actor.GetVisibility()) != rbe_should_show:
                rbe_actor.SetVisibility(rbe_should_show)
                n_changed += 1
        self._category_visible_ntypes = new_visible
        self._using_legacy_mesh = False
        perf_event('fast_path', 'category_toggle',
                   group_kind=group_kind, group_name=group_name,
                   n_actors_changed=n_changed,
                   n_visible_categories=len(new_visible))
        self.plotter.render()
        return True

    # ─── v5.0.0 item 11a: sub-100 ms PID / MID toggles ──────────────────
    def _pid_to_mid_map(self):
        """Build (and cache) the {pid: mid} map used by MID toggles.

        Cleared on every viewer rebuild via ``current_generator`` swap.
        """
        cache = getattr(self, '_pid_to_mid_cache', None)
        cur_gen = getattr(self, 'current_generator', None)
        cache_gen = getattr(self, '_pid_to_mid_cache_gen', None)
        if cache is not None and cache_gen is cur_gen:
            return cache
        out: dict[int, int] = {}
        try:
            model = cur_gen.model if cur_gen is not None else None
            for pid, prop in (model.properties or {}).items():
                mid = getattr(prop, 'mid', None)
                if mid is None:
                    mid = getattr(prop, 'mid1', None)
                if mid is not None:
                    out[int(pid)] = int(mid)
        except Exception:
            pass
        self._pid_to_mid_cache = out
        self._pid_to_mid_cache_gen = cur_gen
        return out

    def _apply_pid_mid_visibility_to_actors(self):
        """Refresh every cat_<ntype> actor's vtkGhostType cell-data from
        the current ``_hidden_pids`` and ``_hidden_mids`` sets.

        Cheap: numpy isin over each actor's cell PID array. Sub-100 ms
        on a deck with 2.4M cells. Marks each grid Modified() so the mapper
        re-picks up the array on next render.
        """
        if not self._category_actor_grids:
            return 0
        hidden_pids = self._hidden_pids
        hidden_mids = self._hidden_mids
        # Resolve hidden MIDs to the corresponding PIDs via property->material.
        if hidden_mids:
            pid_to_mid = self._pid_to_mid_map()
            hidden_via_mid = {
                pid for pid, mid in pid_to_mid.items() if mid in hidden_mids
            }
        else:
            hidden_via_mid = set()
        hidden_combined = hidden_pids | hidden_via_mid
        hidden_arr = (np.fromiter(hidden_combined, dtype=np.int64)
                      if hidden_combined else None)
        n_cells_total = 0
        n_hidden_cells = 0
        for ntype, grid in self._category_actor_grids.items():
            try:
                if 'PID' not in grid.cell_data:
                    continue
                n_cells = int(grid.n_cells)
                if n_cells == 0:
                    continue
                pids = np.asarray(grid.cell_data['PID'])
                if hidden_arr is None or hidden_arr.size == 0:
                    ghost = np.zeros(n_cells, dtype=np.uint8)
                else:
                    mask = np.isin(pids, hidden_arr)
                    ghost = np.where(mask,
                                     np.uint8(0x01),  # HIDDENCELL
                                     np.uint8(0)).astype(np.uint8)
                    n_hidden_cells += int(mask.sum())
                grid.cell_data['vtkGhostType'] = ghost
                # Force VTK to re-read the array.
                try:
                    grid.Modified()
                except Exception:
                    pass
                n_cells_total += n_cells
            except Exception:
                continue
        return (n_cells_total, n_hidden_cells)

    def _fast_path_pid_toggle(self, pid, checked):
        """v5.0.0 item 11a: in-place PID visibility flip on pre-built
        actors. Sub-100 ms on a large production deck.
        """
        from node_runner.profiling import perf_event
        import time as _time
        if not self._category_actor_grids:
            return False
        # Isolate / hidden-groups still force legacy.
        if (getattr(self, '_isolate_mode', None)
                or getattr(self, '_hidden_groups', None)):
            return False
        try:
            pid_int = int(pid)
        except Exception:
            return False
        if checked:
            self._hidden_pids.discard(pid_int)
        else:
            self._hidden_pids.add(pid_int)
        t0 = _time.perf_counter()
        try:
            res = self._apply_pid_mid_visibility_to_actors()
            n_total, n_hidden = (res if isinstance(res, tuple) else (0, 0))
        except Exception:
            n_total, n_hidden = 0, 0
        self.plotter.render()
        perf_event('fast_path', 'pid_toggle',
                   wall_s=round(_time.perf_counter() - t0, 4),
                   pid=pid_int, checked=bool(checked),
                   n_hidden_pids=len(self._hidden_pids),
                   n_cells_total=int(n_total),
                   n_cells_hidden=int(n_hidden))
        return True

    def _fast_path_mid_toggle(self, mid, checked):
        """v5.0.0 item 11a: in-place MID visibility flip on pre-built
        actors. Sub-100 ms on a large production deck. Material visibility propagates to
        every property that references it via _pid_to_mid_map().
        """
        from node_runner.profiling import perf_event
        import time as _time
        if not self._category_actor_grids:
            return False
        if (getattr(self, '_isolate_mode', None)
                or getattr(self, '_hidden_groups', None)):
            return False
        try:
            mid_int = int(mid)
        except Exception:
            return False
        if checked:
            self._hidden_mids.discard(mid_int)
        else:
            self._hidden_mids.add(mid_int)
        t0 = _time.perf_counter()
        try:
            res = self._apply_pid_mid_visibility_to_actors()
            n_total, n_hidden = (res if isinstance(res, tuple) else (0, 0))
        except Exception:
            n_total, n_hidden = 0, 0
        self.plotter.render()
        perf_event('fast_path', 'mid_toggle',
                   wall_s=round(_time.perf_counter() - t0, 4),
                   mid=mid_int, checked=bool(checked),
                   n_hidden_mids=len(self._hidden_mids),
                   n_cells_total=int(n_total),
                   n_cells_hidden=int(n_hidden))
        return True

    def _refresh_visible_sid_sets(self):
        """v4.0.10: rebuild ``_visible_load_sids`` and
        ``_visible_spc_sids`` from the current loads_tab tree state.

        Called at:
          - Top of every ``_update_plot_visibility`` cycle (defensive
            refresh in case the fast-path-maintained sets drift).
          - End of ``_populate_loads_tab`` (because the tree rebuild
            runs under blockSignals, no fast-path events fire during
            it, so the sets need an explicit refresh).

        Cost: O(n_load_sids + n_spc_sids) dict lookups + checkState
        calls. ~5-10 ms on a deck with 1,962 load SIDs. Far cheaper than
        the visibility cycle it enables.
        """
        if not self.current_generator:
            self._visible_load_sids = set()
            self._visible_spc_sids = set()
            return
        model = self.current_generator.model
        li = getattr(self, '_loads_tab_item_index', None) or {}
        new_load = set()
        new_spc = set()
        if model.loads:
            for sid in model.loads.keys():
                items = li.get(('load_set', sid))
                if items and items[0].checkState(0) == QtCore.Qt.Checked:
                    new_load.add(sid)
        if model.spcs:
            for sid in model.spcs.keys():
                items = li.get(('constraint_set', sid))
                if items and items[0].checkState(0) == QtCore.Qt.Checked:
                    new_spc.add(sid)
        self._visible_load_sids = new_load
        self._visible_spc_sids = new_spc

    def _create_all_load_actors(self):
        """v3.4.2: scaling-info pass + cache prep ONLY. Per-SID actors
        are built on demand by ``_ensure_load_actor_for_sid`` when the
        visibility-update path actually needs them.

        On a 1,962-SID production deck the prior
        per-SID `add_mesh` loop hit Qt/VTK at ~50ms each = ~100 s of
        scene build, even though every actor was created with
        opacity=0 (invisible until a checkbox toggled it). Now we
        delay actor creation until the user actually shows a SID,
        which in practice is 1-5 SIDs at a time.
        """
        # Clear any previously-built actors AND invalidate our caches.
        actor_names_to_remove = [
            name for name in self.plotter.actors
            if name.startswith(('force_actors_', 'moment_actors_',
                                'pressure_actors_'))
        ]
        self.plotter.remove_actor(actor_names_to_remove, render=False)
        self.load_scaling_info.clear()
        self._lazy_load_centers = None
        self._lazy_load_normals = None
        self._lazy_load_eid_to_idx = None
        self._lazy_load_built_sids = set()

        if not self.current_generator or not self.current_grid:
            return
        model = self.current_generator.model
        if not model.loads:
            return

        # ---- Scaling pass (cheap; bounded by total load cards). ----
        all_force_mags, all_moment_mags, all_pressure_mags = [], [], []
        for load_list in model.loads.values():
            for load in load_list:
                t = load.type
                if t == 'FORCE':
                    all_force_mags.append(
                        np.linalg.norm(np.array(load.xyz) * load.mag))
                elif t == 'MOMENT':
                    all_moment_mags.append(
                        np.linalg.norm(np.array(load.xyz) * load.mag))
                elif t == 'PLOAD4':
                    all_pressure_mags.append(abs(load.pressures[0]))
                elif t == 'PLOAD2':
                    # v5.2.2 issue 1: same magnitude treatment as PLOAD4.
                    all_pressure_mags.append(
                        abs(float(getattr(load, 'pressure', 0.0))))
        self.load_scaling_info['force_max_mag'] = (
            float(np.max(all_force_mags)) if all_force_mags else 0.0)
        self.load_scaling_info['moment_max_mag'] = (
            float(np.max(all_moment_mags)) if all_moment_mags else 0.0)
        self.load_scaling_info['pressure_max_mag'] = (
            float(np.max(all_pressure_mags)) if all_pressure_mags else 0.0)
        # Deliberately NO per-SID actor creation here. The visibility
        # update path calls _ensure_load_actor_for_sid(sid) for any
        # SID whose checkbox is checked.

    def _ensure_load_actor_for_sid(self, sid):
        """v3.4.2: lazily build the FORCE / MOMENT / PRESSURE actors
        for ONE load SID. Idempotent - returns immediately if any
        actor for this SID is already in the plotter.

        Centers/normals/eid-lookup caches are populated on first call
        and reused for subsequent SIDs in the same grid lifetime.

        v4.0.10: catchers on every branch so the next log proves
        whether the lazy actor build path is firing as expected.
        """
        from node_runner.profiling import perf_event
        if sid in self._lazy_load_built_sids:
            perf_event('lazy_load_actor', 'already_built', sid=sid)
            return
        if not self.current_generator or not self.current_grid:
            perf_event('lazy_load_actor', 'skipped_no_grid', sid=sid)
            return
        grid = self.current_grid
        if 'EID' not in grid.cell_data:
            perf_event('lazy_load_actor', 'skipped_no_eid', sid=sid)
            return
        model = self.current_generator.model
        load_list = (model.loads or {}).get(sid)
        if not load_list:
            self._lazy_load_built_sids.add(sid)
            perf_event('lazy_load_actor', 'skipped_no_loads', sid=sid)
            return

        # First-call cache: compute centers + normals + EID lookup once.
        if self._lazy_load_centers is None:
            (self._lazy_load_centers,
             self._lazy_load_normals) = self._compute_grid_centers_and_normals(grid)
        if self._lazy_load_eid_to_idx is None:
            self._lazy_load_eid_to_idx = {
                int(e): int(i)
                for i, e in enumerate(np.asarray(grid.cell_data['EID']))
            }
        cell_centers = self._lazy_load_centers
        cell_normals = self._lazy_load_normals
        eid_to_idx = self._lazy_load_eid_to_idx

        f_pts, f_vecs = [], []
        m_pts, m_vecs = [], []
        p_pts, p_vecs = [], []
        nodes = model.nodes

        for load in load_list:
            t = load.type
            if t == 'FORCE':
                nid = load.node
                if nid in nodes:
                    f_pts.append(nodes[nid].get_position())
                    f_vecs.append(np.asarray(load.xyz) * load.mag)
            elif t == 'MOMENT':
                nid = load.node
                if nid in nodes:
                    m_pts.append(nodes[nid].get_position())
                    m_vecs.append(np.asarray(load.xyz) * load.mag)
            elif t == 'PLOAD4':
                eids = load.eids
                if not eids:
                    continue
                idxs = [eid_to_idx.get(int(e)) for e in eids]
                idxs = [i for i in idxs if i is not None]
                if not idxs:
                    continue
                idx_arr = np.asarray(idxs, dtype=np.int64)
                p_pts.append(cell_centers[idx_arr])
                p_vecs.append(cell_normals[idx_arr]
                              * float(load.pressures[0]))
            elif t == 'PLOAD2':
                # v5.2.2 issue 1: PLOAD2 is the uniform-pressure-on-2D-
                # element card. Same visual treatment as PLOAD4 but the
                # pressure value is a scalar attribute, not a list.
                eids = getattr(load, 'eids', None) or []
                if not eids:
                    continue
                idxs = [eid_to_idx.get(int(e)) for e in eids]
                idxs = [i for i in idxs if i is not None]
                if not idxs:
                    continue
                idx_arr = np.asarray(idxs, dtype=np.int64)
                pressure = float(getattr(load, 'pressure', 0.0))
                p_pts.append(cell_centers[idx_arr])
                p_vecs.append(cell_normals[idx_arr] * pressure)

        if f_pts:
            actor = pv.PolyData(np.asarray(f_pts))
            actor['vectors'] = np.asarray(f_vecs)
            actor.set_active_vectors('vectors')
            self.plotter.add_mesh(actor, name=f"force_actors_{sid}",
                                  style='wireframe',
                                  render_lines_as_tubes=False,
                                  show_edges=False, pickable=False,
                                  color='red', opacity=0,
                                  reset_camera=False)
        if m_pts:
            actor = pv.PolyData(np.asarray(m_pts))
            actor['vectors'] = np.asarray(m_vecs)
            actor.set_active_vectors('vectors')
            self.plotter.add_mesh(actor, name=f"moment_actors_{sid}",
                                  style='wireframe',
                                  render_lines_as_tubes=False,
                                  show_edges=False, pickable=False,
                                  color='cyan', opacity=0,
                                  reset_camera=False)
        if p_pts:
            pts_arr = np.concatenate(p_pts, axis=0)
            vecs_arr = np.concatenate(p_vecs, axis=0)
            actor = pv.PolyData(pts_arr)
            actor['vectors'] = vecs_arr
            actor.set_active_vectors('vectors')
            self.plotter.add_mesh(actor, name=f"pressure_actors_{sid}",
                                  style='wireframe',
                                  render_lines_as_tubes=False,
                                  show_edges=False, pickable=False,
                                  color='yellow', opacity=0,
                                  reset_camera=False)

        self._lazy_load_built_sids.add(sid)
        perf_event('lazy_load_actor', 'built', sid=sid,
                   n_force=len(f_pts), n_moment=len(m_pts),
                   n_pressure=int(sum(arr.shape[0] for arr in p_pts)) if p_pts else 0)


    def _create_all_constraint_actors(self):
        """Creates all constraint actors for the current model, one per SID."""
        from node_runner.profiling import perf_stage, perf_event
        actor_names_to_remove = [name for name in self.plotter.actors if name.startswith('constraint_actors_')]
        self.plotter.remove_actor(actor_names_to_remove, render=False)

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.spcs: return

        perf_event('constraints', 'rebuild_all', n_sids=len(model.spcs),
                   removed=len(actor_names_to_remove))
        for sid, spc_list in model.spcs.items():
            with perf_stage('constraints', 'per_sid', sid=int(sid),
                            n_spc_cards=len(spc_list)):
                coords = []
                for spc in spc_list:
                    for nid in spc.nodes:
                        if nid in model.nodes:
                            coords.append(model.nodes[nid].get_position())
                if coords:
                    actor = pv.PolyData(np.array(coords))
                    self.plotter.add_mesh(actor, name=f"constraint_actors_{sid}", style='wireframe', render_lines_as_tubes=False, show_edges=False, pickable=False, color='cyan', opacity=0)

    def _create_rbe_actors(self):
        """v4.0.2 (Stage B): vectorized RBE2/RBE3 spider-line builder.

        Pre-v4.0.2 the per-element loop called
        ``model.nodes[nid].get_position()`` for the center + every
        leg (pyNastran-level Python overhead). On production-deck (46,782 RBE/
        RBAR cards with avg ~10 legs each) this took 144 s.

        New path: vectorized via ``np.searchsorted`` against
        ``self.current_node_ids_sorted`` + array indexing into
        ``self.current_grid.points``. The first pass (~46k iterations)
        only gathers NIDs into Python lists. no per-leg `get_position`
        calls. The second pass is pure numpy.

        Legacy fallback retained as
        ``_create_rbe_actors_legacy``; the new path falls back on any
        mismatch detected via a length assertion.
        """
        from node_runner.profiling import perf_stage, perf_event
        self.plotter.remove_actor('rbe_actors', render=False)

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.rigid_elements: return

        # v4.0.2 (Stage B): integrity counters.
        integrity = {
            'center_missing': 0,
            'partial_legs': 0,
            'empty': 0,
            'sample_eids': [],
        }
        # Per-spoke source data, collected in pass 1.
        center_nids = []
        leg_nids = []
        spoke_eids = []
        spoke_types = []  # 0 = RBE2, 1 = RBE3

        with perf_stage('rbe', 'gather_spokes',
                        n_rigid=len(model.rigid_elements)):
            for eid, elem in model.rigid_elements.items():
                try:
                    t = elem.type
                    if t == 'RBE2':
                        center_nid = elem.gn
                        legs = list(elem.Gmi or [])
                        rbe_type = 0
                    elif t == 'RBE3':
                        center_nid = elem.refgrid
                        legs = []
                        for wt_cg in (elem.wt_cg_groups or []):
                            legs.extend(wt_cg[2] or [])
                        rbe_type = 1
                    else:
                        continue

                    if center_nid is None or center_nid not in model.nodes:
                        integrity['center_missing'] += 1
                        if len(integrity['sample_eids']) < 50:
                            integrity['sample_eids'].append(int(eid))
                        continue

                    n_valid = 0
                    n_missing = 0
                    for leg_nid in legs:
                        if leg_nid is None or leg_nid not in model.nodes:
                            n_missing += 1
                            continue
                        center_nids.append(int(center_nid))
                        leg_nids.append(int(leg_nid))
                        spoke_eids.append(int(eid))
                        spoke_types.append(rbe_type)
                        n_valid += 1
                    if n_valid == 0:
                        integrity['empty'] += 1
                        if len(integrity['sample_eids']) < 50:
                            integrity['sample_eids'].append(int(eid))
                    elif n_missing > 0:
                        integrity['partial_legs'] += 1
                        if len(integrity['sample_eids']) < 50:
                            integrity['sample_eids'].append(int(eid))
                except (AttributeError, KeyError, IndexError) as e:
                    print(f"Warning: Skipping RBE {eid} visualization: {e}", flush=True)
                    continue

        try:
            model._rbe_integrity_warnings = integrity
        except Exception:
            pass

        n_spokes = len(center_nids)
        perf_event('rbe', 'count_spokes',
                   n_rigid=len(model.rigid_elements),
                   n_spokes=n_spokes,
                   center_missing=integrity['center_missing'],
                   empty=integrity['empty'],
                   partial=integrity['partial_legs'])

        if n_spokes == 0:
            return

        # v4.0.2 (Stage B): vectorized NID → grid-index lookup.
        # ``current_node_ids_sorted`` is sorted ascending; searchsorted
        # gives O(log N) per query, vectorized over the whole array.
        nids_sorted = np.asarray(self.current_node_ids_sorted, dtype=np.int64)
        n_nodes = nids_sorted.size
        pts_array = np.asarray(self.current_grid.points)

        with perf_stage('rbe', 'searchsorted_lookup', n_spokes=n_spokes,
                        n_nodes=int(n_nodes)):
            centers_arr = np.asarray(center_nids, dtype=np.int64)
            legs_arr = np.asarray(leg_nids, dtype=np.int64)
            center_idx = np.searchsorted(nids_sorted, centers_arr)
            leg_idx = np.searchsorted(nids_sorted, legs_arr)
            # Mask invalid (out-of-range or non-matching) entries.
            valid = (center_idx < n_nodes) & (leg_idx < n_nodes)
            if valid.any():
                # Clamp before indexing to avoid out-of-range.
                ci_clamped = np.clip(center_idx, 0, n_nodes - 1)
                li_clamped = np.clip(leg_idx, 0, n_nodes - 1)
                valid &= (nids_sorted[ci_clamped] == centers_arr)
                valid &= (nids_sorted[li_clamped] == legs_arr)

        n_valid = int(valid.sum())
        perf_event('rbe', 'valid_after_lookup',
                   n_spokes=n_spokes, n_valid=n_valid,
                   n_dropped=n_spokes - n_valid)

        if n_valid == 0:
            return

        with perf_stage('rbe', 'fill_arrays', n_spokes=n_valid):
            center_pts = pts_array[center_idx[valid]]  # (M, 3)
            leg_pts = pts_array[leg_idx[valid]]        # (M, 3)
            spoke_types_arr = np.asarray(spoke_types, dtype=np.int64)[valid]
            spoke_eids_arr = np.asarray(spoke_eids, dtype=np.int64)[valid]

            # Interleave to [center, leg, center, leg, ...].
            points = np.empty((2 * n_valid, 3), dtype=center_pts.dtype)
            points[0::2] = center_pts
            points[1::2] = leg_pts

            # VTK lines array: [2, p0, p1, 2, p2, p3, …]
            lines = np.empty(3 * n_valid, dtype=np.int64)
            lines[0::3] = 2
            lines[1::3] = np.arange(0, 2 * n_valid, 2)
            lines[2::3] = np.arange(1, 2 * n_valid, 2)

            # Vectorized color assignment: pick from a 2-row palette
            # by spoke type.
            from PySide6.QtGui import QColor as _QColor
            rbe2_rgb = np.array(
                [int(x * 255) for x in
                 _QColor(self.type_color_map.get('RBE2', '#ff3131')).getRgbF()[:3]],
                dtype=np.uint8)
            rbe3_rgb = np.array(
                [int(x * 255) for x in
                 _QColor(self.type_color_map.get('RBE3', '#ffd700')).getRgbF()[:3]],
                dtype=np.uint8)
            color_rgb = np.where(spoke_types_arr[:, None] == 0,
                                 rbe2_rgb, rbe3_rgb).astype(np.uint8)
            point_colors = np.repeat(color_rgb, 2, axis=0)

        with perf_stage('rbe', 'build_polydata',
                        n_points=int(points.shape[0]),
                        n_lines=n_valid):
            poly = pv.PolyData(points, lines=lines)
            poly.point_data['colors'] = point_colors
            poly.cell_data['EID'] = spoke_eids_arr

        # Sanity assertion. if anything is off, we want to know in
        # the log rather than producing a silently-wrong actor.
        if poly.n_points != 2 * n_valid or poly.n_cells != n_valid:
            perf_event('rbe', 'mismatch',
                       expected_points=2 * n_valid,
                       expected_cells=n_valid,
                       actual_points=int(poly.n_points),
                       actual_cells=int(poly.n_cells))
            # Fall back to the legacy path on shape mismatch.
            try:
                self._create_rbe_actors_legacy()
            except Exception as exc:
                perf_event('rbe', 'legacy_fallback_failed',
                           exc=str(exc)[:120])
            return

        self.plotter.add_mesh(poly, scalars='colors', rgb=True,
                              render_lines_as_tubes=True, line_width=3,
                              name='rbe_actors', pickable=True)

    def _create_rbe_actors_legacy(self):
        """Pre-v4.0.2 RBE spider-line builder. Kept as a fallback for
        the vectorized path. Slow on big decks (~144 s on a large production deck).
        """
        self.plotter.remove_actor('rbe_actors', render=False)

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.rigid_elements: return

        points = []
        lines = []
        colors = []
        line_eids = []
        pt_idx = 0

        for eid, elem in model.rigid_elements.items():
            try:
                if elem.type == 'RBE2':
                    center_nid = elem.gn
                    leg_nids = list(elem.Gmi or [])
                    color = self.type_color_map.get('RBE2', '#ff3131')
                elif elem.type == 'RBE3':
                    center_nid = elem.refgrid
                    leg_nids = []
                    for wt_cg in (elem.wt_cg_groups or []):
                        leg_nids.extend(wt_cg[2] or [])
                    color = self.type_color_map.get('RBE3', '#ffd700')
                else:
                    continue

                if center_nid is None or center_nid not in model.nodes:
                    continue

                center_pos = model.nodes[center_nid].get_position()

                for leg_nid in leg_nids:
                    if leg_nid is None or leg_nid not in model.nodes:
                        continue
                    leg_pos = model.nodes[leg_nid].get_position()
                    points.append(center_pos)
                    points.append(leg_pos)
                    lines.append([2, pt_idx, pt_idx + 1])
                    colors.append(color)
                    line_eids.append(int(eid))
                    pt_idx += 2
            except (AttributeError, KeyError, IndexError):
                continue

        if not points:
            return

        flat_lines = np.asarray(
            [val for line in lines for val in line], dtype=np.int64)
        poly = pv.PolyData(np.array(points), lines=flat_lines)

        from PySide6.QtGui import QColor
        point_colors = []
        for c in colors:
            rgb = [int(x * 255) for x in QColor(c).getRgbF()[:3]]
            point_colors.append(rgb)
            point_colors.append(rgb)
        poly.point_data['colors'] = np.array(point_colors, dtype=np.uint8)
        poly.cell_data['EID'] = np.asarray(line_eids, dtype=np.int64)

        self.plotter.add_mesh(poly, scalars='colors', rgb=True,
                              render_lines_as_tubes=True, line_width=3,
                              name='rbe_actors', pickable=True)

    def _create_mass_actors(self):
        """Creates cube glyph visualization for CONM2 mass elements.

        v3.3.0: each glyph cube's 6 faces are tagged with the mass
        element's eid in cell_data['EID'] so the cell picker can map
        any clicked face back to the mass element. Actor is now pickable.
        """
        self.plotter.remove_actor('mass_actors', render=False)

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.masses: return

        coords = []
        eids = []
        for eid, mass_elem in model.masses.items():
            try:
                nid = mass_elem.nid
                if nid in model.nodes:
                    coords.append(model.nodes[nid].get_position())
                    eids.append(int(eid))
            except (AttributeError, KeyError):
                continue

        if not coords:
            return

        pts = pv.PolyData(np.array(coords))
        cube = pv.Cube(x_length=1.0, y_length=1.0, z_length=1.0)
        # v3.5.0 item 3: mass glyph scale is now a user preference
        # (Preferences > Entity Sizes > Mass glyph scale).
        # v4.0.7: default lowered to 0.0075 (half of pre-v4.0.7 0.015).
        scale_factor = getattr(self, 'mass_glyph_scale', 0.0075)
        scale = (self.current_grid.length * scale_factor
                 if self.current_grid.length > 0 else 1.0)
        glyphs = pts.glyph(scale=False, factor=scale, geom=cube)
        # v3.3.0: each input point produces cube.n_cells cells in the
        # glyph output. Repeat each mass eid across its 6 cube faces so
        # cell_data['EID'][cell_id] yields the correct mass element id.
        cells_per_glyph = max(1, cube.n_cells)
        n_out = glyphs.n_cells
        expected = len(eids) * cells_per_glyph
        if n_out == expected:
            glyphs.cell_data['EID'] = np.repeat(
                np.asarray(eids, dtype=np.int64), cells_per_glyph)
        else:
            # Defensive: if the glyph filter ever changes its output
            # layout, fall back to a best-effort even split.
            try:
                per = max(1, n_out // max(1, len(eids)))
                glyphs.cell_data['EID'] = np.repeat(
                    np.asarray(eids, dtype=np.int64), per)[:n_out]
            except Exception:
                pass
        self.plotter.add_mesh(glyphs, color=self.type_color_map.get('Masses', '#00cc66'),
                              name='mass_actors', pickable=True)

    def _create_plotel_actors(self):
        """Creates line visualization for PLOTEL elements."""
        self.plotter.remove_actor('plotel_actors', render=False)
        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.plotels: return
        points = []; lines = []; pt_idx = 0
        for eid, elem in model.plotels.items():
            try:
                n1, n2 = elem.nodes
                if n1 in model.nodes and n2 in model.nodes:
                    points.append(model.nodes[n1].get_position())
                    points.append(model.nodes[n2].get_position())
                    lines.append([2, pt_idx, pt_idx + 1])
                    pt_idx += 2
            except (KeyError, AttributeError): continue
        if not points: return
        poly = pv.PolyData(np.array(points))
        poly.lines = np.array([v for l in lines for v in l])
        color = self.type_color_map.get('Plotel', '#ff69b4')
        self.plotter.add_mesh(poly, color=color, render_lines_as_tubes=True,
                              line_width=2, name='plotel_actors', pickable=False)

    def _create_element_normals_actor(self):
        """
        Calculates and creates a renderable actor for shell element normals by
        manually performing a cross-product on each element's nodes. This
        ensures the visualization is a 1:1 match with the model's data.
        """
        self.plotter.remove_actor('element_normals_actor', render=False)

        if not self.current_grid or 'is_shell' not in self.current_grid.cell_data:
            return

        shell_indices = np.where(self.current_grid.cell_data['is_shell'] == 1)[0]
        if shell_indices.size == 0:
            return
        
        centers = []
        vectors = []

        # --- FIX: Manually calculate the normal for each shell element ---
        # This bypasses PyVista's automatic normal "fixing" and guarantees we
        # visualize the true normal based on the raw node order from the model.
        for i in shell_indices:
            cell = self.current_grid.get_cell(i)
            points = cell.points
            
            # Calculate normal using the cross product, respecting Nastran's right-hand rule
            if cell.type == vtk.VTK_QUAD or cell.type == vtk.VTK_PIXEL: # For CQUAD4
                # Use the two diagonals for a robust calculation
                vec1 = points[2] - points[0]
                vec2 = points[3] - points[1]
                normal = np.cross(vec1, vec2)
            elif cell.type == vtk.VTK_TRIANGLE: # For CTRIA3
                vec1 = points[1] - points[0]
                vec2 = points[2] - points[0]
                normal = np.cross(vec1, vec2)
            else:
                continue # Should not be possible due to our filter

            # Normalize the vector to get a unit direction
            norm_len = np.linalg.norm(normal)
            if norm_len > 1e-9:
                normal = normal / norm_len

            centers.append(cell.center)
            vectors.append(normal)

        if not centers:
            return

        try:
            scale_percent = float(self.normal_arrow_scale_input.text())
        except ValueError:
            scale_percent = 2.5
        
        scale = self.plotter.camera.parallel_scale * (scale_percent / 50.0)

        points_pd = pv.PolyData(np.array(centers))
        points_pd['vectors'] = np.array(vectors)
        arrow_geom = pv.Arrow(tip_length=0.2, shaft_radius=0.02)
        glyphs = points_pd.glyph(orient='vectors', scale=False, factor=scale, geom=arrow_geom)
        self.plotter.add_mesh(glyphs, name='element_normals_actor', color='#ffbe0b', reset_camera=False)

    def _create_free_edges_actor(self):
        """Creates a line overlay highlighting free (unshared) shell element edges."""
        self.plotter.remove_actor('free_edges_actor', render=False)

        if not self.current_generator or not self.current_generator.model.elements:
            return

        free_edges = self.current_generator.find_free_edges()
        if not free_edges:
            self._update_status("No free edges found.")
            return

        # Build line PolyData from edge node pairs
        points = []
        lines = []
        node_map = {}  # nid -> index in points list
        model_nodes = self.current_generator.model.nodes

        for nid1, nid2 in free_edges:
            if nid1 not in model_nodes or nid2 not in model_nodes:
                continue  # Skip edges with missing nodes
            for nid in (nid1, nid2):
                if nid not in node_map:
                    node_map[nid] = len(points)
                    points.append(model_nodes[nid].get_position())
            lines.append([2, node_map[nid1], node_map[nid2]])

        if not points:
            self._update_status("No free edges with valid nodes found.")
            return

        poly = pv.PolyData(np.array(points), lines=np.hstack(lines))
        self.plotter.add_mesh(poly, name='free_edges_actor', color='#ff3131',
                              line_width=4, render_lines_as_tubes=True, reset_camera=False)
        self._update_status(f"Showing {len(free_edges)} free edge(s).")

# START: Final replacement for _create_bush_orientation_glyphs in main.py
    def _create_bush_orientation_glyphs(self, visibility_mask):
        """Creates and renders orientation triads for visible CBUSH elements."""
        self.plotter.remove_actor('bush_orientation_glyphs', render=False)
        if not self.show_bush_orient_check.isChecked() or not self.current_grid or self.current_grid.n_cells == 0:
            return

        if 'type' not in self.current_grid.cell_data: return
        all_bush_indices_mask = self.current_grid.cell_data['type'] == 'CBUSH'
        
        visible_bush_mask = all_bush_indices_mask & visibility_mask
        visible_bush_indices = np.where(visible_bush_mask)[0]

        if visible_bush_indices.size == 0:
            return

        visible_eids = self.current_grid.cell_data['EID'][visible_bush_indices]
        centers = self.current_grid.cell_centers().points[visible_bush_indices]
        
        glyph_points = []
        glyph_vectors = []
        glyph_colors = []
        
        # --- MODIFIED: Define colors as numerical RGB tuples ---
        colors = {
            "x": [255, 0, 0],   # Red
            "y": [0, 255, 0],   # Green
            "z": [0, 0, 255]    # Blue
        }
        
        scale = self.plotter.camera.parallel_scale * 0.05

        for i, eid in enumerate(visible_eids):
            matrix = self.current_generator.get_cbush_orientation_matrix(eid)
            center_point = centers[i]

            glyph_points.extend([center_point, center_point, center_point])
            glyph_vectors.extend([matrix[:,0], matrix[:,1], matrix[:,2]])
            glyph_colors.extend([colors['x'], colors['y'], colors['z']])

        if not glyph_points: return

        polydata = pv.PolyData(np.array(glyph_points))
        polydata['vectors'] = np.array(glyph_vectors)
        polydata['colors'] = np.array(glyph_colors, dtype=np.uint8)

        arrows = polydata.glyph(orient='vectors', scale=False, factor=scale, geom=pv.Arrow(tip_length=0.2, shaft_radius=0.02))
        self.plotter.add_mesh(arrows, scalars='colors', rgb=True, name='bush_orientation_glyphs', reset_camera=False)
# END: Final replacement for _create_bush_orientation_glyphs in main.py
    


# START: New glyph creation method in main.py
    def _create_beam_orientation_glyphs(self, visibility_mask):
        """Creates and renders orientation triads for visible CBEAM/CBAR elements."""
        self.plotter.remove_actor('beam_orientation_glyphs', render=False)
        if not self.show_beam_orient_check.isChecked() or not self.current_grid or self.current_grid.n_cells == 0:
            return

        if 'type' not in self.current_grid.cell_data: return
        beam_types = ['CBEAM', 'CBAR']
        all_beam_indices_mask = np.isin(self.current_grid.cell_data['type'], beam_types)
        
        visible_beam_mask = all_beam_indices_mask & visibility_mask
        visible_beam_indices = np.where(visible_beam_mask)[0]

        if visible_beam_indices.size == 0:
            return

        visible_eids = self.current_grid.cell_data['EID'][visible_beam_indices]
        centers = self.current_grid.cell_centers().points[visible_beam_indices]
        
        glyph_points, glyph_vectors, glyph_colors = [], [], []
        colors = {"x": [255, 0, 0], "y": [0, 255, 0], "z": [0, 0, 255]}
        scale = self.plotter.camera.parallel_scale * 0.05

        for i, eid in enumerate(visible_eids):
            matrix = self.current_generator.get_beam_orientation_matrix(eid)
            center_point = centers[i]

            glyph_points.extend([center_point, center_point, center_point])
            glyph_vectors.extend([matrix[:,0], matrix[:,1], matrix[:,2]])
            glyph_colors.extend([colors['x'], colors['y'], colors['z']])

        if not glyph_points: return

        polydata = pv.PolyData(np.array(glyph_points))
        polydata['vectors'] = np.array(glyph_vectors)
        polydata['colors'] = np.array(glyph_colors, dtype=np.uint8)

        arrows = polydata.glyph(orient='vectors', scale=False, factor=scale, geom=pv.Arrow(tip_length=0.2, shaft_radius=0.02))
        self.plotter.add_mesh(arrows, scalars='colors', rgb=True, name='beam_orientation_glyphs', reset_camera=False)

    def _schedule_glyph_refresh(self):
        """Start/restart the debounce timer for orientation glyph refresh."""
        if hasattr(self, '_glyph_refresh_timer'):
            self._glyph_refresh_timer.start()

    def _refresh_orientation_glyphs(self):
        """Recreate visible orientation glyphs using cached visibility mask."""
        if self._last_visibility_mask is None or not self.current_grid:
            return
        mask = self._last_visibility_mask
        if self.show_bush_orient_check.isChecked():
            self._create_bush_orientation_glyphs(mask)
        if self.show_beam_orient_check.isChecked():
            self._create_beam_orientation_glyphs(mask)
        if self.show_normals_check.isChecked():
            self._create_element_normals_actor()
        self.plotter.render()



    def _create_coord_system(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        
        next_cid = max(list(self.current_generator.model.coords.keys()) + [0]) + 1
        
        # Pass the model so the dialog can populate the reference CID list
        dialog = CreateCoordDialog(next_cid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_coord_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _on_coord_creation_accept(self):
        dialog = self.active_creation_dialog
        if not dialog: return

        if params := dialog.get_parameters():
            cmd = AddCoordCommand(params)
            self.command_manager.execute(cmd, self.current_generator.model)
            self._update_status(f"Created/Updated Coordinate System CID {params['cid']}.")
            self._update_viewer(self.current_generator, reset_camera=False)

        self.active_creation_dialog = None


    def _create_coord_system_actors(self, cid):
        """Creates static csys triad arrows + optional label for a coordinate system.

        Reads display-tab settings (axis size, label toggle).
        Safely skips if the coord object cannot provide origin/axes.
        """
        if not self.current_generator or cid not in self.current_generator.model.coords:
            return False

        try:
            coord = self.current_generator.model.coords[cid]
            origin = np.asarray(coord.origin, dtype=float)
            x_axis = np.asarray(coord.i, dtype=float)
            y_axis = np.asarray(coord.j, dtype=float)
            z_axis = np.asarray(coord.k, dtype=float)
        except Exception:
            return False  # coord not fully initialised (e.g. cross-ref missing)

        # Scale from Display-tab slider (‰ of model diagonal)
        model_size = (self.current_grid.length
                      if self.current_grid and self.current_grid.n_points > 1
                      else 10.0)
        if np.isclose(model_size, 0.0):
            model_size = 10.0
        scale = model_size * self.csys_scale_slider.value() / 1000.0

        arrow_kw = dict(shaft_radius=0.015, tip_radius=0.04, tip_length=0.25)
        self.plotter.add_mesh(pv.Arrow(start=origin, direction=x_axis, scale=scale, **arrow_kw),
                              color='red', name=f"coord_{cid}_x", reset_camera=False)
        self.plotter.add_mesh(pv.Arrow(start=origin, direction=y_axis, scale=scale, **arrow_kw),
                              color='green', name=f"coord_{cid}_y", reset_camera=False)
        self.plotter.add_mesh(pv.Arrow(start=origin, direction=z_axis, scale=scale, **arrow_kw),
                              color='blue', name=f"coord_{cid}_z", reset_camera=False)

        # Label - honour the "Show Labels" checkbox
        if self.show_csys_labels_check.isChecked():
            label_map = {0: "CID 0: Global Rect", 1: "CID 1: Global Cyl",
                         2: "CID 2: Global Sph"}
            label_text = label_map.get(cid, f"CID {cid}")
            label_color = 'white' if self.is_dark_theme else 'black'
            self.plotter.add_point_labels(
                origin, [label_text], name=f"coord_label_{cid}",
                font_size=10, text_color=label_color, shape=None,
                always_visible=True, reset_camera=False)
        return True

    def _remove_coord_actors(self, cid):
        """Remove all actors (arrows + label) for a coordinate system."""
        for suffix in ('_x', '_y', '_z'):
            self.plotter.remove_actor(f"coord_{cid}{suffix}", render=False)
        self.plotter.remove_actor(f"coord_label_{cid}", render=False)

    def _refresh_coord_actors(self, _=None, render=True):
        """Rebuild all coord-system actors honouring Display-tab settings + tree state."""
        if not self.current_generator:
            return
        from node_runner.profiling import perf_stage, perf_event
        model = self.current_generator.model

        master_on = self.show_csys_check.isChecked()
        show_global = self.show_global_csys_check.isChecked()

        # v4.0.9: sub-stage timers. The whole refresh costs ~0.22 s per
        # visibility cycle on a large production deck despite only 3 coord systems; sub-stages
        # tell us which step actually dominates.
        with perf_stage('coord_actors', 'remove_all',
                        n_coords=len(model.coords)):
            for cid in list(model.coords.keys()):
                self._remove_coord_actors(cid)

        if not master_on:
            if render:
                self.plotter.render()
            return

        with perf_stage('coord_actors', 'resolve_tree_visibility',
                        n_coords=len(model.coords)):
            tree_visible = {}
            coord_group_items = self._find_tree_items(('group', 'coords'))
            if coord_group_items:
                parent = coord_group_items[0]
                for i in range(parent.childCount()):
                    child = parent.child(i)
                    child_data = child.data(0, QtCore.Qt.UserRole)
                    if isinstance(child_data, tuple) and len(child_data) >= 2:
                        tree_visible[child_data[1]] = (
                            child.checkState(0) == QtCore.Qt.Checked)

        built = 0
        with perf_stage('coord_actors', 'build_visible',
                        n_coords=len(model.coords)):
            for cid in sorted(model.coords.keys()):
                if cid <= 2 and not show_global:
                    continue
                if not tree_visible.get(cid, cid > 2):
                    continue
                self._create_coord_system_actors(cid)
                built += 1
        perf_event('coord_actors', 'build_summary', n_built=built)

        if render:
            self.plotter.render()



    def visible_entity_ids(self, entity_type):
        """v3.3.0: return the set of currently-visible EIDs / NIDs based
        on tree-checkbox state.

        Consumed by EntitySelectionBar._select_visible. The filter
        mirrors the one used by _update_plot_visibility:
          - elem_by_type_group + elem_by_shape_group checkboxes hide
            elements whose type maps to an unchecked group.
          - property / material checkboxes hide elements by PID.
          - Masses group checkbox hides CONM2.
          - Rigid By-Type checkbox hides RBE2/RBE3.
        For 'Node', returns every grid id (no per-node tree checkbox).
        """
        if not self.current_generator:
            return set()
        model = self.current_generator.model
        if entity_type == 'Node':
            return set(model.nodes.keys())
        if entity_type != 'Element':
            return set()

        # Read tree state by group key.
        def _checked(parent_data, child_label):
            for parent in self._find_tree_items(parent_data):
                for i in range(parent.childCount()):
                    child = parent.child(i)
                    cd = child.data(0, QtCore.Qt.UserRole)
                    if isinstance(cd, tuple) and cd[-1] == child_label:
                        return child.checkState(0) == QtCore.Qt.Checked
            return True  # default visible if no tree item exists

        type_to_nastran = {
            'Beams': ('CBEAM',), 'Bars': ('CBAR',), 'Rods': ('CROD',),
            'Bushes': ('CBUSH',),
            'Plates': ('CQUAD4', 'CMEMBRAN', 'CTRIA3'),
            'Rigid': ('RBE2', 'RBE3'),
            'Solids': ('CHEXA', 'CHEXA8', 'CHEXA20',
                       'CTETRA', 'CTETRA4', 'CTETRA10',
                       'CPENTA', 'CPENTA6', 'CPENTA15'),
            'Shear': ('CSHEAR',), 'Gap': ('CGAP',),
        }
        shape_to_nastran = {
            'Line': ('CBEAM', 'CBAR', 'CROD', 'CBUSH', 'CGAP'),
            'Quad': ('CQUAD4', 'CMEMBRAN', 'CSHEAR'),
            'Tria': ('CTRIA3',),
            'Rigid': ('RBE2', 'RBE3'),
            'Hex': ('CHEXA', 'CHEXA8', 'CHEXA20'),
            'Tet': ('CTETRA', 'CTETRA4', 'CTETRA10'),
            'Wedge': ('CPENTA', 'CPENTA6', 'CPENTA15'),
        }
        type_check = {}
        for group, types in type_to_nastran.items():
            ok = True
            # Iterate to find the by-type item with matching group name.
            for parent in self._find_tree_items(('group', 'elem_by_type')):
                for i in range(parent.childCount()):
                    child = parent.child(i)
                    cd = child.data(0, QtCore.Qt.UserRole)
                    if (isinstance(cd, tuple)
                            and cd[0] == 'elem_by_type_group'
                            and cd[1] == group):
                        ok = child.checkState(0) == QtCore.Qt.Checked
            for t in types:
                type_check[t] = ok and type_check.get(t, True)

        for group, types in shape_to_nastran.items():
            ok = True
            for parent in self._find_tree_items(('group', 'elem_by_shape')):
                for i in range(parent.childCount()):
                    child = parent.child(i)
                    cd = child.data(0, QtCore.Qt.UserRole)
                    if (isinstance(cd, tuple)
                            and cd[0] == 'elem_shape_group'
                            and cd[1] == group):
                        ok = child.checkState(0) == QtCore.Qt.Checked
            for t in types:
                type_check[t] = type_check.get(t, True) and ok

        # Property visibility
        hidden_props = set()
        for parent in self._find_tree_items(('group', 'properties')):
            for i in range(parent.childCount()):
                child = parent.child(i)
                if child.checkState(0) != QtCore.Qt.Checked:
                    cd = child.data(0, QtCore.Qt.UserRole)
                    if isinstance(cd, tuple) and cd[0] == 'property':
                        hidden_props.add(cd[1])

        # Mass group
        mass_visible = True
        for parent in self._find_tree_items(('group', 'masses')):
            mass_visible = (parent.checkState(0) == QtCore.Qt.Checked)

        visible = set()
        for eid, elem in (model.elements or {}).items():
            etype = getattr(elem, 'type', None)
            if not type_check.get(etype, True):
                continue
            pid = getattr(elem, 'pid', None)
            if pid in hidden_props:
                continue
            visible.add(int(eid))
        for eid, elem in (model.rigid_elements or {}).items():
            etype = getattr(elem, 'type', None)
            if not type_check.get(etype, True):
                continue
            visible.add(int(eid))
        if mass_visible:
            for eid in (model.masses or {}).keys():
                visible.add(int(eid))
        return visible

    def _update_plot_visibility(self):
        """Fast update of actor visibilities based on tree state. Does not rebuild actors."""
        from node_runner.profiling import perf_stage, perf_event
        _t_total = time.perf_counter()
        # v4.0.2 (Stage A): record the _find_tree_items counter at entry
        # so we know how many lookups this single visibility cycle did.
        self._find_tree_calls_at_entry = getattr(self, '_find_tree_calls_total', 0)
        if not self.current_generator or not self.current_grid:
            with perf_stage('visibility', 'rebuild_plot_empty'):
                self._rebuild_plot()
            return
        # v4.0.10: defensive refresh of the maintained visible-SID sets
        # from current tree state. Cheap (~5-10 ms) and ensures the
        # bounded load_sid_loop iteration below sees the truth even
        # if a fast-path event slipped through unobserved.
        try:
            self._refresh_visible_sid_sets()
        except Exception:
            pass

        n_cells = int(self.current_grid.n_cells)
        visibility_mask = np.ones(n_cells, dtype=bool)
        model = self.current_generator.model
        # v4.0.14: 'Bushes' was missing. Unchecking the Bushes
        # category tree item ran the full visibility cycle but never
        # masked CBUSH cells, so the visual didn't change.
        type_to_nastran_map = { 'Beams': ['CBEAM'], 'Bars': ['CBAR'], 'Rods': ['CROD'], 'Bushes': ['CBUSH'], 'Plates': ['CQUAD4', 'CMEMBRAN', 'CTRIA3'], 'Rigid': ['RBE2', 'RBE3'], 'Solids': ['CHEXA', 'CHEXA8', 'CHEXA20', 'CTETRA', 'CTETRA4', 'CTETRA10', 'CPENTA', 'CPENTA6', 'CPENTA15'], 'Shear': ['CSHEAR'], 'Gap': ['CGAP'] }
        shape_to_nastran_map = { 'Line': ['CBEAM', 'CBAR', 'CROD', 'CBUSH', 'CGAP'], 'Quad': ['CQUAD4', 'CMEMBRAN', 'CSHEAR'], 'Tria': ['CTRIA3'], 'Rigid': ['RBE2', 'RBE3'], 'Hex': ['CHEXA', 'CHEXA8', 'CHEXA20'], 'Tet': ['CTETRA', 'CTETRA4', 'CTETRA10'], 'Wedge': ['CPENTA', 'CPENTA6', 'CPENTA15'] }
        with perf_stage('visibility', 'tree_walk_type_shape', n_cells=n_cells):
            # v4.0.7: use the tree index instead of walking the
            # entire model tree (~23k items on a large production deck, ~1.3 s wall
            # per visibility cycle). With the v4.0.3 tree index in
            # place we can look up each category checkbox by its
            # known UserRole key in O(1). One lookup per category
            # instead of one walk per cycle.
            for cat_name, nastran_types in type_to_nastran_map.items():
                items = self._find_tree_items(('elem_by_type_group', cat_name))
                if items and items[0].checkState(0) != QtCore.Qt.Checked:
                    if nastran_types:
                        visibility_mask[np.isin(
                            self.current_grid.cell_data['type'],
                            nastran_types)] = False
            for shape_name, nastran_types in shape_to_nastran_map.items():
                items = self._find_tree_items(('elem_shape_group', shape_name))
                if items and items[0].checkState(0) != QtCore.Qt.Checked:
                    if nastran_types:
                        visibility_mask[np.isin(
                            self.current_grid.cell_data['type'],
                            nastran_types)] = False
        with perf_stage('visibility', 'tree_walk_pid_mid', n_cells=n_cells):
            props_to_hide = {prop_item.data(0, QtCore.Qt.UserRole)[1] for item in self._find_tree_items(('group', 'properties')) for i in range(item.childCount()) if (prop_item := item.child(i)).checkState(0) != QtCore.Qt.Checked}
            if props_to_hide: visibility_mask[np.isin(self.current_grid.cell_data['PID'], list(props_to_hide))] = False
            mids_to_hide = {mat_item.data(0, QtCore.Qt.UserRole)[1] for item in self._find_tree_items(('group', 'materials')) for i in range(item.childCount()) if (mat_item := item.child(i)).checkState(0) != QtCore.Qt.Checked}
        if mids_to_hide:
            props_using_hidden_mats = set()
            for pid, prop in model.properties.items():
                for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                    val = getattr(prop, attr, None)
                    if val is not None and val in mids_to_hide:
                        props_using_hidden_mats.add(pid)
                        break
            if props_using_hidden_mats: visibility_mask[np.isin(self.current_grid.cell_data['PID'], list(props_using_hidden_mats))] = False
        if (rbe_actor := self.plotter.actors.get('rbe_actors')):
            rbe_actor.SetVisibility(not any(item.checkState(0) != QtCore.Qt.Checked for item in self._find_tree_items(('elem_by_type_group', 'Rigid'))))
        if (mass_actor := self.plotter.actors.get('mass_actors')):
            mass_items = self._find_tree_items(('group', 'masses'))
            mass_actor.SetVisibility(not mass_items or mass_items[0].checkState(0) == QtCore.Qt.Checked)
        if (plotel_actor := self.plotter.actors.get('plotel_actors')):
            plotel_items = self._find_tree_items(('group', 'plotels'))
            plotel_actor.SetVisibility(not plotel_items or plotel_items[0].checkState(0) == QtCore.Qt.Checked)

        # --- Phase 8 / v3.4.0: Group visibility ---
        # v3.4.0 (item 7): a group's effective element set is the union
        # of its explicit elements + any element whose PID is in
        # group['properties'] + any element whose material (via PID) is
        # in group['materials']. That way a group can be defined purely
        # by "every shell using PID 100" or "every element using MAT 2"
        # and isolation works as expected.
        def _effective_eids(grp):
            eids = set(grp.get("elements", []) or [])
            grp_pids = set(grp.get("properties", []) or [])
            grp_mids = set(grp.get("materials", []) or [])
            if grp_pids or grp_mids:
                # Build pid -> eids and mid -> pids maps once.
                if grp_mids:
                    pid_for_mid: set[int] = set()
                    for pid, prop in (model.properties or {}).items():
                        for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                            if getattr(prop, attr, None) in grp_mids:
                                pid_for_mid.add(pid)
                                break
                    grp_pids |= pid_for_mid
                if grp_pids:
                    all_e = {**model.elements, **model.rigid_elements}
                    for eid, elem in all_e.items():
                        pid = getattr(elem, 'pid', None)
                        if pid and pid in grp_pids:
                            eids.add(int(eid))
            return eids

        if self._isolate_mode:
            group_data = self.groups.get(self._isolate_mode, {})
            group_eids = _effective_eids(group_data)
            if group_eids:
                visibility_mask &= np.isin(self.current_grid.cell_data['EID'], list(group_eids))
            else:
                visibility_mask[:] = False
        elif self._hidden_groups:
            hidden_eids: set[int] = set()
            for grp_name in self._hidden_groups:
                grp = self.groups.get(grp_name, {})
                hidden_eids |= _effective_eids(grp)
            if hidden_eids:
                visibility_mask[np.isin(self.current_grid.cell_data['EID'], list(hidden_eids))] = False
                # Phase 4.2: ghost-mode overlay - render hidden groups as a
                # translucent wireframe shadow instead of fully removing them
                # from view, so the user keeps spatial context.
                self.plotter.remove_actor('ghost_actor', render=False)
                if self._ghost_mode_enabled:
                    ghost_idx = np.where(np.isin(
                        self.current_grid.cell_data['EID'], list(hidden_eids),
                    ))[0]
                    if ghost_idx.size > 0:
                        ghost_grid = self.current_grid.extract_cells(ghost_idx)
                        self.plotter.add_mesh(
                            ghost_grid,
                            name='ghost_actor',
                            style='wireframe',
                            color='#6c7086',  # Catppuccin overlay-1: muted gray
                            opacity=0.3,
                            line_width=1,
                            pickable=False,
                            reset_camera=False,
                        )
            else:
                self.plotter.remove_actor('ghost_actor', render=False)
        else:
            self.plotter.remove_actor('ghost_actor', render=False)

        self._last_visibility_mask = visibility_mask.copy()
        with perf_stage('visibility', 'rebuild_plot',
                        n_visible=int(visibility_mask.sum()),
                        n_total=n_cells):
            self._rebuild_plot(visibility_mask=visibility_mask)

        try:
            arrow_size_percent = float(self.arrow_scale_input.text()); use_relative_scaling = self.relative_scaling_check.isChecked()
        except (ValueError, AttributeError):
            arrow_size_percent = 10.0; use_relative_scaling = True
        base_scale = self.current_grid.length * (arrow_size_percent / 100.0)

        max_force_mag = self.load_scaling_info.get('force_max_mag', 0.0)
        max_moment_mag = self.load_scaling_info.get('moment_max_mag', 0.0)
        max_pressure_mag = self.load_scaling_info.get('pressure_max_mag', 0.0)

        # v4.0.2 (Stage A): wrap the 1962-SID load-actor visibility
        # loop. The v4.0.1 log showed 265 s unaccounted for between
        # _rebuild_plot and final_render; this loop is the prime
        # suspect because each iteration calls _find_tree_items which
        # walks the entire tree.
        _t_load_loop = time.perf_counter()
        _n_load_visible = 0
        _find_tree_calls_before = getattr(self, '_find_tree_calls_total', 0)
        # v4.0.9: short-circuit when no load SID is visible AND no stale
        # load glyph exists in plotter.actors. Loads default OFF since
        # v4.0.5 so this is the common case on large production decks (1,962 SIDs
        # all unchecked → no work needed). Saves ~1.0-2.0 s per
        # visibility cycle. Catcher stays in forever so we can spot
        # regressions in the dirty/visible accounting.
        # v4.0.10: short-circuit check now reads the maintained
        # _visible_load_sids set (O(1)) instead of walking 1,962 SIDs
        # to check checkboxes. The set is kept in sync by
        # _fast_path_sid_toggle and the defensive
        # _refresh_visible_sid_sets above.
        _short_circuit_loads = False
        if model.loads:
            _any_stale_glyph = any(
                isinstance(name, str) and (
                    name.startswith('force_glyphs_')
                    or name.startswith('moment_glyphs_')
                    or name.startswith('pressure_glyphs_'))
                for name in self.plotter.actors.keys())
            _short_circuit_loads = (not self._visible_load_sids) and (not _any_stale_glyph)
        if _short_circuit_loads:
            perf_event('load_sid_loop', 'short_circuit',
                       n_sids=len(model.loads),
                       taken=True,
                       wall_s=f"{time.perf_counter()-_t_load_loop:.3f}")
        elif model.loads:
            perf_event('load_sid_loop', 'short_circuit',
                       n_sids=len(model.loads),
                       n_visible=len(self._visible_load_sids),
                       taken=False)
            # v4.0.10: stale-glyph sweep. With bounded iteration below,
            # the per-SID loop only touches SIDs currently in
            # _visible_load_sids. Glyphs for previously-visible-but-
            # now-hidden SIDs would leak. Sweep removes those once
            # per cycle. O(N_actors), bounded by visible+stale count
            # (small in practice).
            _stale_removed = 0
            for actor_name in list(self.plotter.actors.keys()):
                if not isinstance(actor_name, str):
                    continue
                for prefix in ('force_glyphs_', 'moment_glyphs_', 'pressure_glyphs_'):
                    if actor_name.startswith(prefix):
                        try:
                            _stale_sid = int(actor_name[len(prefix):])
                        except ValueError:
                            break
                        if _stale_sid not in self._visible_load_sids:
                            self.plotter.remove_actor(actor_name, render=False)
                            _stale_removed += 1
                        break
            if _stale_removed:
                perf_event('load_sid_loop', 'stale_glyphs_swept',
                           n_removed=_stale_removed)
        # v4.0.10: iterate only the visible SIDs instead of all 1,962.
        # Saves ~3.7 s per click on a large production deck when 1 SID visible (was the
        # case in v4.0.9 profile log line 1647).
        if not _short_circuit_loads and model.loads:
            perf_event('load_sid_loop', 'iter_bounded',
                       n_iter=len(self._visible_load_sids),
                       n_total=len(model.loads))
        for sid in (self._visible_load_sids if not _short_circuit_loads else ()):
            is_visible = False
            if (sid_item_list := self._find_tree_items(('load_set', sid))) and sid_item_list[0].checkState(0) == QtCore.Qt.Checked:
                is_visible = True
                _n_load_visible += 1

            # v3.4.2: lazily build per-SID actors the first time the
            # user shows this SID. Decks with thousands of SIDs (e.g.
            # a production deck we exercised had 1962) used to spend
            # ~50ms per SID at import time creating opacity=0 actors;
            # deferring it caps import time and adds only the cost
            # of the few SIDs the user actually shows.
            if is_visible:
                self._ensure_load_actor_for_sid(sid)

            for type, max_mag in [('force', max_force_mag), ('moment', max_moment_mag), ('pressure', max_pressure_mag)]:
                actor_name = f"{type}_actors_{sid}"; glyph_name = f"{type}_glyphs_{sid}"
                # v4.0.6: skip the no-op remove_actor when the glyph
                # doesn't exist. On production-deck (1,962 SIDs × 3 types = 5,886
                # remove_actor calls per visibility cycle), most are
                # no-ops. but each one still walks PyVista's internal
                # actor list. The membership check is O(1) on the
                # plotter.actors dict.
                if glyph_name in self.plotter.actors:
                    self.plotter.remove_actor(glyph_name, render=False)
                if (actor := self.plotter.actors.get(actor_name)) and is_visible:
                    dataset = actor.mapper.dataset
                    if not dataset.n_points: continue
                    if type == 'force':
                        vectors = dataset['vectors']; mags = np.linalg.norm(vectors, axis=1)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                        if use_relative_scaling and max_mag > 1e-9:
                            scaled_vectors = unit_vectors * (mags / max_mag * base_scale)[:, np.newaxis]
                        else:
                            scaled_vectors = unit_vectors * base_scale
                        self.plotter.add_arrows(dataset.points, scaled_vectors, color='red', name=glyph_name, reset_camera=False)
                    elif type == 'moment':
                        # v4.0.5 (Stage T): pv.DoubleArrow was removed
                        # in PyVista 0.46+. Build a "two-tipped" arrow
                        # by merging two pv.Arrow polydatas pointing in
                        # opposite directions. Fallback to plain
                        # pv.Arrow on any failure so a future PyVista
                        # API change won't crash the visibility loop.
                        # (NOTE: ``type`` is shadowed by the outer
                        # for-loop variable in this block, so use
                        # ``_exc.__class__.__name__`` instead of
                        # ``type(_exc).__name__``.)
                        from node_runner.profiling import perf_event
                        # v4.0.12 (Bug 4 fix): two-tier fallback.
                        # composite first, plain pv.Arrow on
                        # composite-side failure. Same pattern as the
                        # _build_load_glyphs_for_sid path. Keep this
                        # mirror in sync.
                        def _try_build_moment_glyph_inline(geom_pd):
                            if use_relative_scaling and max_mag > 1e-9:
                                _mags = np.linalg.norm(dataset['vectors'], axis=1); dataset['relative_mag'] = _mags / max_mag
                                _g = dataset.glyph(orient='vectors', scale='relative_mag', factor=base_scale, geom=geom_pd)
                            else:
                                _g = dataset.glyph(orient='vectors', scale=False, factor=base_scale, geom=geom_pd)
                            self.plotter.add_mesh(_g, color='cyan', name=glyph_name, reset_camera=False)
                            return _g
                        try:
                            moment_glyph = self._get_moment_glyph_geom()
                            _g = None
                            if moment_glyph is not None:
                                try:
                                    _g = _try_build_moment_glyph_inline(moment_glyph)
                                except Exception as _exc_c:
                                    perf_event('glyph', 'moment_glyph_failed',
                                               sid=sid,
                                               reason='composite_attempt_inline',
                                               exc=f"{_exc_c.__class__.__name__}: {str(_exc_c)[:120]}")
                                    _g = None
                            if _g is None:
                                fb_arrow = pv.Arrow(
                                    start=(0.0, 0.0, 0.0),
                                    direction=(1.0, 0.0, 0.0),
                                    tip_length=0.2, tip_radius=0.08,
                                    shaft_radius=0.03)
                                try:
                                    _g = _try_build_moment_glyph_inline(fb_arrow)
                                    perf_event('glyph', 'moment_fallback_arrow',
                                               sid=sid)
                                except Exception as _exc_f:
                                    perf_event('glyph', 'moment_build_failed',
                                               sid=sid,
                                               reason='inline_fallback_failed',
                                               exc=f"{_exc_f.__class__.__name__}: {str(_exc_f)[:120]}")
                        except Exception as _exc:
                            # Catch-all so a glyph build failure for one
                            # SID doesn't kill the post-parse visibility
                            # for the whole deck.
                            perf_event('glyph', 'moment_build_failed',
                                       sid=sid,
                                       exc=f"{_exc.__class__.__name__}: {str(_exc)[:120]}")
                    elif type == 'pressure':
                        vectors = dataset['vectors']; mags = np.linalg.norm(vectors, axis=1)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                        uniform_size = base_scale * 0.75
                        if use_relative_scaling and max_mag > 1e-9:
                            scaled_vectors = unit_vectors * (mags / max_mag * uniform_size)[:, np.newaxis]
                        else:
                            scaled_vectors = unit_vectors * uniform_size
                        self.plotter.add_arrows(dataset.points, scaled_vectors, color='yellow', name=glyph_name, reset_camera=False)
        perf_event('visibility', 'load_sid_loop',
                   n_sids=len(model.loads),
                   n_visible=_n_load_visible,
                   wall_s=f"{time.perf_counter()-_t_load_loop:.3f}",
                   find_tree_calls=getattr(self, '_find_tree_calls_total', 0) - _find_tree_calls_before)

        _t_spc_loop = time.perf_counter()
        _n_spc_visible = 0
        _find_tree_calls_before = getattr(self, '_find_tree_calls_total', 0)
        # v4.0.9 / v4.0.10: short-circuit + bounded iteration mirror
        # of the load_sid_loop optimization above. Reads
        # self._visible_spc_sids (maintained by _fast_path_sid_toggle
        # and refreshed in _refresh_visible_sid_sets) instead of
        # walking model.spcs.keys() to check checkboxes.
        _short_circuit_spcs = False
        if model.spcs:
            _any_stale_constraint_glyph = any(
                isinstance(name, str) and name.startswith('constraint_glyphs_')
                for name in self.plotter.actors.keys())
            _short_circuit_spcs = (not self._visible_spc_sids) and (not _any_stale_constraint_glyph)
        if _short_circuit_spcs:
            perf_event('spc_sid_loop', 'short_circuit',
                       n_sids=len(model.spcs),
                       taken=True,
                       wall_s=f"{time.perf_counter()-_t_spc_loop:.3f}")
        elif model.spcs:
            perf_event('spc_sid_loop', 'short_circuit',
                       n_sids=len(model.spcs),
                       n_visible=len(self._visible_spc_sids),
                       taken=False)
            # v4.0.10: stale-glyph sweep for constraint glyphs.
            _stale_removed_spc = 0
            for actor_name in list(self.plotter.actors.keys()):
                if not isinstance(actor_name, str):
                    continue
                if actor_name.startswith('constraint_glyphs_'):
                    try:
                        _stale_sid = int(actor_name[len('constraint_glyphs_'):])
                    except ValueError:
                        continue
                    if _stale_sid not in self._visible_spc_sids:
                        self.plotter.remove_actor(actor_name, render=False)
                        _stale_removed_spc += 1
            if _stale_removed_spc:
                perf_event('spc_sid_loop', 'stale_glyphs_swept',
                           n_removed=_stale_removed_spc)
        if not _short_circuit_spcs and model.spcs:
            perf_event('spc_sid_loop', 'iter_bounded',
                       n_iter=len(self._visible_spc_sids),
                       n_total=len(model.spcs))
        for sid in (self._visible_spc_sids if not _short_circuit_spcs else ()):
            is_visible = False
            if (sid_item_list := self._find_tree_items(('constraint_set', sid))) and sid_item_list[0].checkState(0) == QtCore.Qt.Checked:
                is_visible = True
                _n_spc_visible += 1

            glyph_name = f"constraint_glyphs_{sid}"
            # v4.0.6: only call remove_actor if the glyph actually
            # exists.
            if glyph_name in self.plotter.actors:
                self.plotter.remove_actor(glyph_name, render=False)
            if (actor := self.plotter.actors.get(f"constraint_actors_{sid}")) and is_visible:
                dataset = actor.mapper.dataset
                if not dataset.n_points: continue
                # v4.0.7: constraint triangle size -25% (was 0.02 →
                # now 0.015). User preference for default glyph scale.
                glyph_geom = pv.Cone(direction=(0, 0, 1), height=1.0, radius=0.3); scale_factor = self.current_grid.length * 0.015
                glyphs = dataset.glyph(scale=False, factor=scale_factor, orient=False, geom=glyph_geom)
                self.plotter.add_mesh(glyphs, color='cyan', name=glyph_name, reset_camera=False)
        perf_event('visibility', 'spc_sid_loop',
                   n_sids=len(model.spcs),
                   n_visible=_n_spc_visible,
                   wall_s=f"{time.perf_counter()-_t_spc_loop:.3f}",
                   find_tree_calls=getattr(self, '_find_tree_calls_total', 0) - _find_tree_calls_before)

        # --- Control visibility of coordinate system actors ---
        with perf_stage('visibility', 'refresh_coord_actors'):
            self._refresh_coord_actors(render=False)

        with perf_stage('visibility', 'bush_orientation_glyphs'):
            self._create_bush_orientation_glyphs(visibility_mask)
        with perf_stage('visibility', 'beam_orientation_glyphs'):
            self._create_beam_orientation_glyphs(visibility_mask)
        if self.show_free_edges_check.isChecked():
            with perf_stage('visibility', 'free_edges_actor'):
                self._create_free_edges_actor()
        with perf_stage('visibility', 'render_geometry'):
            self._render_geometry()
        # If cross-section is active, push the clipping plane into the
        # mappers of the freshly-rebuilt actors so the section stays on.
        if hasattr(self, "_cross_section") and self._cross_section.is_enabled:
            with perf_stage('visibility', 'cross_section_reapply'):
                self._cross_section.reapply()
        with perf_stage('visibility', 'final_render'):
            self.plotter.render()
        perf_event('visibility', 'total',
                   wall_s=f"{time.perf_counter()-_t_total:.3f}",
                   find_tree_calls_in_this_call=getattr(self, '_find_tree_calls_total', 0) - getattr(self, '_find_tree_calls_at_entry', 0))
        self._emit_mem_event('visibility.cycle')

    # --- Geometry visualization ---

    def _render_geometry(self):
        """Render geometry entities (points, lines, arcs, circles) in the viewer."""
        # Remove old geometry actors
        actors_to_remove = [name for name in self.plotter.actors if name.startswith('geom_')]
        for name in actors_to_remove:
            self.plotter.remove_actor(name, render=False)

        if self.geometry_store.is_empty:
            return

        # Check tree visibility
        geom_items = self._find_tree_items(('group', 'geometry'))
        if geom_items and geom_items[0].checkState(0) == QtCore.Qt.Unchecked:
            return

        # Points group visibility
        pts_visible = True
        pts_group = self._find_tree_items(('geom_group', 'points'))
        if pts_group and pts_group[0].checkState(0) == QtCore.Qt.Unchecked:
            pts_visible = False

        if pts_visible and self.geometry_store.points:
            pts_array = np.array([pt.xyz for pt in self.geometry_store.points.values()])
            self.plotter.add_points(pts_array, color='orange', point_size=8,
                                    render_points_as_spheres=True, name='geom_points')

        # Curves group visibility
        curves_visible = True
        curves_group = self._find_tree_items(('geom_group', 'curves'))
        if curves_group and curves_group[0].checkState(0) == QtCore.Qt.Unchecked:
            curves_visible = False

        if curves_visible:
            for lid, line in self.geometry_store.lines.items():
                pts = line.evaluate_array(2)
                poly = pv.Line(pts[0], pts[1])
                self.plotter.add_mesh(poly, color='magenta', line_width=3,
                                      name=f'geom_line_{lid}', render=False)

            for aid, arc in self.geometry_store.arcs.items():
                pts = arc.evaluate_array(32)
                poly = pv.PolyData(pts)
                poly.lines = np.hstack([[32] + list(range(32))])
                self.plotter.add_mesh(poly, color='cyan', line_width=3,
                                      name=f'geom_arc_{aid}', render=False)

            for cid, circle in self.geometry_store.circles.items():
                pts = circle.evaluate_array(64)
                # Close the circle
                connectivity = list(range(64)) + [0]
                poly = pv.PolyData(pts)
                poly.lines = np.hstack([[65] + connectivity])
                self.plotter.add_mesh(poly, color='cyan', line_width=3,
                                      name=f'geom_circle_{cid}', render=False)

        # Surfaces - just render boundary curves for now (fill deferred)
        surfs_visible = True
        surfs_group = self._find_tree_items(('geom_group', 'surfaces'))
        if surfs_group and surfs_group[0].checkState(0) == QtCore.Qt.Unchecked:
            surfs_visible = False

        if surfs_visible and self.geometry_store.surfaces:
            for sid, surf in self.geometry_store.surfaces.items():
                for curve_id in surf.boundary_curve_ids:
                    try:
                        pts = self.geometry_store.evaluate_curve(curve_id, 32)
                        poly = pv.PolyData(pts)
                        poly.lines = np.hstack([[32] + list(range(32))])
                        self.plotter.add_mesh(poly, color='orange', line_width=2,
                                              opacity=0.7, name=f'geom_surf_{sid}_curve_{curve_id}',
                                              render=False)
                    except KeyError:
                        pass

    # --- Geometry creation handlers ---

    def _create_geometry_point(self):
        next_id = self.geometry_store.next_point_id()
        dialog = CreateGeometryPointDialog(next_id, self)
        if dialog.exec():
            if dialog.is_batch:
                batch = dialog.get_batch_parameters()
                cmds = []
                for pt in batch:
                    pid = pt['id'] if pt['id'] else self.geometry_store.next_point_id()
                    cmd = AddGeometryPointCommand(
                        self.geometry_store, pid,
                        [pt['x'], pt['y'], pt['z']])
                    cmds.append(cmd)
                compound = CompoundCommand(cmds, f"Import {len(cmds)} geometry points from Excel")
                compound.execute(None)
                self.command_manager.undo_stack.append(compound)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Imported {len(cmds)} geometry points from Excel")
            else:
                params = dialog.get_parameters()
                cmd = AddGeometryPointCommand(
                    self.geometry_store, params['id'],
                    [params['x'], params['y'], params['z']])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Created geometry point P{params['id']}")

    def _create_geometry_line(self):
        if len(self.geometry_store.points) < 2:
            QMessageBox.warning(self, "Insufficient Points",
                                "Need at least 2 geometry points to create a line.")
            return
        next_id = self.geometry_store.next_curve_id()
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = CreateGeometryLineDialog(next_id, point_ids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            cmd = AddGeometryLineCommand(
                self.geometry_store, params['id'],
                params['start_point_id'], params['end_point_id'])
            cmd.execute(None)
            self.command_manager.undo_stack.append(cmd)
            self.command_manager.redo_stack.clear()
            self._populate_tree()
            self._update_plot_visibility()
            self._update_status(f"Created geometry line {params['id']}")

    def _create_geometry_arc(self):
        if len(self.geometry_store.points) < 3:
            QMessageBox.warning(self, "Insufficient Points",
                                "Need at least 3 geometry points to create an arc.")
            return
        next_id = self.geometry_store.next_curve_id()
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = CreateGeometryArcDialog(next_id, point_ids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            try:
                cmd = AddGeometryArcCommand(
                    self.geometry_store, params['id'],
                    params['start_point_id'], params['mid_point_id'],
                    params['end_point_id'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Created geometry arc {params['id']}")
            except ValueError as e:
                QMessageBox.critical(self, "Arc Error", str(e))

    def _create_geometry_circle(self):
        if not self.geometry_store.points:
            QMessageBox.warning(self, "No Points",
                                "Need at least 1 geometry point for the circle center.")
            return
        next_id = self.geometry_store.next_curve_id()
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = CreateGeometryCircleDialog(next_id, point_ids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            cmd = AddGeometryCircleCommand(
                self.geometry_store, params['id'],
                params['center_point_id'], params['radius'],
                [params['nx'], params['ny'], params['nz']])
            cmd.execute(None)
            self.command_manager.undo_stack.append(cmd)
            self.command_manager.redo_stack.clear()
            self._populate_tree()
            self._update_plot_visibility()
            self._update_status(f"Created geometry circle {params['id']}")

    def _create_geometry_surface(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.warning(self, "No Curves",
                                "Need at least 1 curve to create a surface.")
            return
        next_id = self.geometry_store.next_surface_id()
        dialog = CreateGeometrySurfaceDialog(next_id, all_curves, self)
        if dialog.exec():
            params = dialog.get_parameters()
            cmd = AddGeometrySurfaceCommand(
                self.geometry_store, params['id'], params['curve_ids'])
            cmd.execute(None)
            self.command_manager.undo_stack.append(cmd)
            self.command_manager.redo_stack.clear()
            self._populate_tree()
            self._update_plot_visibility()
            self._update_status(f"Created geometry surface {params['id']}")

    # --- Geometry modify handlers ---

    def _edit_geometry_point(self):
        if not self.geometry_store.points:
            QMessageBox.information(self, "No Points", "No geometry points to edit.")
            return
        point_ids = sorted(self.geometry_store.points.keys())
        items = [f"P{pid}" for pid in point_ids]
        item, ok = QInputDialog.getItem(self, "Edit Geometry Point",
                                         "Select point:", items, 0, False)
        if ok and item:
            pid = int(item[1:])
            pt = self.geometry_store.points[pid]
            from node_runner.dialogs.geometry import EditGeometryPointDialog
            dialog = EditGeometryPointDialog(pid, pt.xyz, self)
            if dialog.exec():
                params = dialog.get_parameters()
                try:
                    cmd = EditGeometryPointCommand(
                        self.geometry_store, pid,
                        [params['x'], params['y'], params['z']])
                    cmd.execute(None)
                    self.command_manager.undo_stack.append(cmd)
                    self.command_manager.redo_stack.clear()
                    self._populate_tree()
                    self._update_plot_visibility()
                    self._update_status(f"Edited geometry point P{pid}")
                except ValueError as e:
                    QMessageBox.critical(self, "Edit Error", str(e))

    def _edit_geometry_curve(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.information(self, "No Curves", "No geometry curves to edit.")
            return
        items = [f"Curve {cid}" for cid in all_curves]
        item, ok = QInputDialog.getItem(self, "Edit Geometry Curve",
                                         "Select curve:", items, 0, False)
        if ok and item:
            cid = int(item.split()[1])
            point_ids = sorted(self.geometry_store.points.keys())
            if cid in self.geometry_store.lines:
                line = self.geometry_store.lines[cid]
                current_data = {
                    'start_point_id': line.start_point_id,
                    'end_point_id': line.end_point_id,
                }
                curve_type = 'line'
            elif cid in self.geometry_store.arcs:
                arc = self.geometry_store.arcs[cid]
                current_data = {
                    'start_point_id': arc.start_point_id,
                    'mid_point_id': arc.mid_point_id,
                    'end_point_id': arc.end_point_id,
                }
                curve_type = 'arc'
            elif cid in self.geometry_store.circles:
                circle = self.geometry_store.circles[cid]
                current_data = {
                    'center_point_id': circle.center_point_id,
                    'radius': circle.radius,
                    'normal': circle.normal.tolist(),
                }
                curve_type = 'circle'
            else:
                return

            from node_runner.dialogs.geometry import EditGeometryCurveDialog
            dialog = EditGeometryCurveDialog(
                cid, curve_type, point_ids, current_data, self)
            if dialog.exec():
                params = dialog.get_parameters()
                try:
                    cmd = EditGeometryCurveCommand(
                        self.geometry_store, cid, params)
                    cmd.execute(None)
                    self.command_manager.undo_stack.append(cmd)
                    self.command_manager.redo_stack.clear()
                    self._populate_tree()
                    self._update_plot_visibility()
                    self._update_status(f"Edited geometry curve {cid}")
                except ValueError as e:
                    QMessageBox.critical(self, "Edit Error", str(e))

    def _edit_geometry_surface(self):
        if not self.geometry_store.surfaces:
            QMessageBox.information(self, "No Surfaces",
                                    "No geometry surfaces to edit.")
            return
        surface_ids = sorted(self.geometry_store.surfaces.keys())
        items = [f"Surface {sid}" for sid in surface_ids]
        item, ok = QInputDialog.getItem(self, "Edit Geometry Surface",
                                         "Select surface:", items, 0, False)
        if ok and item:
            sid = int(item.split()[1])
            surf = self.geometry_store.surfaces[sid]
            all_curves = self.geometry_store.all_curve_ids()
            from node_runner.dialogs.geometry import EditGeometrySurfaceDialog
            dialog = EditGeometrySurfaceDialog(
                sid, all_curves, surf.boundary_curve_ids, self)
            if dialog.exec():
                params = dialog.get_parameters()
                cmd = EditGeometrySurfaceBoundariesCommand(
                    self.geometry_store, sid, params['curve_ids'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Edited geometry surface {sid} boundaries")

    # --- Geometry transform handlers ---

    def _transform_geometry(self, default_tab='translate'):
        if not self.geometry_store.points:
            QMessageBox.information(self, "No Points",
                                    "No geometry points to transform.")
            return
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = EntitySelectionDialog("Geometry Point", set(point_ids), self)
        if dialog.exec():
            selected = dialog.get_selected_ids()
            if not selected:
                return
            from node_runner.dialogs.geometry import TransformGeometryDialog
            t_dialog = TransformGeometryDialog(selected, self)
            tab_map = {'translate': 0, 'rotate': 1, 'mirror': 2, 'scale': 3}
            t_dialog.tabs.setCurrentIndex(tab_map.get(default_tab, 0))
            if t_dialog.exec():
                params = t_dialog.get_parameters()
                try:
                    cmd = TransformGeometryPointsCommand(
                        self.geometry_store, selected, params)
                    cmd.execute(None)
                    self.command_manager.undo_stack.append(cmd)
                    self.command_manager.redo_stack.clear()
                    self._populate_tree()
                    self._update_plot_visibility()
                    self._update_status(
                        f"{params['type'].capitalize()}d {len(selected)} "
                        f"geometry point(s)")
                except ValueError as e:
                    QMessageBox.critical(self, "Transform Error", str(e))

    def _copy_geometry(self):
        if not self.geometry_store.points:
            QMessageBox.information(self, "No Points",
                                    "No geometry points to copy.")
            return
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = EntitySelectionDialog("Geometry Point", set(point_ids), self)
        if dialog.exec():
            selected = dialog.get_selected_ids()
            if not selected:
                return
            from node_runner.dialogs.geometry import CopyGeometryDialog
            c_dialog = CopyGeometryDialog(selected, self)
            if c_dialog.exec():
                params = c_dialog.get_parameters()
                cmd = CopyGeometryCommand(
                    self.geometry_store, selected, params['delta'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(
                    f"Copied {len(selected)} geometry point(s) "
                    f"and dependent curves")

    # --- CAD-like geometry operation handlers ---

    def _split_curve(self):
        all_curves = self.geometry_store.all_curve_ids()
        # Exclude circles (cannot split full circles)
        splittable = [cid for cid in all_curves
                      if cid not in self.geometry_store.circles]
        if not splittable:
            QMessageBox.information(self, "No Curves",
                                    "No splittable curves (lines/arcs).")
            return
        from node_runner.dialogs.geometry import SplitCurveDialog
        dialog = SplitCurveDialog(splittable, self)
        if dialog.exec():
            params = dialog.get_parameters()
            try:
                cmd = SplitCurveCommand(
                    self.geometry_store,
                    params['curve_id'], params['n_segments'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(
                    f"Split curve {params['curve_id']} into "
                    f"{params['n_segments']} segments")
            except ValueError as e:
                QMessageBox.critical(self, "Split Error", str(e))

    def _offset_curve(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.information(self, "No Curves",
                                    "No geometry curves to offset.")
            return
        from node_runner.dialogs.geometry import OffsetCurveDialog
        dialog = OffsetCurveDialog(all_curves, self)
        if dialog.exec():
            params = dialog.get_parameters()
            try:
                cmd = OffsetCurveCommand(
                    self.geometry_store,
                    params['curve_id'], params['distance'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(
                    f"Offset curve {params['curve_id']} "
                    f"by {params['distance']}")
            except ValueError as e:
                QMessageBox.critical(self, "Offset Error", str(e))

    def _project_point_to_curve(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.information(self, "No Curves",
                                    "No geometry curves to project onto.")
            return
        from node_runner.dialogs.geometry import ProjectPointToCurveDialog
        dialog = ProjectPointToCurveDialog(all_curves, self)
        if dialog.exec():
            params = dialog.get_parameters()
            cmd = ProjectPointToCurveCommand(
                self.geometry_store,
                params['curve_id'], params['source_xyz'])
            cmd.execute(None)
            self.command_manager.undo_stack.append(cmd)
            self.command_manager.redo_stack.clear()
            self._populate_tree()
            self._update_plot_visibility()
            self._update_status(
                f"Projected point onto curve {params['curve_id']}")

    def _fillet_curves(self):
        line_ids = sorted(self.geometry_store.lines.keys())
        if len(line_ids) < 2:
            QMessageBox.information(self, "Insufficient Lines",
                                    "Need at least 2 geometry lines for a fillet.")
            return
        from node_runner.dialogs.geometry import FilletDialog
        dialog = FilletDialog(line_ids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            try:
                cmd = FilletCommand(
                    self.geometry_store,
                    params['line_id_1'], params['line_id_2'],
                    params['radius'])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(
                    f"Created fillet R={params['radius']} between "
                    f"lines {params['line_id_1']} and {params['line_id_2']}")
            except ValueError as e:
                QMessageBox.critical(self, "Fillet Error", str(e))

    # --- Geometry meshing handlers ---

    def _mesh_curve(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.warning(self, "No Curves", "No geometry curves to mesh.")
            return
        if not self.current_generator or not self.current_generator.model.properties:
            QMessageBox.warning(self, "No Properties",
                                "Create a line element property first.")
            return
        model = self.current_generator.model
        line_prop_ids = []
        for pid, prop in model.properties.items():
            if prop.type in ('PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PROD', 'PBUSH'):
                line_prop_ids.append(pid)
        if not line_prop_ids:
            QMessageBox.warning(self, "No Line Properties",
                                "No bar/beam/rod/bush properties found.")
            return
        dialog = MeshCurveDialog(all_curves, line_prop_ids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            cmd = MeshCurveCommand(
                self.geometry_store, params['curve_ids'],
                params['n_elements'], params['elem_type'], params['pid'])
            self.command_manager.execute(cmd, model)
            self._update_viewer(self.current_generator)
            self._update_status(
                f"Meshed {len(params['curve_ids'])} curve(s) with "
                f"{params['n_elements']} elements each")

    def _mesh_surface(self):
        surface_ids = sorted(self.geometry_store.surfaces.keys())
        if not surface_ids:
            QMessageBox.warning(self, "No Surfaces",
                                "No geometry surfaces to mesh.")
            return
        if not self.current_generator:
            return
        model = self.current_generator.model
        shell_pids = []
        for pid, prop in model.properties.items():
            if prop.type in ('PSHELL', 'PCOMP'):
                shell_pids.append(pid)
        if not shell_pids:
            QMessageBox.warning(self, "No Shell Properties",
                                "Create a PSHELL or PCOMP property first.")
            return
        dialog = MeshSurfaceDialog(surface_ids, shell_pids, self)
        if dialog.exec():
            params = dialog.get_parameters()
            try:
                cmd = MeshSurfaceCommand(
                    self.geometry_store, params['surface_id'],
                    params['n_per_edge'], params['pid'],
                    params['elem_preference'])
                self.command_manager.execute(cmd, model)
                self._update_viewer(self.current_generator)
                n_elems = len(cmd._created_eids)
                n_nodes = len(cmd._created_nids)
                self._update_status(
                    f"Meshed surface {params['surface_id']}: "
                    f"{n_nodes} nodes, {n_elems} elements")
            except Exception as e:
                QMessageBox.critical(self, "Mesh Error", str(e))

    def _nodes_at_geometry_points(self):
        if not self.geometry_store.points:
            QMessageBox.warning(self, "No Points", "No geometry points exist.")
            return
        if not self.current_generator:
            return
        model = self.current_generator.model
        point_ids = sorted(self.geometry_store.points.keys())
        dialog = EntitySelectionDialog("Geometry Point", set(point_ids), self)
        if dialog.exec():
            selected = dialog.get_selected_ids()
            if not selected:
                return
            cmd = NodesAtGeometryPointsCommand(self.geometry_store, selected)
            self.command_manager.execute(cmd, model)
            self._update_viewer(self.current_generator)
            self._update_status(f"Created {len(selected)} node(s) at geometry points")

    # --- Geometry deletion handlers ---

    def _delete_geometry_points(self):
        if not self.geometry_store.points:
            QMessageBox.information(self, "No Points", "No geometry points to delete.")
            return
        point_ids = sorted(self.geometry_store.points.keys())
        items = [f"P{pid}" for pid in point_ids]
        item, ok = QInputDialog.getItem(self, "Delete Geometry Point",
                                         "Select point to delete:", items, 0, False)
        if ok and item:
            pid = int(item[1:])
            dep_curves = self.geometry_store.curves_referencing_point(pid)
            dep_surfs = self.geometry_store.surfaces_referencing_curves(set(dep_curves))
            msg = f"Delete geometry point P{pid}?"
            if dep_curves:
                msg += f"\n\nThis will also delete {len(dep_curves)} dependent curve(s)"
            if dep_surfs:
                msg += f" and {len(dep_surfs)} dependent surface(s)"
            msg += "."
            reply = QMessageBox.question(self, "Confirm Deletion", msg,
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteGeometryPointsCommand(self.geometry_store, [pid])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Deleted geometry point P{pid}")

    def _delete_geometry_curves(self):
        all_curves = self.geometry_store.all_curve_ids()
        if not all_curves:
            QMessageBox.information(self, "No Curves", "No geometry curves to delete.")
            return
        items = [f"Curve {cid}" for cid in all_curves]
        item, ok = QInputDialog.getItem(self, "Delete Geometry Curve",
                                         "Select curve to delete:", items, 0, False)
        if ok and item:
            cid = int(item.split()[-1])
            dep_surfs = self.geometry_store.surfaces_referencing_curve(cid)
            msg = f"Delete geometry curve {cid}?"
            if dep_surfs:
                msg += f"\n\nThis will also delete {len(dep_surfs)} dependent surface(s)."
            reply = QMessageBox.question(self, "Confirm Deletion", msg,
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteGeometryCurvesCommand(self.geometry_store, [cid])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Deleted geometry curve {cid}")

    def _delete_geometry_surfaces(self):
        if not self.geometry_store.surfaces:
            QMessageBox.information(self, "No Surfaces", "No geometry surfaces to delete.")
            return
        surface_ids = sorted(self.geometry_store.surfaces.keys())
        items = [f"Surface {sid}" for sid in surface_ids]
        item, ok = QInputDialog.getItem(self, "Delete Geometry Surface",
                                         "Select surface to delete:", items, 0, False)
        if ok and item:
            sid = int(item.split()[-1])
            reply = QMessageBox.question(self, "Confirm Deletion",
                                          f"Delete geometry surface {sid}?",
                                          QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                cmd = DeleteGeometrySurfacesCommand(self.geometry_store, [sid])
                cmd.execute(None)
                self.command_manager.undo_stack.append(cmd)
                self.command_manager.redo_stack.clear()
                self._populate_tree()
                self._update_plot_visibility()
                self._update_status(f"Deleted geometry surface {sid}")

    # --- Analysis (SOL/EIGRL) handlers ---

    def _set_sol_type(self, sol):
        self.sol_type = sol
        label = f"SOL {sol}" if sol else "None (Punch file)"
        self._update_status(f"Solution type set to {label}")

    def _create_eigrl(self):
        from node_runner.dialogs.loads import CreateEigrlDialog
        dialog = CreateEigrlDialog(self)
        if dialog.exec():
            params = dialog.get_parameters()
            self.eigrl_cards.append(params)
            self._update_status(f"Created EIGRL SID {params['sid']} (ND={params['nd']})")

    def _open_analysis_set_manager(self):
        """Open the Analysis Set Manager dialog.

        v5.3.0 item 53: no longer auto-creates a "Default" AnalysisSet
        when none exist. The Manager opens empty; the user clicks
        ``New`` to create one explicitly. This matches common-tool behavior
        (no implicit/punch set) and removes the silent default that
        was leaking into Run Analysis dialogs.
        """
        from node_runner.dialogs.analysis import AnalysisSetManagerDialog, AnalysisSet

        model = self.current_generator.model if self.current_generator else None
        dialog = AnalysisSetManagerDialog(
            self.analysis_sets, self.active_analysis_set_id,
            model=model, groups=self.groups, parent=self)
        if dialog.exec():
            self.analysis_sets, self.active_analysis_set_id = dialog.get_results()
            active = self.analysis_sets.get(self.active_analysis_set_id)
            if active:
                # Sync legacy fields from the active analysis set
                self.sol_type = active.sol_type
                self.subcases = active.subcases
                if active.eigrl:
                    # Replace/add EIGRL card
                    self.eigrl_cards = [active.eigrl]
                self._update_status(
                    f"Analysis Set Manager: {len(self.analysis_sets)} set(s), "
                    f"active = '{active.name}'")
            else:
                self._update_status(
                    f"Analysis Set Manager: {len(self.analysis_sets)} set(s)")
            self._populate_tree()

    def _open_output_requests_dialog(self):
        """Open the quick Output Requests dialog for the active analysis set.

        v5.3.0 item 53: no longer auto-creates a Default set. If the
        user has none defined, point them at the Manager.
        """
        from node_runner.dialogs.analysis import OutputRequestsDialog, AnalysisSet
        if not self.analysis_sets:
            QMessageBox.information(
                self, "No AnalysisSet Defined",
                "No AnalysisSet exists yet. Open the Analysis Set "
                "Manager (Analysis -> Analysis Set Manager…) and click "
                "New to create one before editing output requests.")
            return

        active = self.analysis_sets.get(self.active_analysis_set_id)
        if not active:
            QMessageBox.information(self, "No Active Set",
                                    "No active analysis set. Use the Analysis Set Manager first.")
            return

        dialog = OutputRequestsDialog(active.output_requests, parent=self)
        if dialog.exec():
            active.output_requests = dialog.get_output_requests()
            self._update_status(f"Output requests updated for '{active.name}'")

    # ─── v5.0.0 items 17/18: MYSTRAN run flow ─────────────────────────
    def _configure_mystran(self):
        """Analysis -> Configure MYSTRAN…: open Preferences on the MYSTRAN tab."""
        from node_runner.dialogs.preferences import (
            PreferencesDialog, save_preferences, load_preferences)
        current = load_preferences()
        dlg = PreferencesDialog(current, self)
        # Try to switch to the MYSTRAN tab.
        try:
            for i in range(dlg._tabs.count()):
                if dlg._tabs.tabText(i) == "MYSTRAN":
                    dlg._tabs.setCurrentIndex(i)
                    break
        except Exception:
            pass
        if dlg.exec() == QDialog.Accepted:
            save_preferences(dlg.result_payload())
            self._update_status("MYSTRAN preferences saved.")

    def _open_analysis_history(self):
        """Analysis -> Analysis History…"""
        from node_runner.dialogs.solver_run import AnalysisHistoryDialog
        from node_runner.solve import mystran_settings
        dlg = AnalysisHistoryDialog(mystran_settings.get_scratch_root(),
                                    parent=self)
        dlg.reload_requested.connect(self._on_mystran_history_reload)
        dlg.exec()

    def _on_mystran_history_reload(self, run):
        """User picked a prior run from the history dialog. Re-load
        its results into the current model.
        """
        try:
            from node_runner.solve.mystran_results import load_mystran_results
            bundle = load_mystran_results(run)
            if bundle is None:
                QMessageBox.warning(
                    self, "Analysis History",
                    "Couldn't load OP2 or F06 from that run folder.")
                return
            self._consume_mystran_bundle(bundle, run)
            # v5.0.0 item 20: update the persisted results_source so a
            # prior run that previously read "null" now records the
            # source that actually loaded (OP2 vs F06).
            try:
                if run.bdf_path:
                    run.save_meta(Path(run.bdf_path).parent)
            except Exception:
                pass
        except Exception as exc:
            QMessageBox.critical(
                self, "Analysis History",
                f"Failed to reload that run: {exc}")

    def _run_mystran(self):
        """Analysis -> Run Analysis (MYSTRAN)…: pre-flight -> run-options
        -> solver progress -> auto-load results.

        v5.1.0 item 26: requires an active AnalysisSet. The set drives
        what gets exported (scope/group/load/SPC filtering via
        scope.py) and what shows up in the deck's case-control block.
        common convention -- you can't just "run a model"; you run an
        analysis set.
        """
        if not self.current_generator or not self.current_generator.model:
            QMessageBox.warning(
                self, "Run MYSTRAN",
                "Open or create a model first.")
            return

        # v5.1.0 item 26: block-with-pointer when no AnalysisSet is
        # active. The user goes to AnalysisSet Manager to create or
        # pick one -- we no longer silently run the full model.
        active_set = (self.analysis_sets or {}).get(
            self.active_analysis_set_id)
        if active_set is None:
            ret = QMessageBox.question(
                self, "Run requires an active AnalysisSet",
                "MYSTRAN runs need an active AnalysisSet so the scope "
                "(group target, load/SPC SIDs, solver target) is "
                "explicit.\n\nOpen the AnalysisSet Manager now?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ret == QMessageBox.Yes:
                self._open_analysis_set_manager()
                active_set = (self.analysis_sets or {}).get(
                    self.active_analysis_set_id)
            if active_set is None:
                self._update_status(
                    "MYSTRAN run cancelled -- no active AnalysisSet.")
                return
        # Confirm the AnalysisSet's solver target is MYSTRAN.
        target = getattr(active_set, 'solver_target', 'MYSTRAN')
        if target != 'MYSTRAN':
            ret = QMessageBox.question(
                self, "AnalysisSet target mismatch",
                f"Active AnalysisSet '{active_set.name}' has "
                f"solver_target='{target}', not MYSTRAN. Run anyway?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ret != QMessageBox.Yes:
                return

        # Verify MYSTRAN binary is configured.
        from node_runner.solve import mystran_settings
        from node_runner.solve.mystran_runner import (
            discover_mystran_executable)
        exe = discover_mystran_executable()
        if not exe:
            ret = QMessageBox.question(
                self, "MYSTRAN not configured",
                "I couldn't find the MYSTRAN executable. Configure it now?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ret == QMessageBox.Yes:
                self._configure_mystran()
                exe = discover_mystran_executable()
            if not exe:
                return

        # Open Run-options dialog.
        from node_runner.dialogs.solver_run import (
            RunMystranDialog, MystranPreflightDialog, SolverProgressDialog)
        run_dlg = RunMystranDialog(self, parent=self)
        if run_dlg.exec() != QDialog.Accepted:
            return
        opts = run_dlg.get_options()

        # v5.1.0 item 26: build the scoped model for pre-flight + export.
        from node_runner.scope import scope_model_to_analysis_set
        try:
            scoped_model, scope_report = scope_model_to_analysis_set(
                self.current_generator.model, active_set,
                groups=self.groups)
        except Exception as exc:
            QMessageBox.critical(
                self, "Scoping failed",
                f"Could not scope model to AnalysisSet '{active_set.name}': "
                f"{exc}")
            return

        # Pre-flight (on the scoped model so users see only relevant
        # blocking issues).
        from node_runner.solve.mystran_preflight import scan_for_mystran
        report = scan_for_mystran(scoped_model, sol=opts.get('sol'))
        pf_dlg = MystranPreflightDialog(report, parent=self)
        if pf_dlg.exec() != QDialog.Accepted:
            return

        # Export deck to working dir. Run-folder name carries the
        # AnalysisSet name so Analysis History distinguishes sets.
        safe_setname = ''.join(c if c.isalnum() or c in '-_' else '_'
                                for c in active_set.name)[:32]
        working_dir = Path(opts['working_dir'])
        # Replace the trailing timestamp-only folder with one that
        # also carries the set name.
        if safe_setname and not working_dir.name.endswith(safe_setname):
            working_dir = working_dir.with_name(
                f"{working_dir.name}_{safe_setname}")
        working_dir.mkdir(parents=True, exist_ok=True)
        bdf_stem = Path(self._loaded_filepath).stem if self._loaded_filepath \
            else "untitled"
        bdf_path = working_dir / f"{bdf_stem}.bdf"
        try:
            self.current_generator._write_bdf(
                str(bdf_path),
                field_format=getattr(
                    active_set, 'field_format', 'short'),
                target='mystran',
                analysis_set=active_set,
                groups=self.groups,
                sollib=mystran_settings.get('default_sollib'),
                quad4typ=mystran_settings.get('default_quad4typ'),
                wtmass=mystran_settings.get('default_wtmass'),
            )
        except Exception as exc:
            QMessageBox.critical(
                self, "Export failed",
                f"Could not write MYSTRAN BDF: {exc}")
            return

        # Status line tells the user exactly what got scoped in.
        self._update_status(
            f"Running '{active_set.name}': "
            f"{scope_report.n_elements_kept:,} elements, "
            f"{scope_report.n_loads_kept} load SID(s), "
            f"{scope_report.n_spcs_kept} SPC SID(s)"
            + (f" (group: {scope_report.group_target})"
               if scope_report.group_target else ""))

        # Launch SolverRunWorker.
        from node_runner.solve.mystran_runner import SolverRunWorker
        self._solver_worker = SolverRunWorker(self)
        self._solver_dialog = SolverProgressDialog(bdf_path.name, parent=self)
        self._solver_worker.progress.connect(
            self._solver_dialog.update_progress)
        self._solver_worker.finished.connect(self._on_mystran_finished)
        self._solver_worker.failed.connect(self._on_mystran_failed)
        self._solver_worker.cancelled.connect(self._on_mystran_cancelled)
        self._solver_dialog.cancel_requested.connect(
            self._solver_worker.cancel)
        self._solver_dialog.start_timer()
        self._solver_dialog.show()
        self._solver_worker.start(
            bdf_path=bdf_path,
            exe_path=exe,
            output_dir=working_dir,
            sol=opts.get('sol', 101),
            analysis_set_name=active_set.name,
        )

    def _on_mystran_finished(self, run):
        """Solver completed cleanly. Auto-load results per preference,
        then show a single summary dialog that says where the files
        landed and surfaces any MYSTRAN-level errors from the ERR file.
        """
        from node_runner.solve import mystran_settings
        from node_runner.solve.mystran_results import load_mystran_results
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        self._update_status(
            f"MYSTRAN finished in {run.wall_time:.1f}s "
            f"(exit code {run.return_code}).")

        # Always try to load results -- the user expectation is "click
        # Run, see results". Status of the parse + the run folder are
        # surfaced together in a single post-run dialog below.
        bundle = None
        if mystran_settings.get('auto_load_results', True):
            bundle = load_mystran_results(run)
            if bundle and self._bundle_has_data(bundle):
                self._consume_mystran_bundle(bundle, run)
                self._mystran_results_loaded_msg = (
                    f"Loaded {self._bundle_disp_count(bundle):,} "
                    f"displacements ({run.results_source.upper()} source).")
            elif bundle is not None:
                # Bundle parsed but empty -- treat as failure visually.
                self._mystran_results_loaded_msg = (
                    "Results file parsed but contained NO displacement "
                    "data. MYSTRAN likely aborted with errors -- see "
                    "the ERR file below.")
            else:
                self._mystran_results_loaded_msg = (
                    "No OP2 or F06 results file found in the run folder.")
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished runs before
        # load_mystran_results() sets this field, so the on-disk
        # run_meta.json shows "results_source": null until we re-save
        # here. Affects Analysis History display.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished runs before
        # load_mystran_results() sets this field, so the on-disk
        # run_meta.json shows "results_source": null until we re-save
        # here. Affects Analysis History display.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished runs before
        # load_mystran_results() sets this field, so the on-disk
        # run_meta.json shows "results_source": null until we re-save
        # here. Affects Analysis History display.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished happens before
        # load_mystran_results() can set this field, so the on-disk
        # run_meta.json reads "results_source": null until we re-save
        # here. Affects Analysis History display + downstream tooling.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished happens before
        # load_mystran_results() can set this field, so the on-disk
        # run_meta.json reads "results_source": null until we re-save
        # here. Affects Analysis History display + downstream tooling.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished happens before
        # load_mystran_results() can set this field, so the on-disk
        # run_meta.json reads "results_source": null until we re-save
        # here. Affects Analysis History display + downstream tooling.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass
        else:
            # v5.1.1 item 30: this branch belongs to the outer
            # `if mystran_settings.get('auto_load_results', True):`
            # condition. In v5.0.0/v5.1.0 it was accidentally indented
            # under the `try` block below, which made it fire on every
            # successful save_meta() -- i.e. always -- and silently
            # overwrite the real "Loaded N displacements..." message.
            self._mystran_results_loaded_msg = (
                "Auto-load is disabled in Edit -> Preferences -> "
                "MYSTRAN. Click 'Load Results' below to import.")

        # v5.0.0 item 20: persist the final results_source. The first
        # save_meta() inside SolverRunWorker._on_finished runs before
        # load_mystran_results() sets this field, so the on-disk
        # run_meta.json shows "results_source": null until we re-save
        # here. Affects Analysis History display.
        try:
            if run.bdf_path:
                run.save_meta(Path(run.bdf_path).parent)
        except Exception:
            pass

        self._show_mystran_summary(run, bundle)

    def _bundle_has_data(self, bundle):
        """True iff the bundle contains at least one non-empty
        displacement / eigenvalue / eigenvector entry."""
        for sc in (bundle.get('subcases') or {}).values():
            if sc.get('displacements'):
                return True
            if sc.get('eigenvalues'):
                return True
            if sc.get('eigenvectors'):
                return True
        return False

    def _bundle_disp_count(self, bundle):
        return sum(len(sc.get('displacements', {}))
                   for sc in (bundle.get('subcases') or {}).values())

    def _show_mystran_summary(self, run, bundle):
        """Single post-run dialog: where the files are, MYSTRAN ERR
        contents (truncated), and an Open Folder button.

        Replaces the v5.0.0-first-cut behaviour where a successful exit
        code but no parsed results left the user staring at an
        unchanged viewport with no idea where to look.
        """
        from PySide6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
            QPlainTextEdit, QFrame, QDialogButtonBox,
        )
        from PySide6.QtGui import QFont

        run_dir = Path(run.bdf_path).parent if run.bdf_path else None
        err_path = (Path(run.bdf_path).with_suffix('.ERR')
                    if run.bdf_path else None)
        err_text = ""
        if err_path and err_path.exists():
            try:
                err_text = err_path.read_text(
                    encoding='utf-8', errors='replace')
            except Exception:
                err_text = ""
        # Count MYSTRAN-level errors in the ERR for the banner.
        n_err = err_text.upper().count("*ERROR")
        n_warn = err_text.upper().count("*WARNING")

        dlg = QDialog(self)
        dlg.setWindowTitle("MYSTRAN run summary")
        dlg.setMinimumSize(760, 480)

        layout = QVBoxLayout(dlg)

        # ---- Banner ----
        has_data = bool(bundle) and self._bundle_has_data(bundle)
        if has_data and n_err == 0:
            banner_text = (
                f"<b>Run successful.</b> {self._mystran_results_loaded_msg}")
            banner_style = (
                "background-color: #a6e3a1; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        elif has_data:
            banner_text = (
                f"<b>Run produced results, with {n_err} MYSTRAN error(s).</b> "
                f"{self._mystran_results_loaded_msg} Inspect the ERR "
                f"file below before trusting these numbers.")
            banner_style = (
                "background-color: #f9e2af; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        else:
            banner_text = (
                f"<b>MYSTRAN exited cleanly but produced no result data.</b><br>"
                f"{n_err} error(s), {n_warn} warning(s) in the ERR file. "
                f"{self._mystran_results_loaded_msg}")
            banner_style = (
                "background-color: #f38ba8; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        banner = QLabel(banner_text)
        banner.setStyleSheet(banner_style)
        banner.setWordWrap(True)
        layout.addWidget(banner)

        # ---- File summary ----
        files_box = QLabel(
            f"<b>Run folder:</b> <code>{run_dir}</code><br>"
            f"<b>BDF:</b> <code>{Path(run.bdf_path).name if run.bdf_path else '-'}</code> | "
            f"<b>F06:</b> <code>{Path(run.f06_path).name if run.f06_path else '(none)'}</code> | "
            f"<b>OP2:</b> <code>{Path(run.op2_path).name if run.op2_path else '(none)'}</code>"
        )
        files_box.setStyleSheet(
            "font-family: Consolas, monospace; font-size: 11px; "
            "padding: 6px 0;")
        files_box.setWordWrap(True)
        layout.addWidget(files_box)

        # ---- ERR contents ----
        layout.addWidget(QLabel("<b>MYSTRAN .ERR file:</b>"))
        err_view = QPlainTextEdit()
        err_view.setReadOnly(True)
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.Monospace)
        err_view.setFont(mono)
        if err_text:
            # Truncate very long ERR files so the dialog stays usable.
            if len(err_text) > 20_000:
                err_view.setPlainText(
                    err_text[:20_000]
                    + "\n\n[... truncated; open the file to see the rest ...]")
            else:
                err_view.setPlainText(err_text)
        else:
            err_view.setPlainText("(No ERR file produced.)")
        layout.addWidget(err_view, 1)

        sep = QFrame(); sep.setFrameShape(QFrame.HLine)
        layout.addWidget(sep)

        # ---- Buttons ----
        button_row = QHBoxLayout()
        # v5.1.1 item 33: in-dialog Load Results button so the user
        # always has a one-click path to import what MYSTRAN just
        # produced, regardless of the auto-load Preference.
        load_btn = QPushButton("Load Results")
        load_btn.setToolTip(
            "Parse the OP2 / F06 from this run and feed the results "
            "into the Result Browser. Same as auto-load when that "
            "Preference is on.")
        # Disable when there's nothing to load.
        can_load = bool(
            (run.op2_path and Path(run.op2_path).exists())
            or (run.f06_path and Path(run.f06_path).exists())
        )
        load_btn.setEnabled(can_load)
        load_btn.clicked.connect(
            lambda: self._load_results_from_run_summary(run, dlg))
        button_row.addWidget(load_btn)
        if run_dir is not None:
            open_btn = QPushButton("Open Run Folder")
            open_btn.clicked.connect(lambda: self._open_folder(run_dir))
            button_row.addWidget(open_btn)
            open_err_btn = QPushButton("Open .ERR")
            open_err_btn.setEnabled(bool(err_path and err_path.exists()))
            open_err_btn.clicked.connect(lambda: self._open_folder(err_path))
            button_row.addWidget(open_err_btn)
            open_f06_btn = QPushButton("Open .F06")
            open_f06_btn.setEnabled(
                bool(run.f06_path and Path(run.f06_path).exists()))
            open_f06_btn.clicked.connect(
                lambda: self._open_folder(run.f06_path))
            button_row.addWidget(open_f06_btn)
        button_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        button_row.addWidget(close_btn)
        layout.addLayout(button_row)
        dlg.exec()

    def _load_results_from_run_summary(self, run, dlg):
        """v5.1.1 item 33: handler for the Load Results button on the
        post-run summary dialog. Calls the same MYSTRAN result adapter
        that auto-load uses, then closes the summary on success.
        """
        from node_runner.solve.mystran_results import load_mystran_results
        bundle = load_mystran_results(run)
        if bundle and self._bundle_has_data(bundle):
            self._consume_mystran_bundle(bundle, run)
            self._update_status(
                f"Loaded {self._bundle_disp_count(bundle):,} "
                f"displacements ({(run.results_source or '?').upper()} "
                f"source).")
            try:
                dlg.accept()
            except Exception:
                pass
        else:
            QMessageBox.warning(
                self, "Load Results",
                "No displacement data could be parsed from the OP2/F06 "
                "in this run folder. Inspect the .ERR file for the "
                "MYSTRAN-reported failure.")

    def _open_folder(self, path):
        """Open a folder (or file's parent) in the host file manager."""
        import subprocess, sys
        p = Path(path)
        target = str(p if p.is_dir() else p.parent)
        try:
            if sys.platform == 'win32':
                if p.is_file():
                    # Open the file in its default app.
                    os.startfile(str(p))
                else:
                    os.startfile(target)
            elif sys.platform == 'darwin':
                subprocess.Popen(['open', target])
            else:
                subprocess.Popen(['xdg-open', target])
        except Exception:
            pass

    def _on_mystran_failed(self, msg):
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        # Build a MystranRun-like stub from any partial info we have.
        run = getattr(self._solver_worker, '_run', None)
        if run is not None:
            self._mystran_results_loaded_msg = (
                "MYSTRAN reported a non-zero exit code.")
            self._show_mystran_summary(run, None)
        else:
            QMessageBox.critical(self, "MYSTRAN run failed", msg)
        self._update_status("MYSTRAN run failed.", is_error=True)

    def _on_mystran_cancelled(self):
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        self._update_status("MYSTRAN run cancelled.", is_error=False)

    def _consume_mystran_bundle(self, bundle, run):
        """Push a MYSTRAN result bundle into the existing OP2 plumbing.

        Mirrors the data shape of ``model.load_op2_results`` so all the
        existing Color-By-Results / Animation / Vector-Overlay paths
        work unchanged.
        """
        self.op2_results = bundle
        # v5.2.1 item 49: snapshot the undeformed bbox diagonal now so
        # ``_sync_deform_scale_to_widget`` has a stable reference even
        # if downstream rendering ever leaks a deformed grid back into
        # ``current_grid`` (it shouldn't anymore -- see Item 49 fix --
        # but the snapshot makes the %-of-model math reproducible).
        try:
            if self.current_grid is not None and self.current_grid.n_points > 0:
                b = self.current_grid.bounds
                self._undeformed_bbox_diag = (
                    (b[1] - b[0]) ** 2
                    + (b[3] - b[2]) ** 2
                    + (b[5] - b[4]) ** 2) ** 0.5
            else:
                self._undeformed_bbox_diag = 0.0
        except Exception:
            self._undeformed_bbox_diag = 0.0
        self._update_results_status(run, bundle)
        # Make the legacy results sidebar widget visible + populate its
        # subcase / type / component combos. Without this the v3.0.0
        # color-by-results path doesn't know which scalar to pull.
        # v5.2.0: the legacy results_widget stays HIDDEN; it acts as a
        # state holder only. The Results sidebar tab is the visible UI.
        try:
            self._populate_results_combos()
        except Exception:
            pass
        # v5.1.1 item 31: feed the legacy ResultBrowserDock's Nodes /
        # Elements tables. v5.0.0 + v5.1.0 set ``self.op2_results`` and
        # called ``_on_results_changed`` but never invoked the dock
        # populator, so the dock auto-opened empty. Calling
        # ``_populate_result_browser`` here mirrors the
        # ``_load_results_dialog`` path so MYSTRAN runs feed the dock
        # the same way a manual File -> Load Results would.
        try:
            self._populate_result_browser()
        except Exception:
            pass
        # v5.1.2 item 36: populate the Results sidebar tab tree
        # unconditionally, regardless of color_mode. Without this the
        # tree stays blank until the user manually flips Color By ->
        # Results, which fragments the "results loaded" UX.
        try:
            if hasattr(self, 'results_tab'):
                self.results_tab.populate_from_results(
                    bundle, file_label=self._results_tab_file_label())
        except Exception:
            pass
        # Focus the Results sidebar tab now that it has content.
        try:
            if hasattr(self, 'sidebar_tabs') and hasattr(self, 'results_tab'):
                self.sidebar_tabs.setCurrentWidget(self.results_tab)
        except Exception:
            pass
        # Trigger the same downstream that the OP2 menu item uses.
        try:
            self._on_results_changed()
        except Exception:
            pass
        # v5.1.2 item 37: the right-side Result Browser / Animation /
        # Vector docks were eliminated -- all three panels now live in
        # the Results sidebar tab as sub-tabs and are always visible
        # there. The legacy ``mystran/auto_open_browser`` QSettings key
        # is retained for back-compat but is no longer read here.

    def _update_results_status(self, run, bundle):
        """Refresh the status-bar Results segment."""
        if not hasattr(self, '_status_results_lbl'):
            return
        n_subcases = len(bundle.get('subcases', {}))
        n_disp = sum(
            len(sc.get('displacements', {}))
            for sc in bundle.get('subcases', {}).values())
        src = (run.results_source or "?").upper()
        sol_label = {101: "SOL 1", 103: "SOL 3", 105: "SOL 5"}.get(
            run.sol, f"SOL {run.sol}")
        self._status_results_lbl.setText(
            f"Results: {sol_label} ({src}) | "
            f"{n_subcases} subcase{'s' if n_subcases != 1 else ''} | "
            f"{n_disp:,} disp")

    # ─── v5.0.0 items 17/18: MYSTRAN run flow ─────────────────────────
    def _configure_mystran(self):
        """Analysis -> Configure MYSTRAN…: open Preferences on the MYSTRAN tab."""
        from node_runner.dialogs.preferences import (
            PreferencesDialog, save_preferences, load_preferences)
        current = load_preferences()
        dlg = PreferencesDialog(current, self)
        # Try to switch to the MYSTRAN tab.
        try:
            for i in range(dlg._tabs.count()):
                if dlg._tabs.tabText(i) == "MYSTRAN":
                    dlg._tabs.setCurrentIndex(i)
                    break
        except Exception:
            pass
        if dlg.exec() == QDialog.Accepted:
            save_preferences(dlg.result_payload())
            self._update_status("MYSTRAN preferences saved.")

    def _open_analysis_history(self):
        """Analysis -> Analysis History…"""
        from node_runner.dialogs.solver_run import AnalysisHistoryDialog
        from node_runner.solve import mystran_settings
        dlg = AnalysisHistoryDialog(mystran_settings.get_scratch_root(),
                                    parent=self)
        dlg.reload_requested.connect(self._on_mystran_history_reload)
        dlg.exec()

    def _on_mystran_history_reload(self, run):
        """User picked a prior run from the history dialog. Re-load
        its results into the current model.
        """
        try:
            from node_runner.solve.mystran_results import load_mystran_results
            bundle = load_mystran_results(run)
            if bundle is None:
                QMessageBox.warning(
                    self, "Analysis History",
                    "Couldn't load OP2 or F06 from that run folder.")
                return
            self._consume_mystran_bundle(bundle, run)
            # v5.0.0 item 20: update the persisted results_source so a
            # prior run that previously read "null" now records the
            # source that actually loaded (OP2 vs F06).
            try:
                if run.bdf_path:
                    run.save_meta(Path(run.bdf_path).parent)
            except Exception:
                pass
        except Exception as exc:
            QMessageBox.critical(
                self, "Analysis History",
                f"Failed to reload that run: {exc}")

    def _run_mystran(self):
        """Analysis -> Run Analysis (MYSTRAN)…: pre-flight -> run-options
        -> solver progress -> auto-load results.

        v5.1.0 item 26: requires an active AnalysisSet. The set drives
        what gets exported (scope/group/load/SPC filtering via
        scope.py) and what shows up in the deck's case-control block.
        common convention -- you can't just "run a model"; you run an
        analysis set.
        """
        if not self.current_generator or not self.current_generator.model:
            QMessageBox.warning(
                self, "Run MYSTRAN",
                "Open or create a model first.")
            return

        # v5.1.0 item 26: block-with-pointer when no AnalysisSet is
        # active. The user goes to AnalysisSet Manager to create or
        # pick one -- we no longer silently run the full model.
        active_set = (self.analysis_sets or {}).get(
            self.active_analysis_set_id)
        if active_set is None:
            ret = QMessageBox.question(
                self, "Run requires an active AnalysisSet",
                "MYSTRAN runs need an active AnalysisSet so the scope "
                "(group target, load/SPC SIDs, solver target) is "
                "explicit.\n\nOpen the AnalysisSet Manager now?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ret == QMessageBox.Yes:
                self._open_analysis_set_manager()
                active_set = (self.analysis_sets or {}).get(
                    self.active_analysis_set_id)
            if active_set is None:
                self._update_status(
                    "MYSTRAN run cancelled -- no active AnalysisSet.")
                return
        # Confirm the AnalysisSet's solver target is MYSTRAN.
        target = getattr(active_set, 'solver_target', 'MYSTRAN')
        if target != 'MYSTRAN':
            ret = QMessageBox.question(
                self, "AnalysisSet target mismatch",
                f"Active AnalysisSet '{active_set.name}' has "
                f"solver_target='{target}', not MYSTRAN. Run anyway?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if ret != QMessageBox.Yes:
                return

        # Verify MYSTRAN binary is configured.
        from node_runner.solve import mystran_settings
        from node_runner.solve.mystran_runner import (
            discover_mystran_executable)
        exe = discover_mystran_executable()
        if not exe:
            ret = QMessageBox.question(
                self, "MYSTRAN not configured",
                "I couldn't find the MYSTRAN executable. Configure it now?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ret == QMessageBox.Yes:
                self._configure_mystran()
                exe = discover_mystran_executable()
            if not exe:
                return

        # Open Run-options dialog.
        from node_runner.dialogs.solver_run import (
            RunMystranDialog, MystranPreflightDialog, SolverProgressDialog)
        run_dlg = RunMystranDialog(self, parent=self)
        if run_dlg.exec() != QDialog.Accepted:
            return
        opts = run_dlg.get_options()

        # v5.1.0 item 26: build the scoped model for pre-flight + export.
        from node_runner.scope import scope_model_to_analysis_set
        try:
            scoped_model, scope_report = scope_model_to_analysis_set(
                self.current_generator.model, active_set,
                groups=self.groups)
        except Exception as exc:
            QMessageBox.critical(
                self, "Scoping failed",
                f"Could not scope model to AnalysisSet '{active_set.name}': "
                f"{exc}")
            return

        # Pre-flight (on the scoped model so users see only relevant
        # blocking issues).
        from node_runner.solve.mystran_preflight import scan_for_mystran
        report = scan_for_mystran(scoped_model, sol=opts.get('sol'))
        pf_dlg = MystranPreflightDialog(report, parent=self)
        if pf_dlg.exec() != QDialog.Accepted:
            return

        # Export deck to working dir. Run-folder name carries the
        # AnalysisSet name so Analysis History distinguishes sets.
        safe_setname = ''.join(c if c.isalnum() or c in '-_' else '_'
                                for c in active_set.name)[:32]
        working_dir = Path(opts['working_dir'])
        # Replace the trailing timestamp-only folder with one that
        # also carries the set name.
        if safe_setname and not working_dir.name.endswith(safe_setname):
            working_dir = working_dir.with_name(
                f"{working_dir.name}_{safe_setname}")
        working_dir.mkdir(parents=True, exist_ok=True)
        bdf_stem = Path(self._loaded_filepath).stem if self._loaded_filepath \
            else "untitled"
        bdf_path = working_dir / f"{bdf_stem}.bdf"
        try:
            self.current_generator._write_bdf(
                str(bdf_path),
                field_format=getattr(
                    active_set, 'field_format', 'short'),
                target='mystran',
                analysis_set=active_set,
                groups=self.groups,
                sollib=mystran_settings.get('default_sollib'),
                quad4typ=mystran_settings.get('default_quad4typ'),
                wtmass=mystran_settings.get('default_wtmass'),
            )
        except Exception as exc:
            QMessageBox.critical(
                self, "Export failed",
                f"Could not write MYSTRAN BDF: {exc}")
            return

        # Status line tells the user exactly what got scoped in.
        self._update_status(
            f"Running '{active_set.name}': "
            f"{scope_report.n_elements_kept:,} elements, "
            f"{scope_report.n_loads_kept} load SID(s), "
            f"{scope_report.n_spcs_kept} SPC SID(s)"
            + (f" (group: {scope_report.group_target})"
               if scope_report.group_target else ""))

        # Launch SolverRunWorker.
        from node_runner.solve.mystran_runner import SolverRunWorker
        self._solver_worker = SolverRunWorker(self)
        self._solver_dialog = SolverProgressDialog(bdf_path.name, parent=self)
        self._solver_worker.progress.connect(
            self._solver_dialog.update_progress)
        self._solver_worker.finished.connect(self._on_mystran_finished)
        self._solver_worker.failed.connect(self._on_mystran_failed)
        self._solver_worker.cancelled.connect(self._on_mystran_cancelled)
        self._solver_dialog.cancel_requested.connect(
            self._solver_worker.cancel)
        self._solver_dialog.start_timer()
        self._solver_dialog.show()
        self._solver_worker.start(
            bdf_path=bdf_path,
            exe_path=exe,
            output_dir=working_dir,
            sol=opts.get('sol', 101),
            analysis_set_name=active_set.name,
        )

    def _on_mystran_finished(self, run):
        """Solver completed cleanly. Auto-load results per preference,
        then show a single summary dialog that says where the files
        landed and surfaces any MYSTRAN-level errors from the ERR file.
        """
        from node_runner.solve import mystran_settings
        from node_runner.solve.mystran_results import load_mystran_results
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        self._update_status(
            f"MYSTRAN finished in {run.wall_time:.1f}s "
            f"(exit code {run.return_code}).")

        # Always try to load results -- the user expectation is "click
        # Run, see results". Status of the parse + the run folder are
        # surfaced together in a single post-run dialog below.
        bundle = None
        if mystran_settings.get('auto_load_results', True):
            bundle = load_mystran_results(run)
            if bundle and self._bundle_has_data(bundle):
                self._consume_mystran_bundle(bundle, run)
                self._mystran_results_loaded_msg = (
                    f"Loaded {self._bundle_disp_count(bundle):,} "
                    f"displacements ({run.results_source.upper()} source).")
            elif bundle is not None:
                # Bundle parsed but empty -- treat as failure visually.
                self._mystran_results_loaded_msg = (
                    "Results file parsed but contained NO displacement "
                    "data. MYSTRAN likely aborted with errors -- see "
                    "the ERR file below.")
            else:
                self._mystran_results_loaded_msg = (
                    "No OP2 or F06 results file found in the run folder.")
        else:
            self._mystran_results_loaded_msg = (
                "Auto-load is disabled in Preferences → MYSTRAN.")

        self._show_mystran_summary(run, bundle)

    def _bundle_has_data(self, bundle):
        """True iff the bundle contains at least one non-empty
        displacement / eigenvalue / eigenvector entry."""
        for sc in (bundle.get('subcases') or {}).values():
            if sc.get('displacements'):
                return True
            if sc.get('eigenvalues'):
                return True
            if sc.get('eigenvectors'):
                return True
        return False

    def _bundle_disp_count(self, bundle):
        return sum(len(sc.get('displacements', {}))
                   for sc in (bundle.get('subcases') or {}).values())

    def _show_mystran_summary(self, run, bundle):
        """Single post-run dialog: where the files are, MYSTRAN ERR
        contents (truncated), and an Open Folder button.

        Replaces the v5.0.0-first-cut behaviour where a successful exit
        code but no parsed results left the user staring at an
        unchanged viewport with no idea where to look.
        """
        from PySide6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
            QPlainTextEdit, QFrame, QDialogButtonBox,
        )
        from PySide6.QtGui import QFont

        run_dir = Path(run.bdf_path).parent if run.bdf_path else None
        err_path = (Path(run.bdf_path).with_suffix('.ERR')
                    if run.bdf_path else None)
        err_text = ""
        if err_path and err_path.exists():
            try:
                err_text = err_path.read_text(
                    encoding='utf-8', errors='replace')
            except Exception:
                err_text = ""
        # Count MYSTRAN-level errors in the ERR for the banner.
        n_err = err_text.upper().count("*ERROR")
        n_warn = err_text.upper().count("*WARNING")

        dlg = QDialog(self)
        dlg.setWindowTitle("MYSTRAN run summary")
        dlg.setMinimumSize(760, 480)

        layout = QVBoxLayout(dlg)

        # ---- Banner ----
        has_data = bool(bundle) and self._bundle_has_data(bundle)
        if has_data and n_err == 0:
            banner_text = (
                f"<b>Run successful.</b> {self._mystran_results_loaded_msg}")
            banner_style = (
                "background-color: #a6e3a1; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        elif has_data:
            banner_text = (
                f"<b>Run produced results, with {n_err} MYSTRAN error(s).</b> "
                f"{self._mystran_results_loaded_msg} Inspect the ERR "
                f"file below before trusting these numbers.")
            banner_style = (
                "background-color: #f9e2af; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        else:
            banner_text = (
                f"<b>MYSTRAN exited cleanly but produced no result data.</b><br>"
                f"{n_err} error(s), {n_warn} warning(s) in the ERR file. "
                f"{self._mystran_results_loaded_msg}")
            banner_style = (
                "background-color: #f38ba8; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;")
        banner = QLabel(banner_text)
        banner.setStyleSheet(banner_style)
        banner.setWordWrap(True)
        layout.addWidget(banner)

        # ---- File summary ----
        files_box = QLabel(
            f"<b>Run folder:</b> <code>{run_dir}</code><br>"
            f"<b>BDF:</b> <code>{Path(run.bdf_path).name if run.bdf_path else '-'}</code> | "
            f"<b>F06:</b> <code>{Path(run.f06_path).name if run.f06_path else '(none)'}</code> | "
            f"<b>OP2:</b> <code>{Path(run.op2_path).name if run.op2_path else '(none)'}</code>"
        )
        files_box.setStyleSheet(
            "font-family: Consolas, monospace; font-size: 11px; "
            "padding: 6px 0;")
        files_box.setWordWrap(True)
        layout.addWidget(files_box)

        # ---- ERR contents ----
        layout.addWidget(QLabel("<b>MYSTRAN .ERR file:</b>"))
        err_view = QPlainTextEdit()
        err_view.setReadOnly(True)
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.Monospace)
        err_view.setFont(mono)
        if err_text:
            # Truncate very long ERR files so the dialog stays usable.
            if len(err_text) > 20_000:
                err_view.setPlainText(
                    err_text[:20_000]
                    + "\n\n[... truncated; open the file to see the rest ...]")
            else:
                err_view.setPlainText(err_text)
        else:
            err_view.setPlainText("(No ERR file produced.)")
        layout.addWidget(err_view, 1)

        sep = QFrame(); sep.setFrameShape(QFrame.HLine)
        layout.addWidget(sep)

        # ---- Buttons ----
        button_row = QHBoxLayout()
        # v5.1.1 item 33: in-dialog Load Results button so the user
        # always has a one-click path to import what MYSTRAN just
        # produced, regardless of the auto-load Preference.
        load_btn = QPushButton("Load Results")
        load_btn.setToolTip(
            "Parse the OP2 / F06 from this run and feed the results "
            "into the Result Browser. Same as auto-load when that "
            "Preference is on.")
        # Disable when there's nothing to load.
        can_load = bool(
            (run.op2_path and Path(run.op2_path).exists())
            or (run.f06_path and Path(run.f06_path).exists())
        )
        load_btn.setEnabled(can_load)
        load_btn.clicked.connect(
            lambda: self._load_results_from_run_summary(run, dlg))
        button_row.addWidget(load_btn)
        if run_dir is not None:
            open_btn = QPushButton("Open Run Folder")
            open_btn.clicked.connect(lambda: self._open_folder(run_dir))
            button_row.addWidget(open_btn)
            open_err_btn = QPushButton("Open .ERR")
            open_err_btn.setEnabled(bool(err_path and err_path.exists()))
            open_err_btn.clicked.connect(lambda: self._open_folder(err_path))
            button_row.addWidget(open_err_btn)
            open_f06_btn = QPushButton("Open .F06")
            open_f06_btn.setEnabled(
                bool(run.f06_path and Path(run.f06_path).exists()))
            open_f06_btn.clicked.connect(
                lambda: self._open_folder(run.f06_path))
            button_row.addWidget(open_f06_btn)
        button_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        button_row.addWidget(close_btn)
        layout.addLayout(button_row)
        dlg.exec()

    def _load_results_from_run_summary(self, run, dlg):
        """v5.1.1 item 33: handler for the Load Results button on the
        post-run summary dialog. Calls the same MYSTRAN result adapter
        that auto-load uses, then closes the summary on success.
        """
        from node_runner.solve.mystran_results import load_mystran_results
        bundle = load_mystran_results(run)
        if bundle and self._bundle_has_data(bundle):
            self._consume_mystran_bundle(bundle, run)
            self._update_status(
                f"Loaded {self._bundle_disp_count(bundle):,} "
                f"displacements ({(run.results_source or '?').upper()} "
                f"source).")
            try:
                dlg.accept()
            except Exception:
                pass
        else:
            QMessageBox.warning(
                self, "Load Results",
                "No displacement data could be parsed from the OP2/F06 "
                "in this run folder. Inspect the .ERR file for the "
                "MYSTRAN-reported failure.")

    def _open_folder(self, path):
        """Open a folder (or file's parent) in the host file manager."""
        import subprocess, sys
        p = Path(path)
        target = str(p if p.is_dir() else p.parent)
        try:
            if sys.platform == 'win32':
                if p.is_file():
                    # Open the file in its default app.
                    os.startfile(str(p))
                else:
                    os.startfile(target)
            elif sys.platform == 'darwin':
                subprocess.Popen(['open', target])
            else:
                subprocess.Popen(['xdg-open', target])
        except Exception:
            pass

    def _on_mystran_failed(self, msg):
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        # Build a MystranRun-like stub from any partial info we have.
        run = getattr(self._solver_worker, '_run', None)
        if run is not None:
            self._mystran_results_loaded_msg = (
                "MYSTRAN reported a non-zero exit code.")
            self._show_mystran_summary(run, None)
        else:
            QMessageBox.critical(self, "MYSTRAN run failed", msg)
        self._update_status("MYSTRAN run failed.", is_error=True)

    def _on_mystran_cancelled(self):
        if self._solver_dialog is not None:
            self._solver_dialog.stop_timer()
            self._solver_dialog.close()
        self._update_status("MYSTRAN run cancelled.", is_error=False)

    def _update_results_status(self, run, bundle):
        """Refresh the status-bar Results segment."""
        if not hasattr(self, '_status_results_lbl'):
            return
        n_subcases = len(bundle.get('subcases', {}))
        n_disp = sum(
            len(sc.get('displacements', {}))
            for sc in bundle.get('subcases', {}).values())
        src = (run.results_source or "?").upper()
        sol_label = {101: "SOL 1", 103: "SOL 3", 105: "SOL 5"}.get(
            run.sol, f"SOL {run.sol}")
        self._status_results_lbl.setText(
            f"Results: {sol_label} ({src}) | "
            f"{n_subcases} subcase{'s' if n_subcases != 1 else ''} | "
            f"{n_disp:,} disp")

    def _create_default_analysis_set(self):
        """Create a default analysis set from legacy fields (sol_type, subcases, eigrl_cards)."""
        from node_runner.dialogs.analysis import AnalysisSet
        aset = AnalysisSet(
            id=1,
            name="Default",
            sol_type=self.sol_type,
            subcases=list(self.subcases) if self.subcases else [],
            eigrl=dict(self.eigrl_cards[0]) if self.eigrl_cards else None,
        )
        self.analysis_sets[1] = aset
        self.active_analysis_set_id = 1

    def _set_active_analysis_set(self, set_id):
        """Set the given analysis set as active (from context menu)."""
        if set_id in self.analysis_sets:
            self.active_analysis_set_id = set_id
            aset = self.analysis_sets[set_id]
            # Sync legacy fields
            self.sol_type = aset.sol_type
            self.subcases = aset.subcases
            if aset.eigrl:
                self.eigrl_cards = [aset.eigrl]
            self._update_status(f"Active analysis set: '{aset.name}'")
            self._populate_tree()

    def _delete_analysis_set(self, set_id):
        """Delete an analysis set (from context menu)."""
        if set_id not in self.analysis_sets:
            return
        if len(self.analysis_sets) <= 1:
            QMessageBox.information(self, "Cannot Delete",
                                    "At least one analysis set must remain.")
            return
        reply = QMessageBox.question(
            self, "Confirm Deletion",
            f"Delete analysis set '{self.analysis_sets[set_id].name}'?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            del self.analysis_sets[set_id]
            if self.active_analysis_set_id == set_id:
                self.active_analysis_set_id = next(iter(self.analysis_sets.keys()))
            self._update_status(f"Deleted analysis set {set_id}.")
            self._populate_tree()

    # --- Analysis Tab helpers ---

    def _show_analysis_tab_context_menu(self, pos):
        """Context menu for the Analysis sidebar tab tree."""
        item = self.analysis_tab_tree.itemAt(pos)
        if not item:
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if not data:
            return

        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        kind = data[0]

        if kind == 'analysis_set':
            set_id = data[1]
            aset = (self.analysis_sets or {}).get(set_id)
            target = getattr(aset, 'solver_target', 'MYSTRAN') if aset else 'MYSTRAN'
            # v5.1.0 item 26e: run + export entries at the top of the
            # context menu so the user can fire a set without going
            # through Analysis > Run Analysis (MYSTRAN)...
            run_label = f"Run with {target}"
            run_action = menu.addAction(run_label)
            run_action.setToolTip(
                "Set this analysis set as active and run MYSTRAN.")
            run_action.triggered.connect(
                lambda checked=False, sid=set_id:
                    self._run_analysis_set_from_tab(sid))
            if target != 'MYSTRAN':
                run_action.setEnabled(False)
                run_action.setToolTip(
                    f"v5.1.0 only supports MYSTRAN runs. "
                    f"Change AnalysisSet solver_target to MYSTRAN or "
                    f"use Export... for an MSC/NX deck.")
            export_action = menu.addAction("Export…")
            export_action.setToolTip(
                "Export this AnalysisSet's scoped deck to disk "
                "(solver target / field format are inherited from "
                "the set; override in the Save dialog).")
            export_action.triggered.connect(
                lambda checked=False, sid=set_id:
                    self._export_analysis_set_from_tab(sid))
            menu.addSeparator()
            active_action = menu.addAction("Set as Active")
            active_action.triggered.connect(
                lambda checked=False, sid=set_id:
                    self._set_active_analysis_set(sid))
            edit_action = menu.addAction(f"Edit Analysis Set {set_id}…")
            edit_action.triggered.connect(self._open_analysis_set_manager)
            menu.addSeparator()
            delete_action = menu.addAction(f"Delete Analysis Set {set_id}…")
            delete_action.triggered.connect(
                lambda checked=False, sid=set_id: self._delete_analysis_set(sid))

        if menu.actions():
            menu.exec(self.analysis_tab_tree.viewport().mapToGlobal(pos))

    def _run_analysis_set_from_tab(self, set_id):
        """v5.1.0 item 26e: right-click 'Run with MYSTRAN' on an
        AnalysisSet. Sets active + invokes the existing run flow."""
        if set_id not in (self.analysis_sets or {}):
            return
        self._set_active_analysis_set(set_id)
        self._run_mystran()

    def _export_analysis_set_from_tab(self, set_id):
        """v5.1.0 item 26e: right-click 'Export…' on an AnalysisSet.
        Pre-populates the Save BDF dialog with this set's
        solver_target + field_format and a group-scope-prefilled
        scope picker.
        """
        if set_id not in (self.analysis_sets or {}):
            return
        prior_active = self.active_analysis_set_id
        try:
            self._set_active_analysis_set(set_id)
            # Reuse the canonical Save flow; SaveBdfDialog's Scope
            # combo will list "Active AnalysisSet: <name>" as an option.
            self._save_file_dialog()
        finally:
            # Restore the previously-active set so an Export doesn't
            # silently switch what Run-with-MYSTRAN targets.
            if prior_active != set_id:
                self.active_analysis_set_id = prior_active

    def _new_analysis_set_from_tab(self):
        """Create a new analysis set and open the manager."""
        from node_runner.dialogs.analysis import AnalysisSet
        new_id = max(self.analysis_sets.keys(), default=0) + 1
        aset = AnalysisSet(id=new_id, name=f"Analysis Set {new_id}")
        self.analysis_sets[new_id] = aset
        if self.active_analysis_set_id is None:
            self.active_analysis_set_id = new_id
        self._populate_analysis_tab()
        self._open_analysis_set_manager()

    def _set_active_from_tab(self):
        """Set the currently selected analysis set as active (from tab button)."""
        item = self.analysis_tab_tree.currentItem()
        if not item:
            QMessageBox.information(self, "No Selection",
                                    "Select an analysis set first.")
            return
        data = item.data(0, QtCore.Qt.UserRole)
        if data and data[0] == 'analysis_set':
            self._set_active_analysis_set(data[1])
        elif data and data[0] in ('analysis_subcase', 'analysis_params'):
            # Selected a child - use the parent set's ID
            self._set_active_analysis_set(data[1])
        else:
            QMessageBox.information(self, "No Selection",
                                    "Select an analysis set first.")

    # --- Post-processing (OP2 results) handlers ---

    def _detach_results(self):
        """v5.2.2 issue 5: clear the loaded results bundle while
        keeping the model itself open. Drops back to property coloring
        so the user sees the geometry without the contour."""
        self.op2_results = None
        # Status bar Results segment.
        if hasattr(self, '_status_results_lbl'):
            try:
                self._status_results_lbl.setText("")
            except Exception:
                pass
        # Clear the toolbox dropdowns.
        try:
            if hasattr(self, 'results_tab'):
                self.results_tab.populate_from_results(None)
        except Exception:
            pass
        # v5.3.0 item 58: prefix-sweep so labels-text companion actors
        # also get cleaned up.
        for nm in ('result_min_marker', 'result_max_marker',
                   'result_element_labels', 'result_undeformed_ref'):
            self._remove_actors_by_prefix(nm)
        # Flip the legacy results_widget hidden + switch color mode.
        try:
            self._set_coloring_mode('property')
        except Exception:
            pass
        self._update_plot_visibility()
        try:
            from node_runner.profiling import perf_event
            perf_event('results', 'detached')
        except Exception:
            pass
        self._update_status("Results detached.")

    def _load_results_dialog(self):
        """v5.1.1 item 34: graceful OP2 + F06 loader.

        Accepts both .op2 and .f06 files. On OP2 load failure (typical
        for MYSTRAN OP2s under pyNastran 1.4.x, which lacks a MYSTRAN
        dialect hook), looks for a sibling .F06 next to the .op2 and
        falls through to the minimal F06 parser that the MYSTRAN run
        path uses. Routes the bundle through ``_consume_mystran_bundle``
        so both the new Results sidebar tab AND the legacy Result
        Browser dock get populated identically to the post-run path.
        """
        if not self.current_generator:
            QMessageBox.warning(self, "No Model",
                                "Please load a model first.")
            return
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Results File", "",
            "Results files (*.op2 *.f06 *.F06);;OP2 Files (*.op2);;"
            "F06 Files (*.f06 *.F06);;All Files (*)")
        if not filepath:
            return

        # Build a tiny pseudo-run so _consume_mystran_bundle can do its
        # status-bar update with a sensible label.
        from node_runner.solve import MystranRun
        pseudo_run = MystranRun(
            bdf_path=Path(filepath),
            f06_path=None,
            op2_path=None,
            sol=101,
            analysis_set_name="(manual load)",
        )

        bundle = None
        load_path_label = ""

        if filepath.lower().endswith('.f06'):
            # Direct F06 path -- use the MYSTRAN minimal F06 parser.
            try:
                from node_runner.solve.mystran_results import (
                    _load_f06_minimal)
                bundle = _load_f06_minimal(filepath)
                pseudo_run.f06_path = Path(filepath)
                pseudo_run.results_source = "f06"
                load_path_label = "F06"
            except Exception as e:
                QMessageBox.critical(
                    self, "Error",
                    f"Failed to parse F06: {e}")
                return
        else:
            # OP2 first via pyNastran's load_op2_results.
            from node_runner.model import load_op2_results
            try:
                bundle = load_op2_results(filepath)
                pseudo_run.op2_path = Path(filepath)
                pseudo_run.results_source = "op2"
                load_path_label = "OP2"
            except Exception as exc:
                # MYSTRAN OP2 fallback: pyNastran 1.4.x raises a
                # NoneType / write attribute error on MYSTRAN-dialect
                # OP2s. Look for a sibling .F06 and route to the
                # F06 parser.
                try:
                    from node_runner.profiling import perf_event
                    perf_event('results_load', 'op2_failed',
                               err=str(exc)[:160])
                except Exception:
                    pass
                op2_path = Path(filepath)
                f06_guess = None
                for ext in (".F06", ".f06"):
                    cand = op2_path.with_suffix(ext)
                    if cand.exists():
                        f06_guess = cand
                        break
                if f06_guess is None:
                    # Last-resort: any .f06 file in the same dir.
                    try:
                        siblings = sorted(op2_path.parent.glob("*.f06")) + \
                                   sorted(op2_path.parent.glob("*.F06"))
                        if siblings:
                            f06_guess = siblings[0]
                    except Exception:
                        f06_guess = None
                if f06_guess is None:
                    QMessageBox.critical(
                        self, "OP2 read failed",
                        "pyNastran could not read this OP2 file. If "
                        "it was produced by MYSTRAN, point Node Runner "
                        "at the matching .F06 (in the same folder) "
                        "instead -- the minimal F06 parser handles "
                        "MYSTRAN's native displacement format.\n\n"
                        f"Underlying error:\n{exc}")
                    return
                # Fall through to F06.
                try:
                    from node_runner.solve.mystran_results import (
                        _load_f06_minimal)
                    bundle = _load_f06_minimal(str(f06_guess))
                    pseudo_run.op2_path = Path(filepath)
                    pseudo_run.f06_path = f06_guess
                    pseudo_run.results_source = "f06"
                    load_path_label = (
                        f"F06 fallback (pyNastran rejected the OP2; "
                        f"used sibling {f06_guess.name})")
                except Exception as fexc:
                    QMessageBox.critical(
                        self, "Both OP2 and F06 fallback failed",
                        f"OP2: {exc}\n\nF06 fallback "
                        f"({f06_guess.name}): {fexc}")
                    return

        if not bundle or not bundle.get('subcases'):
            QMessageBox.warning(
                self, "No data",
                "The file loaded but contained no displacement / "
                "eigenvalue tables Node Runner can show.")
            return

        # Funnel through the same consumer as post-run auto-load so
        # both the new Results tab AND the legacy Result Browser dock
        # land populated.
        self._consume_mystran_bundle(bundle, pseudo_run)
        subcases = sorted(bundle['subcases'].keys())
        self._update_status(
            f"Loaded results: {os.path.basename(filepath)} "
            f"({len(subcases)} subcase(s), {load_path_label}).")

    def _populate_results_combos(self):
        """Populate results sidebar combos from loaded OP2 data."""
        if not self.op2_results:
            return
        self.results_subcase_combo.clear()
        for sc_id in sorted(self.op2_results['subcases'].keys()):
            self.results_subcase_combo.addItem(f"Subcase {sc_id}", sc_id)
        self._on_results_type_changed()

    def _populate_mode_combo(self):
        """Populate mode selector from eigenvector data of current subcase."""
        self.results_mode_combo.blockSignals(True)
        self.results_mode_combo.clear()
        sc_id = self.results_subcase_combo.currentData()
        if sc_id is not None and self.op2_results:
            sc_data = self.op2_results['subcases'].get(sc_id, {})
            eig_list = sc_data.get('eigenvectors', [])
            freqs = sc_data.get('frequencies', [])
            for i in range(len(eig_list)):
                freq = freqs[i] if i < len(freqs) and freqs[i] is not None else None
                label = f"Mode {i + 1}"
                if freq is not None:
                    label += f" ({freq:.4g} Hz)"
                self.results_mode_combo.addItem(label, i)
        self.results_mode_combo.blockSignals(False)

    def _on_results_type_changed(self):
        """Update component combo based on selected result type."""
        result_type = self.results_type_combo.currentText()
        self.results_component_combo.clear()
        if result_type in ("Displacement", "Eigenvector", "SPC Forces"):
            self.results_component_combo.addItems([
                "Magnitude", "T1 (X)", "T2 (Y)", "T3 (Z)",
                "R1 (RX)", "R2 (RY)", "R3 (RZ)"
            ])
        elif result_type == "Stress":
            self.results_component_combo.addItems([
                "von Mises", "Max Principal", "Min Principal",
                "XX", "YY", "XY"
            ])
        is_eig = (result_type == "Eigenvector")
        self.results_mode_combo.setVisible(is_eig)
        if is_eig:
            self._populate_mode_combo()
        self._on_results_changed()

    # ─── v5.1.0 item 25: Results sidebar tab handlers ──────────────────
    #
    # The ResultsTab widget emits high-level requests (contour, vector,
    # deform, animate, ...) plus control-panel state changes (color
    # range, levels, min/max markers, element labels, deformation).
    # Each handler below translates the request into operations on
    # the existing v5.0.x results infrastructure -- in particular the
    # display-side results_subcase_combo / results_type_combo /
    # results_component_combo, which drive _update_plot_visibility's
    # results branch via _on_results_changed.

    def _results_tab_select_combos(self, subcase, kind, component):
        """Drive the existing results combos from a ResultsTab request.

        ``kind`` is one of 'displacement', 'stress', 'spc_forces',
        'eigenvector'. ``component`` is a string like 'T3' /
        'von_mises' / 'Mode 2'.
        """
        # Map kind to the existing results_type_combo entries.
        type_map = {
            'displacement': 'Displacement',
            'stress':       'Stress',
            'spc_forces':   'SPC Forces',
            'eigenvector':  'Eigenvector',
        }
        type_label = type_map.get(kind)
        if type_label is None:
            return
        idx = self.results_type_combo.findText(type_label)
        if idx >= 0:
            self.results_type_combo.setCurrentIndex(idx)
        # Subcase combo holds (subcase_id) as currentData.
        for i in range(self.results_subcase_combo.count()):
            if self.results_subcase_combo.itemData(i) == int(subcase):
                self.results_subcase_combo.setCurrentIndex(i)
                break
        # For eigenvectors the component label encodes the mode index;
        # the existing UI shows a separate mode_combo.
        if kind == 'eigenvector' and component.startswith('Mode '):
            try:
                mode_idx = int(component.split()[1]) - 1
                for i in range(self.results_mode_combo.count()):
                    if self.results_mode_combo.itemData(i) == mode_idx:
                        self.results_mode_combo.setCurrentIndex(i)
                        break
            except Exception:
                pass
        else:
            idx = self.results_component_combo.findText(component)
            if idx >= 0:
                self.results_component_combo.setCurrentIndex(idx)

    # ------------------------------------------------------------------
    # v5.2.0 item 42: Post-Processing Toolbox signal handlers
    # ------------------------------------------------------------------

    def _on_results_output_set_changed(self, sid):
        """User picked a different Output Set in the toolbox.

        v5.3.0 item 56: no longer auto-engages Results color mode.
        Results load with the model undeformed + uncoloured
        until the user explicitly picks a Deform Style or Contour
        Style. Picking an Output Set only stages the selection.
        """
        sid = int(sid)
        # Drive the legacy results_subcase_combo so the existing
        # rendering branch (mainwindow.py:11879+) keeps working.
        for i in range(self.results_subcase_combo.count()):
            if self.results_subcase_combo.itemData(i) == sid:
                self.results_subcase_combo.setCurrentIndex(i)
                break
        # Only re-render if we're already in Results mode (i.e. user
        # already opted in). Don't auto-engage.
        if self.color_mode == 'results':
            self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'output_set_changed', subcase=sid)
        except Exception:
            pass

    def _on_results_output_vector_changed(self, kind, component, mode_idx):
        """v5.3.2 item 67: single Output Vector dropdown handler.

        Stashes (kind, component, mode_idx) on self so both the
        contour render path and the deformation block read from it.
        Also drives the legacy results_type / results_component /
        results_mode combos so the existing render pipeline (which
        reads from those) gets the right scalar.
        """
        self._output_vector_kind = kind
        self._output_vector_component = component
        self._output_vector_mode_idx = int(mode_idx)
        type_map = {
            'displacement': 'Displacement',
            'stress':       'Stress',
            'spc_forces':   'SPC Forces',
            'eigenvector':  'Eigenvector',
        }
        type_label = type_map.get(kind, 'Displacement')
        idx = self.results_type_combo.findText(type_label)
        if idx >= 0:
            self.results_type_combo.setCurrentIndex(idx)
        legacy_disp_label = {
            'Magnitude':  'Magnitude', 'T1': 'T1 (X)', 'T2': 'T2 (Y)',
            'T3': 'T3 (Z)', 'R1': 'R1 (RX)', 'R2': 'R2 (RY)',
            'R3': 'R3 (RZ)', 'RMAG': 'Magnitude',
        }
        legacy_stress_label = {
            'von_mises': 'von Mises',
            'max_principal': 'Max Principal',
            'min_principal': 'Min Principal',
            'oxx': 'XX', 'oyy': 'YY', 'txy': 'XY',
        }
        if kind == 'stress':
            comp_label = legacy_stress_label.get(component, component)
        else:
            comp_label = legacy_disp_label.get(component, component)
        idx = self.results_component_combo.findText(comp_label)
        if idx >= 0:
            self.results_component_combo.setCurrentIndex(idx)
        if kind == 'eigenvector' and mode_idx >= 0:
            for i in range(self.results_mode_combo.count()):
                if self.results_mode_combo.itemData(i) == int(mode_idx):
                    self.results_mode_combo.setCurrentIndex(i)
                    break
        # Re-render only if we're already in Results mode (user
        # opted in via a Style toggle). Picking an Output Vector
        # never auto-engages results coloring.
        if self.color_mode == 'results':
            self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'output_vector_changed',
                       kind=kind, component=component, mode_idx=mode_idx)
        except Exception:
            pass

    def _on_results_deform_style_changed(self, style):
        """User picked Undeformed / Deformed / Animate."""
        self._results_deform_style = style
        if style == 'undeformed':
            self.deformation_scale_input.setText("0.0")
            # Stop any running animation.
            if self.anim_play_btn.isChecked():
                self.anim_play_btn.setChecked(False)
        elif style == 'deformed':
            self._sync_deform_scale_to_widget()
            if self.anim_play_btn.isChecked():
                self.anim_play_btn.setChecked(False)
            self._ensure_results_mode_active(trigger='deform')
        elif style == 'animate':
            self._sync_deform_scale_to_widget()
            self._ensure_results_mode_active(trigger='animate')
            if not self.anim_play_btn.isChecked():
                self.anim_play_btn.setChecked(True)
        self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'deform_style', style=style)
        except Exception:
            pass

    def _on_results_deform_scale_changed(self, mode, value):
        """User changed the deform scale mode (% of model | actual)."""
        self._results_scale_mode = mode      # 'pct' or 'actual'
        self._results_scale_value = float(value)
        self._sync_deform_scale_to_widget()
        self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'scale_mode', mode=mode, value=value)
        except Exception:
            pass

    def _sync_deform_scale_to_widget(self):
        """Translate the toolbox's scale mode/value into the legacy
        ``deformation_scale_input`` value that ``_update_plot_visibility``
        reads.

        - ``actual``: pass the spinbox value straight through.
        - ``pct``: scale so that the largest displacement reaches
          ``value / 100 * model_bbox_diagonal``. Falls back to the raw
          value if no displacement data is available yet.
        """
        mode = getattr(self, '_results_scale_mode', 'pct')
        value = float(getattr(self, '_results_scale_value', 10.0))
        if mode == 'actual':
            self.deformation_scale_input.setText(f"{value:.6g}")
            return
        # Percent-of-model mode. Need bbox diagonal + max disp magnitude.
        scale = value   # safe fallback
        try:
            grid = getattr(self, 'current_grid', None)
            bundle = getattr(self, 'op2_results', None)
            # v5.2.0 bug-fix (Round 2): MainWindow initializes
            # ``current_grid`` to an empty UnstructuredGrid with
            # degenerate bounds (1, -1, 1, -1, 1, -1) before any model
            # loads. Guard against that so the % of model math doesn't
            # produce nonsense (e.g. ~0.17 from sqrt(12)*0.05) when no
            # real geometry is present.
            if grid is not None and bundle is not None:
                bnds = grid.bounds
                has_real_bounds = (
                    bnds[1] > bnds[0]
                    and bnds[3] > bnds[2]
                    and bnds[5] > bnds[4])
            else:
                has_real_bounds = False
            if grid is not None and bundle is not None and has_real_bounds:
                # v5.2.1 item 49: prefer the snapshot taken at results
                # load time so subsequent renders never see a diag that
                # has crept upward via an accidental grid mutation.
                snap = float(getattr(self, '_undeformed_bbox_diag', 0.0) or 0.0)
                if snap > 1e-30:
                    diag = snap
                else:
                    diag = (
                        (bnds[1] - bnds[0]) ** 2 +
                        (bnds[3] - bnds[2]) ** 2 +
                        (bnds[5] - bnds[4]) ** 2) ** 0.5
                sid = self.results_subcase_combo.currentData()
                sc = bundle.get('subcases', {}).get(sid, {})
                disps = sc.get('displacements') or {}
                if disps:
                    max_mag = max(
                        (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5
                        for v in disps.values()
                        if v and len(v) >= 3
                    )
                    if max_mag > 1e-30 and diag > 1e-30:
                        scale = (value / 100.0) * diag / max_mag
        except Exception:
            scale = value
        self.deformation_scale_input.setText(f"{scale:.6g}")

    def _on_results_animation_settings_changed(self, mode, n_frames, delay_ms):
        """User changed Sine vs Through Modes, frame count, or delay."""
        self._anim_mode = mode                   # 'sine' or 'modes'
        self._anim_n_frames = int(n_frames)
        self._anim_delay_ms = int(delay_ms)
        # If the timer is running, apply the new interval immediately.
        if self._anim_timer.isActive():
            self._anim_timer.setInterval(max(10, int(delay_ms)))
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'anim_settings',
                       mode=mode, n_frames=n_frames, delay_ms=delay_ms)
        except Exception:
            pass

    def _on_results_show_undeformed_ref_changed(self, on):
        """Toggle the faint undeformed-wireframe overlay."""
        self._results_show_undeformed_ref = bool(on)
        self._refresh_undeformed_ref_actor()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'undeformed_ref', enabled=bool(on))
        except Exception:
            pass

    def _refresh_undeformed_ref_actor(self):
        """Build (or remove) the undeformed wireframe reference actor."""
        actor_name = 'result_undeformed_ref'
        try:
            if actor_name in self.plotter.actors:
                self.plotter.remove_actor(actor_name)
        except Exception:
            pass
        if not getattr(self, '_results_show_undeformed_ref', False):
            try:
                self.plotter.render()
            except Exception:
                pass
            return
        grid = getattr(self, 'current_grid', None)
        if grid is None or grid.n_points == 0:
            return
        # Use the ORIGINAL (un-deformed) points to draw a faint
        # wireframe. v5.2.2 issue 4: dropped the
        # ``hasattr(gen, 'node_id_to_index')`` gate -- some loaded
        # decks don't have that attribute on the generator, which made
        # the wireframe silently skip even when the user toggled it on.
        try:
            import numpy as _np
            # Build a wireframe copy from the source grid. Deep copy
            # so any later mutation of current_grid doesn't leak.
            wire = grid.copy(deep=True)
            self.plotter.add_mesh(
                wire,
                style='wireframe',
                color='#888888',
                opacity=0.35,
                line_width=1,
                name=actor_name,
                    reset_camera=False,
                    show_scalar_bar=False,
                )
            self.plotter.render()
        except Exception:
            pass

    def _on_results_contour_style_changed(self, style):
        """User picked Filled / Filled+Edges / Discrete Bands / No Contour.

        v5.3.2 item 68: explicit, aggressive teardown when style=='none'
        so the contour really goes away:
        1. Sweep every scalar bar from the plotter (the title-keyed
           ``remove_scalar_bar`` is unreliable; also sweep the renderer
           for any leftover vtkScalarBarActor).
        2. Remove the named ``shells`` and ``beams`` actors so the
           next ``add_mesh`` builds fresh ones with property RGB
           (instead of inheriting the old contour mapper).
        3. Remove min/max marker + element-label actors.
        4. Force a re-render via ``_update_plot_visibility`` regardless
           of current color_mode -- if results aren't active we should
           still render the unchanged geometry to flush any stale
           contour pixels.
        """
        self._results_contour_style = style
        if style in ('filled', 'filled_edges', 'bands'):
            self._ensure_results_mode_active(trigger='contour_style')
        elif style == 'none':
            # 1. Scalar bars.
            self._clear_scalar_bars()
            # 2. Force shells/beams actors to be recreated. Without
            #    this PyVista's add_mesh(name='shells', scalars=...)
            #    reuses the existing actor and may keep the old
            #    scalar/cmap binding alive.
            for nm in ('shells', 'beams'):
                try:
                    if nm in self.plotter.actors:
                        self.plotter.remove_actor(nm)
                except Exception:
                    pass
            # 3. Markers + labels (prefix sweep handles companion text
            #    actors PyVista adds for add_point_labels).
            for prefix in ('result_min_marker', 'result_max_marker',
                           'result_element_labels'):
                self._remove_actors_by_prefix(prefix)
        # 4. Re-render. Even if color_mode != 'results' (e.g. user
        #    never engaged any style), we still want the next frame
        #    to reflect the property RGB state.
        try:
            self._update_plot_visibility()
        except Exception:
            pass
        try:
            self.plotter.render()
        except Exception:
            pass
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'contour_style', style=style,
                       color_mode=self.color_mode)
        except Exception:
            pass

    def _on_results_data_conversion_changed(self, conversion):
        """User picked Average / No Averaging / Max at Node / Min at Node."""
        self._results_data_conversion = conversion
        self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'data_conversion', mode=conversion)
        except Exception:
            pass

    def _on_results_palette_changed(self, palette):
        """User picked a colormap (jet/viridis/plasma/coolwarm/gray)."""
        self._results_palette = palette
        self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'palette_changed', palette=palette)
        except Exception:
            pass

    def _on_results_copy_values(self, subcase, kind, component):
        """Legacy 'Copy values to CSV' kept for back-compat with any
        external caller (e.g. tests). v5.2.0 routes the same workflow
        through Tools -> Data Table instead."""
        bundle = getattr(self, 'op2_results', None)
        if not bundle:
            return
        sc = bundle.get('subcases', {}).get(int(subcase))
        if not sc:
            return
        lines = ["id,value"]
        try:
            if kind == 'displacement':
                arr = sc.get('displacements') or {}
                comp_idx = {'T1': 0, 'T2': 1, 'T3': 2,
                            'R1': 3, 'R2': 4, 'R3': 5}.get(component)
                for nid, v in sorted(arr.items()):
                    if comp_idx is not None and comp_idx < len(v):
                        lines.append(f"{nid},{v[comp_idx]:.6e}")
                    elif component == 'Magnitude' and len(v) >= 3:
                        mag = (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5
                        lines.append(f"{nid},{mag:.6e}")
            elif kind == 'stress':
                arr = sc.get('stresses') or {}
                for eid, sdict in sorted(arr.items()):
                    v = sdict.get(component)
                    if v is not None:
                        lines.append(f"{eid},{float(v):.6e}")
        except Exception:
            pass
        from PySide6.QtWidgets import QApplication
        QApplication.clipboard().setText("\n".join(lines))
        self._update_status(
            f"Copied {len(lines) - 1} value(s) to clipboard "
            f"(subcase {subcase} / {kind} / {component}).")

    # --- Control-panel state (v5.1.0 names preserved for the contour pipeline) ---

    def _on_results_color_range_changed(self, auto, vmin, vmax):
        """User toggled Auto OR entered Manual min/max."""
        self._results_user_auto_range = bool(auto)
        self._results_user_vmin = float(vmin)
        self._results_user_vmax = float(vmax)
        # Apply to the active scalar-bar actor if one exists.
        try:
            if not auto and vmax > vmin:
                for actor in self.plotter.actors.values():
                    mapper = getattr(actor, 'mapper', None)
                    if mapper is None:
                        continue
                    try:
                        mapper.SetScalarRange(vmin, vmax)
                    except Exception:
                        continue
                self.plotter.render()
        except Exception:
            pass

    def _on_results_levels_changed(self, n_levels):
        """User changed the contour level count."""
        self._results_n_levels = int(n_levels)
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'levels_changed', n=int(n_levels))
        except Exception:
            pass
        # Trigger a redraw so the new level count takes effect.
        self._on_results_changed()

    def _on_results_markers_changed(self, show_min, show_max):
        """User toggled Min/Max markers."""
        self._results_show_min_marker = bool(show_min)
        self._results_show_max_marker = bool(show_max)
        self._refresh_results_min_max_markers()

    def _on_results_element_labels_changed(self, show, top_n):
        """User toggled Show Element Value Labels."""
        self._results_show_element_labels = bool(show)
        self._results_element_labels_top_n = int(top_n)
        self._refresh_results_element_value_labels()

    # v5.2.0: ``_on_results_deformation_changed`` removed -- the
    # toolbox now drives the deform style via three separate signals
    # (``deform_style_changed`` / ``deform_scale_changed`` /
    # ``animation_settings_changed``). Callers that previously called
    # ``_on_results_deformation_changed(show, scale)`` should route
    # through ``_on_results_deform_style_changed('deformed')`` +
    # ``_on_results_deform_scale_changed('actual', scale)`` instead.

    def _ensure_results_mode_active(self, trigger: str = ''):
        """v5.3.3: switch into Results color mode WITHOUT clobbering the
        user's current Output Vector selection. The previous version
        hardcoded a Displacement-Magnitude default which kept the
        scalar bar stuck on that title regardless of what the user
        picked.

        New behaviour:
        1. If combos haven't been primed yet, prime them.
        2. If the toolbox has a selected Output Vector
           (``_output_vector_*`` attrs), respect it -- drive the legacy
           combos to match.
        3. Otherwise fall back to the first valid scalar in the active
           subcase (Displacement / Magnitude on most decks; Eigenvector
           for SOL 3).
        4. Flip ``color_mode`` to ``'results'`` if not already.
        """
        if not getattr(self, 'op2_results', None):
            return
        bundle = self.op2_results
        subs = bundle.get('subcases', {}) or {}
        if not subs:
            return
        try:
            if self.results_subcase_combo.count() == 0:
                self._populate_results_combos()
        except Exception:
            pass
        try:
            cur_sid = self.results_subcase_combo.currentData()
        except Exception:
            cur_sid = None
        if cur_sid is None:
            first_sid = sorted(subs.keys())[0]
            for i in range(self.results_subcase_combo.count()):
                if self.results_subcase_combo.itemData(i) == first_sid:
                    self.results_subcase_combo.setCurrentIndex(i)
                    break
        try:
            sid_now = self.results_subcase_combo.currentData()
        except Exception:
            sid_now = None
        sc_data = subs.get(sid_now, {}) if sid_now is not None else {}
        # v5.3.3: respect the toolbox's current Output Vector. The
        # _output_vector_kind / _component / _mode_idx attrs are set
        # by _on_results_output_vector_changed; default to displacement
        # / Magnitude only if those attrs aren't meaningful.
        ov_kind = getattr(self, '_output_vector_kind', None)
        ov_component = getattr(self, '_output_vector_component', None)
        ov_mode = int(getattr(self, '_output_vector_mode_idx', -1))
        if not ov_kind or not ov_component:
            # No toolbox selection -- pick a sensible default for this
            # subcase.
            prefer_eig = (
                bool(sc_data.get('eigenvectors')) and
                not sc_data.get('displacements'))
            ov_kind = 'eigenvector' if prefer_eig else 'displacement'
            ov_component = 'Magnitude'
            ov_mode = 0 if prefer_eig else -1
        type_map = {
            'displacement': 'Displacement',
            'stress':       'Stress',
            'spc_forces':   'SPC Forces',
            'eigenvector':  'Eigenvector',
        }
        target_type = type_map.get(ov_kind, 'Displacement')
        try:
            idx = self.results_type_combo.findText(target_type)
            if idx >= 0:
                self.results_type_combo.setCurrentIndex(idx)
        except Exception:
            pass
        legacy_disp_label = {
            'Magnitude': 'Magnitude', 'T1': 'T1 (X)', 'T2': 'T2 (Y)',
            'T3': 'T3 (Z)', 'R1': 'R1 (RX)', 'R2': 'R2 (RY)',
            'R3': 'R3 (RZ)', 'RMAG': 'Magnitude',
        }
        legacy_stress_label = {
            'von_mises': 'von Mises',
            'max_principal': 'Max Principal',
            'min_principal': 'Min Principal',
            'oxx': 'XX', 'oyy': 'YY', 'txy': 'XY',
        }
        if ov_kind == 'stress':
            comp_label = legacy_stress_label.get(ov_component, ov_component)
        else:
            comp_label = legacy_disp_label.get(ov_component, ov_component)
        try:
            comp_idx = self.results_component_combo.findText(comp_label)
            if comp_idx < 0 and ov_kind != 'stress':
                # Fall back to Magnitude/T3 only for displacement-shaped
                # kinds.
                comp_idx = self.results_component_combo.findText("Magnitude")
            if comp_idx < 0 and ov_kind != 'stress':
                comp_idx = self.results_component_combo.findText("T3 (Z)")
            if comp_idx >= 0:
                self.results_component_combo.setCurrentIndex(comp_idx)
        except Exception:
            pass
        if ov_kind == 'eigenvector' and ov_mode >= 0:
            try:
                for i in range(self.results_mode_combo.count()):
                    if self.results_mode_combo.itemData(i) == int(ov_mode):
                        self.results_mode_combo.setCurrentIndex(i)
                        break
            except Exception:
                pass
        # Flip color mode into Results if it isn't already.
        if self.color_mode != 'results':
            try:
                self._set_coloring_mode('results')
            except Exception:
                pass
        try:
            from node_runner.profiling import perf_event
            perf_event(
                'results', 'auto_engage',
                trigger=trigger or 'unspecified',
                subcase=int(sid_now) if sid_now is not None else -1,
                type=self.results_type_combo.currentText(),
                component=self.results_component_combo.currentText())
        except Exception:
            pass

    # v5.2.0: ``_on_results_open_animation_dock`` and
    # ``_on_results_open_vector_dock`` removed -- the sub-tabs they
    # used to focus no longer exist. Vector overlay now lives under
    # View -> Show Displacement Vectors; animation controls live in
    # the Deform section of the Post-Processing Toolbox.

    # --- Actor builders for min/max markers + element value labels ---

    def _active_results_array(self):
        """v5.1.2 item 38: return ``(nids, values, label)`` for the
        currently active subcase/type/component combination, sourced
        from ``self.op2_results`` (the canonical OP2 / F06 bundle)
        rather than scraping VTK grid arrays.

        Returns ``None`` if no bundle, no active selection, or no data
        for the chosen kind/component.

        ``nids`` is an int64 numpy array of node IDs.
        ``values`` is a float numpy array of the same length carrying
        the scalar component for each node.
        ``label`` is a short kind/component string for telemetry.
        """
        import numpy as _np
        bundle = getattr(self, 'op2_results', None)
        if not bundle:
            return None
        # Subcase combo carries the SID via Qt.UserRole (currentData()).
        sc_id = None
        try:
            sc_id = self.results_subcase_combo.currentData()
        except Exception:
            sc_id = None
        if sc_id is None:
            # Fall back to the first subcase in the bundle so the marker
            # is meaningful even before the user touches the combos.
            subs = bundle.get('subcases', {}) or {}
            if not subs:
                return None
            sc_id = sorted(subs.keys())[0]
        sc = bundle.get('subcases', {}).get(int(sc_id))
        if not sc:
            return None
        # Result type
        try:
            result_type = self.results_type_combo.currentText()
        except Exception:
            result_type = "Displacement"
        # Component (display text -> bundle slot)
        try:
            component = self.results_component_combo.currentText()
        except Exception:
            component = "Magnitude"
        comp_idx = {"T1 (X)": 0, "T2 (Y)": 1, "T3 (Z)": 2,
                    "R1 (RX)": 3, "R2 (RY)": 4, "R3 (RZ)": 5,
                    # tolerate the short labels used by the tree
                    "T1": 0, "T2": 1, "T3": 2, "R1": 3, "R2": 4, "R3": 5,
                    }.get(component)
        # Pull the data dict for this kind.
        key_map = {"Displacement": "displacements",
                   "Eigenvector": "eigenvectors",
                   "SPC Forces": "spc_forces"}
        if result_type in key_map:
            raw = sc.get(key_map[result_type], {})
            if result_type == "Eigenvector" and isinstance(raw, list):
                try:
                    mode_idx = max(0, self.results_mode_combo.currentIndex())
                except Exception:
                    mode_idx = 0
                data_dict = raw[mode_idx] if mode_idx < len(raw) else {}
            else:
                data_dict = raw if isinstance(raw, dict) else {}
            if not data_dict:
                return None
            nids = _np.array(sorted(data_dict.keys()), dtype=_np.int64)
            values = _np.zeros(len(nids), dtype=float)
            for i, n in enumerate(nids):
                v = data_dict.get(int(n))
                if not v:
                    continue
                if component == "Magnitude" and len(v) >= 3:
                    values[i] = float((v[0]**2 + v[1]**2 + v[2]**2) ** 0.5)
                elif comp_idx is not None and comp_idx < len(v):
                    values[i] = float(v[comp_idx])
            return nids, values, f"{result_type.lower()}/{component}"
        # Stress is element-centred; signal the caller to handle that
        # via _active_results_element_array() instead.
        return None

    def _active_results_element_array(self):
        """Companion to ``_active_results_array`` for stress / strain.

        Returns ``(eids, values, label)`` for element-centred output
        (currently only Stress) or ``None``.
        """
        import numpy as _np
        bundle = getattr(self, 'op2_results', None)
        if not bundle:
            return None
        try:
            sc_id = self.results_subcase_combo.currentData()
        except Exception:
            sc_id = None
        if sc_id is None:
            subs = bundle.get('subcases', {}) or {}
            if not subs:
                return None
            sc_id = sorted(subs.keys())[0]
        sc = bundle.get('subcases', {}).get(int(sc_id))
        if not sc:
            return None
        try:
            result_type = self.results_type_combo.currentText()
        except Exception:
            result_type = "Stress"
        if result_type != "Stress":
            return None
        try:
            component = self.results_component_combo.currentText()
        except Exception:
            component = "von Mises"
        comp_key_map = {"von Mises": "von_mises",
                        "Max Principal": "max_principal",
                        "Min Principal": "min_principal",
                        "XX": "oxx", "YY": "oyy", "XY": "txy"}
        comp_key = comp_key_map.get(component, "von_mises")
        stresses = sc.get('stresses', {}) or {}
        if not stresses:
            return None
        eids = _np.array(sorted(stresses.keys()), dtype=_np.int64)
        values = _np.array(
            [float(stresses[int(e)].get(comp_key, 0.0)) for e in eids],
            dtype=float)
        return eids, values, f"stress/{comp_key}"

    def _refresh_results_min_max_markers(self):
        """v5.1.2 item 38: build (or remove) the min/max marker
        callouts. Reads scalar values from the canonical OP2 / F06
        bundle (``self.op2_results``) via ``_active_results_array``,
        not from VTK grid arrays. This fixes the v5.1.1 bug where
        VTK's ``vtkOriginalPointIds`` housekeeping array was being
        picked up and reported as if it were a displacement field.
        """
        show_min = getattr(self, '_results_show_min_marker', False)
        show_max = getattr(self, '_results_show_max_marker', False)
        # v5.3.0 item 58: prefix-sweep so the labels-text companion
        # actors get removed too. add_point_labels creates BOTH
        # ``result_min_marker`` AND ``result_min_marker-labels``.
        self._remove_actors_by_prefix('result_min_marker')
        self._remove_actors_by_prefix('result_max_marker')
        if not (show_min or show_max):
            try:
                self.plotter.render()
            except Exception:
                pass
            return
        gen = getattr(self, 'current_generator', None)
        grid = getattr(self, 'current_grid', None)
        if gen is None or grid is None:
            return
        import numpy as _np
        # Try nodal data first (Displacement / Eigenvector / SPC forces).
        active = self._active_results_array()
        if active is not None:
            nids, values, label = active
            if values.size == 0:
                return
            # Map nodes -> grid points via current_node_ids_sorted.
            node_ids_sorted = getattr(self, 'current_node_ids_sorted', None)
            if node_ids_sorted is None:
                return
            id_to_idx = {int(n): i for i, n in enumerate(node_ids_sorted)}
            pts = _np.asarray(grid.points)
            argmin_local = int(_np.argmin(values))
            argmax_local = int(_np.argmax(values))
            self._place_marker_at_node(
                'min', int(nids[argmin_local]), float(values[argmin_local]),
                label, id_to_idx, pts, enabled=show_min)
            self._place_marker_at_node(
                'max', int(nids[argmax_local]), float(values[argmax_local]),
                label, id_to_idx, pts, enabled=show_max)
            try:
                self.plotter.render()
            except Exception:
                pass
            return
        # Element-centred (Stress) fall-through.
        active_e = self._active_results_element_array()
        if active_e is not None:
            eids, values, label = active_e
            if values.size == 0:
                return
            try:
                centers = _np.asarray(grid.cell_centers().points)
            except Exception:
                return
            # Map eid -> grid cell index via grid.cell_data['EID'] if
            # present, else fall back to assuming eids align with cells.
            eid_to_cell = None
            try:
                cell_eids = grid.cell_data.get('EID')
                if cell_eids is not None:
                    cell_eids = _np.asarray(cell_eids).astype(_np.int64)
                    eid_to_cell = {int(e): int(i)
                                   for i, e in enumerate(cell_eids)}
            except Exception:
                eid_to_cell = None
            argmin_local = int(_np.argmin(values))
            argmax_local = int(_np.argmax(values))

            def _place_e(which, eid, value, enabled):
                if not enabled:
                    return
                if eid_to_cell is not None:
                    ci = eid_to_cell.get(int(eid))
                else:
                    ci = int(_np.where(eids == int(eid))[0][0]) \
                        if len(eids) else None
                if ci is None or ci < 0 or ci >= len(centers):
                    return
                p = centers[ci]
                v = float(value)
                self._draw_marker_actor(which, p, v, label)
                try:
                    from node_runner.profiling import perf_event
                    perf_event('results_tab', 'min_max_marker',
                               which=which, scalar=label,
                               value=round(v, 6), eid=int(eid))
                except Exception:
                    pass

            _place_e('min', int(eids[argmin_local]),
                     values[argmin_local], show_min)
            _place_e('max', int(eids[argmax_local]),
                     values[argmax_local], show_max)
            try:
                self.plotter.render()
            except Exception:
                pass

    def _place_marker_at_node(self, which, nid, value, label,
                              id_to_idx, pts, enabled):
        """Place a single min/max marker at a node, or no-op if disabled."""
        if not enabled:
            return
        idx = id_to_idx.get(int(nid))
        if idx is None or idx < 0 or idx >= len(pts):
            return
        self._draw_marker_actor(which, pts[idx], float(value), label)
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'min_max_marker',
                       which=which, scalar=label,
                       value=round(float(value), 6),
                       nid=int(nid))
        except Exception:
            pass

    def _draw_marker_actor(self, which, point, value, label):
        """Place a single point-label callout for the min or max marker."""
        try:
            self.plotter.add_point_labels(
                [point], [f"{which.upper()}: {value:.4g}"],
                name=f"result_{which}_marker",
                point_color='red' if which == 'max' else 'blue',
                point_size=10,
                font_size=12,
                text_color='white',
                shape_color='black',
                shape_opacity=0.6,
                always_visible=True,
                reset_camera=False,
            )
        except Exception:
            pass

    def _build_scalar_bar_kw(self):
        """v5.3.2 item 70: shorter bar + pushed further down so the
        multi-line title above has visible breathing room before the
        top tick value.
        """
        orient = getattr(self, '_results_bar_orientation', 'right')
        base = {
            'color': 'white',
            'title_font_size': 12,
            'label_font_size': 10,
            'fmt': '%.4g',
        }
        if orient == 'bottom':
            base.update({
                'vertical': False,
                'position_x': 0.20,
                'position_y': 0.08,
                'width': 0.60,
                'height': 0.05,
            })
        else:  # 'right'
            # Shorter bar (0.55 height vs 0.65 before) + start lower
            # (0.12 vs 0.18) so the title sits well clear of the top
            # tick. The user reported the title was still too close in
            # v5.3.1 -- this widens the gap by ~50 px on a typical
            # 900-px viewport.
            base.update({
                'vertical': True,
                'position_x': 0.90,
                'position_y': 0.12,
                'width': 0.06,
                'height': 0.55,
            })
        return base

    def _clear_scalar_bars(self):
        """v5.3.1 item 65: sweep every scalar bar from the plotter.

        PyVista's ``add_mesh(scalar_bar_args=...)`` adds a NEW scalar
        bar each time the title differs from existing ones -- it does
        not remove the prior bar. Cycling Contour Vector
        Displacement->Stress accumulates ghost bars. Clear them all
        before each render to keep exactly ONE bar at a time.
        """
        n = 0
        try:
            bars = list(getattr(self.plotter, 'scalar_bars', {}).keys())
            for title in bars:
                try:
                    self.plotter.remove_scalar_bar(title)
                    n += 1
                except Exception:
                    pass
        except Exception:
            pass
        if n > 0:
            try:
                from node_runner.profiling import perf_event
                perf_event('results', 'scalar_bars_cleared', n_removed=n)
            except Exception:
                pass

    def _remove_actors_by_prefix(self, prefix: str):
        """v5.3.0 item 58: prefix-sweep actor removal.

        PyVista's ``add_point_labels(...)`` creates BOTH a points actor
        with the requested ``name`` AND a labels-text actor named
        ``<name>-labels`` (or similar suffix). Removing by the base name
        alone leaves the labels-text actor lingering -- which is why
        v5.2.x's "Show Element Value Labels" checkbox couldn't be
        toggled off. This helper sweeps every actor that matches the
        base name OR begins with ``<base>-`` to wipe both.
        """
        try:
            for nm in list(self.plotter.actors.keys()):
                if nm == prefix or nm.startswith(prefix + '-'):
                    try:
                        self.plotter.remove_actor(nm)
                    except Exception:
                        pass
        except Exception:
            pass

    def _refresh_results_element_value_labels(self):
        """v5.1.2 item 38: build (or remove) the per-element value-label
        callouts. Reads from the canonical OP2 / F06 bundle rather than
        scraping VTK grid arrays. Supports both nodal output (label per
        node) and element output (label per element center).
        """
        show = getattr(self, '_results_show_element_labels', False)
        top_n = int(getattr(self, '_results_element_labels_top_n', 0))
        # v5.3.0 item 58: prefix-sweep so both the points actor AND
        # the labels-text actor get removed.
        self._remove_actors_by_prefix('result_element_labels')
        if not show:
            try:
                self.plotter.render()
            except Exception:
                pass
            return
        grid = getattr(self, 'current_grid', None)
        if grid is None:
            return
        import numpy as _np
        # Prefer element-centred output if active; else fall back to
        # nodal output and label nodes.
        active_e = self._active_results_element_array()
        if active_e is not None:
            eids, values, _label = active_e
            if values.size == 0:
                return
            try:
                centers = _np.asarray(grid.cell_centers().points)
            except Exception:
                return
            try:
                cell_eids = grid.cell_data.get('EID')
                cell_eids = (_np.asarray(cell_eids).astype(_np.int64)
                             if cell_eids is not None else None)
            except Exception:
                cell_eids = None
            if cell_eids is None:
                return
            eid_to_cell = {int(e): int(i)
                           for i, e in enumerate(cell_eids)}
            order = _np.argsort(_np.abs(values))[::-1]
            if top_n > 0:
                order = order[:top_n]
            pick_pts, pick_labels = [], []
            for k in order:
                ci = eid_to_cell.get(int(eids[k]))
                if ci is None or ci < 0 or ci >= len(centers):
                    continue
                pick_pts.append(centers[ci])
                pick_labels.append(f"{float(values[k]):.3g}")
            if not pick_pts:
                return
            self._draw_value_labels(pick_pts, pick_labels)
            try:
                from node_runner.profiling import perf_event
                perf_event('results_tab', 'element_value_labels',
                           n_labels=int(len(pick_pts)), top_n=int(top_n),
                           scope='element')
            except Exception:
                pass
            try:
                self.plotter.render()
            except Exception:
                pass
            return
        active = self._active_results_array()
        if active is None:
            return
        nids, values, _label = active
        if values.size == 0:
            return
        node_ids_sorted = getattr(self, 'current_node_ids_sorted', None)
        if node_ids_sorted is None:
            return
        id_to_idx = {int(n): i for i, n in enumerate(node_ids_sorted)}
        pts = _np.asarray(grid.points)
        order = _np.argsort(_np.abs(values))[::-1]
        if top_n > 0:
            order = order[:top_n]
        pick_pts, pick_labels = [], []
        for k in order:
            pi = id_to_idx.get(int(nids[k]))
            if pi is None or pi < 0 or pi >= len(pts):
                continue
            pick_pts.append(pts[pi])
            pick_labels.append(f"{float(values[k]):.3g}")
        if not pick_pts:
            return
        self._draw_value_labels(pick_pts, pick_labels)
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'element_value_labels',
                       n_labels=int(len(pick_pts)), top_n=int(top_n),
                       scope='node')
        except Exception:
            pass
        try:
            self.plotter.render()
        except Exception:
            pass

    def _draw_value_labels(self, points, labels):
        """Place the value-label actor on the plotter."""
        try:
            self.plotter.add_point_labels(
                points, labels,
                name='result_element_labels',
                font_size=10,
                point_size=4,
                text_color='white',
                shape_opacity=0.0,
                always_visible=False,
                reset_camera=False,
            )
        except Exception:
            pass

    # ─── v5.1.0 item 25: Results sidebar tab handlers ──────────────────
    #
    # The ResultsTab widget emits high-level requests (contour, vector,
    # deform, animate, ...) plus control-panel state changes (color
    # range, levels, min/max markers, element labels, deformation).
    # Each handler below translates the request into operations on
    # the existing v5.0.x results infrastructure -- in particular the
    # display-side results_subcase_combo / results_type_combo /
    # results_component_combo, which drive _update_plot_visibility's
    # results branch via _on_results_changed.

    def _results_tab_select_combos(self, subcase, kind, component):
        """Drive the existing results combos from a ResultsTab request.

        ``kind`` is one of 'displacement', 'stress', 'spc_forces',
        'eigenvector'. ``component`` is a string like 'T3' /
        'von_mises' / 'Mode 2'.
        """
        # Map kind to the existing results_type_combo entries.
        type_map = {
            'displacement': 'Displacement',
            'stress':       'Stress',
            'spc_forces':   'SPC Forces',
            'eigenvector':  'Eigenvector',
        }
        type_label = type_map.get(kind)
        if type_label is None:
            return
        idx = self.results_type_combo.findText(type_label)
        if idx >= 0:
            self.results_type_combo.setCurrentIndex(idx)
        # Subcase combo holds (subcase_id) as currentData.
        for i in range(self.results_subcase_combo.count()):
            if self.results_subcase_combo.itemData(i) == int(subcase):
                self.results_subcase_combo.setCurrentIndex(i)
                break
        # For eigenvectors the component label encodes the mode index;
        # the existing UI shows a separate mode_combo.
        if kind == 'eigenvector' and component.startswith('Mode '):
            try:
                mode_idx = int(component.split()[1]) - 1
                for i in range(self.results_mode_combo.count()):
                    if self.results_mode_combo.itemData(i) == mode_idx:
                        self.results_mode_combo.setCurrentIndex(i)
                        break
            except Exception:
                pass
        else:
            idx = self.results_component_combo.findText(component)
            if idx >= 0:
                self.results_component_combo.setCurrentIndex(idx)

    def _on_results_request_contour(self, subcase, kind, component):
        """User clicked Contour on a Results-tab leaf."""
        if self.color_mode != 'results':
            try:
                self._set_coloring_mode('results')
            except Exception:
                pass
        self._results_tab_select_combos(subcase, kind, component)
        self._on_results_changed()
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'contour_request',
                       subcase=int(subcase), kind=kind, component=component)
        except Exception:
            pass

    def _on_results_request_vector(self, subcase, kind, component):
        """User clicked 'Vector overlay' on a leaf."""
        # Forward to the existing VectorOverlayWidget by toggling its
        # 'disp' / 'reac' checkbox based on the kind.
        try:
            widget = getattr(self, '_vector_overlay_widget', None)
            if widget is not None and hasattr(widget, 'check_disp'):
                if kind in ('displacement', 'eigenvector'):
                    widget.check_disp.setChecked(True)
                if kind == 'spc_forces':
                    widget.check_reac.setChecked(True)
        except Exception:
            pass
        # Ensure the right subcase/component is selected first.
        self._results_tab_select_combos(subcase, kind, component)
        try:
            self._on_vector_overlay_toggled('disp', True)
        except Exception:
            pass

    def _on_results_request_deform(self, subcase, kind, component):
        """User clicked 'Apply as deformed shape'."""
        self._results_tab_select_combos(subcase, kind, component)
        # Set a reasonable default scale if currently zero.
        try:
            if not self.deformation_scale_input.text() or float(
                    self.deformation_scale_input.text()) == 0.0:
                self.deformation_scale_input.setText("1.0")
        except Exception:
            self.deformation_scale_input.setText("1.0")
        self._on_results_changed()

    def _on_results_request_animate(self, subcase, kind, component):
        """User clicked 'Animate'."""
        self._results_tab_select_combos(subcase, kind, component)
        # Start the existing animation timer via anim_play_btn.
        if not self.anim_play_btn.isChecked():
            self.anim_play_btn.setChecked(True)

    def _on_results_copy_values(self, subcase, kind, component):
        """User clicked 'Copy values to CSV'."""
        bundle = getattr(self, 'op2_results', None)
        if not bundle:
            return
        sc = bundle.get('subcases', {}).get(int(subcase))
        if not sc:
            return
        lines = ["id,value"]
        try:
            if kind == 'displacement':
                arr = sc.get('displacements') or {}
                comp_idx = {'T1': 0, 'T2': 1, 'T3': 2,
                            'R1': 3, 'R2': 4, 'R3': 5}.get(component)
                for nid, v in sorted(arr.items()):
                    if comp_idx is not None and comp_idx < len(v):
                        lines.append(f"{nid},{v[comp_idx]:.6e}")
                    elif component == 'TMAG' and len(v) >= 3:
                        mag = (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5
                        lines.append(f"{nid},{mag:.6e}")
                    elif component == 'RMAG' and len(v) >= 6:
                        mag = (v[3]**2 + v[4]**2 + v[5]**2) ** 0.5
                        lines.append(f"{nid},{mag:.6e}")
            elif kind == 'stress':
                arr = sc.get('stresses') or {}
                for eid, sdict in sorted(arr.items()):
                    v = sdict.get(component)
                    if v is not None:
                        lines.append(f"{eid},{float(v):.6e}")
        except Exception:
            pass
        from PySide6.QtWidgets import QApplication
        QApplication.clipboard().setText("\n".join(lines))
        self._update_status(
            f"Copied {len(lines) - 1} value(s) to clipboard "
            f"(subcase {subcase} / {kind} / {component}).")

    # --- Control-panel state ---

    def _on_results_color_range_changed(self, auto, vmin, vmax):
        """User toggled Auto OR entered Manual min/max."""
        self._results_user_auto_range = bool(auto)
        self._results_user_vmin = float(vmin)
        self._results_user_vmax = float(vmax)
        # Apply to the active scalar-bar actor if one exists.
        try:
            if not auto and vmax > vmin:
                for actor in self.plotter.actors.values():
                    mapper = getattr(actor, 'mapper', None)
                    if mapper is None:
                        continue
                    try:
                        mapper.SetScalarRange(vmin, vmax)
                    except Exception:
                        continue
                self.plotter.render()
        except Exception:
            pass

    def _on_results_levels_changed(self, n_levels):
        """User changed the contour level count."""
        self._results_n_levels = int(n_levels)
        try:
            from node_runner.profiling import perf_event
            perf_event('results_tab', 'levels_changed', n=int(n_levels))
        except Exception:
            pass
        # Trigger a redraw so the new level count takes effect.
        self._on_results_changed()

    def _on_results_markers_changed(self, show_min, show_max):
        """User toggled Min/Max markers."""
        self._results_show_min_marker = bool(show_min)
        self._results_show_max_marker = bool(show_max)
        self._refresh_results_min_max_markers()

    def _on_results_element_labels_changed(self, show, top_n):
        """User toggled Show Element Value Labels."""
        self._results_show_element_labels = bool(show)
        self._results_element_labels_top_n = int(top_n)
        self._refresh_results_element_value_labels()

    def _on_results_deformation_changed(self, show, scale):
        """User toggled Show Deformed Shape."""
        if show:
            self.deformation_scale_input.setText(f"{float(scale):.4f}")
        else:
            self.deformation_scale_input.setText("0.0")
        self._on_results_changed()

    def _on_results_open_animation_dock(self):
        """Show the existing Animation dock."""
        dock = getattr(self, '_anim_dock', None)
        if dock is not None:
            dock.show()
            dock.raise_()

    def _on_results_open_vector_dock(self):
        """Show the existing VectorOverlay dock."""
        dock = getattr(self, '_vector_dock', None)
        if dock is not None:
            dock.show()
            dock.raise_()

    # --- Actor builders for min/max markers + element value labels ---

    def _refresh_results_min_max_markers(self):
        """Build (or remove) the min/max marker callouts based on the
        currently active scalar field."""
        show_min = getattr(self, '_results_show_min_marker', False)
        show_max = getattr(self, '_results_show_max_marker', False)
        # Always remove first so toggling off works.
        for name in ('result_min_marker', 'result_max_marker'):
            try:
                if name in self.plotter.actors:
                    self.plotter.remove_actor(name)
            except Exception:
                pass
        if not (show_min or show_max):
            self.plotter.render()
            return
        grid = getattr(self, 'current_grid', None)
        if grid is None:
            return
        # Find the active scalar array on the grid (point_data first,
        # then cell_data).
        scalars = None
        is_point = False
        scalar_name = None
        try:
            for name in grid.point_data.keys():
                if name.lower() in ('eid', 'pid', 'type', 'is_shell',
                                    'cell_rgb', 'cell_visible',
                                    'vtkghosttype'):
                    continue
                scalars = grid.point_data[name]
                scalar_name = name
                is_point = True
                break
        except Exception:
            scalars = None
        if scalars is None:
            try:
                for name in grid.cell_data.keys():
                    if name.lower() in ('eid', 'pid', 'type', 'is_shell',
                                        'cell_rgb', 'cell_visible',
                                        'vtkghosttype'):
                        continue
                    scalars = grid.cell_data[name]
                    scalar_name = name
                    is_point = False
                    break
            except Exception:
                scalars = None
        if scalars is None or len(scalars) == 0:
            return
        import numpy as _np
        arr = _np.asarray(scalars)
        if arr.size == 0:
            return
        argmin = int(_np.argmin(arr))
        argmax = int(_np.argmax(arr))
        if is_point:
            pts = _np.asarray(grid.points)
        else:
            try:
                pts = _np.asarray(grid.cell_centers().points)
            except Exception:
                return
        def _add(name, idx):
            try:
                p = pts[idx]
                v = float(arr[idx])
                self.plotter.add_point_labels(
                    [p], [f"{name.upper()}: {v:.4g}"],
                    name=f"result_{name}_marker",
                    point_color='red' if name == 'max' else 'blue',
                    point_size=10,
                    font_size=12,
                    text_color='white',
                    shape_color='black',
                    shape_opacity=0.6,
                    always_visible=True,
                    reset_camera=False,
                )
                try:
                    from node_runner.profiling import perf_event
                    perf_event('results_tab', 'min_max_marker',
                               which=name, scalar=scalar_name or '',
                               value=round(v, 6))
                except Exception:
                    pass
            except Exception:
                pass
        if show_min:
            _add('min', argmin)
        if show_max:
            _add('max', argmax)
        self.plotter.render()

    def _refresh_results_element_value_labels(self):
        """Build (or remove) the per-element value-label callouts."""
        show = getattr(self, '_results_show_element_labels', False)
        top_n = int(getattr(self, '_results_element_labels_top_n', 0))
        try:
            if 'result_element_labels' in self.plotter.actors:
                self.plotter.remove_actor('result_element_labels')
        except Exception:
            pass
        if not show:
            self.plotter.render()
            return
        grid = getattr(self, 'current_grid', None)
        if grid is None:
            return
        import numpy as _np
        # Use cell_data scalars if any non-housekeeping array exists.
        scalars = None
        try:
            for name in grid.cell_data.keys():
                if name.lower() in ('eid', 'pid', 'type', 'is_shell',
                                    'cell_rgb', 'cell_visible',
                                    'vtkghosttype'):
                    continue
                scalars = _np.asarray(grid.cell_data[name])
                break
        except Exception:
            return
        if scalars is None or scalars.size == 0:
            return
        try:
            centers = _np.asarray(grid.cell_centers().points)
        except Exception:
            return
        if top_n > 0 and top_n < scalars.size:
            order = _np.argsort(_np.abs(scalars))[-top_n:]
            pick = order
        else:
            pick = _np.arange(scalars.size)
        try:
            labels = [f"{float(scalars[i]):.3g}" for i in pick]
            self.plotter.add_point_labels(
                centers[pick], labels,
                name='result_element_labels',
                font_size=10,
                point_size=4,
                text_color='white',
                shape_opacity=0.0,
                always_visible=False,
                reset_camera=False,
            )
            try:
                from node_runner.profiling import perf_event
                perf_event('results_tab', 'element_value_labels',
                           n_labels=int(len(pick)), top_n=int(top_n))
            except Exception:
                pass
        except Exception:
            pass
        self.plotter.render()

    def _on_results_changed(self):
        """Re-render when any results selection changes.

        v5.2.0 Round-2 Bug #4: the v5.1.2 path used to call
        ``results_tab.populate_from_results(...)`` here so the tree-
        widget refreshed on every contour update. That worked when the
        v5.1.2 tab was passive. v5.2.0's toolbox actively drives
        ``output_vector_changed`` on populate, which lands here and
        re-calls populate, producing infinite recursion. Population is
        now done ONCE at load time (``_consume_mystran_bundle`` and
        ``_load_results_dialog``); this method only refreshes the
        rendering side.
        """
        if self.color_mode == "results":
            self._update_plot_visibility()
            try:
                self._refresh_results_min_max_markers()
            except Exception:
                pass
            try:
                self._refresh_results_element_value_labels()
            except Exception:
                pass

    def _results_tab_file_label(self):
        """Compose a human-friendly tree-root label for the Results tab."""
        bundle = getattr(self, 'op2_results', None) or {}
        path_str = bundle.get('filepath', '') or ''
        try:
            stem = Path(path_str).stem if path_str else "Results"
        except Exception:
            stem = "Results"
        src = (bundle.get('results_source') or '').upper()
        if src:
            return f"{stem} ({src} source)"
        return stem or "Results"

    def _toggle_animation(self, checked):
        import math
        if checked:
            # v5.1.2 item 39: auto-engage Results color mode + a default
            # scalar selection so Play actually shows motion instead of
            # silently doing nothing.
            self._ensure_results_mode_active(trigger='animate')
            self.anim_play_btn.setText("Pause")
            self._anim_phase = 0.0
            # v5.2.0 item 44: timer interval comes from the toolbox's
            # delay spinbox (`_anim_delay_ms`) rather than the legacy
            # speed slider. Fall back to the slider value if the new
            # state hasn't been set yet.
            delay = int(getattr(self, '_anim_delay_ms', 0))
            if delay <= 0:
                delay = max(10, 210 - self.anim_speed_slider.value())
            self._anim_timer.start(max(10, delay))
        else:
            self.anim_play_btn.setText("Play")
            self._anim_timer.stop()
            self._anim_override_scale = None
            self._on_results_changed()

    def _on_anim_speed_changed(self, value):
        if self._anim_timer.isActive():
            self._anim_timer.setInterval(max(10, 210 - value))

    def _perf_actor_snapshot(self, label: str, **kw):
        """v5.3.0 item 60: log the current set of plotter actors.

        Used by the animation-tick telemetry to identify whether any
        actor briefly disappears between two ticks (the suspected
        flicker mechanism). Cheap -- only fires when NR_PROFILE is
        active.
        """
        try:
            from node_runner.profiling import perf_event, enabled
            if not enabled():
                return
            names = sorted(self.plotter.actors.keys())
            perf_event('actor_snapshot', label,
                       n=len(names),
                       names=','.join(names)[:300],
                       **kw)
        except Exception:
            pass

    def _on_animation_tick(self):
        """v5.2.0 item 44: dispatch on ``_anim_mode``.

        - ``'sine'``: existing cosine-amplitude path, but stepping with
          a phase increment derived from the configured frame count
          (so changing Frames in the toolbox really changes the period).
        - ``'modes'``: advance one mode index per tick. Wraps at the
          end of the eigenvector list.

        v5.3.0 item 60: wrapped with detailed telemetry to debug the
        persistent animation flicker. Every tick emits:
        * ``actor_snapshot`` events before / after update + render
        * ``animation.tick`` with wall-times for update vs render
        """
        import math
        import time
        t0 = time.perf_counter()
        self._anim_tick_counter = getattr(self, '_anim_tick_counter', 0) + 1
        tick = self._anim_tick_counter
        self._perf_actor_snapshot('anim_tick.before', tick=tick,
                                  phase=round(self._anim_phase, 3))
        mode = getattr(self, '_anim_mode', 'sine')
        n_frames = max(3, int(getattr(self, '_anim_n_frames', 20)))
        if mode == 'modes':
            # Through Modes: cycle the legacy results_mode_combo.
            count = self.results_mode_combo.count()
            if count <= 0:
                # No eigenvector data; fall back to sine.
                mode = 'sine'
            else:
                cur = self.results_mode_combo.currentIndex()
                new = (cur + 1) % count
                self.results_mode_combo.setCurrentIndex(new)
                self._anim_through_modes_current = new
                # Use full deformation scale for each pose (no sine
                # modulation -- we're showing the mode shape itself).
                try:
                    max_scale = float(self.deformation_scale_input.text() or 0)
                except ValueError:
                    max_scale = 1.0
                if max_scale == 0:
                    max_scale = 1.0
                self._anim_override_scale = max_scale
                t_upd_start = time.perf_counter()
                self._update_plot_visibility()
                t_upd_end = time.perf_counter()
                self._perf_actor_snapshot('anim_tick.after_update', tick=tick,
                                          update_ms=round((t_upd_end - t_upd_start) * 1000, 2))
                try:
                    from node_runner.profiling import perf_event
                    perf_event('animation.tick', tick=tick, mode='modes',
                               update_ms=round((t_upd_end - t_upd_start) * 1000, 2),
                               total_ms=round((t_upd_end - t0) * 1000, 2))
                except Exception:
                    pass
                return
        # Default: sine animation.
        phase_step = (2 * math.pi) / float(n_frames)
        self._anim_phase += phase_step
        if self._anim_phase > 2 * math.pi:
            self._anim_phase -= 2 * math.pi
        try:
            max_scale = float(self.deformation_scale_input.text() or 0)
        except ValueError:
            max_scale = 1.0
        if max_scale == 0:
            max_scale = 1.0
        self._anim_override_scale = max_scale * math.sin(self._anim_phase)
        t_upd_start = time.perf_counter()
        self._update_plot_visibility()
        t_upd_end = time.perf_counter()
        self._perf_actor_snapshot('anim_tick.after_update', tick=tick,
                                  update_ms=round((t_upd_end - t_upd_start) * 1000, 2))
        try:
            from node_runner.profiling import perf_event
            perf_event('animation.tick', tick=tick, mode='sine',
                       phase=round(self._anim_phase, 3),
                       override_scale=round(self._anim_override_scale, 6),
                       update_ms=round((t_upd_end - t_upd_start) * 1000, 2),
                       total_ms=round((t_upd_end - t0) * 1000, 2))
        except Exception:
            pass

    def _export_results_gif(self):
        import math
        if not self.op2_results or not self.current_grid:
            QMessageBox.warning(self, "No Results", "Load OP2 results first.")
            return
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Export Animation GIF", "", "GIF Files (*.gif)")
        if not filepath:
            return
        try:
            import imageio
        except ImportError:
            QMessageBox.critical(self, "Missing Package",
                                 "Install imageio: pip install imageio")
            return
        try:
            max_scale = float(self.deformation_scale_input.text() or 1.0)
        except ValueError:
            max_scale = 1.0
        if max_scale == 0:
            max_scale = 1.0
        n_frames = 30
        frames = []
        for i in range(n_frames):
            phase = 2 * math.pi * i / n_frames
            self._anim_override_scale = max_scale * math.sin(phase)
            self._update_plot_visibility()
            self.plotter.render()
            img = self.plotter.screenshot(return_img=True)
            frames.append(img)
        self._anim_override_scale = None
        self._on_results_changed()
        imageio.mimsave(filepath, frames, fps=20, loop=0)
        self._update_status(f"Exported animation GIF: {os.path.basename(filepath)}")

    def _show_free_body_diagram(self):
        if not self.op2_results:
            QMessageBox.warning(self, "No Results",
                                "Load OP2 results first (Results > Load OP2).")
            return
        if not self.current_generator:
            return
        all_nids = list(self.current_generator.model.nodes.keys())
        dlg = FreeBodyDiagramDialog(self.op2_results, all_nids, self)
        dlg.request_render_arrows = self._render_fbd_arrows
        dlg.exec()
        for c in ('green', 'blue', 'red'):
            self.plotter.remove_actor(f'fbd_arrows_{c}', render=False)
        self.plotter.render()

    def _render_fbd_arrows(self, nid, forces, scale, show_applied,
                           show_spc, show_elements):
        """Render colored force arrows at a node for the free body diagram."""
        for c in ('green', 'blue', 'red'):
            self.plotter.remove_actor(f'fbd_arrows_{c}', render=False)
        if not self.current_generator or nid not in self.current_generator.model.nodes:
            return
        node_pos = np.array(
            self.current_generator.model.nodes[nid].get_position())
        color_map = {'APPLIED': 'green', 'SPC': 'blue'}
        all_starts, all_dirs, all_colors = [], [], []
        for f in forces:
            src = f['source']
            if src == 'APPLIED' and not show_applied:
                continue
            if src == 'SPC' and not show_spc:
                continue
            if src.startswith('EID') and not show_elements:
                continue
            fvec = np.array(f['forces'][:3])
            mag = np.linalg.norm(fvec)
            if mag < 1e-20:
                continue
            all_starts.append(node_pos)
            all_dirs.append(fvec * scale)
            all_colors.append(color_map.get(src, 'red'))
        if not all_starts:
            self.plotter.render()
            return
        # Render arrows grouped by color
        for color in ('green', 'blue', 'red'):
            pts = [s for s, c in zip(all_starts, all_colors) if c == color]
            dirs = [d for d, c in zip(all_dirs, all_colors) if c == color]
            if pts:
                import pyvista as pv
                centers = np.array(pts)
                directions = np.array(dirs)
                self.plotter.add_arrows(
                    centers, directions, color=color,
                    name=f'fbd_arrows_{color}', pickable=False,
                    reset_camera=False)
        self.plotter.render()

