"""v5.1.0 item 25: the Femap-killer Results sidebar tab.

Lives as a 5th tab in MainWindow's ``sidebar_tabs`` widget alongside
Model / Loads / Groups / Analysis. The tab has two parts:

1. **Results tree** (top): output sets -> output vectors -> components.
   Right-click a component for: Contour / Vector / Deform / Animate /
   Copy values to CSV.
2. **Contour controls panel** (bottom):
   - Color range Auto/Manual + Min/Max spin boxes
   - Number of color levels (2..32)
   - Show Min Marker / Show Max Marker checkboxes
   - Show Element Value Labels (top-N spin)
   - Coordinate System (Global / Element)
   - Deformed-shape toggle + scale
   - Link buttons that open the Animation and Vector docks

The widget is **presentation-only** -- it emits signals; MainWindow
handles the actual contour rendering / actor builders. This keeps the
file independent of the rendering pipeline and easy to test.

Signals:
- ``request_contour(int subcase, str kind, str component)``
- ``request_vector(int subcase, str kind, str component)``
- ``request_deform(int subcase, str kind, str component)``
- ``request_animate(int subcase, str kind, str component)``
- ``copy_values(int subcase, str kind, str component)``
- ``color_range_changed(bool auto, float vmin, float vmax)``
- ``levels_changed(int n)``
- ``markers_changed(bool show_min, bool show_max)``
- ``element_labels_changed(bool show, int top_n)``
- ``deformation_changed(bool show, float scale)``
- ``open_animation_dock()``
- ``open_vector_dock()``

Where ``kind`` is one of:
``'displacement'``, ``'stress'``, ``'spc_forces'``, ``'eigenvector'``.
``component`` depends on kind (e.g. T1/T2/T3/Mag for displacement;
von_mises/max_principal/oxx/... for stress; mode 1/2/.. for
eigenvectors).
"""

from __future__ import annotations

from typing import Optional

from PySide6 import QtCore
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QTreeWidget, QTreeWidgetItem, QGroupBox, QCheckBox, QSpinBox,
    QDoubleSpinBox, QComboBox, QSizePolicy, QMenu, QFrame,
)


# Component lists per result kind. Used by _populate_subcase_tree.
DISP_COMPONENTS = ['T1', 'T2', 'T3', 'TMAG', 'R1', 'R2', 'R3', 'RMAG']
STRESS_COMPONENTS = ['von_mises', 'max_principal', 'min_principal',
                     'oxx', 'oyy', 'txy']


class ResultsTab(QWidget):
    """The new Results sidebar widget."""

    # ---- Tree-driven render requests ----
    request_contour = Signal(int, str, str)
    request_vector = Signal(int, str, str)
    request_deform = Signal(int, str, str)
    request_animate = Signal(int, str, str)
    copy_values = Signal(int, str, str)

    # ---- Control-panel state changes ----
    color_range_changed = Signal(bool, float, float)   # auto, vmin, vmax
    levels_changed = Signal(int)
    markers_changed = Signal(bool, bool)               # show_min, show_max
    element_labels_changed = Signal(bool, int)          # show, top_n
    deformation_changed = Signal(bool, float)          # show, scale
    open_animation_dock = Signal()
    open_vector_dock = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._current = None  # tuple (subcase_id, kind, component) or None
        self._building_tree = False

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # --- Tree ---
        self._tree = QTreeWidget()
        self._tree.setHeaderHidden(True)
        self._tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self._tree.customContextMenuRequested.connect(
            self._show_context_menu)
        self._tree.itemDoubleClicked.connect(self._on_item_double_clicked)
        self._tree.itemSelectionChanged.connect(self._on_selection_changed)
        outer.addWidget(self._tree, 1)

        # --- Controls panel ---
        controls_outer = QVBoxLayout()
        controls_outer.setSpacing(4)

        # Color range
        cr_box = QGroupBox("Color Range")
        cr_lay = QFormLayout(cr_box)
        cr_lay.setContentsMargins(6, 4, 6, 4)
        self._auto_check = QCheckBox("Auto (from data)")
        self._auto_check.setChecked(True)
        self._auto_check.toggled.connect(self._on_range_changed)
        cr_lay.addRow(self._auto_check)
        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(-1e30, 1e30)
        self._min_spin.setDecimals(6)
        self._min_spin.setEnabled(False)
        self._min_spin.editingFinished.connect(self._on_range_changed)
        cr_lay.addRow("Min:", self._min_spin)
        self._max_spin = QDoubleSpinBox()
        self._max_spin.setRange(-1e30, 1e30)
        self._max_spin.setDecimals(6)
        self._max_spin.setEnabled(False)
        self._max_spin.editingFinished.connect(self._on_range_changed)
        cr_lay.addRow("Max:", self._max_spin)
        self._auto_check.toggled.connect(
            lambda on: self._min_spin.setEnabled(not on))
        self._auto_check.toggled.connect(
            lambda on: self._max_spin.setEnabled(not on))
        self._levels_spin = QSpinBox()
        self._levels_spin.setRange(2, 32)
        self._levels_spin.setValue(9)
        self._levels_spin.valueChanged.connect(self.levels_changed.emit)
        cr_lay.addRow("Levels:", self._levels_spin)
        controls_outer.addWidget(cr_box)

        # Annotations: min/max markers + element labels
        annot_box = QGroupBox("Annotations")
        annot_lay = QFormLayout(annot_box)
        annot_lay.setContentsMargins(6, 4, 6, 4)
        self._min_marker = QCheckBox("Show Min Marker")
        self._max_marker = QCheckBox("Show Max Marker")
        self._min_marker.toggled.connect(self._on_markers_changed)
        self._max_marker.toggled.connect(self._on_markers_changed)
        annot_lay.addRow(self._min_marker)
        annot_lay.addRow(self._max_marker)
        self._labels_check = QCheckBox("Show Element Value Labels")
        self._labels_check.toggled.connect(self._on_labels_changed)
        annot_lay.addRow(self._labels_check)
        self._labels_top_n = QSpinBox()
        self._labels_top_n.setRange(0, 1_000_000)
        self._labels_top_n.setValue(0)
        self._labels_top_n.setSpecialValueText("All visible elements")
        self._labels_top_n.setSuffix(" (top-N)")
        self._labels_top_n.valueChanged.connect(self._on_labels_changed)
        annot_lay.addRow("Limit:", self._labels_top_n)
        controls_outer.addWidget(annot_box)

        # Deformed shape
        deform_box = QGroupBox("Deformed Shape")
        deform_lay = QFormLayout(deform_box)
        deform_lay.setContentsMargins(6, 4, 6, 4)
        self._deform_check = QCheckBox("Show deformed shape")
        self._deform_check.toggled.connect(self._on_deformation_changed)
        deform_lay.addRow(self._deform_check)
        self._deform_scale = QDoubleSpinBox()
        self._deform_scale.setRange(0.0, 1.0e6)
        self._deform_scale.setDecimals(4)
        self._deform_scale.setValue(1.0)
        self._deform_scale.valueChanged.connect(
            self._on_deformation_changed)
        deform_lay.addRow("Scale:", self._deform_scale)
        controls_outer.addWidget(deform_box)

        # Coordinate system (stub for v5.1.x)
        cs_box = QGroupBox("Coordinate System (stub)")
        cs_lay = QFormLayout(cs_box)
        cs_lay.setContentsMargins(6, 4, 6, 4)
        self._cs_combo = QComboBox()
        self._cs_combo.addItems(["Global"])
        self._cs_combo.setEnabled(False)
        self._cs_combo.setToolTip(
            "Element-local rotation lands in v5.1.x; v5.1.0 contours "
            "in the global coordinate system only.")
        cs_lay.addRow("Output in:", self._cs_combo)
        controls_outer.addWidget(cs_box)

        # Dock links
        link_box = QGroupBox("Related panels")
        link_lay = QVBoxLayout(link_box)
        link_lay.setContentsMargins(6, 4, 6, 4)
        open_anim = QPushButton("Open Animation Timeline…")
        open_anim.clicked.connect(self.open_animation_dock.emit)
        link_lay.addWidget(open_anim)
        open_vec = QPushButton("Open Vector Overlay…")
        open_vec.clicked.connect(self.open_vector_dock.emit)
        link_lay.addWidget(open_vec)
        controls_outer.addWidget(link_box)

        controls_outer.addStretch(0)
        outer.addLayout(controls_outer, 0)

    # ------------------------------------------------------------------
    # Tree population
    # ------------------------------------------------------------------

    def populate_from_results(self, op2_results: Optional[dict],
                              file_label: str = ""):
        """Rebuild the tree from an ``op2_results`` dict
        (the same shape produced by node_runner.model.load_op2_results
        and node_runner.solve.mystran_results.load_mystran_results).

        Passes None to clear the tree.
        """
        self._building_tree = True
        try:
            self._tree.clear()
            if not op2_results or not op2_results.get('subcases'):
                placeholder = QTreeWidgetItem(["(no results loaded)"])
                self._tree.addTopLevelItem(placeholder)
                return
            root = QTreeWidgetItem([file_label or "Results"])
            self._tree.addTopLevelItem(root)
            root.setExpanded(True)
            for sc_id, sc_data in sorted(
                    op2_results['subcases'].items(),
                    key=lambda kv: kv[0]):
                sc_item = QTreeWidgetItem(
                    [f"Subcase {sc_id}"])
                sc_item.setData(0, Qt.UserRole, ('subcase', int(sc_id)))
                root.addChild(sc_item)
                sc_item.setExpanded(True)
                # Displacements
                if sc_data.get('displacements'):
                    self._add_components(
                        sc_item, sc_id, 'displacement', 'Displacement',
                        DISP_COMPONENTS)
                # SPC forces
                if sc_data.get('spc_forces'):
                    self._add_components(
                        sc_item, sc_id, 'spc_forces', 'SPC Forces',
                        DISP_COMPONENTS[:7])
                # Stress
                if sc_data.get('stresses'):
                    self._add_components(
                        sc_item, sc_id, 'stress', 'Stress',
                        STRESS_COMPONENTS)
                # Eigenvectors -> one node per mode
                eigs = sc_data.get('eigenvectors') or []
                if eigs:
                    eig_root = QTreeWidgetItem([f"Eigenvectors ({len(eigs)} modes)"])
                    eig_root.setData(0, Qt.UserRole,
                                     ('eigenvectors_group', int(sc_id)))
                    sc_item.addChild(eig_root)
                    for mode_idx in range(len(eigs)):
                        mode_label = f"Mode {mode_idx + 1}"
                        freqs = sc_data.get('frequencies') or []
                        if mode_idx < len(freqs):
                            mode_label += f"  ({freqs[mode_idx]:.3f} Hz)"
                        mode_node = QTreeWidgetItem([mode_label])
                        mode_node.setData(
                            0, Qt.UserRole,
                            ('eigenvector', int(sc_id), int(mode_idx)))
                        eig_root.addChild(mode_node)
        finally:
            self._building_tree = False

    def _add_components(self, parent_item, subcase_id, kind, label,
                        components):
        kind_item = QTreeWidgetItem([label])
        kind_item.setData(0, Qt.UserRole, ('kind', int(subcase_id), kind))
        parent_item.addChild(kind_item)
        kind_item.setExpanded(False)
        for comp in components:
            ci = QTreeWidgetItem([comp])
            ci.setData(0, Qt.UserRole,
                       ('component', int(subcase_id), kind, comp))
            kind_item.addChild(ci)

    # ------------------------------------------------------------------
    # Selection + context menu
    # ------------------------------------------------------------------

    def _on_selection_changed(self):
        if self._building_tree:
            return
        items = self._tree.selectedItems()
        if not items:
            self._current = None
            return
        data = items[0].data(0, Qt.UserRole)
        if isinstance(data, tuple) and data and data[0] == 'component':
            _, sc, kind, comp = data
            self._current = (sc, kind, comp)
        elif isinstance(data, tuple) and data and data[0] == 'eigenvector':
            _, sc, mode_idx = data
            self._current = (sc, 'eigenvector', f"Mode {mode_idx + 1}")
        else:
            self._current = None

    def _on_item_double_clicked(self, item, _col):
        data = item.data(0, Qt.UserRole)
        if not (isinstance(data, tuple) and data):
            return
        if data[0] == 'component':
            _, sc, kind, comp = data
            self.request_contour.emit(int(sc), kind, comp)
        elif data[0] == 'eigenvector':
            _, sc, mode_idx = data
            self.request_contour.emit(
                int(sc), 'eigenvector', f"Mode {mode_idx + 1}")

    def _show_context_menu(self, pos):
        item = self._tree.itemAt(pos)
        if not item:
            return
        data = item.data(0, Qt.UserRole)
        if not (isinstance(data, tuple) and data):
            return
        if data[0] not in ('component', 'eigenvector'):
            return

        if data[0] == 'component':
            _, sc, kind, comp = data
        else:
            _, sc, mode_idx = data
            kind, comp = 'eigenvector', f"Mode {mode_idx + 1}"

        menu = QMenu(self._tree)
        contour_act = menu.addAction("Contour")
        contour_act.triggered.connect(
            lambda: self.request_contour.emit(int(sc), kind, comp))
        vector_act = menu.addAction("Vector overlay")
        vector_act.setEnabled(kind in ('displacement', 'spc_forces',
                                       'eigenvector'))
        vector_act.triggered.connect(
            lambda: self.request_vector.emit(int(sc), kind, comp))
        deform_act = menu.addAction("Apply as deformed shape")
        deform_act.setEnabled(kind in ('displacement', 'eigenvector'))
        deform_act.triggered.connect(
            lambda: self.request_deform.emit(int(sc), kind, comp))
        animate_act = menu.addAction("Animate")
        animate_act.setEnabled(kind in ('displacement', 'eigenvector'))
        animate_act.triggered.connect(
            lambda: self.request_animate.emit(int(sc), kind, comp))
        menu.addSeparator()
        copy_act = menu.addAction("Copy values to clipboard (CSV)")
        copy_act.triggered.connect(
            lambda: self.copy_values.emit(int(sc), kind, comp))
        menu.exec(self._tree.viewport().mapToGlobal(pos))

    # ------------------------------------------------------------------
    # Control signals
    # ------------------------------------------------------------------

    def _on_range_changed(self, *_a):
        auto = self._auto_check.isChecked()
        vmin = float(self._min_spin.value())
        vmax = float(self._max_spin.value())
        self.color_range_changed.emit(auto, vmin, vmax)

    def _on_markers_changed(self, *_a):
        self.markers_changed.emit(
            self._min_marker.isChecked(),
            self._max_marker.isChecked())

    def _on_labels_changed(self, *_a):
        self.element_labels_changed.emit(
            self._labels_check.isChecked(),
            int(self._labels_top_n.value()))

    def _on_deformation_changed(self, *_a):
        self.deformation_changed.emit(
            self._deform_check.isChecked(),
            float(self._deform_scale.value()))

    # ------------------------------------------------------------------
    # Helpers used by MainWindow after a contour lands
    # ------------------------------------------------------------------

    def set_data_range(self, vmin: float, vmax: float):
        """Called by MainWindow after a contour build so the Manual
        spinboxes start at the data range when the user flips Auto OFF."""
        self._min_spin.blockSignals(True)
        self._max_spin.blockSignals(True)
        try:
            self._min_spin.setValue(float(vmin))
            self._max_spin.setValue(float(vmax))
        finally:
            self._min_spin.blockSignals(False)
            self._max_spin.blockSignals(False)
