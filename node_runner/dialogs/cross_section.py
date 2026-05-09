"""Cross-Section dock + Define-Plane dialog.

Femap-style redesign:
  * The old non-modal `CrossSectionDialog` (covers the model, single
    method = global axis only) is replaced with `CrossSectionDock`, a
    `QDockWidget` that tabifies alongside Results / Animation / Vectors.
  * A "Define Plane..." button opens `DefinePlaneDialog`, a small modal
    with a method picker (Global Axis, 3 Points, Normal + Point,
    Coordinate Value, Two Nodes). The dialog returns a `PlaneDefinition`;
    the dock pushes it into the controller, then scrubs along the
    resolved normal with the slider.

`CrossSectionDialog` is kept as a thin back-compat alias pointing to the
dock so any callers/tests that import the old name still work.
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt, QSettings, Signal
from PySide6.QtWidgets import (
    QDockWidget, QWidget, QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
    QGridLayout, QLabel, QComboBox, QSlider, QPushButton, QCheckBox,
    QDoubleSpinBox, QSpinBox, QFrame, QStackedWidget, QGroupBox,
    QDialogButtonBox, QLineEdit, QMessageBox,
)

from node_runner.cross_section import (
    AXIS_PRESETS, PlaneDefinition, PlaneDefinitionError,
    METHOD_AXIS, METHOD_THREE_POINT, METHOD_NORMAL_POINT,
    METHOD_COORD_VALUE, METHOD_TWO_NODES, METHOD_LABELS,
)


SLIDER_RANGE = 1000  # finer slider resolution = smoother scrubbing


# ---------------------------------------------------------------------------
# Define-Plane dialog (method picker)
# ---------------------------------------------------------------------------

class DefinePlaneDialog(QDialog):
    """Modal dialog: pick a method, fill its inputs, return a PlaneDefinition.

    Called from the dock's "Define Plane..." button. After ``exec()`` the
    caller reads ``result_definition`` (or None on cancel/error). The
    selected method is persisted to QSettings so reopening defaults to
    the most recent choice.
    """

    SETTINGS_KEY = "cross_section/last_method"

    def __init__(self, main_window, current: Optional[PlaneDefinition] = None,
                 parent=None):
        super().__init__(parent or main_window)
        self.setWindowTitle("Define Plane")
        self.setModal(True)
        self.setMinimumWidth(420)
        self._main_window = main_window
        self.result_definition: Optional[PlaneDefinition] = None

        layout = QVBoxLayout(self)

        # Method picker.
        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method:"))
        self._method_combo = QComboBox()
        for key in (METHOD_AXIS, METHOD_COORD_VALUE, METHOD_THREE_POINT,
                    METHOD_NORMAL_POINT, METHOD_TWO_NODES):
            self._method_combo.addItem(METHOD_LABELS[key], key)
        method_row.addWidget(self._method_combo, 1)
        layout.addLayout(method_row)

        # Stacked sub-forms - one per method. Index aligns with combo.
        self._stack = QStackedWidget()
        layout.addWidget(self._stack)

        self._page_axis = self._build_axis_page()
        self._page_coord = self._build_coord_page()
        self._page_three = self._build_three_point_page()
        self._page_normal = self._build_normal_point_page()
        self._page_two = self._build_two_nodes_page()

        self._stack.addWidget(self._page_axis)
        self._stack.addWidget(self._page_coord)
        self._stack.addWidget(self._page_three)
        self._stack.addWidget(self._page_normal)
        self._stack.addWidget(self._page_two)

        self._method_combo.currentIndexChanged.connect(self._stack.setCurrentIndex)

        # Initialize with the current definition or a sensible default.
        if current is not None:
            self._prefill(current)
        else:
            settings = QSettings()
            last = settings.value(self.SETTINGS_KEY, METHOD_AXIS)
            self._select_method(last)

        # OK / Cancel.
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._on_accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    # ----- page builders -----

    def _build_axis_page(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self._axis_combo = QComboBox()
        for label, _ in AXIS_PRESETS:
            self._axis_combo.addItem(label)
        form.addRow("Axis:", self._axis_combo)
        info = QLabel("Plane perpendicular to the chosen global axis, passing through the origin. "
                      "Use the dock's slider to scrub the cut along the axis.")
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        form.addRow(info)
        return page

    def _build_coord_page(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self._coord_axis = QComboBox()
        for label, _ in AXIS_PRESETS[::2]:  # +X, +Y, +Z (sign comes from value)
            self._coord_axis.addItem(label)
        form.addRow("Axis:", self._coord_axis)
        self._coord_value = QDoubleSpinBox()
        self._coord_value.setDecimals(4)
        self._coord_value.setRange(-1e9, 1e9)
        self._coord_value.setSingleStep(1.0)
        form.addRow("Value:", self._coord_value)
        info = QLabel("Plane perpendicular to the chosen axis, at the typed coordinate. "
                      "Example: Axis +X, Value 1500 -> plane at X = 1500.")
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        form.addRow(info)
        return page

    def _build_three_point_page(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self._three_n1 = QSpinBox(); self._three_n1.setRange(1, 999_999_999)
        self._three_n2 = QSpinBox(); self._three_n2.setRange(1, 999_999_999)
        self._three_n3 = QSpinBox(); self._three_n3.setRange(1, 999_999_999)
        form.addRow("Node 1:", self._three_n1)
        form.addRow("Node 2:", self._three_n2)
        form.addRow("Node 3:", self._three_n3)
        info = QLabel("Plane passes through the three node IDs. Pick non-colinear nodes.")
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        form.addRow(info)
        return page

    def _build_normal_point_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        n_box = QGroupBox("Normal vector")
        n_form = QFormLayout(n_box)
        self._np_nx = QDoubleSpinBox(); self._np_nx.setDecimals(6); self._np_nx.setRange(-1e6, 1e6); self._np_nx.setValue(1.0)
        self._np_ny = QDoubleSpinBox(); self._np_ny.setDecimals(6); self._np_ny.setRange(-1e6, 1e6)
        self._np_nz = QDoubleSpinBox(); self._np_nz.setDecimals(6); self._np_nz.setRange(-1e6, 1e6)
        n_form.addRow("Nx:", self._np_nx)
        n_form.addRow("Ny:", self._np_ny)
        n_form.addRow("Nz:", self._np_nz)
        layout.addWidget(n_box)

        p_box = QGroupBox("Point on plane")
        p_form = QFormLayout(p_box)
        self._np_use_node = QCheckBox("Use a node as the point")
        self._np_node_id = QSpinBox(); self._np_node_id.setRange(1, 999_999_999); self._np_node_id.setEnabled(False)
        self._np_use_node.toggled.connect(self._np_node_id.setEnabled)
        self._np_use_node.toggled.connect(lambda on: (
            self._np_x.setEnabled(not on),
            self._np_y.setEnabled(not on),
            self._np_z.setEnabled(not on),
        ))
        self._np_x = QDoubleSpinBox(); self._np_x.setDecimals(6); self._np_x.setRange(-1e9, 1e9)
        self._np_y = QDoubleSpinBox(); self._np_y.setDecimals(6); self._np_y.setRange(-1e9, 1e9)
        self._np_z = QDoubleSpinBox(); self._np_z.setDecimals(6); self._np_z.setRange(-1e9, 1e9)
        p_form.addRow(self._np_use_node)
        p_form.addRow("Node ID:", self._np_node_id)
        p_form.addRow("X:", self._np_x)
        p_form.addRow("Y:", self._np_y)
        p_form.addRow("Z:", self._np_z)
        layout.addWidget(p_box)

        info = QLabel("Plane defined by a normal vector and a base point on the plane.")
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)
        return page

    def _build_two_nodes_page(self) -> QWidget:
        page = QWidget()
        form = QFormLayout(page)
        self._two_n1 = QSpinBox(); self._two_n1.setRange(1, 999_999_999)
        self._two_n2 = QSpinBox(); self._two_n2.setRange(1, 999_999_999)
        form.addRow("Node A:", self._two_n1)
        form.addRow("Node B:", self._two_n2)
        info = QLabel("Plane perpendicular to the line A->B, through its midpoint. "
                      "Useful for cuts perpendicular to a beam.")
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        form.addRow(info)
        return page

    # ----- helpers -----

    def _select_method(self, method: str):
        for i in range(self._method_combo.count()):
            if self._method_combo.itemData(i) == method:
                self._method_combo.setCurrentIndex(i)
                self._stack.setCurrentIndex(i)
                return
        self._method_combo.setCurrentIndex(0)
        self._stack.setCurrentIndex(0)

    def _prefill(self, defn: PlaneDefinition):
        self._select_method(defn.method)
        if defn.method == METHOD_AXIS:
            idx = self._axis_combo.findText(defn.axis_label)
            if idx >= 0:
                self._axis_combo.setCurrentIndex(idx)
        elif defn.method == METHOD_COORD_VALUE:
            idx = self._coord_axis.findText(defn.axis_label.lstrip("+-").rjust(2, "+"))
            if idx >= 0:
                self._coord_axis.setCurrentIndex(idx)
            self._coord_value.setValue(float(defn.value))
        elif defn.method == METHOD_THREE_POINT:
            ids = list(defn.node_ids) + [1, 1, 1]
            self._three_n1.setValue(int(ids[0]))
            self._three_n2.setValue(int(ids[1]))
            self._three_n3.setValue(int(ids[2]))
        elif defn.method == METHOD_NORMAL_POINT:
            self._np_nx.setValue(float(defn.normal[0]))
            self._np_ny.setValue(float(defn.normal[1]))
            self._np_nz.setValue(float(defn.normal[2]))
            self._np_x.setValue(float(defn.point[0]))
            self._np_y.setValue(float(defn.point[1]))
            self._np_z.setValue(float(defn.point[2]))
        elif defn.method == METHOD_TWO_NODES:
            ids = list(defn.node_ids) + [1, 1]
            self._two_n1.setValue(int(ids[0]))
            self._two_n2.setValue(int(ids[1]))

    def _model(self):
        gen = getattr(self._main_window, "current_generator", None)
        if gen is None:
            return None
        return getattr(gen, "model", None)

    def _node_xyz(self, nid: int):
        model = self._model()
        if model is None:
            return None
        node = model.nodes.get(int(nid))
        if node is None:
            return None
        try:
            xyz = node.get_position()
        except Exception:
            xyz = node.xyz
        return (float(xyz[0]), float(xyz[1]), float(xyz[2]))

    # ----- accept -----

    def _on_accept(self):
        method = self._method_combo.currentData()
        try:
            defn = self._build_definition(method)
            # Resolve once to validate (raises PlaneDefinitionError on bad input).
            defn.resolve(self._model())
        except PlaneDefinitionError as exc:
            QMessageBox.warning(self, "Define Plane", str(exc))
            return
        except Exception as exc:
            QMessageBox.warning(self, "Define Plane", f"Could not build plane: {exc}")
            return
        QSettings().setValue(self.SETTINGS_KEY, method)
        self.result_definition = defn
        self.accept()

    def _build_definition(self, method: str) -> PlaneDefinition:
        if method == METHOD_AXIS:
            return PlaneDefinition(method=method,
                                   axis_label=self._axis_combo.currentText())
        if method == METHOD_COORD_VALUE:
            return PlaneDefinition(method=method,
                                   axis_label=self._coord_axis.currentText(),
                                   value=float(self._coord_value.value()))
        if method == METHOD_THREE_POINT:
            return PlaneDefinition(method=method,
                                   node_ids=(int(self._three_n1.value()),
                                             int(self._three_n2.value()),
                                             int(self._three_n3.value())))
        if method == METHOD_NORMAL_POINT:
            normal = (float(self._np_nx.value()),
                      float(self._np_ny.value()),
                      float(self._np_nz.value()))
            if self._np_use_node.isChecked():
                xyz = self._node_xyz(int(self._np_node_id.value()))
                if xyz is None:
                    raise PlaneDefinitionError(
                        f"Node {self._np_node_id.value()} not found in model")
                point = xyz
            else:
                point = (float(self._np_x.value()),
                         float(self._np_y.value()),
                         float(self._np_z.value()))
            return PlaneDefinition(method=method, normal=normal, point=point)
        if method == METHOD_TWO_NODES:
            return PlaneDefinition(method=method,
                                   node_ids=(int(self._two_n1.value()),
                                             int(self._two_n2.value())))
        raise PlaneDefinitionError(f"Unknown method: {method}")


# ---------------------------------------------------------------------------
# Cross-Section dock
# ---------------------------------------------------------------------------

class CrossSectionDock(QDockWidget):
    """Compact dock with cross-section controls.

    Replaces the old non-modal dialog. The dock holds the active
    PlaneDefinition; the controller does the GL-side clipping.

    Tabifies alongside Results / Animation / Vectors on the right side
    of the main window (configured by MainWindow._build_cross_section_dock).
    """

    SETTINGS_KEY = "cross_section/last_definition"

    enable_changed = Signal(bool)

    def __init__(self, controller, main_window, parent=None):
        super().__init__("Cross Section", parent or main_window)
        self.setObjectName("CrossSectionDock")
        self.controller = controller
        self.main_window = main_window
        self._definition: Optional[PlaneDefinition] = (
            controller.definition or PlaneDefinition()
        )

        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Enable.
        self.enable_check = QCheckBox("Enable cross section")
        self.enable_check.setChecked(controller.is_enabled)
        self.enable_check.toggled.connect(self._on_enable_toggled)
        layout.addWidget(self.enable_check)

        # Active plane summary.
        summary_box = QGroupBox("Active plane")
        summary_layout = QVBoxLayout(summary_box)
        self._summary_label = QLabel(self._definition.summary())
        self._summary_label.setWordWrap(True)
        self._summary_label.setStyleSheet("font-weight: 500;")
        summary_layout.addWidget(self._summary_label)

        btn_row = QHBoxLayout()
        define_btn = QPushButton("Define Plane...")
        define_btn.clicked.connect(self._on_define_plane)
        flip_btn = QPushButton("Flip Normal")
        flip_btn.clicked.connect(self._on_flip)
        btn_row.addWidget(define_btn)
        btn_row.addWidget(flip_btn)
        summary_layout.addLayout(btn_row)
        layout.addWidget(summary_box)

        # Position slider + spin (along resolved normal).
        pos_box = QGroupBox("Position along normal")
        pos_layout = QVBoxLayout(pos_box)
        slider_row = QHBoxLayout()
        self.position_slider = QSlider(Qt.Horizontal)
        self.position_slider.setRange(0, SLIDER_RANGE)
        self.position_slider.setValue(int(controller.slider_fraction * SLIDER_RANGE))
        self.position_slider.valueChanged.connect(self._on_slider_changed)
        slider_row.addWidget(self.position_slider, 1)
        self.position_value = QDoubleSpinBox()
        self.position_value.setDecimals(3)
        self.position_value.setRange(-1e9, 1e9)
        self.position_value.setSingleStep(0.1)
        self.position_value.setMinimumWidth(110)
        self.position_value.editingFinished.connect(self._on_value_edited)
        slider_row.addWidget(self.position_value)
        pos_layout.addLayout(slider_row)

        center_btn = QPushButton("Center plane")
        center_btn.clicked.connect(self._on_center)
        pos_layout.addWidget(center_btn, alignment=Qt.AlignRight)
        layout.addWidget(pos_box)

        # Display options.
        opts_box = QGroupBox("Display")
        opts_layout = QVBoxLayout(opts_box)
        self.outline_check = QCheckBox("Show plane outline")
        self.outline_check.setChecked(False)
        self.outline_check.toggled.connect(self._on_outline_toggled)
        opts_layout.addWidget(self.outline_check)
        layout.addWidget(opts_box)

        info = QLabel(
            "GPU clipping - even huge models stay smooth. Cross-section "
            "stays on while this dock is open; uncheck Enable to restore "
            "the full view."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)
        layout.addStretch(1)

        self.setWidget(body)

        # Push the initial definition through so the slider matches reality
        # the moment the dock is shown.
        self._apply_definition_to_controller()
        self._apply_slider_to_controller(self.position_slider.value())

    # ----- helpers -----

    def _model(self):
        gen = getattr(self.main_window, "current_generator", None)
        if gen is None:
            return None
        return getattr(gen, "model", None)

    def _model_bounds_along_normal(self):
        """Return (lo, hi) extent of the loaded model projected onto the
        current plane's normal direction. Falls back to (-10, 10) if no
        model is loaded yet."""
        grid = getattr(self.main_window, "current_grid", None)
        if grid is None or grid.n_points == 0:
            return -10.0, 10.0
        b = grid.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
        nx, ny, nz = self.controller.plane.GetNormal()
        corners = [
            (b[0], b[2], b[4]), (b[1], b[2], b[4]),
            (b[0], b[3], b[4]), (b[1], b[3], b[4]),
            (b[0], b[2], b[5]), (b[1], b[2], b[5]),
            (b[0], b[3], b[5]), (b[1], b[3], b[5]),
        ]
        projs = [cx * nx + cy * ny + cz * nz for cx, cy, cz in corners]
        return min(projs), max(projs)

    def _apply_definition_to_controller(self):
        try:
            self.controller.apply_definition(self._definition, self._model())
        except PlaneDefinitionError as exc:
            QMessageBox.warning(self, "Cross Section", str(exc))
            # Fall back to a safe default.
            self._definition = PlaneDefinition()
            self.controller.apply_definition(self._definition, None)
        self._summary_label.setText(self._definition.summary())

    def _apply_slider_to_controller(self, slider_value):
        frac = slider_value / float(SLIDER_RANGE)
        self.controller.set_slider_fraction(frac)
        lo, hi = self._model_bounds_along_normal()
        # Translate fraction into a world-space distance along the normal
        # and place the plane origin at base + that displacement.
        # For axis/coord-value methods the base origin is on-axis already;
        # for everything else we offset along the resolved normal from the
        # definition's resolved origin (which sits ON the plane).
        nx, ny, nz = self.controller.plane.GetNormal()
        # Use the resolved origin captured by the last apply_definition().
        if self._definition is not None and self._definition.method in (
                METHOD_THREE_POINT, METHOD_NORMAL_POINT, METHOD_TWO_NODES):
            base = self._definition._resolved_origin
            # Fraction 0.5 means "no offset from the picked plane".
            offset = (frac - 0.5) * (hi - lo)
            ox, oy, oz = (base[0] + nx * offset,
                          base[1] + ny * offset,
                          base[2] + nz * offset)
            display_value = offset
        else:
            dist = lo + frac * (hi - lo)
            ox, oy, oz = (nx * dist, ny * dist, nz * dist)
            display_value = dist
        self.controller.set_origin(ox, oy, oz)
        self.position_value.blockSignals(True)
        self.position_value.setValue(display_value)
        self.position_value.blockSignals(False)

    # ----- handlers -----

    def _on_enable_toggled(self, on):
        if on:
            self.controller.enable()
        else:
            self.controller.disable()
        self.enable_changed.emit(bool(on))

    def _on_define_plane(self):
        dlg = DefinePlaneDialog(self.main_window, current=self._definition,
                                parent=self)
        if dlg.exec() == QDialog.Accepted and dlg.result_definition is not None:
            self._definition = dlg.result_definition
            self._apply_definition_to_controller()
            # Reset slider to center for picked-plane methods, leave for
            # axis/coord (where 0.5 = model midpoint along axis).
            self.position_slider.blockSignals(True)
            self.position_slider.setValue(SLIDER_RANGE // 2)
            self.position_slider.blockSignals(False)
            self._apply_slider_to_controller(self.position_slider.value())
            # Auto-enable when the user goes through the trouble of
            # picking a plane (matches user expectation).
            if not self.enable_check.isChecked():
                self.enable_check.setChecked(True)

    def _on_flip(self):
        self.controller.flip()

    def _on_slider_changed(self, value):
        self._apply_slider_to_controller(value)

    def _on_value_edited(self):
        lo, hi = self._model_bounds_along_normal()
        v = self.position_value.value()
        if self._definition is not None and self._definition.method in (
                METHOD_THREE_POINT, METHOD_NORMAL_POINT, METHOD_TWO_NODES):
            # Treat the spinbox as offset from the picked plane.
            span = (hi - lo) if hi > lo else 1.0
            frac = 0.5 + v / span
        else:
            frac = (v - lo) / (hi - lo) if hi > lo else 0.5
        frac = max(0.0, min(1.0, frac))
        self.position_slider.blockSignals(True)
        self.position_slider.setValue(int(frac * SLIDER_RANGE))
        self.position_slider.blockSignals(False)
        self._apply_slider_to_controller(self.position_slider.value())

    def _on_center(self):
        self.position_slider.setValue(SLIDER_RANGE // 2)

    def _on_outline_toggled(self, on):
        self.controller.set_outline_visible(bool(on))

    # ----- public -----

    def show_and_raise(self):
        self.show()
        self.raise_()
        # If tabified, bring this tab to front.
        try:
            self.main_window.tabifyDockWidget(self, self)  # no-op but safe
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Backwards-compat shim
# ---------------------------------------------------------------------------

# Several tests/callers still import CrossSectionDialog. Keep the symbol
# pointing at the new dock so old code paths that just need a "controller
# UI" find something they can construct. Anything that genuinely needs the
# old QDialog API has been migrated to MainWindow._open_cross_section_dock.
CrossSectionDialog = CrossSectionDock
