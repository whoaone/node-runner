"""Theme B: Result browser dock and supporting widgets.

Hosts:
  * `ResultTableModel` - QAbstractTableModel backed by numpy arrays for
    fast scrolling on huge OP2 results.
  * `ResultBrowserDock` - the right-side dock with sortable table tabs
    for nodal and element results plus a built-in expression bar.
  * `AnimationTimelineWidget` - play/pause/scrub timeline for mode
    shapes, replacing the old GIF-export-only flow.
  * `VectorOverlayWidget` - checkboxes to toggle displacement /
    principal stress / reaction-force arrows on the 3D view.

These are pure Qt widgets - all model mutation goes through the
MainWindow which still owns the plotter and pyNastran objects.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import (
    Qt, QAbstractTableModel, QModelIndex, QSortFilterProxyModel, QTimer,
    Signal,
)
from PySide6.QtWidgets import (
    QDockWidget, QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QTableView, QHeaderView, QLabel, QLineEdit, QPushButton, QCheckBox,
    QSlider, QSpinBox, QDoubleSpinBox, QComboBox,
)


# ---------------------------------------------------------------------------
# Result table model
# ---------------------------------------------------------------------------

class ResultTableModel(QAbstractTableModel):
    """Read-only table backed by parallel numpy arrays.

    Construct with a list of column names and a dict {col_name: ndarray}.
    All arrays must be the same length. Sorting is delegated to
    QSortFilterProxyModel for free O(n log n) sorts on numpy data.
    """

    def __init__(self, columns: list[str], arrays: dict[str, np.ndarray],
                 id_column: str = 'ID', parent=None):
        super().__init__(parent)
        self._columns = list(columns)
        self._arrays = {c: np.asarray(arrays.get(c, np.array([]))) for c in columns}
        self._id_column = id_column
        self._n_rows = len(self._arrays[columns[0]]) if columns else 0

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else self._n_rows

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._columns)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return self._columns[section]
        return str(section + 1)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        if role not in (Qt.DisplayRole, Qt.UserRole, Qt.ToolTipRole):
            return None
        col = self._columns[index.column()]
        v = self._arrays[col][index.row()]
        if role == Qt.UserRole:
            return v
        if isinstance(v, (np.integer, int)):
            return str(int(v))
        try:
            return f"{float(v):.6g}"
        except (ValueError, TypeError):
            return str(v)

    def id_of_row(self, row: int):
        if 0 <= row < self._n_rows:
            arr = self._arrays.get(self._id_column)
            if arr is not None and row < len(arr):
                return int(arr[row])
        return None

    def row_of_id(self, eid: int):
        arr = self._arrays.get(self._id_column)
        if arr is None:
            return None
        match = np.where(arr == int(eid))[0]
        return int(match[0]) if match.size else None


# ---------------------------------------------------------------------------
# Result browser dock
# ---------------------------------------------------------------------------

class ResultBrowserDock(QDockWidget):
    """Sortable result table + expression bar, dockable on the main window.

    Signals
    -------
    entity_picked : (entity_type, entity_id)
        Emitted when a row's "Zoom to" button is clicked or a row is
        double-clicked. ``entity_type`` is 'Node' or 'Element'.
    expression_changed : str
        Emitted when the expression input loses focus or the user
        presses Enter.
    """

    entity_picked = Signal(str, int)
    expression_changed = Signal(str)

    def __init__(self, parent=None):
        super().__init__("Results", parent)
        self.setObjectName("ResultBrowserDock")
        self.setAllowedAreas(Qt.RightDockWidgetArea | Qt.LeftDockWidgetArea)

        wrapper = QWidget(); layout = QVBoxLayout(wrapper)
        layout.setContentsMargins(6, 6, 6, 6); layout.setSpacing(6)

        # Expression bar
        expr_row = QHBoxLayout()
        expr_row.addWidget(QLabel("Expression:"))
        self.expression_edit = QLineEdit()
        self.expression_edit.setPlaceholderText(
            "e.g. 1 - stress_vm/450e6   |   sqrt(disp_x**2 + disp_y**2)"
        )
        self.expression_edit.editingFinished.connect(
            lambda: self.expression_changed.emit(self.expression_edit.text().strip()),
        )
        expr_row.addWidget(self.expression_edit, 1)
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self._clear_expression)
        expr_row.addWidget(clear_btn)
        layout.addLayout(expr_row)

        # Tabs: nodal + element results
        self.tabs = QTabWidget()
        self.node_view = QTableView(); self._configure_view(self.node_view)
        self.elem_view = QTableView(); self._configure_view(self.elem_view)
        self.tabs.addTab(self.node_view, "Nodes")
        self.tabs.addTab(self.elem_view, "Elements")
        layout.addWidget(self.tabs, 1)

        self.setWidget(wrapper)
        self.hide()  # only show when results loaded

        # Models
        self._node_model = ResultTableModel(['NID'], {'NID': np.array([])}, id_column='NID', parent=self)
        self._elem_model = ResultTableModel(['EID'], {'EID': np.array([])}, id_column='EID', parent=self)
        self._set_model(self.node_view, self._node_model)
        self._set_model(self.elem_view, self._elem_model)

        self.node_view.doubleClicked.connect(self._on_node_double_clicked)
        self.elem_view.doubleClicked.connect(self._on_elem_double_clicked)

    def _configure_view(self, view: QTableView):
        view.setSortingEnabled(True)
        view.setSelectionBehavior(QTableView.SelectRows)
        view.setSelectionMode(QTableView.SingleSelection)
        view.setEditTriggers(QTableView.NoEditTriggers)
        view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        view.horizontalHeader().setStretchLastSection(True)

    def _set_model(self, view, source_model):
        proxy = QSortFilterProxyModel(self)
        proxy.setSourceModel(source_model)
        view.setModel(proxy)

    def _clear_expression(self):
        self.expression_edit.clear()
        self.expression_changed.emit("")

    def _on_node_double_clicked(self, idx):
        proxy = self.node_view.model()
        if proxy is None:
            return
        row = proxy.mapToSource(idx).row()
        nid = self._node_model.id_of_row(row)
        if nid is not None:
            self.entity_picked.emit('Node', nid)

    def _on_elem_double_clicked(self, idx):
        proxy = self.elem_view.model()
        if proxy is None:
            return
        row = proxy.mapToSource(idx).row()
        eid = self._elem_model.id_of_row(row)
        if eid is not None:
            self.entity_picked.emit('Element', eid)

    # --- public API used by MainWindow ---

    def update_nodal_results(self, columns: list[str], arrays: dict):
        self._node_model = ResultTableModel(columns, arrays, id_column='NID', parent=self)
        self._set_model(self.node_view, self._node_model)

    def update_element_results(self, columns: list[str], arrays: dict):
        self._elem_model = ResultTableModel(columns, arrays, id_column='EID', parent=self)
        self._set_model(self.elem_view, self._elem_model)

    def select_node(self, nid: int):
        row = self._node_model.row_of_id(nid)
        if row is None:
            return
        self.tabs.setCurrentWidget(self.node_view)
        proxy = self.node_view.model()
        if proxy is None:
            return
        proxy_idx = proxy.mapFromSource(self._node_model.index(row, 0))
        self.node_view.selectRow(proxy_idx.row())
        self.node_view.scrollTo(proxy_idx)

    def select_element(self, eid: int):
        row = self._elem_model.row_of_id(eid)
        if row is None:
            return
        self.tabs.setCurrentWidget(self.elem_view)
        proxy = self.elem_view.model()
        if proxy is None:
            return
        proxy_idx = proxy.mapFromSource(self._elem_model.index(row, 0))
        self.elem_view.selectRow(proxy_idx.row())
        self.elem_view.scrollTo(proxy_idx)


# ---------------------------------------------------------------------------
# Animation timeline
# ---------------------------------------------------------------------------

class AnimationTimelineWidget(QWidget):
    """Play / pause / scrub timeline for mode-shape animation.

    Signals
    -------
    phase_changed : float
        Current phase in [0, 1]. The MainWindow translates this into a
        deformation scale for the active mode.
    play_state_changed : bool
        True when playing, False when paused.
    """

    phase_changed = Signal(float)
    play_state_changed = Signal(bool)

    SLIDER_RES = 1000

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 2, 4, 2)

        self.play_btn = QPushButton("Play")
        self.play_btn.setCheckable(True)
        self.play_btn.toggled.connect(self._on_play_toggled)
        layout.addWidget(self.play_btn)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, self.SLIDER_RES)
        self.slider.valueChanged.connect(self._on_slider)
        layout.addWidget(self.slider, 1)

        self.speed_combo = QComboBox()
        self.speed_combo.addItems(["0.25x", "0.5x", "1.0x", "2.0x", "4.0x"])
        self.speed_combo.setCurrentIndex(2)
        self.speed_combo.currentIndexChanged.connect(self._on_speed_changed)
        layout.addWidget(self.speed_combo)

        self.loop_check = QCheckBox("Loop")
        self.loop_check.setChecked(True)
        layout.addWidget(self.loop_check)

        # Internal driving timer (60fps target).
        self._timer = QTimer(self)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self._tick)
        self._direction = +1
        self._phase = 0.0

    def _speed_factor(self):
        m = {0: 0.25, 1: 0.5, 2: 1.0, 3: 2.0, 4: 4.0}
        return m.get(self.speed_combo.currentIndex(), 1.0)

    def _on_speed_changed(self, _idx):
        # No-op: speed is read at each tick from _speed_factor()
        pass

    def _on_play_toggled(self, on):
        if on:
            self.play_btn.setText("Pause")
            self._timer.start()
        else:
            self.play_btn.setText("Play")
            self._timer.stop()
        self.play_state_changed.emit(on)

    def _on_slider(self, v):
        self._phase = v / float(self.SLIDER_RES)
        self.phase_changed.emit(self._phase)

    def _tick(self):
        # Step phase by ~speed * (16ms / 1000ms) over a 1-second loop.
        step = self._speed_factor() * 16.0 / 1000.0
        new_phase = self._phase + self._direction * step
        if new_phase >= 1.0:
            if self.loop_check.isChecked():
                new_phase = 0.0
            else:
                new_phase = 1.0
                self.play_btn.setChecked(False)
        elif new_phase <= 0.0:
            new_phase = 0.0
        self._phase = new_phase
        self.slider.blockSignals(True)
        self.slider.setValue(int(new_phase * self.SLIDER_RES))
        self.slider.blockSignals(False)
        self.phase_changed.emit(new_phase)

    @property
    def phase(self) -> float:
        return self._phase


# ---------------------------------------------------------------------------
# Vector overlay panel
# ---------------------------------------------------------------------------

class VectorOverlayWidget(QWidget):
    """Three checkboxes for displacement / principal stress / reaction
    force vector overlays plus a length-scale slider."""

    overlay_toggled = Signal(str, bool)   # ('disp' / 'prin' / 'reac', on)
    scale_changed = Signal(float)         # length scale fraction in [0.01, 0.5]

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self); layout.setContentsMargins(2, 2, 2, 2)

        row = QHBoxLayout()
        for key, label in (('disp', 'Disp'), ('prin', 'Principal'), ('reac', 'Reaction')):
            cb = QCheckBox(label)
            cb.toggled.connect(lambda on, k=key: self.overlay_toggled.emit(k, on))
            row.addWidget(cb)
            setattr(self, f"check_{key}", cb)
        layout.addLayout(row)

        scale_row = QHBoxLayout()
        scale_row.addWidget(QLabel("Length scale:"))
        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.01, 0.5)
        self.scale_spin.setSingleStep(0.01)
        self.scale_spin.setDecimals(2)
        self.scale_spin.setValue(0.05)
        self.scale_spin.valueChanged.connect(self.scale_changed.emit)
        scale_row.addWidget(self.scale_spin)
        scale_row.addStretch(1)
        layout.addLayout(scale_row)
