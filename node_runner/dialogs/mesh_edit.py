"""Dialogs for the Theme A mesh-editing tools.

Smaller and more focused than the entity creators / editors - each one
just gathers a few parameters before invoking a command on the current
selection. The bigger UI (selection bar, undo/redo wiring) lives in
MainWindow; these are pure parameter forms.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
    QSpinBox, QDoubleSpinBox, QCheckBox, QComboBox, QGroupBox,
)


class SmoothNodesDialog(QDialog):
    """Iterations, factor, pin free edges for Laplacian smoothing."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Smooth Nodes")
        self.setMinimumWidth(320)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Laplacian smoothing moves each selected node toward the "
            "centroid of its mesh neighbors. Boundary (free-edge) nodes "
            "are pinned by default."
        ))

        grid = QGridLayout()
        grid.addWidget(QLabel("Iterations:"), 0, 0)
        self.iter_spin = QSpinBox()
        self.iter_spin.setRange(1, 100)
        self.iter_spin.setValue(5)
        grid.addWidget(self.iter_spin, 0, 1)

        grid.addWidget(QLabel("Step factor (0..1):"), 1, 0)
        self.factor_spin = QDoubleSpinBox()
        self.factor_spin.setRange(0.0, 1.0)
        self.factor_spin.setSingleStep(0.05)
        self.factor_spin.setValue(0.5)
        grid.addWidget(self.factor_spin, 1, 1)

        layout.addLayout(grid)

        self.pin_check = QCheckBox("Pin nodes on free edges")
        self.pin_check.setChecked(True)
        layout.addWidget(self.pin_check)

        row = QHBoxLayout()
        row.addStretch(1)
        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        ok = QPushButton("Smooth"); ok.setDefault(True); ok.clicked.connect(self.accept)
        row.addWidget(cancel); row.addWidget(ok)
        layout.addLayout(row)

    @property
    def iterations(self):
        return int(self.iter_spin.value())

    @property
    def factor(self):
        return float(self.factor_spin.value())

    @property
    def pin_free_edges(self):
        return self.pin_check.isChecked()


class MirrorElementsDialog(QDialog):
    """Plane normal (X/Y/Z) and offset for element mirroring."""

    def __init__(self, suggested_plane='Y', suggested_value=0.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mirror Elements")
        self.setMinimumWidth(320)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Reflect the selected elements across an axis-aligned plane. "
            "Connectivity is reversed so shell normals stay outward."
        ))

        grid = QGridLayout()
        grid.addWidget(QLabel("Plane normal:"), 0, 0)
        self.plane_combo = QComboBox()
        self.plane_combo.addItems(["X", "Y", "Z"])
        idx = self.plane_combo.findText(suggested_plane.upper())
        if idx >= 0:
            self.plane_combo.setCurrentIndex(idx)
        grid.addWidget(self.plane_combo, 0, 1)

        grid.addWidget(QLabel("Plane offset:"), 1, 0)
        self.value_spin = QDoubleSpinBox()
        self.value_spin.setRange(-1e9, 1e9)
        self.value_spin.setDecimals(6)
        self.value_spin.setValue(float(suggested_value))
        grid.addWidget(self.value_spin, 1, 1)

        grid.addWidget(QLabel("Weld tolerance:"), 2, 0)
        self.weld_tol_spin = QDoubleSpinBox()
        self.weld_tol_spin.setRange(0.0, 1e3)
        self.weld_tol_spin.setDecimals(8)
        self.weld_tol_spin.setValue(1e-6)
        grid.addWidget(self.weld_tol_spin, 2, 1)

        layout.addLayout(grid)

        row = QHBoxLayout()
        row.addStretch(1)
        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        ok = QPushButton("Mirror"); ok.setDefault(True); ok.clicked.connect(self.accept)
        row.addWidget(cancel); row.addWidget(ok)
        layout.addLayout(row)

    @property
    def plane(self):
        return self.plane_combo.currentText()

    @property
    def value(self):
        return float(self.value_spin.value())

    @property
    def weld_tol(self):
        return float(self.weld_tol_spin.value())


class CopyElementsDialog(QDialog):
    """Translation vector (dx, dy, dz) for the element-copy command."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Copy Elements")
        self.setMinimumWidth(320)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Duplicate the selected elements (with their nodes) and "
            "translate the copies by the offset below."
        ))

        grid = QGridLayout()
        labels = ("dx", "dy", "dz")
        self.spins = []
        for i, lbl in enumerate(labels):
            grid.addWidget(QLabel(lbl + ":"), i, 0)
            sp = QDoubleSpinBox()
            sp.setRange(-1e9, 1e9)
            sp.setDecimals(6)
            sp.setValue(0.0)
            grid.addWidget(sp, i, 1)
            self.spins.append(sp)
        layout.addLayout(grid)

        row = QHBoxLayout()
        row.addStretch(1)
        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        ok = QPushButton("Copy"); ok.setDefault(True); ok.clicked.connect(self.accept)
        row.addWidget(cancel); row.addWidget(ok)
        layout.addLayout(row)

    @property
    def translation(self):
        return tuple(float(s.value()) for s in self.spins)


class CombineTriasDialog(QDialog):
    """Tolerance angle for combine operation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Combine TRIA3 -> QUAD4")
        self.setMinimumWidth(320)
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Pair selected CTRIA3s that share an edge and are nearly "
            "coplanar. Each successful pair becomes one CQUAD4."
        ))

        row = QHBoxLayout()
        row.addWidget(QLabel("Coplanarity tolerance (deg):"))
        self.angle_spin = QDoubleSpinBox()
        self.angle_spin.setRange(0.0, 90.0)
        self.angle_spin.setDecimals(2)
        self.angle_spin.setValue(5.0)
        row.addWidget(self.angle_spin)
        layout.addLayout(row)

        btn_row = QHBoxLayout()
        btn_row.addStretch(1)
        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        ok = QPushButton("Combine"); ok.setDefault(True); ok.clicked.connect(self.accept)
        btn_row.addWidget(cancel); btn_row.addWidget(ok)
        layout.addLayout(btn_row)

    @property
    def angle_tol_deg(self):
        return float(self.angle_spin.value())
