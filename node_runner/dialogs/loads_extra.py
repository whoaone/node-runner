"""Theme C dialogs: load and BC types not previously creatable in the UI.

PLOAD1 (distributed line load on bars/beams), PLOAD2 (uniform pressure
on shells via element list), SPCD (enforced displacement), MPC
(general multi-point constraint), RBAR / RBE1 / RSPLINE (rigid
elements not covered by the existing spider creator), and BOLT
(preload card).

Each dialog is a small focused form. The underlying card creation goes
through new Command subclasses in commands.py so undo/redo works.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QPushButton, QSpinBox, QDoubleSpinBox, QLineEdit, QComboBox,
    QCheckBox, QGroupBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QDialogButtonBox, QMessageBox,
)


def _parse_id_list(text: str) -> list[int]:
    """Parse '1, 2-5, 7' style ID lists. Tolerates spaces, mixed separators."""
    out = []
    if not text:
        return out
    for chunk in text.replace(';', ',').split(','):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '-' in chunk and not chunk.startswith('-'):
            try:
                a, b = chunk.split('-', 1)
                out.extend(range(int(a), int(b) + 1))
            except ValueError:
                pass
        else:
            try:
                out.append(int(chunk))
            except ValueError:
                pass
    return out


# ---------------------------------------------------------------------------
# C1. PLOAD1 - distributed line load
# ---------------------------------------------------------------------------

class CreatePload1Dialog(QDialog):
    """Distributed load along a bar or beam.

    Fields: load SID, element list, type (FX/FY/FZ/MX/MY/MZ), scale (LE/FR),
    X1 (start position), P1 (start magnitude), X2, P2.
    """

    LOAD_TYPES = ('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ',
                  'FXE', 'FYE', 'FZE', 'MXE', 'MYE', 'MZE')

    def __init__(self, parent=None, existing_sid=None):
        super().__init__(parent)
        self.setWindowTitle("Create Distributed Load (PLOAD1)")
        self.setMinimumWidth(420)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.sid_edit = QSpinBox(); self.sid_edit.setRange(1, 99_999_999)
        self.sid_edit.setValue(int(existing_sid) if existing_sid else 1)
        form.addRow("Load SID:", self.sid_edit)

        self.eids_edit = QLineEdit()
        self.eids_edit.setPlaceholderText("e.g. 101, 105-120")
        form.addRow("Bar/beam EIDs:", self.eids_edit)

        self.load_type = QComboBox()
        self.load_type.addItems(self.LOAD_TYPES)
        form.addRow("Load type:", self.load_type)

        self.scale = QComboBox()
        self.scale.addItems(["LE (length)", "FR (fraction)"])
        form.addRow("Scale:", self.scale)

        self.x1 = QDoubleSpinBox(); self.x1.setRange(-1e9, 1e9); self.x1.setDecimals(6); self.x1.setValue(0.0)
        self.p1 = QDoubleSpinBox(); self.p1.setRange(-1e12, 1e12); self.p1.setDecimals(6); self.p1.setValue(1.0)
        self.x2 = QDoubleSpinBox(); self.x2.setRange(-1e9, 1e9); self.x2.setDecimals(6); self.x2.setValue(1.0)
        self.p2 = QDoubleSpinBox(); self.p2.setRange(-1e12, 1e12); self.p2.setDecimals(6); self.p2.setValue(1.0)
        form.addRow("X1 (start):", self.x1)
        form.addRow("P1 (magnitude at X1):", self.p1)
        form.addRow("X2 (end, blank=duplicate of X1):", self.x2)
        form.addRow("P2 (magnitude at X2):", self.p2)

        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def sid(self):
        return int(self.sid_edit.value())

    @property
    def eids(self):
        return _parse_id_list(self.eids_edit.text())

    @property
    def values(self):
        return {
            'load_type': self.load_type.currentText(),
            'scale': 'LE' if self.scale.currentIndex() == 0 else 'FR',
            'x1': float(self.x1.value()),
            'p1': float(self.p1.value()),
            'x2': float(self.x2.value()),
            'p2': float(self.p2.value()),
        }


# ---------------------------------------------------------------------------
# C2. PLOAD2 - uniform pressure on shells via element list
# ---------------------------------------------------------------------------

class CreatePload2Dialog(QDialog):
    def __init__(self, parent=None, existing_sid=None):
        super().__init__(parent)
        self.setWindowTitle("Create Uniform Pressure (PLOAD2)")
        self.setMinimumWidth(380)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.sid_edit = QSpinBox(); self.sid_edit.setRange(1, 99_999_999)
        self.sid_edit.setValue(int(existing_sid) if existing_sid else 1)
        form.addRow("Load SID:", self.sid_edit)

        self.eids_edit = QLineEdit()
        self.eids_edit.setPlaceholderText("e.g. 200, 210-220")
        form.addRow("Shell EIDs:", self.eids_edit)

        self.pressure = QDoubleSpinBox()
        self.pressure.setRange(-1e12, 1e12); self.pressure.setDecimals(6); self.pressure.setValue(1.0)
        form.addRow("Pressure (force per unit area):", self.pressure)

        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def sid(self):
        return int(self.sid_edit.value())

    @property
    def eids(self):
        return _parse_id_list(self.eids_edit.text())

    @property
    def pressure_value(self):
        return float(self.pressure.value())


# ---------------------------------------------------------------------------
# C3. SPCD - enforced displacement
# ---------------------------------------------------------------------------

class CreateSpcdDialog(QDialog):
    def __init__(self, parent=None, existing_sid=None):
        super().__init__(parent)
        self.setWindowTitle("Create Enforced Displacement (SPCD)")
        self.setMinimumWidth(380)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.sid_edit = QSpinBox(); self.sid_edit.setRange(1, 99_999_999)
        self.sid_edit.setValue(int(existing_sid) if existing_sid else 1)
        form.addRow("Load SID:", self.sid_edit)

        self.nids_edit = QLineEdit()
        self.nids_edit.setPlaceholderText("e.g. 1, 5-10, 99")
        form.addRow("Node IDs:", self.nids_edit)

        self.dof_edit = QLineEdit("1")
        self.dof_edit.setPlaceholderText("e.g. 123 (TX TY TZ) or 13 (TX TZ)")
        form.addRow("Constrained DOFs:", self.dof_edit)

        self.value = QDoubleSpinBox()
        self.value.setRange(-1e9, 1e9); self.value.setDecimals(8); self.value.setValue(0.001)
        form.addRow("Enforced value:", self.value)

        layout.addLayout(form)

        info = QLabel(
            "SPCD goes into the LOAD case-control entry, not SPC. Make sure "
            "the analysis subcase references this SID under LOAD."
        )
        info.setWordWrap(True); info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def sid(self):
        return int(self.sid_edit.value())

    @property
    def nids(self):
        return _parse_id_list(self.nids_edit.text())

    @property
    def dof_string(self):
        return self.dof_edit.text().strip() or '1'

    @property
    def enforced_value(self):
        return float(self.value.value())


# ---------------------------------------------------------------------------
# C4. MPC - general multi-point constraint
# ---------------------------------------------------------------------------

class CreateMpcDialog(QDialog):
    """List of (node, DOF, coefficient) terms; first row is the dependent."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Multi-Point Constraint (MPC)")
        self.setMinimumWidth(480)

        layout = QVBoxLayout(self)

        form = QFormLayout()
        self.sid_edit = QSpinBox(); self.sid_edit.setRange(1, 99_999_999); self.sid_edit.setValue(1)
        form.addRow("Constraint SID:", self.sid_edit)
        layout.addLayout(form)

        layout.addWidget(QLabel("Terms (first row = dependent):"))
        self.table = QTableWidget(2, 3)
        self.table.setHorizontalHeaderLabels(["Node", "DOF", "Coefficient"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # Two seed rows (dependent + first independent)
        for r in range(2):
            for c, default in enumerate(("", "", "")):
                self.table.setItem(r, c, QTableWidgetItem(default))
        layout.addWidget(self.table)

        row = QHBoxLayout()
        add_btn = QPushButton("+ Add term"); add_btn.clicked.connect(self._add_row)
        del_btn = QPushButton("- Remove last"); del_btn.clicked.connect(self._del_row)
        row.addWidget(add_btn); row.addWidget(del_btn); row.addStretch(1)
        layout.addLayout(row)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._on_accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self._sid = None
        self._terms: list[tuple[int, int, float]] = []

    def _add_row(self):
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c in range(3):
            self.table.setItem(r, c, QTableWidgetItem(""))

    def _del_row(self):
        r = self.table.rowCount()
        if r > 2:
            self.table.removeRow(r - 1)

    def _on_accept(self):
        terms: list[tuple[int, int, float]] = []
        for r in range(self.table.rowCount()):
            try:
                nid = int(self.table.item(r, 0).text().strip())
                dof = int(self.table.item(r, 1).text().strip())
                coef = float(self.table.item(r, 2).text().strip())
            except (ValueError, AttributeError):
                continue
            terms.append((nid, dof, coef))
        if len(terms) < 2:
            QMessageBox.warning(self, "MPC needs at least 2 terms",
                                "Provide a dependent term and at least one independent.")
            return
        if all(abs(t[2]) < 1e-30 for t in terms):
            QMessageBox.warning(self, "Coefficients all zero",
                                "MPC coefficients cannot all be zero.")
            return
        self._sid = int(self.sid_edit.value())
        self._terms = terms
        self.accept()

    @property
    def sid(self):
        return self._sid

    @property
    def terms(self):
        return list(self._terms)


# ---------------------------------------------------------------------------
# C5. RBAR - rigid bar
# ---------------------------------------------------------------------------

class CreateRbarDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Rigid Bar (RBAR)")
        self.setMinimumWidth(380)
        layout = QVBoxLayout(self); form = QFormLayout()
        self.eid = QSpinBox(); self.eid.setRange(1, 99_999_999); self.eid.setValue(1)
        self.ga = QSpinBox(); self.ga.setRange(1, 99_999_999); self.ga.setValue(1)
        self.gb = QSpinBox(); self.gb.setRange(1, 99_999_999); self.gb.setValue(2)
        self.cna = QLineEdit("123456"); self.cnb = QLineEdit("")
        self.cma = QLineEdit(""); self.cmb = QLineEdit("123456")
        for w in (self.cna, self.cnb, self.cma, self.cmb):
            w.setPlaceholderText("DOFs as digits, e.g. 123456")
        form.addRow("EID:", self.eid)
        form.addRow("Grid A (GA):", self.ga)
        form.addRow("Grid B (GB):", self.gb)
        form.addRow("CNA (independent DOFs at A):", self.cna)
        form.addRow("CNB (independent DOFs at B):", self.cnb)
        form.addRow("CMA (dependent DOFs at A):", self.cma)
        form.addRow("CMB (dependent DOFs at B):", self.cmb)
        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def values(self):
        return {
            'eid': int(self.eid.value()),
            'ga': int(self.ga.value()),
            'gb': int(self.gb.value()),
            'cna': self.cna.text().strip(),
            'cnb': self.cnb.text().strip(),
            'cma': self.cma.text().strip(),
            'cmb': self.cmb.text().strip(),
        }


# ---------------------------------------------------------------------------
# C6. RBE1 - general rigid element
# ---------------------------------------------------------------------------

class CreateRbe1Dialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Rigid Element (RBE1)")
        self.setMinimumWidth(440)
        layout = QVBoxLayout(self); form = QFormLayout()
        self.eid = QSpinBox(); self.eid.setRange(1, 99_999_999); self.eid.setValue(1)
        self.indep_nodes = QLineEdit(); self.indep_nodes.setPlaceholderText("e.g. 1, 5, 10")
        self.indep_dofs = QLineEdit("123456")
        self.dep_nodes = QLineEdit(); self.dep_nodes.setPlaceholderText("e.g. 11, 12")
        self.dep_dofs = QLineEdit("123")
        self.alpha = QDoubleSpinBox(); self.alpha.setRange(0.0, 1.0); self.alpha.setDecimals(6); self.alpha.setValue(0.0)
        form.addRow("EID:", self.eid)
        form.addRow("Independent nodes:", self.indep_nodes)
        form.addRow("Independent DOFs:", self.indep_dofs)
        form.addRow("Dependent nodes:", self.dep_nodes)
        form.addRow("Dependent DOFs:", self.dep_dofs)
        form.addRow("Alpha (thermal coef.):", self.alpha)
        layout.addLayout(form)
        info = QLabel(
            "RBE1 lets you spread loads from many independent grid points "
            "to many dependent ones. Use the simpler Spider creator if "
            "you only have one center node."
        )
        info.setWordWrap(True); info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def values(self):
        return {
            'eid': int(self.eid.value()),
            'indep_nodes': _parse_id_list(self.indep_nodes.text()),
            'indep_dofs': self.indep_dofs.text().strip(),
            'dep_nodes': _parse_id_list(self.dep_nodes.text()),
            'dep_dofs': self.dep_dofs.text().strip(),
            'alpha': float(self.alpha.value()),
        }


# ---------------------------------------------------------------------------
# C7. RSPLINE - rigid spline
# ---------------------------------------------------------------------------

class CreateRsplineDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Rigid Spline (RSPLINE)")
        self.setMinimumWidth(420)
        layout = QVBoxLayout(self); form = QFormLayout()
        self.eid = QSpinBox(); self.eid.setRange(1, 99_999_999); self.eid.setValue(1)
        self.nodes_edit = QLineEdit(); self.nodes_edit.setPlaceholderText("e.g. 1, 5, 10, 15")
        self.dofs_edit = QLineEdit("123")
        self.diam = QDoubleSpinBox(); self.diam.setRange(0.0, 1e6); self.diam.setDecimals(6); self.diam.setValue(0.1)
        form.addRow("EID:", self.eid)
        form.addRow("Spline nodes (in order):", self.nodes_edit)
        form.addRow("Released DOFs:", self.dofs_edit)
        form.addRow("Diameter / blend factor:", self.diam)
        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def values(self):
        return {
            'eid': int(self.eid.value()),
            'nodes': _parse_id_list(self.nodes_edit.text()),
            'dofs': self.dofs_edit.text().strip(),
            'diam': float(self.diam.value()),
        }


# ---------------------------------------------------------------------------
# C8. BOLT - preload
# ---------------------------------------------------------------------------

class CreateBoltDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Bolt Preload (BOLT)")
        self.setMinimumWidth(380)
        layout = QVBoxLayout(self); form = QFormLayout()
        self.bid = QSpinBox(); self.bid.setRange(1, 99_999_999); self.bid.setValue(1)
        self.preload = QDoubleSpinBox()
        self.preload.setRange(-1e12, 1e12); self.preload.setDecimals(3); self.preload.setValue(1000.0)
        self.element_edit = QLineEdit()
        self.element_edit.setPlaceholderText("CBAR/CBEAM EID(s) representing the bolt shank")
        form.addRow("BOLT ID:", self.bid)
        form.addRow("Preload force:", self.preload)
        form.addRow("Bolt element EIDs:", self.element_edit)
        layout.addLayout(form)
        info = QLabel(
            "BOLT is solver-version specific (NX vs MSC). This dialog "
            "captures the basic preload; refine for your solver."
        )
        info.setWordWrap(True); info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @property
    def values(self):
        return {
            'bid': int(self.bid.value()),
            'preload': float(self.preload.value()),
            'eids': _parse_id_list(self.element_edit.text()),
        }
