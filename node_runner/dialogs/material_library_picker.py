"""v5.0.0 item 8: MMPDS-traceable material library picker.

Opened from ``CreateMaterialDialog`` via the new "Load..." button. Lets
the user filter and select a preset material from
``node_runner.material_library.MATERIALS_US``. On accept, returns the
selected material dict; the caller is expected to (a) show the unit
confirmation banner and (b) populate the MAT1 / MAT8 form fields.

Convention: the source units are surfaced both in the table header
AND in the confirmation banner, because the app itself is unitless and
US aerospace work almost always uses psi / lbm / in / degF rather than
the SI defaults that many tools assume.
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtGui import QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit,
    QTableView, QHeaderView, QAbstractItemView, QFrame, QMessageBox,
    QSortFilterProxyModel,
)

from node_runner.material_library import (
    MATERIALS_US, SOURCE_UNITS_TEXT, MATERIAL_LIBRARY_VERSION,
)


_COLUMNS = ("Name", "Type", "Category", "Basis", "Source")


class MaterialLibraryPicker(QDialog):
    """Modal picker over MATERIALS_US."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Load Material from Library")
        self.setModal(True)
        self.setMinimumSize(720, 460)

        self._selected: Optional[dict] = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 10)
        outer.setSpacing(6)

        intro = QLabel(
            "Pick a preset material. On OK, you'll be asked to confirm "
            "the source units match your model before the fields are "
            "populated. Values are inserted exactly as listed - no "
            "conversion is done."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        outer.addWidget(intro)

        units_label = QLabel(f"<b>Source units:</b> {SOURCE_UNITS_TEXT}")
        units_label.setStyleSheet(
            "color: #f9e2af; font-size: 11px; padding: 4px 0;"
        )
        outer.addWidget(units_label)

        filter_row = QHBoxLayout()
        filter_row.setSpacing(6)
        filter_row.addWidget(QLabel("Filter:"))
        self._filter_edit = QLineEdit()
        self._filter_edit.setPlaceholderText(
            "Type to filter by name / category / source...")
        filter_row.addWidget(self._filter_edit, 1)
        outer.addLayout(filter_row)

        # Model + proxy
        self._model = QStandardItemModel(0, len(_COLUMNS), self)
        self._model.setHorizontalHeaderLabels(list(_COLUMNS))
        for m in MATERIALS_US:
            row = [
                QStandardItem(str(m.get("name", ""))),
                QStandardItem(str(m.get("type", ""))),
                QStandardItem(str(m.get("category", ""))),
                QStandardItem(str(m.get("basis", ""))),
                QStandardItem(str(m.get("source", ""))),
            ]
            for it in row:
                it.setEditable(False)
            # Stash the material dict on the first cell for retrieval.
            row[0].setData(m, Qt.UserRole + 1)
            self._model.appendRow(row)

        self._proxy = QSortFilterProxyModel(self)
        self._proxy.setSourceModel(self._model)
        self._proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self._proxy.setFilterKeyColumn(-1)  # all columns
        self._filter_edit.textChanged.connect(self._proxy.setFilterFixedString)

        self._view = QTableView()
        self._view.setModel(self._proxy)
        self._view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._view.setSelectionMode(QAbstractItemView.SingleSelection)
        self._view.setSortingEnabled(True)
        self._view.setAlternatingRowColors(True)
        self._view.horizontalHeader().setStretchLastSection(True)
        self._view.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self._view.verticalHeader().setVisible(False)
        self._view.doubleClicked.connect(self._on_accept)
        self._view.selectionModel().selectionChanged.connect(
            self._on_selection_changed)
        outer.addWidget(self._view, 1)

        # Selection preview block.
        self._preview_label = QLabel("Select a material to preview values.")
        self._preview_label.setStyleSheet(
            "color: #94e2d5; font-family: Consolas, monospace; "
            "font-size: 11px; padding: 6px 4px;")
        self._preview_label.setWordWrap(True)
        outer.addWidget(self._preview_label)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        outer.addWidget(sep)

        button_row = QHBoxLayout()
        button_row.addWidget(QLabel(
            f"<i>Library v{MATERIAL_LIBRARY_VERSION} - "
            f"{len(MATERIALS_US)} entries</i>"))
        button_row.addStretch(1)
        self._ok_btn = QPushButton("Load…")
        self._ok_btn.setDefault(True)
        self._ok_btn.setEnabled(False)
        self._ok_btn.clicked.connect(self._on_accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_row.addWidget(self._ok_btn)
        button_row.addWidget(cancel_btn)
        outer.addLayout(button_row)

    # ------------------------------------------------------------------
    # Behavior
    # ------------------------------------------------------------------

    def _current_material(self) -> Optional[dict]:
        idx = self._view.currentIndex()
        if not idx.isValid():
            return None
        src_idx = self._proxy.mapToSource(idx)
        first_cell = self._model.item(src_idx.row(), 0)
        if first_cell is None:
            return None
        return first_cell.data(Qt.UserRole + 1)

    def _on_selection_changed(self, *_a):
        mat = self._current_material()
        if not mat:
            self._ok_btn.setEnabled(False)
            self._preview_label.setText(
                "Select a material to preview values.")
            return
        self._ok_btn.setEnabled(True)
        # Build a compact preview line.
        if mat.get("type") == "MAT1":
            text = (
                f"<b>{mat['name']}</b><br/>"
                f"E={mat['E']:.3g}, G={mat['G']:.3g}, nu={mat['nu']:.3g}, "
                f"rho={mat['rho']:.4g}, "
                f"alpha={mat.get('alpha', 0):.3g} 1/degF, "
                f"Tref={mat.get('Tref', 70.0):.1f} degF<br/>"
                f"<i>{mat.get('source', '')}</i>"
            )
        elif mat.get("type") == "MAT8":
            text = (
                f"<b>{mat['name']}</b><br/>"
                f"E1={mat['E1']:.3g}, E2={mat['E2']:.3g}, "
                f"nu12={mat['nu12']:.3g}, G12={mat['G12']:.3g}, "
                f"rho={mat['rho']:.4g}<br/>"
                f"<i>{mat.get('source', '')}</i>"
            )
        else:
            text = f"<b>{mat['name']}</b><br/><i>{mat.get('source', '')}</i>"
        self._preview_label.setText(text)

    def _on_accept(self, *_a):
        mat = self._current_material()
        if not mat:
            return
        # Unit-safety banner: explicit confirmation before fields populate.
        banner = (
            f"Loading MMPDS-traceable values for:\n\n"
            f"    {mat['name']}\n\n"
            f"Source units:\n"
            f"    • Stress  : psi\n"
            f"    • Density : lbm/in^3\n"
            f"    • Thermal : 1/degF  (Tref = {mat.get('Tref', 70.0):.0f} degF)\n\n"
            f"Values insert EXACTLY as listed - no conversion is done.\n"
            f"Make sure this matches your model's unit system before continuing."
        )
        ret = QMessageBox.question(
            self, "Confirm material values", banner,
            QMessageBox.Ok | QMessageBox.Cancel,
            QMessageBox.Ok,
        )
        if ret != QMessageBox.Ok:
            try:
                from node_runner.profiling import perf_event
                perf_event('materials', 'picker_cancel_unit_banner',
                           name=mat['name'])
            except Exception:
                pass
            return
        try:
            from node_runner.profiling import perf_event
            perf_event('materials', 'picker_choose',
                       name=mat['name'], mat_type=mat.get('type', ''))
        except Exception:
            pass
        self._selected = mat
        self.accept()

    def selected_material(self) -> Optional[dict]:
        return self._selected
