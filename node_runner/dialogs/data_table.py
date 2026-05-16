"""v5.2.0 item 45: professional Data Table dialog.

Modeless dialog opened from Tools -> Data Table. Live-syncs with the
main viewport selection -- whenever the user picks new nodes or
elements the table updates immediately. The user picks columns from
a dropdown listing available output vectors; copy-to-clipboard and
Excel/CSV export buttons sit at the bottom.

The dialog talks to MainWindow via:

  * ``main_window.selection_changed`` signal (str kind, list ids)
  * ``main_window.op2_results`` dict (for value lookups)
  * ``main_window._current_selection_nids`` / ``_current_selection_eids``
    (initial population when the dialog opens)

It owns no rendering state and never mutates the model.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Optional

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QRadioButton, QButtonGroup, QComboBox, QPushButton, QTableView,
    QHeaderView, QFileDialog, QMessageBox, QApplication,
)

from node_runner.dialogs.result_browser import ResultTableModel


# Columns the user can add. (kind, comp_key, label_template)
# label_template uses {mode_idx} placeholder for eigenvector entries.
DISP_COLUMNS = [
    ('displacement', 'Magnitude',  'Disp Mag'),
    ('displacement', 'T1',         'Disp T1'),
    ('displacement', 'T2',         'Disp T2'),
    ('displacement', 'T3',         'Disp T3'),
    ('displacement', 'R1',         'Disp R1'),
    ('displacement', 'R2',         'Disp R2'),
    ('displacement', 'R3',         'Disp R3'),
]
SPC_COLUMNS = [
    ('spc_forces', 'Magnitude', 'SPC Mag'),
    ('spc_forces', 'T1',        'SPC T1'),
    ('spc_forces', 'T2',        'SPC T2'),
    ('spc_forces', 'T3',        'SPC T3'),
]
STRESS_COLUMNS = [
    ('stress', 'von_mises',      'σ vM'),
    ('stress', 'max_principal',  'σ MaxPrin'),
    ('stress', 'min_principal',  'σ MinPrin'),
    ('stress', 'oxx',            'σ xx'),
    ('stress', 'oyy',            'σ yy'),
    ('stress', 'txy',            'τ xy'),
]


class DataTableDialog(QDialog):
    """Live-syncing post-processing data table."""

    def __init__(self, main_window):
        super().__init__(main_window)
        self.setWindowTitle("Data Table")
        self.setModal(False)
        self.resize(720, 480)
        self._mw = main_window
        # User-configured column list. Each entry is a tuple
        # (kind, component_key, mode_idx, label).
        self._columns: list[tuple[str, str, int, str]] = []
        # Cache of the last selection (kind + ids).
        self._source_kind = 'Node'        # 'Node' or 'Element'

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.setContentsMargins(8, 8, 8, 8)

        # ----- Top: source rows + output set -----
        top_form = QFormLayout()
        top_form.setContentsMargins(0, 0, 0, 0)
        source_row = QHBoxLayout()
        self._radio_nodes = QRadioButton("Selected Nodes")
        self._radio_elems = QRadioButton("Selected Elements")
        self._radio_nodes.setChecked(True)
        self._radio_group = QButtonGroup(self)
        self._radio_group.addButton(self._radio_nodes)
        self._radio_group.addButton(self._radio_elems)
        self._radio_nodes.toggled.connect(self._on_source_changed)
        self._radio_elems.toggled.connect(self._on_source_changed)
        source_row.addWidget(self._radio_nodes)
        source_row.addWidget(self._radio_elems)
        source_row.addStretch(1)
        top_form.addRow("Source rows:", source_row)

        self._output_set_combo = QComboBox()
        self._output_set_combo.currentIndexChanged.connect(
            self._rebuild_table)
        top_form.addRow("Output Set:", self._output_set_combo)
        layout.addLayout(top_form)

        # ----- Column controls -----
        col_row = QHBoxLayout()
        col_row.addWidget(QLabel("Columns:"))
        self._add_column_btn = QPushButton("+ Add Column…")
        self._add_column_btn.clicked.connect(self._on_add_column)
        col_row.addWidget(self._add_column_btn)
        self._remove_column_btn = QPushButton("- Remove Last")
        self._remove_column_btn.clicked.connect(self._on_remove_column)
        col_row.addWidget(self._remove_column_btn)
        col_row.addStretch(1)
        layout.addLayout(col_row)

        # ----- Table view -----
        self._view = QTableView()
        self._view.setSortingEnabled(True)
        self._view.setSelectionBehavior(QTableView.SelectRows)
        self._view.setEditTriggers(QTableView.NoEditTriggers)
        self._view.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self._view.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self._view, 1)

        # ----- Status + buttons -----
        bottom_row = QHBoxLayout()
        self._status_lbl = QLabel("No selection.")
        bottom_row.addWidget(self._status_lbl, 1)
        self._copy_btn = QPushButton("Copy to Clipboard")
        self._copy_btn.clicked.connect(self._on_copy)
        bottom_row.addWidget(self._copy_btn)
        self._export_btn = QPushButton("Export…")
        self._export_btn.clicked.connect(self._on_export)
        bottom_row.addWidget(self._export_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        bottom_row.addWidget(close_btn)
        layout.addLayout(bottom_row)

        # Live-sync.
        try:
            self._mw.selection_changed.connect(self._on_selection_changed)
        except Exception:
            pass

        # Initial population.
        self._populate_output_set_combo()
        self._rebuild_table()

    # ------------------------------------------------------------------
    # Population helpers
    # ------------------------------------------------------------------

    def _populate_output_set_combo(self):
        self._output_set_combo.blockSignals(True)
        try:
            self._output_set_combo.clear()
            bundle = getattr(self._mw, 'op2_results', None)
            if not bundle or not bundle.get('subcases'):
                self._output_set_combo.addItem("(no results)", None)
                self._output_set_combo.setEnabled(False)
                return
            self._output_set_combo.setEnabled(True)
            for sid in sorted(bundle['subcases'].keys()):
                title = (bundle['subcases'][sid].get('title')
                         or bundle['subcases'][sid].get('subtitle')
                         or '')
                label = f"Subcase {sid}: {title}" if title else f"Subcase {sid}"
                self._output_set_combo.addItem(label, int(sid))
        finally:
            self._output_set_combo.blockSignals(False)

    def _available_columns(self) -> list[tuple[str, str, int, str]]:
        """Return the column choices for the active subcase + source kind."""
        bundle = getattr(self._mw, 'op2_results', None)
        if not bundle:
            return []
        sid = self._output_set_combo.currentData()
        if sid is None:
            return []
        sc = bundle.get('subcases', {}).get(int(sid), {})
        cols: list[tuple[str, str, int, str]] = []
        if self._source_kind == 'Node':
            if sc.get('displacements'):
                for kind, comp, label in DISP_COLUMNS:
                    cols.append((kind, comp, -1, label))
            if sc.get('spc_forces'):
                for kind, comp, label in SPC_COLUMNS:
                    cols.append((kind, comp, -1, label))
            eigs = sc.get('eigenvectors') or []
            for mode_idx in range(len(eigs)):
                for kind, comp, label in DISP_COLUMNS[:4]:   # Mag,T1,T2,T3
                    cols.append(('eigenvector', comp, mode_idx,
                                 f"Mode {mode_idx+1} {label}"))
        else:  # Element
            if sc.get('stresses'):
                for kind, comp, label in STRESS_COLUMNS:
                    cols.append((kind, comp, -1, label))
        return cols

    # ------------------------------------------------------------------
    # Live selection
    # ------------------------------------------------------------------

    def _on_selection_changed(self, kind, ids):
        # We only need to react if the changed kind matches our current
        # source filter; toggling the radio also rebuilds.
        if kind == self._source_kind:
            self._rebuild_table()

    def _on_source_changed(self):
        self._source_kind = ('Node'
                             if self._radio_nodes.isChecked()
                             else 'Element')
        self._rebuild_table()

    # ------------------------------------------------------------------
    # Column add / remove
    # ------------------------------------------------------------------

    def _on_add_column(self):
        choices = self._available_columns()
        if not choices:
            QMessageBox.information(
                self, "No columns available",
                "Load results and pick an Output Set first.")
            return
        # Show a small popup with the available columns.
        from PySide6.QtWidgets import QInputDialog
        labels = [c[3] for c in choices]
        text, ok = QInputDialog.getItem(
            self, "Add Column", "Output vector:", labels, 0, False)
        if not ok or text not in labels:
            return
        idx = labels.index(text)
        chosen = choices[idx]
        # Avoid duplicates.
        if any(c == chosen for c in self._columns):
            return
        self._columns.append(chosen)
        self._rebuild_table()

    def _on_remove_column(self):
        if self._columns:
            self._columns.pop()
            self._rebuild_table()

    # ------------------------------------------------------------------
    # Build the table model from current state
    # ------------------------------------------------------------------

    def _selected_ids(self) -> list[int]:
        if self._source_kind == 'Node':
            return list(getattr(self._mw, '_current_selection_nids', []) or [])
        return list(getattr(self._mw, '_current_selection_eids', []) or [])

    def _rebuild_table(self):
        ids = self._selected_ids()
        id_col_name = 'NID' if self._source_kind == 'Node' else 'EID'
        column_names = [id_col_name] + [c[3] for c in self._columns]
        arrays = {id_col_name: np.array(ids, dtype=int)}
        bundle = getattr(self._mw, 'op2_results', None)
        sid = self._output_set_combo.currentData()
        sc = (bundle['subcases'].get(int(sid))
              if (bundle and sid is not None) else None)
        for (kind, comp, mode_idx, label) in self._columns:
            arrays[label] = self._compute_column_values(
                ids, kind, comp, mode_idx, sc)
        model = ResultTableModel(column_names, arrays,
                                 id_column=id_col_name, parent=self)
        self._view.setModel(model)
        self._update_status_label(len(ids))

    def _compute_column_values(self, ids, kind, comp, mode_idx,
                                sc) -> np.ndarray:
        """Extract a numpy array of values for the given column."""
        out = np.zeros(len(ids), dtype=float)
        if sc is None or not ids:
            return out
        if kind == 'displacement':
            d = sc.get('displacements') or {}
            for i, nid in enumerate(ids):
                v = d.get(int(nid))
                out[i] = self._comp_from_vec(v, comp)
        elif kind == 'spc_forces':
            d = sc.get('spc_forces') or {}
            for i, nid in enumerate(ids):
                v = d.get(int(nid))
                out[i] = self._comp_from_vec(v, comp)
        elif kind == 'eigenvector':
            eigs = sc.get('eigenvectors') or []
            if 0 <= mode_idx < len(eigs):
                d = eigs[mode_idx]
                for i, nid in enumerate(ids):
                    v = d.get(int(nid))
                    out[i] = self._comp_from_vec(v, comp)
        elif kind == 'stress':
            d = sc.get('stresses') or {}
            for i, eid in enumerate(ids):
                s = d.get(int(eid)) or {}
                out[i] = float(s.get(comp, 0.0))
        return out

    @staticmethod
    def _comp_from_vec(v, comp: str) -> float:
        if not v:
            return 0.0
        idx_map = {'T1': 0, 'T2': 1, 'T3': 2,
                   'R1': 3, 'R2': 4, 'R3': 5}
        if comp == 'Magnitude' and len(v) >= 3:
            return float((v[0]**2 + v[1]**2 + v[2]**2) ** 0.5)
        i = idx_map.get(comp)
        if i is not None and i < len(v):
            return float(v[i])
        return 0.0

    def _update_status_label(self, n: int):
        kind_word = 'node' if self._source_kind == 'Node' else 'element'
        plural = '' if n == 1 else 's'
        self._status_lbl.setText(
            f"Selection: {n} {kind_word}{plural} "
            f"x {len(self._columns)} column(s)")

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _gather_rows(self) -> tuple[list[str], list[list]]:
        """Pull (headers, rows) from the current table model."""
        model = self._view.model()
        if model is None:
            return [], []
        n_cols = model.columnCount()
        headers = [model.headerData(c, Qt.Horizontal)
                   for c in range(n_cols)]
        rows: list[list] = []
        for r in range(model.rowCount()):
            row = []
            for c in range(n_cols):
                idx = model.index(r, c)
                row.append(model.data(idx, Qt.DisplayRole))
            rows.append(row)
        return headers, rows

    def _on_copy(self):
        headers, rows = self._gather_rows()
        if not rows:
            QMessageBox.information(self, "Empty", "Nothing to copy.")
            return
        # Tab-separated values: Excel pastes these directly.
        buf = io.StringIO()
        buf.write("\t".join(str(h) for h in headers) + "\n")
        for row in rows:
            buf.write("\t".join("" if v is None else str(v)
                                for v in row) + "\n")
        QApplication.clipboard().setText(buf.getvalue())
        try:
            from node_runner.profiling import perf_event
            perf_event('data_table', 'copy',
                       n_rows=len(rows), n_cols=len(headers))
        except Exception:
            pass
        self._status_lbl.setText(
            f"Copied {len(rows)} row(s) x {len(headers)} column(s).")

    def _on_export(self):
        headers, rows = self._gather_rows()
        if not rows:
            QMessageBox.information(self, "Empty", "Nothing to export.")
            return
        filepath, selected_filter = QFileDialog.getSaveFileName(
            self, "Export Data Table", "",
            "Excel (*.xlsx);;CSV (*.csv);;All Files (*)")
        if not filepath:
            return
        path = Path(filepath)
        fmt = 'xlsx' if path.suffix.lower() == '.xlsx' else 'csv'
        if fmt == 'csv':
            try:
                with open(filepath, 'w', newline='', encoding='utf-8') as f:
                    w = csv.writer(f)
                    w.writerow(headers)
                    for row in rows:
                        w.writerow(['' if v is None else v for v in row])
            except Exception as e:
                QMessageBox.critical(self, "Export failed", str(e))
                return
        else:
            try:
                import openpyxl
            except ImportError:
                QMessageBox.critical(
                    self, "openpyxl not installed",
                    "Install openpyxl (pip install openpyxl) to export "
                    "to .xlsx, or use the CSV option instead.")
                return
            try:
                wb = openpyxl.Workbook()
                ws = wb.active
                ws.title = "Data Table"
                ws.append([str(h) for h in headers])
                for row in rows:
                    ws.append([None if v is None else v for v in row])
                wb.save(filepath)
            except Exception as e:
                QMessageBox.critical(self, "Export failed", str(e))
                return
        try:
            from node_runner.profiling import perf_event
            perf_event('data_table', 'export',
                       format=fmt, n_rows=len(rows), n_cols=len(headers))
        except Exception:
            pass
        self._status_lbl.setText(
            f"Exported {len(rows)} rows to {path.name}.")

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def showEvent(self, event):
        # Refresh in case results were loaded while the dialog was closed.
        # v5.2.0 bug-fix (Round 2): use UniqueConnection so re-showing
        # the dialog doesn't end up with N duplicate slot invocations.
        try:
            self._mw.selection_changed.connect(
                self._on_selection_changed, Qt.UniqueConnection)
        except (RuntimeError, TypeError):
            # Already connected -- no harm done.
            pass
        self._populate_output_set_combo()
        self._rebuild_table()
        super().showEvent(event)

    def closeEvent(self, event):
        # v5.2.0 bug-fix (Round 2): leave the signal connected. The
        # dialog instance is reused across open/close; disconnecting
        # here would silently break live-sync on every re-open. Qt
        # tears down the connection automatically when MainWindow
        # destroys this dialog.
        super().closeEvent(event)
