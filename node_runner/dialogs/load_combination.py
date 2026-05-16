"""v3.5.0 (item 2): professional Load Combination edit dialog.

A Nastran LOAD card combines multiple load sets with scale factors:

    LOAD,sid,S,S1,L1,S2,L2,...

where `S` is the overall scale and each `(Si, Li)` pair multiplies
load set `Li` by `Si`. Common tools expose this through "Model > Load >
Combine" with a (scale, SID) row editor.

This dialog mirrors that UI:

    - Combo SID (read-only when editing; pre-filled with max+1 for new)
    - Title (free text; stored on the combo if supported)
    - Overall scale factor (QDoubleSpinBox)
    - Member table: rows of (scale, load-set SID, summary), +/- buttons
    - Live preview of resolved member count
    - Restore / OK / Cancel

Modes: 'create' / 'edit' / 'copy'. Storage roundtrips through
``model.load_combinations[sid]``. After a real ``read_bdf`` parse this
slot is ``list[LOAD]`` (one card object per overlapping SID); in-session
edits write back as a dict ``{scale, scale_factors, load_ids, title}``.
The ``read_combo_payload`` adapter below normalizes both shapes so the
rest of the UI can stay dict-shaped.
"""

from __future__ import annotations

from typing import Optional

from PySide6 import QtCore
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QLineEdit, QDoubleSpinBox, QSpinBox, QFrame, QTableWidget, QHeaderView,
    QTableWidgetItem, QComboBox, QMessageBox, QAbstractItemView,
)


def read_combo_payload(combo_data) -> dict:
    """Normalize ``model.load_combinations[sid]`` to dict shape.

    After ``read_bdf`` the value is ``list[LOAD]`` (one pyNastran LOAD
    card per overlapping SID). In-session edits write back a dict. This
    adapter returns ``{scale, scale_factors, load_ids, title}`` for both.
    """
    if isinstance(combo_data, list) and combo_data:
        L = combo_data[0]
        comment = getattr(L, 'comment', '') or ''
        return {
            'scale': float(getattr(L, 'scale', 1.0)),
            'scale_factors': [float(s) for s in getattr(L, 'scale_factors', [])],
            'load_ids': [int(i) for i in getattr(L, 'load_ids', [])],
            'title': comment.strip(),
        }
    if isinstance(combo_data, dict):
        return {
            'scale': float(combo_data.get('scale', 1.0)),
            'scale_factors': [float(s) for s in combo_data.get('scale_factors', [])],
            'load_ids': [int(i) for i in combo_data.get('load_ids', [])],
            'title': str(combo_data.get('title', '') or ''),
        }
    return {'scale': 1.0, 'scale_factors': [], 'load_ids': [], 'title': ''}


def _summarize_load_set(model, sid) -> str:
    """Compact 'N FORCE, K PLOAD4, ...' summary for a load-set SID."""
    if model is None or sid not in (model.loads or {}):
        return "(unknown)"
    counts = {}
    for L in model.loads[sid]:
        t = getattr(L, 'type', '?')
        counts[t] = counts.get(t, 0) + 1
    return ", ".join(f"{c} {t}" for t, c in counts.items()) or "(empty)"


def _load_set_count(model, sid) -> int:
    if model is None or sid not in (model.loads or {}):
        return 0
    return len(model.loads[sid])


class LoadCombinationDialog(QDialog):
    """Edit / Create / Copy a Nastran LOAD card.

    Result is read via ``result_payload()`` after ``exec()`` returns
    ``Accepted``. Returns a tuple ``(sid, payload_dict)`` where the
    dict has keys ``scale`` (float), ``scale_factors`` (list[float]),
    ``load_ids`` (list[int]). Caller is responsible for writing it
    back to ``model.load_combinations[sid]``.
    """

    MODE_CREATE = 'create'
    MODE_EDIT = 'edit'
    MODE_COPY = 'copy'

    def __init__(self, model, mode: str = MODE_CREATE,
                 combo_sid: Optional[int] = None, parent=None):
        super().__init__(parent)
        self.model = model
        self.mode = mode
        self.original_sid = combo_sid

        all_combo_sids = set((model.load_combinations or {}).keys()) if model else set()
        all_load_sids = sorted((model.loads or {}).keys()) if model else []
        self._available_load_sids = all_load_sids

        # Compute initial state from mode. read_combo_payload handles
        # both list[LOAD] (real parse) and dict (in-session) shapes.
        if mode == self.MODE_EDIT and combo_sid in (model.load_combinations or {}):
            init_data = read_combo_payload(model.load_combinations[combo_sid])
            self._initial_sid = combo_sid
            self._sid_editable = False
            title = f"Edit Load Combination: SID {combo_sid}"
        elif mode == self.MODE_COPY and combo_sid in (model.load_combinations or {}):
            init_data = read_combo_payload(model.load_combinations[combo_sid])
            self._initial_sid = self._next_sid(all_combo_sids | set(all_load_sids))
            self._sid_editable = True
            title = f"Copy Load Combination (from SID {combo_sid})"
        else:
            init_data = read_combo_payload(None)
            self._initial_sid = self._next_sid(all_combo_sids | set(all_load_sids))
            self._sid_editable = True
            title = "New Load Combination"

        self._initial_payload = {
            'scale': init_data['scale'],
            'scale_factors': list(init_data['scale_factors']),
            'load_ids': list(init_data['load_ids']),
        }
        # Title text is not part of the Nastran LOAD card; we keep it
        # only for the dialog's display and a comment annotation.
        self._initial_title = init_data['title']

        self.setWindowTitle(title)
        self.setMinimumWidth(620)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 8)
        outer.setSpacing(8)

        title_label = QLabel(title)
        title_label.setStyleSheet(
            "font-weight: bold; font-size: 13px; color: #f9e2af;")
        outer.addWidget(title_label)

        # --- header form ---
        form = QFormLayout()
        form.setSpacing(6)

        self._sid_input = QSpinBox()
        self._sid_input.setRange(1, 99999999)
        self._sid_input.setValue(self._initial_sid)
        self._sid_input.setEnabled(self._sid_editable)
        if not self._sid_editable:
            self._sid_input.setSuffix("  (read-only when editing)")
        form.addRow("Combo SID:", self._sid_input)

        self._title_input = QLineEdit(self._initial_title)
        form.addRow("Title:", self._title_input)

        self._scale_input = QDoubleSpinBox()
        self._scale_input.setRange(-1e9, 1e9)
        self._scale_input.setDecimals(4)
        self._scale_input.setValue(self._initial_payload['scale'])
        self._scale_input.valueChanged.connect(self._refresh_preview)
        form.addRow("Overall scale (S):", self._scale_input)

        outer.addLayout(form)

        outer.addWidget(QLabel("Member load sets:"))

        # --- member table ---
        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(
            ["Scale", "Load Set SID", "Summary"])
        h = self._table.horizontalHeader()
        h.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        h.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        h.setSectionResizeMode(2, QHeaderView.Stretch)
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.setMinimumHeight(160)
        outer.addWidget(self._table)

        for sf, lid in zip(self._initial_payload['scale_factors'],
                           self._initial_payload['load_ids']):
            self._add_row(sf, lid)

        # +/- buttons
        pm_row = QHBoxLayout()
        pm_row.setSpacing(4)
        pm_row.addStretch(1)
        add_btn = QPushButton("+ Add member")
        add_btn.clicked.connect(lambda: self._add_row())
        rem_btn = QPushButton("- Remove selected")
        rem_btn.clicked.connect(self._remove_selected_rows)
        pm_row.addWidget(add_btn)
        pm_row.addWidget(rem_btn)
        outer.addLayout(pm_row)

        # preview
        self._preview_label = QLabel("")
        self._preview_label.setStyleSheet(
            "color: #94e2d5; font-style: italic;")
        outer.addWidget(self._preview_label)
        self._refresh_preview()

        # separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        outer.addWidget(sep)

        # bottom buttons
        bot = QHBoxLayout()
        self._restore_btn = QPushButton("Restore")
        self._restore_btn.setToolTip("Reset all fields to the on-open state.")
        self._restore_btn.clicked.connect(self._restore)
        bot.addWidget(self._restore_btn)
        bot.addStretch(1)
        ok = QPushButton("OK")
        ok.setStyleSheet("font-weight: bold;")
        ok.setMinimumWidth(90)
        ok.clicked.connect(self._on_accept)
        cancel = QPushButton("Cancel")
        cancel.setMinimumWidth(80)
        cancel.clicked.connect(self.reject)
        bot.addWidget(ok)
        bot.addWidget(cancel)
        outer.addLayout(bot)

    # ------------------------------------------------------------------
    # Row management
    # ------------------------------------------------------------------

    def _add_row(self, scale_factor: float = 1.0,
                 load_id: Optional[int] = None):
        row = self._table.rowCount()
        self._table.insertRow(row)

        sb = QDoubleSpinBox()
        sb.setRange(-1e9, 1e9)
        sb.setDecimals(4)
        sb.setValue(float(scale_factor))
        sb.valueChanged.connect(self._refresh_preview)
        self._table.setCellWidget(row, 0, sb)

        cb = QComboBox()
        for sid in self._available_load_sids:
            cb.addItem(str(sid), int(sid))
        if load_id is not None:
            idx = cb.findData(int(load_id))
            if idx >= 0:
                cb.setCurrentIndex(idx)
            else:
                # The combo doesn't include this SID. Add it.
                cb.addItem(str(load_id), int(load_id))
                cb.setCurrentIndex(cb.count() - 1)
        cb.currentIndexChanged.connect(
            lambda _r=row: self._refresh_summary_cell(_r))
        cb.currentIndexChanged.connect(self._refresh_preview)
        self._table.setCellWidget(row, 1, cb)

        self._refresh_summary_cell(row)
        self._refresh_preview()

    def _refresh_summary_cell(self, row: int):
        cb = self._table.cellWidget(row, 1)
        if cb is None:
            return
        sid = cb.currentData()
        summary = _summarize_load_set(self.model, sid)
        item = QTableWidgetItem(summary)
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
        self._table.setItem(row, 2, item)

    def _remove_selected_rows(self):
        rows = sorted(
            {idx.row() for idx in self._table.selectedIndexes()},
            reverse=True)
        for r in rows:
            self._table.removeRow(r)
        self._refresh_preview()

    def _refresh_preview(self):
        # _add_row fires during __init__ before _preview_label exists;
        # skip silently in that case (the final init step also calls
        # _refresh_preview after the label is in place).
        if not hasattr(self, '_preview_label'):
            return
        n_members = self._table.rowCount()
        total = 0
        for r in range(n_members):
            cb = self._table.cellWidget(r, 1)
            if cb is None:
                continue
            sid = cb.currentData()
            total += _load_set_count(self.model, sid)
        self._preview_label.setText(
            f"Resolves to {n_members} member set(s), "
            f"totalling {total:,} load entries")

    # ------------------------------------------------------------------
    # OK / Restore handlers
    # ------------------------------------------------------------------

    def _restore(self):
        self._sid_input.setValue(self._initial_sid)
        self._title_input.setText(self._initial_title)
        self._scale_input.setValue(self._initial_payload['scale'])
        # Clear table, re-populate.
        while self._table.rowCount():
            self._table.removeRow(0)
        for sf, lid in zip(self._initial_payload['scale_factors'],
                           self._initial_payload['load_ids']):
            self._add_row(sf, lid)
        self._refresh_preview()

    def _on_accept(self):
        # Validate: SID not colliding with another combo (or load) when
        # creating / copying.
        new_sid = int(self._sid_input.value())
        if self._sid_editable:
            existing_combos = set(
                (self.model.load_combinations or {}).keys()) if self.model else set()
            existing_loads = set((self.model.loads or {}).keys()) if self.model else set()
            collide = existing_combos | existing_loads
            if new_sid in collide:
                QMessageBox.warning(
                    self, "SID conflict",
                    f"SID {new_sid} is already used by another load set "
                    f"or combination. Please choose a different SID.")
                return
        if self._table.rowCount() == 0:
            confirm = QMessageBox.question(
                self, "Empty combination",
                "This combination has no member load sets. "
                "Save it anyway?")
            if confirm != QMessageBox.Yes:
                return
        self.accept()

    # ------------------------------------------------------------------
    # Result API
    # ------------------------------------------------------------------

    def result_payload(self):
        """Return (sid, payload_dict) describing the dialog's state.

        payload_dict matches the pyNastran ``model.load_combinations``
        shape with an extra ``title`` key the caller can persist
        separately (e.g. as a comment)."""
        sid = int(self._sid_input.value())
        scale = float(self._scale_input.value())
        scale_factors = []
        load_ids = []
        for r in range(self._table.rowCount()):
            sb = self._table.cellWidget(r, 0)
            cb = self._table.cellWidget(r, 1)
            if sb is None or cb is None:
                continue
            scale_factors.append(float(sb.value()))
            load_ids.append(int(cb.currentData()))
        return sid, {
            'scale': scale,
            'scale_factors': scale_factors,
            'load_ids': load_ids,
            'title': self._title_input.text().strip(),
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _next_sid(used: set) -> int:
        return (max(used) + 1) if used else 1
