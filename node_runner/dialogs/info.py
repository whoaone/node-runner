import numpy as np

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTextEdit, QDialogButtonBox,
    QLabel, QTableWidget, QTableWidgetItem, QHeaderView, QPushButton,
    QApplication, QFormLayout, QComboBox, QLineEdit, QGroupBox, QCheckBox,
    QTreeWidget, QTreeWidgetItem, QSpinBox, QDoubleSpinBox, QTabWidget, QWidget,
)
from PySide6.QtGui import QColor, QBrush
from PySide6 import QtCore


class InfoDialog(QDialog):
    def __init__(self, title, content, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(400, 300)
        layout = QVBoxLayout(self)
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setText(content)
        text_edit.setFontFamily("monospace")
        layout.addWidget(text_edit)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box)


class NodeElementInfoDialog(QDialog):
    """Tabular info dialog for nodes or elements with labeled columns."""

    def __init__(self, title, data_rows, columns, parent=None, column_tooltips=None):
        """
        Args:
            title: Window title.
            data_rows: list of dicts, each dict keyed by column name.
            columns: list of column header strings in display order.
            column_tooltips: optional dict mapping column name to tooltip string.
        """
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(600, 350)
        self.resize(750, 450)
        self._columns = columns
        self._data_rows = data_rows

        layout = QVBoxLayout(self)

        self._table = QTableWidget(len(data_rows), len(columns))
        self._table.setHorizontalHeaderLabels(columns)
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        self._table.verticalHeader().setVisible(False)
        self._table.setSortingEnabled(False)  # must be OFF during population

        if column_tooltips:
            for col, key in enumerate(columns):
                if key in column_tooltips:
                    header_item = self._table.horizontalHeaderItem(col)
                    if header_item:
                        header_item.setToolTip(column_tooltips[key])

        for row, data in enumerate(data_rows):
            for col, key in enumerate(columns):
                val = data.get(key, '')
                item = QTableWidgetItem()
                if isinstance(val, (int, float)):
                    item.setData(QtCore.Qt.DisplayRole, val)
                else:
                    item.setText(str(val))
                self._table.setItem(row, col, item)

        self._table.setSortingEnabled(True)  # enable after all items populated
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self._table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self._table)

        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(self._copy_table)
        btn_layout.addWidget(copy_btn)
        btn_layout.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        btn_layout.addWidget(ok_btn)
        layout.addLayout(btn_layout)

    def _copy_table(self):
        lines = ['\t'.join(self._columns)]
        for data in self._data_rows:
            lines.append('\t'.join(str(data.get(c, '')) for c in self._columns))
        QApplication.clipboard().setText('\n'.join(lines))


class GroupEditorDialog(QDialog):
    """Edit group contents — view and remove nodes/elements."""

    def __init__(self, group_name, group_data, model, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Group: {group_name}")
        self.setMinimumSize(500, 400)
        self.resize(600, 450)
        self._group_name = group_name
        self._nodes = list(group_data.get("nodes", []))
        self._elements = list(group_data.get("elements", []))
        self._model = model
        self._removed_nodes = []
        self._removed_elements = []
        self._added_nodes = []
        self._added_elements = []

        layout = QVBoxLayout(self)

        tabs = QTabWidget()
        layout.addWidget(tabs)

        # Nodes tab
        nodes_tab = QWidget()
        nodes_layout = QVBoxLayout(nodes_tab)
        self._nodes_table = QTableWidget(len(self._nodes), 4)
        self._nodes_table.setHorizontalHeaderLabels(["NID", "X", "Y", "Z"])
        self._nodes_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._nodes_table.setSelectionBehavior(QTableWidget.SelectRows)
        self._nodes_table.setSelectionMode(QTableWidget.ExtendedSelection)
        self._nodes_table.verticalHeader().setVisible(False)
        self._nodes_table.setSortingEnabled(True)
        for row, nid in enumerate(sorted(self._nodes)):
            item = QTableWidgetItem()
            item.setData(QtCore.Qt.DisplayRole, nid)
            self._nodes_table.setItem(row, 0, item)
            if model and nid in model.nodes:
                node = model.nodes[nid]
                for col, val in enumerate(node.xyz):
                    vi = QTableWidgetItem()
                    vi.setData(QtCore.Qt.DisplayRole, round(float(val), 6))
                    self._nodes_table.setItem(row, col + 1, vi)
        self._nodes_table.horizontalHeader().setStretchLastSection(True)
        nodes_layout.addWidget(self._nodes_table)
        nodes_btn_layout = QHBoxLayout()
        remove_nodes_btn = QPushButton("Remove Selected")
        remove_nodes_btn.clicked.connect(self._remove_selected_nodes)
        add_nodes_btn = QPushButton("Add by ID...")
        add_nodes_btn.clicked.connect(self._add_nodes_by_id)
        nodes_btn_layout.addWidget(remove_nodes_btn)
        nodes_btn_layout.addWidget(add_nodes_btn)
        nodes_btn_layout.addStretch()
        nodes_layout.addLayout(nodes_btn_layout)
        tabs.addTab(nodes_tab, f"Nodes ({len(self._nodes)})")

        # Elements tab
        elems_tab = QWidget()
        elems_layout = QVBoxLayout(elems_tab)
        self._elems_table = QTableWidget(len(self._elements), 3)
        self._elems_table.setHorizontalHeaderLabels(["EID", "Type", "PID"])
        self._elems_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._elems_table.setSelectionBehavior(QTableWidget.SelectRows)
        self._elems_table.setSelectionMode(QTableWidget.ExtendedSelection)
        self._elems_table.verticalHeader().setVisible(False)
        self._elems_table.setSortingEnabled(True)
        all_elements = {}
        if model:
            all_elements = {**model.elements, **model.rigid_elements}
        for row, eid in enumerate(sorted(self._elements)):
            item = QTableWidgetItem()
            item.setData(QtCore.Qt.DisplayRole, eid)
            self._elems_table.setItem(row, 0, item)
            if eid in all_elements:
                elem = all_elements[eid]
                self._elems_table.setItem(row, 1, QTableWidgetItem(elem.type))
                pid_item = QTableWidgetItem()
                pid_val = getattr(elem, 'pid', '')
                if isinstance(pid_val, int):
                    pid_item.setData(QtCore.Qt.DisplayRole, pid_val)
                else:
                    pid_item.setText(str(pid_val))
                self._elems_table.setItem(row, 2, pid_item)
        self._elems_table.horizontalHeader().setStretchLastSection(True)
        elems_layout.addWidget(self._elems_table)
        elems_btn_layout = QHBoxLayout()
        remove_elems_btn = QPushButton("Remove Selected")
        remove_elems_btn.clicked.connect(self._remove_selected_elements)
        add_elems_btn = QPushButton("Add by ID...")
        add_elems_btn.clicked.connect(self._add_elements_by_id)
        elems_btn_layout.addWidget(remove_elems_btn)
        elems_btn_layout.addWidget(add_elems_btn)
        elems_btn_layout.addStretch()
        elems_layout.addLayout(elems_btn_layout)
        tabs.addTab(elems_tab, f"Elements ({len(self._elements)})")
        self._tabs = tabs

        # Dialog buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

    def _remove_selected_nodes(self):
        rows = sorted(set(idx.row() for idx in self._nodes_table.selectedIndexes()), reverse=True)
        for row in rows:
            nid_item = self._nodes_table.item(row, 0)
            if nid_item:
                nid = nid_item.data(QtCore.Qt.DisplayRole)
                if nid in self._nodes:
                    self._nodes.remove(nid)
                    self._removed_nodes.append(nid)
            self._nodes_table.removeRow(row)
        self._tabs.setTabText(0, f"Nodes ({self._nodes_table.rowCount()})")

    def _remove_selected_elements(self):
        rows = sorted(set(idx.row() for idx in self._elems_table.selectedIndexes()), reverse=True)
        for row in rows:
            eid_item = self._elems_table.item(row, 0)
            if eid_item:
                eid = eid_item.data(QtCore.Qt.DisplayRole)
                if eid in self._elements:
                    self._elements.remove(eid)
                    self._removed_elements.append(eid)
            self._elems_table.removeRow(row)
        self._tabs.setTabText(1, f"Elements ({self._elems_table.rowCount()})")

    def _add_nodes_by_id(self):
        from PySide6.QtWidgets import QInputDialog
        text, ok = QInputDialog.getText(self, "Add Nodes", "Enter node IDs (comma-separated):")
        if ok and text.strip():
            new_nids = []
            for part in text.split(','):
                part = part.strip()
                if part.isdigit():
                    nid = int(part)
                    if nid not in self._nodes:
                        self._nodes.append(nid)
                        self._added_nodes.append(nid)
                        new_nids.append(nid)
            for nid in new_nids:
                row = self._nodes_table.rowCount()
                self._nodes_table.insertRow(row)
                item = QTableWidgetItem()
                item.setData(QtCore.Qt.DisplayRole, nid)
                self._nodes_table.setItem(row, 0, item)
                if self._model and nid in self._model.nodes:
                    node = self._model.nodes[nid]
                    for col, val in enumerate(node.xyz):
                        vi = QTableWidgetItem()
                        vi.setData(QtCore.Qt.DisplayRole, round(float(val), 6))
                        self._nodes_table.setItem(row, col + 1, vi)
            self._tabs.setTabText(0, f"Nodes ({self._nodes_table.rowCount()})")

    def _add_elements_by_id(self):
        from PySide6.QtWidgets import QInputDialog
        text, ok = QInputDialog.getText(self, "Add Elements", "Enter element IDs (comma-separated):")
        if ok and text.strip():
            all_elements = {}
            if self._model:
                all_elements = {**self._model.elements, **self._model.rigid_elements}
            new_eids = []
            for part in text.split(','):
                part = part.strip()
                if part.isdigit():
                    eid = int(part)
                    if eid not in self._elements:
                        self._elements.append(eid)
                        self._added_elements.append(eid)
                        new_eids.append(eid)
            for eid in new_eids:
                row = self._elems_table.rowCount()
                self._elems_table.insertRow(row)
                item = QTableWidgetItem()
                item.setData(QtCore.Qt.DisplayRole, eid)
                self._elems_table.setItem(row, 0, item)
                if eid in all_elements:
                    elem = all_elements[eid]
                    self._elems_table.setItem(row, 1, QTableWidgetItem(elem.type))
                    pid_item = QTableWidgetItem()
                    pid_val = getattr(elem, 'pid', '')
                    if isinstance(pid_val, int):
                        pid_item.setData(QtCore.Qt.DisplayRole, pid_val)
                    else:
                        pid_item.setText(str(pid_val))
                    self._elems_table.setItem(row, 2, pid_item)
            self._tabs.setTabText(1, f"Elements ({self._elems_table.rowCount()})")

    def get_changes(self):
        """Returns (added_nodes, added_elements, removed_nodes, removed_elements)."""
        return self._added_nodes, self._added_elements, self._removed_nodes, self._removed_elements


class LenientImportReportDialog(QDialog):
    """Report dialog shown after a lenient BDF import.

    Displays a summary of successfully imported card counts and
    a table of skipped/failed cards with line numbers, card type,
    and error messages.
    """

    def __init__(self, filename, lenient_result, parent=None):
        super().__init__(parent)
        self._filename = filename
        self._result = lenient_result

        self.setWindowTitle("Import Report")
        self.setMinimumSize(720, 480)
        self.resize(780, 560)

        layout = QVBoxLayout(self)

        # ---- Banner ----
        n_skipped = len(lenient_result.skipped)
        total = sum(lenient_result.counts.values())
        if n_skipped:
            banner_text = (
                f"<b>Partial Import: {filename}</b><br>"
                f"{total} cards imported, "
                f"{n_skipped} card{'s' if n_skipped != 1 else ''} skipped."
            )
            banner_style = (
                "background-color: #f9e2af; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;"
            )
        else:
            banner_text = (
                f"<b>Import Successful: {filename}</b><br>"
                f"All {total} cards imported via lenient parsing."
            )
            banner_style = (
                "background-color: #a6e3a1; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;"
            )
        banner = QLabel(banner_text)
        banner.setStyleSheet(banner_style)
        banner.setWordWrap(True)
        layout.addWidget(banner)

        # ---- Summary of imported cards ----
        layout.addWidget(QLabel("<b>Successfully Imported:</b>"))
        sorted_types = sorted(
            lenient_result.counts.items(), key=lambda x: -x[1])

        n_rows = len(sorted_types) + 1  # +1 for total row
        summary_table = QTableWidget(n_rows, 2)
        summary_table.setHorizontalHeaderLabels(["Card Type", "Count"])
        summary_table.setEditTriggers(QTableWidget.NoEditTriggers)
        summary_table.setAlternatingRowColors(True)
        summary_table.verticalHeader().setVisible(False)
        summary_table.setSelectionBehavior(QTableWidget.SelectRows)
        summary_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.Stretch)
        summary_table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeToContents)

        for row, (name, count) in enumerate(sorted_types):
            summary_table.setItem(row, 0, QTableWidgetItem(name))
            count_item = QTableWidgetItem(f"{count:,d}")
            count_item.setTextAlignment(
                QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            summary_table.setItem(row, 1, count_item)

        # Total row (bold)
        total_row = len(sorted_types)
        total_name = QTableWidgetItem("Total")
        total_count = QTableWidgetItem(f"{total:,d}")
        total_count.setTextAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        bold_font = total_name.font()
        bold_font.setBold(True)
        total_name.setFont(bold_font)
        total_count.setFont(bold_font)
        summary_table.setItem(total_row, 0, total_name)
        summary_table.setItem(total_row, 1, total_count)

        summary_table.setMaximumHeight(min(200, 30 + n_rows * 28))
        layout.addWidget(summary_table)

        # ---- Skipped cards table ----
        if n_skipped:
            layout.addWidget(QLabel(
                f"<b>Skipped Cards ({n_skipped}):</b>"))

            table = QTableWidget(n_skipped, 3)
            table.setHorizontalHeaderLabels(
                ["Line", "Card Type", "Error"])
            table.horizontalHeader().setStretchLastSection(True)
            table.horizontalHeader().setSectionResizeMode(
                0, QHeaderView.ResizeToContents)
            table.horizontalHeader().setSectionResizeMode(
                1, QHeaderView.ResizeToContents)
            table.setEditTriggers(QTableWidget.NoEditTriggers)
            table.setAlternatingRowColors(True)
            table.verticalHeader().setVisible(False)

            for row, sc in enumerate(lenient_result.skipped):
                line_str = str(sc.line) if sc.line > 0 else "?"
                line_item = QTableWidgetItem(line_str)
                line_item.setTextAlignment(
                    QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                table.setItem(row, 0, line_item)
                table.setItem(row, 1, QTableWidgetItem(sc.card))
                err_display = sc.error
                if len(err_display) > 150:
                    err_display = err_display[:147] + '...'
                err_item = QTableWidgetItem(err_display)
                err_item.setToolTip(sc.error)
                table.setItem(row, 2, err_item)

            layout.addWidget(table)

        # ---- Buttons ----
        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        btn_layout.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        btn_layout.addWidget(ok_btn)
        layout.addLayout(btn_layout)

    def _copy_report(self):
        """Copy a plain-text version of the report to clipboard."""
        lines = [f"Lenient Import Report: {self._filename}", ""]
        lines.append("Imported:")
        for ctype, count in sorted(self._result.counts.items()):
            lines.append(f"  {ctype}: {count}")
        lines.append(f"  Total: {sum(self._result.counts.values())}")
        lines.append("")
        if self._result.skipped:
            lines.append(
                f"Skipped ({len(self._result.skipped)}):")
            for sc in self._result.skipped:
                lines.append(
                    f"  Line {sc.line}: [{sc.card}] {sc.error}")
                lines.append(f"    Raw: {sc.raw}")
        QApplication.clipboard().setText('\n'.join(lines))


class MaterialImportReportDialog(QDialog):
    """Report dialog shown after importing materials from an external BDF.

    Displays a table of successfully imported materials and a table of
    skipped/failed materials with clear reasons.
    """

    def __init__(self, filename, imported_params, skipped_dup, skipped_type,
                 parse_failures, parent=None):
        """
        Args:
            filename: Name of the imported BDF file.
            imported_params: list of param dicts for successfully imported materials.
            skipped_dup: list of (mid, mat_type) tuples skipped due to MID conflict.
            skipped_type: list of (mid, mat_type) tuples skipped due to unsupported type.
            parse_failures: list of SkippedCard namedtuples for MAT cards that failed parsing.
        """
        super().__init__(parent)
        self._filename = filename
        self._imported = imported_params
        self._skipped_dup = skipped_dup
        self._skipped_type = skipped_type
        self._parse_failures = parse_failures

        self.setWindowTitle("Material Import Report")
        self.setMinimumSize(600, 360)
        self.resize(700, 480)

        layout = QVBoxLayout(self)

        total_imported = len(imported_params)
        n_problems = len(skipped_dup) + len(skipped_type) + len(parse_failures)

        # ---- Banner ----
        if n_problems:
            banner_text = (
                f"<b>Partial Import: {filename}</b><br>"
                f"{total_imported} material{'s' if total_imported != 1 else ''} imported, "
                f"{n_problems} skipped."
            )
            banner_style = (
                "background-color: #f9e2af; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;"
            )
        else:
            banner_text = (
                f"<b>Import Successful: {filename}</b><br>"
                f"{total_imported} material{'s' if total_imported != 1 else ''} imported."
            )
            banner_style = (
                "background-color: #a6e3a1; color: #1e1e2e; "
                "padding: 10px; border-radius: 4px; font-size: 13px;"
            )
        banner = QLabel(banner_text)
        banner.setStyleSheet(banner_style)
        banner.setWordWrap(True)
        layout.addWidget(banner)

        # ---- Imported materials table ----
        if imported_params:
            layout.addWidget(QLabel(f"<b>Imported ({total_imported}):</b>"))
            imp_table = QTableWidget(total_imported, 3)
            imp_table.setHorizontalHeaderLabels(["MID", "Type", "Title"])
            imp_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            imp_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
            imp_table.horizontalHeader().setStretchLastSection(True)
            imp_table.setEditTriggers(QTableWidget.NoEditTriggers)
            imp_table.setAlternatingRowColors(True)
            imp_table.verticalHeader().setVisible(False)
            for row, p in enumerate(imported_params):
                imp_table.setItem(row, 0, QTableWidgetItem(str(p['mid'])))
                imp_table.setItem(row, 1, QTableWidgetItem(p['type']))
                title = p.get('comment', '').lstrip('$ ').strip()
                imp_table.setItem(row, 2, QTableWidgetItem(title))
            layout.addWidget(imp_table)

        # ---- Skipped table ----
        all_skipped = []
        for mid, mtype in skipped_dup:
            all_skipped.append((str(mid), mtype, "MID already exists in model"))
        for mid, mtype in skipped_type:
            all_skipped.append((str(mid), mtype, f"Unsupported type: {mtype}"))
        for sc in parse_failures:
            all_skipped.append((str(sc.line), sc.card, sc.error))

        if all_skipped:
            layout.addWidget(QLabel(f"<b>Skipped ({len(all_skipped)}):</b>"))
            skip_table = QTableWidget(len(all_skipped), 3)
            skip_table.setHorizontalHeaderLabels(["MID / Line", "Type", "Reason"])
            skip_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            skip_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
            skip_table.horizontalHeader().setStretchLastSection(True)
            skip_table.setEditTriggers(QTableWidget.NoEditTriggers)
            skip_table.setAlternatingRowColors(True)
            skip_table.verticalHeader().setVisible(False)
            for row, (id_str, mtype, reason) in enumerate(all_skipped):
                skip_table.setItem(row, 0, QTableWidgetItem(id_str))
                skip_table.setItem(row, 1, QTableWidgetItem(mtype))
                reason_item = QTableWidgetItem(reason if len(reason) <= 120 else reason[:117] + '...')
                reason_item.setForeground(QBrush(QColor('#f38ba8')))
                skip_table.setItem(row, 2, reason_item)
            layout.addWidget(skip_table)

        # ---- Buttons ----
        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        btn_layout.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        btn_layout.addWidget(ok_btn)
        layout.addLayout(btn_layout)

    def _copy_report(self):
        """Copy a plain-text version of the report to clipboard."""
        lines = [f"Material Import Report: {self._filename}", ""]
        if self._imported:
            lines.append(f"Imported ({len(self._imported)}):")
            lines.append(f"  {'MID':<8s} {'Type':<8s} Title")
            lines.append(f"  {'-'*8} {'-'*8} {'-'*20}")
            for p in self._imported:
                title = p.get('comment', '').lstrip('$ ').strip()
                lines.append(f"  {p['mid']:<8d} {p['type']:<8s} {title}")
            lines.append("")
        if self._skipped_dup:
            lines.append(f"Skipped - MID conflict ({len(self._skipped_dup)}):")
            for mid, mtype in self._skipped_dup:
                lines.append(f"  MID {mid} ({mtype})")
        if self._skipped_type:
            lines.append(f"Skipped - unsupported type ({len(self._skipped_type)}):")
            for mid, mtype in self._skipped_type:
                lines.append(f"  MID {mid}: {mtype}")
        if self._parse_failures:
            lines.append(f"Skipped - parse failures ({len(self._parse_failures)}):")
            for sc in self._parse_failures:
                lines.append(f"  Line {sc.line}: [{sc.card}] {sc.error}")
        QApplication.clipboard().setText('\n'.join(lines))


class FreeBodyDiagramDialog(QDialog):
    """Displays grid point force balance for selected nodes.

    Shows applied loads, SPC reactions, and element internal force
    contributions at each node, with a residual (equilibrium) check.
    """

    def __init__(self, op2_results, all_node_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Free Body Diagram")
        self.setMinimumSize(760, 520)
        self.resize(820, 580)
        self._op2 = op2_results
        self._all_nids = sorted(all_node_ids)
        self._arrow_scale = 1.0
        self.request_render_arrows = None

        layout = QVBoxLayout(self)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Subcase:"))
        self.sc_combo = QComboBox()
        for sc_id in sorted(op2_results['subcases'].keys()):
            self.sc_combo.addItem(f"Subcase {sc_id}", sc_id)
        ctrl.addWidget(self.sc_combo)

        ctrl.addWidget(QLabel("Node ID:"))
        self.nid_input = QLineEdit()
        self.nid_input.setPlaceholderText("Enter node ID")
        self.nid_input.setMaximumWidth(120)
        ctrl.addWidget(self.nid_input)
        query_btn = QPushButton("Query")
        query_btn.clicked.connect(self._query_node)
        ctrl.addWidget(query_btn)

        ctrl.addStretch()
        ctrl.addWidget(QLabel("Arrow Scale:"))
        self.scale_input = QLineEdit("1.0")
        self.scale_input.setMaximumWidth(80)
        ctrl.addWidget(self.scale_input)
        layout.addLayout(ctrl)

        filter_layout = QHBoxLayout()
        self.show_applied = QCheckBox("Applied"); self.show_applied.setChecked(True)
        self.show_spc = QCheckBox("SPC Reactions"); self.show_spc.setChecked(True)
        self.show_elements = QCheckBox("Element Forces"); self.show_elements.setChecked(True)
        self.show_applied.toggled.connect(self._refresh_table)
        self.show_spc.toggled.connect(self._refresh_table)
        self.show_elements.toggled.connect(self._refresh_table)
        filter_layout.addWidget(self.show_applied)
        filter_layout.addWidget(self.show_spc)
        filter_layout.addWidget(self.show_elements)
        filter_layout.addStretch()
        layout.addLayout(filter_layout)

        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(
            ["Source", "Fx", "Fy", "Fz", "Mx", "My", "Mz", "Type"])
        self.table.horizontalHeader().setStretchLastSection(True)
        for col in range(1, 7):
            self.table.horizontalHeader().setSectionResizeMode(
                col, QHeaderView.ResizeToContents)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        layout.addWidget(self.table)

        self.residual_label = QLabel("")
        self.residual_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        layout.addWidget(self.residual_label)

        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        btn_layout.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

        self._current_forces = []
        self._current_nid = None

    def _query_node(self):
        try:
            nid = int(self.nid_input.text())
        except ValueError:
            self.residual_label.setText("Enter a valid node ID.")
            return

        sc_id = self.sc_combo.currentData()
        if sc_id is None:
            return

        sc_data = self._op2['subcases'].get(sc_id, {})
        gpf = sc_data.get('grid_point_forces', {})
        forces = gpf.get(nid, [])

        if not forces:
            self.residual_label.setText(
                f"No grid point force data for node {nid} in subcase {sc_id}.")
            self._current_forces = []
            self._current_nid = nid
            self.table.setRowCount(0)
            return

        self._current_forces = forces
        self._current_nid = nid
        self._refresh_table()

        if self.request_render_arrows:
            try:
                scale = float(self.scale_input.text() or 1.0)
            except ValueError:
                scale = 1.0
            self.request_render_arrows(nid, forces, scale,
                                        self.show_applied.isChecked(),
                                        self.show_spc.isChecked(),
                                        self.show_elements.isChecked())

    def _refresh_table(self):
        forces = self._current_forces
        if not forces:
            return

        filtered = []
        for f in forces:
            src = f['source']
            if src == 'APPLIED' and not self.show_applied.isChecked():
                continue
            if src == 'SPC' and not self.show_spc.isChecked():
                continue
            if src.startswith('EID') and not self.show_elements.isChecked():
                continue
            filtered.append(f)

        self.table.setRowCount(len(filtered) + 1)

        totals = np.zeros(6)
        type_colors = {'APPLIED': '#a6e3a1', 'SPC': '#89b4fa'}

        for row, f in enumerate(filtered):
            vals = f['forces']
            totals += np.array(vals)
            src = f['source']
            ftype = ('Applied' if src == 'APPLIED'
                     else ('Reaction' if src == 'SPC' else 'Element'))

            self.table.setItem(row, 0, QTableWidgetItem(src))
            for col in range(6):
                item = QTableWidgetItem(f"{vals[col]:.4g}")
                self.table.setItem(row, col + 1, item)
            type_item = QTableWidgetItem(ftype)
            color = type_colors.get(src, '#f38ba8')
            type_item.setForeground(QBrush(QColor(color)))
            self.table.setItem(row, 7, type_item)

        total_row = len(filtered)
        total_src = QTableWidgetItem("TOTAL")
        total_src.setForeground(QBrush(QColor('#f9e2af')))
        self.table.setItem(total_row, 0, total_src)
        for col in range(6):
            item = QTableWidgetItem(f"{totals[col]:.4g}")
            item.setForeground(QBrush(QColor('#f9e2af')))
            self.table.setItem(total_row, col + 1, item)
        self.table.setItem(total_row, 7, QTableWidgetItem(""))

        residual = np.linalg.norm(totals)
        max_force = max(abs(totals).max(), 1e-30)
        ratio = residual / max_force
        if ratio < 0.01:
            self.residual_label.setText(
                f"Node {self._current_nid}: Equilibrium OK "
                f"(residual = {residual:.4g}, ratio = {ratio:.2e})")
            self.residual_label.setStyleSheet(
                "font-weight: bold; font-size: 13px; color: #a6e3a1;")
        else:
            self.residual_label.setText(
                f"Node {self._current_nid}: Equilibrium WARNING "
                f"(residual = {residual:.4g}, ratio = {ratio:.2e})")
            self.residual_label.setStyleSheet(
                "font-weight: bold; font-size: 13px; color: #f38ba8;")

    def _copy_report(self):
        if not self._current_forces:
            return
        lines = [f"Free Body Diagram: Node {self._current_nid}", ""]
        lines.append(f"{'Source':<16s} {'Fx':>12s} {'Fy':>12s} {'Fz':>12s} "
                     f"{'Mx':>12s} {'My':>12s} {'Mz':>12s}")
        lines.append("-" * 100)
        totals = np.zeros(6)
        for f in self._current_forces:
            v = f['forces']
            totals += np.array(v)
            lines.append(f"{f['source']:<16s} {v[0]:>12.4g} {v[1]:>12.4g} "
                         f"{v[2]:>12.4g} {v[3]:>12.4g} {v[4]:>12.4g} "
                         f"{v[5]:>12.4g}")
        lines.append("-" * 100)
        lines.append(f"{'TOTAL':<16s} {totals[0]:>12.4g} {totals[1]:>12.4g} "
                     f"{totals[2]:>12.4g} {totals[3]:>12.4g} "
                     f"{totals[4]:>12.4g} {totals[5]:>12.4g}")
        lines.append(f"\nResidual: {np.linalg.norm(totals):.4g}")
        QApplication.clipboard().setText('\n'.join(lines))


# ---------------------------------------------------------------------------
# Model Checking Dialogs
# ---------------------------------------------------------------------------

class MassPropertiesReportDialog(QDialog):
    """Mass properties report with per-PID breakdown, CG, and inertia tensor."""

    def __init__(self, mass_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mass Properties Report")
        self.setMinimumSize(600, 500)
        self._mass_data = mass_data
        layout = QVBoxLayout(self)

        # Summary
        summary_group = QGroupBox("Summary")
        summary_layout = QFormLayout(summary_group)
        summary_layout.addRow("Total Mass:", QLabel(f"{mass_data['total_mass']:.6g}"))
        summary_layout.addRow("Structural Mass:", QLabel(f"{mass_data['structural_mass']:.6g}"))
        summary_layout.addRow("CONM2 Mass:", QLabel(f"{mass_data['conm2_mass']:.6g}"))
        cg = mass_data['cg']
        summary_layout.addRow("Center of Gravity:",
                              QLabel(f"({cg[0]:.6g}, {cg[1]:.6g}, {cg[2]:.6g})"))
        layout.addWidget(summary_group)

        # Inertia tensor
        inertia = mass_data.get('inertia', [0]*6)
        inertia_group = QGroupBox("Inertia Tensor (about CG)")
        inertia_layout = QFormLayout(inertia_group)
        labels = ['Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz']
        for label, val in zip(labels, inertia):
            inertia_layout.addRow(f"{label}:", QLabel(f"{val:.6g}"))
        layout.addWidget(inertia_group)

        # Per-PID table
        by_pid = mass_data.get('by_pid', {})
        if by_pid:
            pid_group = QGroupBox("Mass by Property")
            pid_layout = QVBoxLayout(pid_group)
            table = QTableWidget(len(by_pid), 2)
            table.setHorizontalHeaderLabels(["PID", "Mass"])
            table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            for row, (pid, m) in enumerate(sorted(by_pid.items())):
                table.setItem(row, 0, QTableWidgetItem(str(pid)))
                table.setItem(row, 1, QTableWidgetItem(f"{m:.6g}"))
            pid_layout.addWidget(table)
            layout.addWidget(pid_group)

        # Notes
        for note in mass_data.get('notes', []):
            layout.addWidget(QLabel(f"Note: {note}"))

        # Buttons
        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

    def _copy_report(self):
        d = self._mass_data
        lines = ["Mass Properties Report", "=" * 40]
        lines.append(f"Total Mass:      {d['total_mass']:.6g}")
        lines.append(f"Structural Mass: {d['structural_mass']:.6g}")
        lines.append(f"CONM2 Mass:      {d['conm2_mass']:.6g}")
        cg = d['cg']
        lines.append(f"CG:              ({cg[0]:.6g}, {cg[1]:.6g}, {cg[2]:.6g})")
        inertia = d.get('inertia', [0]*6)
        labels = ['Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz']
        lines.append("")
        lines.append("Inertia Tensor (about CG):")
        for label, val in zip(labels, inertia):
            lines.append(f"  {label}: {val:.6g}")
        lines.append("")
        lines.append("Mass by Property:")
        for pid, m in sorted(d.get('by_pid', {}).items()):
            lines.append(f"  PID {pid}: {m:.6g}")
        for note in d.get('notes', []):
            lines.append(f"Note: {note}")
        QApplication.clipboard().setText('\n'.join(lines))


class FreeEdgeReportDialog(QDialog):
    """Free edge report with edge table and Zoom To button."""

    zoom_to_edge = None  # callback set by caller: (nid1, nid2) -> None

    def __init__(self, free_edges, node_coords, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Free Edge Report")
        self.setMinimumSize(500, 400)
        self._edges = free_edges
        self._node_coords = node_coords
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(f"Found {len(free_edges)} free edges"))

        self.table = QTableWidget(len(free_edges), 3)
        self.table.setHorizontalHeaderLabels(["Node 1", "Node 2", "Length"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)

        for row, (n1, n2) in enumerate(free_edges):
            self.table.setItem(row, 0, QTableWidgetItem(str(n1)))
            self.table.setItem(row, 1, QTableWidgetItem(str(n2)))
            try:
                p1 = np.array(node_coords[n1])
                p2 = np.array(node_coords[n2])
                length = float(np.linalg.norm(p2 - p1))
                self.table.setItem(row, 2, QTableWidgetItem(f"{length:.4g}"))
            except (KeyError, TypeError):
                self.table.setItem(row, 2, QTableWidgetItem("N/A"))

        layout.addWidget(self.table)

        btn_layout = QHBoxLayout()
        zoom_btn = QPushButton("Zoom To Selected Edge")
        zoom_btn.clicked.connect(self._zoom_to_edge)
        btn_layout.addWidget(zoom_btn)
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

    def _zoom_to_edge(self):
        row = self.table.currentRow()
        if row < 0 or row >= len(self._edges):
            return
        if self.zoom_to_edge:
            self.zoom_to_edge(self._edges[row][0], self._edges[row][1])

    def _copy_report(self):
        lines = [f"Free Edge Report: {len(self._edges)} edges", ""]
        lines.append(f"{'Node 1':>8s} {'Node 2':>8s} {'Length':>12s}")
        lines.append("-" * 32)
        for n1, n2 in self._edges:
            try:
                p1 = np.array(self._node_coords[n1])
                p2 = np.array(self._node_coords[n2])
                length = float(np.linalg.norm(p2 - p1))
                lines.append(f"{n1:>8d} {n2:>8d} {length:>12.4g}")
            except (KeyError, TypeError):
                lines.append(f"{n1:>8d} {n2:>8d} {'N/A':>12s}")
        QApplication.clipboard().setText('\n'.join(lines))


class QualitySummaryReportDialog(QDialog):
    """Element quality summary with per-metric statistics and failing element selector."""

    select_failing_callback = None  # callback: (list_of_eids) -> None

    def __init__(self, quality_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Element Quality Summary")
        self.setMinimumSize(650, 450)
        self._quality_data = quality_data
        layout = QVBoxLayout(self)

        # Compute statistics per metric
        all_metrics = {}
        for eid, metrics in quality_data.items():
            for name, val in metrics.items():
                all_metrics.setdefault(name, []).append((eid, val))

        self._metric_data = all_metrics

        # Statistics table
        stat_headers = ["Metric", "Count", "Min", "Max", "Mean", "Threshold", "Failing"]
        self.stat_table = QTableWidget(len(all_metrics), len(stat_headers))
        self.stat_table.setHorizontalHeaderLabels(stat_headers)
        self.stat_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # Default thresholds
        defaults = {
            'warp': 10.0, 'aspect': 5.0, 'skew': 60.0,
            'jacobian': 0.5, 'taper': 0.5, 'min_angle': 30.0,
            'max_angle': 150.0,
        }

        for row, (metric, pairs) in enumerate(sorted(all_metrics.items())):
            vals = [v for _, v in pairs]
            self.stat_table.setItem(row, 0, QTableWidgetItem(metric))
            self.stat_table.setItem(row, 1, QTableWidgetItem(str(len(vals))))
            self.stat_table.setItem(row, 2, QTableWidgetItem(f"{min(vals):.4g}"))
            self.stat_table.setItem(row, 3, QTableWidgetItem(f"{max(vals):.4g}"))
            self.stat_table.setItem(row, 4,
                                    QTableWidgetItem(f"{sum(vals)/len(vals):.4g}"))

            thresh_item = QTableWidgetItem(str(defaults.get(metric, "")))
            thresh_item.setFlags(thresh_item.flags() | QtCore.Qt.ItemIsEditable)
            self.stat_table.setItem(row, 5, thresh_item)

            # Count failing
            threshold = defaults.get(metric)
            if threshold is not None:
                if metric in ('min_angle', 'jacobian'):
                    failing = sum(1 for v in vals if v < threshold)
                else:
                    failing = sum(1 for v in vals if v > threshold)
            else:
                failing = 0
            self.stat_table.setItem(row, 6, QTableWidgetItem(str(failing)))

        layout.addWidget(self.stat_table)

        btn_layout = QHBoxLayout()
        refresh_btn = QPushButton("Refresh Failing Counts")
        refresh_btn.clicked.connect(self._refresh_failing)
        btn_layout.addWidget(refresh_btn)
        select_btn = QPushButton("Select Failing Elements")
        select_btn.clicked.connect(self._select_failing)
        btn_layout.addWidget(select_btn)
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(self._copy_report)
        btn_layout.addWidget(copy_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

    def _get_thresholds(self):
        """Read thresholds from the table."""
        thresholds = {}
        for row in range(self.stat_table.rowCount()):
            metric = self.stat_table.item(row, 0).text()
            try:
                val = float(self.stat_table.item(row, 5).text())
                thresholds[metric] = val
            except (ValueError, AttributeError):
                pass
        return thresholds

    def _refresh_failing(self):
        thresholds = self._get_thresholds()
        for row in range(self.stat_table.rowCount()):
            metric = self.stat_table.item(row, 0).text()
            threshold = thresholds.get(metric)
            if threshold is None:
                self.stat_table.setItem(row, 6, QTableWidgetItem("0"))
                continue
            pairs = self._metric_data.get(metric, [])
            vals = [v for _, v in pairs]
            if metric in ('min_angle', 'jacobian'):
                failing = sum(1 for v in vals if v < threshold)
            else:
                failing = sum(1 for v in vals if v > threshold)
            self.stat_table.setItem(row, 6, QTableWidgetItem(str(failing)))

    def _select_failing(self):
        if not self.select_failing_callback:
            return
        thresholds = self._get_thresholds()
        failing_eids = set()
        for metric, threshold in thresholds.items():
            pairs = self._metric_data.get(metric, [])
            if metric in ('min_angle', 'jacobian'):
                failing_eids.update(eid for eid, v in pairs if v < threshold)
            else:
                failing_eids.update(eid for eid, v in pairs if v > threshold)
        self.select_failing_callback(list(failing_eids))

    def _copy_report(self):
        lines = ["Element Quality Summary", "=" * 60]
        lines.append(f"{'Metric':<14s} {'Count':>6s} {'Min':>10s} {'Max':>10s} "
                     f"{'Mean':>10s}")
        lines.append("-" * 60)
        for metric, pairs in sorted(self._metric_data.items()):
            vals = [v for _, v in pairs]
            lines.append(f"{metric:<14s} {len(vals):>6d} {min(vals):>10.4g} "
                         f"{max(vals):>10.4g} {sum(vals)/len(vals):>10.4g}")
        QApplication.clipboard().setText('\n'.join(lines))


class OrphanCheckDialog(QDialog):
    """Shows unused properties, materials, and unreferenced nodes."""

    def __init__(self, orphan_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Orphan Entity Check")
        self.setMinimumSize(450, 350)
        layout = QVBoxLayout(self)

        unused_pids = orphan_data.get('unused_pids', [])
        unused_mids = orphan_data.get('unused_mids', [])
        orphan_nids = orphan_data.get('orphan_nids', [])

        total = len(unused_pids) + len(unused_mids) + len(orphan_nids)
        if total == 0:
            layout.addWidget(QLabel("No orphan entities found. Model is clean."))
        else:
            layout.addWidget(QLabel(f"Found {total} orphan entities:"))

        tree = QTreeWidget()
        tree.setHeaderLabels(["Category", "Count"])
        tree.setColumnCount(2)

        if unused_pids:
            pid_item = QTreeWidgetItem(["Unused Properties", str(len(unused_pids))])
            for pid in unused_pids:
                QTreeWidgetItem(pid_item, [f"PID {pid}", ""])
            tree.addTopLevelItem(pid_item)
            pid_item.setExpanded(True)

        if unused_mids:
            mid_item = QTreeWidgetItem(["Unused Materials", str(len(unused_mids))])
            for mid in unused_mids:
                QTreeWidgetItem(mid_item, [f"MID {mid}", ""])
            tree.addTopLevelItem(mid_item)
            mid_item.setExpanded(True)

        if orphan_nids:
            nid_item = QTreeWidgetItem(["Unreferenced Nodes", str(len(orphan_nids))])
            for nid in orphan_nids[:500]:
                QTreeWidgetItem(nid_item, [f"Node {nid}", ""])
            if len(orphan_nids) > 500:
                QTreeWidgetItem(nid_item,
                                [f"... and {len(orphan_nids) - 500} more", ""])
            tree.addTopLevelItem(nid_item)
            nid_item.setExpanded(True)

        tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        layout.addWidget(tree)

        btn_layout = QHBoxLayout()
        copy_btn = QPushButton("Copy Report")
        copy_btn.clicked.connect(
            lambda: QApplication.clipboard().setText(self._build_report()))
        btn_layout.addWidget(copy_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)

        self._orphan_data = orphan_data

    def _build_report(self):
        d = self._orphan_data
        lines = ["Orphan Entity Check", "=" * 40]
        lines.append(f"Unused Properties ({len(d['unused_pids'])}): "
                     + ", ".join(str(p) for p in d['unused_pids']))
        lines.append(f"Unused Materials ({len(d['unused_mids'])}): "
                     + ", ".join(str(m) for m in d['unused_mids']))
        lines.append(f"Unreferenced Nodes ({len(d['orphan_nids'])}): "
                     + ", ".join(str(n) for n in d['orphan_nids'][:100]))
        if len(d['orphan_nids']) > 100:
            lines.append(f"  ... and {len(d['orphan_nids']) - 100} more")
        return '\n'.join(lines)
