"""Dialog classes for creating and editing loads and constraints."""

import numpy as np

from PySide6 import QtCore
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget, QMessageBox,
    QPushButton, QRadioButton, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)
from PySide6.QtCore import Qt


class CreateLoadDialog(QDialog):
    selection_requested = QtCore.Signal(str)

    def __init__(self, model, parent=None, existing_sid=None):
        super().__init__(parent)
        self.model = model

        self.setWindowTitle("Create Load")
        self.setMinimumWidth(500)

        main_layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._create_nodal_tab(), "Nodal (Force/Moment)")
        self.tabs.addTab(self._create_pressure_tab(), "Elemental (Pressure)")
        self.tabs.addTab(self._create_temp_tab(), "Body (Temperature)")
        self.tabs.addTab(self._create_grav_tab(), "Body (Gravity)")

        main_layout.addWidget(self.tabs)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        # Pre-fill SID if provided (for "add to existing set" workflow)
        if existing_sid:
            self.nodal_sid.setText(str(existing_sid))
            self.pres_sid.setText(str(existing_sid))
            self.temp_sid.setText(str(existing_sid))
            self.grav_sid.setText(str(existing_sid))

    def _populate_from_existing(self, sid):
        """Fills the dialog with data from an existing load set."""
        if sid not in self.model.loads:
            return

        load_cards = self.model.loads[sid]
        if not load_cards:
            return

        # Assumption: all cards in a set are of the same type
        first_card = load_cards[0]
        load_type = first_card.type

        if load_type in ['FORCE', 'MOMENT']:
            self.tabs.setCurrentWidget(self.nodal_tab)
            self.pres_tab.setEnabled(False)
            self.temp_tab.setEnabled(False)
            self.grav_tab.setEnabled(False)

            self.nodal_sid.setText(str(sid))
            self.nodal_sid.setEnabled(False)

            if load_type == 'FORCE':
                self.force_rb.setChecked(True)
            else: # MOMENT
                self.moment_rb.setChecked(True)

            # Get components from first card and scale by magnitude
            mag = first_card.mag
            xyz = np.array(first_card.xyz) * mag
            self.nodal_fx.setText(str(xyz[0]))
            self.nodal_fy.setText(str(xyz[1]))
            self.nodal_fz.setText(str(xyz[2]))

            all_nodes = sorted([card.node_id for card in load_cards])
            self.nodal_node_list.addItems([str(nid) for nid in all_nodes])

        elif load_type == 'PLOAD4':
            self.tabs.setCurrentWidget(self.pres_tab)
            self.nodal_tab.setEnabled(False)
            self.temp_tab.setEnabled(False)
            self.grav_tab.setEnabled(False)

            self.pres_sid.setText(str(sid))
            self.pres_sid.setEnabled(False)
            self.pres_val.setText(str(first_card.pressures[0]))

            all_elems = sorted([card.eids[0] for card in load_cards])
            self.pres_elem_list.addItems([str(eid) for eid in all_elems])

        elif load_type == 'TEMPD':
            self.tabs.setCurrentWidget(self.temp_tab)
            self.nodal_tab.setEnabled(False)
            self.pres_tab.setEnabled(False)
            self.grav_tab.setEnabled(False)

            self.temp_sid.setText(str(sid))
            self.temp_sid.setEnabled(False)
            self.temp_val.setText(str(first_card.temperature))

        elif load_type == 'GRAV':
            self.tabs.setCurrentWidget(self.grav_tab)
            self.nodal_tab.setEnabled(False)
            self.pres_tab.setEnabled(False)
            self.temp_tab.setEnabled(False)

            self.grav_sid.setText(str(sid))
            self.grav_sid.setEnabled(False)
            self.grav_cid.setText(str(first_card.cid))
            self.grav_scale.setText(str(first_card.scale))
            N = first_card.N
            self.grav_nx.setText(str(N[0]))
            self.grav_ny.setText(str(N[1]))
            self.grav_nz.setText(str(N[2]))

    def _create_nodal_tab(self):
        self.nodal_tab = QWidget() # Store reference for easy access
        layout = QVBoxLayout(self.nodal_tab)
        form = QFormLayout()

        self.nodal_sid = QLineEdit("1")
        self.nodal_fx, self.nodal_fy, self.nodal_fz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        self.force_rb = QRadioButton("FORCE", checked=True)
        self.moment_rb = QRadioButton("MOMENT")
        type_layout = QHBoxLayout(); type_layout.addWidget(self.force_rb); type_layout.addWidget(self.moment_rb)

        form.addRow("Load Set ID (SID):", self.nodal_sid)
        form.addRow("Load Type:", type_layout)
        form.addRow("Component Fx / Mx:", self.nodal_fx)
        form.addRow("Component Fy / My:", self.nodal_fy)
        form.addRow("Component Fz / Mz:", self.nodal_fz)
        layout.addLayout(form)

        node_group = QGroupBox("Target Nodes")
        node_layout = QVBoxLayout(node_group)
        select_nodes_btn = QPushButton("Select Nodes...")
        self.nodal_node_list = QListWidget(); self.nodal_node_list.setMaximumHeight(100)
        node_layout.addWidget(select_nodes_btn)
        node_layout.addWidget(self.nodal_node_list)
        layout.addWidget(node_group)

        select_nodes_btn.clicked.connect(lambda: self.selection_requested.emit('Node'))
        return self.nodal_tab

    def _create_pressure_tab(self):
        self.pres_tab = QWidget()
        layout = QVBoxLayout(self.pres_tab)
        form = QFormLayout()

        self.pres_sid = QLineEdit("1")
        self.pres_val = QLineEdit("0.0")
        form.addRow("Load Set ID (SID):", self.pres_sid)
        form.addRow("Pressure (Force/Area):", self.pres_val)
        layout.addLayout(form)

        elem_group = QGroupBox("Target Elements (CQUAD4, CTRIA3)")
        elem_layout = QVBoxLayout(elem_group)
        select_elems_btn = QPushButton("Select Elements...")
        self.pres_elem_list = QListWidget(); self.pres_elem_list.setMaximumHeight(100)
        elem_layout.addWidget(select_elems_btn)
        elem_layout.addWidget(self.pres_elem_list)
        layout.addWidget(elem_group)

        select_elems_btn.clicked.connect(lambda: self.selection_requested.emit('Element'))
        return self.pres_tab

    def _create_temp_tab(self):
        self.temp_tab = QWidget()
        form = QFormLayout(self.temp_tab)
        self.temp_sid = QLineEdit("1")
        self.temp_val = QLineEdit("0.0")
        form.addRow("Temperature Set ID (SID):", self.temp_sid)
        form.addRow("Default Temperature:", self.temp_val)
        return self.temp_tab

    def _create_grav_tab(self):
        self.grav_tab = QWidget()
        layout = QVBoxLayout(self.grav_tab)
        form = QFormLayout()

        self.grav_sid = QLineEdit("1")
        self.grav_cid = QLineEdit("0")
        self.grav_scale = QLineEdit("9810.0")
        self.grav_nx = QLineEdit("0.0")
        self.grav_ny = QLineEdit("0.0")
        self.grav_nz = QLineEdit("-1.0")

        form.addRow("Load Set ID (SID):", self.grav_sid)
        form.addRow("Coordinate System (CID):", self.grav_cid)
        form.addRow("Acceleration Magnitude:", self.grav_scale)
        form.addRow("Direction Nx:", self.grav_nx)
        form.addRow("Direction Ny:", self.grav_ny)
        form.addRow("Direction Nz:", self.grav_nz)
        layout.addLayout(form)

        info_label = QLabel(
            "<b>Note:</b> The gravity vector is defined as N × Scale. "
            "For example, N=(0,0,-1) and Scale=386.1 gives standard "
            "gravity in the -Z direction (in/s² units)."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("font-size: 11px; color: #a6adc8; padding-top: 5px;")
        layout.addWidget(info_label)
        layout.addStretch()
        return self.grav_tab

    def update_selection_list(self, entity_type, ids):
        if entity_type == 'Node':
            self.nodal_node_list.clear()
            self.nodal_node_list.addItems([str(nid) for nid in ids])
        elif entity_type == 'Element':
            self.pres_elem_list.clear()
            self.pres_elem_list.addItems([str(eid) for eid in ids])

    def get_parameters(self):
        try:
            current_tab_widget = self.tabs.currentWidget()
            if current_tab_widget == self.nodal_tab:
                sid = int(self.nodal_sid.text())
                load_type = 'FORCE' if self.force_rb.isChecked() else 'MOMENT'
                components = [float(self.nodal_fx.text()), float(self.nodal_fy.text()), float(self.nodal_fz.text())]
                node_ids = [int(self.nodal_node_list.item(i).text()) for i in range(self.nodal_node_list.count())]
                if not node_ids: QMessageBox.warning(self, "Input Error", "Please select at least one node."); return None
                return {'type': 'nodal', 'sid': sid, 'load_type': load_type, 'components': components, 'node_ids': node_ids}

            elif current_tab_widget == self.pres_tab:
                sid = int(self.pres_sid.text())
                pressure = float(self.pres_val.text())
                element_ids = [int(self.pres_elem_list.item(i).text()) for i in range(self.pres_elem_list.count())]
                if not element_ids: QMessageBox.warning(self, "Input Error", "Please select at least one element."); return None
                return {'type': 'pressure', 'sid': sid, 'pressure': pressure, 'element_ids': element_ids}

            elif current_tab_widget == self.temp_tab:
                sid = int(self.temp_sid.text())
                temperature = float(self.temp_val.text())
                return {'type': 'temperature', 'sid': sid, 'temperature': temperature}

            elif current_tab_widget == self.grav_tab:
                sid = int(self.grav_sid.text())
                cid = int(self.grav_cid.text())
                scale = float(self.grav_scale.text())
                N = [float(self.grav_nx.text()), float(self.grav_ny.text()), float(self.grav_nz.text())]
                if np.linalg.norm(N) == 0:
                    QMessageBox.warning(self, "Input Error", "Direction vector cannot be zero.")
                    return None
                return {'type': 'gravity', 'sid': sid, 'cid': cid, 'scale': scale, 'N': N}

        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please provide valid numerical inputs for all fields.")
            return None


class CreateConstraintDialog(QDialog):
    selection_requested = QtCore.Signal(str)

    def __init__(self, model, parent=None, existing_sid=None):
        super().__init__(parent)
        self.model = model

        title = f"Edit Constraint (SID: {existing_sid})" if existing_sid else "Create Nodal Constraint (SPC)"
        self.setWindowTitle(title)

        main_layout = QVBoxLayout(self)
        form_layout = QFormLayout()

        self.sid_input = QLineEdit("1")
        form_layout.addRow("Constraint Set ID (SID):", self.sid_input)
        main_layout.addLayout(form_layout)

        dof_group = QGroupBox("Degrees of Freedom to Constrain")
        dof_layout = QHBoxLayout(dof_group)
        self.dof_checks = {
            '1': QCheckBox("TX"), '2': QCheckBox("TY"), '3': QCheckBox("TZ"),
            '4': QCheckBox("RX"), '5': QCheckBox("RY"), '6': QCheckBox("RZ")
        }
        for check in self.dof_checks.values():
            dof_layout.addWidget(check)
        main_layout.addWidget(dof_group)

        node_group = QGroupBox("Target Nodes")
        node_layout = QVBoxLayout(node_group)
        select_nodes_btn = QPushButton("Select Nodes...")
        self.node_list_widget = QListWidget()
        self.node_list_widget.setMaximumHeight(100)
        node_layout.addWidget(select_nodes_btn)
        node_layout.addWidget(self.node_list_widget)
        main_layout.addWidget(node_group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        select_nodes_btn.clicked.connect(lambda: self.selection_requested.emit('Node'))
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        if existing_sid:
            self.sid_input.setText(str(existing_sid))
            self.sid_input.setEnabled(False)
            self._populate_from_existing(existing_sid)

    def _populate_from_existing(self, sid):
        """Fills the dialog with data from an existing constraint set."""
        if sid not in self.model.spcs:
            return

        spc_cards = self.model.spcs[sid]
        if not spc_cards:
            return

        # Assume all cards in a set have the same DOFs
        first_card = spc_cards[0]
        dofs = str(first_card.components)
        for dof_char, checkbox in self.dof_checks.items():
            if dof_char in dofs:
                checkbox.setChecked(True)

        # Collect all unique nodes from all cards in the set
        all_nodes = set()
        for card in spc_cards:
            all_nodes.update(card.nodes)

        self.node_list_widget.addItems([str(nid) for nid in sorted(list(all_nodes))])

    def update_selection_list(self, entity_type, ids):
        self.node_list_widget.clear()
        self.node_list_widget.addItems([str(nid) for nid in ids])

    def get_parameters(self):
        try:
            sid = int(self.sid_input.text())
            nodes = [int(self.node_list_widget.item(i).text()) for i in range(self.node_list_widget.count())]
            dof = "".join(key for key, check in self.dof_checks.items() if check.isChecked())

            if not nodes:
                QMessageBox.warning(self, "Input Error", "Please select at least one node.")
                return None
            if not dof:
                QMessageBox.warning(self, "Input Error", "Please select at least one degree of freedom.")
                return None

            return {'sid': sid, 'node_ids': nodes, 'dof': dof}
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please provide a valid integer for the SID.")
            return None


class CreateLoadCombinationDialog(QDialog):
    """Dialog for creating a LOAD combination card."""

    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        self.setWindowTitle("Create LOAD Combination")
        self.setMinimumWidth(500)

        main_layout = QVBoxLayout(self)

        form = QFormLayout()
        self.sid_input = QLineEdit("")
        self.overall_scale = QLineEdit("1.0")
        form.addRow("Combination SID:", self.sid_input)
        form.addRow("Overall Scale (S):", self.overall_scale)
        main_layout.addLayout(form)

        # Suggest next available SID
        existing_sids = set(model.loads.keys())
        if hasattr(model, 'load_combinations'):
            existing_sids |= set(model.load_combinations.keys())
        next_sid = max(existing_sids, default=0) + 1
        self.sid_input.setText(str(next_sid))

        # Table of available load sets
        group = QGroupBox("Component Load Sets")
        group_layout = QVBoxLayout(group)
        info_label = QLabel("Check the load sets to include and set their scale factors.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("font-size: 11px; color: #a6adc8;")
        group_layout.addWidget(info_label)

        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Include", "Load SID", "Scale Factor", "Type"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)

        load_sids = sorted(model.loads.keys())
        self.table.setRowCount(len(load_sids))
        for row, sid in enumerate(load_sids):
            # Include checkbox
            check_item = QTableWidgetItem()
            check_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            check_item.setCheckState(Qt.Unchecked)
            self.table.setItem(row, 0, check_item)

            # Load SID (read-only)
            sid_item = QTableWidgetItem(str(sid))
            sid_item.setFlags(Qt.ItemIsEnabled)
            self.table.setItem(row, 1, sid_item)

            # Scale factor (editable)
            scale_item = QTableWidgetItem("1.0")
            self.table.setItem(row, 2, scale_item)

            # Type summary
            cards = model.loads[sid]
            if cards:
                types = set(c.type for c in cards)
                type_str = ", ".join(sorted(types))
            else:
                type_str = "-"
            type_item = QTableWidgetItem(type_str)
            type_item.setFlags(Qt.ItemIsEnabled)
            self.table.setItem(row, 3, type_item)

        group_layout.addWidget(self.table)
        main_layout.addWidget(group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

    def get_parameters(self):
        try:
            sid = int(self.sid_input.text())
            overall_scale = float(self.overall_scale.text())
        except ValueError:
            QMessageBox.warning(self, "Input Error", "SID must be integer and scale must be numeric.")
            return None

        components = []
        for row in range(self.table.rowCount()):
            check = self.table.item(row, 0)
            if check.checkState() == Qt.Checked:
                try:
                    load_sid = int(self.table.item(row, 1).text())
                    scale = float(self.table.item(row, 2).text())
                    components.append((scale, load_sid))
                except ValueError:
                    QMessageBox.warning(self, "Input Error",
                                        f"Invalid scale factor in row {row + 1}.")
                    return None

        if not components:
            QMessageBox.warning(self, "Input Error",
                                "Please check at least one load set to include.")
            return None

        return {'sid': sid, 'overall_scale': overall_scale, 'components': components}


class SubcaseEditorDialog(QDialog):
    """Dialog for editing the Case Control Deck (subcases)."""

    def __init__(self, model, subcases=None, eigrl_sids=None, parent=None):
        super().__init__(parent)
        self.model = model
        self.setWindowTitle("Subcase Editor (Case Control)")
        self.setMinimumWidth(800)
        self.setMinimumHeight(400)

        main_layout = QVBoxLayout(self)

        info_label = QLabel(
            "Define subcases for the Case Control Deck. Each subcase references "
            "a load set, SPC set, and output requests. For SOL 105 buckling, "
            "set METHOD to the EIGRL SID and STATSUB to the preload subcase ID."
        )
        info_label.setWordWrap(True)
        info_label.setStyleSheet("font-size: 11px; color: #a6adc8; padding-bottom: 5px;")
        main_layout.addWidget(info_label)

        # Table with 9 columns (added METHOD and STATSUB)
        self.table = QTableWidget()
        self.table.setColumnCount(9)
        self.table.setHorizontalHeaderLabels([
            "Subcase ID", "Load SID", "SPC SID",
            "METHOD", "STATSUB",
            "DISP", "STRESS", "FORCE", "STRAIN"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        main_layout.addWidget(self.table)

        # Buttons for add/remove
        btn_layout = QHBoxLayout()
        add_btn = QPushButton("Add Subcase")
        remove_btn = QPushButton("Remove Selected")
        add_btn.clicked.connect(self._add_subcase_row)
        remove_btn.clicked.connect(self._remove_selected_row)
        btn_layout.addWidget(add_btn)
        btn_layout.addWidget(remove_btn)
        btn_layout.addStretch()
        main_layout.addLayout(btn_layout)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

        # Gather available load and SPC SIDs
        self._load_sids = sorted(model.loads.keys())
        if hasattr(model, 'load_combinations'):
            self._load_sids = sorted(set(self._load_sids) | set(model.load_combinations.keys()))
        self._spc_sids = sorted(model.spcs.keys())
        self._eigrl_sids = eigrl_sids or []

        # Populate with existing subcases
        if subcases:
            for sc in subcases:
                self._add_subcase_row(sc)
        else:
            # Start with one empty row
            self._add_subcase_row()

    def _make_output_combo(self, default="NONE"):
        combo = QComboBox()
        combo.addItems(["NONE", "ALL"])
        combo.setCurrentText(default)
        return combo

    def _make_sid_combo(self, sid_list, default=None):
        combo = QComboBox()
        combo.addItem("-")
        for sid in sid_list:
            combo.addItem(str(sid))
        if default is not None:
            combo.setCurrentText(str(default))
        return combo

    def _add_subcase_row(self, data=None):
        row = self.table.rowCount()
        self.table.setRowCount(row + 1)

        if isinstance(data, dict):
            sc_id = data.get('id', row + 1)
            load_sid = data.get('load_sid')
            spc_sid = data.get('spc_sid')
            method_sid = data.get('method_sid')
            statsub_id = data.get('statsub_id')
            disp = data.get('disp', 'NONE')
            stress = data.get('stress', 'NONE')
            force = data.get('force', 'NONE')
            strain = data.get('strain', 'NONE')
        else:
            sc_id = row + 1
            load_sid = spc_sid = method_sid = statsub_id = None
            disp = stress = force = strain = "NONE"

        # Subcase ID
        self.table.setItem(row, 0, QTableWidgetItem(str(sc_id)))
        # Load SID combo
        self.table.setCellWidget(row, 1, self._make_sid_combo(self._load_sids, load_sid))
        # SPC SID combo
        self.table.setCellWidget(row, 2, self._make_sid_combo(self._spc_sids, spc_sid))
        # METHOD combo (EIGRL SIDs)
        self.table.setCellWidget(row, 3, self._make_sid_combo(self._eigrl_sids, method_sid))
        # STATSUB combo (subcase IDs - use generic editable field)
        statsub_combo = QComboBox()
        statsub_combo.setEditable(True)
        statsub_combo.addItem("-")
        # Add existing subcase IDs as options
        for r in range(self.table.rowCount() - 1):
            try:
                existing_id = self.table.item(r, 0).text()
                statsub_combo.addItem(existing_id)
            except (AttributeError, ValueError):
                pass
        if statsub_id is not None:
            statsub_combo.setCurrentText(str(statsub_id))
        self.table.setCellWidget(row, 4, statsub_combo)
        # Output request combos
        self.table.setCellWidget(row, 5, self._make_output_combo(disp))
        self.table.setCellWidget(row, 6, self._make_output_combo(stress))
        self.table.setCellWidget(row, 7, self._make_output_combo(force))
        self.table.setCellWidget(row, 8, self._make_output_combo(strain))

    def _remove_selected_row(self):
        row = self.table.currentRow()
        if row >= 0:
            self.table.removeRow(row)

    def get_subcases(self):
        subcases = []
        for row in range(self.table.rowCount()):
            try:
                sc_id = int(self.table.item(row, 0).text())
            except (ValueError, AttributeError):
                QMessageBox.warning(self, "Input Error",
                                    f"Invalid Subcase ID in row {row + 1}.")
                return None

            load_combo = self.table.cellWidget(row, 1)
            spc_combo = self.table.cellWidget(row, 2)
            method_combo = self.table.cellWidget(row, 3)
            statsub_combo = self.table.cellWidget(row, 4)

            load_text = load_combo.currentText() if load_combo else "-"
            spc_text = spc_combo.currentText() if spc_combo else "-"
            method_text = method_combo.currentText() if method_combo else "-"
            statsub_text = statsub_combo.currentText() if statsub_combo else "-"

            load_sid = int(load_text) if load_text != "-" else None
            spc_sid = int(spc_text) if spc_text != "-" else None
            try:
                method_sid = int(method_text) if method_text != "-" else None
            except ValueError:
                method_sid = None
            try:
                statsub_id = int(statsub_text) if statsub_text != "-" else None
            except ValueError:
                statsub_id = None

            disp = self.table.cellWidget(row, 5).currentText()
            stress = self.table.cellWidget(row, 6).currentText()
            force = self.table.cellWidget(row, 7).currentText()
            strain = self.table.cellWidget(row, 8).currentText()

            subcases.append({
                'id': sc_id, 'load_sid': load_sid, 'spc_sid': spc_sid,
                'method_sid': method_sid, 'statsub_id': statsub_id,
                'disp': disp, 'stress': stress, 'force': force, 'strain': strain,
            })
        return subcases


class CreateEigrlDialog(QDialog):
    """Dialog for creating an EIGRL eigenvalue extraction card."""

    def __init__(self, parent=None, existing=None):
        super().__init__(parent)
        self.setWindowTitle("EIGRL - Eigenvalue Extraction")
        self.setMinimumWidth(400)

        layout = QVBoxLayout(self)
        info = QLabel("Define eigenvalue extraction parameters for buckling/modal analysis.")
        info.setWordWrap(True)
        layout.addWidget(info)

        form = QFormLayout()
        self.sid_input = QLineEdit("100")
        self.sid_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("SID (METHOD reference):", self.sid_input)

        self.v1_input = QLineEdit("")
        self.v1_input.setPlaceholderText("Lower bound (optional)")
        self.v1_input.setValidator(QDoubleValidator())
        form.addRow("V1 (Lower Freq):", self.v1_input)

        self.v2_input = QLineEdit("")
        self.v2_input.setPlaceholderText("Upper bound (optional)")
        self.v2_input.setValidator(QDoubleValidator())
        form.addRow("V2 (Upper Freq):", self.v2_input)

        self.nd_input = QLineEdit("10")
        self.nd_input.setValidator(QIntValidator(1, 10000))
        form.addRow("ND (# Eigenvalues):", self.nd_input)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

        if existing:
            self.sid_input.setText(str(existing.get('sid', 100)))
            if existing.get('v1') is not None:
                self.v1_input.setText(str(existing['v1']))
            if existing.get('v2') is not None:
                self.v2_input.setText(str(existing['v2']))
            self.nd_input.setText(str(existing.get('nd', 10)))

    def get_parameters(self):
        v1 = float(self.v1_input.text()) if self.v1_input.text() else None
        v2 = float(self.v2_input.text()) if self.v2_input.text() else None
        return {
            'sid': int(self.sid_input.text()),
            'v1': v1,
            'v2': v2,
            'nd': int(self.nd_input.text()),
        }


class LoadSetManagerDialog(QDialog):
    """Table-based manager for viewing, adding, and removing load entries in a set."""

    def __init__(self, model, sid, parent=None):
        super().__init__(parent)
        self.model = model
        self.sid = sid
        self.setWindowTitle(f"Load Set Manager - SID {sid}")
        self.setMinimumWidth(750)
        self.setMinimumHeight(400)

        import copy
        self._working_entries = copy.deepcopy(model.loads.get(sid, []))
        self._tempd_entry = copy.deepcopy(model.tempds.get(sid))
        self._tempd_removed = False
        self._modified = False

        main_layout = QVBoxLayout(self)
        total = len(self._working_entries) + (1 if self._tempd_entry else 0)
        main_layout.addWidget(QLabel(f"<b>Load Set SID {sid}</b> - "
                                     f"{total} entries"))

        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["#", "Type", "Target", "Value"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(0,
            QHeaderView.ResizeToContents)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        main_layout.addWidget(self.table)

        btn_layout = QHBoxLayout()
        remove_btn = QPushButton("Remove Selected")
        remove_btn.clicked.connect(self._remove_entry)
        btn_layout.addWidget(remove_btn)
        btn_layout.addStretch()
        main_layout.addLayout(btn_layout)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

        self._populate_table()

    def _populate_table(self):
        total_rows = len(self._working_entries) + (1 if self._tempd_entry else 0)
        self.table.setRowCount(total_rows)
        for i, card in enumerate(self._working_entries):
            self.table.setItem(i, 0, QTableWidgetItem(str(i)))
            card_type = card.type if hasattr(card, 'type') else '?'
            self.table.setItem(i, 1, QTableWidgetItem(card_type))

            target = ""
            value = ""
            try:
                if card_type in ('FORCE', 'FORCE1', 'FORCE2'):
                    target = f"Node {card.node_id}"
                    mag = card.mag
                    xyz = np.array(card.xyz) * mag
                    value = f"[{xyz[0]:.4g}, {xyz[1]:.4g}, {xyz[2]:.4g}]"
                elif card_type in ('MOMENT', 'MOMENT1', 'MOMENT2'):
                    target = f"Node {card.node_id}"
                    mag = card.mag
                    xyz = np.array(card.xyz) * mag
                    value = f"[{xyz[0]:.4g}, {xyz[1]:.4g}, {xyz[2]:.4g}]"
                elif card_type == 'PLOAD4':
                    target = f"Elem {card.eids[0]}"
                    value = f"P={card.pressures[0]:.4g}"
                elif card_type == 'GRAV':
                    target = "Global"
                    N = card.N
                    value = (f"Scale={card.scale:.4g} "
                             f"[{N[0]:.4g}, {N[1]:.4g}, {N[2]:.4g}]")
                elif card_type == 'TEMPD':
                    target = "Global"
                    value = f"T={card.temperature:.4g}"
                else:
                    target = str(card_type)
                    value = "(details unavailable)"
            except Exception:
                target = "?"
                value = "?"

            self.table.setItem(i, 2, QTableWidgetItem(target))
            self.table.setItem(i, 3, QTableWidgetItem(value))

        # Append TEMPD entry at end (from model.tempds, not model.loads)
        if self._tempd_entry:
            i = len(self._working_entries)
            self.table.setItem(i, 0, QTableWidgetItem("T"))
            self.table.setItem(i, 1, QTableWidgetItem("TEMPD"))
            self.table.setItem(i, 2, QTableWidgetItem("Global"))
            try:
                self.table.setItem(i, 3, QTableWidgetItem(
                    f"T={self._tempd_entry.temperature:.4g}"))
            except Exception:
                self.table.setItem(i, 3, QTableWidgetItem("?"))

    def _remove_entry(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.information(self, "No Selection",
                                    "Select a row to remove.")
            return
        reply = QMessageBox.question(
            self, "Confirm",
            f"Remove entry #{row} ({self.table.item(row, 1).text()} - "
            f"{self.table.item(row, 2).text()})?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            if row < len(self._working_entries):
                del self._working_entries[row]
            elif self._tempd_entry and row == len(self._working_entries):
                # Removing the TEMPD entry
                self._tempd_entry = None
                self._tempd_removed = True
            self._modified = True
            self._populate_table()

    def get_modified_entries(self):
        """Return the modified list if changed, or None if unchanged."""
        if self._modified:
            return self._working_entries
        return None

    def tempd_was_removed(self):
        """Return True if the TEMPD entry was removed by the user."""
        return self._tempd_removed
