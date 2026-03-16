import copy
import numpy as np
from PySide6 import QtCore
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QComboBox, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QDialogButtonBox,
)

from node_runner.utils import get_entity_title_from_comment
from node_runner.dialogs.creators import CreateMaterialDialog, CreatePropertyDialog
from node_runner.dialogs.widgets import OrientationWidget
from node_runner.commands import DeleteMaterialCommand, DeletePropertyCommand


class MaterialEditorDialog(QDialog):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        self.parent_window = parent
        self.setWindowTitle("Material Editor")
        self.setMinimumSize(500, 400)

        layout = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["ID", "Title", "Type"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers) # Disable in-place editing
        self.table.setSelectionBehavior(QTableWidget.SelectRows)

        self.populate_table()
        layout.addWidget(self.table)

        button_layout = QHBoxLayout()
        delete_btn = QPushButton("Delete Selected")
        delete_btn.clicked.connect(self._delete_selected_material)
        button_layout.addWidget(delete_btn)
        button_layout.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        button_layout.addWidget(close_btn)
        layout.addLayout(button_layout)

        self.table.itemDoubleClicked.connect(self.edit_material)

    def populate_table(self):
        self.table.setRowCount(len(self.model.materials))
        for row, (mid, mat) in enumerate(sorted(self.model.materials.items())):
            title = get_entity_title_from_comment(mat.comment, mat.type, mid)
            self.table.setItem(row, 0, QTableWidgetItem(str(mid)))
            self.table.setItem(row, 1, QTableWidgetItem(title))
            self.table.setItem(row, 2, QTableWidgetItem(mat.type))

    def edit_material(self, item):
        mid = int(self.table.item(item.row(), 0).text())
        material_to_edit = self.model.materials.get(mid)
        if not material_to_edit:
            return

        dialog = CreateMaterialDialog(mid, self.parent_window, material_to_edit)
        if dialog.exec():
            if params := dialog.get_parameters():
                self.parent_window._edit_material(mid, params)
                self.populate_table()
                self.parent_window._populate_tree()

    def _delete_selected_material(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.information(self, "No Selection", "Please select a material to delete.")
            return
        mid = int(self.table.item(row, 0).text())
        using_pids = [pid for pid, prop in self.model.properties.items() if hasattr(prop, 'mid') and prop.mid == mid]
        msg = f"Are you sure you want to delete Material {mid}?"
        if using_pids:
            msg += f"\n\nWarning: {len(using_pids)} property/properties reference this material."
        reply = QMessageBox.question(self, "Confirm Deletion", msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            cmd = DeleteMaterialCommand(mid)
            self.parent_window.command_manager.execute(cmd, self.model)
            self.populate_table()
            self.parent_window._populate_tree()
            self.parent_window._update_status(f"Deleted Material {mid}.")


class PropertyEditorDialog(QDialog):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        self.parent_window = parent
        self.setWindowTitle("Property Editor")
        self.setMinimumSize(500, 400)

        layout = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["ID", "Title", "Type"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)

        self.populate_table()
        layout.addWidget(self.table)

        button_layout = QHBoxLayout()
        delete_btn = QPushButton("Delete Selected")
        delete_btn.clicked.connect(self._delete_selected_property)
        button_layout.addWidget(delete_btn)
        button_layout.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.reject)
        button_layout.addWidget(close_btn)
        layout.addLayout(button_layout)

        self.table.itemDoubleClicked.connect(self.edit_property)

    def populate_table(self):
        self.table.setRowCount(len(self.model.properties))
        for row, (pid, prop) in enumerate(sorted(self.model.properties.items())):
            title = get_entity_title_from_comment(prop.comment, prop.type, pid)
            self.table.setItem(row, 0, QTableWidgetItem(str(pid)))
            self.table.setItem(row, 1, QTableWidgetItem(title))
            self.table.setItem(row, 2, QTableWidgetItem(prop.type))

    def edit_property(self, item):
        pid = int(self.table.item(item.row(), 0).text())
        prop_to_edit = self.model.properties.get(pid)
        if not prop_to_edit:
            return

        dialog = CreatePropertyDialog(pid, self.model, self.parent_window, prop_to_edit)
        if dialog.exec():
            if params := dialog.get_parameters():
                self.parent_window._edit_property(pid, params)
                self.populate_table()
                self.parent_window._populate_tree()
                self.parent_window._update_plot_visibility()

    def _delete_selected_property(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.information(self, "No Selection", "Please select a property to delete.")
            return
        pid = int(self.table.item(row, 0).text())
        using_eids = [eid for eid, elem in self.model.elements.items() if hasattr(elem, 'pid') and elem.pid == pid]
        msg = f"Are you sure you want to delete Property {pid}?"
        if using_eids:
            msg += f"\n\nWarning: {len(using_eids)} element(s) reference this property."
        reply = QMessageBox.question(self, "Confirm Deletion", msg, QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            cmd = DeletePropertyCommand(pid)
            self.parent_window.command_manager.execute(cmd, self.model)
            self.populate_table()
            self.parent_window._populate_tree()
            self.parent_window._update_viewer(self.parent_window.current_generator, reset_camera=False)
            self.parent_window._update_status(f"Deleted Property {pid}.")


    def accept(self):
        try:
            for row in range(self.table.rowCount()):
                pid = int(self.table.item(row, 0).text())
                prop = self.model.properties[pid]
                prop.comment = f"$ {self.table.item(row, 1).text()}"
                if hasattr(prop, 'mid'):
                    prop.mid = int(self.table.cellWidget(row, 3).currentText())
            super().accept()
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Invalid input for property data: {e}")


class ElementEditorDialog(QDialog):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent
        self.setWindowTitle("Edit Element")
        self.setMinimumWidth(450)
        self.model = model
        self.current_eid = None
        self.node_inputs = []
        self.orientation_widget = None

        main_layout = QVBoxLayout(self)

        load_group = QGroupBox("Load Element")
        load_layout = QHBoxLayout(load_group)
        self.eid_input = QLineEdit()
        self.eid_input.setPlaceholderText("Enter Element ID (EID)")
        load_button = QPushButton("Load")
        load_layout.addWidget(QLabel("EID:"))
        load_layout.addWidget(self.eid_input)
        load_layout.addWidget(load_button)
        main_layout.addWidget(load_group)

        self.props_group = QGroupBox("Element Properties")
        self.props_group.setVisible(False)
        self.props_layout = QFormLayout(self.props_group)
        main_layout.addWidget(self.props_group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

        load_button.clicked.connect(self._load_element_data)

    def _load_element_data(self):
        try:
            eid = int(self.eid_input.text())
            if eid not in self.model.elements: raise KeyError
            self.current_eid = eid
            element = self.model.elements[eid]
            self._populate_ui_for_element(element)
            self.props_group.setVisible(True)
        except (ValueError, KeyError):
            QMessageBox.warning(self, "Error", f"Element ID '{self.eid_input.text()}' not found.")
            self.props_group.setVisible(False)
            self.current_eid = None

    def _clear_layout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if widget := child.widget():
                widget.deleteLater()

    def _populate_ui_for_element(self, element):
        self._clear_layout(self.props_layout)
        self.node_inputs = []
        self.orientation_widget = None

        self.props_layout.addRow(QLabel(f"<b>Type: {element.type}</b>"))

        self.pid_combo = QComboBox()
        prop_ids = [str(pid) for pid in self.model.properties.keys()]
        self.pid_combo.addItems(prop_ids)
        if hasattr(element, 'pid') and (current_pid := str(element.pid)) in prop_ids:
            self.pid_combo.setCurrentText(current_pid)
        self.props_layout.addRow("Property ID (PID):", self.pid_combo)

        if hasattr(element, 'nodes'):
            for i, nid in enumerate(element.nodes):
                node_edit = QLineEdit(str(nid))
                self.node_inputs.append(node_edit)
                self.props_layout.addRow(f"Node {i+1} ID:", node_edit)

        # --- MODIFIED: Add CBEAM and CBAR to this check ---
        if element.type in ['CBUSH', 'CBEAM', 'CBAR']:
            self.orientation_widget = OrientationWidget(self.model, self)
            self.orientation_widget.set_orientation(element)
            self.props_layout.addRow(self.orientation_widget)
            self.orientation_widget.pick_orientation_node_requested.connect(self._pick_orientation_node_for_edit)

        # --- Beam offset fields for CBEAM ---
        self.offset_fields = None
        if element.type == 'CBEAM':
            offset_group = QGroupBox("Beam Offsets")
            offset_form = QFormLayout(offset_group)

            self.offt_combo = QComboBox()
            offt_options = ["GGG", "BGG", "GBG", "GGB", "BGG", "BBG", "BGB", "BBB"]
            self.offt_combo.addItems(offt_options)
            offt_val = getattr(element, 'offt', 'GGG') or 'GGG'
            if offt_val in offt_options:
                self.offt_combo.setCurrentText(offt_val)
            offset_form.addRow("Offset Type (OFFT):", self.offt_combo)

            wa = getattr(element, 'wa', None)
            if wa is None:
                wa = [0.0, 0.0, 0.0]
            wb = getattr(element, 'wb', None)
            if wb is None:
                wb = [0.0, 0.0, 0.0]

            self.wa_x = QLineEdit(str(float(wa[0])))
            self.wa_y = QLineEdit(str(float(wa[1])))
            self.wa_z = QLineEdit(str(float(wa[2])))
            self.wb_x = QLineEdit(str(float(wb[0])))
            self.wb_y = QLineEdit(str(float(wb[1])))
            self.wb_z = QLineEdit(str(float(wb[2])))

            offset_form.addRow("WA x:", self.wa_x)
            offset_form.addRow("WA y:", self.wa_y)
            offset_form.addRow("WA z:", self.wa_z)
            offset_form.addRow("WB x:", self.wb_x)
            offset_form.addRow("WB y:", self.wb_y)
            offset_form.addRow("WB z:", self.wb_z)

            self.props_layout.addRow(offset_group)
            self.offset_fields = True

    def _pick_orientation_node_for_edit(self):
        if self.parent_window and self.orientation_widget:
            self.hide()
            self.parent_window._activate_single_node_picker(self.orientation_widget.orient_node_id.setText, calling_dialog=self)

    def accept(self):
        if self.current_eid is None:
            super().reject(); return

        try:
            element = self.model.elements[self.current_eid]
            element.pid = int(self.pid_combo.currentText())

            # Recreate node list, handling potential grounding (e.g., CBUSH)
            new_nodes = [int(le.text()) for le in self.node_inputs]
            if len(element.nodes) == 1 and len(new_nodes) == 1:
                element.nodes = new_nodes
            elif len(new_nodes) == 2:
                 element.nodes = new_nodes

            if self.orientation_widget:
                new_orient = self.orientation_widget.get_orientation()
                element.g0 = element.x = None # CID is not on CBEAM/CBAR
                if hasattr(element, 'cid'): element.cid = None

                method = new_orient.get('method')
                if method == 'vector': element.x = new_orient['values']
                elif method == 'node': element.g0 = new_orient['values'][0]
                elif method == 'cid' and hasattr(element, 'cid'):
                    element.cid = new_orient['values'][0]

            # Save beam offset data
            if self.offset_fields and element.type == 'CBEAM':
                element.offt = self.offt_combo.currentText()
                element.wa = np.array([
                    float(self.wa_x.text()),
                    float(self.wa_y.text()),
                    float(self.wa_z.text()),
                ])
                element.wb = np.array([
                    float(self.wb_x.text()),
                    float(self.wb_y.text()),
                    float(self.wb_z.text()),
                ])

            super().accept()
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Failed to save element: {e}")
