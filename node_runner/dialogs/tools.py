import os
import math
import random
import numpy as np

from PySide6 import QtCore, QtGui
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGridLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QComboBox, QCheckBox, QRadioButton,
    QDialogButtonBox, QWidget, QListWidget, QListWidgetItem, QMessageBox,
    QColorDialog, QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea,
    QTabWidget, QStackedLayout, QFrame
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPalette, QColor, QPixmap, QDoubleValidator, QAction


class NodeTransformDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Transform Nodes"); self.setMinimumWidth(500); main_layout = QVBoxLayout(self); self.tabs = QTabWidget()
        self.tabs.addTab(self._create_translate_tab(), "Translate (Delta)"); self.tabs.addTab(self._create_move_to_point_tab(), "Move to Point")
        self.tabs.addTab(self._create_scale_tab(), "Scale"); self.tabs.addTab(self._create_rotate_tab(), "Rotate")
        main_layout.addWidget(self.tabs); button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box); button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
    def _create_translate_tab(self):
        widget = QWidget(); layout = QFormLayout(widget); self.trans_dx, self.trans_dy, self.trans_dz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        layout.addRow("ΔX:", self.trans_dx); layout.addRow("ΔY:", self.trans_dy); layout.addRow("ΔZ:", self.trans_dz); return widget
    def _create_move_to_point_tab(self):
        widget = QWidget(); layout = QFormLayout(widget); self.move_x, self.move_y, self.move_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        layout.addRow("Target X:", self.move_x); layout.addRow("Target Y:", self.move_y); layout.addRow("Target Z:", self.move_z); return widget
    def _create_scale_tab(self):
        widget = QWidget(); main_layout = QVBoxLayout(widget); form_layout = QFormLayout()
        self.scale_fx, self.scale_fy, self.scale_fz = QLineEdit("1.0"), QLineEdit("1.0"), QLineEdit("1.0")
        form_layout.addRow("Factor X:", self.scale_fx); form_layout.addRow("Factor Y:", self.scale_fy); form_layout.addRow("Factor Z:", self.scale_fz); main_layout.addLayout(form_layout)
        center_group = QGroupBox("Scale Center"); center_layout = QVBoxLayout(center_group)
        self.scale_origin_rb, self.scale_centroid_rb, self.scale_custom_rb = QRadioButton("Global Origin (0,0,0)", checked=True), QRadioButton("Selection Centroid"), QRadioButton("Custom Point")
        self.scale_cx, self.scale_cy, self.scale_cz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        for w in [self.scale_cx, self.scale_cy, self.scale_cz]: w.setEnabled(False)
        custom_layout = QFormLayout(); custom_layout.addRow("Center X:", self.scale_cx); custom_layout.addRow("Center Y:", self.scale_cy); custom_layout.addRow("Center Z:", self.scale_cz)
        center_layout.addWidget(self.scale_origin_rb); center_layout.addWidget(self.scale_centroid_rb); center_layout.addWidget(self.scale_custom_rb); center_layout.addLayout(custom_layout); main_layout.addWidget(center_group)
        self.scale_custom_rb.toggled.connect(lambda c: [w.setEnabled(c) for w in [self.scale_cx, self.scale_cy, self.scale_cz]]); return widget
    def _create_rotate_tab(self):
        widget = QWidget(); main_layout = QVBoxLayout(widget); form_layout = QFormLayout()
        self.rot_angle = QLineEdit("0.0"); form_layout.addRow("Angle (degrees):", self.rot_angle)
        axis_group = QGroupBox("Axis of Rotation"); axis_layout = QHBoxLayout(axis_group)
        self.rot_x_rb, self.rot_y_rb, self.rot_z_rb = QRadioButton("X-Axis"), QRadioButton("Y-Axis"), QRadioButton("Z-Axis", checked=True)
        axis_layout.addWidget(self.rot_x_rb); axis_layout.addWidget(self.rot_y_rb); axis_layout.addWidget(self.rot_z_rb); main_layout.addLayout(form_layout); main_layout.addWidget(axis_group)
        center_group = QGroupBox("Center of Rotation"); center_layout = QVBoxLayout(center_group)
        self.rot_origin_rb, self.rot_centroid_rb, self.rot_custom_rb = QRadioButton("Global Origin (0,0,0)", checked=True), QRadioButton("Selection Centroid"), QRadioButton("Custom Point")
        self.rot_cx, self.rot_cy, self.rot_cz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        for w in [self.rot_cx, self.rot_cy, self.rot_cz]: w.setEnabled(False)
        custom_layout = QFormLayout(); custom_layout.addRow("Center X:", self.rot_cx); custom_layout.addRow("Center Y:", self.rot_cy); custom_layout.addRow("Center Z:", self.rot_cz)
        center_layout.addWidget(self.rot_origin_rb); center_layout.addWidget(self.rot_centroid_rb); center_layout.addWidget(self.rot_custom_rb); center_layout.addLayout(custom_layout); main_layout.addWidget(center_group)
        self.rot_custom_rb.toggled.connect(lambda c: [w.setEnabled(c) for w in [self.rot_cx, self.rot_cy, self.rot_cz]]); return widget
    def get_transform_parameters(self):
        params = {}
        try:
            if (idx := self.tabs.currentIndex()) == 0: params.update({'type': 'translate', 'delta': [float(self.trans_dx.text()), float(self.trans_dy.text()), float(self.trans_dz.text())]})
            elif idx == 1: params.update({'type': 'move_to_point', 'target': [float(self.move_x.text()), float(self.move_y.text()), float(self.move_z.text())]})
            elif idx == 2:
                params.update({'type': 'scale', 'factors': [float(self.scale_fx.text()), float(self.scale_fy.text()), float(self.scale_fz.text())]})
                if self.scale_centroid_rb.isChecked(): params['center_type'] = 'centroid'
                elif self.scale_custom_rb.isChecked(): params.update({'center_type': 'custom', 'custom_center': [float(self.scale_cx.text()), float(self.scale_cy.text()), float(self.scale_cz.text())]})
                else: params['center_type'] = 'origin'
            elif idx == 3:
                params.update({'type': 'rotate', 'angle': float(self.rot_angle.text())})
                if self.rot_x_rb.isChecked(): params['axis'] = 'x'
                elif self.rot_y_rb.isChecked(): params['axis'] = 'y'
                else: params['axis'] = 'z'
                if self.rot_centroid_rb.isChecked(): params['center_type'] = 'centroid'
                elif self.rot_custom_rb.isChecked(): params.update({'center_type': 'custom', 'custom_center': [float(self.rot_cx.text()), float(self.rot_cy.text()), float(self.rot_cz.text())]})
                else: params['center_type'] = 'origin'
        except ValueError: QMessageBox.warning(self, "Input Error", "All numerical inputs must be valid floats."); return None
        return params


class ColorManagerDialog(QDialog):
    def __init__(self, color_map, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Color Manager"); self.setMinimumWidth(300); self.color_map, self.color_buttons = color_map.copy(), {}
        main_layout = QVBoxLayout(self); scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll_content = QWidget(); self.grid_layout = QGridLayout(scroll_content); scroll.setWidget(scroll_content)
        for row, (name, color) in enumerate(self.color_map.items()):
            label, btn = QLabel(str(name)), QPushButton(); btn.setFixedSize(24, 24); btn.setStyleSheet(f"background-color: {color}; border: 1px solid #45475a;")
            btn.clicked.connect(lambda c=False, n=name, b=btn: self.pick_color(n, b)); self.color_buttons[name] = btn; self.grid_layout.addWidget(label, row, 0); self.grid_layout.addWidget(btn, row, 1)
        main_layout.addWidget(scroll); button_box = QHBoxLayout(); randomize_btn, ok_btn, cancel_btn = QPushButton("Randomize"), QPushButton("OK"), QPushButton("Cancel")
        button_box.addWidget(randomize_btn); button_box.addStretch(); button_box.addWidget(ok_btn); button_box.addWidget(cancel_btn); main_layout.addLayout(button_box)
        randomize_btn.clicked.connect(self.randomize_colors); ok_btn.clicked.connect(self.accept); cancel_btn.clicked.connect(self.reject)
    def pick_color(self, name, btn):
        if (color := QColorDialog.getColor(QtGui.QColor(self.color_map[name]), self)).isValid(): self.color_map[name] = color.name(); btn.setStyleSheet(f"background-color: {self.color_map[name]}; border: 1px solid #45475a;")
    def randomize_colors(self):
        for name, btn in self.color_buttons.items():
            color = QtGui.QColor(random.randint(50, 220), random.randint(50, 220), random.randint(50, 220)); self.color_map[name] = color.name(); btn.setStyleSheet(f"background-color: {self.color_map[name]}; border: 1px solid #45475a;")


class GeneratorDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fuselage Generator")
        self.setMinimumWidth(500)
        self.dim_inputs = {}
        main_layout = QVBoxLayout(self)
        top_row_layout = QHBoxLayout()
        global_group = QGroupBox("Global Geometry")
        form = QFormLayout(global_group)
        self.radius_input = QLineEdit("99.0")
        form.addRow("Fuselage Radius (in):", self.radius_input)
        self.bays_input = QLineEdit("6")
        form.addRow("Number of Bays:", self.bays_input)
        self.spacing_input = QLineEdit("20.0")
        form.addRow("Frame Spacing (in):", self.spacing_input)
        top_row_layout.addWidget(global_group)
        skin_group = QGroupBox("Skin Properties")
        form = QFormLayout(skin_group)
        self.skin_elem_combo = QComboBox()
        self.skin_elem_combo.addItems(["CQUAD4", "CMEMBRAN"])
        form.addRow("Element Type:", self.skin_elem_combo)
        self.skin_thick_input = QLineEdit("0.063")
        form.addRow("Thickness (in):", self.skin_thick_input)
        self.skin_mat_combo = QComboBox()
        form.addRow("Material:", self.skin_mat_combo)
        top_row_layout.addWidget(skin_group)
        main_layout.addLayout(top_row_layout)
        props_group = QGroupBox("Stringer & Frame Properties")
        props_grid_layout = QGridLayout(props_group)
        props_grid_layout.setColumnStretch(1, 0)
        props_grid_layout.setColumnStretch(2, 1)
        self._setup_component_ui(
            prefix="stringer",
            grid=props_grid_layout,
            start_row=0,
            defaults={
                "num": "90", "elem_type": "CBEAM", "section_type": "J",
                "dims": {
                    "L": {"H": "2.5", "W": "1.0", "t": "0.06"},
                    "T": {"W": "1.5", "tf": "0.05", "H": "2.0", "tw": "0.05"},
                    "J": {"H": "2.0", "W1": "1.0", "W2": "0.8", "t": "0.05"}
                }
            }
        )
        separator = QFrame(); separator.setFrameShape(QFrame.HLine); separator.setFrameShadow(QFrame.Sunken)
        props_grid_layout.addWidget(separator, 5, 0, 1, 4)
        self._setup_component_ui(
            prefix="frame",
            grid=props_grid_layout,
            start_row=6,
            defaults={
                "elem_type": "CBEAM", "section_type": "Z",
                "dims": {
                    "L": {"H": "5.0", "W": "1.5", "t": "0.08"},
                    "T": {"W": "2.0", "tf": "0.06", "H": "2.5", "tw": "0.06"},
                    "J": {"H": "2.0", "W1": "1.0", "W2": "0.8", "t": "0.05"}
                }
            }
        )
        main_layout.addWidget(props_group)
        self.floor_check = QCheckBox("Add Floor Beams and Stanchions")
        main_layout.addWidget(self.floor_check)
        self.floor_group = QGroupBox("Floor Structure Details")
        floor_layout = QFormLayout(self.floor_group)
        self.floor_z_input = QLineEdit("0.0")
        self.floor_elem_count_input = QLineEdit("5")
        int_validator = QtGui.QIntValidator(self); int_validator.setBottom(3)
        self.floor_elem_count_input.setValidator(int_validator)
        self.stanchion_y_input = QLineEdit("-40.0, 40.0")
        self.stanchion_elem_count_input = QLineEdit("3")
        floor_beam_widget = QWidget(); floor_beam_h_layout = QHBoxLayout(floor_beam_widget)
        floor_beam_h_layout.setContentsMargins(0, 0, 0, 0); floor_beam_h_layout.addWidget(self.floor_z_input)
        floor_beam_h_layout.addWidget(QLabel("Num. Elements:")); floor_beam_h_layout.addWidget(self.floor_elem_count_input)
        floor_layout.addRow("Floor Waterline (Z-coord):", floor_beam_widget)
        stanchion_widget = QWidget(); stanchion_h_layout = QHBoxLayout(stanchion_widget)
        stanchion_h_layout.setContentsMargins(0, 0, 0, 0); stanchion_h_layout.addWidget(self.stanchion_y_input)
        stanchion_h_layout.addWidget(QLabel("Num. Elements:")); stanchion_h_layout.addWidget(self.stanchion_elem_count_input)
        floor_layout.addRow("Stanchion Butt Lines (Y, comma-sep):", stanchion_widget)
        connection_group = QGroupBox("Connection Method"); connection_layout = QVBoxLayout(connection_group)
        self.snap_method_rb = QRadioButton("Snap to Closest Skin Node (Simple)"); self.modify_method_rb = QRadioButton("Modify Fuselage Mesh (Accurate)")
        self.snap_method_rb.setChecked(True); self.vertical_stanchions_check = QCheckBox("Make Stanchions Vertical")
        snap_layout = QHBoxLayout(); snap_layout.addWidget(self.snap_method_rb); snap_layout.addWidget(self.vertical_stanchions_check)
        connection_layout.addLayout(snap_layout); connection_layout.addWidget(self.modify_method_rb)
        floor_layout.addRow(connection_group)
        self.snap_method_rb.toggled.connect(self.vertical_stanchions_check.setEnabled)
        self.floor_group.setVisible(False)
        main_layout.addWidget(self.floor_group)
        self.floor_check.toggled.connect(self.floor_group.setVisible)
        main_layout.addStretch()
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)
        self._update_dims_ui("frame")
        self._update_dims_ui("stringer")
        self._update_image("frame")
        self._update_image("stringer")

    def _setup_component_ui(self, prefix, grid, start_row, defaults):
        is_stringer = (prefix == "stringer")
        header = QLabel(f"<b>{prefix.capitalize()}s</b>")
        grid.addWidget(header, start_row, 0, 1, 2)
        if is_stringer:
            num_input = QLineEdit(defaults.get("num", "0"))
            grid.addWidget(QLabel("Number of Stringers:"), start_row + 1, 0)
            grid.addWidget(num_input, start_row + 1, 1)
            setattr(self, "num_stringers_input", num_input)
        row_offset = 2 if is_stringer else 1
        elem_combo = QComboBox(); elem_combo.addItems(["CBEAM", "CBAR", "CROD"]); elem_combo.setCurrentText(defaults.get("elem_type"))
        grid.addWidget(QLabel("Element Type:"), start_row + row_offset, 0)
        grid.addWidget(elem_combo, start_row + row_offset, 1)
        type_combo = QComboBox(); type_combo.addItems(["L", "T", "Z", "J"]); type_combo.setCurrentText(defaults.get("section_type"))
        grid.addWidget(QLabel("Section Type:"), start_row + row_offset + 1, 0)
        grid.addWidget(type_combo, start_row + row_offset + 1, 1)
        mat_combo = QComboBox()
        grid.addWidget(QLabel("Material:"), start_row + row_offset + 2, 0)
        grid.addWidget(mat_combo, start_row + row_offset + 2, 1)
        dims_widgets, stacked_layout = self._create_dims_widgets(prefix, defaults.get("dims", {}))
        grid.addWidget(dims_widgets, start_row + 1, 2, 4, 1)
        image_label = QLabel()
        grid.addWidget(image_label, start_row, 3, 5, 1)
        setattr(self, f"{prefix}_elem_combo", elem_combo)
        setattr(self, f"{prefix}_type_combo", type_combo)
        setattr(self, f"{prefix}_mat_combo", mat_combo)
        setattr(self, f"{prefix}_image_label", image_label)
        setattr(self, f"{prefix}_dims_stack", stacked_layout)
        type_combo.currentIndexChanged.connect(lambda: self._update_dims_ui(prefix))
        type_combo.currentIndexChanged.connect(lambda: self._update_image(prefix))

    def _create_dims_widgets(self, prefix, defaults):
        container = QWidget()
        stacked_layout = QStackedLayout(container)
        section_defs = {
            "L": [("Height H (in):", "H"), ("Width W (in):", "W"), ("Thickness t (in):", "t")],
            "T": [("Flange W (in):", "W"), ("Flange t (in):", "tf"), ("Web H (in):", "H"), ("Web t (in):", "tw")],
            "J": [("Height H (in):", "H"), ("Top Flange W1 (in):", "W1"), ("Bot Flange W2 (in):", "W2"), ("Thickness t (in):", "t")]
        }
        for sec_type, fields in section_defs.items():
            form_widget = QWidget()
            layout = QFormLayout(form_widget); layout.setContentsMargins(0, 0, 0, 0)
            for label, key_suffix in fields:
                widget_key = f"{prefix}_{key_suffix}"
                default_val = defaults.get(sec_type, {}).get(key_suffix, "0.0")
                line_edit = QLineEdit(default_val)
                self.dim_inputs[widget_key] = line_edit
                layout.addRow(label, line_edit)
            stacked_layout.addWidget(form_widget)
        return container, stacked_layout

    def _update_dims_ui(self, prefix):
        type_combo = getattr(self, f"{prefix}_type_combo")
        stack = getattr(self, f"{prefix}_dims_stack")
        sec_type = type_combo.currentText()
        if sec_type in ["L", "Z"]: stack.setCurrentIndex(0)
        elif sec_type == "T": stack.setCurrentIndex(1)
        elif sec_type == "J": stack.setCurrentIndex(2)
        else: stack.setCurrentWidget(QWidget())

    def _update_image(self, prefix):
        combo_box = getattr(self, f"{prefix}_type_combo")
        image_label = getattr(self, f"{prefix}_image_label")
        img_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "images", f"{combo_box.currentText()}_section.png")
        if os.path.exists(img_path):
            pixmap = QtGui.QPixmap(img_path)
            scaled_pixmap = pixmap.scaled(150, 150, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
            image_label.setPixmap(scaled_pixmap)
        else:
            image_label.clear()

    def set_materials(self, model):
        material_ids = [str(mid) for mid in model.materials.keys()]
        for combo in [self.skin_mat_combo, self.frame_mat_combo, self.stringer_mat_combo]:
            combo.clear()
            combo.addItems(material_ids)

    def _get_dims_from_ui(self, prefix, section_type):
        dims = {}
        section_keys = {
            "L": ["H", "W", "t"], "Z": ["H", "W", "t"],
            "T": ["W", "tf", "H", "tw"],
            "J": ["H", "W1", "W2", "t"]
        }
        if section_type in section_keys:
            for key_suffix in section_keys[section_type]:
                widget_key = f"{prefix}_{key_suffix}"
                if widget_key in self.dim_inputs:
                    dims[key_suffix] = float(self.dim_inputs[widget_key].text())
        return dims

    def get_parameters(self):
        params = {
            "fuselage_radius": float(self.radius_input.text()), "num_bays": int(self.bays_input.text()), "frame_spacing": float(self.spacing_input.text()),
            "skin_element_type": self.skin_elem_combo.currentText(), "skin_thickness": float(self.skin_thick_input.text()), "skin_material_id": int(self.skin_mat_combo.currentText()),
            "frame_element_type": self.frame_elem_combo.currentText(), "frame_section_type": self.frame_type_combo.currentText(), "frame_material_id": int(self.frame_mat_combo.currentText()),
            "frame_dims": self._get_dims_from_ui("frame", self.frame_type_combo.currentText()),
            "num_stringers": int(self.num_stringers_input.text()), "stringer_element_type": self.stringer_elem_combo.currentText(), "stringer_section_type": self.stringer_type_combo.currentText(), "stringer_material_id": int(self.stringer_mat_combo.currentText()),
            "stringer_dims": self._get_dims_from_ui("stringer", self.stringer_type_combo.currentText()),
            "add_floor": self.floor_check.isChecked()
        }
        if params["add_floor"]:
            stanchion_y_coords = [float(val.strip()) for val in self.stanchion_y_input.text().split(',') if val.strip()]
            params.update({
                "floor_z": float(self.floor_z_input.text()), "floor_elem_count": int(self.floor_elem_count_input.text()),
                "stanchion_y_coords": stanchion_y_coords, "stanchion_elem_count": int(self.stanchion_elem_count_input.text()),
                "connection_method": 'modify' if self.modify_method_rb.isChecked() else 'snap',
                "force_vertical": self.vertical_stanchions_check.isChecked() and self.snap_method_rb.isChecked()
            })
        return params


class CoincidentNodeDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Coincident Node Checker")
        self.setMinimumSize(450, 400) # Increased size slightly for new buttons

        main_layout = QVBoxLayout(self)

        # --- Top panel for user input ---
        input_group = QGroupBox("Input Parameters")
        form_layout = QFormLayout(input_group)

        self.tolerance_input = QLineEdit("0.001")
        validator = QDoubleValidator(0.0, 1.0, 8, self)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.tolerance_input.setValidator(validator)
        form_layout.addRow("Merge Tolerance:", self.tolerance_input)

        self.find_button = QPushButton("Find Coincident Nodes")
        form_layout.addWidget(self.find_button)
        main_layout.addWidget(input_group)

        # --- Results panel for displaying and editing found sets ---
        results_group = QGroupBox("Found Sets")
        results_layout = QVBoxLayout(results_group)

        self.results_list = QListWidget()
        # --- MODIFIED: Enable selection of multiple items ---
        self.results_list.setSelectionMode(QListWidget.ExtendedSelection)
        results_layout.addWidget(self.results_list)

        # --- NEW: Add buttons for removing items from the list ---
        edit_buttons_layout = QHBoxLayout()
        remove_button = QPushButton("Remove Selected")
        clear_button = QPushButton("Clear All")
        edit_buttons_layout.addStretch() # Pushes buttons to the right
        edit_buttons_layout.addWidget(remove_button)
        edit_buttons_layout.addWidget(clear_button)
        results_layout.addLayout(edit_buttons_layout)

        main_layout.addWidget(results_group)

        # --- Bottom action buttons ---
        button_box = QDialogButtonBox()
        self.merge_button = button_box.addButton("Merge Selected", QDialogButtonBox.AcceptRole)
        self.close_button = button_box.addButton(QDialogButtonBox.Close)
        self.merge_button.setEnabled(False)
        main_layout.addWidget(button_box)

        # --- Connect signals to slots ---
        self.find_button.clicked.connect(self._find_nodes)
        remove_button.clicked.connect(self._remove_selected_items)
        clear_button.clicked.connect(self._clear_all_items)
        self.merge_button.clicked.connect(self.accept)
        self.close_button.clicked.connect(self.reject)

    def _find_nodes(self):
        """Calls the backend find method and populates the results list."""
        self._clear_all_items()
        try:
            tolerance = float(self.tolerance_input.text())
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Please enter a valid tolerance value.")
            return

        if not self.main_window.current_generator:
            QMessageBox.warning(self, "No Model", "No model is loaded.")
            return

        duplicate_groups = self.main_window.current_generator.find_duplicate_nodes(tolerance)

        if not duplicate_groups:
            self.results_list.addItem("No coincident nodes found with this tolerance.")
            return

        for i, group in enumerate(duplicate_groups):
            master_nid = min(group)
            others = [str(nid) for nid in group if nid != master_nid]
            item_text = f"Group {i+1}: Merge {', '.join(others)} -> into {master_nid}"

            # --- MODIFIED: Store raw data with the list item ---
            item = QListWidgetItem(item_text)
            item.setData(QtCore.Qt.UserRole, group) # Store the raw list e.g., [11, 12, 10]
            self.results_list.addItem(item)

        self._update_merge_button_state()

    # --- NEW: Method to remove only the selected items ---
    def _remove_selected_items(self):
        selected_items = self.results_list.selectedItems()
        if not selected_items:
            return
        for item in selected_items:
            # Taking the item from the list widget also deletes it
            self.results_list.takeItem(self.results_list.row(item))
        self._update_merge_button_state()

    # --- NEW: Method to clear all items from the list ---
    def _clear_all_items(self):
        self.results_list.clear()
        self._update_merge_button_state()

    # --- NEW: Helper to dynamically set the merge button state ---
    def _update_merge_button_state(self):
        # Enable the merge button only if there are items in the list.
        is_enabled = self.results_list.count() > 0 and \
                     "No coincident nodes found" not in self.results_list.item(0).text()
        self.merge_button.setEnabled(is_enabled)
        # Also update the button text to be more accurate
        self.merge_button.setText(f"Merge {self.results_list.count()} Groups") if is_enabled else "Merge Selected"

    # --- NEW: Handle the delete key press for convenience ---
    def keyPressEvent(self, event):
        """Handles the Delete key to remove selected items from the list."""
        if event.key() == QtCore.Qt.Key_Delete and self.results_list.hasFocus():
            self._remove_selected_items()
        else:
            super().keyPressEvent(event) # Pass on other key events

    # --- NEW: Getter method to retrieve results (Revised Phase 4) ---
    def get_groups_to_merge(self):
        """Retrieves the raw node group data from the items remaining in the list."""
        groups = []
        for i in range(self.results_list.count()):
            item = self.results_list.item(i)
            # Ensure we don't process the "No nodes found" informational message
            if item.data(QtCore.Qt.UserRole):
                groups.append(item.data(QtCore.Qt.UserRole))
        return groups


class DuplicateElementDialog(QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.main_window = main_window
        self.setWindowTitle("Duplicate Element Checker")
        self.setMinimumSize(450, 400)

        main_layout = QVBoxLayout(self)

        # --- Top panel ---
        input_group = QGroupBox("Find Duplicate Elements")
        input_layout = QVBoxLayout(input_group)
        info_label = QLabel("Finds elements that share the exact same set of nodes.\n"
                            "The lowest EID in each group is kept; others are marked for deletion.")
        info_label.setWordWrap(True)
        input_layout.addWidget(info_label)
        self.find_button = QPushButton("Find Duplicates")
        input_layout.addWidget(self.find_button)
        main_layout.addWidget(input_group)

        # --- Results panel ---
        results_group = QGroupBox("Found Groups")
        results_layout = QVBoxLayout(results_group)
        self.results_list = QListWidget()
        self.results_list.setSelectionMode(QListWidget.ExtendedSelection)
        results_layout.addWidget(self.results_list)

        edit_buttons_layout = QHBoxLayout()
        remove_button = QPushButton("Remove Selected")
        clear_button = QPushButton("Clear All")
        edit_buttons_layout.addStretch()
        edit_buttons_layout.addWidget(remove_button)
        edit_buttons_layout.addWidget(clear_button)
        results_layout.addLayout(edit_buttons_layout)
        main_layout.addWidget(results_group)

        # --- Bottom action buttons ---
        button_box = QDialogButtonBox()
        self.delete_button = button_box.addButton("Delete Duplicates", QDialogButtonBox.AcceptRole)
        self.close_button = button_box.addButton(QDialogButtonBox.Close)
        self.delete_button.setEnabled(False)
        main_layout.addWidget(button_box)

        # --- Signals ---
        self.find_button.clicked.connect(self._find_duplicates)
        remove_button.clicked.connect(self._remove_selected_items)
        clear_button.clicked.connect(self._clear_all_items)
        self.delete_button.clicked.connect(self.accept)
        self.close_button.clicked.connect(self.reject)

    def _find_duplicates(self):
        self._clear_all_items()
        if not self.main_window.current_generator:
            QMessageBox.warning(self, "No Model", "No model is loaded.")
            return

        model = self.main_window.current_generator.model
        from collections import defaultdict
        sig_map = defaultdict(list)

        for eid, elem in model.elements.items():
            sig = tuple(sorted(elem.node_ids))
            sig_map[sig].append(eid)
        for eid, elem in model.rigid_elements.items():
            sig = tuple(sorted(elem.node_ids))
            sig_map[sig].append(eid)

        dup_groups = {sig: sorted(eids) for sig, eids in sig_map.items() if len(eids) > 1}

        if not dup_groups:
            self.results_list.addItem("No duplicate elements found.")
            return

        for i, (sig, eids) in enumerate(sorted(dup_groups.items(), key=lambda x: x[1][0])):
            keep_eid = eids[0]
            delete_eids = eids[1:]
            item_text = f"Group {i+1}: EIDs {eids} — keep {keep_eid}, delete {', '.join(str(e) for e in delete_eids)}"
            item = QListWidgetItem(item_text)
            item.setData(Qt.UserRole, {'keep': keep_eid, 'delete': delete_eids})
            self.results_list.addItem(item)

        self._update_delete_button_state()

    def _remove_selected_items(self):
        for item in self.results_list.selectedItems():
            self.results_list.takeItem(self.results_list.row(item))
        self._update_delete_button_state()

    def _clear_all_items(self):
        self.results_list.clear()
        self._update_delete_button_state()

    def _update_delete_button_state(self):
        has_items = self.results_list.count() > 0 and \
                    self.results_list.item(0).data(Qt.UserRole) is not None
        self.delete_button.setEnabled(has_items)
        if has_items:
            self.delete_button.setText(f"Delete Duplicates ({self.results_list.count()} groups)")
        else:
            self.delete_button.setText("Delete Duplicates")

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Delete and self.results_list.hasFocus():
            self._remove_selected_items()
        else:
            super().keyPressEvent(event)

    def get_eids_to_delete(self):
        eids = []
        for i in range(self.results_list.count()):
            item = self.results_list.item(i)
            data = item.data(Qt.UserRole)
            if data:
                eids.extend(data['delete'])
        return eids


class FindReplaceDialog(QDialog):
    """Dialog to reassign property or material across elements."""
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        self.setWindowTitle("Find/Replace Property or Material")
        self.setMinimumWidth(400)

        main_layout = QVBoxLayout(self)

        # Mode selection
        mode_group = QGroupBox("Mode")
        mode_layout = QHBoxLayout(mode_group)
        self.prop_rb = QRadioButton("Replace Property", checked=True)
        self.mat_rb = QRadioButton("Replace Material")
        mode_layout.addWidget(self.prop_rb)
        mode_layout.addWidget(self.mat_rb)
        main_layout.addWidget(mode_group)

        # From / To combos
        form = QFormLayout()
        self.from_combo = QComboBox()
        self.to_combo = QComboBox()
        form.addRow("From:", self.from_combo)
        form.addRow("To:", self.to_combo)
        main_layout.addLayout(form)

        # Preview
        self.preview_label = QLabel("Select source and target to see affected count.")
        self.preview_label.setWordWrap(True)
        self.preview_label.setStyleSheet("font-size: 11px; color: #a6adc8; padding: 5px;")
        main_layout.addWidget(self.preview_label)

        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

        # Signals
        self.prop_rb.toggled.connect(self._refresh_combos)
        self.from_combo.currentIndexChanged.connect(self._update_preview)
        self.to_combo.currentIndexChanged.connect(self._update_preview)

        self._refresh_combos()

    def _refresh_combos(self):
        self.from_combo.clear()
        self.to_combo.clear()
        if self.prop_rb.isChecked():
            for pid in sorted(self.model.properties.keys()):
                prop = self.model.properties[pid]
                label = f"PID {pid} ({prop.type})"
                self.from_combo.addItem(label, pid)
                self.to_combo.addItem(label, pid)
        else:
            for mid in sorted(self.model.materials.keys()):
                mat = self.model.materials[mid]
                label = f"MID {mid} ({mat.type})"
                self.from_combo.addItem(label, mid)
                self.to_combo.addItem(label, mid)
        self._update_preview()

    def _update_preview(self):
        from_id = self.from_combo.currentData()
        to_id = self.to_combo.currentData()
        if from_id is None or to_id is None:
            self.preview_label.setText("No items to replace.")
            return
        if from_id == to_id:
            self.preview_label.setText("Source and target are the same.")
            return

        if self.prop_rb.isChecked():
            count = sum(1 for e in self.model.elements.values()
                        if hasattr(e, 'pid') and e.pid == from_id)
            self.preview_label.setText(f"{count} element(s) will change from PID {from_id} to PID {to_id}.")
        else:
            count = sum(1 for p in self.model.properties.values()
                        if hasattr(p, 'mid') and p.mid == from_id)
            self.preview_label.setText(f"{count} property(ies) will change from MID {from_id} to MID {to_id}.")

    def get_parameters(self):
        from_id = self.from_combo.currentData()
        to_id = self.to_combo.currentData()
        if from_id is None or to_id is None or from_id == to_id:
            QMessageBox.warning(self, "Input Error", "Please select different source and target.")
            return None
        mode = 'property' if self.prop_rb.isChecked() else 'material'
        return {'mode': mode, 'from_id': from_id, 'to_id': to_id}


class RenumberDialog(QDialog):
    """Dialog for renumbering node and element IDs."""
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model
        self.setWindowTitle("Renumber Nodes/Elements")
        self.setMinimumWidth(380)

        main_layout = QVBoxLayout(self)
        form = QFormLayout()

        self.renumber_nodes_check = QCheckBox("Renumber Nodes")
        self.renumber_nodes_check.setChecked(True)
        self.renumber_elems_check = QCheckBox("Renumber Elements")
        self.renumber_elems_check.setChecked(True)
        main_layout.addWidget(self.renumber_nodes_check)
        main_layout.addWidget(self.renumber_elems_check)

        self.start_nid = QLineEdit("1")
        self.start_eid = QLineEdit("1")
        form.addRow("Start Node ID:", self.start_nid)
        form.addRow("Start Element ID:", self.start_eid)
        main_layout.addLayout(form)

        # Info labels
        nids = sorted(model.nodes.keys()) if model.nodes else []
        eids = sorted(list(model.elements.keys()) + list(model.rigid_elements.keys())) if model.elements else []
        nid_range = f"{nids[0]}-{nids[-1]}" if nids else "none"
        eid_range = f"{eids[0]}-{eids[-1]}" if eids else "none"
        info = QLabel(f"Current nodes: {nid_range} ({len(nids)} total)\n"
                      f"Current elements: {eid_range} ({len(eids)} total)")
        info.setStyleSheet("font-size: 11px; color: #a6adc8; padding: 5px;")
        main_layout.addWidget(info)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

    def get_parameters(self):
        try:
            return {
                'start_nid': int(self.start_nid.text()),
                'start_eid': int(self.start_eid.text()),
                'do_nodes': self.renumber_nodes_check.isChecked(),
                'do_elements': self.renumber_elems_check.isChecked(),
            }
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Start IDs must be valid integers.")
            return None


class ImportOptionsDialog(QDialog):
    """A simple dialog to ask the user how to import a BDF file."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import Options")

        main_layout = QVBoxLayout(self)
        main_layout.addWidget(QLabel("How would you like to import this file?"))

        self.new_model_rb = QRadioButton("Start a new model from this file")
        self.append_model_rb = QRadioButton("Append cards to the current model")
        self.append_model_rb.setChecked(True)

        main_layout.addWidget(self.new_model_rb)
        main_layout.addWidget(self.append_model_rb)

        note_label = QLabel(
            "<b>Note:</b> Appending will overwrite any existing entities (nodes, elements, etc.) "
            "that have the same ID as an entity in the imported file."
        )
        note_label.setWordWrap(True)
        note_label.setStyleSheet("font-size: 11px; color: #a6adc8; padding-top: 5px;")
        main_layout.addWidget(note_label)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

    def get_selected_option(self):
        if self.new_model_rb.isChecked():
            return 'new'
        return 'append'


class ImportCADDialog(QDialog):
    """Dialog for CAD file import settings (STEP/IGES/STL)."""

    def __init__(self, filepath, available_pids=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import CAD File")
        self.setMinimumWidth(400)
        layout = QVBoxLayout(self)

        # File info
        file_group = QGroupBox("File")
        file_layout = QFormLayout(file_group)
        self.filepath = filepath
        file_layout.addRow("File:", QLabel(os.path.basename(filepath)))
        ext = os.path.splitext(filepath)[1].lower()
        file_layout.addRow("Type:", QLabel(ext.upper().lstrip('.')))
        layout.addWidget(file_group)

        # Mesh settings
        mesh_group = QGroupBox("Mesh Settings")
        mesh_layout = QFormLayout(mesh_group)
        self.mesh_size_input = QLineEdit("10.0")
        self.mesh_size_input.setValidator(QDoubleValidator(0.001, 1e6, 4, self))
        mesh_layout.addRow("Mesh Size:", self.mesh_size_input)

        self.elem_pref_combo = QComboBox()
        self.elem_pref_combo.addItems(["Quad (CQUAD4)", "Tri (CTRIA3)"])
        mesh_layout.addRow("Element Type:", self.elem_pref_combo)

        # Disable mesh settings for STL (already meshed)
        if ext == '.stl':
            self.mesh_size_input.setEnabled(False)
            self.elem_pref_combo.setEnabled(False)

        layout.addWidget(mesh_group)

        # Property
        prop_group = QGroupBox("Property")
        prop_layout = QFormLayout(prop_group)
        self.pid_input = QLineEdit("1")
        if available_pids:
            self.pid_input.setText(str(min(available_pids)))
        prop_layout.addRow("Property ID (PID):", self.pid_input)
        layout.addWidget(prop_group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_settings(self):
        elem_pref = 'quad' if self.elem_pref_combo.currentIndex() == 0 else 'tri'
        try:
            mesh_size = float(self.mesh_size_input.text())
        except ValueError:
            mesh_size = 10.0
        try:
            pid = int(self.pid_input.text())
        except ValueError:
            pid = 1
        return {
            'filepath': self.filepath,
            'mesh_size': mesh_size,
            'elem_preference': elem_pref,
            'pid': pid,
        }


class DisplaySettingsDialog(QDialog):
    """Dialog for configuring display settings: background, display options, and screenshot."""

    def __init__(self, settings=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Display Settings")
        self.setMinimumWidth(520)
        self.setMinimumHeight(450)

        self.settings = settings or {}

        main_layout = QVBoxLayout(self)
        tabs = QTabWidget()

        tabs.addTab(self._create_background_tab(), "Background")
        tabs.addTab(self._create_display_tab(), "Display Options")
        tabs.addTab(self._create_screenshot_tab(), "Screenshot")

        main_layout.addWidget(tabs)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

    def _create_color_button(self, initial_color="#1e1e2e"):
        """Create a clickable color button."""
        btn = QPushButton()
        btn.setFixedSize(60, 28)
        btn.setProperty("hex_color", initial_color)
        btn.setStyleSheet(
            f"background-color: {initial_color}; border: 1px solid #555; border-radius: 3px;")
        btn.clicked.connect(lambda: self._pick_color(btn))
        return btn

    def _pick_color(self, btn):
        """Open color picker for a color button."""
        current = QColor(btn.property("hex_color"))
        color = QColorDialog.getColor(current, self, "Pick Color")
        if color.isValid():
            hex_color = color.name()
            btn.setProperty("hex_color", hex_color)
            btn.setStyleSheet(
                f"background-color: {hex_color}; border: 1px solid #555; border-radius: 3px;")

    def _create_background_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Mode selector
        mode_group = QGroupBox("Background Mode")
        mode_layout = QVBoxLayout(mode_group)

        self.bg_solid_radio = QRadioButton("Solid Color")
        self.bg_gradient_radio = QRadioButton("Gradient (Top / Bottom)")
        self.bg_image_radio = QRadioButton("Image")

        mode_layout.addWidget(self.bg_solid_radio)
        mode_layout.addWidget(self.bg_gradient_radio)
        mode_layout.addWidget(self.bg_image_radio)

        # Solid color options
        self.solid_color_widget = QWidget()
        solid_layout = QHBoxLayout(self.solid_color_widget)
        solid_layout.setContentsMargins(20, 0, 0, 0)
        solid_layout.addWidget(QLabel("Color:"))
        self.solid_color_btn = self._create_color_button(
            self.settings.get('bg_color', '#1e1e2e'))
        solid_layout.addWidget(self.solid_color_btn)
        solid_layout.addStretch()
        mode_layout.addWidget(self.solid_color_widget)

        # Gradient options
        self.gradient_widget = QWidget()
        grad_layout = QHBoxLayout(self.gradient_widget)
        grad_layout.setContentsMargins(20, 0, 0, 0)
        grad_layout.addWidget(QLabel("Top:"))
        self.grad_top_btn = self._create_color_button(
            self.settings.get('bg_top_color', '#1a2a4a'))
        grad_layout.addWidget(self.grad_top_btn)
        grad_layout.addWidget(QLabel("Bottom:"))
        self.grad_bottom_btn = self._create_color_button(
            self.settings.get('bg_bottom_color', '#0a0a1e'))
        grad_layout.addWidget(self.grad_bottom_btn)
        grad_layout.addStretch()
        mode_layout.addWidget(self.gradient_widget)

        # Image options
        self.image_widget = QWidget()
        img_layout = QFormLayout(self.image_widget)
        img_layout.setContentsMargins(20, 0, 0, 0)
        self.bg_image_path = QLineEdit(self.settings.get('bg_image_path', ''))
        self.bg_image_path.setPlaceholderText("Path to image (JPG/PNG)")
        browse_layout = QHBoxLayout()
        browse_layout.addWidget(self.bg_image_path)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_image)
        browse_layout.addWidget(browse_btn)
        img_layout.addRow("Image:", browse_layout)
        mode_layout.addWidget(self.image_widget)

        layout.addWidget(mode_group)

        # Presets
        presets_group = QGroupBox("Presets")
        presets_layout = QHBoxLayout(presets_group)
        self.preset_combo = QComboBox()
        from node_runner.theme import BACKGROUND_PRESETS
        self.preset_combo.addItems(list(BACKGROUND_PRESETS.keys()))
        self.preset_combo.currentTextChanged.connect(self._apply_preset)
        apply_preset_btn = QPushButton("Apply Preset")
        apply_preset_btn.clicked.connect(
            lambda: self._apply_preset(self.preset_combo.currentText()))
        presets_layout.addWidget(self.preset_combo, 1)
        presets_layout.addWidget(apply_preset_btn)
        layout.addWidget(presets_group)

        layout.addStretch()

        # Set initial mode
        bg_mode = self.settings.get('bg_mode', 'solid')
        if bg_mode == 'gradient':
            self.bg_gradient_radio.setChecked(True)
        elif bg_mode == 'image':
            self.bg_image_radio.setChecked(True)
        else:
            self.bg_solid_radio.setChecked(True)

        # Connect mode changes to update visibility
        self.bg_solid_radio.toggled.connect(self._update_bg_visibility)
        self.bg_gradient_radio.toggled.connect(self._update_bg_visibility)
        self.bg_image_radio.toggled.connect(self._update_bg_visibility)
        self._update_bg_visibility()

        return widget

    def _update_bg_visibility(self):
        """Show/hide background option widgets based on selected mode."""
        self.solid_color_widget.setVisible(self.bg_solid_radio.isChecked())
        self.gradient_widget.setVisible(self.bg_gradient_radio.isChecked())
        self.image_widget.setVisible(self.bg_image_radio.isChecked())

    def _apply_preset(self, preset_name):
        """Apply a background preset."""
        from node_runner.theme import BACKGROUND_PRESETS
        preset = BACKGROUND_PRESETS.get(preset_name)
        if not preset:
            return
        mode = preset['mode']
        if mode == 'solid':
            self.bg_solid_radio.setChecked(True)
            color = preset['color']
            self.solid_color_btn.setProperty("hex_color", color)
            self.solid_color_btn.setStyleSheet(
                f"background-color: {color}; border: 1px solid #555; border-radius: 3px;")
        elif mode == 'gradient':
            self.bg_gradient_radio.setChecked(True)
            top = preset['top']
            bottom = preset['bottom']
            self.grad_top_btn.setProperty("hex_color", top)
            self.grad_top_btn.setStyleSheet(
                f"background-color: {top}; border: 1px solid #555; border-radius: 3px;")
            self.grad_bottom_btn.setProperty("hex_color", bottom)
            self.grad_bottom_btn.setStyleSheet(
                f"background-color: {bottom}; border: 1px solid #555; border-radius: 3px;")

    def _browse_image(self):
        """Open file dialog to select background image."""
        from PySide6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Background Image", "",
            "Images (*.png *.jpg *.jpeg *.bmp *.tiff)")
        if path:
            self.bg_image_path.setText(path)

    def _create_display_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Axes & general
        axes_group = QGroupBox("Axes && General")
        axes_layout = QFormLayout(axes_group)
        self.show_axes_check = QCheckBox()
        self.show_axes_check.setChecked(self.settings.get('show_axes', True))
        axes_layout.addRow("Show Axes:", self.show_axes_check)
        self.axes_color_btn = self._create_color_button(
            self.settings.get('axes_color', '#cdd6f4'))
        axes_layout.addRow("Axis Label Color:", self.axes_color_btn)
        layout.addWidget(axes_group)

        # Element rendering
        render_group = QGroupBox("Rendering")
        render_layout = QFormLayout(render_group)
        self.edge_color_btn = self._create_color_button(
            self.settings.get('edge_color', '#000000'))
        render_layout.addRow("Edge Color:", self.edge_color_btn)

        from PySide6.QtWidgets import QSpinBox
        self.node_size_spin = QSpinBox()
        self.node_size_spin.setRange(1, 20)
        self.node_size_spin.setValue(self.settings.get('node_size', 8))
        render_layout.addRow("Node Render Size:", self.node_size_spin)

        self.line_width_spin = QSpinBox()
        self.line_width_spin.setRange(1, 10)
        self.line_width_spin.setValue(self.settings.get('line_width', 2))
        render_layout.addRow("Line Element Width:", self.line_width_spin)

        self.selection_color_btn = self._create_color_button(
            self.settings.get('selection_color', '#ffff00'))
        render_layout.addRow("Selection Highlight:", self.selection_color_btn)
        layout.addWidget(render_group)

        # Watermark
        watermark_group = QGroupBox("Watermark / Overlay Text")
        wm_layout = QFormLayout(watermark_group)
        self.watermark_text = QLineEdit(self.settings.get('watermark_text', ''))
        self.watermark_text.setPlaceholderText("Optional text overlay (e.g. company name)")
        wm_layout.addRow("Text:", self.watermark_text)
        self.watermark_color_btn = self._create_color_button(
            self.settings.get('watermark_color', '#ffffff'))
        wm_layout.addRow("Color:", self.watermark_color_btn)
        self.watermark_position = QComboBox()
        self.watermark_position.addItems(
            ["Top-Left", "Top-Right", "Bottom-Left", "Bottom-Right"])
        self.watermark_position.setCurrentText(
            self.settings.get('watermark_position', 'Bottom-Right'))
        wm_layout.addRow("Position:", self.watermark_position)
        layout.addWidget(watermark_group)

        layout.addStretch()
        return widget

    def _create_screenshot_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        res_group = QGroupBox("Resolution")
        res_layout = QFormLayout(res_group)
        self.screenshot_res_combo = QComboBox()
        self.screenshot_res_combo.addItems([
            "Current Viewport",
            "1920 x 1080 (Full HD)",
            "2560 x 1440 (2K)",
            "3840 x 2160 (4K)",
            "Custom",
        ])
        self.screenshot_res_combo.setCurrentText(
            self.settings.get('screenshot_resolution', 'Current Viewport'))
        res_layout.addRow("Resolution:", self.screenshot_res_combo)

        from PySide6.QtWidgets import QSpinBox
        custom_widget = QWidget()
        custom_layout = QHBoxLayout(custom_widget)
        custom_layout.setContentsMargins(0, 0, 0, 0)
        self.custom_width = QSpinBox()
        self.custom_width.setRange(100, 7680)
        self.custom_width.setValue(self.settings.get('custom_width', 1920))
        self.custom_height = QSpinBox()
        self.custom_height.setRange(100, 4320)
        self.custom_height.setValue(self.settings.get('custom_height', 1080))
        custom_layout.addWidget(QLabel("W:"))
        custom_layout.addWidget(self.custom_width)
        custom_layout.addWidget(QLabel("H:"))
        custom_layout.addWidget(self.custom_height)
        res_layout.addRow("Custom Size:", custom_widget)
        layout.addWidget(res_group)

        opts_group = QGroupBox("Options")
        opts_layout = QFormLayout(opts_group)
        self.transparent_bg_check = QCheckBox()
        self.transparent_bg_check.setChecked(
            self.settings.get('transparent_bg', False))
        opts_layout.addRow("Transparent Background:", self.transparent_bg_check)

        self.aa_check = QCheckBox()
        self.aa_check.setChecked(self.settings.get('antialiasing', True))
        opts_layout.addRow("Anti-aliasing:", self.aa_check)

        self.file_format_combo = QComboBox()
        self.file_format_combo.addItems(["PNG", "JPG", "BMP", "TIFF"])
        self.file_format_combo.setCurrentText(
            self.settings.get('file_format', 'PNG'))
        opts_layout.addRow("File Format:", self.file_format_combo)
        layout.addWidget(opts_group)

        # Capture button
        capture_btn = QPushButton("Capture Screenshot Now...")
        capture_btn.setStyleSheet(
            "padding: 8px; font-weight: bold; background-color: #89b4fa; "
            "color: #1e1e2e; border-radius: 4px;")
        capture_btn.clicked.connect(self._capture_now)
        layout.addWidget(capture_btn)

        layout.addStretch()
        return widget

    def _capture_now(self):
        """Signal parent to take a screenshot with current settings."""
        self._capture_requested = True
        self.accept()

    def get_settings(self):
        """Return all display settings as a dict."""
        # Background mode
        if self.bg_gradient_radio.isChecked():
            bg_mode = 'gradient'
        elif self.bg_image_radio.isChecked():
            bg_mode = 'image'
        else:
            bg_mode = 'solid'

        return {
            'bg_mode': bg_mode,
            'bg_color': self.solid_color_btn.property("hex_color"),
            'bg_top_color': self.grad_top_btn.property("hex_color"),
            'bg_bottom_color': self.grad_bottom_btn.property("hex_color"),
            'bg_image_path': self.bg_image_path.text(),
            'show_axes': self.show_axes_check.isChecked(),
            'axes_color': self.axes_color_btn.property("hex_color"),
            'edge_color': self.edge_color_btn.property("hex_color"),
            'node_size': self.node_size_spin.value(),
            'line_width': self.line_width_spin.value(),
            'selection_color': self.selection_color_btn.property("hex_color"),
            'watermark_text': self.watermark_text.text(),
            'watermark_color': self.watermark_color_btn.property("hex_color"),
            'watermark_position': self.watermark_position.currentText(),
            'screenshot_resolution': self.screenshot_res_combo.currentText(),
            'custom_width': self.custom_width.value(),
            'custom_height': self.custom_height.value(),
            'transparent_bg': self.transparent_bg_check.isChecked(),
            'antialiasing': self.aa_check.isChecked(),
            'file_format': self.file_format_combo.currentText(),
            'capture_requested': getattr(self, '_capture_requested', False),
        }
