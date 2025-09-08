# main.py

import sys
import os
import json
import random
import re
import numpy as np
import pyvista as pv
import vtk
import copy
import time
from pyvistaqt import QtInteractor

from PySide6 import QtCore, QtGui
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QFileDialog, QMessageBox, QFormLayout, QGroupBox, QFrame,
    QLineEdit, QComboBox, QLabel, QCheckBox, QColorDialog, QDialog,
    QScrollArea, QGridLayout, QTreeWidget, QTreeWidgetItem, QSplitter,
    QDialogButtonBox, QTreeWidgetItemIterator, QTableWidget, QTableWidgetItem,
    QHeaderView, QListWidget, QListWidgetItem, QTextEdit, QTabWidget, QRadioButton,
    QStackedLayout, QInputDialog, QMenu
)
from PySide6.QtGui import QPalette, QColor, QAction, QActionGroup, QDoubleValidator
from collections import Counter

from pyNastran.bdf.bdf import BDF
from nas import NastranModelGenerator

# --- Define Color Palettes for Theming ---
dark_palette = QPalette()
dark_palette.setColor(QPalette.Window, QColor(53, 53, 53)); dark_palette.setColor(QPalette.WindowText, QtCore.Qt.white); dark_palette.setColor(QPalette.Base, QColor(35, 35, 35)); dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53)); dark_palette.setColor(QPalette.Text, QtCore.Qt.white); dark_palette.setColor(QPalette.Button, QColor(53, 53, 53)); dark_palette.setColor(QPalette.ButtonText, QtCore.Qt.white); dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218)); dark_palette.setColor(QPalette.HighlightedText, QtCore.Qt.black)
light_palette = QPalette()

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

        button_box = QDialogButtonBox(QDialogButtonBox.Close)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.table.itemDoubleClicked.connect(self.edit_material)

    def populate_table(self):
        self.table.setRowCount(len(self.model.materials))
        for row, (mid, mat) in enumerate(sorted(self.model.materials.items())):
            title = mat.comment.strip().lstrip('$').strip() or f'{mat.type} {mid}'
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
                del self.model.materials[mid] 
                self.parent_window._add_material_card_from_params(params)
                self.populate_table()
                self.parent_window._populate_tree()


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

        button_box = QDialogButtonBox(QDialogButtonBox.Close)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        
        self.table.itemDoubleClicked.connect(self.edit_property)

    def populate_table(self):
        self.table.setRowCount(len(self.model.properties))
        for row, (pid, prop) in enumerate(sorted(self.model.properties.items())):
            title = prop.comment.strip().lstrip('$').strip() or f'{prop.type} {pid}'
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
                del self.model.properties[pid]
                self.parent_window._add_property_card_from_params(params)
                self.populate_table()
                self.parent_window._populate_tree()
                self.parent_window._update_plot_visibility()

            
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


class CreateMaterialDialog(QDialog):
    def __init__(self, next_mid, parent=None, existing_material=None):
        super().__init__(parent)
        self.existing_material = existing_material
        title = "Edit Material" if existing_material else "Create Material"
        self.setWindowTitle(title)
        self.setMinimumWidth(550)

        main_layout = QVBoxLayout(self)
        form_layout = QFormLayout()

        self.mid_input = QLineEdit(str(next_mid))
        self.mid_input.setValidator(QDoubleValidator())
        if existing_material:
            self.mid_input.setEnabled(False)
        self.title_input = QLineEdit()

        self.type_combo = QComboBox()
        self.type_combo.addItems(["MAT1 (Isotropic)", "MAT8 (Orthotropic 2D)", "MAT9 (Anisotropic 3D)"])
        
        form_layout.addRow("Material ID (MID):", self.mid_input)
        form_layout.addRow("Title:", self.title_input)
        form_layout.addRow("Material Type:", self.type_combo)
        main_layout.addLayout(form_layout)
        
        self.stacked_layout = QStackedLayout()
        self.mat1_widget = self._create_mat1_ui()
        self.mat8_widget = self._create_mat8_ui()
        self.mat9_widget = self._create_mat9_ui()
        self.stacked_layout.addWidget(self.mat1_widget)
        self.stacked_layout.addWidget(self.mat8_widget)
        self.stacked_layout.addWidget(self.mat9_widget)
        main_layout.addLayout(self.stacked_layout)
        
        self.type_combo.currentIndexChanged.connect(self.stacked_layout.setCurrentIndex)
        
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)
        
        if existing_material:
            self._populate_from_existing(existing_material)

    def _create_mat_qlineedit(self, default_val="0.0"):
        le = QLineEdit(default_val)
        le.setValidator(QDoubleValidator())
        return le

    def _create_mat1_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.mat1_e = self._create_mat_qlineedit("10.0E6")
        self.mat1_g = self._create_mat_qlineedit("3.8E6")
        self.mat1_nu = self._create_mat_qlineedit("0.33")
        self.mat1_rho = self._create_mat_qlineedit("0.1")
        self.mat1_a = self._create_mat_qlineedit("1.2E-5")
        self.mat1_tref = self._create_mat_qlineedit("70.0")
        self.mat1_ge = self._create_mat_qlineedit("0.0")
        layout.addRow("Young's Modulus (E):", self.mat1_e)
        layout.addRow("Shear Modulus (G):", self.mat1_g)
        layout.addRow("Poisson's Ratio (nu):", self.mat1_nu)
        layout.addRow("Density (rho):", self.mat1_rho)
        layout.addRow("Thermal Exp. (a):", self.mat1_a)
        layout.addRow("Ref. Temp (Tref):", self.mat1_tref)
        layout.addRow("Damping (ge):", self.mat1_ge)
        return widget

    def _create_mat8_ui(self):
        widget = QWidget()
        layout = QGridLayout(widget)
        self.mat8_e1 = self._create_mat_qlineedit(); self.mat8_e2 = self._create_mat_qlineedit()
        self.mat8_nu12 = self._create_mat_qlineedit(); self.mat8_g12 = self._create_mat_qlineedit()
        self.mat8_g1z = self._create_mat_qlineedit(); self.mat8_g2z = self._create_mat_qlineedit()
        self.mat8_rho = self._create_mat_qlineedit(); self.mat8_a1 = self._create_mat_qlineedit()
        self.mat8_a2 = self._create_mat_qlineedit(); self.mat8_tref = self._create_mat_qlineedit()
        layout.addWidget(QLabel("E1:"), 0, 0); layout.addWidget(self.mat8_e1, 0, 1)
        layout.addWidget(QLabel("E2:"), 0, 2); layout.addWidget(self.mat8_e2, 0, 3)
        layout.addWidget(QLabel("nu12:"), 1, 0); layout.addWidget(self.mat8_nu12, 1, 1)
        layout.addWidget(QLabel("G12:"), 1, 2); layout.addWidget(self.mat8_g12, 1, 3)
        layout.addWidget(QLabel("G1z:"), 2, 0); layout.addWidget(self.mat8_g1z, 2, 1)
        layout.addWidget(QLabel("G2z:"), 2, 2); layout.addWidget(self.mat8_g2z, 2, 3)
        layout.addWidget(QLabel("rho:"), 3, 0); layout.addWidget(self.mat8_rho, 3, 1)
        layout.addWidget(QLabel("a1:"), 3, 2); layout.addWidget(self.mat8_a1, 3, 3)
        layout.addWidget(QLabel("a2:"), 4, 0); layout.addWidget(self.mat8_a2, 4, 1)
        layout.addWidget(QLabel("Tref:"), 4, 2); layout.addWidget(self.mat8_tref, 4, 3)
        return widget

    def _create_mat9_ui(self):
        widget = QWidget()
        main_layout = QVBoxLayout(widget)
        main_layout.addWidget(QLabel("Anisotropic Stiffness Matrix [Gij]"))
        grid = QGridLayout()
        self.mat9_gij = {}
        labels = [
            "G11", "G12", "G13", "G14", "G15", "G16",
            "G22", "G23", "G24", "G25", "G26",
            "G33", "G34", "G35", "G36",
            "G44", "G45", "G46",
            "G55", "G56", "G66"
        ]
        row, col = 0, 0
        for label in labels:
            grid.addWidget(QLabel(f"{label}:"), row, col)
            le = self._create_mat_qlineedit()
            self.mat9_gij[label] = le
            grid.addWidget(le, row, col + 1)
            col += 2
            if col > 4:
                col = 0
                row += 1
        main_layout.addLayout(grid)
        
        other_props = QFormLayout()
        self.mat9_rho = self._create_mat_qlineedit(); self.mat9_tref = self._create_mat_qlineedit()
        other_props.addRow("rho:", self.mat9_rho); other_props.addRow("Tref:", self.mat9_tref)
        main_layout.addLayout(other_props)
        return widget
        
    def _populate_from_existing(self, mat):
        self.title_input.setText(mat.comment.strip().lstrip('$').strip())
        if mat.type == 'MAT1':
            self.type_combo.setCurrentIndex(0)
            self.mat1_e.setText(str(mat.e)); self.mat1_g.setText(str(mat.g)); self.mat1_nu.setText(str(mat.nu))
            self.mat1_rho.setText(str(mat.rho)); self.mat1_a.setText(str(mat.a)); self.mat1_tref.setText(str(mat.tref)); self.mat1_ge.setText(str(mat.ge))
        elif mat.type == 'MAT8':
            self.type_combo.setCurrentIndex(1)
            self.mat8_e1.setText(str(mat.e1)); self.mat8_e2.setText(str(mat.e2)); self.mat8_nu12.setText(str(mat.nu12))
            self.mat8_g12.setText(str(mat.g12)); self.mat8_g1z.setText(str(mat.g1z)); self.mat8_g2z.setText(str(mat.g2z))
            self.mat8_rho.setText(str(mat.rho)); self.mat8_a1.setText(str(mat.a1)); self.mat8_a2.setText(str(mat.a2)); self.mat8_tref.setText(str(mat.tref))
        elif mat.type == 'MAT9':
            self.type_combo.setCurrentIndex(2)
            for key, le in self.mat9_gij.items():
                le.setText(str(getattr(mat, key.lower(), 0.0)))
            self.mat9_rho.setText(str(mat.rho)); self.mat9_tref.setText(str(mat.tref))
        self.type_combo.setEnabled(False)

    def get_parameters(self):
        try:
            base_params = {
                'mid': int(self.mid_input.text()),
                'comment': f"$ {self.title_input.text()}"
            }
            if self.type_combo.currentIndex() == 0: # MAT1
                base_params['type'] = 'MAT1'
                base_params.update({
                    'E': float(self.mat1_e.text()), 'G': float(self.mat1_g.text()), 'nu': float(self.mat1_nu.text()),
                    'rho': float(self.mat1_rho.text()), 'a': float(self.mat1_a.text()), 'tref': float(self.mat1_tref.text()), 'ge': float(self.mat1_ge.text())
                })
            elif self.type_combo.currentIndex() == 1: # MAT8
                base_params['type'] = 'MAT8'
                base_params.update({
                    'E1': float(self.mat8_e1.text()), 'E2': float(self.mat8_e2.text()), 'nu12': float(self.mat8_nu12.text()),
                    'G12': float(self.mat8_g12.text()), 'G1z': float(self.mat8_g1z.text()), 'G2z': float(self.mat8_g2z.text()),
                    'rho': float(self.mat8_rho.text()), 'a1': float(self.mat8_a1.text()), 'a2': float(self.mat8_a2.text()), 'tref': float(self.mat8_tref.text())
                })
            elif self.type_combo.currentIndex() == 2: # MAT9
                base_params['type'] = 'MAT9'
                gij_dict = {key: float(le.text()) for key, le in self.mat9_gij.items()}
                base_params.update(gij_dict)
                base_params.update({
                    'rho': float(self.mat9_rho.text()), 'tref': float(self.mat9_tref.text())
                })
            return base_params
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Please provide valid numerical inputs. Error: {e}")
            return None

class CreatePropertyDialog(QDialog):
    def __init__(self, next_pid, model, parent=None, existing_prop=None):
        super().__init__(parent)
        self.model = model
        self.existing_prop = existing_prop
        title = "Edit Property" if existing_prop else "Create Property"
        self.setWindowTitle(title)
        self.setMinimumWidth(550)

        main_layout = QVBoxLayout(self)
        form_layout = QFormLayout()

        self.pid_input = QLineEdit(str(next_pid))
        self.pid_input.setValidator(QDoubleValidator())
        if existing_prop:
            self.pid_input.setEnabled(False)
        self.title_input = QLineEdit()

        self.type_combo = QComboBox()
        self.type_combo.addItems(["PSHELL", "PCOMP", "PBAR", "PBEAM"])
        
        form_layout.addRow("Property ID (PID):", self.pid_input)
        form_layout.addRow("Title:", self.title_input)
        form_layout.addRow("Property Type:", self.type_combo)
        main_layout.addLayout(form_layout)
        
        self.stacked_layout = QStackedLayout()
        self.pshell_widget = self._create_pshell_ui()
        self.pcomp_widget = self._create_pcomp_ui()
        self.pbar_widget = self._create_pbar_ui()
        self.pbeam_widget = self._create_pbeam_ui()

        self.stacked_layout.addWidget(self.pshell_widget)
        self.stacked_layout.addWidget(self.pcomp_widget)
        self.stacked_layout.addWidget(self.pbar_widget)
        self.stacked_layout.addWidget(self.pbeam_widget)
        main_layout.addLayout(self.stacked_layout)
        
        self.type_combo.currentIndexChanged.connect(self.stacked_layout.setCurrentIndex)
        
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)
        
        if existing_prop:
            self._populate_from_existing(existing_prop)
    
    def _create_prop_qlineedit(self, default_val="0.0"):
        le = QLineEdit(default_val)
        le.setValidator(QDoubleValidator())
        return le

    def _get_material_combo(self):
        combo = QComboBox()
        combo.addItems([str(mid) for mid in self.model.materials.keys()])
        return combo

    def _create_pshell_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.pshell_mid1 = self._get_material_combo()
        self.pshell_t = self._create_prop_qlineedit("0.1")
        self.pshell_nsm = self._create_prop_qlineedit("0.0")
        layout.addRow("Material ID (MID1):", self.pshell_mid1)
        layout.addRow("Thickness (t):", self.pshell_t)
        layout.addRow("Non-Structural Mass (nsm):", self.pshell_nsm)
        return widget

    def _create_pcomp_ui(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        form = QFormLayout()
        self.pcomp_nsm = self._create_prop_qlineedit("0.0")
        self.pcomp_ft = QComboBox(); self.pcomp_ft.addItems(['HILL', 'HOFF', 'TSAI', 'STRN'])
        form.addRow("Non-Structural Mass (nsm):", self.pcomp_nsm)
        form.addRow("Failure Theory (ft):", self.pcomp_ft)
        layout.addLayout(form)
        
        ply_group = QGroupBox("Ply Definition")
        ply_layout = QVBoxLayout(ply_group)
        self.pcomp_ply_table = QTableWidget(0, 4)
        self.pcomp_ply_table.setHorizontalHeaderLabels(["Material ID", "Thickness", "Theta", "SOUT"])
        self.pcomp_ply_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        ply_layout.addWidget(self.pcomp_ply_table)

        ply_button_layout = QHBoxLayout()
        add_ply_btn = QPushButton("Add Ply")
        remove_ply_btn = QPushButton("Remove Ply")
        ply_button_layout.addWidget(add_ply_btn)
        ply_button_layout.addWidget(remove_ply_btn)
        ply_layout.addLayout(ply_button_layout)
        
        add_ply_btn.clicked.connect(self._add_pcomp_ply)
        remove_ply_btn.clicked.connect(self._remove_pcomp_ply)

        layout.addWidget(ply_group)
        return widget

    def _add_pcomp_ply(self):
        row = self.pcomp_ply_table.rowCount()
        self.pcomp_ply_table.insertRow(row)
        
        mid_combo = self._get_material_combo()
        t_edit = self._create_prop_qlineedit("0.05")
        theta_edit = self._create_prop_qlineedit("0.0")
        sout_combo = QComboBox(); sout_combo.addItems(["NO", "YES"])

        self.pcomp_ply_table.setCellWidget(row, 0, mid_combo)
        self.pcomp_ply_table.setCellWidget(row, 1, t_edit)
        self.pcomp_ply_table.setCellWidget(row, 2, theta_edit)
        self.pcomp_ply_table.setCellWidget(row, 3, sout_combo)

    def _remove_pcomp_ply(self):
        current_row = self.pcomp_ply_table.currentRow()
        if current_row >= 0:
            self.pcomp_ply_table.removeRow(current_row)

    def _create_pbar_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.pbar_mid = self._get_material_combo()
        self.pbar_A = self._create_prop_qlineedit("1.0")
        self.pbar_i1 = self._create_prop_qlineedit("1.0")
        self.pbar_i2 = self._create_prop_qlineedit("1.0")
        self.pbar_j = self._create_prop_qlineedit("1.0")
        layout.addRow("Material ID (MID):", self.pbar_mid)
        layout.addRow("Area (A):", self.pbar_A)
        layout.addRow("Moment of Inertia (I1):", self.pbar_i1)
        layout.addRow("Moment of Inertia (I2):", self.pbar_i2)
        layout.addRow("Torsional Constant (J):", self.pbar_j)
        return widget

    def _create_pbeam_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.pbeam_mid = self._get_material_combo()
        self.pbeam_A = self._create_prop_qlineedit("1.0")
        self.pbeam_i1 = self._create_prop_qlineedit("1.0")
        self.pbeam_i2 = self._create_prop_qlineedit("1.0")
        self.pbeam_j = self._create_prop_qlineedit("1.0")
        layout.addRow("Material ID (MID):", self.pbeam_mid)
        layout.addRow("Area (A):", self.pbeam_A)
        layout.addRow("Moment of Inertia (I1):", self.pbeam_i1)
        layout.addRow("Moment of Inertia (I2):", self.pbeam_i2)
        layout.addRow("Torsional Constant (J):", self.pbeam_j)
        return widget

    def _populate_from_existing(self, prop):
        self.title_input.setText(prop.comment.strip().lstrip('$').strip())
        prop_type_map = {"PSHELL": 0, "PCOMP": 1, "PBAR": 2, "PBEAM": 3}
        self.type_combo.setCurrentIndex(prop_type_map.get(prop.type, 0))
        self.type_combo.setEnabled(False)
        
        if prop.type == 'PSHELL':
            self.pshell_mid1.setCurrentText(str(prop.mid1))
            self.pshell_t.setText(str(prop.t))
            self.pshell_nsm.setText(str(prop.nsm))
        elif prop.type == 'PCOMP':
            self.pcomp_nsm.setText(str(prop.nsm))
            self.pcomp_ft.setCurrentText(prop.ft)
            for i in range(prop.nplies):
                self._add_pcomp_ply()
                self.pcomp_ply_table.cellWidget(i, 0).setCurrentText(str(prop.mids[i]))
                self.pcomp_ply_table.cellWidget(i, 1).setText(str(prop.thicknesses[i]))
                self.pcomp_ply_table.cellWidget(i, 2).setText(str(prop.thetas[i]))
                self.pcomp_ply_table.cellWidget(i, 3).setCurrentText("YES" if prop.souts[i] == "YES" else "NO")
        elif prop.type == 'PBAR':
            self.pbar_mid.setCurrentText(str(prop.mid))
            self.pbar_A.setText(str(prop.A)); self.pbar_i1.setText(str(prop.i1));
            self.pbar_i2.setText(str(prop.i2)); self.pbar_j.setText(str(prop.j));
        elif prop.type == 'PBEAM':
            self.pbeam_mid.setCurrentText(str(prop.mid))
            self.pbeam_A.setText(str(prop.A[0])); self.pbeam_i1.setText(str(prop.i1[0]));
            self.pbeam_i2.setText(str(prop.i2[0])); self.pbeam_j.setText(str(prop.j[0]));

    def get_parameters(self):
        try:
            prop_type = self.type_combo.currentText()
            base_params = {
                'pid': int(self.pid_input.text()),
                'comment': f"$ {self.title_input.text()}",
                'type': prop_type
            }
            if prop_type == 'PSHELL':
                base_params.update({
                    'mid1': int(self.pshell_mid1.currentText()), 't': float(self.pshell_t.text()), 'nsm': float(self.pshell_nsm.text())
                })
            elif prop_type == 'PCOMP':
                plies = []
                for row in range(self.pcomp_ply_table.rowCount()):
                    mid = int(self.pcomp_ply_table.cellWidget(row, 0).currentText())
                    t = float(self.pcomp_ply_table.cellWidget(row, 1).text())
                    theta = float(self.pcomp_ply_table.cellWidget(row, 2).text())
                    sout = self.pcomp_ply_table.cellWidget(row, 3).currentText()
                    plies.append([mid, t, theta, sout])
                base_params.update({
                    'nsm': float(self.pcomp_nsm.text()), 'ft': self.pcomp_ft.currentText(), 'plies': plies
                })
            elif prop_type == 'PBAR':
                base_params.update({
                    'mid': int(self.pbar_mid.currentText()), 'A': float(self.pbar_A.text()),
                    'i1': float(self.pbar_i1.text()), 'i2': float(self.pbar_i2.text()), 'j': float(self.pbar_j.text())
                })
            elif prop_type == 'PBEAM':
                 base_params.update({
                    'mid': int(self.pbeam_mid.currentText()), 'A': float(self.pbeam_A.text()),
                    'i1': float(self.pbeam_i1.text()), 'i2': float(self.pbeam_i2.text()), 'j': float(self.pbeam_j.text()),
                })
            return base_params
        except (ValueError, KeyError, AttributeError) as e:
            QMessageBox.warning(self, "Input Error", f"Invalid input for property data. Error: {e}")
            return None



class ElementEditorDialog(QDialog):
    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Edit Element")
        self.setMinimumWidth(400)
        self.model = model
        self.current_eid = None
        self.node_inputs = []

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
            if eid not in self.model.elements:
                raise KeyError
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
            if child.widget():
                child.widget().deleteLater()

    def _populate_ui_for_element(self, element):
        self._clear_layout(self.props_layout)
        self.node_inputs = []

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

    def accept(self):
        if self.current_eid is None:
            super().reject() # Nothing to save
            return
        
        try:
            element = self.model.elements[self.current_eid]
            
            new_pid = int(self.pid_combo.currentText())
            element.pid = new_pid

            new_node_ids = [int(le.text()) for le in self.node_inputs]
            if len(new_node_ids) != len(element.nodes):
                raise ValueError("Number of nodes does not match element type.")
            
            element.nodes = new_node_ids

            super().accept()
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Failed to save element: {e}")




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



class EntitySelectionDialog(QDialog):
    request_show_selection = QtCore.Signal(str, list)
    request_picking_mode = QtCore.Signal(str, bool)

    def __init__(self, entity_type, all_entity_ids, parent=None, single_selection_mode=False):
        super().__init__(parent)
        self.entity_type = entity_type
        self.setWindowTitle(f"{self.entity_type} Selection")
        
        self.all_entity_ids, self.selected_ids = set(all_entity_ids), set()
        self.single_selection_mode = single_selection_mode
        
        self.list_widget = QListWidget()
        if self.single_selection_mode:
            self.list_widget.setSelectionMode(QListWidget.SingleSelection)
        self.list_widget.setMaximumHeight(100)

        main_layout = QVBoxLayout(self)
        panels_layout = QHBoxLayout()

        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.addWidget(QLabel(f"Selected {self.entity_type}s:"))

        selection_area_layout = QHBoxLayout()
        selection_area_layout.addWidget(self.list_widget)

        button_vbox_layout = QVBoxLayout()
        show_button, remove_button, select_all_button, clear_button = QPushButton("Show"), QPushButton("Remove"), QPushButton("Select All"), QPushButton("Clear All")
        button_vbox_layout.addWidget(show_button)
        button_vbox_layout.addWidget(remove_button)
        button_vbox_layout.addWidget(select_all_button)
        button_vbox_layout.addWidget(clear_button)
        button_vbox_layout.addStretch(1)
        selection_area_layout.addLayout(button_vbox_layout)
        left_layout.addLayout(selection_area_layout)

        pick_info_text = "Press ENTER to accept selection.\nPress ESC to cancel."
        pick_button = QPushButton("Pick from Viewer")
        pick_info = QLabel(pick_info_text)
        pick_info.setStyleSheet("font-size: 11px; color: #9aa3b2;")
        
        graphical_selection_layout = QHBoxLayout()
        graphical_selection_layout.addWidget(pick_button)
        graphical_selection_layout.addWidget(pick_info, 1)
        left_layout.addLayout(graphical_selection_layout)
        
        panels_layout.addWidget(left_widget, 2, QtCore.Qt.AlignTop)

        manual_group = QGroupBox("Manual Input")
        manual_form = QFormLayout(manual_group)
        self.single_id_input, add_single_button = QLineEdit(), QPushButton("Add")
        
        single_id_widget = QWidget()
        single_id_layout = QHBoxLayout(single_id_widget)
        single_id_layout.setContentsMargins(0,0,0,0)
        single_id_layout.addWidget(self.single_id_input)
        single_id_layout.addWidget(add_single_button)
        manual_form.addRow(f"{self.entity_type} ID:", single_id_widget)

        self.start_id_input, self.end_id_input, add_range_button = QLineEdit(), QLineEdit(), QPushButton("Add Range")
        manual_form.addRow("Start ID:", self.start_id_input)
        manual_form.addRow("End ID:", self.end_id_input)
        manual_form.addWidget(add_range_button)

        paste_group = QGroupBox("Paste List")
        paste_layout = QVBoxLayout(paste_group)
        self.paste_edit = QTextEdit()
        self.paste_edit.setPlaceholderText("Paste comma, space, or newline separated IDs")
        self.paste_edit.setMaximumHeight(100)
        add_list_button = QPushButton("Add from List")
        paste_layout.addWidget(self.paste_edit)
        paste_layout.addWidget(add_list_button)
        
        right_panel_widget = QWidget()
        right_panel_layout = QHBoxLayout(right_panel_widget)
        right_panel_layout.addWidget(manual_group)
        right_panel_layout.addWidget(paste_group, 0, QtCore.Qt.AlignTop)
        panels_layout.addWidget(right_panel_widget, 3)

        if self.single_selection_mode:
            select_all_button.setEnabled(False)
            self.start_id_input.setEnabled(False)
            self.end_id_input.setEnabled(False)
            add_range_button.setEnabled(False)
            self.paste_edit.setEnabled(False)
            add_list_button.setEnabled(False)
        
        main_layout.addLayout(panels_layout)
        main_layout.addStretch(1)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        show_button.clicked.connect(self._show_selection)
        pick_button.clicked.connect(self._activate_picking)
        add_single_button.clicked.connect(self._add_single_id)
        add_range_button.clicked.connect(self._add_range_ids)
        add_list_button.clicked.connect(self._add_list_ids)
        select_all_button.clicked.connect(self._select_all)
        remove_button.clicked.connect(self._remove_selected)
        clear_button.clicked.connect(self._clear_all)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        self.adjustSize()

        if self.parent():
            parent_geom = self.parent().frameGeometry()
            dialog_geom = self.frameGeometry()
            
            new_x = parent_geom.x() + (parent_geom.width() - dialog_geom.width()) / 2
            new_y = parent_geom.y() + (parent_geom.height() - dialog_geom.height())
            
            self.move(new_x, new_y)

    def add_selection(self, ids_to_add):
        if self.single_selection_mode:
            self.selected_ids.clear()
            if ids_to_add:
                ids_to_add = {list(ids_to_add)[0]}
        
        valid_ids = self.all_entity_ids.intersection(set(ids_to_add))
        self.selected_ids.update(valid_ids)
        self._refresh_list_widget()
        self._show_selection()
        if self.single_selection_mode and self.selected_ids:
            self.accept()

    def toggle_selection(self, id_to_toggle):
        if self.single_selection_mode:
            self.selected_ids.clear()
            if id_to_toggle in self.all_entity_ids:
                self.selected_ids.add(id_to_toggle)
            self._refresh_list_widget()
            self._show_selection()
            if self.selected_ids:
                self.accept()
            return

        if id_to_toggle in self.selected_ids:
            self.selected_ids.remove(id_to_toggle)
        elif id_to_toggle in self.all_entity_ids:
            self.selected_ids.add(id_to_toggle)
        self._refresh_list_widget()
        self._show_selection()

    def _select_all(self):
        self.selected_ids = self.all_entity_ids.copy()
        self._refresh_list_widget()
        self._show_selection()

    def _refresh_list_widget(self):
        self.list_widget.clear()
        self.list_widget.addItems([str(i) for i in sorted(list(self.selected_ids))])

    def _add_single_id(self):
        try:
            self.add_selection({int(self.single_id_input.text())})
            self.single_id_input.clear()
        except ValueError:
            pass

    def _add_range_ids(self):
        try:
            start_id, end_id = int(self.start_id_input.text()), int(self.end_id_input.text())
            self.add_selection(set(range(min(start_id, end_id), max(start_id, end_id) + 1)))
            self.start_id_input.clear()
            self.end_id_input.clear()
        except ValueError:
            pass

    def _add_list_ids(self):
        try:
            self.add_selection({int(i) for i in re.findall(r'\d+', self.paste_edit.toPlainText())})
            self.paste_edit.clear()
        except ValueError:
            pass

    def _remove_selected(self):
        for item in self.list_widget.selectedItems():
            self.selected_ids.remove(int(item.text()))
        self._refresh_list_widget()
        self._show_selection()

    def _clear_all(self):
        self.selected_ids.clear()
        self._refresh_list_widget()
        self.request_show_selection.emit(self.entity_type, [])

    def _show_selection(self):
        self.request_show_selection.emit(self.entity_type, self.get_selected_ids())

    def _activate_picking(self):
        self.hide()
        self.request_picking_mode.emit(self.entity_type, True)

    def get_selected_ids(self):
        return sorted(list(self.selected_ids))


class CreateNodesDialog(QDialog):
    def __init__(self, next_id, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Nodes"); self.setMinimumWidth(400); self.next_id = next_id
        main_layout = QVBoxLayout(self); self.tabs = QTabWidget()
        self.tabs.addTab(self._create_single_node_tab(), "Single Node"); self.tabs.addTab(self._create_nodes_between_tab(), "Nodes Between"); main_layout.addWidget(self.tabs)
        self.button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept); self.button_box.rejected.connect(self.reject); main_layout.addWidget(self.button_box)
    def _create_single_node_tab(self):
        widget = QWidget(); layout = QFormLayout(widget)
        self.node_id_input = QLineEdit(str(self.next_id)); self.node_x_input = QLineEdit("0.0"); self.node_y_input = QLineEdit("0.0"); self.node_z_input = QLineEdit("0.0")
        for w in [self.node_id_input, self.node_x_input, self.node_y_input, self.node_z_input]: w.setValidator(QDoubleValidator())
        layout.addRow("Node ID (0 for auto):", self.node_id_input); layout.addRow("X Coordinate:", self.node_x_input); layout.addRow("Y Coordinate:", self.node_y_input); layout.addRow("Z Coordinate:", self.node_z_input)
        return widget
    def _create_nodes_between_tab(self):
        widget = QWidget(); layout = QFormLayout(widget)
        self.start_nid_input = QLineEdit(); self.end_nid_input = QLineEdit(); self.num_nodes_input = QLineEdit("1")
        layout.addRow("Start Node ID:", self.start_nid_input); layout.addRow("End Node ID:", self.end_nid_input); layout.addRow("Number of Nodes to Add:", self.num_nodes_input)
        return widget
    def get_creation_parameters(self):
        params = {}
        try:
            if self.tabs.currentIndex() == 0:
                params.update({'type': 'single', 'id': int(float(self.node_id_input.text())), 'coords': [float(self.node_x_input.text()), float(self.node_y_input.text()), float(self.node_z_input.text())]})
            else:
                params.update({'type': 'between', 'start_nid': int(self.start_nid_input.text()), 'end_nid': int(self.end_nid_input.text()), 'num_nodes': int(self.num_nodes_input.text())})
        except ValueError: QMessageBox.warning(self, "Input Error", "Please provide valid numerical inputs."); return None
        return params

class CreateLineElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Line Element"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid)); self.type_combo = QComboBox(); self.type_combo.addItems(["CBEAM", "CBAR", "CROD"])
        self.pid_combo = QComboBox(); self.pid_combo.addItems([str(pid) for pid in model.properties.keys()])
        self.n1_input, self.n2_input = QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn = QPushButton("Pick..."), QPushButton("Pick...")
        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn)
        top_form.addRow("Element ID (0 for auto):", self.eid_input); top_form.addRow("Element Type:", self.type_combo); top_form.addRow("Property ID (PID):", self.pid_combo)
        top_form.addRow("Node 1 ID:", n1_widget); top_form.addRow("Node 2 ID:", n2_widget); main_layout.addLayout(top_form)
        self.orient_group = QGroupBox("Orientation"); orient_layout = QVBoxLayout(self.orient_group)
        self.orient_global_z_rb = QRadioButton("Global Z-Axis (0,0,1)", checked=True); self.orient_global_x_rb = QRadioButton("Global X-Axis (1,0,0)")
        self.orient_global_y_rb = QRadioButton("Global Y-Axis (0,1,0)"); self.orient_vector_rb = QRadioButton("Custom Vector"); self.orient_node_rb = QRadioButton("By Node")
        self.orient_vector_widget = QWidget()
        vec_form = QFormLayout(self.orient_vector_widget); vec_form.setContentsMargins(20, 5, 5, 5)
        self.base_x, self.base_y, self.base_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        self.dir_x, self.dir_y, self.dir_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        base_layout = QHBoxLayout(); base_layout.addWidget(self.base_x); base_layout.addWidget(self.base_y); base_layout.addWidget(self.base_z)
        dir_layout = QHBoxLayout(); dir_layout.addWidget(self.dir_x); dir_layout.addWidget(self.dir_y); dir_layout.addWidget(self.dir_z)
        vec_form.addRow("Base (X,Y,Z):", base_layout); vec_form.addRow("Direction (X,Y,Z):", dir_layout)
        self.orient_node_widget = QWidget()
        node_layout = QHBoxLayout(self.orient_node_widget); node_layout.setContentsMargins(20, 0, 0, 0)
        self.orient_node_id = QLineEdit(); pick_orient_node_btn = QPushButton("Pick...")
        node_layout.addWidget(QLabel("Node ID:")); node_layout.addWidget(self.orient_node_id); node_layout.addWidget(pick_orient_node_btn)
        orient_layout.addWidget(self.orient_global_z_rb); orient_layout.addWidget(self.orient_global_x_rb); orient_layout.addWidget(self.orient_global_y_rb)
        orient_layout.addWidget(self.orient_vector_rb); orient_layout.addWidget(self.orient_vector_widget)
        orient_layout.addWidget(self.orient_node_rb); orient_layout.addWidget(self.orient_node_widget)
        main_layout.addWidget(self.orient_group)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel); main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        self.type_combo.currentTextChanged.connect(self._update_ui_for_type); self.orient_vector_rb.toggled.connect(self._update_orientation_ui)
        self.orient_node_rb.toggled.connect(self._update_orientation_ui);
        pick_n1_btn.clicked.connect(self._pick_node_1); pick_n2_btn.clicked.connect(self._pick_node_2)
        pick_orient_node_btn.clicked.connect(self._pick_orientation_node)
        self.n1_input.textChanged.connect(self._update_highlight); self.n2_input.textChanged.connect(self._update_highlight)
        self._update_ui_for_type(self.type_combo.currentText()); self._update_orientation_ui()
    def _pick_node_1(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n1_input.setText)
    def _pick_node_2(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n2_input.setText)
    def _pick_orientation_node(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.orient_node_id.setText)
    def _update_highlight(self):
        nodes_to_highlight = []
        try: nodes_to_highlight.append(int(self.n1_input.text()))
        except ValueError: pass
        try: nodes_to_highlight.append(int(self.n2_input.text()))
        except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def _update_ui_for_type(self, text): self.orient_group.setEnabled(text != "CROD")
    def _update_orientation_ui(self):
        self.orient_vector_widget.setVisible(self.orient_vector_rb.isChecked()); self.orient_node_widget.setVisible(self.orient_node_rb.isChecked())
    def get_parameters(self):
        try:
            params = {'eid': int(self.eid_input.text()), 'pid': int(self.pid_combo.currentText()), 'n1': int(self.n1_input.text()), 'n2': int(self.n2_input.text()), 'type': self.type_combo.currentText()}
            orientation = {'method': 'default'}
            if self.orient_group.isEnabled():
                if self.orient_global_x_rb.isChecked(): orientation = {'method': 'vector', 'values': [1.0, 0.0, 0.0]}
                elif self.orient_global_y_rb.isChecked(): orientation = {'method': 'vector', 'values': [0.0, 1.0, 0.0]}
                elif self.orient_global_z_rb.isChecked(): orientation = {'method': 'vector', 'values': [0.0, 0.0, 1.0]}
                elif self.orient_vector_rb.isChecked():
                    base = np.array([float(self.base_x.text()), float(self.base_y.text()), float(self.base_z.text())])
                    direction = np.array([float(self.dir_x.text()), float(self.dir_y.text()), float(self.dir_z.text())])
                    orientation = {'method': 'vector', 'values': (direction - base).tolist()}
                elif self.orient_node_rb.isChecked(): orientation = {'method': 'node', 'values': [int(self.orient_node_id.text())]}
            params['orientation'] = orientation; return params
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None

class CreatePlateElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Plate Element"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid)); self.type_combo = QComboBox(); self.type_combo.addItems(["CQUAD4", "CTRIA3"])
        self.pid_combo = QComboBox(); self.pid_combo.addItems([str(pid) for pid in model.properties.keys()])
        self.n1_input, self.n2_input, self.n3_input, self.n4_input = QLineEdit(), QLineEdit(), QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn, pick_n3_btn, pick_n4_btn = QPushButton("Pick..."), QPushButton("Pick..."), QPushButton("Pick..."), QPushButton("Pick...")
        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn)
        n3_widget = QWidget(); n3_layout = QHBoxLayout(n3_widget); n3_layout.setContentsMargins(0,0,0,0); n3_layout.addWidget(self.n3_input); n3_layout.addWidget(pick_n3_btn)
        self.n4_widget = QWidget(); n4_layout = QHBoxLayout(self.n4_widget); n4_layout.setContentsMargins(0,0,0,0); n4_layout.addWidget(self.n4_input); n4_layout.addWidget(pick_n4_btn)
        top_form.addRow("Element ID (0 for auto):", self.eid_input); top_form.addRow("Element Type:", self.type_combo); top_form.addRow("Property ID (PID):", self.pid_combo)
        top_form.addRow("Node 1 ID:", n1_widget); top_form.addRow("Node 2 ID:", n2_widget); top_form.addRow("Node 3 ID:", n3_widget)
        self.n4_row_label = QLabel("Node 4 ID:"); top_form.addRow(self.n4_row_label, self.n4_widget); main_layout.addLayout(top_form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        self.type_combo.currentTextChanged.connect(self._toggle_n4)
        pick_n1_btn.clicked.connect(self._pick_node_1); pick_n2_btn.clicked.connect(self._pick_node_2)
        pick_n3_btn.clicked.connect(self._pick_node_3); pick_n4_btn.clicked.connect(self._pick_node_4)
        self.n1_input.textChanged.connect(self._update_highlight); self.n2_input.textChanged.connect(self._update_highlight)
        self.n3_input.textChanged.connect(self._update_highlight); self.n4_input.textChanged.connect(self._update_highlight)
        self._toggle_n4(self.type_combo.currentText())
    def _pick_node_1(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n1_input.setText)
    def _pick_node_2(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n2_input.setText)
    def _pick_node_3(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n3_input.setText)
    def _pick_node_4(self):
        if self.parent_window: self.parent_window._activate_single_node_picker(self.n4_input.setText)
    def _toggle_n4(self, text):
        is_quad = (text == "CQUAD4")
        self.n4_widget.setVisible(is_quad); self.n4_row_label.setVisible(is_quad)
        self._update_highlight()
    def _update_highlight(self):
        nodes_to_highlight = []
        try: nodes_to_highlight.append(int(self.n1_input.text()))
        except ValueError: pass
        try: nodes_to_highlight.append(int(self.n2_input.text()))
        except ValueError: pass
        try: nodes_to_highlight.append(int(self.n3_input.text()))
        except ValueError: pass
        if self.type_combo.currentText() == "CQUAD4":
            try: nodes_to_highlight.append(int(self.n4_input.text()))
            except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def get_parameters(self):
        try:
            nodes = [int(self.n1_input.text()), int(self.n2_input.text()), int(self.n3_input.text())]
            if self.type_combo.currentText() == "CQUAD4": nodes.append(int(self.n4_input.text()))
            return {'eid': int(self.eid_input.text()), 'pid': int(self.pid_combo.currentText()), 'nodes': nodes, 'type': self.type_combo.currentText()}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None

class CreateRbeDialog(QDialog):
    def __init__(self, next_eid, all_node_ids, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.all_node_ids = all_node_ids; self.setWindowTitle("Create RBE Element"); self.setMinimumSize(500, 400)
        main_layout = QVBoxLayout(self); form_layout = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid)); self.type_combo = QComboBox(); self.type_combo.addItems(["RBE2", "RBE3"])
        self.dof_input = QLineEdit("123456"); form_layout.addRow("Element ID (0 for auto):", self.eid_input); form_layout.addRow("Element Type:", self.type_combo)
        form_layout.addRow("DOF (CM):", self.dof_input); main_layout.addLayout(form_layout); splitter = QSplitter(QtCore.Qt.Horizontal)
        indep_box = QGroupBox("Independent Node(s)"); indep_layout = QVBoxLayout(indep_box); self.indep_list = QListWidget()
        self.indep_list.setSelectionMode(QListWidget.ExtendedSelection); select_indep_btn = QPushButton("Select...")
        indep_layout.addWidget(self.indep_list); indep_layout.addWidget(select_indep_btn)
        dep_box = QGroupBox("Dependent Node(s)"); dep_layout = QVBoxLayout(dep_box); self.dep_list = QListWidget()
        self.dep_list.setSelectionMode(QListWidget.ExtendedSelection); select_dep_btn = QPushButton("Select...")
        dep_layout.addWidget(self.dep_list); dep_layout.addWidget(select_dep_btn)
        splitter.addWidget(indep_box); splitter.addWidget(dep_box); main_layout.addWidget(splitter); button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box); select_indep_btn.clicked.connect(lambda: self.select_nodes_for_list(self.indep_list))
        select_dep_btn.clicked.connect(lambda: self.select_nodes_for_list(self.dep_list)); button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)

    def select_nodes_for_list(self, target_list):
        dialog = EntitySelectionDialog('Node', self.all_node_ids, self)
        dialog.request_show_selection.connect(self.parent_window._highlight_entities)
        dialog.request_picking_mode.connect(self.parent_window._set_picking_mode)
        if dialog.exec():
            target_list.clear()
            target_list.addItems([str(nid) for nid in dialog.get_selected_ids()])
        dialog.request_show_selection.disconnect(self.parent_window._highlight_entities)
        dialog.request_picking_mode.disconnect(self.parent_window._set_picking_mode)
        self.parent_window._highlight_entities('Node', [])
    def get_parameters(self):
        try:
            indep_nodes = [int(self.indep_list.item(i).text()) for i in range(self.indep_list.count())]; dep_nodes = [int(self.dep_list.item(i).text()) for i in range(self.dep_list.count())]
            if not indep_nodes or not dep_nodes: QMessageBox.warning(self, "Input Error", "Must select at least one independent and one dependent node."); return None
            return {'eid': int(self.eid_input.text()), 'type': self.type_combo.currentText(), 'dof': self.dof_input.text(), 'indep': indep_nodes, 'dep': dep_nodes}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None

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

class CreateLoadDialog(QDialog):
    selection_requested = QtCore.Signal(str)

    def __init__(self, model, parent=None, existing_sid=None):
        super().__init__(parent)
        self.model = model
        
        title = f"Edit Load (SID: {existing_sid})" if existing_sid else "Create Load"
        self.setWindowTitle(title)
        self.setMinimumWidth(500)

        main_layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._create_nodal_tab(), "Nodal (Force/Moment)")
        self.tabs.addTab(self._create_pressure_tab(), "Elemental (Pressure)")
        self.tabs.addTab(self._create_temp_tab(), "Body (Temperature)")
        
        main_layout.addWidget(self.tabs)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        
        if existing_sid:
            self._populate_from_existing(existing_sid)

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

            self.pres_sid.setText(str(sid))
            self.pres_sid.setEnabled(False)
            self.pres_val.setText(str(first_card.pressure))

            all_elems = sorted([card.eid for card in load_cards])
            self.pres_elem_list.addItems([str(eid) for eid in all_elems])

        elif load_type == 'TEMPD':
            self.tabs.setCurrentWidget(self.temp_tab)
            self.nodal_tab.setEnabled(False)
            self.pres_tab.setEnabled(False)

            self.temp_sid.setText(str(sid))
            self.temp_sid.setEnabled(False)
            self.temp_val.setText(str(first_card.temperature))

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



class ColorManagerDialog(QDialog):
    def __init__(self, color_map, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Color Manager"); self.setMinimumWidth(300); self.color_map, self.color_buttons = color_map.copy(), {}
        main_layout = QVBoxLayout(self); scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll_content = QWidget(); self.grid_layout = QGridLayout(scroll_content); scroll.setWidget(scroll_content)
        for row, (name, color) in enumerate(self.color_map.items()):
            label, btn = QLabel(str(name)), QPushButton(); btn.setFixedSize(24, 24); btn.setStyleSheet(f"background-color: {color}; border: 1px solid #555;")
            btn.clicked.connect(lambda c=False, n=name, b=btn: self.pick_color(n, b)); self.color_buttons[name] = btn; self.grid_layout.addWidget(label, row, 0); self.grid_layout.addWidget(btn, row, 1)
        main_layout.addWidget(scroll); button_box = QHBoxLayout(); randomize_btn, ok_btn, cancel_btn = QPushButton("Randomize"), QPushButton("OK"), QPushButton("Cancel")
        button_box.addWidget(randomize_btn); button_box.addStretch(); button_box.addWidget(ok_btn); button_box.addWidget(cancel_btn); main_layout.addLayout(button_box)
        randomize_btn.clicked.connect(self.randomize_colors); ok_btn.clicked.connect(self.accept); cancel_btn.clicked.connect(self.reject)
    def pick_color(self, name, btn):
        if (color := QColorDialog.getColor(QtGui.QColor(self.color_map[name]), self)).isValid(): self.color_map[name] = color.name(); btn.setStyleSheet(f"background-color: {self.color_map[name]}; border: 1px solid #555;")
    def randomize_colors(self):
        for name, btn in self.color_buttons.items():
            color = QtGui.QColor(random.randint(50, 220), random.randint(50, 220), random.randint(50, 220)); self.color_map[name] = color.name(); btn.setStyleSheet(f"background-color: {self.color_map[name]}; border: 1px solid #555;")

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
        img_path = os.path.join(os.path.dirname(__file__), "images", f"{combo_box.currentText()}_section.png")
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

class ClickAndDragInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, main_window):
        super().__init__()
        self.main_window = main_window
        self.mouse_press_time = 0
        self.mouse_press_pos = (0, 0)
        self.click_duration_threshold = 0.2
        self.drag_distance_threshold = 5
        self.AddObserver(vtk.vtkCommand.LeftButtonPressEvent, self.left_button_press_event)
        self.AddObserver(vtk.vtkCommand.LeftButtonReleaseEvent, self.left_button_release_event)
        self.AddObserver(vtk.vtkCommand.KeyPressEvent, self.key_press_event)
    
    def key_press_event(self, obj, event):
        key = self.GetInteractor().GetKeySym()
        if key == "Return":
            self.main_window._request_accept_picking_mode()
        elif key == "Escape":
            self.main_window._request_disable_picking_mode()
        else:
            self.OnKeyPress()

    def left_button_press_event(self, obj, event): 
        self.mouse_press_time = time.time()
        self.mouse_press_pos = self.GetInteractor().GetEventPosition()
        self.OnLeftButtonDown()

    def left_button_release_event(self, obj, event):
        is_click = (time.time() - self.mouse_press_time < self.click_duration_threshold and 
                    np.linalg.norm(np.array(self.GetInteractor().GetEventPosition()) - self.mouse_press_pos) < self.drag_distance_threshold)
        
        self.OnLeftButtonUp()

        if is_click:
            self.perform_pick(self.GetInteractor().GetEventPosition())

    def perform_pick(self, click_pos):
        print(f"DEBUG: Performing pick at screen position {click_pos}")
        if not (renderer := self.main_window.plotter.renderer): return

        entity_type = self.main_window.current_selection_type
        if not entity_type:
            entity_type = "Node"

        try:
            entity_id = -1
            if entity_type == 'Node':
                picker = vtk.vtkPointPicker()
                picker.Pick(click_pos[0], click_pos[1], 0, renderer)

                point_id_local = picker.GetPointId()
                picked_actor = picker.GetActor()

                if picked_actor and point_id_local != -1:
                    dataset = picked_actor.GetMapper().GetInput()
                    if dataset and 'vtkOriginalPointIds' in dataset.point_data:
                        original_point_id = int(dataset.point_data['vtkOriginalPointIds'][point_id_local])
                        if 0 <= original_point_id < len(self.main_window.current_node_ids_sorted):
                            entity_id = self.main_window.current_node_ids_sorted[original_point_id]

            elif entity_type == 'Element':
                picker = vtk.vtkCellPicker()
                picker.Pick(click_pos[0], click_pos[1], 0, renderer)
                cell_id_local = picker.GetCellId()
                picked_actor = picker.GetActor()

                if picked_actor and cell_id_local != -1:
                    dataset = picked_actor.GetMapper().GetInput()
                    if dataset and 'EID' in dataset.cell_data:
                        entity_id = int(dataset.cell_data["EID"][cell_id_local])

            if entity_id == -1:
                 self.main_window._update_status(f"No {entity_type.lower()} found at this location.")
                 return

            if self.main_window.picking_target_callback:
                self.main_window.picking_target_callback(str(entity_id))
                self.main_window._set_picking_mode("Node", False)
                self.main_window.picking_target_callback = None
                return

            if self.main_window.active_selection_dialog:
                dialog = self.main_window.active_selection_dialog
                dialog.toggle_selection(entity_id)
                status_verb = 'Selected' if entity_id in dialog.selected_ids else 'Deselected'
                self.main_window._update_status(f"{status_verb} {entity_type}: {entity_id}")
                self.main_window._highlight_entities(dialog.entity_type, dialog.get_selected_ids())

        except (IndexError, TypeError, KeyError) as e:
            self.main_window._update_status(f"Could not select {entity_type.lower()}: {e}", is_error=True)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Node Runner v1.0"); self.setGeometry(100, 100, 1200, 800)
        self.is_dark_theme, self.current_generator, self.current_grid = True, None, None
        self.shell_opacity, self.color_mode, self.render_style = 1.0, "property", "surface"
        self.load_scaling_info = {}
        self.type_color_map, self.pid_color_map = {"Shells": "#0077be", "Beams": "#f85a40"}, {}
        self.active_selection_dialog, self.current_node_ids_sorted = None, []; self.default_interactor_style = None
        self.picking_target_callback = None
        self.active_creation_dialog = None
        self.current_selection_type = None
        self.initUI(); self._auto_load_materials()
        self.tree_widget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tree_widget.customContextMenuRequested.connect(self._show_tree_context_menu)

    def initUI(self):
        self._create_menu_bar()
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        left_panel_container = QWidget()
        left_panel_layout = QVBoxLayout(left_panel_container)
        self.tree_widget = QTreeWidget()
        self.tree_widget.setHeaderLabels(["Model Tree"])
        left_panel_layout.addWidget(self.tree_widget)
        display_group = QGroupBox("Display Options")
        display_layout = QFormLayout(display_group)
        self.show_node_labels_check, self.show_elem_labels_check = QCheckBox(), QCheckBox()
        display_layout.addRow("Show Node IDs:", self.show_node_labels_check)
        display_layout.addRow("Show Element IDs:", self.show_elem_labels_check)

        # --- NEW: Add arrow scaling controls ---
        self.arrow_scale_input = QLineEdit("10.0")
        self.arrow_scale_input.setValidator(QDoubleValidator(0.1, 100.0, 2, self))
        display_layout.addRow("Arrow Size (% of Model):", self.arrow_scale_input)

        self.relative_scaling_check = QCheckBox()
        self.relative_scaling_check.setChecked(True)
        display_layout.addRow("Use Relative Scaling:", self.relative_scaling_check)

        left_panel_layout.addWidget(display_group)
        entity_group = QGroupBox("View Entity")
        entity_layout = QFormLayout(entity_group)
        self.entity_type_combo = QComboBox()
        self.entity_type_combo.addItems(["Node", "Element"])
        self.entity_id_input = QLineEdit()
        find_button, clear_button = QPushButton("Find"), QPushButton("Clear")
        button_layout = QHBoxLayout()
        button_layout.addWidget(find_button)
        button_layout.addWidget(clear_button)
        entity_layout.addRow("Entity Type:", self.entity_type_combo)
        entity_layout.addRow("Entity ID:", self.entity_id_input)
        entity_layout.addRow(button_layout)
        left_panel_layout.addWidget(entity_group)

        self.plotter = QtInteractor(self)
        self.plotter.set_background('#353535')
        
        self.plotter.add_axes()
        self.axes_actor = self.plotter.renderer.axes_actor
        self._set_axes_label_color((1, 1, 1))

        splitter = QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_panel_container)
        splitter.addWidget(self.plotter.interactor)
        splitter.setSizes([300, 900])
        main_layout.addWidget(splitter)
        self.statusBar().showMessage("Ready.")
        self.theme_button.clicked.connect(self._toggle_theme)
        self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)
        self.show_node_labels_check.stateChanged.connect(self._toggle_node_labels)
        self.show_elem_labels_check.stateChanged.connect(self._toggle_element_labels)

        # --- NEW: Connect scaling controls to the redraw function ---
        self.arrow_scale_input.textChanged.connect(self._update_plot_visibility)
        self.relative_scaling_check.stateChanged.connect(self._update_plot_visibility)

        find_button.clicked.connect(self._find_and_zoom_to_entity)
        clear_button.clicked.connect(self._clear_entity_highlight)

    def _set_axes_label_color(self, color_tuple):
        if not hasattr(self, 'axes_actor') or not self.axes_actor:
            return
        self.axes_actor.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)
        self.axes_actor.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)
        self.axes_actor.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(color_tuple)

    def _create_menu_bar(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("&File")
        open_action = QAction("&Open BDF/DAT...", self)
        open_action.triggered.connect(self._open_file_dialog)
        file_menu.addAction(open_action)
        
        # --- NEW: Added Save BDF action ---
        save_action = QAction("&Save BDF File...", self)
        save_action.triggered.connect(self._save_file_dialog)
        file_menu.addAction(save_action)

        file_menu.addSeparator()
        save_conf_action = QAction("Save Fuselage Config...", self)
        save_conf_action.triggered.connect(self._save_configuration)
        file_menu.addAction(save_conf_action)
        load_conf_action = QAction("Load Fuselage Config...", self)
        load_conf_action.triggered.connect(self._load_configuration)
        file_menu.addAction(load_conf_action)
        file_menu.addSeparator()
        quit_action = QAction("&Quit", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        edit_menu = menu_bar.addMenu("&Edit")
        move_nodes_action = QAction("Transform Nodes...", self)
        move_nodes_action.triggered.connect(self._open_node_move_tool)
        edit_menu.addAction(move_nodes_action)
        edit_element_action = QAction("Element...", self)
        edit_element_action.triggered.connect(self._open_element_editor)
        edit_menu.addAction(edit_element_action)
        edit_menu.addSeparator()
        edit_props_action = QAction("Properties...", self)
        edit_props_action.triggered.connect(self._open_property_editor)
        edit_menu.addAction(edit_props_action)
        edit_mats_action = QAction("Materials...", self)
        edit_mats_action.triggered.connect(self._open_material_editor)
        edit_menu.addAction(edit_mats_action)

        model_menu = menu_bar.addMenu("&Model")
        create_prop_action = QAction("Property...", self)
        create_prop_action.triggered.connect(self._create_property)
        model_menu.addAction(create_prop_action)
        create_mat_action = QAction("Material...", self)
        create_mat_action.triggered.connect(self._create_material)
        model_menu.addAction(create_mat_action)
        model_menu.addSeparator()
        create_load_action = QAction("Load...", self)
        create_load_action.triggered.connect(self._create_load)
        model_menu.addAction(create_load_action)
        create_constraint_action = QAction("Constraint...", self)
        create_constraint_action.triggered.connect(self._create_constraint)
        model_menu.addAction(create_constraint_action)

        mesh_menu = menu_bar.addMenu("M&esh")
        create_nodes_action = QAction("Nodes...", self)
        create_nodes_action.triggered.connect(self._create_nodes)
        mesh_menu.addAction(create_nodes_action)
        elements_menu = mesh_menu.addMenu("Elements")
        create_lines_action = QAction("Line...", self)
        create_lines_action.triggered.connect(self._create_line_elements)
        elements_menu.addAction(create_lines_action)
        create_plates_action = QAction("Plate...", self)
        create_plates_action.triggered.connect(self._create_plate_elements)
        elements_menu.addAction(create_plates_action)
        create_rbes_action = QAction("Rigid...", self)
        create_rbes_action.triggered.connect(self._create_rbes)
        elements_menu.addAction(create_rbes_action)

        delete_menu = menu_bar.addMenu("&Delete")
        delete_nodes_action = QAction("Nodes...", self)
        delete_nodes_action.triggered.connect(self._delete_nodes)
        delete_menu.addAction(delete_nodes_action)
        delete_elements_action = QAction("Elements...", self)
        delete_elements_action.triggered.connect(self._delete_elements)
        delete_menu.addAction(delete_elements_action)

        list_menu = menu_bar.addMenu("&List")
        list_node_action = QAction("Node Information...", self)
        list_node_action.triggered.connect(self._list_node_info)
        list_menu.addAction(list_node_action)
        list_element_action = QAction("Element Information...", self)
        list_element_action.triggered.connect(self._list_element_info)
        list_menu.addAction(list_element_action)

        tools_menu = menu_bar.addMenu("&Tools")
        fuselage_gen_action = QAction("Fuselage Generator...", self)
        fuselage_gen_action.triggered.connect(self._open_fuselage_generator)
        tools_menu.addAction(fuselage_gen_action)

        view_menu = menu_bar.addMenu("&View")
        self.show_origin_action = QAction("Show Origin", self, checkable=True)
        self.show_origin_action.toggled.connect(self._toggle_origin)
        view_menu.addAction(self.show_origin_action)
        self.shading_action = QAction("Shading", self, checkable=True, checked=True)
        self.shading_action.toggled.connect(self._toggle_shading)
        view_menu.addAction(self.shading_action)
        view_menu.addSeparator()
        views_menu = view_menu.addMenu("Standard Views")
        view_top = QAction("Top (+Z)", self); view_top.triggered.connect(lambda: self._set_standard_view('top')); views_menu.addAction(view_top)
        view_bottom = QAction("Bottom (-Z)", self); view_bottom.triggered.connect(lambda: self._set_standard_view('bottom')); views_menu.addAction(view_bottom)
        view_front = QAction("Front (+X)", self); view_front.triggered.connect(lambda: self._set_standard_view('front')); views_menu.addAction(view_front)
        view_back = QAction("Back (-X)", self); view_back.triggered.connect(lambda: self._set_standard_view('back')); views_menu.addAction(view_back)
        view_right = QAction("Right (+Y)", self); view_right.triggered.connect(lambda: self._set_standard_view('right')); views_menu.addAction(view_right)
        view_left = QAction("Left (-Y)", self); view_left.triggered.connect(lambda: self._set_standard_view('left')); views_menu.addAction(view_left)
        camera_menu = view_menu.addMenu("Camera")
        self.perspective_action = QAction("Perspective View", self, checkable=True, checked=True)
        self.perspective_action.toggled.connect(self._toggle_perspective_view)
        camera_menu.addAction(self.perspective_action)
        zoom_menu = view_menu.addMenu("Zoom")
        zoom_model_action = QAction("Zoom to Model", self)
        zoom_model_action.triggered.connect(self._zoom_to_model)
        zoom_menu.addAction(zoom_model_action)
        zoom_all_action = QAction("Zoom to All", self)
        zoom_all_action.triggered.connect(self._zoom_to_all)
        zoom_menu.addAction(zoom_all_action)
        style_menu = view_menu.addMenu("Style")
        style_group = QActionGroup(self)
        self.style_wire_action = QAction("Wireframe", self, checkable=True); self.style_wire_action.triggered.connect(lambda: self._set_render_style("wireframe")); style_group.addAction(self.style_wire_action)
        self.style_surf_action = QAction("Surface with Edges", self, checkable=True, checked=True); self.style_surf_action.triggered.connect(lambda: self._set_render_style("surface")); style_group.addAction(self.style_surf_action)
        style_menu.addAction(self.style_wire_action); style_menu.addAction(self.style_surf_action)
        color_menu = view_menu.addMenu("Color By")
        color_group = QActionGroup(self)
        color_type = QAction("Element Type", self, checkable=True); color_type.triggered.connect(lambda: self._set_coloring_mode("type")); color_group.addAction(color_type)
        color_prop = QAction("Property ID", self, checkable=True, checked=True); color_prop.triggered.connect(lambda: self._set_coloring_mode("property")); color_group.addAction(color_prop)
        color_menu.addAction(color_type); color_menu.addAction(color_prop)
        self.transparent_shells_action = QAction("Transparent Shells", self, checkable=True)
        self.transparent_shells_action.toggled.connect(self._toggle_shell_transparency)
        view_menu.addAction(self.transparent_shells_action)
        view_menu.addSeparator()
        edit_colors_action = QAction("Edit Colors...", self)
        edit_colors_action.triggered.connect(self._open_color_manager)
        view_menu.addAction(edit_colors_action)

        # --- NEW: Added Help menu ---
        help_menu = menu_bar.addMenu("&Help")
        about_action = QAction("&About...", self)
        about_action.triggered.connect(self._show_about_dialog)
        help_menu.addAction(about_action)
        license_action = QAction("&License...", self)
        license_action.triggered.connect(self._show_license_dialog)
        help_menu.addAction(license_action)

        self.theme_button = QPushButton("Switch to Light Mode")
        self.theme_button.setFlat(True)
        menu_bar.setCornerWidget(self.theme_button, QtCore.Qt.TopRightCorner)
    
    
# --- NEW: Handler for the "Save BDF File..." menu action ---
    def _save_file_dialog(self):
        if not self.current_generator:
            self._update_status("No model to save.", is_error=True)
            return

        filepath, _ = QFileDialog.getSaveFileName(self, "Save Nastran File", "", "Nastran Files (*.bdf *.dat);;All Files (*)")
        if not filepath:
            self._update_status("Save cancelled.")
            return
        
        try:
            self.current_generator.model.write_bdf(filepath, size=8, is_double=False)
            self._update_status(f"Model saved to {os.path.basename(filepath)}")
        except Exception as e:
            self._update_status(f"File save failed: {e}", is_error=True)
            QMessageBox.critical(self, "Error", f"Could not save file: {e}")

    # --- NEW: Handler for the "License..." menu action ---
    def _show_license_dialog(self):
        license_text = """
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent

      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS
        """
        dialog = QDialog(self)
        dialog.setWindowTitle("License Information")
        dialog.setMinimumSize(600, 500)
        layout = QVBoxLayout(dialog)
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setText(license_text)
        text_edit.setFontFamily("monospace")
        layout.addWidget(text_edit)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(dialog.accept)
        layout.addWidget(button_box)
        dialog.exec()

    # --- BUG FIX: Implement the missing method to add material cards ---
    def _add_material_card_from_params(self, params):
        """Adds a MAT1, MAT8, or MAT9 card to the model from a parameters dict."""
        params_copy = params.copy()
        mat_type = params_copy.pop('type')
        mid = params_copy.pop('mid')
        comment = params_copy.pop('comment')
        model = self.current_generator.model

        try:
            if mat_type == 'MAT1':
                model.add_mat1(mid, params_copy['E'], params_copy['G'], params_copy['nu'],
                               rho=params_copy.get('rho', 0.0), a=params_copy.get('a', 0.0),
                               tref=params_copy.get('tref', 0.0), ge=params_copy.get('ge', 0.0), comment=comment)
            elif mat_type == 'MAT8':
                model.add_mat8(mid, params_copy['E1'], params_copy['E2'], params_copy['nu12'],
                               g12=params_copy.get('G12', 0.0), g1z=params_copy.get('G1z', 0.0), g2z=params_copy.get('G2z', 0.0),
                               rho=params_copy.get('rho', 0.0), a1=params_copy.get('a1', 0.0), a2=params_copy.get('a2', 0.0),
                               tref=params_copy.get('tref', 0.0), comment=comment)
            elif mat_type == 'MAT9':
                g_values = {key.lower(): val for key, val in params_copy.items() if key.startswith('G')}
                model.add_mat9(mid,
                               rho=params_copy.get('rho', 0.0), 
                               tref=params_copy.get('tref', 0.0),
                               comment=comment, **g_values)
        except KeyError as e:
            QMessageBox.warning(self, "Creation Error", f"Could not create material. Missing parameter: {e}")    


    def _show_about_dialog(self):
        title = "About Node Runner"
        text = """
        <p><b>Node Runner v1.0</b></p>
        <p>Created by Angel Linares, September 2025</p>
        <hr>
        <p>A lightweight pre-processor for creating, editing, and visualizing Nastran models. 
        Includes a parametric fuselage generator and tools for general finite element modeling.</p>
        
        <p><b>Capabilities:</b></p>
        <p style="margin-left: 10px;">✅ Open, view, and <b>save</b> standard Nastran BDF files.</p>
        <p style="margin-left: 10px;">✅ Create & Edit: Nodes, Elements, Properties, and Materials.</p>
        <p style="margin-left: 10px;">✅ Create, Edit & Delete: <b>Loads</b> (Force, Moment, Pressure) and <b>Constraints</b> (SPC).</p>
        <p style="margin-left: 10px;">✅ Interactively select and transform nodes (move, scale, rotate).</p>
        <p style="margin-left: 10px;">✅ Visualize loads and constraints with an advanced, user-controlled scaling system.</p>
        <p style="margin-left: 10px;">✅ Includes a tool for parametric fuselage and floor structure generation.</p>

        <p><b>Current Limitations:</b></p>
        <p style="margin-left: 10px;">⚠️ Does not create <b>coordinate systems</b> (CORDs).</p>
        <p style="margin-left: 10px;">⚠️ The fuselage generator does not create tapered or non-cylindrical sections.</p>
        <p style="margin-left: 10px;">⚠️ Does not perform analysis or post-processing.</p>
        """
        QMessageBox.about(self, title, text)

    def _show_tree_context_menu(self, position):
        item = self.tree_widget.itemAt(position)
        if not item:
            return

        item_data = item.data(0, QtCore.Qt.UserRole)
        if not isinstance(item_data, tuple) or len(item_data) < 2:
            return

        entity_type, entity_id = item_data
        menu = QMenu()

        # Add "Edit..." action for applicable types
        if entity_type in ['property', 'material', 'load_set', 'constraint_set']:
            action_text = f"Edit {entity_type.replace('_', ' ').capitalize()} {entity_id}..."
            edit_action = QAction(action_text, self)
            
            if entity_type == 'property':
                edit_action.triggered.connect(self._open_property_editor)
            elif entity_type == 'material':
                edit_action.triggered.connect(self._open_material_editor)
            elif entity_type in ['load_set', 'constraint_set']:
                edit_action.triggered.connect(lambda: self._handle_edit_from_tree(entity_type, entity_id))
            menu.addAction(edit_action)

        # Add "Delete" action for applicable types
        if entity_type in ['load_set', 'constraint_set']:
            if menu.actions():
                menu.addSeparator()
            delete_action = QAction(f"Delete Set {entity_id}...", self)
            delete_action.triggered.connect(lambda: self._handle_delete_from_tree(entity_type, entity_id))
            menu.addAction(delete_action)

        # Only show the menu if actions were added
        if menu.actions():
            menu.exec(self.tree_widget.viewport().mapToGlobal(position))

    def _handle_delete_from_tree(self, entity_type, entity_id):
        if not self.current_generator: return
        type_name = "Load Set" if entity_type == 'load_set' else "Constraint Set"
        reply = QMessageBox.question(self, "Confirm Deletion", f"Are you sure you want to delete {type_name} {entity_id}?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            deleted_count = 0
            if entity_type == 'load_set':
                deleted_count = self.current_generator.delete_loads_by_sid(entity_id)
                # --- MODIFIED: Rebuild only the necessary actors after deletion ---
                if deleted_count > 0: self._create_all_load_actors()
            elif entity_type == 'constraint_set':
                deleted_count = self.current_generator.delete_constraints_by_sid(entity_id)
                if deleted_count > 0: self._create_all_constraint_actors()
            
            if deleted_count > 0:
                self._update_status(f"Deleted {type_name} {entity_id}.")
                self._populate_tree()
                self._update_plot_visibility()
        
    def _handle_edit_from_tree(self, entity_type, entity_id):
        """Launches the correct dialog in 'edit' mode for a given entity."""
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator: return
        
        dialog = None
        if entity_type == 'load_set':
            dialog = CreateLoadDialog(self.current_generator.model, self, existing_sid=entity_id)
            dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
            dialog.accepted.connect(self._on_load_creation_accept)
            dialog.rejected.connect(self._on_creation_reject)
        elif entity_type == 'constraint_set':
            dialog = CreateConstraintDialog(self.current_generator.model, self, existing_sid=entity_id)
            dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
            dialog.accepted.connect(self._on_constraint_creation_accept)
            dialog.rejected.connect(self._on_creation_reject)
        
        if dialog:
            dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            dialog.show()
            self.active_creation_dialog = dialog

    def _create_material(self):
        if not self.current_generator:
            self.current_generator = NastranModelGenerator()
        
        next_mid = max(self.current_generator.model.materials.keys()) + 1 if self.current_generator.model.materials else 1
        
        dialog = CreateMaterialDialog(next_mid, self)
        if dialog.exec():
            if params := dialog.get_parameters():
                mat_type = params['type']
                self._add_material_card_from_params(params)
                self._update_status(f"Created Material MID {params['mid']} ({mat_type}).")
                self._populate_tree()

    def _create_property(self):
        if not self.current_generator:
            self.current_generator = NastranModelGenerator()

        if not self.current_generator.model.materials:
            QMessageBox.warning(self, "No Materials", "You must create a material before creating a property.")
            return

        next_pid = max(self.current_generator.model.properties.keys()) + 1 if self.current_generator.model.properties else 1
        
        dialog = CreatePropertyDialog(next_pid, self.current_generator.model, self)
        if dialog.exec():
            if params := dialog.get_parameters():
                self._add_property_card_from_params(params)
                self._update_status(f"Created Property PID {params['pid']} ({params['type']}).")
                self._populate_tree()
                self._update_viewer(self.current_generator, reset_camera=False)

    def _create_coord(self):
        QMessageBox.information(self, "Not Implemented", "Coordinate system creation will be added in a future update.")

    def _add_property_card_from_params(self, params):
        params_copy = params.copy()
        prop_type = params_copy.pop('type')
        model = self.current_generator.model
        pid = params_copy.pop('pid')
        comment = params_copy.pop('comment')

        if prop_type == 'PSHELL':
            model.add_pshell(pid, mid1=params_copy['mid1'], t=params_copy['t'], nsm=params_copy['nsm'], comment=comment)
        elif prop_type == 'PCOMP':
            plies = params_copy['plies']
            if plies:
                mids, thicknesses, thetas, souts = zip(*plies)
                model.add_pcomp(pid, list(mids), list(thicknesses), list(thetas), souts=list(souts),
                                nsm=params_copy['nsm'], ft=params_copy['ft'], comment=comment)
            else:
                model.add_pcomp(pid, [], [], [], nsm=params_copy['nsm'], ft=params_copy['ft'], comment=comment)

        elif prop_type == 'PBAR':
            model.add_pbar(pid, params_copy['mid'], params_copy['A'], params_copy['i1'], params_copy['i2'], params_copy['j'], comment=comment)
        elif prop_type == 'PBEAM':
            model.add_pbeam(pid, params_copy['mid'], [0.0], ['C'], [params_copy['A']], [params_copy['i1']],
                            [params_copy['i2']], [0.0], [params_copy['j']], comment=comment)

    def _set_standard_view(self, view):
        if view == 'top': self.plotter.view_xy()
        elif view == 'bottom': self.plotter.view_xy(negative=True)
        elif view == 'front': self.plotter.view_yz()
        elif view == 'back': self.plotter.view_yz(negative=True)
        elif view == 'right': self.plotter.view_xz()
        elif view == 'left': self.plotter.view_xz(negative=True)
        self._update_status(f"View set to: {view.capitalize()}")

    def _toggle_origin(self, state):
        if state: self.plotter.add_mesh(pv.Sphere(center=(0,0,0), radius=0.01 * self.current_grid.length if self.current_grid else 0.05), color='cyan', name='origin_actor')
        else: self.plotter.remove_actor('origin_actor')

    def _toggle_shading(self, state):
        for actor in self.plotter.renderer.actors.values():
            if hasattr(actor, 'prop') and actor.prop is not None:
                actor.prop.lighting = state
                
        self._update_status(f"Shading {'ON' if state else 'OFF'}.")
        self.plotter.render()
        
    def _zoom_to_model(self):
        if self.current_grid: self.plotter.reset_camera(bounds=self.current_grid.bounds)
        
    def _zoom_to_all(self):
        if self.current_grid:
            # --- FIX: Manually calculate the expanded bounding box ---
            # Get the original bounds (xmin, xmax, ymin, ymax, zmin, zmax)
            bounds = np.array(self.current_grid.bounds).reshape(3, 2)
            
            # Calculate the center and size (extent) of the box
            center = np.mean(bounds, axis=1)
            extents = bounds[:, 1] - bounds[:, 0]
            
            # Increase the extents by a 1.1x factor
            new_extents = extents * 1.1
            
            # Calculate the new min/max points from the new center and extents
            new_bounds_arr = np.vstack([
                center - new_extents / 2,
                center + new_extents / 2
            ]).T.ravel()
            
            # Apply the new, expanded bounds to the camera
            self.plotter.reset_camera(bounds=new_bounds_arr)
        
    def _create_nodes(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator: self.current_generator = NastranModelGenerator()
        next_id = self.current_generator.get_next_available_id('node'); dialog = CreateNodesDialog(next_id, self)
        dialog.accepted.connect(self._on_node_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        
    def _on_node_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_creation_parameters()
        if params:
            if params['type'] == 'single':
                new_id = self.current_generator.add_node(params['id'], *params['coords'])
                self._update_status(f"Created Node: {new_id}")
            elif params['type'] == 'between':
                created_ids = self.current_generator.add_nodes_between(params['start_nid'], params['end_nid'], params['num_nodes'])
                self._update_status(f"Created {len(created_ids)} nodes.")
            self._update_viewer(self.current_generator, reset_camera=False)
        self.active_creation_dialog = None
        
    def _create_line_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateLineElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_line_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        
    def _on_line_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_parameters()
        if params:
            new_eid = self.current_generator.add_line_element(**params)
            self._update_status(f"Created {params['type']}: {new_eid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        self._highlight_nodes([])
        self.active_creation_dialog = None
        
    def _create_plate_elements(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "Create properties first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreatePlateElementDialog(next_eid, self.current_generator.model, self)
        dialog.accepted.connect(self._on_plate_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog
        

    def _on_plate_creation_accept(self):
        if not self.active_creation_dialog: return
        try:
            params = self.active_creation_dialog.get_parameters()
            if params:
                new_eid = self.current_generator.add_plate_element(**params)
                self._update_status(f"Created {params['type']}: {new_eid}")
                self._update_viewer(self.current_generator, reset_camera=False)
        finally:
            self._highlight_entities('Node', [])
            self.active_creation_dialog = None
        
    def _create_rbes(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator or not self.current_generator.model.nodes: QMessageBox.warning(self, "No Nodes", "Create nodes first."); return
        next_eid = self.current_generator.get_next_available_id('element'); dialog = CreateRbeDialog(next_eid, self.current_generator.model.nodes.keys(), self)
        dialog.accepted.connect(self._on_rbe_creation_accept); dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose); dialog.show(); self.active_creation_dialog = dialog

    def _delete_nodes(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "No nodes to delete.")
            return
        
        if self.active_selection_dialog:
            self.active_selection_dialog.activateWindow()
            return

        self.current_selection_type = 'Node'
        all_node_ids = self.current_generator.model.nodes.keys()
        
        self.active_selection_dialog = EntitySelectionDialog('Node', all_node_ids, self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)
        self.active_selection_dialog.accepted.connect(self._on_delete_nodes_accept)
        self.active_selection_dialog.rejected.connect(self._on_selection_dialog_reject)
        
        self.active_selection_dialog.show()

    def _on_delete_nodes_accept(self):
        dialog = self.active_selection_dialog
        if not dialog: return

        node_ids = dialog.get_selected_ids()
        
        self._on_selection_dialog_accept() 

        if not node_ids: return

        reply = QMessageBox.question(self, "Confirm Deletion",
                                   f"Are you sure you want to delete {len(node_ids)} selected node(s)? "
                                   "All connected elements will also be deleted.",
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            result = self.current_generator.delete_nodes(node_ids)
            self._update_status(f"Deleted {result['nodes']} nodes and {result['elements']} elements.")
            self._update_viewer(self.current_generator, reset_camera=False)

    def _delete_elements(self):
        if not self.current_generator or not self.current_grid or self.current_grid.n_cells == 0:
            QMessageBox.warning(self, "No Model", "No elements to delete.")
            return

        if self.active_selection_dialog:
            self.active_selection_dialog.activateWindow()
            return
            
        self.current_selection_type = 'Element'
        all_eids = self.current_grid.cell_data['EID']

        self.active_selection_dialog = EntitySelectionDialog('Element', all_eids, self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)
        self.active_selection_dialog.accepted.connect(self._on_delete_elements_accept)
        self.active_selection_dialog.rejected.connect(self._on_selection_dialog_reject)
        
        self.active_selection_dialog.show()

    def _on_delete_elements_accept(self):
        dialog = self.active_selection_dialog
        if not dialog: return

        eids_to_delete = dialog.get_selected_ids()
        
        self._on_selection_dialog_accept()

        if not eids_to_delete: return
        
        reply = QMessageBox.question(self, "Confirm Deletion",
                                f"Are you sure you want to delete {len(eids_to_delete)} element(s)?",
                                QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            deleted_count = self.current_generator.delete_elements(eids_to_delete)
            self._update_status(f"Deleted {deleted_count} elements.")
            self._update_viewer(self.current_generator, reset_camera=False)

    def _delete_mesh(self):
        if not self.current_generator: return
        reply = QMessageBox.question(self, "Confirm Delete Mesh",
                                   "Are you sure you want to delete all nodes and elements?",
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.current_generator.model = BDF(debug=False)
            self._update_viewer(self.current_generator)
            self._update_status("Mesh deleted.")

    def _list_node_info(self):
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "No nodes to list.")
            return
        if self.active_selection_dialog:
            return self.active_selection_dialog.activateWindow()

        self.current_selection_type = 'Node'
        all_node_ids = self.current_generator.model.nodes.keys()
        
        self.active_selection_dialog = EntitySelectionDialog('Node', all_node_ids, self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)
        self.active_selection_dialog.accepted.connect(self._on_list_nodes_accept)
        self.active_selection_dialog.rejected.connect(self._on_selection_dialog_reject)
        
        self.active_selection_dialog.show()

    def _list_element_info(self):
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Model", "No elements to list.")
            return
        if self.active_selection_dialog:
            return self.active_selection_dialog.activateWindow()
            
        all_elements = {**self.current_generator.model.elements, **self.current_generator.model.rigid_elements}
        if not all_elements:
            QMessageBox.warning(self, "No Model", "No elements to list.")
            return

        self.current_selection_type = 'Element'
        self.active_selection_dialog = EntitySelectionDialog('Element', all_elements.keys(), self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)
        self.active_selection_dialog.accepted.connect(self._on_list_elements_accept)
        self.active_selection_dialog.rejected.connect(self._on_selection_dialog_reject)

        self.active_selection_dialog.show()

    def _on_rbe_creation_accept(self):
        if not self.active_creation_dialog: return
        params = self.active_creation_dialog.get_parameters()
        if params:
            new_eid = self.current_generator.add_rbe_element(**params)
            self._update_status(f"Created {params['type']}: {new_eid}")
            self._update_viewer(self.current_generator, reset_camera=False)
        self.active_creation_dialog = None
        
    def _on_creation_reject(self):
        self._update_status("Creation cancelled.")
        self._highlight_entities('Node', [])
        self.active_creation_dialog = None

    def _open_node_move_tool(self):
        if self.active_selection_dialog:
            self.active_selection_dialog.activateWindow()
            return
        if not self.current_generator or not self.current_generator.model.nodes:
            QMessageBox.warning(self, "No Model", "Load a model first.")
            return
        
        self.current_selection_type = 'Node'
        all_node_ids = self.current_generator.model.nodes.keys()
        self.active_selection_dialog = EntitySelectionDialog('Node', all_node_ids, self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)
        self.active_selection_dialog.accepted.connect(self._on_selection_dialog_accept)
        self.active_selection_dialog.rejected.connect(self._on_selection_dialog_reject)
        self.active_selection_dialog.show()

    def _end_selection_mode(self):
        dialog = self.active_selection_dialog
        if dialog:
            try:
                dialog.request_show_selection.disconnect(self._highlight_entities)
                dialog.request_picking_mode.disconnect(self._set_picking_mode)
                dialog.accepted.disconnect()
                dialog.rejected.disconnect()
            except (TypeError, RuntimeError):
                pass
        self.active_selection_dialog = None

        if self.current_selection_type:
             self._highlight_entities(self.current_selection_type, [])
             self._set_picking_mode(self.current_selection_type, False)

        self.current_selection_type = None

    def _on_list_nodes_accept(self):
        dialog = self.active_selection_dialog
        if not dialog:
            return

        selected_ids = dialog.get_selected_ids()
        if selected_ids:
            info_content = []
            for nid in selected_ids:
                if nid in self.current_generator.model.nodes:
                    node = self.current_generator.model.nodes[nid]
                    info_content.append(f"--- Node {nid} ---\n{str(node)}")
            
            full_content = "\n\n".join(info_content)
            title = f"Information for {len(selected_ids)} Node(s)"
            info_dialog = InfoDialog(title, full_content, self)
            info_dialog.exec()
        
        self._end_selection_mode()

    def _on_list_elements_accept(self):
        dialog = self.active_selection_dialog
        if not dialog:
            return

        selected_ids = dialog.get_selected_ids()
        if selected_ids:
            all_elements = {**self.current_generator.model.elements, **self.current_generator.model.rigid_elements}
            info_content = []
            for eid in selected_ids:
                if eid in all_elements:
                    element = all_elements[eid]
                    info_content.append(f"--- Element {eid} ({element.type}) ---\n{str(element)}")
            
            full_content = "\n\n".join(info_content)
            title = f"Information for {len(selected_ids)} Element(s)"
            info_dialog = InfoDialog(title, full_content, self)
            info_dialog.exec()
        
        self._end_selection_mode()

    def _on_selection_dialog_accept(self):
        dialog = self.active_selection_dialog
        if not dialog:
            return

        is_transform_op = self.current_selection_type == 'Node' and dialog.sender() != self 
        
        if is_transform_op:
            selected_ids = dialog.get_selected_ids()
            if selected_ids:
              transform_dialog = NodeTransformDialog(self)
              if transform_dialog.exec() and (params := transform_dialog.get_transform_parameters()):
                  self.current_generator.transform_nodes(selected_ids, params)
                  self._update_viewer(self.current_generator)
                  self._update_status(f"Moved {len(selected_ids)} nodes.")
              else:
                  self._update_status("Node transform cancelled.")
        
        self._highlight_entities(self.current_selection_type, [])
        
        try:
            dialog.request_show_selection.disconnect(self._highlight_entities)
            dialog.request_picking_mode.disconnect(self._set_picking_mode)
        except (TypeError, RuntimeError): pass
        self.active_selection_dialog = None 
        
        self._set_picking_mode(self.current_selection_type, False)
        self.current_selection_type = None

    def _on_selection_dialog_reject(self):
        self._update_status("Selection cancelled.")
        self._end_selection_mode()

    def _open_element_editor(self):
        if not self.current_generator or not self.current_generator.model.elements:
            QMessageBox.warning(self, "No Elements", "The model contains no elements to edit.")
            return

        dialog = ElementEditorDialog(self.current_generator.model, self)
        if dialog.exec():
            self._update_status(f"Element {dialog.current_eid} modified.")
            self._update_viewer(self.current_generator)

    def _activate_single_node_picker(self, target_callback):
        self.picking_target_callback = target_callback
        self._set_picking_mode("Node", True)
        self._update_status("PICKING MODE: Click a single node.")

    def _highlight_entities(self, entity_type, entity_ids):
        camera_before = self.plotter.camera.copy()

        self.plotter.remove_actor('selection_highlight', render=False)

        if not entity_ids or not self.current_generator or not self.current_grid:
            self.plotter.render()
            return

        if entity_type == 'Node':
            coords = []
            for nid in entity_ids:
                if nid in self.current_generator.model.nodes:
                    coords.append(self.current_generator.model.nodes[nid].get_position())
            
            if coords:
                self.plotter.add_points(
                    np.array(coords),
                    color='yellow',
                    point_size=12,
                    render_points_as_spheres=True,
                    name='selection_highlight',
                )

        elif entity_type == 'Element':
            indices = np.isin(self.current_grid.cell_data['EID'], entity_ids)
            if np.any(indices):
                highlight_grid = self.current_grid.extract_cells(indices)
                self.plotter.add_mesh(
                    highlight_grid, style='wireframe', color='yellow',
                    line_width=5, name='selection_highlight', pickable=False
                )

        self.plotter.camera = camera_before
        self.plotter.render()

    def _set_picking_mode(self, entity_type, enabled):
        if enabled:
            if not self.default_interactor_style:
                self.default_interactor_style = self.plotter.interactor_style
            self._update_status(f"PICKING MODE: Left-click to select. ENTER to accept, ESC to cancel.")
            self.plotter.setCursor(QtCore.Qt.CrossCursor)
            self.plotter.add_text("PICKING MODE (ENTER to accept, ESC to cancel)", position='upper_right', color='yellow', font_size=8, name='picking_mode_text')
            self.plotter.interactor.SetInteractorStyle(ClickAndDragInteractor(main_window=self))
        else:
            self.plotter.setCursor(QtCore.Qt.ArrowCursor)
            self.plotter.remove_actor('picking_mode_text')
            if self.default_interactor_style:
                self.plotter.interactor.SetInteractorStyle(self.default_interactor_style)
            if self.active_selection_dialog and not self.active_selection_dialog.isVisible():
                self.active_selection_dialog.show()
                self.active_selection_dialog.activateWindow()
            self.picking_target_callback = None
            if self.current_selection_type is None:
                self._update_status("Picking mode disabled.")

    def _request_disable_picking_mode(self):
        QtCore.QTimer.singleShot(0, self._disable_picking_mode)

    def _request_accept_picking_mode(self):
        QtCore.QTimer.singleShot(0, self._accept_picking_mode)

    def _accept_picking_mode(self):
        if self.current_selection_type:
            self._set_picking_mode(self.current_selection_type, False)
        
        if self.active_selection_dialog:
            self.active_selection_dialog.show()
            self.active_selection_dialog.activateWindow()
            self._update_status("Selection complete. Click OK to proceed.")

    def _disable_picking_mode(self):
        if self.current_selection_type:
            self._set_picking_mode(self.current_selection_type, False)
        else:
            self._set_picking_mode('Node', False)
        
        if self.active_selection_dialog:
            self.active_selection_dialog.reject()

    def _update_viewer(self, generator, reset_camera=True):
        # --- MODIFIED: Explicitly remove all actors, including old load/constraint actors ---
        self.plotter.clear()
        self.plotter.add_axes()
        self.axes_actor = self.plotter.renderer.axes_actor
        self._set_axes_label_color((1, 1, 1) if self.is_dark_theme else (0, 0, 0))

        self.current_generator, self.current_grid = generator, None
        model = generator.model

        if not model.nodes:
            self.tree_widget.blockSignals(True)
            self._populate_tree()
            self.tree_widget.blockSignals(False)
            return

        self.current_node_ids_sorted = sorted(model.nodes.keys())
        node_map = {nid: i for i, nid in enumerate(self.current_node_ids_sorted)}
        node_coords = np.array([model.nodes[nid].get_position() for nid in self.current_node_ids_sorted])

        shells_conn, trias_conn, beams_conn, rods_conn = [], [], [], []
        elem_data = {'PID': [], 'EID': [], 'type': [], 'is_shell': []}

        all_elements = {**model.elements, **model.rigid_elements}
        for eid, elem in all_elements.items():
            try:
                if elem.type in ['CQUAD4', 'CMEMBRAN']:
                    shells_conn.append([4] + [node_map[nid] for nid in elem.nodes])
                    elem_data['PID'].append(elem.pid); elem_data['EID'].append(eid); elem_data['type'].append(elem.type); elem_data['is_shell'].append(1)
                elif elem.type == 'CTRIA3':
                    trias_conn.append([3] + [node_map[nid] for nid in elem.nodes])
                    elem_data['PID'].append(elem.pid); elem_data['EID'].append(eid); elem_data['type'].append(elem.type); elem_data['is_shell'].append(1)
                elif elem.type in ['CBEAM', 'CBAR']:
                    beams_conn.append([2] + [node_map[nid] for nid in elem.nodes])
                    elem_data['PID'].append(elem.pid); elem_data['EID'].append(eid); elem_data['type'].append(elem.type); elem_data['is_shell'].append(0)
                elif elem.type == 'CROD':
                    rods_conn.append([2] + [node_map[nid] for nid in elem.nodes])
                    elem_data['PID'].append(elem.pid); elem_data['EID'].append(eid); elem_data['type'].append(elem.type); elem_data['is_shell'].append(0)
            except KeyError as e:
                print(f"Warning: Skipping element {eid} due to missing node {e}.")

        all_conns = shells_conn + trias_conn + beams_conn + rods_conn
        if not all_conns:
            self.current_grid = pv.UnstructuredGrid(np.array([]), np.array([]), node_coords)
        else:
            flat_cells_list = [value for conn_list in all_conns for value in conn_list]
            cells = np.array(flat_cells_list)
            cell_types = np.concatenate([
                np.full(len(shells_conn), pv.CellType.QUAD), np.full(len(trias_conn), pv.CellType.TRIANGLE),
                np.full(len(beams_conn) + len(rods_conn), pv.CellType.LINE)
            ])
            self.current_grid = pv.UnstructuredGrid(cells, cell_types, node_coords)

        self.current_grid.point_data['vtkOriginalPointIds'] = np.arange(self.current_grid.n_points)
        if elem_data['EID']:
            for key, value in elem_data.items(): self.current_grid.cell_data[key] = np.array(value)

        for pid in model.properties.keys():
            if pid not in self.pid_color_map:
                self.pid_color_map[pid] = QColor(random.randint(50, 220), random.randint(50, 220), random.randint(50, 220)).name()

        self.tree_widget.blockSignals(True)
        self._populate_tree()
        self.tree_widget.blockSignals(False)

        # --- MODIFIED: Call new functions to build the actors ONCE per model load ---
        self._create_all_load_actors()
        self._create_all_constraint_actors()
        self._update_plot_visibility() # This now just toggles visibility

        if reset_camera and self.current_grid and self.current_grid.n_points > 0:
            self.plotter.reset_camera()
            self.plotter.view_isometric()
        else:
            self.plotter.render()

    def _find_tree_items(self, data_tuple):
        return [it.value() for it in QTreeWidgetItemIterator(self.tree_widget) if it.value().data(0, QtCore.Qt.UserRole) == data_tuple]
        
    def _open_property_editor(self):
        if not self.current_generator or not self.current_generator.model.properties: QMessageBox.warning(self, "No Properties", "No properties to edit."); return
        dialog = PropertyEditorDialog(self.current_generator.model, self)
        if dialog.exec(): self._update_plot_visibility(); self._populate_tree(); self._update_status("Properties updated.")
        
    def _open_material_editor(self):
        if not self.current_generator or not self.current_generator.model.materials: QMessageBox.warning(self, "No Materials", "No materials to edit."); return
        dialog = MaterialEditorDialog(self.current_generator.model, self)
        if dialog.exec(): self._update_status("Materials updated.")


    def _populate_tree(self):
        try:
            self.tree_widget.itemChanged.disconnect(self._handle_tree_item_changed)
        except (RuntimeError, TypeError): pass
        self.tree_widget.clear()

        if not self.current_generator:
            self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)
            return
        
        model = self.current_generator.model
        def create_item(parent, text, data=None, checked=True):
            item = QTreeWidgetItem(parent, [text])
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0, QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            item.setData(0, QtCore.Qt.UserRole, data)
            return item
        
        nodes_item = create_item(self.tree_widget, "Nodes", data=('group', 'nodes'))
        create_item(nodes_item, f"Total ({len(model.nodes)})", data=('nodes_total', 'all'))

        all_elements = {**model.elements, **model.rigid_elements}
        if all_elements:
            elems_item = create_item(self.tree_widget, f"Elements ({len(all_elements)})", data=('group', 'elements'))
            type_item = create_item(elems_item, "By Type", data=('group', 'elem_by_type'))
            shape_item = create_item(elems_item, "By Shape", data=('group', 'elem_by_shape'))

            by_type_map = {
                'CBEAM': 'Beams', 'CBAR': 'Bars', 'CROD': 'Rods',
                'CQUAD4': 'Plates', 'CMEMBRAN': 'Plates', 'CTRIA3': 'Plates',
                'RBE2': 'Rigid', 'RBE3': 'Rigid',
            }
            by_type_counts = Counter(by_type_map.get(e.type, 'Other') for e in all_elements.values())
            for type_name, count in sorted(by_type_counts.items()):
                create_item(type_item, f"{type_name} ({count})", data=('elem_by_type_group', type_name))

            shape_map = { 'Line': ['CBEAM', 'CBAR', 'CROD'], 'Quad': ['CQUAD4', 'CMEMBRAN'], 'Tria': ['CTRIA3'], 'Rigid': ['RBE2', 'RBE3'] }
            type_to_shape = {elem_type: shape for shape, types in shape_map.items() for elem_type in types}
            shape_counts = Counter(type_to_shape.get(e.type, 'Other') for e in all_elements.values())
            for shape, count in sorted(shape_counts.items()):
                create_item(shape_item, f"{shape} ({count})", data=('elem_shape_group', shape))

        if model.properties:
            props_item = create_item(self.tree_widget, "Properties", data=('group', 'properties'))
            for pid, p in sorted(model.properties.items()):
                title = p.comment.strip().lstrip('$').strip() or f'{p.type} {pid}'
                create_item(props_item, f"{pid}: {title}", data=('property', pid))
        
        if model.materials:
            mats_item = create_item(self.tree_widget, "Materials", data=('group', 'materials'))
            for mid, m in sorted(model.materials.items()):
                title = m.comment.strip().lstrip('$').strip() or f'Material {mid}'
                create_item(mats_item, f"{mid}: {title}", data=('material', mid))
        
        # --- MODIFICATION: Add Loads and Constraints to the tree ---
        if model.loads:
            loads_item = create_item(self.tree_widget, "Loads", data=('group', 'loads'))
            for sid, load_list in sorted(model.loads.items()):
                counts = Counter(load.type for load in load_list)
                summary = ", ".join(f"{count} {ltype}" for ltype, count in counts.items())
                create_item(loads_item, f"SID {sid}: {summary}", data=('load_set', sid))

        if model.spcs:
            constraints_item = create_item(self.tree_widget, "Constraints", data=('group', 'constraints'))
            for sid, spc_list in sorted(model.spcs.items()):
                # Assuming SPC1 or SPC, count unique nodes
                all_nodes = set()
                for spc in spc_list:
                    all_nodes.update(spc.nodes)
                summary = f"{len(all_nodes)} Nodes"
                create_item(constraints_item, f"SID {sid}: {spc_list[0].type} ({summary})", data=('constraint_set', sid))
        
        self.tree_widget.expandAll()
        self.tree_widget.itemChanged.connect(self._handle_tree_item_changed)




    def _handle_tree_item_changed(self, item, column):
        self.tree_widget.blockSignals(True)
        def set_child_states(parent_item):
            if parent_item.checkState(0) == QtCore.Qt.PartiallyChecked: return
            for i in range(parent_item.childCount()):
                child = parent_item.child(i); child.setCheckState(0, parent_item.checkState(0)); set_child_states(child)
        set_child_states(item)
        parent = item.parent()
        while parent:
            child_states = [parent.child(i).checkState(0) for i in range(parent.childCount())]
            if all(s == QtCore.Qt.Checked for s in child_states): parent.setCheckState(0, QtCore.Qt.Checked)
            elif all(s == QtCore.Qt.Unchecked for s in child_states): parent.setCheckState(0, QtCore.Qt.Unchecked)
            else: parent.setCheckState(0, QtCore.Qt.PartiallyChecked)
            parent = parent.parent()
        self.tree_widget.blockSignals(False); self._update_plot_visibility()
        

    def _open_fuselage_generator(self):
        mat_model = None
        if self.current_generator:
            mat_model = self.current_generator.model
        else:
            mat_model = self._auto_load_materials()

        if not mat_model or not mat_model.materials:
            QMessageBox.warning(self, "No Materials", "Load a model with materials or place 'materials.bdf' in the application folder."); return

        dialog = GeneratorDialog(self); dialog.set_materials(mat_model)
        if dialog.exec():
            try:
                params = dialog.get_parameters()
                generator = NastranModelGenerator(params=params)
                for mid, mat in mat_model.materials.items():
                    generator.model.add_mat1(mid, mat.e, mat.g, mat.nu, comment=mat.comment)
                
                generator._generate_nodes()
                generator._generate_elements()
                generator.generate_floor_structure()
                
                self._update_viewer(generator); self._update_status("Fuselage model generated.")
            except Exception as e: self._update_status(f"Generation failed: {e}", is_error=True); QMessageBox.critical(self, "Error", f"Failed to generate model: {e}")

    def _open_file_dialog(self):
        if not (filepath := QFileDialog.getOpenFileName(self, "Open Nastran File", "", "Nastran Files (*.bdf *.dat);;All Files (*)")[0]): return
        try:
            generator = NastranModelGenerator();
            with open(filepath, 'r') as f:
                generator.parse_from_text(f.read())
            
            self._update_viewer(generator);
            self._update_status(f"Displayed {os.path.basename(filepath)}")
        except Exception as e: self._update_status(f"File open failed.", is_error=True); QMessageBox.critical(self, "Error", f"Could not open/parse file: {e}")
    
    def _auto_load_materials(self):
        if os.path.exists("materials.bdf"):
            try:
                mat_model = BDF(debug=False)
                mat_model.read_bdf("materials.bdf", punch=True)
                self._update_status(f"Auto-loaded {len(mat_model.materials)} materials.")
                return mat_model
            except Exception as e:
                self._update_status(f"Could not parse default materials: {e}", is_error=True)
        else:
            self._update_status("'materials.bdf' not found.", is_error=True)
        return None

    def _save_configuration(self):
        if not (self.current_generator and self.current_generator.params): QMessageBox.warning(self, "No Data", "Only generator parameters can be saved."); return
        if not (save_path := QFileDialog.getSaveFileName(self, "Save Configuration", "", "JSON Files (*.json)")[0]): return
        try:
            with open(save_path, 'w') as f: json.dump(self.current_generator.params, f, indent=4)
            self._update_status(f"Config saved to {os.path.basename(save_path)}")
        except Exception as e: self._update_status(f"Error saving config: {e}", is_error=True)

    def _load_configuration(self):
        if not (load_path := QFileDialog.getOpenFileName(self, "Load Configuration", "", "JSON Files (*.json)")[0]): return
        try:
            with open(load_path, 'r') as f: config = json.load(f)
            dialog = GeneratorDialog(self);
            
            mat_model = self.current_generator.model if self.current_generator else self._auto_load_materials()
            if not mat_model:
                 QMessageBox.warning(self, "No Materials", "Could not find materials to use for generator."); return
            dialog.set_materials(mat_model)

            for key, w in [('fuselage_radius',dialog.radius_input),('num_bays',dialog.bays_input),('frame_spacing',dialog.spacing_input),('skin_thickness',dialog.skin_thick_input),('skin_material_id',dialog.skin_mat_combo),('num_stringers',dialog.num_stringers_input),('stringer_section_type',dialog.stringer_type_combo),('stringer_material_id',dialog.stringer_mat_combo),('frame_section_type',dialog.frame_type_combo),('frame_material_id',dialog.frame_mat_combo)]:
                if key in config: getattr(w, 'setText' if isinstance(w, QLineEdit) else 'setCurrentText')(str(config[key]))
            self._update_status(f"Config loaded from {os.path.basename(load_path)}.")
            if dialog.exec():
                params = dialog.get_parameters()
                generator = NastranModelGenerator(params=params)
                for mid, mat in mat_model.materials.items():
                    generator.model.add_mat1(mid, mat.e, mat.g, mat.nu, comment=mat.comment)
                
                generator._generate_nodes(); generator._generate_elements()
                self._update_viewer(generator); self._update_status("Model generated from loaded configuration.")
        except Exception as e: self._update_status(f"Error loading config: {e}", is_error=True); QMessageBox.critical(self, "Error", f"Could not load config: {e}")
        
    def _open_color_manager(self):
        if not self.current_grid: self._update_status("No model loaded.", is_error=True); return
        current_map = self.type_color_map if self.color_mode == "type" else self.pid_color_map; dialog = ColorManagerDialog(current_map, self)
        if dialog.exec():
            (self.type_color_map if self.color_mode == "type" else self.pid_color_map).update(dialog.color_map)
            self._update_plot_visibility(); self._update_status("Colors updated.")
            
    def _set_render_style(self, style): 
        self.render_style = style; self.transparent_shells_action.setEnabled(style=="surface"); self._update_plot_visibility(); self._refresh_labels(); self._update_status(f"Render style: {style}.")
        
    def _toggle_shell_transparency(self, state): 
        self.shell_opacity = 0.4 if state else 1.0; self._update_plot_visibility(); self._refresh_labels(); self._update_status(f"Shell transparency {'ON' if state else 'OFF'}.")
        
    def _toggle_perspective_view(self, state): 
        self.plotter.camera.SetParallelProjection(not state); self.plotter.render(); self._update_status(f"Perspective view {'ON' if state else 'OFF'}.")
        
    def _set_coloring_mode(self, mode): 
        self.color_mode = mode; self._update_plot_visibility(); self._update_status(f"Color mode: {mode}.")
        
    def _rebuild_plot(self, visibility_mask=None):
        self.plotter.remove_actor(['nodes_actor', 'shells', 'beams'])

        nodes_visible = False
        if (nodes_item_list := self._find_tree_items(('group', 'nodes'))):
            if nodes_item_list[0].checkState(0) == QtCore.Qt.Checked:
                nodes_visible = True

        if nodes_visible and self.current_grid and self.current_grid.n_points > 0:
            self.plotter.add_points(self.current_grid, color='lime', point_size=5,
                                    render_points_as_spheres=True, name='nodes_actor')

        grid_to_render = self.current_grid
        if visibility_mask is not None:
            visible_indices = np.where(visibility_mask)[0]
            if visible_indices.size > 0:
                grid_to_render = self.current_grid.extract_cells(visible_indices)
            else:
                self.plotter.remove_actor(['shells', 'beams'])
                self.plotter.render()
                return

        if grid_to_render and grid_to_render.n_cells > 0:
            shells = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 1)
            beams = grid_to_render.extract_cells(grid_to_render.cell_data['is_shell'] == 0)
            show_edges = self.render_style == "surface"

            if self.color_mode == "type":
                if shells.n_cells > 0:
                    self.plotter.add_mesh(shells, style=self.render_style, show_edges=show_edges,
                                          opacity=self.shell_opacity, color=self.type_color_map["Shells"], name="shells")
                if beams.n_cells > 0:
                    self.plotter.add_mesh(beams, color=self.type_color_map["Beams"],
                                          render_lines_as_tubes=True, line_width=2, name="beams")
            elif self.color_mode == "property":
                if shells.n_cells > 0:
                    pids = np.unique(shells.cell_data['PID'])
                    cmap = [self.pid_color_map.get(p, "#FFFFFF") for p in pids]
                    self.plotter.add_mesh(shells, style=self.render_style, show_edges=show_edges, scalars="PID",
                                          opacity=self.shell_opacity, cmap=cmap, categories=True, name="shells", show_scalar_bar=False)
                if beams.n_cells > 0:
                    pids = np.unique(beams.cell_data['PID'])
                    cmap = [self.pid_color_map.get(p, "#FFFFFF") for p in pids]
                    self.plotter.add_mesh(beams, scalars="PID", cmap=cmap, categories=True,
                                          render_lines_as_tubes=True, line_width=2, name="beams", show_scalar_bar=False)
        else:
             self.plotter.remove_actor(['shells', 'beams'])


        self.plotter.render()
        
    def _refresh_labels(self):
        if self.show_node_labels_check.isChecked(): self._toggle_node_labels(False); self._toggle_node_labels(True)
        if self.show_elem_labels_check.isChecked(): self._toggle_element_labels(False); self._toggle_element_labels(True)
        
    def _toggle_node_labels(self, state):
        self.plotter.remove_actor('node_labels', render=False)
        if state and self.current_grid:
            color = 'white' if self.is_dark_theme else 'black'
            occlude = self.render_style == 'surface' and self.shell_opacity == 1.0
            self.plotter.add_point_labels(self.current_grid.points, self.current_node_ids_sorted, name='node_labels', font_size=12, text_color=color, shape=None, always_visible=(not occlude), render=True)
        self._update_status(f"Node labels {'ON' if state else 'OFF'}")

    def _toggle_element_labels(self, state):
        self.plotter.remove_actor('elem_labels', render=False)
        if state and self.current_grid:
            color = 'white' if self.is_dark_theme else 'black'
            occlude = self.render_style == 'surface' and self.shell_opacity == 1.0
            self.plotter.add_point_labels(self.current_grid.cell_centers().points, self.current_grid.cell_data["EID"], name='elem_labels', font_size=12, text_color=color, shape=None, always_visible=(not occlude), render=True)
        self._update_status(f"Element labels {'ON' if state else 'OFF'}")
    
    def _find_and_zoom_to_entity(self):
        if not self.current_generator or not self.current_grid: self._update_status("No model.", is_error=True); return
        try: entity_id = int(self.entity_id_input.text()); etype = self.entity_type_combo.currentText()
        except ValueError: self._update_status("Invalid ID.", is_error=True); return
        
        coords, found = None, False
        model = self.current_generator.model

        if etype == "Node":
            if entity_id in model.nodes:
                coords = model.nodes[entity_id].get_position()
                found = True
        elif etype == "Element":
            if entity_id in self.current_grid.cell_data["EID"]:
                idx = np.where(self.current_grid.cell_data["EID"] == entity_id)[0][0]
                coords = self.current_grid.get_cell(idx).center
                found = True

        self.plotter.remove_actor('highlight_actor', render=False)
        if found:
            self.plotter.add_points(np.array(coords), color='yellow', point_size=15, render_points_as_spheres=True, name='highlight_actor'); self.plotter.fly_to(coords)
            self._update_status(f"Found {etype} {entity_id}.")
        else: self._update_status(f"{etype} {entity_id} not found.", is_error=True)

    def _clear_entity_highlight(self): 
        self.plotter.remove_actor('highlight_actor', render=True); self.entity_id_input.clear(); self._update_status("Highlight cleared.")
        
    def _update_status(self, msg, is_error=False): 
        self.statusBar().showMessage(msg); self.statusBar().setStyleSheet(f"color: {'#ff6b6b' if is_error else '#9aa3b2'};")
        
    def _toggle_theme(self):
        app = QApplication.instance()
        self.is_dark_theme = not self.is_dark_theme
        app.setPalette(dark_palette if self.is_dark_theme else light_palette)

        bg_color = '#353535' if self.is_dark_theme else 'white'
        label_color_tuple = (1, 1, 1) if self.is_dark_theme else (0, 0, 0)
        
        self.plotter.set_background(bg_color)
        
        self.plotter.remove_axes()
        self.plotter.add_axes()
        self.axes_actor = self.plotter.renderer.axes_actor
        self._set_axes_label_color(label_color_tuple)

        self.theme_button.setText(f"Switch to {'Light' if self.is_dark_theme else 'Dark'} Mode")
        self._update_plot_visibility()


    def _create_load(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please open or generate a model first.")
            return

        dialog = CreateLoadDialog(self.current_generator.model, self)
        dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
        dialog.accepted.connect(self._on_load_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _create_constraint(self):
        if self.active_creation_dialog: return self.active_creation_dialog.activateWindow()
        if not self.current_generator:
            QMessageBox.warning(self, "No Model", "Please open or generate a model first.")
            return

        dialog = CreateConstraintDialog(self.current_generator.model, self)
        # --- FIX: Connect directly to the handler, since the signal now carries the entity type ---
        dialog.selection_requested.connect(self._handle_sub_dialog_selection_request)
        dialog.accepted.connect(self._on_constraint_creation_accept)
        dialog.rejected.connect(self._on_creation_reject)
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.show()
        self.active_creation_dialog = dialog

    def _on_load_creation_accept(self):
        dialog = self.active_creation_dialog
        if not dialog: return
        
        params = dialog.get_parameters()
        if params:
            gen = self.current_generator
            sid = params['sid']
            
            if not dialog.tabs.currentWidget().findChild(QLineEdit).isEnabled():
                gen.delete_loads_by_sid(sid)
                self._update_status(f"Updating Load Set {sid}...")
            
            load_type = params.pop('type'); count = 0
            if load_type == 'nodal': count = gen.add_nodal_load(**params); self._update_status(f"Created/Updated {count} {params['load_type']} cards in SID {sid}.")
            elif load_type == 'pressure': count = gen.add_pressure_load(**params); self._update_status(f"Created/Updated PLOAD4 on {count} elements in SID {sid}.")
            elif load_type == 'temperature': count = gen.add_default_temperature(**params); self._update_status(f"Created/Updated TEMPD card in SID {sid}.")
            
            # --- MODIFIED: Rebuild tree, then rebuild all actors, then update visibility ---
            self._populate_tree()
            self._create_all_load_actors()
            self._update_plot_visibility()
        
        self.active_creation_dialog = None

    def _on_constraint_creation_accept(self):
        dialog = self.active_creation_dialog
        if not dialog: return
        
        params = dialog.get_parameters()
        if params:
            sid = params['sid']
            if not dialog.sid_input.isEnabled():
                self.current_generator.delete_constraints_by_sid(sid)
                self._update_status(f"Updating Constraint Set {sid}...")

            count = self.current_generator.add_nodal_constraint(**params)
            self._update_status(f"Created/Updated SPC on {count} nodes in SID {sid}.")
            
            # --- MODIFIED: Rebuild tree, then rebuild all actors, then update visibility ---
            self._populate_tree()
            self._create_all_constraint_actors()
            self._update_plot_visibility()

        self.active_creation_dialog = None


    # In class MainWindow, ADD this new method.
    def _handle_sub_dialog_selection_request(self, entity_type):
        if self.active_creation_dialog:
            self.active_creation_dialog.hide()

        if entity_type == 'Node':
            all_ids = self.current_generator.model.nodes.keys()
        else:
            all_ids = [eid for eid, elem in self.current_generator.model.elements.items() if elem.type in ['CQUAD4', 'CTRIA3']]
        
        self.current_selection_type = entity_type
        self.active_selection_dialog = EntitySelectionDialog(entity_type, all_ids, self)
        self.active_selection_dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.active_selection_dialog.request_show_selection.connect(self._highlight_entities)
        self.active_selection_dialog.request_picking_mode.connect(self._set_picking_mode)

        def on_accept():
            selected_ids = self.active_selection_dialog.get_selected_ids()
            if self.active_creation_dialog:
                self.active_creation_dialog.update_selection_list(entity_type, selected_ids)
                self.active_creation_dialog.show()
            self._end_selection_mode()

        def on_reject():
            if self.active_creation_dialog:
                self.active_creation_dialog.show()
            self._end_selection_mode()

        self.active_selection_dialog.accepted.connect(on_accept)
        self.active_selection_dialog.rejected.connect(on_reject)
        self.active_selection_dialog.show()

    def _create_all_load_actors(self):
        """
        Creates all load actors for the current model, one per SID.
        This is slow and should only be called when loads are created/edited/deleted.
        """
        actor_names_to_remove = [name for name in self.plotter.actors if name.startswith(('force_actors_', 'moment_actors_', 'pressure_actors_'))]
        self.plotter.remove_actor(actor_names_to_remove, render=False)

        # --- MODIFIED: Clear and prepare the scaling info dictionary ---
        self.load_scaling_info.clear()

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.loads: return
        
        all_force_mags, all_moment_mags, all_pressure_mags = [], [], []
        for load_list in model.loads.values():
            for load in load_list:
                if load.type == 'FORCE': all_force_mags.append(np.linalg.norm(np.array(load.xyz) * load.mag))
                elif load.type == 'MOMENT': all_moment_mags.append(np.linalg.norm(np.array(load.xyz) * load.mag))
                elif load.type == 'PLOAD4': all_pressure_mags.append(abs(load.pressure))

        # --- MODIFIED: Store max magnitudes in the new dictionary ---
        self.load_scaling_info['force_max_mag'] = np.max(all_force_mags) if all_force_mags else 0.0
        self.load_scaling_info['moment_max_mag'] = np.max(all_moment_mags) if all_moment_mags else 0.0
        self.load_scaling_info['pressure_max_mag'] = np.max(all_pressure_mags) if all_pressure_mags else 0.0
        
        eid_to_cell_idx_map = {eid: i for i, eid in enumerate(self.current_grid.cell_data['EID'])}

        for sid, load_list in model.loads.items():
            force_points, force_vectors = [], []
            moment_points, moment_vectors = [], []
            pressure_points, pressure_vectors = [], []

            for load in load_list:
                if load.type == 'FORCE' and load.node in model.nodes:
                    force_points.append(model.nodes[load.node].get_position())
                    force_vectors.append(np.array(load.xyz) * load.mag)
                elif load.type == 'MOMENT' and load.node in model.nodes:
                    moment_points.append(model.nodes[load.node].get_position())
                    moment_vectors.append(np.array(load.xyz) * load.mag)
                elif load.type == 'PLOAD4':
                    if (cell_idx := eid_to_cell_idx_map.get(load.eid)) is not None:
                        cell = self.current_grid.get_cell(cell_idx)
                        pressure_points.append(cell.center)
                        pressure_vectors.append(cell.normal * load.pressure)

            if force_points:
                actor = pv.PolyData(np.array(force_points))
                actor['vectors'] = np.array(force_vectors)
                actor.set_active_vectors('vectors')
                # --- FIX: The buggy actor.userData line has been removed ---
                self.plotter.add_mesh(actor, name=f"force_actors_{sid}", style='wireframe', render_lines_as_tubes=False, show_edges=False, pickable=False, color='red', opacity=0)

            if moment_points:
                actor = pv.PolyData(np.array(moment_points))
                actor['vectors'] = np.array(moment_vectors)
                actor.set_active_vectors('vectors')
                self.plotter.add_mesh(actor, name=f"moment_actors_{sid}", style='wireframe', render_lines_as_tubes=False, show_edges=False, pickable=False, color='cyan', opacity=0)
            
            if pressure_points:
                actor = pv.PolyData(np.array(pressure_points))
                actor['vectors'] = np.array(pressure_vectors)
                actor.set_active_vectors('vectors')
                self.plotter.add_mesh(actor, name=f"pressure_actors_{sid}", style='wireframe', render_lines_as_tubes=False, show_edges=False, pickable=False, color='yellow', opacity=0)
    
    
    def _create_all_constraint_actors(self):
        """Creates all constraint actors for the current model, one per SID."""
        actor_names_to_remove = [name for name in self.plotter.actors if name.startswith('constraint_actors_')]
        self.plotter.remove_actor(actor_names_to_remove, render=False)

        if not self.current_generator or not self.current_grid: return
        model = self.current_generator.model
        if not model.spcs: return

        for sid, spc_list in model.spcs.items():
            coords = []
            for spc in spc_list:
                for nid in spc.nodes:
                    if nid in model.nodes:
                        coords.append(model.nodes[nid].get_position())
            if coords:
                actor = pv.PolyData(np.array(coords))
                self.plotter.add_mesh(actor, name=f"constraint_actors_{sid}", style='wireframe', render_lines_as_tubes=False, show_edges=False, pickable=False, color='cyan', opacity=0)

    def _update_plot_visibility(self):
        """Fast update of actor visibilities based on tree state. Does not rebuild actors."""
        if not self.current_generator or not self.current_grid:
            self._rebuild_plot()
            return

        visibility_mask = np.ones(self.current_grid.n_cells, dtype=bool)
        model = self.current_generator.model
        type_to_nastran_map = { 'Beams': ['CBEAM'], 'Bars': ['CBAR'], 'Rods': ['CROD'], 'Plates': ['CQUAD4', 'CMEMBRAN', 'CTRIA3'], 'Rigid': ['RBE2', 'RBE3'] }
        iterator = QTreeWidgetItemIterator(self.tree_widget)
        while iterator.value():
            item = iterator.value()
            item_data = item.data(0, QtCore.Qt.UserRole)
            if isinstance(item_data, tuple) and item_data[0] == 'elem_by_type_group':
                if item.checkState(0) != QtCore.Qt.Checked:
                    nastran_types_to_hide = type_to_nastran_map.get(item_data[1], [])
                    if nastran_types_to_hide: visibility_mask[np.isin(self.current_grid.cell_data['type'], nastran_types_to_hide)] = False
            iterator += 1
        props_to_hide = {prop_item.data(0, QtCore.Qt.UserRole)[1] for item in self._find_tree_items(('group', 'properties')) for i in range(item.childCount()) if (prop_item := item.child(i)).checkState(0) != QtCore.Qt.Checked}
        if props_to_hide: visibility_mask[np.isin(self.current_grid.cell_data['PID'], list(props_to_hide))] = False
        mids_to_hide = {mat_item.data(0, QtCore.Qt.UserRole)[1] for item in self._find_tree_items(('group', 'materials')) for i in range(item.childCount()) if (mat_item := item.child(i)).checkState(0) != QtCore.Qt.Checked}
        if mids_to_hide:
            props_using_hidden_mats = {pid for pid, prop in model.properties.items() if (hasattr(prop, 'mid') and prop.mid in mids_to_hide) or (hasattr(prop, 'mids') and not mids_to_hide.isdisjoint(prop.mids))}
            if props_using_hidden_mats: visibility_mask[np.isin(self.current_grid.cell_data['PID'], list(props_using_hidden_mats))] = False
        if (rbe_actor := self.plotter.actors.get('rbe_actors')):
            rbe_actor.SetVisibility(not any(item.checkState(0) != QtCore.Qt.Checked for item in self._find_tree_items(('elem_by_type_group', 'Rigid'))))
        
        self._rebuild_plot(visibility_mask=visibility_mask if np.any(visibility_mask) else None)

        try:
            arrow_size_percent = float(self.arrow_scale_input.text()); use_relative_scaling = self.relative_scaling_check.isChecked()
        except (ValueError, AttributeError):
            arrow_size_percent = 10.0; use_relative_scaling = True
        base_scale = self.current_grid.length * (arrow_size_percent / 100.0)

        max_force_mag = self.load_scaling_info.get('force_max_mag', 0.0)
        max_moment_mag = self.load_scaling_info.get('moment_max_mag', 0.0)
        max_pressure_mag = self.load_scaling_info.get('pressure_max_mag', 0.0)

        for sid in model.loads.keys():
            is_visible = False
            if (sid_item_list := self._find_tree_items(('load_set', sid))) and sid_item_list[0].checkState(0) == QtCore.Qt.Checked:
                is_visible = True
            
            for type, max_mag in [('force', max_force_mag), ('moment', max_moment_mag), ('pressure', max_pressure_mag)]:
                actor_name = f"{type}_actors_{sid}"; glyph_name = f"{type}_glyphs_{sid}"
                self.plotter.remove_actor(glyph_name, render=False)
                if (actor := self.plotter.actors.get(actor_name)) and is_visible:
                    dataset = actor.mapper.dataset
                    if not dataset.n_points: continue

                    # --- FIX: All `add_arrows` calls now use pre-scaled vectors ---
                    if type == 'force':
                        vectors = dataset['vectors']; mags = np.linalg.norm(vectors, axis=1)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                        if use_relative_scaling and max_mag > 1e-9:
                            scaled_vectors = unit_vectors * (mags / max_mag * base_scale)[:, np.newaxis]
                        else:
                            scaled_vectors = unit_vectors * base_scale
                        self.plotter.add_arrows(dataset.points, scaled_vectors, color='red', name=glyph_name)

                    elif type == 'moment':
                        moment_glyph = pv.DoubleArrow(tip_length=0.2, tip_radius=0.08, shaft_radius=0.03)
                        if use_relative_scaling and max_mag > 1e-9:
                            mags = np.linalg.norm(dataset['vectors'], axis=1); dataset['relative_mag'] = mags / max_mag
                            glyphs = dataset.glyph(orient='vectors', scale='relative_mag', factor=base_scale, geom=moment_glyph)
                        else:
                            glyphs = dataset.glyph(orient='vectors', scale=False, factor=base_scale, geom=moment_glyph)
                        self.plotter.add_mesh(glyphs, color='cyan', name=glyph_name)

                    elif type == 'pressure':
                        vectors = dataset['vectors']; mags = np.linalg.norm(vectors, axis=1)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            unit_vectors = np.nan_to_num(vectors / mags[:, np.newaxis])
                        uniform_size = base_scale * 0.75
                        if use_relative_scaling and max_mag > 1e-9:
                            scaled_vectors = unit_vectors * (mags / max_mag * uniform_size)[:, np.newaxis]
                        else:
                            scaled_vectors = unit_vectors * uniform_size
                        self.plotter.add_arrows(dataset.points, scaled_vectors, color='yellow', name=glyph_name)
        
        for sid in model.spcs.keys():
            is_visible = False
            if (sid_item_list := self._find_tree_items(('constraint_set', sid))) and sid_item_list[0].checkState(0) == QtCore.Qt.Checked:
                is_visible = True
            
            glyph_name = f"constraint_glyphs_{sid}"
            self.plotter.remove_actor(glyph_name, render=False)
            if (actor := self.plotter.actors.get(f"constraint_actors_{sid}")) and is_visible:
                dataset = actor.mapper.dataset
                if not dataset.n_points: continue
                glyph_geom = pv.Cone(direction=(0, 0, 1), height=1.0, radius=0.3); scale_factor = self.current_grid.length * 0.02
                glyphs = dataset.glyph(scale=False, factor=scale_factor, orient=False, geom=glyph_geom)
                self.plotter.add_mesh(glyphs, color='cyan', name=glyph_name)
                
        self.plotter.render()
    

if __name__ == '__main__':
    app = QApplication(sys.argv); app.setStyle("Fusion"); app.setPalette(dark_palette)
    window = MainWindow(); window.show()
    sys.exit(app.exec())