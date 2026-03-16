import numpy as np

from PySide6 import QtCore, QtGui
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QMessageBox, QFormLayout, QGroupBox,
    QLineEdit, QComboBox, QLabel, QCheckBox, QDialog,
    QGridLayout, QSplitter,
    QDialogButtonBox, QTableWidget,
    QHeaderView, QListWidget, QTabWidget, QRadioButton,
    QStackedLayout
)
from PySide6.QtGui import QDoubleValidator

from node_runner.dialogs.widgets import OrientationWidget
from node_runner.dialogs.selection import EntitySelectionDialog
from node_runner.utils import get_entity_title_from_comment


def _select_combo_by_data(combo, value):
    """Select a combo box item by its userData value."""
    for i in range(combo.count()):
        if combo.itemData(i) == value:
            combo.setCurrentIndex(i)
            return


class CreateMaterialDialog(QDialog):
    import_requested = QtCore.Signal()

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

        btn_row = QHBoxLayout()
        if not existing_material:
            import_btn = QPushButton("Import from BDF...")
            import_btn.clicked.connect(self.import_requested.emit)
            btn_row.addWidget(import_btn)
        btn_row.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(ok_btn)
        btn_row.addWidget(cancel_btn)
        main_layout.addLayout(btn_row)

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

# START: Corrected replacement for get_parameters method in CreateMaterialDialog
    def get_parameters(self):
        """Retrieves the material parameters from the UI fields."""
        try:
            mat_type_full = self.type_combo.currentText()
            mat_type = mat_type_full.split(' ')[0] # e.g., "MAT1" from "MAT1 (Isotropic)"

            base_params = {
                'mid': int(self.mid_input.text()),
                'comment': f"$ {self.title_input.text()}",
                'type': mat_type
            }

            if mat_type == 'MAT1':
                base_params.update({
                    'E': float(self.mat1_e.text()), 'G': float(self.mat1_g.text()),
                    'nu': float(self.mat1_nu.text()), 'rho': float(self.mat1_rho.text()),
                    'a': float(self.mat1_a.text()), 'tref': float(self.mat1_tref.text()),
                    'ge': float(self.mat1_ge.text())
                })
            elif mat_type == 'MAT8':
                base_params.update({
                    'E1': float(self.mat8_e1.text()), 'E2': float(self.mat8_e2.text()),
                    'nu12': float(self.mat8_nu12.text()), 'G12': float(self.mat8_g12.text()),
                    'G1z': float(self.mat8_g1z.text()), 'G2z': float(self.mat8_g2z.text()),
                    'rho': float(self.mat8_rho.text()), 'a1': float(self.mat8_a1.text()),
                    'a2': float(self.mat8_a2.text()), 'tref': float(self.mat8_tref.text())
                })
            elif mat_type == 'MAT9':
                g_values = {key: float(le.text()) for key, le in self.mat9_gij.items()}
                base_params.update(g_values)
                base_params.update({
                    'rho': float(self.mat9_rho.text()), 'tref': float(self.mat9_tref.text())
                })

            return base_params
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Please provide valid numerical inputs for all material fields. Error: {e}")
            return None
# END: Corrected replacement for get_parameters method in CreateMaterialDialog

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

        self.type_combo.addItems(["PSHELL", "PCOMP", "PBAR", "PBEAM", "PBUSH", "PSOLID", "PROD", "PSHEAR", "PGAP"])

        form_layout.addRow("Property ID (PID):", self.pid_input)
        form_layout.addRow("Title:", self.title_input)
        form_layout.addRow("Property Type:", self.type_combo)
        main_layout.addLayout(form_layout)

        self.stacked_layout = QStackedLayout()
        self.pshell_widget = self._create_pshell_ui()
        self.pcomp_widget = self._create_pcomp_ui()
        self.pbar_widget = self._create_pbar_ui()

        self.pbeam_widget = self._create_pbeam_ui()

        self.pbush_widget = self._create_pbush_ui()
        self.psolid_widget = self._create_psolid_ui()
        self.prod_widget = self._create_prod_ui()
        self.pshear_widget = self._create_pshear_ui()
        self.pgap_widget = self._create_pgap_ui()

        self.stacked_layout.addWidget(self.pshell_widget)
        self.stacked_layout.addWidget(self.pcomp_widget)
        self.stacked_layout.addWidget(self.pbar_widget)
        self.stacked_layout.addWidget(self.pbeam_widget)
        self.stacked_layout.addWidget(self.pbush_widget)
        self.stacked_layout.addWidget(self.psolid_widget)
        self.stacked_layout.addWidget(self.prod_widget)
        self.stacked_layout.addWidget(self.pshear_widget)
        self.stacked_layout.addWidget(self.pgap_widget)
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
        for mid, mat in sorted(self.model.materials.items()):
            title = get_entity_title_from_comment(mat.comment, "Material", mid)
            combo.addItem(f"{mid}: {title}", mid)
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
        pbar_lib_btn = QPushButton("Section Library...")
        pbar_lib_btn.clicked.connect(lambda: self._open_section_library('pbar'))
        layout.addRow(pbar_lib_btn)
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
        pbeam_lib_btn = QPushButton("Section Library...")
        pbeam_lib_btn.clicked.connect(lambda: self._open_section_library('pbeam'))
        layout.addRow(pbeam_lib_btn)
        return widget

    def _open_section_library(self, target):
        dlg = BeamSectionLibraryDialog(self)
        if dlg.exec() == QDialog.Accepted:
            props = dlg.get_properties()
            if target == 'pbar':
                self.pbar_A.setText(f"{props['A']:.6g}")
                self.pbar_i1.setText(f"{props['I1']:.6g}")
                self.pbar_i2.setText(f"{props['I2']:.6g}")
                self.pbar_j.setText(f"{props['J']:.6g}")
            else:
                self.pbeam_A.setText(f"{props['A']:.6g}")
                self.pbeam_i1.setText(f"{props['I1']:.6g}")
                self.pbeam_i2.setText(f"{props['I2']:.6g}")
                self.pbeam_j.setText(f"{props['J']:.6g}")

    def _populate_from_existing(self, prop):
        self.title_input.setText(prop.comment.strip().lstrip('$').strip())
        prop_type_map = {"PSHELL": 0, "PCOMP": 1, "PBAR": 2, "PBEAM": 3, "PBUSH": 4,
                         "PSOLID": 5, "PROD": 6, "PSHEAR": 7, "PGAP": 8}
        self.type_combo.setCurrentIndex(prop_type_map.get(prop.type, 0))
        self.type_combo.setEnabled(False)

        if prop.type == 'PSHELL':
            _select_combo_by_data(self.pshell_mid1, prop.mid1)
            self.pshell_t.setText(str(prop.t))
            self.pshell_nsm.setText(str(prop.nsm))
        elif prop.type == 'PCOMP':
            self.pcomp_nsm.setText(str(prop.nsm))
            self.pcomp_ft.setCurrentText(prop.ft)
            for i in range(prop.nplies):
                self._add_pcomp_ply()
                _select_combo_by_data(self.pcomp_ply_table.cellWidget(i, 0), prop.mids[i])
                self.pcomp_ply_table.cellWidget(i, 1).setText(str(prop.thicknesses[i]))
                self.pcomp_ply_table.cellWidget(i, 2).setText(str(prop.thetas[i]))
                self.pcomp_ply_table.cellWidget(i, 3).setCurrentText("YES" if prop.souts[i] == "YES" else "NO")
        elif prop.type == 'PBAR':
            _select_combo_by_data(self.pbar_mid, prop.mid)
            self.pbar_A.setText(str(prop.A)); self.pbar_i1.setText(str(prop.i1));
            self.pbar_i2.setText(str(prop.i2)); self.pbar_j.setText(str(prop.j));
        elif prop.type == 'PBEAM':
            _select_combo_by_data(self.pbeam_mid, prop.mid)
            self.pbeam_A.setText(str(prop.A[0])); self.pbeam_i1.setText(str(prop.i1[0]));
            self.pbeam_i2.setText(str(prop.i2[0])); self.pbeam_j.setText(str(prop.j[0]));
        elif prop.type == 'PBUSH':
            k_values = prop.k
            self.pbush_k1.setText(str(k_values[0] if len(k_values) > 0 else 0.0))
            self.pbush_k2.setText(str(k_values[1] if len(k_values) > 1 else 0.0))
            self.pbush_k3.setText(str(k_values[2] if len(k_values) > 2 else 0.0))
            self.pbush_k4.setText(str(k_values[3] if len(k_values) > 3 else 0.0))
            self.pbush_k5.setText(str(k_values[4] if len(k_values) > 4 else 0.0))
            self.pbush_k6.setText(str(k_values[5] if len(k_values) > 5 else 0.0))
        elif prop.type == 'PSOLID':
            _select_combo_by_data(self.psolid_mid, prop.mid)
        elif prop.type == 'PROD':
            _select_combo_by_data(self.prod_mid, prop.mid)
            self.prod_A.setText(str(prop.A))
            self.prod_j.setText(str(prop.j))
            self.prod_c.setText(str(prop.c))
            self.prod_nsm.setText(str(prop.nsm))
        elif prop.type == 'PSHEAR':
            _select_combo_by_data(self.pshear_mid, prop.mid)
            self.pshear_t.setText(str(prop.t))
            self.pshear_nsm.setText(str(prop.nsm))
        elif prop.type == 'PGAP':
            self.pgap_u0.setText(str(prop.u0))
            self.pgap_f0.setText(str(prop.f0))
            self.pgap_ka.setText(str(prop.ka))
            if prop.kb is not None:
                self.pgap_kb.setText(str(prop.kb))

    # START: New UI method in CreatePropertyDialog class in main.py
    def _create_pbush_ui(self):
        """Creates the UI panel for PBUSH stiffness inputs."""
        widget = QWidget()
        layout = QGridLayout(widget)

        self.pbush_k1 = self._create_prop_qlineedit("0.0")
        self.pbush_k2 = self._create_prop_qlineedit("0.0")
        self.pbush_k3 = self._create_prop_qlineedit("0.0")
        self.pbush_k4 = self._create_prop_qlineedit("0.0")
        self.pbush_k5 = self._create_prop_qlineedit("0.0")
        self.pbush_k6 = self._create_prop_qlineedit("0.0")

        layout.addWidget(QLabel("K1 (TX Stiffness):"), 0, 0)
        layout.addWidget(self.pbush_k1, 0, 1)
        layout.addWidget(QLabel("K2 (TY Stiffness):"), 1, 0)
        layout.addWidget(self.pbush_k2, 1, 1)
        layout.addWidget(QLabel("K3 (TZ Stiffness):"), 2, 0)
        layout.addWidget(self.pbush_k3, 2, 1)

        layout.addWidget(QLabel("K4 (RX Stiffness):"), 0, 2)
        layout.addWidget(self.pbush_k4, 0, 3)
        layout.addWidget(QLabel("K5 (RY Stiffness):"), 1, 2)
        layout.addWidget(self.pbush_k5, 1, 3)
        layout.addWidget(QLabel("K6 (RZ Stiffness):"), 2, 2)
        layout.addWidget(self.pbush_k6, 2, 3)

        return widget
# END: New UI method in CreatePropertyDialog class in main.py

    def _create_psolid_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.psolid_mid = self._get_material_combo()
        layout.addRow("Material ID (MID):", self.psolid_mid)
        return widget

    def _create_prod_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.prod_mid = self._get_material_combo()
        self.prod_A = self._create_prop_qlineedit("1.0")
        self.prod_j = self._create_prop_qlineedit("0.0")
        self.prod_c = self._create_prop_qlineedit("0.0")
        self.prod_nsm = self._create_prop_qlineedit("0.0")
        layout.addRow("Material ID (MID):", self.prod_mid)
        layout.addRow("Area (A):", self.prod_A)
        layout.addRow("Torsional Constant (J):", self.prod_j)
        layout.addRow("Stress Recovery Coeff (c):", self.prod_c)
        layout.addRow("Non-Structural Mass (nsm):", self.prod_nsm)
        return widget

    def _create_pshear_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.pshear_mid = self._get_material_combo()
        self.pshear_t = self._create_prop_qlineedit("0.1")
        self.pshear_nsm = self._create_prop_qlineedit("0.0")
        layout.addRow("Material ID (MID):", self.pshear_mid)
        layout.addRow("Thickness (t):", self.pshear_t)
        layout.addRow("Non-Structural Mass (nsm):", self.pshear_nsm)
        return widget

    def _create_pgap_ui(self):
        widget = QWidget()
        layout = QFormLayout(widget)
        self.pgap_u0 = self._create_prop_qlineedit("0.0")
        self.pgap_f0 = self._create_prop_qlineedit("0.0")
        self.pgap_ka = self._create_prop_qlineedit("1.0e8")
        self.pgap_kb = self._create_prop_qlineedit("")
        self.pgap_kb.setPlaceholderText("Optional (blank = default)")
        layout.addRow("Initial Gap (u0):", self.pgap_u0)
        layout.addRow("Preload (f0):", self.pgap_f0)
        layout.addRow("Axial Stiffness (ka):", self.pgap_ka)
        layout.addRow("Contact Stiffness (kb):", self.pgap_kb)
        return widget

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
                    'mid1': int(self.pshell_mid1.currentData()), 't': float(self.pshell_t.text()), 'nsm': float(self.pshell_nsm.text())
                })
            elif prop_type == 'PCOMP':
                plies = []
                for row in range(self.pcomp_ply_table.rowCount()):
                    mid = int(self.pcomp_ply_table.cellWidget(row, 0).currentData())
                    t = float(self.pcomp_ply_table.cellWidget(row, 1).text())
                    theta = float(self.pcomp_ply_table.cellWidget(row, 2).text())
                    sout = self.pcomp_ply_table.cellWidget(row, 3).currentText()
                    plies.append([mid, t, theta, sout])
                base_params.update({
                    'nsm': float(self.pcomp_nsm.text()), 'ft': self.pcomp_ft.currentText(), 'plies': plies
                })
            elif prop_type == 'PBAR':
                base_params.update({
                    'mid': int(self.pbar_mid.currentData()), 'A': float(self.pbar_A.text()),
                    'i1': float(self.pbar_i1.text()), 'i2': float(self.pbar_i2.text()), 'j': float(self.pbar_j.text())
                })
            elif prop_type == 'PBEAM':
                 base_params.update({
                    'mid': int(self.pbeam_mid.currentData()), 'A': float(self.pbeam_A.text()),
                    'i1': float(self.pbeam_i1.text()), 'i2': float(self.pbeam_i2.text()), 'j': float(self.pbeam_j.text()),
                })
            elif prop_type == 'PBUSH':
                k_values = [
                    float(self.pbush_k1.text()), float(self.pbush_k2.text()),
                    float(self.pbush_k3.text()), float(self.pbush_k4.text()),
                    float(self.pbush_k5.text()), float(self.pbush_k6.text())
                ]
                base_params.update({'k': k_values})
            elif prop_type == 'PSOLID':
                base_params.update({'mid': int(self.psolid_mid.currentData())})
            elif prop_type == 'PROD':
                base_params.update({
                    'mid': int(self.prod_mid.currentData()), 'A': float(self.prod_A.text()),
                    'j': float(self.prod_j.text()), 'c': float(self.prod_c.text()),
                    'nsm': float(self.prod_nsm.text())
                })
            elif prop_type == 'PSHEAR':
                base_params.update({
                    'mid': int(self.pshear_mid.currentData()), 't': float(self.pshear_t.text()),
                    'nsm': float(self.pshear_nsm.text())
                })
            elif prop_type == 'PGAP':
                kb_text = self.pgap_kb.text().strip()
                base_params.update({
                    'u0': float(self.pgap_u0.text()), 'f0': float(self.pgap_f0.text()),
                    'ka': float(self.pgap_ka.text()),
                    'kb': float(kb_text) if kb_text else None
                })
            return base_params
        except (ValueError, KeyError, AttributeError) as e:
            QMessageBox.warning(self, "Input Error", f"Invalid input for property data. Error: {e}")
            return None

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

class CreateBushElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent
        self.model = model
        self.setWindowTitle("Create Bush Element")
        main_layout = QVBoxLayout(self)

        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid))
        self.pid_combo = QComboBox()
        pbush_ids = [str(pid) for pid, prop in model.properties.items() if prop.type == 'PBUSH']
        self.pid_combo.addItems(pbush_ids)

        self.n1_input, self.n2_input = QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn = QPushButton("Pick..."), QPushButton("Pick...")
        self.ground_n2_check = QCheckBox("Ground Node 2")

        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn); n2_layout.addWidget(self.ground_n2_check)

        top_form.addRow("Element ID (0 for auto):", self.eid_input)
        top_form.addRow("Property ID (PID):", self.pid_combo)
        top_form.addRow("Node 1 ID:", n1_widget)
        top_form.addRow("Node 2 ID:", n2_widget)
        main_layout.addLayout(top_form)

        self.orientation_widget = OrientationWidget(model, self)
        main_layout.addWidget(self.orientation_widget)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        self.ground_n2_check.toggled.connect(self.n2_input.setDisabled)
        self.ground_n2_check.toggled.connect(pick_n2_btn.setDisabled)
        pick_n1_btn.clicked.connect(self._pick_node_1)
        pick_n2_btn.clicked.connect(self._pick_node_2)
        # --- CONNECT THE WIDGET'S SIGNAL TO OUR PICKER METHOD ---
        self.orientation_widget.pick_orientation_node_requested.connect(self._pick_orientation_node)

    def _pick_node_1(self):
        if self.parent_window: self.hide(); self.parent_window._activate_single_node_picker(self.n1_input.setText, calling_dialog=self)
    def _pick_node_2(self):
        if self.parent_window: self.hide(); self.parent_window._activate_single_node_picker(self.n2_input.setText, calling_dialog=self)

    # --- ADD THIS METHOD BACK ---
    def _pick_orientation_node(self):
        if self.parent_window: self.hide(); self.parent_window._activate_single_node_picker(self.orientation_widget.orient_node_id.setText, calling_dialog=self)

    def get_parameters(self):
        try:
            params = {
                'eid': int(self.eid_input.text()),
                'pid': int(self.pid_combo.currentText()),
                'n1': int(self.n1_input.text()),
                'n2': None if self.ground_n2_check.isChecked() else int(self.n2_input.text()),
                'orientation': self.orientation_widget.get_orientation()
            }
            return params
        except (ValueError, KeyError, AttributeError):
            QMessageBox.warning(self, "Input Error", "Invalid input. Ensure all required fields are filled correctly.")
            return None
# END: Updated CreateBushElementDialog in main.py


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
        """Open EntitySelectionDialog non-blocking for node picking."""
        self._pending_target_list = target_list
        dialog = EntitySelectionDialog('Node', self.all_node_ids, self)
        dialog.request_show_selection.connect(self.parent_window._highlight_entities)
        dialog.request_picking_mode.connect(self._on_picking_mode_requested)
        dialog.accepted.connect(lambda: self._on_selection_accepted(dialog))
        dialog.rejected.connect(lambda: self._on_selection_rejected(dialog))
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
        # Store reference so mainwindow can manage it
        self.parent_window.active_selection_dialog = dialog
        # Hide this RBE dialog and show selection dialog
        self.hide()
        dialog.show()

    def _on_picking_mode_requested(self, entity_type, enabled):
        """When EntitySelectionDialog requests picking, ensure RBE dialog is hidden."""
        self.parent_window._set_picking_mode(entity_type, enabled)

    def _on_selection_accepted(self, dialog):
        """Process accepted selection and restore RBE dialog."""
        self._pending_target_list.clear()
        self._pending_target_list.addItems([str(nid) for nid in dialog.get_selected_ids()])
        self._cleanup_selection(dialog)

    def _on_selection_rejected(self, dialog):
        """Handle cancelled selection and restore RBE dialog."""
        self._cleanup_selection(dialog)

    def _cleanup_selection(self, dialog):
        """Disconnect signals, clear highlights, restore RBE dialog."""
        try:
            dialog.request_show_selection.disconnect(self.parent_window._highlight_entities)
            dialog.request_picking_mode.disconnect(self._on_picking_mode_requested)
        except RuntimeError:
            pass
        self.parent_window._highlight_entities('Node', [])
        self.parent_window.active_selection_dialog = None
        self.show()
        self.activateWindow()
        dialog.deleteLater()
    def get_parameters(self):
        try:
            indep_nodes = [int(self.indep_list.item(i).text()) for i in range(self.indep_list.count())]; dep_nodes = [int(self.dep_list.item(i).text()) for i in range(self.dep_list.count())]
            if not indep_nodes or not dep_nodes: QMessageBox.warning(self, "Input Error", "Must select at least one independent and one dependent node."); return None
            return {'eid': int(self.eid_input.text()), 'elem_type': self.type_combo.currentText(), 'dof': self.dof_input.text(), 'indep_nodes': indep_nodes, 'dep_nodes': dep_nodes}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None

class CreateConm2Dialog(QDialog):
    """Dialog for creating a CONM2 concentrated mass element."""
    def __init__(self, next_eid, all_node_ids, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.all_node_ids = all_node_ids
        self.setWindowTitle("Create Concentrated Mass (CONM2)"); self.setMinimumWidth(400)
        main_layout = QVBoxLayout(self)

        form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid))
        form.addRow("Element ID (0 for auto):", self.eid_input)

        # Node selection row
        node_row = QHBoxLayout()
        self.node_input = QLineEdit(); self.node_input.setPlaceholderText("Click Select...")
        self.node_input.setReadOnly(True)
        select_node_btn = QPushButton("Select...")
        select_node_btn.clicked.connect(self._select_node)
        node_row.addWidget(self.node_input); node_row.addWidget(select_node_btn)
        form.addRow("Node:", node_row)

        self.mass_input = QLineEdit(); self.mass_input.setPlaceholderText("e.g. 100.0")
        form.addRow("Mass:", self.mass_input)

        self.cid_input = QLineEdit("0")
        form.addRow("Coord System (CID):", self.cid_input)

        main_layout.addLayout(form)

        # Offset group
        offset_group = QGroupBox("Offset (optional)")
        offset_layout = QFormLayout(offset_group)
        self.x1_input = QLineEdit("0.0"); self.x2_input = QLineEdit("0.0"); self.x3_input = QLineEdit("0.0")
        offset_layout.addRow("X1:", self.x1_input); offset_layout.addRow("X2:", self.x2_input); offset_layout.addRow("X3:", self.x3_input)
        main_layout.addWidget(offset_group)

        # Inertia group
        inertia_group = QGroupBox("Inertia (optional)")
        inertia_layout = QFormLayout(inertia_group)
        self.i11_input = QLineEdit("0.0"); self.i21_input = QLineEdit("0.0"); self.i22_input = QLineEdit("0.0")
        self.i31_input = QLineEdit("0.0"); self.i32_input = QLineEdit("0.0"); self.i33_input = QLineEdit("0.0")
        inertia_layout.addRow("I11:", self.i11_input); inertia_layout.addRow("I21:", self.i21_input); inertia_layout.addRow("I22:", self.i22_input)
        inertia_layout.addRow("I31:", self.i31_input); inertia_layout.addRow("I32:", self.i32_input); inertia_layout.addRow("I33:", self.i33_input)
        inertia_group.setCheckable(True); inertia_group.setChecked(False)
        main_layout.addWidget(inertia_group)
        self.inertia_group = inertia_group

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        main_layout.addWidget(button_box)

    def _select_node(self):
        """Open EntitySelectionDialog non-blocking for node picking."""
        dialog = EntitySelectionDialog('Node', self.all_node_ids, self)
        dialog.request_show_selection.connect(self.parent_window._highlight_entities)
        dialog.request_picking_mode.connect(self._on_picking_mode_requested)
        dialog.accepted.connect(lambda: self._on_selection_accepted(dialog))
        dialog.rejected.connect(lambda: self._on_selection_rejected(dialog))
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
        self.parent_window.active_selection_dialog = dialog
        self.hide()
        dialog.show()

    def _on_picking_mode_requested(self, entity_type, enabled):
        self.parent_window._set_picking_mode(entity_type, enabled)

    def _on_selection_accepted(self, dialog):
        selected = dialog.get_selected_ids()
        if selected:
            self.node_input.setText(str(selected[0]))
        self._cleanup_selection(dialog)

    def _on_selection_rejected(self, dialog):
        self._cleanup_selection(dialog)

    def _cleanup_selection(self, dialog):
        try:
            dialog.request_show_selection.disconnect(self.parent_window._highlight_entities)
            dialog.request_picking_mode.disconnect(self._on_picking_mode_requested)
        except RuntimeError:
            pass
        self.parent_window._highlight_entities('Node', [])
        self.parent_window.active_selection_dialog = None
        self.show()
        self.activateWindow()
        dialog.deleteLater()

    def get_parameters(self):
        try:
            nid = int(self.node_input.text())
            mass = float(self.mass_input.text())
            if mass <= 0:
                QMessageBox.warning(self, "Input Error", "Mass must be greater than 0."); return None
            params = {
                'eid': int(self.eid_input.text()),
                'nid': nid,
                'mass': mass,
                'cid': int(self.cid_input.text()),
                'X': [float(self.x1_input.text()), float(self.x2_input.text()), float(self.x3_input.text())],
            }
            if self.inertia_group.isChecked():
                params['I'] = [float(self.i11_input.text()), float(self.i21_input.text()), float(self.i22_input.text()),
                               float(self.i31_input.text()), float(self.i32_input.text()), float(self.i33_input.text())]
            return params
        except (ValueError, KeyError):
            QMessageBox.warning(self, "Input Error", "Invalid input. Please check all fields."); return None


class CreateSolidElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Solid Element"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid)); self.type_combo = QComboBox(); self.type_combo.addItems(["CHEXA", "CTETRA", "CPENTA"])
        self.pid_combo = QComboBox(); self.pid_combo.addItems([str(pid) for pid in model.properties.keys()])
        self.node_inputs = []
        self.node_labels = []
        self.node_widgets = []
        for i in range(8):
            ni_input = QLineEdit(); pick_btn = QPushButton("Pick...")
            widget = QWidget(); layout = QHBoxLayout(widget); layout.setContentsMargins(0,0,0,0); layout.addWidget(ni_input); layout.addWidget(pick_btn)
            label = QLabel(f"Node {i+1} ID:")
            self.node_inputs.append(ni_input); self.node_labels.append(label); self.node_widgets.append(widget)
            pick_btn.clicked.connect(lambda checked, inp=ni_input: self._pick_node(inp))
            ni_input.textChanged.connect(self._update_highlight)
        top_form.addRow("Element ID (0 for auto):", self.eid_input); top_form.addRow("Element Type:", self.type_combo); top_form.addRow("Property ID (PID):", self.pid_combo)
        for i in range(8): top_form.addRow(self.node_labels[i], self.node_widgets[i])
        main_layout.addLayout(top_form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        self.type_combo.currentTextChanged.connect(self._toggle_nodes)
        self._toggle_nodes(self.type_combo.currentText())
    def _pick_node(self, input_widget):
        if self.parent_window: self.parent_window._activate_single_node_picker(input_widget.setText)
    def _toggle_nodes(self, text):
        node_count = {'CHEXA': 8, 'CTETRA': 4, 'CPENTA': 6}.get(text, 8)
        for i in range(8):
            visible = (i < node_count)
            self.node_widgets[i].setVisible(visible); self.node_labels[i].setVisible(visible)
        self._update_highlight()
    def _update_highlight(self):
        node_count = {'CHEXA': 8, 'CTETRA': 4, 'CPENTA': 6}.get(self.type_combo.currentText(), 8)
        nodes_to_highlight = []
        for i in range(node_count):
            try: nodes_to_highlight.append(int(self.node_inputs[i].text()))
            except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def get_parameters(self):
        try:
            node_count = {'CHEXA': 8, 'CTETRA': 4, 'CPENTA': 6}.get(self.type_combo.currentText(), 8)
            nodes = [int(self.node_inputs[i].text()) for i in range(node_count)]
            return {'eid': int(self.eid_input.text()), 'pid': int(self.pid_combo.currentText()), 'nodes': nodes, 'type': self.type_combo.currentText()}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None


class CreateShearElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Shear Panel (CSHEAR)"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid))
        self.pid_combo = QComboBox(); self.pid_combo.addItems([str(pid) for pid in model.properties.keys()])
        self.n1_input, self.n2_input, self.n3_input, self.n4_input = QLineEdit(), QLineEdit(), QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn, pick_n3_btn, pick_n4_btn = QPushButton("Pick..."), QPushButton("Pick..."), QPushButton("Pick..."), QPushButton("Pick...")
        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn)
        n3_widget = QWidget(); n3_layout = QHBoxLayout(n3_widget); n3_layout.setContentsMargins(0,0,0,0); n3_layout.addWidget(self.n3_input); n3_layout.addWidget(pick_n3_btn)
        n4_widget = QWidget(); n4_layout = QHBoxLayout(n4_widget); n4_layout.setContentsMargins(0,0,0,0); n4_layout.addWidget(self.n4_input); n4_layout.addWidget(pick_n4_btn)
        top_form.addRow("Element ID (0 for auto):", self.eid_input); top_form.addRow("Property ID (PID):", self.pid_combo)
        top_form.addRow("Node 1 ID:", n1_widget); top_form.addRow("Node 2 ID:", n2_widget); top_form.addRow("Node 3 ID:", n3_widget); top_form.addRow("Node 4 ID:", n4_widget)
        main_layout.addLayout(top_form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        pick_n1_btn.clicked.connect(lambda: self._pick_node(self.n1_input)); pick_n2_btn.clicked.connect(lambda: self._pick_node(self.n2_input))
        pick_n3_btn.clicked.connect(lambda: self._pick_node(self.n3_input)); pick_n4_btn.clicked.connect(lambda: self._pick_node(self.n4_input))
        self.n1_input.textChanged.connect(self._update_highlight); self.n2_input.textChanged.connect(self._update_highlight)
        self.n3_input.textChanged.connect(self._update_highlight); self.n4_input.textChanged.connect(self._update_highlight)
    def _pick_node(self, input_widget):
        if self.parent_window: self.parent_window._activate_single_node_picker(input_widget.setText)
    def _update_highlight(self):
        nodes_to_highlight = []
        for inp in [self.n1_input, self.n2_input, self.n3_input, self.n4_input]:
            try: nodes_to_highlight.append(int(inp.text()))
            except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def get_parameters(self):
        try:
            nodes = [int(self.n1_input.text()), int(self.n2_input.text()), int(self.n3_input.text()), int(self.n4_input.text())]
            return {'eid': int(self.eid_input.text()), 'pid': int(self.pid_combo.currentText()), 'nodes': nodes}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None


class CreateGapElementDialog(QDialog):
    def __init__(self, next_eid, model, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Gap Element (CGAP)"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid))
        self.pid_combo = QComboBox(); self.pid_combo.addItems([str(pid) for pid in model.properties.keys()])
        self.n1_input, self.n2_input = QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn = QPushButton("Pick..."), QPushButton("Pick...")
        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn)
        top_form.addRow("Element ID (0 for auto):", self.eid_input); top_form.addRow("Property ID (PID):", self.pid_combo)
        top_form.addRow("Node A ID:", n1_widget); top_form.addRow("Node B ID:", n2_widget)
        # Orientation
        orient_group = QGroupBox("Orientation Vector")
        orient_layout = QFormLayout(orient_group)
        self.orient_method = QComboBox(); self.orient_method.addItems(["Vector", "Node"])
        self.orient_x, self.orient_y, self.orient_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("1.0")
        for w in [self.orient_x, self.orient_y, self.orient_z]: w.setValidator(QDoubleValidator())
        orient_vec_widget = QWidget(); orient_vec_layout = QHBoxLayout(orient_vec_widget); orient_vec_layout.setContentsMargins(0,0,0,0)
        orient_vec_layout.addWidget(self.orient_x); orient_vec_layout.addWidget(self.orient_y); orient_vec_layout.addWidget(self.orient_z)
        self.orient_node_input = QLineEdit(); self.orient_node_label = QLabel("Orientation Node ID:")
        orient_layout.addRow("Method:", self.orient_method)
        orient_layout.addRow("X, Y, Z:", orient_vec_widget)
        self.orient_vec_widget = orient_vec_widget; self.orient_vec_label = orient_layout.labelForField(orient_vec_widget)
        orient_layout.addRow(self.orient_node_label, self.orient_node_input)
        top_form.addRow(orient_group)
        main_layout.addLayout(top_form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        pick_n1_btn.clicked.connect(lambda: self._pick_node(self.n1_input)); pick_n2_btn.clicked.connect(lambda: self._pick_node(self.n2_input))
        self.n1_input.textChanged.connect(self._update_highlight); self.n2_input.textChanged.connect(self._update_highlight)
        self.orient_method.currentTextChanged.connect(self._toggle_orient)
        self._toggle_orient(self.orient_method.currentText())
    def _pick_node(self, input_widget):
        if self.parent_window: self.parent_window._activate_single_node_picker(input_widget.setText)
    def _toggle_orient(self, text):
        is_vector = (text == "Vector")
        self.orient_vec_widget.setVisible(is_vector)
        if self.orient_vec_label: self.orient_vec_label.setVisible(is_vector)
        self.orient_node_input.setVisible(not is_vector); self.orient_node_label.setVisible(not is_vector)
    def _update_highlight(self):
        nodes_to_highlight = []
        try: nodes_to_highlight.append(int(self.n1_input.text()))
        except ValueError: pass
        try: nodes_to_highlight.append(int(self.n2_input.text()))
        except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def get_parameters(self):
        try:
            orientation = None
            if self.orient_method.currentText() == "Vector":
                orientation = {'method': 'vector', 'values': [float(self.orient_x.text()), float(self.orient_y.text()), float(self.orient_z.text())]}
            else:
                orientation = {'method': 'node', 'values': [int(self.orient_node_input.text())]}
            return {'eid': int(self.eid_input.text()), 'pid': int(self.pid_combo.currentText()), 'n1': int(self.n1_input.text()), 'n2': int(self.n2_input.text()), 'orientation': orientation}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None


class CreatePlotelDialog(QDialog):
    def __init__(self, next_eid, all_node_ids, parent=None):
        super().__init__(parent)
        self.parent_window = parent; self.setWindowTitle("Create Plot Element (PLOTEL)"); main_layout = QVBoxLayout(self)
        top_form = QFormLayout()
        self.eid_input = QLineEdit(str(next_eid))
        self.n1_input, self.n2_input = QLineEdit(), QLineEdit()
        pick_n1_btn, pick_n2_btn = QPushButton("Pick..."), QPushButton("Pick...")
        n1_widget = QWidget(); n1_layout = QHBoxLayout(n1_widget); n1_layout.setContentsMargins(0,0,0,0); n1_layout.addWidget(self.n1_input); n1_layout.addWidget(pick_n1_btn)
        n2_widget = QWidget(); n2_layout = QHBoxLayout(n2_widget); n2_layout.setContentsMargins(0,0,0,0); n2_layout.addWidget(self.n2_input); n2_layout.addWidget(pick_n2_btn)
        top_form.addRow("Element ID (0 for auto):", self.eid_input)
        top_form.addRow("Node 1 ID:", n1_widget); top_form.addRow("Node 2 ID:", n2_widget)
        main_layout.addLayout(top_form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)
        button_box.accepted.connect(self.accept); button_box.rejected.connect(self.reject)
        pick_n1_btn.clicked.connect(lambda: self._pick_node(self.n1_input)); pick_n2_btn.clicked.connect(lambda: self._pick_node(self.n2_input))
        self.n1_input.textChanged.connect(self._update_highlight); self.n2_input.textChanged.connect(self._update_highlight)
    def _pick_node(self, input_widget):
        if self.parent_window: self.parent_window._activate_single_node_picker(input_widget.setText)
    def _update_highlight(self):
        nodes_to_highlight = []
        try: nodes_to_highlight.append(int(self.n1_input.text()))
        except ValueError: pass
        try: nodes_to_highlight.append(int(self.n2_input.text()))
        except ValueError: pass
        if self.parent_window: self.parent_window._highlight_entities('Node', nodes_to_highlight)
    def get_parameters(self):
        try:
            return {'eid': int(self.eid_input.text()), 'nodes': [int(self.n1_input.text()), int(self.n2_input.text())]}
        except (ValueError, KeyError): QMessageBox.warning(self, "Input Error", "Invalid input."); return None


class CreateCoordDialog(QDialog):
    def __init__(self, next_cid, model, parent=None, existing_coord=None):
        super().__init__(parent)
        self.model = model
        title = "Edit Coordinate System" if existing_coord else "Create Coordinate System"
        self.setWindowTitle(title)
        self.setMinimumWidth(450)

        # --- Main Layout ---
        main_layout = QVBoxLayout(self)
        top_form = QFormLayout()

        # --- Top-Level Controls (always visible) ---
        self.cid_input = QLineEdit(str(next_cid))
        self.cid_input.setValidator(QtGui.QIntValidator(3, 99999999, self))

        self.type_combo = QComboBox()
        self.type_combo.addItems(["Rectangular", "Cylindrical", "Spherical"])

        self.method_combo = QComboBox()
        self.method_combo.addItems(["3 Points", "Rotate and Translate"])

        self.title_input = QLineEdit()

        top_form.addRow("Coordinate ID (CID):", self.cid_input)
        top_form.addRow("Title:", self.title_input)
        top_form.addRow("System Type:", self.type_combo)
        top_form.addRow("Definition Method:", self.method_combo)
        main_layout.addLayout(top_form)

        # --- Stacked Layout for Different Method Panels ---
        self.stacked_layout = QStackedLayout()
        self.panel_3_points = self._create_3_points_panel()
        self.panel_rotate_translate = self._create_rotate_translate_panel()
        self.stacked_layout.addWidget(self.panel_3_points)
        self.stacked_layout.addWidget(self.panel_rotate_translate)
        main_layout.addLayout(self.stacked_layout)

        # --- Bottom Button Box ---
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        # --- Connections ---
        self.method_combo.currentIndexChanged.connect(self.stacked_layout.setCurrentIndex)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        # --- Populate for Edit Mode ---
        if existing_coord:
            self.cid_input.setText(str(existing_coord.cid))
            self.cid_input.setEnabled(False)
            type_map = {'CORD2R': 0, 'CORD2C': 1, 'CORD2S': 2}
            self.type_combo.setCurrentIndex(type_map.get(existing_coord.type, 0))

            self.method_combo.setCurrentIndex(0)
            self.method_combo.setEnabled(False)

            # Populate title from comment
            if hasattr(existing_coord, 'comment') and existing_coord.comment:
                self.title_input.setText(existing_coord.comment.strip().lstrip('$').strip())

            # --- FIX: Changed existing_coord.Origin() to existing_coord.origin ---
            origin = existing_coord.origin
            self.p3_origin_x.setText(str(origin[0])); self.p3_origin_y.setText(str(origin[1])); self.p3_origin_z.setText(str(origin[2]))

            z_axis_pt = existing_coord.e2; self.p3_z_axis_x.setText(str(z_axis_pt[0])); self.p3_z_axis_y.setText(str(z_axis_pt[1])); self.p3_z_axis_z.setText(str(z_axis_pt[2]))
            xz_plane_pt = existing_coord.e3; self.p3_xz_plane_x.setText(str(xz_plane_pt[0])); self.p3_xz_plane_y.setText(str(xz_plane_pt[1])); self.p3_xz_plane_z.setText(str(xz_plane_pt[2]))

    def _create_3_points_panel(self):
        """Creates the widget with inputs for the '3 Points' definition method."""
        panel = QWidget()
        form = QFormLayout(panel)
        form.setContentsMargins(0, 10, 0, 0)

        self.p3_origin_x, self.p3_origin_y, self.p3_origin_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        self.p3_z_axis_x, self.p3_z_axis_y, self.p3_z_axis_z = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("1.0")
        self.p3_xz_plane_x, self.p3_xz_plane_y, self.p3_xz_plane_z = QLineEdit("1.0"), QLineEdit("0.0"), QLineEdit("0.0")

        for w in [self.p3_origin_x, self.p3_origin_y, self.p3_origin_z, self.p3_z_axis_x, self.p3_z_axis_y, self.p3_z_axis_z, self.p3_xz_plane_x, self.p3_xz_plane_y, self.p3_xz_plane_z]:
            w.setValidator(QDoubleValidator())

        p3_origin_layout = QHBoxLayout(); p3_origin_layout.addWidget(self.p3_origin_x); p3_origin_layout.addWidget(self.p3_origin_y); p3_origin_layout.addWidget(self.p3_origin_z)
        p3_z_axis_layout = QHBoxLayout(); p3_z_axis_layout.addWidget(self.p3_z_axis_x); p3_z_axis_layout.addWidget(self.p3_z_axis_y); p3_z_axis_layout.addWidget(self.p3_z_axis_z)
        p3_xz_plane_layout = QHBoxLayout(); p3_xz_plane_layout.addWidget(self.p3_xz_plane_x); p3_xz_plane_layout.addWidget(self.p3_xz_plane_y); p3_xz_plane_layout.addWidget(self.p3_xz_plane_z)

        form.addRow("Origin (X,Y,Z):", p3_origin_layout)
        form.addRow("Point on Z-Axis (X,Y,Z):", p3_z_axis_layout)
        form.addRow("Point on XZ-Plane (X,Y,Z):", p3_xz_plane_layout)
        return panel

    def _create_rotate_translate_panel(self):
        """Creates the widget with inputs for the 'Rotate and Translate' method."""
        panel = QWidget()
        form = QFormLayout(panel)
        form.setContentsMargins(0, 10, 0, 0)

        self.ref_cid_combo = QComboBox()
        for cid in sorted(self.model.coords.keys()):
            self.ref_cid_combo.addItem(str(cid))

        self.trans_dx, self.trans_dy, self.trans_dz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")
        self.rot_rx, self.rot_ry, self.rot_rz = QLineEdit("0.0"), QLineEdit("0.0"), QLineEdit("0.0")

        for w in [self.trans_dx, self.trans_dy, self.trans_dz, self.rot_rx, self.rot_ry, self.rot_rz]:
            w.setValidator(QDoubleValidator())

        trans_layout = QHBoxLayout(); trans_layout.addWidget(self.trans_dx); trans_layout.addWidget(self.trans_dy); trans_layout.addWidget(self.trans_dz)
        rot_layout = QHBoxLayout(); rot_layout.addWidget(self.rot_rx); rot_layout.addWidget(self.rot_ry); rot_layout.addWidget(self.rot_rz)

        form.addRow("Reference System (CID):", self.ref_cid_combo)
        form.addRow("Translate (dX,dY,dZ):", trans_layout)
        form.addRow("Rotate (rX,rY,rZ) deg:", rot_layout)
        return panel

    def get_parameters(self):
        try:
            cid = int(self.cid_input.text())
            if cid in (0, 1, 2):
                QMessageBox.warning(self, "Reserved CID",
                                    "CIDs 0, 1, and 2 are predefined global coordinate systems and cannot be used.")
                return None
            params = {
                'cid': cid,
                'comment': f"$ {self.title_input.text()}" if self.title_input.text().strip() else '',
                'type': self.type_combo.currentText().lower(),
                'method': self.method_combo.currentText()
            }
            if self.method_combo.currentText() == "3 Points":
                params['origin'] = [float(self.p3_origin_x.text()), float(self.p3_origin_y.text()), float(self.p3_origin_z.text())]
                params['z_axis_point'] = [float(self.p3_z_axis_x.text()), float(self.p3_z_axis_y.text()), float(self.p3_z_axis_z.text())]
                params['xz_plane_point'] = [float(self.p3_xz_plane_x.text()), float(self.p3_xz_plane_y.text()), float(self.p3_xz_plane_z.text())]
                origin_v, z_axis_v, xz_plane_v = np.array(params['origin']), np.array(params['z_axis_point']), np.array(params['xz_plane_point'])
                if np.allclose(origin_v, z_axis_v) or np.linalg.norm(np.cross(z_axis_v - origin_v, xz_plane_v - origin_v)) < 1e-9:
                    QMessageBox.warning(self, "Invalid Input", "The three defining points for the coordinate system cannot be collinear.")
                    return None

            elif self.method_combo.currentText() == "Rotate and Translate":
                params['ref_cid'] = int(self.ref_cid_combo.currentText())
                params['translations'] = [float(self.trans_dx.text()), float(self.trans_dy.text()), float(self.trans_dz.text())]
                params['rotations'] = [float(self.rot_rx.text()), float(self.rot_ry.text()), float(self.rot_rz.text())]

            return params
        except (ValueError, KeyError) as e:
            QMessageBox.warning(self, "Input Error", f"Please provide valid numerical inputs. Error: {e}")
            return None


class BeamSectionLibraryDialog(QDialog):
    """Dialog for selecting standard beam cross-section shapes and computing A/I1/I2/J."""

    SHAPES = [
        'I / H Beam', 'C Channel', 'Box (Rect Tube)', 'Round Tube',
        'Solid Rectangle', 'Solid Circle',
        'L Section', 'T Section', 'Z Section',
    ]

    # (label, key, default) tuples per shape
    SHAPE_DIMS = {
        'I / H Beam': [('Height (H)', 'H', '1.0'), ('Top Flange Width', 'W_top', '0.5'),
                       ('Bottom Flange Width', 'W_bot', '0.5'), ('Web Thickness', 't_web', '0.05'),
                       ('Top Flange Thickness', 't_ftop', '0.08'), ('Bottom Flange Thickness', 't_fbot', '0.08')],
        'C Channel': [('Height (H)', 'H', '1.0'), ('Flange Width (W)', 'W', '0.4'),
                      ('Web Thickness', 't_web', '0.05'), ('Flange Thickness', 't_flange', '0.08')],
        'Box (Rect Tube)': [('Height (H)', 'H', '1.0'), ('Width (W)', 'W', '0.5'),
                            ('Web Thickness', 't_web', '0.05'), ('Flange Thickness', 't_flange', '0.05')],
        'Round Tube': [('Outer Radius (R)', 'R', '0.5'), ('Wall Thickness (t)', 't', '0.05')],
        'Solid Rectangle': [('Height (H)', 'H', '1.0'), ('Width (W)', 'W', '0.5')],
        'Solid Circle': [('Radius (R)', 'R', '0.5')],
        'L Section': [('Height (H)', 'H', '1.0'), ('Width (W)', 'W', '0.5'), ('Thickness (t)', 't', '0.05')],
        'T Section': [('Flange Width (W)', 'W', '0.5'), ('Flange Thickness (tf)', 'tf', '0.08'),
                      ('Web Height (H)', 'H', '1.0'), ('Web Thickness (tw)', 'tw', '0.05')],
        'Z Section': [('Height (H)', 'H', '1.0'), ('Width (W)', 'W', '0.4'), ('Thickness (t)', 't', '0.05')],
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Beam Section Library")
        self.setMinimumSize(640, 480)
        self._computed = {'A': 0., 'I1': 0., 'I2': 0., 'J': 0.}
        self._dim_inputs = {}

        layout = QVBoxLayout(self)

        top = QHBoxLayout()
        top.addWidget(QLabel("Section Shape:"))
        self.shape_combo = QComboBox()
        self.shape_combo.addItems(self.SHAPES)
        self.shape_combo.currentTextChanged.connect(self._on_shape_changed)
        top.addWidget(self.shape_combo, 1)
        layout.addLayout(top)

        body = QHBoxLayout()

        # Left: dimension inputs (stacked)
        self._dims_stack = QStackedLayout()
        dims_container = QWidget()
        dims_container.setLayout(self._dims_stack)
        for shape in self.SHAPES:
            page = QWidget()
            form = QFormLayout(page)
            self._dim_inputs[shape] = {}
            for label, key, default in self.SHAPE_DIMS[shape]:
                le = QLineEdit(default)
                le.setValidator(QDoubleValidator(0, 1e20, 8))
                le.textChanged.connect(self._recalculate)
                form.addRow(label + ":", le)
                self._dim_inputs[shape][key] = le
            self._dims_stack.addWidget(page)
        body.addWidget(dims_container, 1)

        # Right: matplotlib preview
        try:
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
            from matplotlib.figure import Figure
            self._fig = Figure(figsize=(3, 3))
            self._ax = self._fig.add_subplot(111)
            self._canvas = FigureCanvasQTAgg(self._fig)
            body.addWidget(self._canvas, 1)
            self._has_mpl = True
        except ImportError:
            self._has_mpl = False
            body.addWidget(QLabel("(matplotlib not available for preview)"), 1)

        layout.addLayout(body)

        # Results display
        results_group = QGroupBox("Computed Properties")
        rform = QFormLayout(results_group)
        self._result_A = QLineEdit(); self._result_A.setReadOnly(True)
        self._result_I1 = QLineEdit(); self._result_I1.setReadOnly(True)
        self._result_I2 = QLineEdit(); self._result_I2.setReadOnly(True)
        self._result_J = QLineEdit(); self._result_J.setReadOnly(True)
        rform.addRow("Area (A):", self._result_A)
        rform.addRow("I1 (strong axis):", self._result_I1)
        rform.addRow("I2 (weak axis):", self._result_I2)
        rform.addRow("J (torsion):", self._result_J)
        layout.addWidget(results_group)

        btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        layout.addWidget(btn_box)

        self._on_shape_changed(self.shape_combo.currentText())

    def _on_shape_changed(self, shape_text):
        idx = self.SHAPES.index(shape_text)
        self._dims_stack.setCurrentIndex(idx)
        self._recalculate()

    def _recalculate(self):
        from node_runner.model import NastranModelGenerator as NMG
        shape = self.shape_combo.currentText()
        dims = {}
        try:
            for key, le in self._dim_inputs[shape].items():
                val = float(le.text()) if le.text() else 0.
                if val <= 0:
                    raise ValueError
                dims[key] = val
        except (ValueError, KeyError):
            self._result_A.setText("---")
            self._result_I1.setText("---")
            self._result_I2.setText("---")
            self._result_J.setText("---")
            return

        calc_map = {
            'I / H Beam': NMG.calc_i_section_props,
            'C Channel': NMG.calc_c_section_props,
            'Box (Rect Tube)': NMG.calc_box_section_props,
            'Round Tube': NMG.calc_tube_section_props,
            'Solid Rectangle': NMG.calc_solid_rect_props,
            'Solid Circle': NMG.calc_solid_circle_props,
        }

        if shape in calc_map:
            props = calc_map[shape](dims)
        else:
            gen = NMG.__new__(NMG)
            props = gen._calculate_beam_props(shape[0], dims)

        self._computed = props
        self._result_A.setText(f"{props['A']:.6g}")
        self._result_I1.setText(f"{props['I1']:.6g}")
        self._result_I2.setText(f"{props['I2']:.6g}")
        self._result_J.setText(f"{props['J']:.6g}")

        if self._has_mpl:
            self._draw_preview(shape, dims)

    def _draw_preview(self, shape, dims):
        ax = self._ax
        ax.clear()
        ax.set_aspect('equal')
        ax.set_facecolor('#313244')
        self._fig.set_facecolor('#1e1e2e')
        color = '#89b4fa'

        if shape == 'I / H Beam':
            H = dims['H']; Wt = dims['W_top']; Wb = dims['W_bot']
            tw = dims['t_web']; tft = dims['t_ftop']; tfb = dims['t_fbot']
            # Bottom flange
            ax.fill([-Wb/2, Wb/2, Wb/2, -Wb/2], [0, 0, tfb, tfb], color=color)
            # Web
            ax.fill([-tw/2, tw/2, tw/2, -tw/2], [tfb, tfb, H-tft, H-tft], color=color)
            # Top flange
            ax.fill([-Wt/2, Wt/2, Wt/2, -Wt/2], [H-tft, H-tft, H, H], color=color)

        elif shape == 'C Channel':
            H = dims['H']; W = dims['W']; tw = dims['t_web']; tf = dims['t_flange']
            from matplotlib.patches import Polygon
            pts = [(0,0),(W,0),(W,tf),(tw,tf),(tw,H-tf),(W,H-tf),(W,H),(0,H)]
            ax.add_patch(Polygon(pts, closed=True, fc=color, ec=color))
            ax.set_xlim(-0.1*W, W*1.2); ax.set_ylim(-0.1*H, H*1.1)

        elif shape == 'Box (Rect Tube)':
            H = dims['H']; W = dims['W']; tw = dims['t_web']; tf = dims['t_flange']
            from matplotlib.patches import Rectangle as MplRect
            ax.add_patch(MplRect((-W/2, 0), W, H, fc=color, ec=color))
            ax.add_patch(MplRect((-W/2+tw, tf), W-2*tw, H-2*tf, fc='#313244', ec='#313244'))
            ax.set_xlim(-W, W); ax.set_ylim(-0.2*H, H*1.2)

        elif shape == 'Round Tube':
            R = dims['R']; t = dims['t']
            theta = np.linspace(0, 2*np.pi, 80)
            ax.fill(R*np.cos(theta), R*np.sin(theta), color=color)
            ax.fill((R-t)*np.cos(theta), (R-t)*np.sin(theta), color='#313244')

        elif shape == 'Solid Rectangle':
            H = dims['H']; W = dims['W']
            ax.fill([-W/2, W/2, W/2, -W/2], [0, 0, H, H], color=color)

        elif shape == 'Solid Circle':
            R = dims['R']
            theta = np.linspace(0, 2*np.pi, 80)
            ax.fill(R*np.cos(theta), R*np.sin(theta), color=color)

        elif shape == 'L Section':
            H = dims['H']; W = dims['W']; t = dims['t']
            pts = [(0,0),(W,0),(W,t),(t,t),(t,H),(0,H)]
            from matplotlib.patches import Polygon
            ax.add_patch(Polygon(pts, closed=True, fc=color, ec=color))
            ax.set_xlim(-0.1*W, W*1.2); ax.set_ylim(-0.1*H, H*1.1)

        elif shape == 'T Section':
            W = dims['W']; tf = dims['tf']; H = dims['H']; tw = dims['tw']
            ax.fill([-W/2, W/2, W/2, -W/2], [H, H, H+tf, H+tf], color=color)
            ax.fill([-tw/2, tw/2, tw/2, -tw/2], [0, 0, H, H], color=color)

        elif shape == 'Z Section':
            H = dims['H']; W = dims['W']; t = dims['t']
            from matplotlib.patches import Polygon
            pts = [(0,0),(W,0),(W,t),(t,t),(t,H-t),(W,H-t),(W,H),(0,H),(0,H-t),(-W+t,H-t),(-W+t,H-2*t),(0,H-2*t)]
            # Simplified Z shape
            ax.add_patch(Polygon([(0,0),(W,0),(W,t),(t,t),(t,H),(0,H)], closed=True, fc=color, ec=color))
            ax.add_patch(Polygon([(-W+t,H-t),(0,H-t),(0,H),(-W+t,H)], closed=True, fc=color, ec=color))
            ax.set_xlim(-W*1.2, W*1.2); ax.set_ylim(-0.1*H, H*1.1)

        ax.autoscale_view()
        ax.tick_params(colors='#a6adc8', labelsize=7)
        for spine in ax.spines.values():
            spine.set_color('#585b70')
        self._canvas.draw_idle()

    def get_properties(self):
        return dict(self._computed)


class CreateSpiderDialog(QDialog):
    """Dialog for creating a spider connection (RBE2/RBE3).

    The user specifies ring node IDs, RBE type, and DOF components.
    A center node is auto-computed at the centroid of ring nodes.
    """

    request_picking_mode = QtCore.Signal(str, bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Spider Connection")
        self.setMinimumWidth(400)
        layout = QVBoxLayout(self)

        # RBE Type
        type_group = QGroupBox("Connection Type")
        type_layout = QFormLayout(type_group)
        self.rbe_type_combo = QComboBox()
        self.rbe_type_combo.addItems(["RBE2 (Rigid)", "RBE3 (Interpolation)"])
        type_layout.addRow("Type:", self.rbe_type_combo)
        layout.addWidget(type_group)

        # DOF
        dof_group = QGroupBox("DOF Components")
        dof_layout = QFormLayout(dof_group)
        self.dof_input = QLineEdit("123456")
        dof_layout.addRow("DOF:", self.dof_input)
        layout.addWidget(dof_group)

        # Ring nodes
        nodes_group = QGroupBox("Ring Nodes")
        nodes_layout = QVBoxLayout(nodes_group)
        self.nodes_input = QLineEdit()
        self.nodes_input.setPlaceholderText("Enter node IDs (comma separated)")
        nodes_layout.addWidget(self.nodes_input)
        pick_btn = QPushButton("Pick from Viewer")
        pick_btn.clicked.connect(self._pick_nodes)
        nodes_layout.addWidget(pick_btn)
        self.node_count_label = QLabel("0 nodes selected")
        nodes_layout.addWidget(self.node_count_label)
        layout.addWidget(nodes_group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _pick_nodes(self):
        self.hide()
        self.request_picking_mode.emit('Node', True)

    def set_picked_nodes(self, node_ids):
        """Called by mainwindow when nodes are picked."""
        self.nodes_input.setText(", ".join(str(n) for n in node_ids))
        self.node_count_label.setText(f"{len(node_ids)} nodes selected")

    def get_settings(self):
        import re
        rbe_type = 'RBE2' if self.rbe_type_combo.currentIndex() == 0 else 'RBE3'
        dof = self.dof_input.text().strip() or '123456'
        nids = [int(x) for x in re.findall(r'\d+', self.nodes_input.text())]
        return {'rbe_type': rbe_type, 'dof': dof, 'ring_nids': nids}


class CreateWeldDialog(QDialog):
    """Dialog for creating a weld/fastener connection (CWELD-like).

    The user specifies two node IDs (one on each surface), diameter,
    and property ID.
    """

    picking_target = None  # callback: set by mainwindow

    def __init__(self, available_pids=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Weld/Fastener Connection")
        self.setMinimumWidth(400)
        layout = QVBoxLayout(self)

        # Nodes
        conn_group = QGroupBox("Connection Nodes")
        conn_layout = QFormLayout(conn_group)
        self.nid_a_input = QLineEdit()
        self.nid_a_input.setPlaceholderText("Node on surface A")
        self.nid_b_input = QLineEdit()
        self.nid_b_input.setPlaceholderText("Node on surface B")
        conn_layout.addRow("Node A:", self.nid_a_input)
        conn_layout.addRow("Node B:", self.nid_b_input)
        layout.addWidget(conn_group)

        # Parameters
        param_group = QGroupBox("Parameters")
        param_layout = QFormLayout(param_group)
        self.diameter_input = QLineEdit("6.0")
        self.diameter_input.setValidator(QDoubleValidator(0.001, 1e6, 3, self))
        param_layout.addRow("Diameter:", self.diameter_input)
        self.pid_input = QLineEdit("1")
        if available_pids:
            next_pid = max(available_pids) + 1
            self.pid_input.setText(str(next_pid))
        param_layout.addRow("Property ID:", self.pid_input)
        layout.addWidget(param_group)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_settings(self):
        try:
            nid_a = int(self.nid_a_input.text())
        except ValueError:
            nid_a = 0
        try:
            nid_b = int(self.nid_b_input.text())
        except ValueError:
            nid_b = 0
        try:
            diameter = float(self.diameter_input.text())
        except ValueError:
            diameter = 6.0
        try:
            pid = int(self.pid_input.text())
        except ValueError:
            pid = 1
        return {'nid_a': nid_a, 'nid_b': nid_b, 'diameter': diameter, 'pid': pid}
