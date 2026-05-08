"""Geometry creation dialogs and mesh-on-geometry dialogs."""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QComboBox,
    QDialogButtonBox, QListWidget, QListWidgetItem,
    QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView,
    QMessageBox, QHBoxLayout,
)
from PySide6.QtGui import QDoubleValidator, QIntValidator


class CreateGeometryPointDialog(QDialog):
    """Dialog for creating a geometry point, with optional Excel batch import."""

    def __init__(self, next_id, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Geometry Point")
        self.setMinimumWidth(450)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.id_input = QLineEdit(str(next_id))
        self.id_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("Point ID:", self.id_input)

        self.x_input = QLineEdit("0.0")
        self.x_input.setValidator(QDoubleValidator())
        form.addRow("X:", self.x_input)

        self.y_input = QLineEdit("0.0")
        self.y_input.setValidator(QDoubleValidator())
        form.addRow("Y:", self.y_input)

        self.z_input = QLineEdit("0.0")
        self.z_input.setValidator(QDoubleValidator())
        form.addRow("Z:", self.z_input)

        layout.addLayout(form)

        # Preview table (hidden until Excel loaded)
        self.preview_table = QTableWidget()
        self.preview_table.setColumnCount(4)
        self.preview_table.setHorizontalHeaderLabels(["ID", "X", "Y", "Z"])
        self.preview_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.preview_table.setVisible(False)
        layout.addWidget(self.preview_table)

        self._batch_points = []

        btn_row = QHBoxLayout()
        excel_btn = QPushButton("From Excel...")
        excel_btn.clicked.connect(self._load_excel)
        btn_row.addWidget(excel_btn)
        btn_row.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(ok_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

    def _load_excel(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Open Excel File", "", "Excel Files (*.xlsx);;All Files (*)")
        if not filepath:
            return
        try:
            import openpyxl
            wb = openpyxl.load_workbook(filepath, read_only=True, data_only=True)
            ws = wb.active
            self._batch_points.clear()
            rows = list(ws.iter_rows(min_row=2, values_only=True))  # Skip header
            for row in rows:
                if not row or all(v is None for v in row):
                    continue
                # Columns: ID (optional), CSys (ignored), X, Y, Z
                point_id = int(row[0]) if row[0] is not None else None
                # row[1] is CSys - ignored for now
                x = float(row[2]) if len(row) > 2 and row[2] is not None else 0.0
                y = float(row[3]) if len(row) > 3 and row[3] is not None else 0.0
                z = float(row[4]) if len(row) > 4 and row[4] is not None else 0.0
                self._batch_points.append({'id': point_id, 'x': x, 'y': y, 'z': z})
            wb.close()

            # Populate preview table
            self.preview_table.setRowCount(len(self._batch_points))
            for i, pt in enumerate(self._batch_points):
                self.preview_table.setItem(i, 0, QTableWidgetItem(
                    str(pt['id']) if pt['id'] else "Auto"))
                self.preview_table.setItem(i, 1, QTableWidgetItem(str(pt['x'])))
                self.preview_table.setItem(i, 2, QTableWidgetItem(str(pt['y'])))
                self.preview_table.setItem(i, 3, QTableWidgetItem(str(pt['z'])))
            self.preview_table.setVisible(True)
            # Disable single-point form fields
            self.id_input.setEnabled(False)
            self.x_input.setEnabled(False)
            self.y_input.setEnabled(False)
            self.z_input.setEnabled(False)
        except ImportError:
            QMessageBox.critical(self, "Missing Dependency",
                                 "openpyxl is required for Excel import.\n"
                                 "Install it with: pip install openpyxl")
        except Exception as e:
            QMessageBox.critical(self, "Excel Error", f"Failed to read file: {e}")

    @property
    def is_batch(self):
        return bool(self._batch_points)

    def get_batch_parameters(self):
        return list(self._batch_points)

    def get_parameters(self):
        return {
            'id': int(self.id_input.text()),
            'x': float(self.x_input.text()),
            'y': float(self.y_input.text()),
            'z': float(self.z_input.text()),
        }


class CreateGeometryLineDialog(QDialog):
    """Dialog for creating a geometry line between two points."""

    def __init__(self, next_id, point_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Geometry Line")

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.id_input = QLineEdit(str(next_id))
        self.id_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("Line ID:", self.id_input)

        self.start_combo = QComboBox()
        self.end_combo = QComboBox()
        for pid in point_ids:
            self.start_combo.addItem(f"P{pid}", pid)
            self.end_combo.addItem(f"P{pid}", pid)
        if len(point_ids) > 1:
            self.end_combo.setCurrentIndex(1)

        form.addRow("Start Point:", self.start_combo)
        form.addRow("End Point:", self.end_combo)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'id': int(self.id_input.text()),
            'start_point_id': self.start_combo.currentData(),
            'end_point_id': self.end_combo.currentData(),
        }


class CreateGeometryArcDialog(QDialog):
    """Dialog for creating a geometry arc through three points."""

    def __init__(self, next_id, point_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Geometry Arc")

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.id_input = QLineEdit(str(next_id))
        self.id_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("Arc ID:", self.id_input)

        self.start_combo = QComboBox()
        self.mid_combo = QComboBox()
        self.end_combo = QComboBox()
        for pid in point_ids:
            self.start_combo.addItem(f"P{pid}", pid)
            self.mid_combo.addItem(f"P{pid}", pid)
            self.end_combo.addItem(f"P{pid}", pid)
        if len(point_ids) > 1:
            self.mid_combo.setCurrentIndex(1)
        if len(point_ids) > 2:
            self.end_combo.setCurrentIndex(2)

        form.addRow("Start Point:", self.start_combo)
        form.addRow("Mid Point:", self.mid_combo)
        form.addRow("End Point:", self.end_combo)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'id': int(self.id_input.text()),
            'start_point_id': self.start_combo.currentData(),
            'mid_point_id': self.mid_combo.currentData(),
            'end_point_id': self.end_combo.currentData(),
        }


class CreateGeometryCircleDialog(QDialog):
    """Dialog for creating a geometry circle."""

    def __init__(self, next_id, point_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Geometry Circle")

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.id_input = QLineEdit(str(next_id))
        self.id_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("Circle ID:", self.id_input)

        self.center_combo = QComboBox()
        for pid in point_ids:
            self.center_combo.addItem(f"P{pid}", pid)
        form.addRow("Center Point:", self.center_combo)

        self.radius_input = QLineEdit("1.0")
        self.radius_input.setValidator(QDoubleValidator(0.0001, 1e12, 6))
        form.addRow("Radius:", self.radius_input)

        normal_group = QGroupBox("Plane Normal")
        normal_form = QFormLayout(normal_group)
        self.nx_input = QLineEdit("0.0")
        self.nx_input.setValidator(QDoubleValidator())
        self.ny_input = QLineEdit("0.0")
        self.ny_input.setValidator(QDoubleValidator())
        self.nz_input = QLineEdit("1.0")
        self.nz_input.setValidator(QDoubleValidator())
        normal_form.addRow("NX:", self.nx_input)
        normal_form.addRow("NY:", self.ny_input)
        normal_form.addRow("NZ:", self.nz_input)

        layout.addLayout(form)
        layout.addWidget(normal_group)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'id': int(self.id_input.text()),
            'center_point_id': self.center_combo.currentData(),
            'radius': float(self.radius_input.text()),
            'nx': float(self.nx_input.text()),
            'ny': float(self.ny_input.text()),
            'nz': float(self.nz_input.text()),
        }


class CreateGeometrySurfaceDialog(QDialog):
    """Dialog for creating a geometry surface from boundary curves."""

    def __init__(self, next_id, curve_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create Geometry Surface")

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.id_input = QLineEdit(str(next_id))
        self.id_input.setValidator(QIntValidator(1, 999999999))
        form.addRow("Surface ID:", self.id_input)
        layout.addLayout(form)

        layout.addWidget(QLabel("Select boundary curves (in order):"))
        self.curve_list = QListWidget()
        self.curve_list.setSelectionMode(QListWidget.MultiSelection)
        for cid in curve_ids:
            item = QListWidgetItem(f"Curve {cid}")
            item.setData(256, cid)  # Qt.UserRole = 256
            self.curve_list.addItem(item)
        layout.addWidget(self.curve_list)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        selected_ids = [item.data(256) for item in self.curve_list.selectedItems()]
        return {
            'id': int(self.id_input.text()),
            'curve_ids': selected_ids,
        }


class MeshCurveDialog(QDialog):
    """Dialog for meshing along geometry curves."""

    def __init__(self, curve_ids, property_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mesh Curve")

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Select curves to mesh:"))
        self.curve_list = QListWidget()
        self.curve_list.setSelectionMode(QListWidget.MultiSelection)
        for cid in curve_ids:
            item = QListWidgetItem(f"Curve {cid}")
            item.setData(256, cid)
            self.curve_list.addItem(item)
        layout.addWidget(self.curve_list)

        form = QFormLayout()
        self.n_elements_input = QLineEdit("10")
        self.n_elements_input.setValidator(QIntValidator(1, 10000))
        form.addRow("Elements per curve:", self.n_elements_input)

        self.elem_type_combo = QComboBox()
        self.elem_type_combo.addItems(["CBAR", "CBEAM", "CROD"])
        form.addRow("Element Type:", self.elem_type_combo)

        self.property_combo = QComboBox()
        for pid in property_ids:
            self.property_combo.addItem(f"PID {pid}", pid)
        form.addRow("Property:", self.property_combo)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        selected_ids = [item.data(256) for item in self.curve_list.selectedItems()]
        return {
            'curve_ids': selected_ids,
            'n_elements': int(self.n_elements_input.text()),
            'elem_type': self.elem_type_combo.currentText(),
            'pid': self.property_combo.currentData(),
        }


class MeshSurfaceDialog(QDialog):
    """Dialog for meshing a geometry surface with gmsh."""

    def __init__(self, surface_ids, property_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mesh Surface")

        layout = QVBoxLayout(self)

        form = QFormLayout()
        self.surface_combo = QComboBox()
        for sid in surface_ids:
            self.surface_combo.addItem(f"Surface {sid}", sid)
        form.addRow("Surface:", self.surface_combo)

        self.n_per_edge_input = QLineEdit("10")
        self.n_per_edge_input.setValidator(QIntValidator(1, 1000))
        form.addRow("Elements per edge:", self.n_per_edge_input)

        self.elem_pref_combo = QComboBox()
        self.elem_pref_combo.addItems(["Quad-Dominant", "All Triangles"])
        form.addRow("Element type:", self.elem_pref_combo)

        self.property_combo = QComboBox()
        for pid in property_ids:
            self.property_combo.addItem(f"PID {pid}", pid)
        form.addRow("Property:", self.property_combo)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        pref = 'quad' if self.elem_pref_combo.currentIndex() == 0 else 'tri'
        return {
            'surface_id': self.surface_combo.currentData(),
            'n_per_edge': int(self.n_per_edge_input.text()),
            'elem_preference': pref,
            'pid': self.property_combo.currentData(),
        }


# ---------------------------------------------------------------------------
# Geometry edit dialogs
# ---------------------------------------------------------------------------

class EditGeometryPointDialog(QDialog):
    """Dialog for editing a geometry point's coordinates."""

    def __init__(self, point_id, current_xyz, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Geometry Point P{point_id}")
        self.setMinimumWidth(350)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        id_label = QLineEdit(str(point_id))
        id_label.setEnabled(False)
        form.addRow("Point ID:", id_label)

        self.x_input = QLineEdit(str(current_xyz[0]))
        self.x_input.setValidator(QDoubleValidator())
        form.addRow("X:", self.x_input)

        self.y_input = QLineEdit(str(current_xyz[1]))
        self.y_input.setValidator(QDoubleValidator())
        form.addRow("Y:", self.y_input)

        self.z_input = QLineEdit(str(current_xyz[2]))
        self.z_input.setValidator(QDoubleValidator())
        form.addRow("Z:", self.z_input)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'x': float(self.x_input.text()),
            'y': float(self.y_input.text()),
            'z': float(self.z_input.text()),
        }


class EditGeometryCurveDialog(QDialog):
    """Dialog for editing a geometry curve's references/parameters."""

    def __init__(self, curve_id, curve_type, point_ids, current_data, parent=None):
        """
        Parameters
        ----------
        curve_id : int
        curve_type : str  ('line', 'arc', 'circle')
        point_ids : list[int]  available point IDs
        current_data : dict  current curve parameters
        """
        super().__init__(parent)
        self.setWindowTitle(f"Edit Geometry Curve {curve_id}")
        self.setMinimumWidth(400)
        self._curve_type = curve_type

        layout = QVBoxLayout(self)
        form = QFormLayout()

        id_label = QLineEdit(str(curve_id))
        id_label.setEnabled(False)
        form.addRow("Curve ID:", id_label)

        type_label = QLineEdit(curve_type.capitalize())
        type_label.setEnabled(False)
        form.addRow("Type:", type_label)

        if curve_type == 'line':
            self.start_combo = QComboBox()
            self.end_combo = QComboBox()
            for pid in point_ids:
                self.start_combo.addItem(f"P{pid}", pid)
                self.end_combo.addItem(f"P{pid}", pid)
            _set_combo_by_data(self.start_combo, current_data['start_point_id'])
            _set_combo_by_data(self.end_combo, current_data['end_point_id'])
            form.addRow("Start Point:", self.start_combo)
            form.addRow("End Point:", self.end_combo)

        elif curve_type == 'arc':
            self.start_combo = QComboBox()
            self.mid_combo = QComboBox()
            self.end_combo = QComboBox()
            for pid in point_ids:
                self.start_combo.addItem(f"P{pid}", pid)
                self.mid_combo.addItem(f"P{pid}", pid)
                self.end_combo.addItem(f"P{pid}", pid)
            _set_combo_by_data(self.start_combo, current_data['start_point_id'])
            _set_combo_by_data(self.mid_combo, current_data['mid_point_id'])
            _set_combo_by_data(self.end_combo, current_data['end_point_id'])
            form.addRow("Start Point:", self.start_combo)
            form.addRow("Mid Point:", self.mid_combo)
            form.addRow("End Point:", self.end_combo)

        elif curve_type == 'circle':
            self.center_combo = QComboBox()
            for pid in point_ids:
                self.center_combo.addItem(f"P{pid}", pid)
            _set_combo_by_data(self.center_combo, current_data['center_point_id'])
            form.addRow("Center Point:", self.center_combo)

            self.radius_input = QLineEdit(str(current_data['radius']))
            self.radius_input.setValidator(QDoubleValidator(0.0001, 1e12, 6))
            form.addRow("Radius:", self.radius_input)

            normal = current_data['normal']
            self.nx_input = QLineEdit(str(normal[0]))
            self.nx_input.setValidator(QDoubleValidator())
            self.ny_input = QLineEdit(str(normal[1]))
            self.ny_input.setValidator(QDoubleValidator())
            self.nz_input = QLineEdit(str(normal[2]))
            self.nz_input.setValidator(QDoubleValidator())
            form.addRow("Normal X:", self.nx_input)
            form.addRow("Normal Y:", self.ny_input)
            form.addRow("Normal Z:", self.nz_input)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        if self._curve_type == 'line':
            return {
                'start_point_id': self.start_combo.currentData(),
                'end_point_id': self.end_combo.currentData(),
            }
        elif self._curve_type == 'arc':
            return {
                'start_point_id': self.start_combo.currentData(),
                'mid_point_id': self.mid_combo.currentData(),
                'end_point_id': self.end_combo.currentData(),
            }
        elif self._curve_type == 'circle':
            return {
                'center_point_id': self.center_combo.currentData(),
                'radius': float(self.radius_input.text()),
                'normal': [float(self.nx_input.text()),
                           float(self.ny_input.text()),
                           float(self.nz_input.text())],
            }


class EditGeometrySurfaceDialog(QDialog):
    """Dialog for editing a geometry surface's boundary curves."""

    def __init__(self, surface_id, curve_ids, current_boundary_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Geometry Surface {surface_id}")
        self.setMinimumWidth(400)

        layout = QVBoxLayout(self)

        id_label = QLineEdit(str(surface_id))
        id_label.setEnabled(False)
        form = QFormLayout()
        form.addRow("Surface ID:", id_label)
        layout.addLayout(form)

        layout.addWidget(QLabel("Select boundary curves (in order):"))
        self.curve_list = QListWidget()
        self.curve_list.setSelectionMode(QListWidget.MultiSelection)
        for cid in curve_ids:
            item = QListWidgetItem(f"Curve {cid}")
            item.setData(256, cid)
            self.curve_list.addItem(item)
            if cid in current_boundary_ids:
                item.setSelected(True)
        layout.addWidget(self.curve_list)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'curve_ids': [item.data(256)
                          for item in self.curve_list.selectedItems()],
        }


# ---------------------------------------------------------------------------
# Geometry transform dialog
# ---------------------------------------------------------------------------

class TransformGeometryDialog(QDialog):
    """Dialog for translating, rotating, mirroring, or scaling geometry points."""

    def __init__(self, point_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Transform Geometry")
        self.setMinimumWidth(400)

        from PySide6.QtWidgets import QTabWidget, QSpinBox, QWidget, QRadioButton

        layout = QVBoxLayout(self)

        # Point selection summary
        layout.addWidget(QLabel(f"Transforming {len(point_ids)} point(s)"))

        self.tabs = QTabWidget()

        # --- Translate tab ---
        translate_tab = QWidget()
        t_layout = QFormLayout(translate_tab)
        self.dx_input = QLineEdit("0.0")
        self.dx_input.setValidator(QDoubleValidator())
        self.dy_input = QLineEdit("0.0")
        self.dy_input.setValidator(QDoubleValidator())
        self.dz_input = QLineEdit("0.0")
        self.dz_input.setValidator(QDoubleValidator())
        t_layout.addRow("Delta X:", self.dx_input)
        t_layout.addRow("Delta Y:", self.dy_input)
        t_layout.addRow("Delta Z:", self.dz_input)
        self.tabs.addTab(translate_tab, "Translate")

        # --- Rotate tab ---
        rotate_tab = QWidget()
        r_layout = QFormLayout(rotate_tab)
        self.rot_angle_input = QLineEdit("0.0")
        self.rot_angle_input.setValidator(QDoubleValidator())
        r_layout.addRow("Angle (deg):", self.rot_angle_input)
        self.rot_axis_combo = QComboBox()
        self.rot_axis_combo.addItems(["X", "Y", "Z"])
        self.rot_axis_combo.setCurrentIndex(2)
        r_layout.addRow("Axis:", self.rot_axis_combo)
        self.rot_cx = QLineEdit("0.0")
        self.rot_cx.setValidator(QDoubleValidator())
        self.rot_cy = QLineEdit("0.0")
        self.rot_cy.setValidator(QDoubleValidator())
        self.rot_cz = QLineEdit("0.0")
        self.rot_cz.setValidator(QDoubleValidator())
        r_layout.addRow("Center X:", self.rot_cx)
        r_layout.addRow("Center Y:", self.rot_cy)
        r_layout.addRow("Center Z:", self.rot_cz)
        self.tabs.addTab(rotate_tab, "Rotate")

        # --- Mirror tab ---
        mirror_tab = QWidget()
        m_layout = QFormLayout(mirror_tab)
        self.mirror_plane_combo = QComboBox()
        self.mirror_plane_combo.addItems(["XY Plane", "XZ Plane", "YZ Plane", "Custom"])
        m_layout.addRow("Plane:", self.mirror_plane_combo)
        self.mirror_px = QLineEdit("0.0")
        self.mirror_px.setValidator(QDoubleValidator())
        self.mirror_py = QLineEdit("0.0")
        self.mirror_py.setValidator(QDoubleValidator())
        self.mirror_pz = QLineEdit("0.0")
        self.mirror_pz.setValidator(QDoubleValidator())
        m_layout.addRow("Plane Point X:", self.mirror_px)
        m_layout.addRow("Plane Point Y:", self.mirror_py)
        m_layout.addRow("Plane Point Z:", self.mirror_pz)
        self.mirror_nx = QLineEdit("0.0")
        self.mirror_nx.setValidator(QDoubleValidator())
        self.mirror_ny = QLineEdit("0.0")
        self.mirror_ny.setValidator(QDoubleValidator())
        self.mirror_nz = QLineEdit("1.0")
        self.mirror_nz.setValidator(QDoubleValidator())
        m_layout.addRow("Normal X:", self.mirror_nx)
        m_layout.addRow("Normal Y:", self.mirror_ny)
        m_layout.addRow("Normal Z:", self.mirror_nz)
        self.tabs.addTab(mirror_tab, "Mirror")

        # --- Scale tab ---
        scale_tab = QWidget()
        s_layout = QFormLayout(scale_tab)
        self.scale_factor_input = QLineEdit("1.0")
        self.scale_factor_input.setValidator(QDoubleValidator(0.0001, 1e12, 6))
        s_layout.addRow("Scale Factor:", self.scale_factor_input)
        self.scale_cx = QLineEdit("0.0")
        self.scale_cx.setValidator(QDoubleValidator())
        self.scale_cy = QLineEdit("0.0")
        self.scale_cy.setValidator(QDoubleValidator())
        self.scale_cz = QLineEdit("0.0")
        self.scale_cz.setValidator(QDoubleValidator())
        s_layout.addRow("Center X:", self.scale_cx)
        s_layout.addRow("Center Y:", self.scale_cy)
        s_layout.addRow("Center Z:", self.scale_cz)
        self.tabs.addTab(scale_tab, "Scale")

        layout.addWidget(self.tabs)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        tab_idx = self.tabs.currentIndex()
        if tab_idx == 0:  # Translate
            return {
                'type': 'translate',
                'delta': [float(self.dx_input.text()),
                          float(self.dy_input.text()),
                          float(self.dz_input.text())],
            }
        elif tab_idx == 1:  # Rotate
            return {
                'type': 'rotate',
                'angle': float(self.rot_angle_input.text()),
                'axis': self.rot_axis_combo.currentText().lower(),
                'center': [float(self.rot_cx.text()),
                           float(self.rot_cy.text()),
                           float(self.rot_cz.text())],
            }
        elif tab_idx == 2:  # Mirror
            preset = self.mirror_plane_combo.currentIndex()
            if preset == 0:  # XY
                pp = [0, 0, 0]
                pn = [0, 0, 1]
            elif preset == 1:  # XZ
                pp = [0, 0, 0]
                pn = [0, 1, 0]
            elif preset == 2:  # YZ
                pp = [0, 0, 0]
                pn = [1, 0, 0]
            else:  # Custom
                pp = [float(self.mirror_px.text()),
                      float(self.mirror_py.text()),
                      float(self.mirror_pz.text())]
                pn = [float(self.mirror_nx.text()),
                      float(self.mirror_ny.text()),
                      float(self.mirror_nz.text())]
            return {
                'type': 'mirror',
                'plane_point': pp,
                'plane_normal': pn,
            }
        else:  # Scale
            return {
                'type': 'scale',
                'factor': float(self.scale_factor_input.text()),
                'center': [float(self.scale_cx.text()),
                           float(self.scale_cy.text()),
                           float(self.scale_cz.text())],
            }


class CopyGeometryDialog(QDialog):
    """Dialog for copying geometry with a translation offset."""

    def __init__(self, point_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Copy Geometry")
        self.setMinimumWidth(350)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(f"Copying {len(point_ids)} point(s) and dependent curves"))

        form = QFormLayout()
        self.dx_input = QLineEdit("0.0")
        self.dx_input.setValidator(QDoubleValidator())
        self.dy_input = QLineEdit("0.0")
        self.dy_input.setValidator(QDoubleValidator())
        self.dz_input = QLineEdit("0.0")
        self.dz_input.setValidator(QDoubleValidator())
        form.addRow("Offset X:", self.dx_input)
        form.addRow("Offset Y:", self.dy_input)
        form.addRow("Offset Z:", self.dz_input)
        layout.addLayout(form)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'delta': [float(self.dx_input.text()),
                      float(self.dy_input.text()),
                      float(self.dz_input.text())],
        }


# ---------------------------------------------------------------------------
# CAD-like operation dialogs
# ---------------------------------------------------------------------------

class SplitCurveDialog(QDialog):
    """Dialog for splitting a curve into N segments."""

    def __init__(self, curve_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Split Curve")
        self.setMinimumWidth(350)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.curve_combo = QComboBox()
        for cid in curve_ids:
            self.curve_combo.addItem(f"Curve {cid}", cid)
        form.addRow("Curve:", self.curve_combo)

        self.n_segments_input = QLineEdit("2")
        self.n_segments_input.setValidator(QIntValidator(2, 100))
        form.addRow("Number of Segments:", self.n_segments_input)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'curve_id': self.curve_combo.currentData(),
            'n_segments': int(self.n_segments_input.text()),
        }


class OffsetCurveDialog(QDialog):
    """Dialog for creating a parallel curve at a given distance."""

    def __init__(self, curve_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Offset Curve")
        self.setMinimumWidth(350)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.curve_combo = QComboBox()
        for cid in curve_ids:
            self.curve_combo.addItem(f"Curve {cid}", cid)
        form.addRow("Curve:", self.curve_combo)

        self.distance_input = QLineEdit("1.0")
        self.distance_input.setValidator(QDoubleValidator())
        form.addRow("Offset Distance:", self.distance_input)

        layout.addLayout(form)
        layout.addWidget(QLabel("(Positive = outward, Negative = inward for arcs/circles)"))

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'curve_id': self.curve_combo.currentData(),
            'distance': float(self.distance_input.text()),
        }


class ProjectPointToCurveDialog(QDialog):
    """Dialog for projecting a point onto a curve."""

    def __init__(self, curve_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Project Point to Curve")
        self.setMinimumWidth(400)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.curve_combo = QComboBox()
        for cid in curve_ids:
            self.curve_combo.addItem(f"Curve {cid}", cid)
        form.addRow("Target Curve:", self.curve_combo)

        self.sx_input = QLineEdit("0.0")
        self.sx_input.setValidator(QDoubleValidator())
        self.sy_input = QLineEdit("0.0")
        self.sy_input.setValidator(QDoubleValidator())
        self.sz_input = QLineEdit("0.0")
        self.sz_input.setValidator(QDoubleValidator())
        form.addRow("Source X:", self.sx_input)
        form.addRow("Source Y:", self.sy_input)
        form.addRow("Source Z:", self.sz_input)

        layout.addLayout(form)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'curve_id': self.curve_combo.currentData(),
            'source_xyz': [float(self.sx_input.text()),
                           float(self.sy_input.text()),
                           float(self.sz_input.text())],
        }


class FilletDialog(QDialog):
    """Dialog for creating a fillet arc between two lines."""

    def __init__(self, line_ids, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fillet")
        self.setMinimumWidth(350)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.line1_combo = QComboBox()
        self.line2_combo = QComboBox()
        for lid in line_ids:
            self.line1_combo.addItem(f"Line {lid}", lid)
            self.line2_combo.addItem(f"Line {lid}", lid)
        if len(line_ids) > 1:
            self.line2_combo.setCurrentIndex(1)
        form.addRow("Line 1:", self.line1_combo)
        form.addRow("Line 2:", self.line2_combo)

        self.radius_input = QLineEdit("1.0")
        self.radius_input.setValidator(QDoubleValidator(0.0001, 1e12, 6))
        form.addRow("Fillet Radius:", self.radius_input)

        layout.addLayout(form)
        layout.addWidget(QLabel("Both lines must share a common endpoint."))

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_parameters(self):
        return {
            'line_id_1': self.line1_combo.currentData(),
            'line_id_2': self.line2_combo.currentData(),
            'radius': float(self.radius_input.text()),
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_combo_by_data(combo, value):
    """Set a QComboBox's current index to the item matching the given data value."""
    for i in range(combo.count()):
        if combo.itemData(i) == value:
            combo.setCurrentIndex(i)
            return
