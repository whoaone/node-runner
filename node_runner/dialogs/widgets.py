from PySide6 import QtCore
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QRadioButton, QLineEdit, QComboBox, QLabel, QPushButton,
)


class OrientationWidget(QWidget):
    """A reusable widget for defining element orientation."""
    # Add a signal to notify the parent dialog when the pick button is clicked
    pick_orientation_node_requested = QtCore.Signal()

    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model = model

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        orient_group = QGroupBox("Orientation")
        orient_layout = QVBoxLayout(orient_group)

        self.orient_default_rb = QRadioButton("Align with Global Coordinate System", checked=True)
        self.orient_vector_rb = QRadioButton("Custom Vector (defines XY plane)")
        self.orient_node_rb = QRadioButton("Orientation Node (defines XY plane)")
        self.orient_cid_rb = QRadioButton("By Coordinate System (CID)")

        self.orient_vector_widget = QWidget()
        vec_form = QFormLayout(self.orient_vector_widget); vec_form.setContentsMargins(20, 5, 5, 5)
        self.vec_x, self.vec_y, self.vec_z = QLineEdit("1.0"), QLineEdit("0.0"), QLineEdit("0.0")
        vec_layout = QHBoxLayout(); vec_layout.addWidget(self.vec_x); vec_layout.addWidget(self.vec_y); vec_layout.addWidget(self.vec_z)
        vec_form.addRow("Vector (X,Y,Z):", vec_layout)

        self.orient_node_widget = QWidget()
        node_layout = QHBoxLayout(self.orient_node_widget); node_layout.setContentsMargins(20, 0, 0, 0)
        self.orient_node_id = QLineEdit()
        pick_orient_node_btn = QPushButton("Pick...")
        node_layout.addWidget(QLabel("Node ID (G0):")); node_layout.addWidget(self.orient_node_id); node_layout.addWidget(pick_orient_node_btn)

        self.orient_cid_widget = QWidget()
        cid_layout = QHBoxLayout(self.orient_cid_widget); cid_layout.setContentsMargins(20,0,0,0)
        self.cid_combo = QComboBox()
        self.cid_combo.addItems([str(cid) for cid in model.coords.keys()])
        cid_layout.addWidget(QLabel("Coord ID:")); cid_layout.addWidget(self.cid_combo, 1)

        orient_layout.addWidget(self.orient_default_rb)
        orient_layout.addWidget(self.orient_vector_rb); orient_layout.addWidget(self.orient_vector_widget)
        orient_layout.addWidget(self.orient_node_rb); orient_layout.addWidget(self.orient_node_widget)
        orient_layout.addWidget(self.orient_cid_rb); orient_layout.addWidget(self.orient_cid_widget)
        layout.addWidget(orient_group)

        pick_orient_node_btn.clicked.connect(self.pick_orientation_node_requested.emit)
        self.orient_vector_rb.toggled.connect(self._update_ui_visibility)
        self.orient_node_rb.toggled.connect(self._update_ui_visibility)
        self.orient_cid_rb.toggled.connect(self._update_ui_visibility)
        self._update_ui_visibility()

    def _update_ui_visibility(self):
        self.orient_vector_widget.setVisible(self.orient_vector_rb.isChecked())
        self.orient_node_widget.setVisible(self.orient_node_rb.isChecked())
        self.orient_cid_widget.setVisible(self.orient_cid_rb.isChecked())

    def get_orientation(self):
        """Returns the orientation dictionary based on the UI state."""
        orientation = {'method': 'default'}
        if self.orient_vector_rb.isChecked():
            vec = [float(self.vec_x.text()), float(self.vec_y.text()), float(self.vec_z.text())]
            orientation = {'method': 'vector', 'values': vec}
        elif self.orient_node_rb.isChecked():
            orientation = {'method': 'node', 'values': [int(self.orient_node_id.text())]}
        elif self.orient_cid_rb.isChecked():
            orientation = {'method': 'cid', 'values': [int(self.cid_combo.currentText())]}
        return orientation

    def set_orientation(self, element):
        """Sets the UI state based on a pyNastran element object."""
        # Check if the element CAN have a CID and if it's set
        if hasattr(element, 'cid') and element.cid is not None:
            self.orient_cid_rb.setChecked(True)
            self.cid_combo.setCurrentText(str(element.cid))
        # Then check for G0
        elif hasattr(element, 'g0') and element.g0 is not None:
            self.orient_node_rb.setChecked(True)
            self.orient_node_id.setText(str(element.g0))
        # Then check for the vector
        elif hasattr(element, 'x') and element.x is not None and not all(v is None for v in element.x):
            self.orient_vector_rb.setChecked(True)
            self.vec_x.setText(str(element.x[0]))
            self.vec_y.setText(str(element.x[1]))
            self.vec_z.setText(str(element.x[2]))
        # If none of the above, it's a default orientation
        else:
            self.orient_default_rb.setChecked(True)
