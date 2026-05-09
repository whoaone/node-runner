"""Cross-Section dialog: toggle, choose plane normal, slide position.

Non-modal so the user can scrub the section while interacting with the
3D view. Lives over the main window; closing it disables cross-section
(matching Femap's behavior of "section view is on while the dialog is up").
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, QSlider,
    QPushButton, QCheckBox, QDoubleSpinBox, QFrame,
)

from node_runner.cross_section import AXIS_PRESETS


class CrossSectionDialog(QDialog):
    """Section-plane controls: enable, normal direction, position slider.

    Position is stored as a 0..1 fraction of the model's extent along the
    chosen normal, so the slider stays meaningful when the model changes
    size or when the user picks a different normal.
    """

    SLIDER_RANGE = 1000  # finer slider resolution = smoother scrubbing

    def __init__(self, controller, main_window, parent=None):
        super().__init__(parent or main_window)
        self.setWindowTitle("Cross Section")
        self.setModal(False)
        self.setMinimumWidth(360)

        self.controller = controller
        self.main_window = main_window

        layout = QVBoxLayout(self)

        self.enable_check = QCheckBox("Enable cross section")
        self.enable_check.setChecked(controller.is_enabled)
        self.enable_check.toggled.connect(self._on_enable_toggled)
        layout.addWidget(self.enable_check)

        sep1 = QFrame(); sep1.setFrameShape(QFrame.HLine); sep1.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep1)

        layout.addWidget(QLabel("Plane normal:"))
        self.normal_combo = QComboBox()
        for label, _ in AXIS_PRESETS:
            self.normal_combo.addItem(label)
        # Initialize to controller's stored choice
        idx = self.normal_combo.findText(controller.normal_label)
        if idx >= 0:
            self.normal_combo.setCurrentIndex(idx)
        self.normal_combo.currentTextChanged.connect(self._on_normal_changed)
        layout.addWidget(self.normal_combo)

        layout.addWidget(QLabel("Plane position:"))
        slider_row = QHBoxLayout()
        self.position_slider = QSlider(Qt.Horizontal)
        self.position_slider.setRange(0, self.SLIDER_RANGE)
        self.position_slider.setValue(int(controller.slider_fraction * self.SLIDER_RANGE))
        self.position_slider.valueChanged.connect(self._on_slider_changed)
        slider_row.addWidget(self.position_slider, 1)

        self.position_value = QDoubleSpinBox()
        self.position_value.setDecimals(3)
        self.position_value.setRange(-1e9, 1e9)
        self.position_value.setSingleStep(0.1)
        self.position_value.setMinimumWidth(110)
        self.position_value.editingFinished.connect(self._on_value_edited)
        slider_row.addWidget(self.position_value)
        layout.addLayout(slider_row)

        info = QLabel(
            "Tip: With a model loaded, slide to scrub the section through "
            "the plane normal. Disable to restore the full view. The clip "
            "is GPU-side, so even huge models stay smooth."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(info)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        reset_btn = QPushButton("Center plane")
        reset_btn.clicked.connect(self._on_reset)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        button_row.addWidget(reset_btn)
        button_row.addWidget(close_btn)
        layout.addLayout(button_row)

        # Apply initial state to plane immediately so the slider matches reality
        self._apply_normal_to_controller()
        self._apply_slider_to_controller(self.position_slider.value())

    # ----- helpers -----

    def _model_bounds_along_normal(self):
        """Return (lo, hi) extent of the loaded model projected onto the
        current normal direction. Falls back to (-10, 10) if no model
        is loaded yet."""
        grid = getattr(self.main_window, "current_grid", None)
        if grid is None or grid.n_points == 0:
            return -10.0, 10.0
        b = grid.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
        nx, ny, nz = self.controller.plane.GetNormal()
        # Project the 8 corners of the bounding box onto the normal axis,
        # then take min and max.
        corners = [
            (b[0], b[2], b[4]), (b[1], b[2], b[4]),
            (b[0], b[3], b[4]), (b[1], b[3], b[4]),
            (b[0], b[2], b[5]), (b[1], b[2], b[5]),
            (b[0], b[3], b[5]), (b[1], b[3], b[5]),
        ]
        projs = [cx * nx + cy * ny + cz * nz for cx, cy, cz in corners]
        return min(projs), max(projs)

    def _apply_normal_to_controller(self):
        label = self.normal_combo.currentText()
        for lab, vec in AXIS_PRESETS:
            if lab == label:
                self.controller.set_normal(vec[0], vec[1], vec[2], label=lab)
                break

    def _apply_slider_to_controller(self, slider_value):
        frac = slider_value / float(self.SLIDER_RANGE)
        self.controller.set_slider_fraction(frac)
        lo, hi = self._model_bounds_along_normal()
        # Translate fraction into a world-space distance along the normal
        # and place the plane origin at that point.
        dist = lo + frac * (hi - lo)
        nx, ny, nz = self.controller.plane.GetNormal()
        self.controller.set_origin(nx * dist, ny * dist, nz * dist)
        self.position_value.blockSignals(True)
        self.position_value.setValue(dist)
        self.position_value.blockSignals(False)

    # ----- handlers -----

    def _on_enable_toggled(self, on):
        if on:
            self.controller.enable()
        else:
            self.controller.disable()

    def _on_normal_changed(self, _label):
        self._apply_normal_to_controller()
        # Re-evaluate slider against new bounds so the section visually
        # reappears at the same fraction along the new normal.
        self._apply_slider_to_controller(self.position_slider.value())

    def _on_slider_changed(self, value):
        self._apply_slider_to_controller(value)

    def _on_value_edited(self):
        # User typed a value into the spin box - convert back to slider fraction
        lo, hi = self._model_bounds_along_normal()
        if hi > lo:
            frac = (self.position_value.value() - lo) / (hi - lo)
            frac = max(0.0, min(1.0, frac))
            self.position_slider.blockSignals(True)
            self.position_slider.setValue(int(frac * self.SLIDER_RANGE))
            self.position_slider.blockSignals(False)
        self._apply_slider_to_controller(self.position_slider.value())

    def _on_reset(self):
        self.position_slider.setValue(self.SLIDER_RANGE // 2)

    # ----- lifecycle -----

    def closeEvent(self, ev):
        # Closing the dialog disables the section (Femap behavior). The
        # controller stays around so reopening picks up the same normal
        # and fraction.
        self.controller.disable()
        if self.enable_check.isChecked():
            self.enable_check.blockSignals(True)
            self.enable_check.setChecked(False)
            self.enable_check.blockSignals(False)
        super().closeEvent(ev)
