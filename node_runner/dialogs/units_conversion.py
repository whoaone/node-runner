"""Tools > Convert Units... dialog.

Femap-style: the model is unitless, but you can ask the tool to rescale
every relevant card by a set of factors (length, force, mass).
Stress, density, and moment factors are derived from those three:

    stress  = force / length**2
    density = mass / length**3
    moment  = force * length

What gets scaled when you hit Apply:
  GRID xyz                   x length
  PSHELL t                   x length
  PBEAM/PBAR area            x length**2
  PBEAM/PBAR I1, I2, I12, J  x length**4
  MAT1 E, G                  x stress
  MAT1 rho                   x density
  FORCE magnitudes           x force
  MOMENT magnitudes          x moment
  PLOAD2/4 pressures         x stress
  PLOAD1 forces (FX/FY/FZ)   x force; moments (MX/MY/MZ) x moment
  GRAV magnitude             x length / time**2 (skipped here; gravity is
                              usually a unit-system constant the user
                              re-types after a conversion)
  SPCD enforced values       x length

Anything not recognized is left alone. Coord systems' offsets are
scaled by length.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
    QDoubleSpinBox, QGroupBox, QComboBox, QDialogButtonBox, QCheckBox,
    QFormLayout,
)


# Factor presets keyed by description. Each entry is
# (length_factor, force_factor, mass_factor).
PRESETS = {
    "(custom)":                          (1.0,    1.0,    1.0),
    "m -> mm  (length only)":            (1000.0, 1.0,    1.0),
    "mm -> m":                           (0.001,  1.0,    1.0),
    "in -> mm":                          (25.4,   1.0,    1.0),
    "mm -> in":                          (1.0/25.4, 1.0,  1.0),
    "m,kg,N  -> mm,t,N":                 (1000.0, 1.0,    0.001),
    "in,lbf,lbf-s2/in -> mm,N,t":        (25.4,   4.448222, 175.1268),
    "ft -> m":                           (0.3048, 1.0,    1.0),
    "ft -> in":                          (12.0,   1.0,    1.0),
}


class UnitConversionDialog(QDialog):
    """Three editable factors (length / force / mass) plus a preset combo."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Convert Units")
        self.setMinimumWidth(440)

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(
            "Scale the entire model by independent length, force, and mass "
            "factors. Stress, density, and moment factors are derived. "
            "Use a preset or enter custom factors below."
        ))

        # Preset combo
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Preset:"))
        self._preset_combo = QComboBox()
        for label in PRESETS:
            self._preset_combo.addItem(label)
        self._preset_combo.currentTextChanged.connect(self._on_preset)
        preset_row.addWidget(self._preset_combo, 1)
        layout.addLayout(preset_row)

        # Factor inputs
        factors_box = QGroupBox("Factors")
        form = QFormLayout(factors_box)
        self._length_spin = QDoubleSpinBox()
        self._length_spin.setRange(1e-12, 1e12); self._length_spin.setDecimals(8)
        self._length_spin.setValue(1.0)
        self._length_spin.valueChanged.connect(self._on_factor_changed)
        form.addRow("Length factor:", self._length_spin)

        self._force_spin = QDoubleSpinBox()
        self._force_spin.setRange(1e-12, 1e12); self._force_spin.setDecimals(8)
        self._force_spin.setValue(1.0)
        self._force_spin.valueChanged.connect(self._on_factor_changed)
        form.addRow("Force factor:", self._force_spin)

        self._mass_spin = QDoubleSpinBox()
        self._mass_spin.setRange(1e-12, 1e12); self._mass_spin.setDecimals(8)
        self._mass_spin.setValue(1.0)
        self._mass_spin.valueChanged.connect(self._on_factor_changed)
        form.addRow("Mass factor:", self._mass_spin)

        layout.addWidget(factors_box)

        # Derived factors readout
        self._derived_label = QLabel()
        self._derived_label.setStyleSheet("color: #888; font-family: monospace;")
        layout.addWidget(self._derived_label)
        self._refresh_derived()

        # Optional: also rewrite the unit-system hint after conversion
        self._update_label_check = QCheckBox(
            "Also update the unit-label hint in the status bar"
        )
        self._update_label_check.setChecked(False)
        layout.addWidget(self._update_label_check)

        info = QLabel(
            "Apply runs as a single undoable command. Original values are "
            "snapshotted so you can hit Edit > Undo if the conversion was "
            "wrong."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #666; font-size: 10px;")
        layout.addWidget(info)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Apply")
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def _on_preset(self, label):
        if label == "(custom)" or label not in PRESETS:
            return
        L, F, M = PRESETS[label]
        for w, v in ((self._length_spin, L), (self._force_spin, F), (self._mass_spin, M)):
            w.blockSignals(True); w.setValue(v); w.blockSignals(False)
        self._refresh_derived()

    def _on_factor_changed(self, _v=None):
        # Any manual edit reverts the preset to (custom).
        if self._preset_combo.currentText() != "(custom)":
            self._preset_combo.blockSignals(True)
            self._preset_combo.setCurrentText("(custom)")
            self._preset_combo.blockSignals(False)
        self._refresh_derived()

    def _refresh_derived(self):
        L = self._length_spin.value()
        F = self._force_spin.value()
        M = self._mass_spin.value()
        if L <= 0 or F <= 0 or M <= 0:
            self._derived_label.setText("Derived: <invalid>")
            return
        stress = F / (L * L)
        density = M / (L ** 3)
        moment = F * L
        area = L * L
        moi = L ** 4
        self._derived_label.setText(
            f"Derived:  stress * {stress:.6g}   density * {density:.6g}   "
            f"moment * {moment:.6g}\n"
            f"          area * {area:.6g}        moi/J * {moi:.6g}"
        )

    @property
    def factors(self):
        return {
            'length': float(self._length_spin.value()),
            'force':  float(self._force_spin.value()),
            'mass':   float(self._mass_spin.value()),
        }

    @property
    def update_label_hint(self):
        return self._update_label_check.isChecked()
