"""Export options dialog for choosing Nastran field-width on save."""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QRadioButton, QButtonGroup,
    QLabel, QPushButton, QGroupBox, QCheckBox,
)
from PySide6.QtCore import Qt


FORMAT_LABELS = {
    "short": "Short Field (8-character)",
    "long":  "Long Field (16-character)",
    "free":  "Free Field (comma-separated)",
}

FORMAT_DESCRIPTIONS = {
    "short": "Industry standard for portability and readability.",
    "long":  "Higher precision for coordinates and large IDs.",
    "free":  "Comma-separated; for modern solvers and rapid manual edits.",
}


class ExportOptionsDialog(QDialog):
    """Modal dialog asking the user which Nastran field-width to use.

    The dialog can be primed with a `default_format` (one of 'short', 'long',
    'free') so the user's last choice or QSettings preference is preselected.
    """

    def __init__(self, default_format="short", parent=None,
                 allow_save_as_default=True):
        super().__init__(parent)
        self.setWindowTitle("Export Options")
        self.setModal(True)
        self.setMinimumWidth(380)

        self._selected_format = default_format if default_format in FORMAT_LABELS else "short"
        self._save_as_default = False

        layout = QVBoxLayout(self)

        header = QLabel("Choose the Nastran field format for this export:")
        header.setWordWrap(True)
        layout.addWidget(header)

        group_box = QGroupBox()
        group_layout = QVBoxLayout(group_box)
        self._button_group = QButtonGroup(self)
        self._radios = {}

        for key, label in FORMAT_LABELS.items():
            radio = QRadioButton(label)
            radio.setChecked(key == self._selected_format)
            radio.toggled.connect(lambda checked, k=key: self._on_radio_toggled(k, checked))
            desc = QLabel(FORMAT_DESCRIPTIONS[key])
            desc.setStyleSheet("color: #888; font-size: 11px; margin-left: 22px; margin-bottom: 4px;")
            desc.setWordWrap(True)
            group_layout.addWidget(radio)
            group_layout.addWidget(desc)
            self._button_group.addButton(radio)
            self._radios[key] = radio

        layout.addWidget(group_box)

        if allow_save_as_default:
            self._default_check = QCheckBox("Remember this choice for future exports")
            layout.addWidget(self._default_check)
        else:
            self._default_check = None

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        ok = QPushButton("Export")
        ok.setDefault(True)
        ok.clicked.connect(self.accept)
        button_row.addWidget(cancel)
        button_row.addWidget(ok)
        layout.addLayout(button_row)

    def _on_radio_toggled(self, key, checked):
        if checked:
            self._selected_format = key

    @property
    def selected_format(self):
        return self._selected_format

    @property
    def save_as_default(self):
        if self._default_check is None:
            return False
        return self._default_check.isChecked()


class ExportDefaultsDialog(QDialog):
    """Settings-menu dialog that just sets the default export format.

    Differs from ExportOptionsDialog in intent only: this is reached from
    Settings, not Save, and it always writes to QSettings on Accept.
    """

    def __init__(self, current_format="short", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export Defaults")
        self.setModal(True)
        self.setMinimumWidth(380)

        self._selected_format = current_format if current_format in FORMAT_LABELS else "short"

        layout = QVBoxLayout(self)

        header = QLabel("Default Nastran field format used when saving / exporting:")
        header.setWordWrap(True)
        layout.addWidget(header)

        group_box = QGroupBox()
        group_layout = QVBoxLayout(group_box)
        self._button_group = QButtonGroup(self)
        for key, label in FORMAT_LABELS.items():
            radio = QRadioButton(label)
            radio.setChecked(key == self._selected_format)
            radio.toggled.connect(lambda checked, k=key: self._on_radio_toggled(k, checked))
            desc = QLabel(FORMAT_DESCRIPTIONS[key])
            desc.setStyleSheet("color: #888; font-size: 11px; margin-left: 22px; margin-bottom: 4px;")
            desc.setWordWrap(True)
            group_layout.addWidget(radio)
            group_layout.addWidget(desc)
            self._button_group.addButton(radio)
        layout.addWidget(group_box)

        info = QLabel(
            "The save dialog will offer a per-file override. "
            "If a file's format is auto-detected on import, that detection "
            "wins for round-tripping unless you change the format manually."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #666; font-size: 10px;")
        layout.addWidget(info)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        ok = QPushButton("Save")
        ok.setDefault(True)
        ok.clicked.connect(self.accept)
        button_row.addWidget(cancel)
        button_row.addWidget(ok)
        layout.addLayout(button_row)

    def _on_radio_toggled(self, key, checked):
        if checked:
            self._selected_format = key

    @property
    def selected_format(self):
        return self._selected_format


class UnitsDialog(QDialog):
    """Settings-menu dialog for choosing display units (SI / English).

    Cosmetic for now: drives the status-bar Units label only. No engineering
    unit conversion is performed on the model.
    """

    def __init__(self, current_units="SI", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Units")
        self.setModal(True)
        self.setMinimumWidth(320)

        self._selected_units = current_units if current_units in ("SI", "English") else "SI"

        layout = QVBoxLayout(self)
        header = QLabel("Display units shown in the status bar:")
        layout.addWidget(header)

        self._button_group = QButtonGroup(self)
        for opt in ("SI", "English"):
            radio = QRadioButton(opt)
            radio.setChecked(opt == self._selected_units)
            radio.toggled.connect(lambda checked, o=opt: self._on_radio_toggled(o, checked))
            self._button_group.addButton(radio)
            layout.addWidget(radio)

        info = QLabel(
            "This setting is cosmetic for now: it labels the status bar but "
            "does not convert model values."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #666; font-size: 10px;")
        layout.addWidget(info)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        ok = QPushButton("Save")
        ok.setDefault(True)
        ok.clicked.connect(self.accept)
        button_row.addWidget(cancel)
        button_row.addWidget(ok)
        layout.addLayout(button_row)

    def _on_radio_toggled(self, opt, checked):
        if checked:
            self._selected_units = opt

    @property
    def selected_units(self):
        return self._selected_units
