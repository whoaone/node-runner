"""Export options dialog for choosing Nastran field-width on save."""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QRadioButton, QButtonGroup,
    QLabel, QPushButton, QGroupBox, QCheckBox, QComboBox, QLineEdit,
    QFileDialog, QFormLayout, QFrame, QDialogButtonBox,
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

        # v5.0.0 item 6: surface the current default prominently inside the
        # dialog now that the status-bar "Export: ..." label is gone.
        current_label_text = FORMAT_LABELS.get(self._selected_format, self._selected_format)
        self._current_label = QLabel(
            f"<b>Current default:</b> {current_label_text}"
        )
        self._current_label.setStyleSheet(
            "color: #a6e3a1; font-size: 12px; padding: 4px 0;"
        )
        layout.addWidget(self._current_label)

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
            if hasattr(self, '_current_label') and self._current_label is not None:
                self._current_label.setText(
                    f"<b>Current default:</b> {FORMAT_LABELS.get(key, key)}"
                )

    @property
    def selected_format(self):
        return self._selected_format


class UnitsDialog(QDialog):
    """Settings-menu dialog for setting a free-text unit-system hint.

    Node Runner is intentionally unitless (Femap-style): the model is just
    numbers. This dialog only sets a label that appears in the status bar
    so you remember which unit system you're working in. To actually
    rescale the model values, use **Tools > Convert Units...**.
    """

    SUGGESTIONS = ("", "mm / N / t", "mm / N / kg", "m / N / kg",
                   "in / lbf / slug", "in / lbf / lbf-s2/in")

    def __init__(self, current_units="", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Unit Label")
        self.setModal(True)
        self.setMinimumWidth(380)

        self._selected_units = str(current_units or "")

        layout = QVBoxLayout(self)
        header = QLabel(
            "Unit-system hint shown in the status bar. Free text - no "
            "validation, no enforcement, no automatic conversion. Leave "
            "blank to hide the field."
        )
        header.setWordWrap(True)
        layout.addWidget(header)

        from PySide6.QtWidgets import QLineEdit, QComboBox
        self._line = QLineEdit(self._selected_units)
        self._line.setPlaceholderText("e.g. mm / N / t  or  in / lbf / slug")
        layout.addWidget(self._line)

        # Quick-pick combo seeded with common conventions
        layout.addWidget(QLabel("Quick pick:"))
        self._suggest_combo = QComboBox()
        self._suggest_combo.addItems(["(custom)"] + list(self.SUGGESTIONS[1:]))
        self._suggest_combo.currentTextChanged.connect(self._on_suggestion)
        layout.addWidget(self._suggest_combo)

        info = QLabel(
            "To actually rescale node coordinates, properties, materials, "
            "and loads, use Tools > Convert Units..."
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

    def _on_suggestion(self, text):
        if text and text != "(custom)":
            self._line.setText(text)

    @property
    def selected_units(self):
        return self._line.text().strip()


# ---------------------------------------------------------------------------
# v5.1.0 item 28: SaveBdfDialog
# ---------------------------------------------------------------------------

SOLVER_TARGETS = {
    "generic": "Generic (MSC / NX compatible)",
    "MSC":     "MSC Nastran",
    "NX":      "Siemens NX Nastran",
    "MYSTRAN": "MYSTRAN (open-source)",
}


class SaveBdfDialog(QDialog):
    """v5.1.0 item 28: the Save BDF dialog.

    Replaces the v5.0.x flow of "QFileDialog -> ExportOptionsDialog
    (field-width only)" with a single richer dialog that lets the user
    pick:

    - Output path
    - Solver target (Generic / MSC / NX / MYSTRAN) -- drives the
      ``_write_bdf(target=...)`` translator path
    - Field format (Short / Long / Free)
    - Scope (Full Model / Active AnalysisSet / Group: <name>)
      [the AnalysisSet entry is shown only when one is active; group
      entries are populated from MainWindow.groups]
    - Include $ NR-META block (default ON)
    - Write to sidecar .nrmeta instead of inline (default OFF)
    - Remember choices

    The dialog is presentation-only: it returns its selections via
    :attr:`result_payload`; the caller (MainWindow._save_file_dialog)
    does the actual write.
    """

    def __init__(self, *, default_path="",
                 default_format="short",
                 default_target="generic",
                 active_analysis_set_name=None,
                 group_names=None,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Save BDF — Export Options")
        self.setModal(True)
        self.setMinimumWidth(540)
        group_names = list(group_names or [])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 8)
        outer.setSpacing(8)

        intro = QLabel(
            "Choose where to save and how to translate the deck. "
            "Source models on disk are never modified -- the deck is "
            "written to the path below.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        outer.addWidget(intro)

        form = QFormLayout()
        form.setSpacing(6)

        # ---- Output path ----
        path_row = QHBoxLayout()
        self._path_edit = QLineEdit(default_path)
        self._path_edit.setPlaceholderText("C:\\path\\to\\model.bdf")
        path_row.addWidget(self._path_edit, 1)
        browse_btn = QPushButton("Browse…")
        browse_btn.clicked.connect(self._pick_path)
        path_row.addWidget(browse_btn)
        form.addRow(QLabel("Output path:"), path_row)

        # ---- Solver target ----
        self._target_combo = QComboBox()
        for key, label in SOLVER_TARGETS.items():
            self._target_combo.addItem(label, key)
        # Select the default target.
        for i in range(self._target_combo.count()):
            if self._target_combo.itemData(i) == default_target:
                self._target_combo.setCurrentIndex(i)
                break
        self._target_combo.setToolTip(
            "<b>Solver target</b><br>"
            "<i>Generic</i>: vanilla BDF, MSC/NX-compatible. No "
            "solver-specific PARAMs injected.<br>"
            "<i>MYSTRAN</i>: applies the v5.0.0 MYSTRAN translation "
            "(drops MSC-only PARAMs; injects SOLLIB / QUAD4TYP).<br>"
            "<i>MSC / NX</i>: same as Generic in v5.1.0 -- placeholder "
            "for future per-vendor translators.")
        form.addRow(QLabel("Solver target:"), self._target_combo)

        # ---- Field format ----
        self._format_combo = QComboBox()
        for key, label in FORMAT_LABELS.items():
            self._format_combo.addItem(label, key)
        for i in range(self._format_combo.count()):
            if self._format_combo.itemData(i) == default_format:
                self._format_combo.setCurrentIndex(i)
                break
        form.addRow(QLabel("Field format:"), self._format_combo)

        # ---- Scope ----
        self._scope_combo = QComboBox()
        self._scope_combo.addItem("Full model", ("full", None))
        if active_analysis_set_name:
            self._scope_combo.addItem(
                f"Active AnalysisSet: {active_analysis_set_name}",
                ("analysis_set", None))
        for gname in group_names:
            self._scope_combo.addItem(f"Group: {gname}", ("group", gname))
        self._scope_combo.setToolTip(
            "<b>Scope</b><br>"
            "Choose what to write. <i>Full model</i> writes everything.<br>"
            "<i>Active AnalysisSet</i> writes the scoped sub-deck "
            "defined by the active set (group_target + load/SPC SIDs). "
            "<i>Group: NAME</i> writes only entities in that group plus "
            "auto-collected dependencies.")
        form.addRow(QLabel("Scope:"), self._scope_combo)

        outer.addLayout(form)

        # ---- NR-META options ----
        meta_box = QGroupBox("Node Runner metadata")
        meta_lay = QVBoxLayout(meta_box)
        self._include_nr_meta = QCheckBox(
            "Include Node Runner metadata (groups + tree state)")
        self._include_nr_meta.setChecked(True)
        self._include_nr_meta.setToolTip(
            "Embeds group definitions and hidden-state as ordinary "
            "Nastran $ NR-META comments. MSC / NX / MYSTRAN ignore "
            "them. Node Runner restores them on re-import.")
        meta_lay.addWidget(self._include_nr_meta)
        self._sidecar_nrmeta = QCheckBox(
            "Write to sidecar .nrmeta file instead of inline")
        self._sidecar_nrmeta.setToolTip(
            "When ON, the metadata is written to <bdf>.nrmeta and the "
            ".bdf stays free of Node Runner comments. Useful when "
            "handing the deck to a coworker who doesn't want $ NR-META "
            "lines in their bulk data.")
        meta_lay.addWidget(self._sidecar_nrmeta)
        # When include is unchecked, sidecar should grey out.
        self._include_nr_meta.toggled.connect(
            lambda on: self._sidecar_nrmeta.setEnabled(on))
        outer.addWidget(meta_box)

        # ---- Remember choices ----
        self._remember = QCheckBox(
            "Remember solver target + field format for future saves")
        outer.addWidget(self._remember)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        outer.addWidget(sep)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        save = QPushButton("Save")
        save.setDefault(True)
        save.clicked.connect(self.accept)
        button_row.addWidget(cancel)
        button_row.addWidget(save)
        outer.addLayout(button_row)

    def _pick_path(self):
        start = self._path_edit.text() or ""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save BDF", start,
            "Nastran Bulk Data (*.bdf *.dat *.nas);;All files (*.*)")
        if path:
            self._path_edit.setText(path)

    @property
    def result_payload(self) -> dict:
        scope_kind, scope_arg = self._scope_combo.currentData() or (
            "full", None)
        return {
            "path": self._path_edit.text().strip(),
            "target": self._target_combo.currentData() or "generic",
            "field_format": self._format_combo.currentData() or "short",
            "scope_kind": scope_kind,           # 'full' | 'analysis_set' | 'group'
            "scope_arg": scope_arg,             # group name when scope_kind=='group'
            "include_nr_meta": bool(self._include_nr_meta.isChecked()),
            "sidecar_nrmeta": bool(
                self._sidecar_nrmeta.isChecked()
                and self._include_nr_meta.isChecked()),
            "remember": bool(self._remember.isChecked()),
        }
