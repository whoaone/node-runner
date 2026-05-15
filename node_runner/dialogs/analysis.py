"""Analysis Set Manager - AnalysisSet dataclass, manager dialog, and output requests dialog."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from PySide6 import QtCore
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QComboBox, QSpinBox,
    QDialogButtonBox, QTabWidget, QWidget, QListWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QMessageBox,
    QCheckBox,
)
from PySide6.QtGui import QDoubleValidator, QIntValidator


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class AnalysisSet:
    id: int
    name: str
    sol_type: int | None = None  # None, 101, 103, 105, 106

    # Output requests - global defaults for all subcases
    output_requests: dict = field(default_factory=lambda: {
        'displacement': 'NONE',
        'velocity': 'NONE',
        'acceleration': 'NONE',
        'stress': 'NONE',
        'strain': 'NONE',
        'force': 'NONE',
        'spcforce': 'NONE',
        'mpcforce': 'NONE',
        'gpforce': 'NONE',
        'oload': 'NONE',
    })

    # Eigenvalue params (for SOL 103/105)
    eigrl: dict | None = None  # {'sid': int, 'v1': float|None, 'v2': float|None, 'nd': int|None}

    # PARAM cards
    params: dict = field(default_factory=lambda: {
        'POST': -1,
        'MAXRATIO': 1e7,
    })

    # Subcases
    subcases: list = field(default_factory=list)
    # Each: {'id', 'load_sid', 'spc_sid', 'method_sid', 'statsub_id',
    #         'disp', 'stress', 'force', 'strain'}

    # ----- v5.1.0 item 26: scope + solver target -----
    #
    # solver_target drives which translator path _write_bdf takes when
    # this AnalysisSet is the export scope. 'generic' = vanilla
    # MSC/NX-compatible BDF (no translator); 'MYSTRAN' applies the
    # v5.0.0 mystran_export translator. 'MSC' / 'NX' are placeholders
    # for future vendor-specific translators -- v5.1.0 treats them as
    # 'generic'.
    solver_target: str = 'MYSTRAN'      # 'MYSTRAN' | 'MSC' | 'NX' | 'generic'

    # group_target: None = full model; otherwise the name of a group
    # in MainWindow.groups whose nodes/elements/properties/materials/
    # coords define the scoped sub-deck. Run + export operations on
    # this AnalysisSet use scope.py to build the sub-model.
    group_target: str | None = None

    # Explicit white-lists of LOAD / SPC SIDs that this AnalysisSet
    # exports. Empty lists mean "include every SID present in
    # model.loads / model.spcs"; non-empty lists drop SIDs not in the
    # list. Subcases still reference SIDs by integer; if a subcase
    # references a SID not in these lists, scope.py drops the subcase
    # rather than emit a broken case-control reference.
    load_sids: list = field(default_factory=list)
    spc_sids: list = field(default_factory=list)

    # Per-AnalysisSet default field format. Drives the export pipeline
    # when the user picks "Run / Export with this AnalysisSet" from
    # the right-click context menu without going through the Save BDF
    # dialog. The Save BDF dialog still lets the user override.
    field_format: str = 'short'         # 'short' | 'long' | 'free'


# Preset output request configurations
OUTPUT_PRESETS = {
    'Minimal': {
        'displacement': 'ALL', 'velocity': 'NONE', 'acceleration': 'NONE',
        'stress': 'NONE', 'strain': 'NONE', 'force': 'NONE',
        'spcforce': 'NONE', 'mpcforce': 'NONE', 'gpforce': 'NONE', 'oload': 'NONE',
    },
    'Standard': {
        'displacement': 'ALL', 'velocity': 'NONE', 'acceleration': 'NONE',
        'stress': 'ALL', 'strain': 'NONE', 'force': 'ALL',
        'spcforce': 'ALL', 'mpcforce': 'NONE', 'gpforce': 'NONE', 'oload': 'NONE',
    },
    'Full': {k: 'ALL' for k in [
        'displacement', 'velocity', 'acceleration', 'stress', 'strain',
        'force', 'spcforce', 'mpcforce', 'gpforce', 'oload',
    ]},
}

SOL_TYPES = [
    (None, "None (Punch)"),
    (101, "SOL 101 - Linear Statics"),
    (103, "SOL 103 - Normal Modes"),
    (105, "SOL 105 - Buckling"),
    (106, "SOL 106 - Nonlinear Static"),
]

# Common Nastran PARAM cards with widget types and defaults
# (param_name, label, widget_type, default_value, options_or_range)
COMMON_PARAMS = [
    ('POST',     'POST (OP2 Output)',        'combo',  -1,   [(-1, '-1 (Default)'), (0, '0 (No OP2)'), (1, '1 (OP2)')]),
    ('AUTOSPC',  'AUTOSPC',                  'yesno',  'YES', None),
    ('BAILOUT',  'BAILOUT',                  'combo',  0,    [(0, '0 (Stop on singularity)'), (-1, '-1 (Continue)')]),
    ('GRDPNT',   'GRDPNT (Mass Summary)',    'int',    0,    (0, 999999999)),
    ('MAXRATIO', 'MAXRATIO',                 'float',  1e7,  (1.0, 1e15)),
    ('KROT',     'KROT (Drill Stiffness)',   'combo',  0,    [(0, '0 (Off)'), (100, '100 (Small)'), (1000, '1000 (Moderate)')]),
    ('PRGPST',   'PRGPST (Singularity Tbl)', 'yesno', 'NO',  None),
    ('SRCOMPS',  'SRCOMPS (Stress Recovery)', 'yesno', 'NO',  None),
    ('NOCOMPS',  'NOCOMPS (Composite Out)',  'combo',  1,    [(-1, '-1 (Ply+Total)'), (0, '0 (Ply only)'), (1, '1 (Total only)')]),
    ('WTMASS',   'WTMASS (Mass Units)',      'float',  1.0,  (1e-15, 1e15)),
]


# ---------------------------------------------------------------------------
# Analysis Set Manager Dialog
# ---------------------------------------------------------------------------

class AnalysisSetManagerDialog(QDialog):
    """Two-panel dialog for managing analysis sets."""

    def __init__(self, analysis_sets: dict[int, AnalysisSet],
                 active_set_id: int | None,
                 model=None, groups=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Analysis Set Manager")
        self.setMinimumSize(900, 600)

        self._analysis_sets = {k: self._copy_set(v) for k, v in analysis_sets.items()}
        self._active_set_id = active_set_id
        self._model = model
        # v5.1.0 item 26b: groups dict from MainWindow (name -> data
        # with nodes/elements/properties/materials/coords). Drives the
        # group_target combo on the new Scope tab. None when called
        # from contexts that don't have a model loaded.
        self._groups = groups or {}

        main_layout = QHBoxLayout(self)

        # --- Left panel: set list ---
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)

        left_layout.addWidget(QLabel("Analysis Sets:"))
        self.set_list = QListWidget()
        self.set_list.currentRowChanged.connect(self._on_set_selected)
        left_layout.addWidget(self.set_list)

        btn_layout = QHBoxLayout()
        new_btn = QPushButton("New")
        dup_btn = QPushButton("Duplicate")
        del_btn = QPushButton("Delete")
        active_btn = QPushButton("Set Active")
        new_btn.clicked.connect(self._new_set)
        dup_btn.clicked.connect(self._duplicate_set)
        del_btn.clicked.connect(self._delete_set)
        active_btn.clicked.connect(self._set_active)
        btn_layout.addWidget(new_btn)
        btn_layout.addWidget(dup_btn)
        btn_layout.addWidget(del_btn)
        btn_layout.addWidget(active_btn)
        left_layout.addLayout(btn_layout)
        main_layout.addWidget(left_widget, 1)

        # --- Right panel: tabs ---
        self.tabs = QTabWidget()
        self._build_general_tab()
        self._build_scope_tab()        # v5.1.0 item 26b
        self._build_output_tab()
        self._build_eigrl_tab()
        self._build_params_tab()
        self._build_subcases_tab()
        main_layout.addWidget(self.tabs, 3)

        # --- Bottom buttons ---
        right_vbox = QVBoxLayout()
        right_vbox.addWidget(self.tabs)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self._on_accept)
        button_box.rejected.connect(self.reject)
        right_vbox.addWidget(button_box)

        main_layout.addLayout(right_vbox, 3)
        # Remove tabs from the left side since we just re-added it with button box
        main_layout.removeWidget(self.tabs)

        self._populate_set_list()
        if self.set_list.count() > 0:
            self.set_list.setCurrentRow(0)

    # --- Tab builders ---

    def _build_general_tab(self):
        tab = QWidget()
        layout = QFormLayout(tab)

        self.name_input = QLineEdit()
        layout.addRow("Name:", self.name_input)

        self.sol_combo = QComboBox()
        for sol_val, sol_label in SOL_TYPES:
            self.sol_combo.addItem(sol_label, sol_val)
        self.sol_combo.currentIndexChanged.connect(self._on_sol_changed)
        layout.addRow("Solution Type:", self.sol_combo)

        self.tabs.addTab(tab, "General")

    def _build_scope_tab(self):
        """v5.1.0 item 26b: Scope tab.

        Lets the user pick the solver target (MYSTRAN / MSC / NX /
        Generic), per-set field format, geometric scope (Full model /
        a specific group), and explicit LOAD / SPC SID white-lists.

        Smart options: when solver_target == 'MYSTRAN', the dialog
        shows a banner reminding the user that MSC-only PARAMs will be
        dropped on export.
        """
        from PySide6.QtWidgets import QListWidget, QListWidgetItem
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setSpacing(6)

        intro = QLabel(
            "Configure which solver this AnalysisSet targets and what "
            "subset of the model gets exported when this set is "
            "active.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        outer.addWidget(intro)

        form = QFormLayout()

        self.solver_target_combo = QComboBox()
        for key, label in (
                ('MYSTRAN', 'MYSTRAN (open-source)'),
                ('MSC',     'MSC Nastran'),
                ('NX',      'Siemens NX Nastran'),
                ('generic', 'Generic (MSC/NX compatible)')):
            self.solver_target_combo.addItem(label, key)
        self.solver_target_combo.currentIndexChanged.connect(
            self._on_solver_target_changed)
        form.addRow("Solver target:", self.solver_target_combo)

        self.field_format_combo = QComboBox()
        for key, label in (('short', 'Short Field (8-character)'),
                           ('long', 'Long Field (16-character)'),
                           ('free', 'Free Field (comma-separated)')):
            self.field_format_combo.addItem(label, key)
        form.addRow("Field format:", self.field_format_combo)

        self.group_target_combo = QComboBox()
        self.group_target_combo.addItem("Full Model", None)
        for gname in sorted(self._groups.keys()):
            self.group_target_combo.addItem(f"Group: {gname}", gname)
        self.group_target_combo.setToolTip(
            "When set to a Group, the export only includes elements "
            "in that group plus auto-collected property / material / "
            "coord-system dependencies. 'Full Model' writes everything.")
        form.addRow("Geometric scope:", self.group_target_combo)

        outer.addLayout(form)

        # MYSTRAN smart banner (visible only when solver_target=MYSTRAN)
        self._mystran_banner = QLabel(
            "<b>MYSTRAN mode:</b> MSC-only PARAMs (POST, COUPMASS, "
            "AUTOSPC, etc.) will be dropped on export. Pre-flight "
            "blocks aero / nonlinear / contact / optimization decks.")
        self._mystran_banner.setWordWrap(True)
        self._mystran_banner.setStyleSheet(
            "background-color: #f9e2af; color: #1e1e2e; "
            "padding: 8px; border-radius: 4px;")
        outer.addWidget(self._mystran_banner)

        # Load + SPC SID pickers
        loads_box = QGroupBox("Loads to include (empty = all)")
        loads_lay = QVBoxLayout(loads_box)
        self.load_sids_list = QListWidget()
        self.load_sids_list.setSelectionMode(QListWidget.NoSelection)
        loads_lay.addWidget(self.load_sids_list)
        outer.addWidget(loads_box, 1)

        spcs_box = QGroupBox("SPC sets to include (empty = all)")
        spcs_lay = QVBoxLayout(spcs_box)
        self.spc_sids_list = QListWidget()
        self.spc_sids_list.setSelectionMode(QListWidget.NoSelection)
        spcs_lay.addWidget(self.spc_sids_list)
        outer.addWidget(spcs_box, 1)

        # Populate from the loaded model (these refresh on each set
        # selection in _populate_scope_tab below).
        self._populate_scope_sid_pickers()

        self.tabs.addTab(tab, "Scope")

    def _populate_scope_sid_pickers(self):
        """(Re)populate the LOAD / SPC SID checkable lists from the
        current model. Called once per set load."""
        from PySide6.QtWidgets import QListWidgetItem
        self.load_sids_list.clear()
        self.spc_sids_list.clear()
        if not self._model:
            return
        for sid in sorted((self._model.loads or {}).keys()):
            n_cards = len((self._model.loads or {}).get(sid, []))
            item = QListWidgetItem(f"LOAD SID {sid}  ({n_cards} card(s))")
            item.setData(QtCore.Qt.UserRole, int(sid))
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Checked)
            self.load_sids_list.addItem(item)
        for sid in sorted((self._model.spcs or {}).keys()):
            n_cards = len((self._model.spcs or {}).get(sid, []))
            item = QListWidgetItem(f"SPC SID {sid}  ({n_cards} card(s))")
            item.setData(QtCore.Qt.UserRole, int(sid))
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Checked)
            self.spc_sids_list.addItem(item)

    def _on_solver_target_changed(self, _idx):
        """Toggle the MYSTRAN-mode banner visibility."""
        if not hasattr(self, '_mystran_banner'):
            return
        target = self.solver_target_combo.currentData() or 'MYSTRAN'
        self._mystran_banner.setVisible(target == 'MYSTRAN')

    def _build_output_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Preset buttons
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Presets:"))
        for name in OUTPUT_PRESETS:
            btn = QPushButton(name)
            btn.clicked.connect(lambda checked, n=name: self._apply_output_preset(n))
            preset_layout.addWidget(btn)
        preset_layout.addStretch()
        layout.addLayout(preset_layout)

        # Output request grid
        form = QFormLayout()
        self._output_combos = {}
        output_labels = {
            'displacement': 'Displacement', 'velocity': 'Velocity',
            'acceleration': 'Acceleration', 'stress': 'Stress',
            'strain': 'Strain', 'force': 'Element Force',
            'spcforce': 'SPC Forces', 'mpcforce': 'MPC Forces',
            'gpforce': 'Grid Point Forces', 'oload': 'Applied Loads (OLOAD)',
        }
        for key, label in output_labels.items():
            combo = QComboBox()
            combo.addItems(["NONE", "ALL"])
            self._output_combos[key] = combo
            form.addRow(f"{label}:", combo)

        layout.addLayout(form)
        layout.addStretch()
        self.tabs.addTab(tab, "Output Requests")

    def _build_eigrl_tab(self):
        tab = QWidget()
        layout = QFormLayout(tab)

        self.eigrl_sid_input = QSpinBox()
        self.eigrl_sid_input.setRange(1, 999999)
        self.eigrl_sid_input.setValue(100)
        layout.addRow("SID:", self.eigrl_sid_input)

        self.eigrl_v1_input = QLineEdit("")
        self.eigrl_v1_input.setPlaceholderText("Optional - lower frequency bound")
        self.eigrl_v1_input.setValidator(QDoubleValidator())
        layout.addRow("V1 (Lower Freq):", self.eigrl_v1_input)

        self.eigrl_v2_input = QLineEdit("")
        self.eigrl_v2_input.setPlaceholderText("Optional - upper frequency bound")
        self.eigrl_v2_input.setValidator(QDoubleValidator())
        layout.addRow("V2 (Upper Freq):", self.eigrl_v2_input)

        self.eigrl_nd_input = QSpinBox()
        self.eigrl_nd_input.setRange(1, 10000)
        self.eigrl_nd_input.setValue(10)
        layout.addRow("ND (Number of Modes):", self.eigrl_nd_input)

        self.eigrl_tab_widget = tab
        self.tabs.addTab(tab, "Eigenvalue (EIGRL)")

    def _build_params_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # --- Common Parameters ---
        common_group = QGroupBox("Common Parameters")
        common_form = QFormLayout(common_group)
        common_form.setSpacing(6)

        self._common_param_widgets = {}  # {param_name: (wtype, widget)}

        for param_name, label, wtype, default, options in COMMON_PARAMS:
            if wtype == 'combo':
                widget = QComboBox()
                for val, text in options:
                    widget.addItem(text, val)
                idx = widget.findData(default)
                if idx >= 0:
                    widget.setCurrentIndex(idx)
            elif wtype == 'yesno':
                widget = QComboBox()
                widget.addItem("YES", "YES")
                widget.addItem("NO", "NO")
                idx = widget.findData(default)
                if idx >= 0:
                    widget.setCurrentIndex(idx)
            elif wtype == 'int':
                widget = QSpinBox()
                widget.setRange(options[0], options[1])
                widget.setValue(default)
            elif wtype == 'float':
                widget = QLineEdit(str(default))
                widget.setValidator(QDoubleValidator())
            else:
                continue

            self._common_param_widgets[param_name] = (wtype, widget)
            common_form.addRow(f"{label}:", widget)

        layout.addWidget(common_group)

        # --- Custom Parameters (free-form) ---
        custom_group = QGroupBox("Custom Parameters")
        custom_layout = QVBoxLayout(custom_group)

        self.params_table = QTableWidget()
        self.params_table.setColumnCount(2)
        self.params_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self.params_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        custom_layout.addWidget(self.params_table)

        btn_layout = QHBoxLayout()
        add_btn = QPushButton("Add Parameter")
        remove_btn = QPushButton("Remove Selected")
        add_btn.clicked.connect(self._add_param_row)
        remove_btn.clicked.connect(self._remove_param_row)
        btn_layout.addWidget(add_btn)
        btn_layout.addWidget(remove_btn)
        btn_layout.addStretch()
        custom_layout.addLayout(btn_layout)

        layout.addWidget(custom_group)

        self.tabs.addTab(tab, "Parameters")

    def _build_subcases_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        info_label = QLabel(
            "Define subcases: each references a load set, SPC set, and output requests. "
            "For SOL 105 buckling, set METHOD to the EIGRL SID and STATSUB to the preload subcase.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("font-size: 11px; color: #a6adc8;")
        layout.addWidget(info_label)

        self.subcases_table = QTableWidget()
        self.subcases_table.setColumnCount(9)
        self.subcases_table.setHorizontalHeaderLabels([
            "ID", "Load SID", "SPC SID", "METHOD", "STATSUB",
            "DISP", "STRESS", "FORCE", "STRAIN"
        ])
        self.subcases_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addWidget(self.subcases_table)

        btn_layout = QHBoxLayout()
        add_sc_btn = QPushButton("Add Subcase")
        remove_sc_btn = QPushButton("Remove Selected")
        gen_combos_btn = QPushButton("Generate All Combos")
        add_sc_btn.clicked.connect(self._add_subcase_row)
        remove_sc_btn.clicked.connect(self._remove_subcase_row)
        gen_combos_btn.clicked.connect(self._generate_all_combos)
        select_combos_btn = QPushButton("Select Combos...")
        select_combos_btn.clicked.connect(self._select_combo_matrix)
        btn_layout.addWidget(add_sc_btn)
        btn_layout.addWidget(remove_sc_btn)
        btn_layout.addWidget(gen_combos_btn)
        btn_layout.addWidget(select_combos_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        self.tabs.addTab(tab, "Subcases")

    # --- Set list management ---

    def _populate_set_list(self):
        self.set_list.blockSignals(True)
        self.set_list.clear()
        for sid, aset in sorted(self._analysis_sets.items()):
            prefix = "\u2605 " if sid == self._active_set_id else "  "
            sol_str = f"SOL {aset.sol_type}" if aset.sol_type else "Punch"
            self.set_list.addItem(f"{prefix}{sid}: {aset.name} ({sol_str})")
        self.set_list.blockSignals(False)

    def _on_set_selected(self, row):
        """Load the selected analysis set into the right panel tabs."""
        if row < 0:
            return
        self._save_current_to_set()  # Save any pending edits first
        sid = sorted(self._analysis_sets.keys())[row]
        self._current_set_id = sid
        aset = self._analysis_sets[sid]

        # General tab
        self.name_input.setText(aset.name)
        sol_idx = next((i for i, (v, _) in enumerate(SOL_TYPES) if v == aset.sol_type), 0)
        self.sol_combo.setCurrentIndex(sol_idx)

        # v5.1.0 item 26b: Scope tab
        target = getattr(aset, 'solver_target', 'MYSTRAN') or 'MYSTRAN'
        idx = self.solver_target_combo.findData(target)
        self.solver_target_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self._on_solver_target_changed(idx)
        ff = getattr(aset, 'field_format', 'short') or 'short'
        idx = self.field_format_combo.findData(ff)
        self.field_format_combo.setCurrentIndex(idx if idx >= 0 else 0)
        gt = getattr(aset, 'group_target', None)
        idx = self.group_target_combo.findData(gt)
        self.group_target_combo.setCurrentIndex(idx if idx >= 0 else 0)
        load_keep = set(int(s) for s in (getattr(aset, 'load_sids', []) or []))
        spc_keep = set(int(s) for s in (getattr(aset, 'spc_sids', []) or []))
        # Empty list -> "include all" -> every checkbox checked.
        for i in range(self.load_sids_list.count()):
            it = self.load_sids_list.item(i)
            sid = int(it.data(QtCore.Qt.UserRole))
            it.setCheckState(
                QtCore.Qt.Checked if (not load_keep or sid in load_keep)
                else QtCore.Qt.Unchecked)
        for i in range(self.spc_sids_list.count()):
            it = self.spc_sids_list.item(i)
            sid = int(it.data(QtCore.Qt.UserRole))
            it.setCheckState(
                QtCore.Qt.Checked if (not spc_keep or sid in spc_keep)
                else QtCore.Qt.Unchecked)

        # Output tab
        for key, combo in self._output_combos.items():
            combo.setCurrentText(aset.output_requests.get(key, 'NONE'))

        # EIGRL tab
        if aset.eigrl:
            self.eigrl_sid_input.setValue(aset.eigrl.get('sid', 100))
            self.eigrl_v1_input.setText(str(aset.eigrl['v1']) if aset.eigrl.get('v1') is not None else "")
            self.eigrl_v2_input.setText(str(aset.eigrl['v2']) if aset.eigrl.get('v2') is not None else "")
            self.eigrl_nd_input.setValue(aset.eigrl.get('nd', 10))
        else:
            self.eigrl_sid_input.setValue(100)
            self.eigrl_v1_input.clear()
            self.eigrl_v2_input.clear()
            self.eigrl_nd_input.setValue(10)

        # EIGRL tab enabled only for SOL 103/105
        eigrl_enabled = aset.sol_type in (103, 105)
        self.eigrl_tab_widget.setEnabled(eigrl_enabled)

        # Params tab - Common widgets
        common_param_names = {name for name, _, _, _, _ in COMMON_PARAMS}
        for param_name, (wtype, widget) in self._common_param_widgets.items():
            value = aset.params.get(param_name)
            default = next(
                (d for n, _, _, d, _ in COMMON_PARAMS if n == param_name),
                None)
            effective = value if value is not None else default
            if wtype in ('combo', 'yesno'):
                idx = widget.findData(effective)
                if idx >= 0:
                    widget.setCurrentIndex(idx)
                else:
                    widget.setCurrentText(str(effective))
            elif wtype == 'int':
                try:
                    widget.setValue(int(effective))
                except (ValueError, TypeError):
                    widget.setValue(0)
            elif wtype == 'float':
                widget.setText(str(effective) if effective is not None else "")

        # Params tab - Custom table (non-common params only)
        self.params_table.setRowCount(0)
        for pname, pval in aset.params.items():
            if pname not in common_param_names:
                row_idx = self.params_table.rowCount()
                self.params_table.setRowCount(row_idx + 1)
                self.params_table.setItem(row_idx, 0,
                                          QTableWidgetItem(str(pname)))
                self.params_table.setItem(row_idx, 1,
                                          QTableWidgetItem(str(pval)))

        # Subcases tab
        self.subcases_table.setRowCount(0)
        for sc in aset.subcases:
            self._add_subcase_row(sc)

    def _save_current_to_set(self):
        """Save the current UI state back to the current analysis set."""
        if not hasattr(self, '_current_set_id') or self._current_set_id not in self._analysis_sets:
            return
        aset = self._analysis_sets[self._current_set_id]
        aset.name = self.name_input.text() or f"Set {aset.id}"
        aset.sol_type = self.sol_combo.currentData()

        # v5.1.0 item 26b: Scope tab fields
        if hasattr(self, 'solver_target_combo'):
            aset.solver_target = (
                self.solver_target_combo.currentData() or 'MYSTRAN')
        if hasattr(self, 'field_format_combo'):
            aset.field_format = (
                self.field_format_combo.currentData() or 'short')
        if hasattr(self, 'group_target_combo'):
            aset.group_target = self.group_target_combo.currentData()
        if hasattr(self, 'load_sids_list'):
            checked = []
            unchecked = []
            for i in range(self.load_sids_list.count()):
                it = self.load_sids_list.item(i)
                sid = int(it.data(QtCore.Qt.UserRole))
                if it.checkState() == QtCore.Qt.Checked:
                    checked.append(sid)
                else:
                    unchecked.append(sid)
            # Empty list semantically = "include all"; only persist a
            # non-empty list if the user actually narrowed the set.
            aset.load_sids = checked if unchecked else []
        if hasattr(self, 'spc_sids_list'):
            checked = []
            unchecked = []
            for i in range(self.spc_sids_list.count()):
                it = self.spc_sids_list.item(i)
                sid = int(it.data(QtCore.Qt.UserRole))
                if it.checkState() == QtCore.Qt.Checked:
                    checked.append(sid)
                else:
                    unchecked.append(sid)
            aset.spc_sids = checked if unchecked else []

        # Output requests
        for key, combo in self._output_combos.items():
            aset.output_requests[key] = combo.currentText()

        # EIGRL
        if aset.sol_type in (103, 105):
            v1_text = self.eigrl_v1_input.text().strip()
            v2_text = self.eigrl_v2_input.text().strip()
            aset.eigrl = {
                'sid': self.eigrl_sid_input.value(),
                'v1': float(v1_text) if v1_text else None,
                'v2': float(v2_text) if v2_text else None,
                'nd': self.eigrl_nd_input.value(),
            }
        else:
            aset.eigrl = None

        # Params - Common widgets
        aset.params = {}
        for param_name, (wtype, widget) in self._common_param_widgets.items():
            if wtype == 'combo':
                aset.params[param_name] = widget.currentData()
            elif wtype == 'yesno':
                aset.params[param_name] = widget.currentData()
            elif wtype == 'int':
                aset.params[param_name] = widget.value()
            elif wtype == 'float':
                text = widget.text().strip()
                if text:
                    try:
                        val = float(text)
                        if val == int(val) and 'e' not in text.lower():
                            val = int(val)
                        aset.params[param_name] = val
                    except ValueError:
                        aset.params[param_name] = text

        # Params - Custom table
        for r in range(self.params_table.rowCount()):
            pname_item = self.params_table.item(r, 0)
            pval_item = self.params_table.item(r, 1)
            if pname_item and pval_item:
                pname = pname_item.text().strip()
                pval_str = pval_item.text().strip()
                if pname:
                    try:
                        pval = int(pval_str)
                    except ValueError:
                        try:
                            pval = float(pval_str)
                        except ValueError:
                            pval = pval_str
                    aset.params[pname] = pval

        # Subcases
        aset.subcases = self._read_subcases_from_table()

    def _new_set(self):
        """Create a new analysis set."""
        new_id = max(self._analysis_sets.keys(), default=0) + 1
        aset = AnalysisSet(id=new_id, name=f"Analysis Set {new_id}")
        self._analysis_sets[new_id] = aset
        self._populate_set_list()
        self.set_list.setCurrentRow(self.set_list.count() - 1)

    def _duplicate_set(self):
        """Duplicate the currently selected analysis set."""
        if not hasattr(self, '_current_set_id') or self._current_set_id not in self._analysis_sets:
            return
        self._save_current_to_set()
        src = self._analysis_sets[self._current_set_id]
        new_id = max(self._analysis_sets.keys(), default=0) + 1
        new_set = self._copy_set(src)
        new_set.id = new_id
        new_set.name = f"{src.name} (Copy)"
        self._analysis_sets[new_id] = new_set
        self._populate_set_list()
        self.set_list.setCurrentRow(self.set_list.count() - 1)

    def _delete_set(self):
        """Delete the currently selected analysis set."""
        if not hasattr(self, '_current_set_id') or self._current_set_id not in self._analysis_sets:
            return
        if len(self._analysis_sets) <= 1:
            QMessageBox.information(self, "Cannot Delete", "At least one analysis set must remain.")
            return
        del self._analysis_sets[self._current_set_id]
        if self._active_set_id == self._current_set_id:
            self._active_set_id = next(iter(self._analysis_sets.keys()))
        self._populate_set_list()
        if self.set_list.count() > 0:
            self.set_list.setCurrentRow(0)

    def _set_active(self):
        """Mark the current set as active."""
        if hasattr(self, '_current_set_id') and self._current_set_id in self._analysis_sets:
            self._active_set_id = self._current_set_id
            self._populate_set_list()
            # Re-select same row
            keys = sorted(self._analysis_sets.keys())
            idx = keys.index(self._current_set_id)
            self.set_list.setCurrentRow(idx)

    # --- Output presets ---

    def _apply_output_preset(self, preset_name):
        preset = OUTPUT_PRESETS.get(preset_name, {})
        for key, combo in self._output_combos.items():
            combo.setCurrentText(preset.get(key, 'NONE'))

    def _on_sol_changed(self, index):
        sol_val = self.sol_combo.currentData()
        eigrl_enabled = sol_val in (103, 105)
        self.eigrl_tab_widget.setEnabled(eigrl_enabled)

    # --- Params table ---

    def _add_param_row(self):
        row = self.params_table.rowCount()
        self.params_table.setRowCount(row + 1)
        self.params_table.setItem(row, 0, QTableWidgetItem(""))
        self.params_table.setItem(row, 1, QTableWidgetItem(""))

    def _remove_param_row(self):
        row = self.params_table.currentRow()
        if row >= 0:
            self.params_table.removeRow(row)

    # --- Subcases table ---

    def _make_output_combo(self, default="NONE"):
        combo = QComboBox()
        combo.addItems(["NONE", "ALL"])
        combo.setCurrentText(default)
        return combo

    def _make_sid_combo(self, sid_list, default=None):
        combo = QComboBox()
        combo.addItem("\u2014")
        for sid in sid_list:
            combo.addItem(str(sid))
        if default is not None:
            combo.setCurrentText(str(default))
        return combo

    def _get_load_sids(self):
        if self._model:
            sids = sorted(self._model.loads.keys())
            if hasattr(self._model, 'load_combinations'):
                sids = sorted(set(sids) | set(self._model.load_combinations.keys()))
            return sids
        return []

    def _get_spc_sids(self):
        if self._model:
            return sorted(self._model.spcs.keys())
        return []

    def _add_subcase_row(self, data=None):
        row = self.subcases_table.rowCount()
        self.subcases_table.setRowCount(row + 1)

        load_sids = self._get_load_sids()
        spc_sids = self._get_spc_sids()

        if isinstance(data, dict):
            sc_id = data.get('id', row + 1)
            load_sid = data.get('load_sid')
            spc_sid = data.get('spc_sid')
            method_sid = data.get('method_sid')
            statsub_id = data.get('statsub_id')
            disp = data.get('disp', 'NONE')
            stress = data.get('stress', 'NONE')
            force = data.get('force', 'NONE')
            strain = data.get('strain', 'NONE')
        else:
            sc_id = row + 1
            load_sid = spc_sid = method_sid = statsub_id = None
            disp = stress = force = strain = "NONE"

        self.subcases_table.setItem(row, 0, QTableWidgetItem(str(sc_id)))
        self.subcases_table.setCellWidget(row, 1, self._make_sid_combo(load_sids, load_sid))
        self.subcases_table.setCellWidget(row, 2, self._make_sid_combo(spc_sids, spc_sid))

        # METHOD - editable combo
        method_combo = QComboBox()
        method_combo.setEditable(True)
        method_combo.addItem("\u2014")
        if method_sid is not None:
            method_combo.setCurrentText(str(method_sid))
        self.subcases_table.setCellWidget(row, 3, method_combo)

        # STATSUB - editable combo
        statsub_combo = QComboBox()
        statsub_combo.setEditable(True)
        statsub_combo.addItem("\u2014")
        if statsub_id is not None:
            statsub_combo.setCurrentText(str(statsub_id))
        self.subcases_table.setCellWidget(row, 4, statsub_combo)

        self.subcases_table.setCellWidget(row, 5, self._make_output_combo(disp))
        self.subcases_table.setCellWidget(row, 6, self._make_output_combo(stress))
        self.subcases_table.setCellWidget(row, 7, self._make_output_combo(force))
        self.subcases_table.setCellWidget(row, 8, self._make_output_combo(strain))

    def _remove_subcase_row(self):
        row = self.subcases_table.currentRow()
        if row >= 0:
            self.subcases_table.removeRow(row)

    def _generate_all_combos(self):
        """Generate one subcase per (Load SID × SPC SID) combination."""
        load_sids = self._get_load_sids()
        spc_sids = self._get_spc_sids()
        if not load_sids or not spc_sids:
            QMessageBox.information(self, "No Data",
                                    "Need both load sets and SPC sets to generate combinations.")
            return

        self.subcases_table.setRowCount(0)
        sc_id = 1
        for load_sid in load_sids:
            for spc_sid in spc_sids:
                self._add_subcase_row({
                    'id': sc_id,
                    'load_sid': load_sid,
                    'spc_sid': spc_sid,
                    'disp': 'ALL', 'stress': 'ALL', 'force': 'ALL', 'strain': 'NONE',
                })
                sc_id += 1

    def _select_combo_matrix(self):
        """Open checkbox matrix to select specific Load×SPC combinations."""
        load_sids = self._get_load_sids()
        spc_sids = self._get_spc_sids()
        if not load_sids or not spc_sids:
            QMessageBox.information(
                self, "No Data",
                "Need both load sets and SPC sets to select combinations.")
            return

        dialog = SubcaseMatrixDialog(load_sids, spc_sids, parent=self)
        if dialog.exec():
            combos = dialog.get_selected_combinations()
            if not combos:
                QMessageBox.information(
                    self, "No Selection",
                    "No combinations were selected.")
                return

            self.subcases_table.setRowCount(0)
            sc_id = 1
            for load_sid, spc_sid in combos:
                self._add_subcase_row({
                    'id': sc_id,
                    'load_sid': load_sid,
                    'spc_sid': spc_sid,
                    'disp': 'ALL', 'stress': 'ALL',
                    'force': 'ALL', 'strain': 'NONE',
                })
                sc_id += 1

    def _read_subcases_from_table(self):
        """Read subcases from the table widget."""
        subcases = []
        for row in range(self.subcases_table.rowCount()):
            try:
                sc_id = int(self.subcases_table.item(row, 0).text())
            except (ValueError, AttributeError):
                continue

            def _combo_int(col):
                w = self.subcases_table.cellWidget(row, col)
                if not w:
                    return None
                text = w.currentText()
                if text == "\u2014" or not text.strip():
                    return None
                try:
                    return int(text)
                except ValueError:
                    return None

            def _combo_text(col):
                w = self.subcases_table.cellWidget(row, col)
                return w.currentText() if w else "NONE"

            subcases.append({
                'id': sc_id,
                'load_sid': _combo_int(1),
                'spc_sid': _combo_int(2),
                'method_sid': _combo_int(3),
                'statsub_id': _combo_int(4),
                'disp': _combo_text(5),
                'stress': _combo_text(6),
                'force': _combo_text(7),
                'strain': _combo_text(8),
            })
        return subcases

    # --- Accept / Result ---

    def _on_accept(self):
        self._save_current_to_set()
        self.accept()

    def get_results(self):
        """Return the analysis sets dict and active set ID."""
        return self._analysis_sets, self._active_set_id

    @staticmethod
    def _copy_set(aset: AnalysisSet) -> AnalysisSet:
        """Deep copy an AnalysisSet."""
        return AnalysisSet(
            id=aset.id,
            name=aset.name,
            sol_type=aset.sol_type,
            output_requests=dict(aset.output_requests),
            eigrl=dict(aset.eigrl) if aset.eigrl else None,
            params=dict(aset.params),
            subcases=[dict(sc) for sc in aset.subcases],
        )


# ---------------------------------------------------------------------------
# Output Requests Quick Dialog
# ---------------------------------------------------------------------------

class OutputRequestsDialog(QDialog):
    """Quick dialog for editing global output requests of the active analysis set."""

    def __init__(self, output_requests: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Output Requests")
        self.setMinimumWidth(400)

        layout = QVBoxLayout(self)

        # Preset buttons
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Presets:"))
        for name in OUTPUT_PRESETS:
            btn = QPushButton(name)
            btn.clicked.connect(lambda checked, n=name: self._apply_preset(n))
            preset_layout.addWidget(btn)
        preset_layout.addStretch()
        layout.addLayout(preset_layout)

        # Output combos
        form = QFormLayout()
        self._combos = {}
        output_labels = {
            'displacement': 'Displacement', 'velocity': 'Velocity',
            'acceleration': 'Acceleration', 'stress': 'Stress',
            'strain': 'Strain', 'force': 'Element Force',
            'spcforce': 'SPC Forces', 'mpcforce': 'MPC Forces',
            'gpforce': 'Grid Point Forces', 'oload': 'Applied Loads (OLOAD)',
        }
        for key, label in output_labels.items():
            combo = QComboBox()
            combo.addItems(["NONE", "ALL"])
            combo.setCurrentText(output_requests.get(key, 'NONE'))
            self._combos[key] = combo
            form.addRow(f"{label}:", combo)
        layout.addLayout(form)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _apply_preset(self, preset_name):
        preset = OUTPUT_PRESETS.get(preset_name, {})
        for key, combo in self._combos.items():
            combo.setCurrentText(preset.get(key, 'NONE'))

    def get_output_requests(self):
        return {key: combo.currentText() for key, combo in self._combos.items()}


# ---------------------------------------------------------------------------
# Subcase Combination Matrix Dialog
# ---------------------------------------------------------------------------

class SubcaseMatrixDialog(QDialog):
    """Checkbox matrix for selecting Load SID × SPC SID combinations."""

    def __init__(self, load_sids: list[int], spc_sids: list[int],
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Load / SPC Combinations")
        self.setMinimumSize(500, 400)

        self._load_sids = load_sids
        self._spc_sids = spc_sids

        layout = QVBoxLayout(self)

        info = QLabel(
            "Check the Load SID / SPC SID intersections to generate "
            "subcases for. Each checked cell becomes one subcase.")
        info.setWordWrap(True)
        layout.addWidget(info)

        # --- Matrix table ---
        self.matrix_table = QTableWidget()
        self.matrix_table.setRowCount(len(load_sids))
        self.matrix_table.setColumnCount(len(spc_sids))
        self.matrix_table.setHorizontalHeaderLabels(
            [f"SPC {s}" for s in spc_sids])
        self.matrix_table.setVerticalHeaderLabels(
            [f"LOAD {s}" for s in load_sids])

        for r in range(len(load_sids)):
            for c in range(len(spc_sids)):
                cb = QCheckBox()
                self.matrix_table.setCellWidget(r, c, cb)

        self.matrix_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.matrix_table.verticalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        layout.addWidget(self.matrix_table)

        # --- Action buttons ---
        action_layout = QHBoxLayout()
        select_all_btn = QPushButton("Select All")
        clear_all_btn = QPushButton("Clear All")
        invert_btn = QPushButton("Invert")
        select_all_btn.clicked.connect(self._select_all)
        clear_all_btn.clicked.connect(self._clear_all)
        invert_btn.clicked.connect(self._invert)
        action_layout.addWidget(select_all_btn)
        action_layout.addWidget(clear_all_btn)
        action_layout.addWidget(invert_btn)
        action_layout.addStretch()
        layout.addLayout(action_layout)

        # --- OK / Cancel ---
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def _select_all(self):
        for r in range(self.matrix_table.rowCount()):
            for c in range(self.matrix_table.columnCount()):
                self.matrix_table.cellWidget(r, c).setChecked(True)

    def _clear_all(self):
        for r in range(self.matrix_table.rowCount()):
            for c in range(self.matrix_table.columnCount()):
                self.matrix_table.cellWidget(r, c).setChecked(False)

    def _invert(self):
        for r in range(self.matrix_table.rowCount()):
            for c in range(self.matrix_table.columnCount()):
                cb = self.matrix_table.cellWidget(r, c)
                cb.setChecked(not cb.isChecked())

    def get_selected_combinations(self) -> list[tuple[int, int]]:
        """Return list of (load_sid, spc_sid) for checked cells."""
        combos = []
        for r in range(self.matrix_table.rowCount()):
            for c in range(self.matrix_table.columnCount()):
                if self.matrix_table.cellWidget(r, c).isChecked():
                    combos.append((self._load_sids[r], self._spc_sids[c]))
        return combos
