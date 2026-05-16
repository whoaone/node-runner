"""v5.2.0 item 42: professional Post-Processing Toolbox.

Replaces the v5.1.x tree-and-sub-tabs design with three collapsible
sections (Results / Deform / Contour) and dropdown-driven output set
+ output vector picking. Tables / Animation / Vectors sub-tabs have
been removed; tabular results live under Tools -> Data Table now.

This widget is **presentation-only** -- it emits signals; MainWindow
handles the actual rendering / actor builders / animation timer.

Signals
-------
output_set_changed(int subcase_id)
output_vector_changed(str kind, str component, int mode_idx)
    ``kind`` in {'displacement','stress','spc_forces','eigenvector'}.
    ``component`` is a canonical key (e.g. 'T3', 'von_mises', 'Magnitude').
    ``mode_idx`` is 0-based for eigenvector kinds, else -1.

deform_style_changed(str)               'undeformed' | 'deformed' | 'animate'
deform_scale_changed(str mode, float value)   mode in {'pct','actual'}
animation_settings_changed(str mode, int n_frames, int delay_ms)
                                        mode in {'sine','modes'}
show_undeformed_ref_changed(bool)

contour_style_changed(str)              'filled' | 'filled_edges' | 'bands'
data_conversion_changed(str)            'average' | 'no_avg' | 'max_node' | 'min_node'
levels_changed(int)
palette_changed(str)                    'jet' | 'viridis' | 'plasma' | 'coolwarm' | 'gray'
color_range_changed(bool auto, float vmin, float vmax)
markers_changed(bool show_min, bool show_max)
element_labels_changed(bool show, int top_n)

Public methods
--------------
populate_from_results(op2_results, file_label='')
    Rebuild Output Set + Output Vector dropdowns from a results bundle.
current_subcase_id() -> int | None
current_output_vector() -> tuple[str, str, int]
set_data_range(vmin, vmax)
    Called by MainWindow after a contour build so the Manual min/max
    spinboxes start at the data range when the user toggles Auto OFF.
"""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QComboBox,
    QGroupBox, QCheckBox, QSpinBox, QDoubleSpinBox, QRadioButton,
    QButtonGroup, QSizePolicy, QFrame,
)

from node_runner.dialogs.collapsible_section import CollapsibleSection


# ---------------------------------------------------------------------------
# Output Vector dropdown helpers
# ---------------------------------------------------------------------------

# Component canonical keys + display labels used in the Output Vector combo.
DISP_COMPONENTS = [
    ("Magnitude",     "Magnitude"),
    ("T1 (X)",        "T1"),
    ("T2 (Y)",        "T2"),
    ("T3 (Z)",        "T3"),
    ("Rotation Mag",  "RMAG"),
    ("R1 (RX)",       "R1"),
    ("R2 (RY)",       "R2"),
    ("R3 (RZ)",       "R3"),
]
STRESS_COMPONENTS = [
    ("von Mises",       "von_mises"),
    ("Max Principal",   "max_principal"),
    ("Min Principal",   "min_principal"),
    ("σxx",        "oxx"),
    ("σyy",        "oyy"),
    ("τxy",        "txy"),
]
SPC_COMPONENTS = [
    ("Magnitude",     "Magnitude"),
    ("T1 (X)",        "T1"),
    ("T2 (Y)",        "T2"),
    ("T3 (Z)",        "T3"),
]


PALETTES = [
    ("Jet",        "jet"),
    ("Viridis",    "viridis"),
    ("Plasma",     "plasma"),
    ("Coolwarm",   "coolwarm"),
    ("Grayscale",  "gray"),
]
DEFORM_STYLES = [
    ("Undeformed", "undeformed"),
    ("Deformed",   "deformed"),
    ("Animate",    "animate"),
]
CONTOUR_STYLES = [
    # v5.3.0 item 56: No Contour first so it's the default. Common tools load
    # results with "None - Model Only" contour by default; the user
    # explicitly turns on Filled / Filled+Edges / Discrete Bands.
    ("No Contour",     "none"),
    ("Filled",         "filled"),
    ("Filled + Edges", "filled_edges"),
    ("Discrete Bands", "bands"),
]
DATA_CONVERSIONS = [
    ("Average (nodal)",         "average"),
    ("No Averaging (centroid)", "no_avg"),
    ("Max at Node",             "max_node"),
    ("Min at Node",             "min_node"),
]
ANIM_MODES = [
    ("Sine wave",    "sine"),
    ("Through Modes", "modes"),
]


# ---------------------------------------------------------------------------
# ResultsTab
# ---------------------------------------------------------------------------

class ResultsTab(QWidget):
    """The Post-Processing Toolbox sidebar widget.

    Layout: three collapsible sections (Results, Deform, Contour)
    stacked vertically. The widget owns no rendering state -- it
    just emits signals; MainWindow drives the plotter.
    """

    # ---- Selection ----
    output_set_changed = Signal(int)
    output_vector_changed = Signal(str, str, int)   # kind, component, mode_idx
    # v5.3.2 item 67: v5.3.0's dual deform/contour vector design was
    # over-engineered. Reverted to a single Output Vector. Style
    # toggles in the Deform / Contour sections decide what renders;
    # deformation always uses the displacement triplet from the same
    # subcase regardless of the selected Output Vector.
    # v5.3.2 item 69: bar_orientation_changed signal removed -- the
    # setting moved to File > Preferences.

    # ---- Deform ----
    deform_style_changed = Signal(str)              # 'undeformed'|'deformed'|'animate'
    deform_scale_changed = Signal(str, float)       # mode, value
    animation_settings_changed = Signal(str, int, int)  # mode, n_frames, delay_ms
    show_undeformed_ref_changed = Signal(bool)

    # ---- Contour ----
    contour_style_changed = Signal(str)
    data_conversion_changed = Signal(str)
    levels_changed = Signal(int)
    palette_changed = Signal(str)
    color_range_changed = Signal(bool, float, float)
    markers_changed = Signal(bool, bool)
    element_labels_changed = Signal(bool, int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._loading = False
        self._has_eigenvectors = False
        self._current_subcase_id: Optional[int] = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # ------------------------------------------------------------
        # Section 1: RESULTS (output set + output vector)
        # ------------------------------------------------------------
        # v5.3.2 item 67: ONE Output Vector dropdown lives here -- it
        # drives BOTH the contour AND the deformation (deformation
        # always interprets the displacement triplet from the same
        # subcase, regardless of which vector is selected). The Deform
        # and Contour sections each have their own Style toggle that
        # turn rendering on/off INDEPENDENTLY of each other.
        self._results_section = CollapsibleSection("Results", self)
        results_body = QWidget()
        rb = QFormLayout(results_body)
        rb.setContentsMargins(0, 4, 0, 4)
        rb.setSpacing(4)
        self.output_set_combo = QComboBox()
        self.output_set_combo.setSizePolicy(QSizePolicy.Expanding,
                                            QSizePolicy.Fixed)
        self.output_set_combo.currentIndexChanged.connect(
            self._on_output_set_changed)
        rb.addRow("Output Set:", self.output_set_combo)
        self.output_vector_combo = QComboBox()
        self.output_vector_combo.setSizePolicy(QSizePolicy.Expanding,
                                               QSizePolicy.Fixed)
        self.output_vector_combo.currentIndexChanged.connect(
            self._on_output_vector_changed)
        rb.addRow("Output Vector:", self.output_vector_combo)
        self._results_section.set_body(results_body)
        outer.addWidget(self._results_section)

        # ------------------------------------------------------------
        # Section 2: DEFORM
        # ------------------------------------------------------------
        # v5.3.0 item 57: collapsed by default. The Results section is
        # the entry point and stays expanded; user expands Deform when
        # they want a deformed shape.
        self._deform_section = CollapsibleSection("Deform", self, expanded=False)
        deform_body = QWidget()
        db = QVBoxLayout(deform_body)
        db.setContentsMargins(0, 4, 0, 4)
        db.setSpacing(6)

        # Style
        style_row = QFormLayout()
        style_row.setContentsMargins(0, 0, 0, 0)
        self.deform_style_combo = QComboBox()
        for label, key in DEFORM_STYLES:
            self.deform_style_combo.addItem(label, key)
        # v5.3.0 item 56: default to Undeformed. Matches the common
        # "None - Model Only" convention so results load without
        # auto-engaging a deformed render. The user explicitly turns
        # on Deformed or Animate when they want to see displacement.
        self.deform_style_combo.setCurrentIndex(0)
        self.deform_style_combo.currentIndexChanged.connect(
            self._on_deform_style_changed)
        style_row.addRow("Style:", self.deform_style_combo)
        # v5.3.2 item 67: removed the per-section Deform Vector
        # dropdown. The deformation always uses the displacement
        # triplet (T1/T2/T3) from the current subcase. The single
        # Output Vector in the Results section is what gets contoured.
        db.addLayout(style_row)

        # Scale group
        scale_group = QGroupBox("Deform Options")
        sg = QFormLayout(scale_group)
        sg.setContentsMargins(6, 4, 6, 4)
        scale_mode_row = QHBoxLayout()
        self._scale_mode_pct = QRadioButton("% of model")
        self._scale_mode_actual = QRadioButton("Actual")
        self._scale_mode_pct.setChecked(True)
        self._scale_mode_group = QButtonGroup(self)
        self._scale_mode_group.addButton(self._scale_mode_pct)
        self._scale_mode_group.addButton(self._scale_mode_actual)
        self._scale_mode_pct.toggled.connect(self._on_scale_changed)
        self._scale_mode_actual.toggled.connect(self._on_scale_changed)
        scale_mode_row.addWidget(self._scale_mode_pct)
        scale_mode_row.addWidget(self._scale_mode_actual)
        scale_mode_row.addStretch(1)
        sg.addRow("Scale mode:", scale_mode_row)
        self.scale_value_spin = QDoubleSpinBox()
        self.scale_value_spin.setDecimals(4)
        self.scale_value_spin.setRange(0.0, 1.0e9)
        self.scale_value_spin.setValue(10.0)
        self.scale_value_spin.setSingleStep(1.0)
        self.scale_value_spin.editingFinished.connect(self._on_scale_changed)
        sg.addRow("Scale value:", self.scale_value_spin)
        db.addWidget(scale_group)

        # Animation options group
        self._anim_group = QGroupBox("Animation Options")
        ag = QFormLayout(self._anim_group)
        ag.setContentsMargins(6, 4, 6, 4)
        anim_mode_row = QHBoxLayout()
        self._anim_mode_sine = QRadioButton("Sine wave")
        self._anim_mode_modes = QRadioButton("Through Modes")
        self._anim_mode_sine.setChecked(True)
        self._anim_mode_group = QButtonGroup(self)
        self._anim_mode_group.addButton(self._anim_mode_sine)
        self._anim_mode_group.addButton(self._anim_mode_modes)
        self._anim_mode_sine.toggled.connect(self._on_anim_settings_changed)
        self._anim_mode_modes.toggled.connect(self._on_anim_settings_changed)
        anim_mode_row.addWidget(self._anim_mode_sine)
        anim_mode_row.addWidget(self._anim_mode_modes)
        anim_mode_row.addStretch(1)
        ag.addRow("Mode:", anim_mode_row)
        self.anim_frames_spin = QSpinBox()
        self.anim_frames_spin.setRange(3, 240)
        self.anim_frames_spin.setValue(20)
        self.anim_frames_spin.editingFinished.connect(
            self._on_anim_settings_changed)
        ag.addRow("Frames:", self.anim_frames_spin)
        self.anim_delay_spin = QSpinBox()
        self.anim_delay_spin.setRange(10, 1000)
        self.anim_delay_spin.setValue(40)
        self.anim_delay_spin.setSuffix(" ms")
        self.anim_delay_spin.editingFinished.connect(
            self._on_anim_settings_changed)
        ag.addRow("Delay:", self.anim_delay_spin)
        db.addWidget(self._anim_group)

        # Undeformed reference checkbox (always visible)
        self._show_undef_ref = QCheckBox("Show undeformed reference")
        self._show_undef_ref.toggled.connect(
            lambda on: self.show_undeformed_ref_changed.emit(bool(on)))
        db.addWidget(self._show_undef_ref)
        db.addStretch(0)

        self._deform_section.set_body(deform_body)
        outer.addWidget(self._deform_section)
        self._update_animation_options_enabled()

        # ------------------------------------------------------------
        # Section 3: CONTOUR
        # ------------------------------------------------------------
        # v5.3.0 item 57: collapsed by default. User expands when they
        # want a gradient overlay.
        self._contour_section = CollapsibleSection("Contour", self, expanded=False)
        contour_body = QWidget()
        cb = QVBoxLayout(contour_body)
        cb.setContentsMargins(0, 4, 0, 4)
        cb.setSpacing(6)

        # Style
        style_form = QFormLayout()
        style_form.setContentsMargins(0, 0, 0, 0)
        self.contour_style_combo = QComboBox()
        for label, key in CONTOUR_STYLES:
            self.contour_style_combo.addItem(label, key)
        self.contour_style_combo.currentIndexChanged.connect(
            lambda *_a: self.contour_style_changed.emit(
                self.contour_style_combo.currentData()))
        style_form.addRow("Style:", self.contour_style_combo)
        # v5.3.2 item 67: removed the per-section Contour Vector
        # dropdown. The Output Vector picked in the Results section
        # drives the contour.
        # v5.3.2 item 69: Bar Position moved to File > Preferences.
        cb.addLayout(style_form)

        # Data Conversion group
        conv_group = QGroupBox("Data Conversion")
        cg = QVBoxLayout(conv_group)
        cg.setContentsMargins(6, 4, 6, 4)
        self._conv_buttons: dict[str, QRadioButton] = {}
        self._conv_group = QButtonGroup(self)
        for label, key in DATA_CONVERSIONS:
            rb_ = QRadioButton(label)
            rb_.setProperty("conv_key", key)
            self._conv_buttons[key] = rb_
            self._conv_group.addButton(rb_)
            cg.addWidget(rb_)
            rb_.toggled.connect(self._on_data_conversion_changed)
        self._conv_buttons['average'].setChecked(True)
        cb.addWidget(conv_group)

        # Levels & Range group
        lr_group = QGroupBox("Levels && Range")
        lg = QFormLayout(lr_group)
        lg.setContentsMargins(6, 4, 6, 4)
        self.levels_spin = QSpinBox()
        self.levels_spin.setRange(2, 32)
        self.levels_spin.setValue(9)
        self.levels_spin.valueChanged.connect(self.levels_changed.emit)
        lg.addRow("Levels:", self.levels_spin)
        self.palette_combo = QComboBox()
        for label, key in PALETTES:
            self.palette_combo.addItem(label, key)
        self.palette_combo.currentIndexChanged.connect(
            lambda *_a: self.palette_changed.emit(
                self.palette_combo.currentData()))
        lg.addRow("Palette:", self.palette_combo)
        self._auto_range = QCheckBox("Auto range")
        self._auto_range.setChecked(True)
        self._auto_range.toggled.connect(self._on_range_changed)
        lg.addRow(self._auto_range)
        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(-1e30, 1e30)
        self._min_spin.setDecimals(6)
        self._min_spin.setEnabled(False)
        self._min_spin.editingFinished.connect(self._on_range_changed)
        lg.addRow("Min:", self._min_spin)
        self._max_spin = QDoubleSpinBox()
        self._max_spin.setRange(-1e30, 1e30)
        self._max_spin.setDecimals(6)
        self._max_spin.setEnabled(False)
        self._max_spin.editingFinished.connect(self._on_range_changed)
        lg.addRow("Max:", self._max_spin)
        self._auto_range.toggled.connect(
            lambda on: self._min_spin.setEnabled(not on))
        self._auto_range.toggled.connect(
            lambda on: self._max_spin.setEnabled(not on))
        cb.addWidget(lr_group)

        # Annotations group
        annot_group = QGroupBox("Annotations")
        ag2 = QFormLayout(annot_group)
        ag2.setContentsMargins(6, 4, 6, 4)
        self._min_marker = QCheckBox("Show Min Marker")
        self._max_marker = QCheckBox("Show Max Marker")
        self._min_marker.toggled.connect(self._on_markers_changed)
        self._max_marker.toggled.connect(self._on_markers_changed)
        ag2.addRow(self._min_marker)
        ag2.addRow(self._max_marker)
        self._labels_check = QCheckBox("Show Element Value Labels")
        self._labels_check.toggled.connect(self._on_labels_changed)
        ag2.addRow(self._labels_check)
        self._labels_top_n = QSpinBox()
        self._labels_top_n.setRange(0, 1_000_000)
        self._labels_top_n.setValue(0)
        self._labels_top_n.setSpecialValueText("All visible elements")
        self._labels_top_n.setSuffix(" (top-N)")
        self._labels_top_n.valueChanged.connect(self._on_labels_changed)
        ag2.addRow("Limit:", self._labels_top_n)
        cb.addWidget(annot_group)
        cb.addStretch(0)

        self._contour_section.set_body(contour_body)
        outer.addWidget(self._contour_section)

        outer.addStretch(1)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def sections(self) -> list[CollapsibleSection]:
        return [self._results_section, self._deform_section,
                self._contour_section]

    def populate_from_results(self, op2_results: Optional[dict],
                              file_label: str = ""):
        """Rebuild the Output Set + Output Vector dropdowns from a
        results bundle (same shape produced by
        :func:`node_runner.model.load_op2_results` and
        :func:`node_runner.solve.mystran_results.load_mystran_results`).

        Pass ``None`` to clear the dropdowns.
        """
        self._loading = True
        try:
            self.output_set_combo.clear()
            self.output_vector_combo.clear()
            self._has_eigenvectors = False
            self._current_subcase_id = None
            if not op2_results or not op2_results.get('subcases'):
                self.output_set_combo.addItem("(no results loaded)", None)
                self.output_set_combo.setEnabled(False)
                self.output_vector_combo.setEnabled(False)
                return
            self.output_set_combo.setEnabled(True)
            self.output_vector_combo.setEnabled(True)
            for sc_id in sorted(op2_results['subcases'].keys()):
                sc = op2_results['subcases'][sc_id]
                title = sc.get('title') or sc.get('subtitle') or ''
                if title:
                    label = f"Subcase {sc_id}: {title}"
                elif sc.get('eigenvectors'):
                    label = f"Subcase {sc_id}: Normal Modes"
                else:
                    label = f"Subcase {sc_id}"
                self.output_set_combo.addItem(label, int(sc_id))
        finally:
            self._loading = False
        # v5.2.0 bug-fix (Round 2): the first addItem above flipped the
        # combo from -1 -> 0 while ``_loading=True``, so the
        # currentIndexChanged handler returned early. After unlock we
        # need to drive the handler manually so the Output Vector combo
        # gets populated and the rest of the toolbox sees a real
        # current_subcase_id. We invoke it explicitly rather than try
        # a -1 -> 0 dance (which would emit an extra "no selection"
        # change first and could confuse MainWindow's handler).
        if self.output_set_combo.count() > 0:
            self._on_output_set_changed(self.output_set_combo.currentIndex())

    def current_subcase_id(self) -> Optional[int]:
        return self._current_subcase_id

    def current_output_vector(self) -> tuple[str, str, int]:
        """Return (kind, component, mode_idx) for the active output
        vector. v5.3.2 item 67: single dropdown lives in the Results
        section now.
        """
        data = self.output_vector_combo.currentData()
        if isinstance(data, tuple) and len(data) == 3:
            return data
        return ('displacement', 'Magnitude', -1)

    def set_data_range(self, vmin: float, vmax: float):
        """Pre-load the Manual min/max spinboxes from the data range."""
        self._min_spin.blockSignals(True)
        self._max_spin.blockSignals(True)
        try:
            self._min_spin.setValue(float(vmin))
            self._max_spin.setValue(float(vmax))
        finally:
            self._min_spin.blockSignals(False)
            self._max_spin.blockSignals(False)

    # ------------------------------------------------------------------
    # Internal handlers
    # ------------------------------------------------------------------

    def _on_output_set_changed(self, _index):
        if self._loading:
            return
        sid = self.output_set_combo.currentData()
        if sid is None:
            return
        self._current_subcase_id = int(sid)
        self._rebuild_output_vector_combo()
        self.output_set_changed.emit(int(sid))

    def _rebuild_output_vector_combo(self):
        """v5.3.2 item 67: single Output Vector combo populated from
        the active subcase. Every renderable scalar shows up here in
        category order:
          1. Displacements (T1/T2/T3/Magnitude, R1/R2/R3/RMag)
          2. SPC Forces
          3. Stress components (von Mises, principals, σxx/σyy/τxy)
          4. Eigenvector modes (one group per mode)
        """
        self._loading = True
        try:
            self.output_vector_combo.clear()
            sid = self._current_subcase_id
            bundle = self._find_bundle_from_parent()
            sc = (bundle.get('subcases', {}).get(sid)
                  if (bundle and sid is not None) else None)
            if sc is None:
                self.output_vector_combo.addItem("(no data)", None)
                return
            self._has_eigenvectors = bool(sc.get('eigenvectors'))
            first = True
            if sc.get('displacements'):
                for label, key in DISP_COMPONENTS:
                    self.output_vector_combo.addItem(
                        f"Displacement - {label}",
                        ('displacement', key, -1))
                first = False
            if sc.get('spc_forces'):
                if not first:
                    self.output_vector_combo.insertSeparator(
                        self.output_vector_combo.count())
                for label, key in SPC_COMPONENTS:
                    self.output_vector_combo.addItem(
                        f"SPC Force - {label}",
                        ('spc_forces', key, -1))
                first = False
            if sc.get('stresses'):
                if not first:
                    self.output_vector_combo.insertSeparator(
                        self.output_vector_combo.count())
                for label, key in STRESS_COMPONENTS:
                    self.output_vector_combo.addItem(
                        f"Stress - {label}",
                        ('stress', key, -1))
                first = False
            eigs = sc.get('eigenvectors') or []
            if eigs:
                freqs = sc.get('frequencies') or []
                for mode_idx in range(len(eigs)):
                    if not first:
                        self.output_vector_combo.insertSeparator(
                            self.output_vector_combo.count())
                    freq = (freqs[mode_idx]
                            if mode_idx < len(freqs)
                            and freqs[mode_idx] is not None else None)
                    suffix = f" ({freq:.3g} Hz)" if freq else ""
                    for label, key in DISP_COMPONENTS:
                        self.output_vector_combo.addItem(
                            f"Mode {mode_idx + 1}{suffix} - {label}",
                            ('eigenvector', key, mode_idx))
                    first = False
            if self.output_vector_combo.count() == 0:
                self.output_vector_combo.addItem("(no data)", None)
        finally:
            self._loading = False
        # Pick the first valid item (Displacement - Magnitude on a
        # SOL 101 deck).
        for i in range(self.output_vector_combo.count()):
            data = self.output_vector_combo.itemData(i)
            if isinstance(data, tuple):
                self.output_vector_combo.setCurrentIndex(i)
                self.output_vector_changed.emit(*data)
                return

    def _on_output_vector_changed(self, _index):
        if self._loading:
            return
        data = self.output_vector_combo.currentData()
        if isinstance(data, tuple) and len(data) == 3:
            self.output_vector_changed.emit(*data)

    def _find_bundle_from_parent(self):
        """Locate the live op2_results dict on MainWindow.

        Avoids storing a reference here; ``populate_from_results``
        is always called with the current bundle, but the Output
        Vector combo also rebuilds whenever the user picks a new
        subcase, so we need to read the bundle on demand.
        """
        p = self.parent()
        while p is not None:
            if hasattr(p, 'op2_results'):
                return getattr(p, 'op2_results', None)
            p = p.parent()
        return None

    # ----- Deform handlers -----

    def _on_deform_style_changed(self, _index):
        style = self.deform_style_combo.currentData()
        self.deform_style_changed.emit(style or 'deformed')
        self._update_animation_options_enabled()

    def _on_scale_changed(self, *_a):
        mode = 'pct' if self._scale_mode_pct.isChecked() else 'actual'
        self.deform_scale_changed.emit(mode, float(self.scale_value_spin.value()))

    def _on_anim_settings_changed(self, *_a):
        mode = 'modes' if self._anim_mode_modes.isChecked() else 'sine'
        self.animation_settings_changed.emit(
            mode,
            int(self.anim_frames_spin.value()),
            int(self.anim_delay_spin.value()))

    def _update_animation_options_enabled(self):
        """v5.2 fix: keep Frames + Delay spinboxes editable regardless
        of the current Deform Style so the user can pre-configure
        animation parameters before picking Animate. Only "Through
        Modes" stays gated on eigenvector data availability.
        """
        # Frames + Delay spinboxes always editable.
        self._anim_group.setEnabled(True)
        if not self._has_eigenvectors and self._anim_mode_modes.isChecked():
            self._anim_mode_sine.setChecked(True)
        self._anim_mode_modes.setEnabled(self._has_eigenvectors)

    # ----- Contour handlers -----

    def _on_data_conversion_changed(self, *_a):
        for key, btn in self._conv_buttons.items():
            if btn.isChecked():
                self.data_conversion_changed.emit(key)
                return

    def _on_range_changed(self, *_a):
        auto = self._auto_range.isChecked()
        self.color_range_changed.emit(
            bool(auto), float(self._min_spin.value()),
            float(self._max_spin.value()))

    def _on_markers_changed(self, *_a):
        self.markers_changed.emit(
            self._min_marker.isChecked(),
            self._max_marker.isChecked())

    def _on_labels_changed(self, *_a):
        self.element_labels_changed.emit(
            self._labels_check.isChecked(),
            int(self._labels_top_n.value()))
