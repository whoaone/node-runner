"""v3.5.0 (item 3): Preferences dialog for entity colors, sizes, and
highlight settings.

Persists via QSettings("NodeRunner", "NodeRunner"). On OK the caller
applies the new values to ``self.type_color_map`` /
``self.highlight_color`` / ``self.mass_glyph_scale`` etc. and rebuilds
the viewer.

Tabs:
  Entity Colors   - color swatch per type
  Entity Sizes    - mass glyph %, node size, beam width, edge width, ...
  Highlight       - highlight color + outline width
"""

from __future__ import annotations

from typing import Optional

from PySide6 import QtCore
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QDoubleSpinBox, QSpinBox, QFrame, QTabWidget,
    QColorDialog, QSizePolicy, QComboBox, QLineEdit, QFileDialog,
    QCheckBox,
)


# ----------------------------------------------------------------------
# Defaults
# ----------------------------------------------------------------------

DEFAULT_TYPE_COLORS = {
    "Shells":     "#0077be",
    "Beams":      "#f85a40",
    "Bars":       "#f85a40",
    "Rods":       "#f85a40",
    "Bushes":     "#1f6feb",
    "Plates":     "#0077be",
    "RBE2":       "#ff3131",
    "RBE3":       "#ffd700",
    "Masses":     "#00cc66",
    "Solids":     "#8b5cf6",
    "Shear":      "#ff8c00",
    "Gap":        "#20b2aa",
    "Plotel":     "#ff69b4",
    "Free Nodes": "#89b4fa",
}

DEFAULT_HIGHLIGHT_COLOR = "#fab387"

# v5.0.0 item 9: brightness multiplier presets used by
# `_install_light_rig` and `_apply_shading_to_actor`.
DEFAULT_BRIGHTNESS_MODE = "Normal"
BRIGHTNESS_MODES = ["Subdued", "Normal", "Bright"]

# v5.0.0 item 7: BDF material library path.
DEFAULT_MATERIALS_LIBRARY_PATH = ""

# v5.0.0 items 12-18: MYSTRAN integration settings live in their own
# mystran/* QSettings prefix (see node_runner.solve.mystran_settings).
# The Preferences tab is a thin UI over that module's get/set helpers.

DEFAULT_SIZES = {
    'mass_glyph_scale_pct': 1.5,   # % of model length; was hardcoded 1.5
    'node_size':            3,
    'beam_width':           2,
    'edge_width':           1,
    'rbe_line_width':       3,
    'free_edge_width':      4,
    'highlight_outline_width': 5,
}


# ----------------------------------------------------------------------
# QSettings load/save helpers (module-level so app startup can read
# them BEFORE the dialog is ever constructed).
# ----------------------------------------------------------------------

def load_preferences() -> dict:
    """Read all preferences from QSettings, falling back to the defaults
    defined above. Returns a dict with keys ``colors``, ``highlight_color``,
    ``sizes``, ``brightness_mode``, ``materials_library_path``."""
    qs = QtCore.QSettings("NodeRunner", "NodeRunner")
    colors = {}
    for key, default in DEFAULT_TYPE_COLORS.items():
        colors[key] = qs.value(f"colors/{key}", default, type=str)
    highlight = qs.value(
        "colors/highlight", DEFAULT_HIGHLIGHT_COLOR, type=str)
    sizes = {}
    for key, default in DEFAULT_SIZES.items():
        cls = float if isinstance(default, float) else int
        try:
            sizes[key] = cls(qs.value(f"sizes/{key}", default, type=cls))
        except Exception:
            sizes[key] = default
    brightness_mode = qs.value(
        "render/brightness_mode", DEFAULT_BRIGHTNESS_MODE, type=str)
    if brightness_mode not in BRIGHTNESS_MODES:
        brightness_mode = DEFAULT_BRIGHTNESS_MODE
    materials_library_path = qs.value(
        "materials/library_path", DEFAULT_MATERIALS_LIBRARY_PATH, type=str)
    try:
        from node_runner.solve import mystran_settings
        mystran = mystran_settings.get_all()
    except Exception:
        mystran = {}
    return {
        'colors': colors,
        'highlight_color': highlight,
        'sizes': sizes,
        'brightness_mode': brightness_mode,
        'materials_library_path': materials_library_path,
        'mystran': mystran,
    }


def save_preferences(payload: dict):
    """Write all preferences to QSettings."""
    qs = QtCore.QSettings("NodeRunner", "NodeRunner")
    for key, value in (payload.get('colors') or {}).items():
        qs.setValue(f"colors/{key}", str(value))
    if 'highlight_color' in payload:
        qs.setValue("colors/highlight", str(payload['highlight_color']))
    for key, value in (payload.get('sizes') or {}).items():
        qs.setValue(f"sizes/{key}", value)
    if 'brightness_mode' in payload:
        qs.setValue("render/brightness_mode", str(payload['brightness_mode']))
    if 'materials_library_path' in payload:
        qs.setValue("materials/library_path",
                    str(payload['materials_library_path']))
    if 'mystran' in payload and isinstance(payload['mystran'], dict):
        try:
            from node_runner.solve import mystran_settings
            mystran_settings.save_all(payload['mystran'])
        except Exception:
            pass
    qs.sync()


# ----------------------------------------------------------------------
# Dialog
# ----------------------------------------------------------------------


class _ColorSwatch(QWidget):
    """Swatch button + hex label + reset, sharing one signal."""

    color_changed = QtCore.Signal(str)

    def __init__(self, initial_hex: str, default_hex: str, parent=None):
        super().__init__(parent)
        self._default = default_hex
        self._current = initial_hex
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        self._swatch = QPushButton()
        self._swatch.setFixedSize(28, 18)
        self._label = QLabel(initial_hex)
        self._label.setStyleSheet(
            "color: #94e2d5; font-family: Consolas, monospace;")
        reset = QPushButton("↺")  # restore arrow
        reset.setFixedWidth(24)
        reset.setToolTip("Reset to default")
        lay.addWidget(self._swatch)
        lay.addWidget(self._label)
        lay.addWidget(reset)
        lay.addStretch(1)
        self._apply_swatch_style()
        self._swatch.clicked.connect(self._pick_color)
        reset.clicked.connect(lambda: self._set_color(self._default))

    def _apply_swatch_style(self):
        self._swatch.setStyleSheet(
            f"background: {self._current}; border: 1px solid #585b70; "
            f"border-radius: 3px;")

    def _pick_color(self):
        col = QColorDialog.getColor(
            QColor(self._current), self, "Pick color")
        if col.isValid():
            self._set_color(col.name())

    def _set_color(self, hex_color: str):
        self._current = hex_color
        self._label.setText(hex_color)
        self._apply_swatch_style()
        self.color_changed.emit(hex_color)

    def current_color(self) -> str:
        return self._current

    def set_color(self, hex_color: str):
        self._set_color(hex_color)


class PreferencesDialog(QDialog):
    """Three-tab preferences dialog for v3.5.0."""

    def __init__(self, current: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setMinimumWidth(480)
        self._initial = current
        self._swatches = {}        # type-name -> _ColorSwatch
        self._size_spins = {}      # size key -> Q*SpinBox
        self._highlight_swatch: Optional[_ColorSwatch] = None
        self._highlight_outline_spin: Optional[QSpinBox] = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(10, 10, 10, 8)
        outer.setSpacing(8)

        title = QLabel("Node Runner Preferences")
        title.setStyleSheet(
            "font-weight: bold; font-size: 13px; color: #f9e2af;")
        outer.addWidget(title)

        self._brightness_combo: Optional[QComboBox] = None
        self._library_path_edit: Optional[QLineEdit] = None
        self._library_status_label: Optional[QLabel] = None
        # MYSTRAN tab widgets (v5.0.0 items 12-18)
        self._mystran_exe_edit: Optional[QLineEdit] = None
        self._mystran_version_label: Optional[QLabel] = None
        self._mystran_scratch_edit: Optional[QLineEdit] = None
        self._mystran_keep_check: Optional[QCheckBox] = None
        self._mystran_sollib_combo: Optional[QComboBox] = None
        self._mystran_quad4typ_combo: Optional[QComboBox] = None
        self._mystran_wtmass_spin: Optional[QDoubleSpinBox] = None
        self._mystran_autoload_check: Optional[QCheckBox] = None
        self._mystran_autoopen_check: Optional[QCheckBox] = None
        self._mystran_cleanup_spin: Optional[QSpinBox] = None

        self._tabs = QTabWidget()
        self._tabs.addTab(self._build_colors_tab(), "Entity Colors")
        self._tabs.addTab(self._build_sizes_tab(), "Entity Sizes")
        self._tabs.addTab(self._build_rendering_tab(), "Rendering")
        self._tabs.addTab(self._build_library_tab(), "Library")
        self._tabs.addTab(self._build_mystran_tab(), "MYSTRAN")
        self._tabs.addTab(self._build_highlight_tab(), "Highlight")
        outer.addWidget(self._tabs)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        outer.addWidget(sep)

        bot = QHBoxLayout()
        restore = QPushButton("Restore Defaults")
        restore.setToolTip("Reset every setting on every tab to its default.")
        restore.clicked.connect(self._restore_all_defaults)
        bot.addWidget(restore)
        bot.addStretch(1)
        ok = QPushButton("OK")
        ok.setStyleSheet("font-weight: bold;")
        ok.setMinimumWidth(90)
        ok.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.setMinimumWidth(80)
        cancel.clicked.connect(self.reject)
        bot.addWidget(ok)
        bot.addWidget(cancel)
        outer.addLayout(bot)

    # ------------------------------------------------------------------
    # Tabs
    # ------------------------------------------------------------------

    def _build_colors_tab(self):
        w = QWidget()
        lay = QFormLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)
        cur_colors = self._initial.get('colors') or DEFAULT_TYPE_COLORS
        for key, default in DEFAULT_TYPE_COLORS.items():
            cur = cur_colors.get(key, default)
            swatch = _ColorSwatch(cur, default)
            self._swatches[key] = swatch
            lay.addRow(QLabel(f"{key}:"), swatch)
        return w

    def _build_sizes_tab(self):
        w = QWidget()
        lay = QFormLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)
        cur_sizes = self._initial.get('sizes') or DEFAULT_SIZES
        rows = [
            ('mass_glyph_scale_pct', 'Mass glyph scale (% of model length):',
             0.1, 20.0, 0.1, 1),
            ('node_size', 'Node point size (px):', 1, 30, 1, 0),
            ('beam_width', 'Beam line width (px):', 1, 10, 1, 0),
            ('edge_width', 'Edge width (px):', 1, 8, 1, 0),
            ('rbe_line_width', 'RBE line width (px):', 1, 10, 1, 0),
            ('free_edge_width', 'Free-edge line width (px):', 1, 10, 1, 0),
        ]
        for key, label, lo, hi, step, decimals in rows:
            if decimals:
                sb = QDoubleSpinBox()
                sb.setDecimals(decimals)
            else:
                sb = QSpinBox()
            sb.setRange(lo, hi)
            sb.setSingleStep(step)
            sb.setValue(cur_sizes.get(key, DEFAULT_SIZES[key]))
            sb.setMinimumWidth(110)
            self._size_spins[key] = sb
            lay.addRow(QLabel(label), sb)
        return w

    def _build_rendering_tab(self):
        """v5.0.0 item 9: brightness mode selector."""
        w = QWidget()
        lay = QFormLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)

        cur_mode = self._initial.get('brightness_mode', DEFAULT_BRIGHTNESS_MODE)
        if cur_mode not in BRIGHTNESS_MODES:
            cur_mode = DEFAULT_BRIGHTNESS_MODE
        self._brightness_combo = QComboBox()
        self._brightness_combo.addItems(BRIGHTNESS_MODES)
        self._brightness_combo.setCurrentText(cur_mode)
        self._brightness_combo.setMinimumWidth(140)
        lay.addRow(QLabel("Mesh brightness:"), self._brightness_combo)

        note = QLabel(
            "Multiplies the scene lights and shading ambient. "
            "Subdued = 0.85×, Normal = 1.0×, Bright = 1.15×.\n"
            "Takes effect on the next viewer refresh."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: #888; font-style: italic; margin-top: 8px;")
        lay.addRow(note)
        return w

    def _build_library_tab(self):
        """v5.0.0 item 7: materials library file path."""
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)

        intro = QLabel(
            "Path to your team's MAT-card library (BDF). Materials in this "
            "file are auto-loaded when Node Runner starts. Leave empty to "
            "skip auto-loading.\n\n"
            "For built-in MMPDS preset values, use Materials → Load… in the "
            "Create Material dialog instead."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        lay.addWidget(intro)

        row = QHBoxLayout()
        cur_path = self._initial.get(
            'materials_library_path', DEFAULT_MATERIALS_LIBRARY_PATH)
        self._library_path_edit = QLineEdit(cur_path)
        self._library_path_edit.setPlaceholderText(
            "C:\\path\\to\\materials.bdf")
        self._library_path_edit.textChanged.connect(self._refresh_library_status)
        row.addWidget(self._library_path_edit, 1)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._pick_library_path)
        row.addWidget(browse)
        clear = QPushButton("Clear")
        clear.clicked.connect(lambda: self._library_path_edit.setText(""))
        row.addWidget(clear)
        lay.addLayout(row)

        self._library_status_label = QLabel("")
        self._library_status_label.setStyleSheet(
            "font-family: Consolas, monospace; font-size: 11px;")
        lay.addWidget(self._library_status_label)
        self._refresh_library_status()
        lay.addStretch(1)
        return w

    def _pick_library_path(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select materials library (BDF)",
            self._library_path_edit.text() or "",
            "Nastran decks (*.bdf *.dat *.nas);;All files (*.*)"
        )
        if path:
            self._library_path_edit.setText(path)

    def _refresh_library_status(self):
        if self._library_status_label is None:
            return
        path = self._library_path_edit.text().strip() if self._library_path_edit else ""
        if not path:
            self._library_status_label.setText(
                "<i>No library configured</i>")
            self._library_status_label.setStyleSheet(
                "color: #6c7086; font-family: Consolas, monospace; font-size: 11px;")
            return
        import os
        if os.path.isfile(path):
            self._library_status_label.setText(
                f"<b>✓</b> File exists: <code>{path}</code>")
            self._library_status_label.setStyleSheet(
                "color: #a6e3a1; font-family: Consolas, monospace; font-size: 11px;")
        else:
            self._library_status_label.setText(
                f"<b>✗</b> File not found: <code>{path}</code>")
            self._library_status_label.setStyleSheet(
                "color: #f38ba8; font-family: Consolas, monospace; font-size: 11px;")

    def _build_mystran_tab(self):
        """v5.0.0 items 12-18: MYSTRAN solver integration settings.

        Reads / writes via :mod:`node_runner.solve.mystran_settings` so
        the solve package owns the key prefix (``mystran/*``).
        Tooltips throughout assume the user knows MSC Nastran or Femap
        but is new to MYSTRAN -- they explicitly call out where the
        MYSTRAN dialect diverges.
        """
        try:
            from node_runner.solve import mystran_settings as _ms
            from node_runner.solve import mystran_tips as _tips
            current = _ms.get_all()
        except Exception:
            from node_runner.solve.mystran_settings import DEFAULTS
            current = dict(DEFAULTS)
            _tips = None

        def _tip(attr_name):
            return getattr(_tips, attr_name, "") if _tips else ""

        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)

        intro = QLabel(
            "Configure the open-source MYSTRAN solver. Node Runner drives "
            "the binary as a subprocess; we do NOT bundle MYSTRAN itself.\n"
            "Supported solutions: SOL 1 (linear static), SOL 3 (normal modes), "
            "SOL 5 (linear buckling)."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #cdd6f4; font-size: 11px;")
        lay.addWidget(intro)

        form = QFormLayout()
        form.setSpacing(6)

        # --- Executable path + Detect ---
        exe_row = QHBoxLayout()
        self._mystran_exe_edit = QLineEdit(current.get("executable_path", ""))
        self._mystran_exe_edit.setPlaceholderText(
            "C:\\MYSTRAN\\mystran.exe or /usr/local/bin/mystran")
        self._mystran_exe_edit.setToolTip(_tip("EXE_PATH"))
        self._mystran_exe_edit.textChanged.connect(
            self._refresh_mystran_version)
        exe_row.addWidget(self._mystran_exe_edit, 1)
        browse_btn = QPushButton("Browse…")
        browse_btn.setToolTip("Pick the MYSTRAN binary from disk.")
        browse_btn.clicked.connect(self._pick_mystran_exe)
        exe_row.addWidget(browse_btn)
        detect_btn = QPushButton("Detect")
        detect_btn.setToolTip(
            "Search PATH and common install locations "
            "(C:\\MYSTRAN\\, /usr/local/bin/mystran, MYSTRAN_EXE env var).")
        detect_btn.clicked.connect(self._auto_detect_mystran)
        exe_row.addWidget(detect_btn)
        exe_label = QLabel("MYSTRAN executable:")
        exe_label.setToolTip(_tip("EXE_PATH"))
        form.addRow(exe_label, exe_row)

        self._mystran_version_label = QLabel("")
        self._mystran_version_label.setStyleSheet(
            "font-family: Consolas, monospace; font-size: 11px;")
        form.addRow(QLabel(""), self._mystran_version_label)
        self._refresh_mystran_version()

        # --- Scratch dir ---
        scratch_row = QHBoxLayout()
        self._mystran_scratch_edit = QLineEdit(current.get("scratch_root", ""))
        self._mystran_scratch_edit.setPlaceholderText(
            "Where MYSTRAN run folders go")
        self._mystran_scratch_edit.setToolTip(_tip("SCRATCH_ROOT"))
        scratch_row.addWidget(self._mystran_scratch_edit, 1)
        scratch_browse = QPushButton("Browse…")
        scratch_browse.clicked.connect(self._pick_mystran_scratch)
        scratch_row.addWidget(scratch_browse)
        scratch_label = QLabel("Scratch directory:")
        scratch_label.setToolTip(_tip("SCRATCH_ROOT"))
        form.addRow(scratch_label, scratch_row)

        # --- Solver defaults ---
        self._mystran_sollib_combo = QComboBox()
        # MYSTRAN v18+ accepts SPARSE or BANDED for SOLLIB; IntMKL /
        # LAPACK names are not valid (MYSTRAN warns and falls back to
        # SPARSE).
        self._mystran_sollib_combo.addItems(["SPARSE", "BANDED"])
        cur_sollib = str(current.get("default_sollib", "SPARSE")).upper()
        if cur_sollib not in ("SPARSE", "BANDED"):
            cur_sollib = "SPARSE"
        self._mystran_sollib_combo.setCurrentText(cur_sollib)
        self._mystran_sollib_combo.setToolTip(_tip("SOLLIB"))
        sollib_label = QLabel("Default SOLLIB:")
        sollib_label.setToolTip(_tip("SOLLIB"))
        form.addRow(sollib_label, self._mystran_sollib_combo)

        self._mystran_quad4typ_combo = QComboBox()
        self._mystran_quad4typ_combo.addItems(["MIN4T", "MIN4"])
        self._mystran_quad4typ_combo.setCurrentText(
            str(current.get("default_quad4typ", "MIN4T")))
        self._mystran_quad4typ_combo.setToolTip(_tip("QUAD4TYP"))
        quad4_label = QLabel("Default QUAD4TYP:")
        quad4_label.setToolTip(_tip("QUAD4TYP"))
        form.addRow(quad4_label, self._mystran_quad4typ_combo)

        self._mystran_wtmass_spin = QDoubleSpinBox()
        self._mystran_wtmass_spin.setDecimals(6)
        self._mystran_wtmass_spin.setRange(1e-12, 1e6)
        self._mystran_wtmass_spin.setSingleStep(0.01)
        self._mystran_wtmass_spin.setValue(
            float(current.get("default_wtmass", 1.0)))
        self._mystran_wtmass_spin.setToolTip(_tip("WTMASS"))
        wtmass_label = QLabel("Default WTMASS:")
        wtmass_label.setToolTip(_tip("WTMASS"))
        form.addRow(wtmass_label, self._mystran_wtmass_spin)

        # --- Behaviour ---
        self._mystran_keep_check = QCheckBox(
            "Keep intermediate run files (BDF / F06 / OP2 / log)")
        self._mystran_keep_check.setChecked(
            bool(current.get("keep_intermediates", True)))
        self._mystran_keep_check.setToolTip(_tip("KEEP_INTERMEDIATES"))
        form.addRow(QLabel(""), self._mystran_keep_check)

        self._mystran_autoload_check = QCheckBox(
            "Auto-load results after a successful run")
        self._mystran_autoload_check.setChecked(
            bool(current.get("auto_load_results", True)))
        self._mystran_autoload_check.setToolTip(_tip("AUTO_LOAD_RESULTS"))
        form.addRow(QLabel(""), self._mystran_autoload_check)

        self._mystran_autoopen_check = QCheckBox(
            "Auto-open the Result Browser dock when results load")
        self._mystran_autoopen_check.setChecked(
            bool(current.get("auto_open_browser", True)))
        self._mystran_autoopen_check.setToolTip(_tip("AUTO_OPEN_BROWSER"))
        form.addRow(QLabel(""), self._mystran_autoopen_check)

        self._mystran_cleanup_spin = QSpinBox()
        self._mystran_cleanup_spin.setRange(0, 365)
        self._mystran_cleanup_spin.setValue(
            int(current.get("cleanup_after_days", 30)))
        self._mystran_cleanup_spin.setSuffix(" days  (0 = never)")
        self._mystran_cleanup_spin.setToolTip(_tip("CLEANUP_AFTER_DAYS"))
        cleanup_label = QLabel("Auto-cleanup runs older than:")
        cleanup_label.setToolTip(_tip("CLEANUP_AFTER_DAYS"))
        form.addRow(cleanup_label, self._mystran_cleanup_spin)

        lay.addLayout(form)
        lay.addStretch(1)
        return w

    def _pick_mystran_exe(self):
        from PySide6.QtCore import QStandardPaths
        start_dir = (self._mystran_exe_edit.text()
                     or QStandardPaths.writableLocation(
                         QStandardPaths.HomeLocation))
        path, _ = QFileDialog.getOpenFileName(
            self, "Select MYSTRAN executable", start_dir,
            "Executables (*.exe);;All files (*.*)"
        )
        if path:
            self._mystran_exe_edit.setText(path)

    def _pick_mystran_scratch(self):
        from PySide6.QtCore import QStandardPaths
        start_dir = (self._mystran_scratch_edit.text()
                     or QStandardPaths.writableLocation(
                         QStandardPaths.HomeLocation))
        path = QFileDialog.getExistingDirectory(
            self, "Select MYSTRAN scratch directory", start_dir)
        if path:
            self._mystran_scratch_edit.setText(path)

    def _auto_detect_mystran(self):
        try:
            from node_runner.solve.mystran_runner import (
                discover_mystran_executable)
            found = discover_mystran_executable()
        except Exception:
            found = None
        if found:
            self._mystran_exe_edit.setText(str(found))
            self._refresh_mystran_version()
        else:
            if self._mystran_version_label is not None:
                self._mystran_version_label.setText(
                    "<i>Not found on PATH or in common install locations</i>")
                self._mystran_version_label.setStyleSheet(
                    "color: #f38ba8; font-family: Consolas, monospace; "
                    "font-size: 11px;")

    def _refresh_mystran_version(self):
        if self._mystran_version_label is None:
            return
        path = self._mystran_exe_edit.text().strip()
        if not path:
            self._mystran_version_label.setText(
                "<i>Detect or browse to your MYSTRAN binary</i>")
            self._mystran_version_label.setStyleSheet(
                "color: #6c7086; font-family: Consolas, monospace; "
                "font-size: 11px;")
            return
        import os
        if not os.path.isfile(path):
            self._mystran_version_label.setText(
                "<b>✗</b> File not found")
            self._mystran_version_label.setStyleSheet(
                "color: #f38ba8; font-family: Consolas, monospace; "
                "font-size: 11px;")
            return
        try:
            from node_runner.solve.mystran_runner import detect_mystran_version
            ver = detect_mystran_version(path)
        except Exception as e:
            ver = f"(version probe failed: {e})"
        self._mystran_version_label.setText(
            f"<b>✓</b> {ver}" if ver else "<b>✓</b> Executable found")
        self._mystran_version_label.setStyleSheet(
            "color: #a6e3a1; font-family: Consolas, monospace; "
            "font-size: 11px;")

    def _build_highlight_tab(self):
        w = QWidget()
        lay = QFormLayout(w)
        lay.setSpacing(6)
        lay.setContentsMargins(10, 10, 10, 10)
        cur_highlight = self._initial.get(
            'highlight_color', DEFAULT_HIGHLIGHT_COLOR)
        self._highlight_swatch = _ColorSwatch(
            cur_highlight, DEFAULT_HIGHLIGHT_COLOR)
        lay.addRow(QLabel("Highlight color:"), self._highlight_swatch)
        cur_outline = (self._initial.get('sizes') or {}).get(
            'highlight_outline_width',
            DEFAULT_SIZES['highlight_outline_width'])
        self._highlight_outline_spin = QSpinBox()
        self._highlight_outline_spin.setRange(1, 12)
        self._highlight_outline_spin.setValue(cur_outline)
        self._highlight_outline_spin.setMinimumWidth(110)
        lay.addRow(QLabel("Selection outline width (px):"),
                   self._highlight_outline_spin)
        note = QLabel(
            "Highlight color is shown when entities are selected in\n"
            "the viewport. Outline width affects the wireframe overlay.")
        note.setStyleSheet("color: #888; font-style: italic; margin-top: 8px;")
        lay.addRow(note)
        return w

    # ------------------------------------------------------------------
    # Behavior
    # ------------------------------------------------------------------

    def _restore_all_defaults(self):
        for key, default in DEFAULT_TYPE_COLORS.items():
            if key in self._swatches:
                self._swatches[key].set_color(default)
        if self._highlight_swatch is not None:
            self._highlight_swatch.set_color(DEFAULT_HIGHLIGHT_COLOR)
        for key, default in DEFAULT_SIZES.items():
            sb = self._size_spins.get(key)
            if sb is not None:
                sb.setValue(default)
        if self._highlight_outline_spin is not None:
            self._highlight_outline_spin.setValue(
                DEFAULT_SIZES['highlight_outline_width'])

    # ------------------------------------------------------------------
    # Result
    # ------------------------------------------------------------------

    def result_payload(self) -> dict:
        colors = {k: s.current_color() for k, s in self._swatches.items()}
        sizes = {k: sb.value() for k, sb in self._size_spins.items()}
        if self._highlight_outline_spin is not None:
            sizes['highlight_outline_width'] = (
                self._highlight_outline_spin.value())
        highlight = (self._highlight_swatch.current_color()
                     if self._highlight_swatch else DEFAULT_HIGHLIGHT_COLOR)
        brightness_mode = (self._brightness_combo.currentText()
                           if self._brightness_combo is not None
                           else DEFAULT_BRIGHTNESS_MODE)
        materials_library_path = (
            self._library_path_edit.text().strip()
            if self._library_path_edit is not None
            else DEFAULT_MATERIALS_LIBRARY_PATH
        )
        mystran = {}
        if self._mystran_exe_edit is not None:
            mystran['executable_path'] = self._mystran_exe_edit.text().strip()
        if self._mystran_scratch_edit is not None:
            mystran['scratch_root'] = self._mystran_scratch_edit.text().strip()
        if self._mystran_keep_check is not None:
            mystran['keep_intermediates'] = bool(
                self._mystran_keep_check.isChecked())
        if self._mystran_sollib_combo is not None:
            mystran['default_sollib'] = self._mystran_sollib_combo.currentText()
        if self._mystran_quad4typ_combo is not None:
            mystran['default_quad4typ'] = (
                self._mystran_quad4typ_combo.currentText())
        if self._mystran_wtmass_spin is not None:
            mystran['default_wtmass'] = float(
                self._mystran_wtmass_spin.value())
        if self._mystran_autoload_check is not None:
            mystran['auto_load_results'] = bool(
                self._mystran_autoload_check.isChecked())
        if self._mystran_autoopen_check is not None:
            mystran['auto_open_browser'] = bool(
                self._mystran_autoopen_check.isChecked())
        if self._mystran_cleanup_spin is not None:
            mystran['cleanup_after_days'] = int(
                self._mystran_cleanup_spin.value())
        return {
            'colors': colors,
            'highlight_color': highlight,
            'sizes': sizes,
            'brightness_mode': brightness_mode,
            'materials_library_path': materials_library_path,
            'mystran': mystran,
        }
