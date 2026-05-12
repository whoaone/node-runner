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
    QColorDialog, QSizePolicy,
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
    """Read all v3.5.0 preferences from QSettings, falling back to the
    defaults defined above. Returns a single dict with keys
    ``colors``, ``highlight_color``, ``sizes``."""
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
    return {'colors': colors, 'highlight_color': highlight, 'sizes': sizes}


def save_preferences(payload: dict):
    """Write all v3.5.0 preferences to QSettings."""
    qs = QtCore.QSettings("NodeRunner", "NodeRunner")
    for key, value in (payload.get('colors') or {}).items():
        qs.setValue(f"colors/{key}", str(value))
    if 'highlight_color' in payload:
        qs.setValue("colors/highlight", str(payload['highlight_color']))
    for key, value in (payload.get('sizes') or {}).items():
        qs.setValue(f"sizes/{key}", value)
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

        self._tabs = QTabWidget()
        self._tabs.addTab(self._build_colors_tab(), "Entity Colors")
        self._tabs.addTab(self._build_sizes_tab(), "Entity Sizes")
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
        return {
            'colors': colors,
            'highlight_color': highlight,
            'sizes': sizes,
        }
