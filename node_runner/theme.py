"""Node Runner theme system - premium dark and light palettes with QSS."""

from PySide6 import QtCore
from PySide6.QtGui import QPalette, QColor

# ---------------------------------------------------------------------------
# Color constants  (Catppuccin Mocha inspired)
# ---------------------------------------------------------------------------
BG_PRIMARY = "#1e1e2e"
BG_SECONDARY = "#181825"
SURFACE = "#313244"
BORDER = "#45475a"
TEXT = "#cdd6f4"
TEXT_SECONDARY = "#a6adc8"
ACCENT = "#89b4fa"
ACCENT_HOVER = "#b4d0fb"
SELECTION = "#45475a"
SUCCESS = "#a6e3a1"
WARNING = "#f9e2af"
ERROR = "#f38ba8"
SCROLLBAR_HANDLE = "#585b70"
SCROLLBAR_BG = "#181825"
# Phase 3: selection accent. Catppuccin Peach - contrasts with ACCENT (blue)
# without straying from the existing palette.
SELECTION_ACCENT = "#fab387"

# ---------------------------------------------------------------------------
# Dark palette
# ---------------------------------------------------------------------------
dark_palette = QPalette()
dark_palette.setColor(QPalette.Window, QColor(BG_PRIMARY))
dark_palette.setColor(QPalette.WindowText, QColor(TEXT))
dark_palette.setColor(QPalette.Base, QColor(BG_SECONDARY))
dark_palette.setColor(QPalette.AlternateBase, QColor(SURFACE))
dark_palette.setColor(QPalette.Text, QColor(TEXT))
dark_palette.setColor(QPalette.Button, QColor(SURFACE))
dark_palette.setColor(QPalette.ButtonText, QColor(TEXT))
dark_palette.setColor(QPalette.Highlight, QColor(ACCENT))
dark_palette.setColor(QPalette.HighlightedText, QColor(BG_PRIMARY))
dark_palette.setColor(QPalette.ToolTipBase, QColor(SURFACE))
dark_palette.setColor(QPalette.ToolTipText, QColor(TEXT))
dark_palette.setColor(QPalette.PlaceholderText, QColor(TEXT_SECONDARY))
dark_palette.setColor(QPalette.BrightText, QColor(ERROR))
dark_palette.setColor(QPalette.Link, QColor(ACCENT))
dark_palette.setColor(QPalette.LinkVisited, QColor("#cba6f7"))
# Disabled group
dark_palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor("#6c7086"))
dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor("#6c7086"))
dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor("#6c7086"))
dark_palette.setColor(QPalette.Disabled, QPalette.Highlight, QColor("#45475a"))
dark_palette.setColor(QPalette.Disabled, QPalette.HighlightedText, QColor("#6c7086"))

# ---------------------------------------------------------------------------
# Light palette  (kept minimal - system defaults + overrides)
# ---------------------------------------------------------------------------
light_palette = QPalette()

# ---------------------------------------------------------------------------
# QSS Stylesheet
# ---------------------------------------------------------------------------
DARK_STYLESHEET = f"""
/* ===== Global ===== */
QMainWindow {{
    background-color: {BG_PRIMARY};
}}
QWidget {{
    background-color: {BG_PRIMARY};
    color: {TEXT};
    font-family: "Segoe UI", "Helvetica Neue", sans-serif;
    font-size: 13px;
}}

/* ===== Menu bar ===== */
QMenuBar {{
    background-color: {BG_SECONDARY};
    color: {TEXT};
    border-bottom: 1px solid {SURFACE};
    padding: 2px 0px;
    spacing: 0px;
}}
QMenuBar::item {{
    padding: 5px 10px;
    border-radius: 4px;
    margin: 2px 1px;
}}
QMenuBar::item:selected {{
    background-color: {SURFACE};
}}

/* ===== Menus ===== */
QMenu {{
    background-color: {BG_PRIMARY};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 4px 0px;
}}
QMenu::item {{
    padding: 6px 28px 6px 14px;
    border-radius: 3px;
    margin: 1px 4px;
}}
QMenu::item:selected {{
    background-color: {SURFACE};
    color: {TEXT};
}}
QMenu::separator {{
    height: 1px;
    background: {BORDER};
    margin: 4px 10px;
}}
QMenu::indicator {{
    width: 14px;
    height: 14px;
    margin-left: 6px;
}}

/* ===== Tabs ===== */
QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-top: none;
    background: {BG_PRIMARY};
}}
QTabBar::tab {{
    background: {BG_SECONDARY};
    color: {TEXT_SECONDARY};
    padding: 7px 16px;
    border: 1px solid {BORDER};
    border-bottom: none;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    margin-right: 2px;
}}
QTabBar::tab:selected {{
    background: {BG_PRIMARY};
    color: {TEXT};
    border-bottom: 2px solid {ACCENT};
}}
QTabBar::tab:hover:!selected {{
    background: {SURFACE};
    color: {TEXT};
}}

/* ===== Tree widgets ===== */
QTreeWidget {{
    background-color: {BG_SECONDARY};
    border: 1px solid {BORDER};
    border-radius: 4px;
    outline: none;
    padding: 2px;
}}
QTreeWidget::item {{
    padding: 3px 4px;
    border-radius: 3px;
}}
QTreeWidget::item:hover {{
    background-color: {SURFACE};
}}
QTreeWidget::item:selected {{
    background-color: {SELECTION};
    color: {TEXT};
}}
QTreeWidget QHeaderView::section {{
    background: {BG_SECONDARY};
    color: {TEXT_SECONDARY};
    border: none;
    border-bottom: 1px solid {BORDER};
    padding: 5px 8px;
    font-weight: 600;
}}

/* ===== Buttons ===== */
QPushButton {{
    background-color: {SURFACE};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 6px 14px;
    min-height: 14px;
}}
QPushButton:hover {{
    background-color: {BORDER};
    border-color: {ACCENT};
}}
QPushButton:pressed {{
    background-color: {BG_SECONDARY};
}}
QPushButton:disabled {{
    color: #6c7086;
    border-color: {SURFACE};
}}
QPushButton:flat {{
    background-color: transparent;
    border: none;
    color: {ACCENT};
    padding: 4px 10px;
}}
QPushButton:flat:hover {{
    color: {ACCENT_HOVER};
    background-color: rgba(137, 180, 250, 0.08);
    border-radius: 4px;
}}

/* ===== Group boxes ===== */
QGroupBox {{
    border: 1px solid {BORDER};
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 18px;
    font-weight: 600;
    color: {TEXT};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 2px 8px;
    color: {ACCENT};
}}

/* ===== Checkboxes ===== */
QCheckBox {{
    spacing: 8px;
    color: {TEXT};
}}
QCheckBox::indicator {{
    width: 16px;
    height: 16px;
    border: 2px solid {BORDER};
    border-radius: 4px;
    background: {BG_SECONDARY};
}}
QCheckBox::indicator:hover {{
    border-color: {ACCENT};
}}
QCheckBox::indicator:checked {{
    background-color: {ACCENT};
    border-color: {ACCENT};
}}
QCheckBox::indicator:checked:hover {{
    background-color: {ACCENT_HOVER};
    border-color: {ACCENT_HOVER};
}}
QCheckBox::indicator:indeterminate {{
    background-color: {BORDER};
    border-color: {BORDER};
}}

/* ===== Radio buttons ===== */
QRadioButton {{
    spacing: 8px;
    color: {TEXT};
}}
QRadioButton::indicator {{
    width: 16px;
    height: 16px;
    border: 2px solid {BORDER};
    border-radius: 9px;
    background: {BG_SECONDARY};
}}
QRadioButton::indicator:hover {{
    border-color: {ACCENT};
}}
QRadioButton::indicator:checked {{
    background-color: {ACCENT};
    border-color: {ACCENT};
}}

/* ===== Combo boxes ===== */
QComboBox {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 5px;
    padding: 5px 10px;
    color: {TEXT};
    min-height: 16px;
}}
QComboBox:hover {{
    border-color: {ACCENT};
}}
QComboBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{
    border: none;
    width: 24px;
    border-left: 1px solid {BORDER};
}}
QComboBox QAbstractItemView {{
    background-color: {BG_PRIMARY};
    border: 1px solid {BORDER};
    border-radius: 4px;
    selection-background-color: {SURFACE};
    selection-color: {TEXT};
    padding: 2px;
    outline: none;
}}

/* ===== Line edits ===== */
QLineEdit {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 5px;
    padding: 5px 8px;
    color: {TEXT};
    selection-background-color: {ACCENT};
    selection-color: {BG_PRIMARY};
}}
QLineEdit:focus {{
    border-color: {ACCENT};
}}
QLineEdit:disabled {{
    color: #6c7086;
    background-color: {BG_SECONDARY};
}}

/* ===== Spin boxes ===== */
QSpinBox, QDoubleSpinBox {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 5px;
    padding: 4px 8px;
    color: {TEXT};
}}
QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {ACCENT};
}}

/* ===== Scrollbars  (VS Code thin style) ===== */
QScrollBar:vertical {{
    background: {SCROLLBAR_BG};
    width: 10px;
    border: none;
    border-radius: 5px;
}}
QScrollBar::handle:vertical {{
    background: {BORDER};
    min-height: 24px;
    border-radius: 5px;
    margin: 2px;
}}
QScrollBar::handle:vertical:hover {{
    background: {SCROLLBAR_HANDLE};
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0px;
}}
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {{
    background: transparent;
}}
QScrollBar:horizontal {{
    background: {SCROLLBAR_BG};
    height: 10px;
    border: none;
    border-radius: 5px;
}}
QScrollBar::handle:horizontal {{
    background: {BORDER};
    min-width: 24px;
    border-radius: 5px;
    margin: 2px;
}}
QScrollBar::handle:horizontal:hover {{
    background: {SCROLLBAR_HANDLE};
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0px;
}}
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal {{
    background: transparent;
}}

/* ===== Splitter ===== */
QSplitter::handle {{
    background-color: {SURFACE};
}}
QSplitter::handle:horizontal {{
    width: 3px;
}}
QSplitter::handle:vertical {{
    height: 3px;
}}
QSplitter::handle:hover {{
    background-color: {ACCENT};
}}

/* ===== Status bar ===== */
QStatusBar {{
    background-color: {BG_SECONDARY};
    border-top: 1px solid {SURFACE};
    color: {TEXT_SECONDARY};
    padding: 2px 8px;
}}

/* ===== Tooltips ===== */
QToolTip {{
    background-color: {SURFACE};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 4px;
    padding: 5px 8px;
}}

/* ===== Dialogs ===== */
QDialog {{
    background-color: {BG_PRIMARY};
}}

/* ===== Tables ===== */
QTableWidget {{
    background-color: {BG_SECONDARY};
    gridline-color: {BORDER};
    border: 1px solid {BORDER};
    border-radius: 4px;
    outline: none;
}}
QTableWidget::item {{
    padding: 4px 6px;
    color: {TEXT};
}}
QTableWidget::item:selected {{
    background-color: {SELECTION};
    color: {TEXT};
}}
QHeaderView::section {{
    background-color: {SURFACE};
    color: {TEXT};
    border: none;
    border-right: 1px solid {BORDER};
    border-bottom: 1px solid {BORDER};
    padding: 5px 8px;
    font-weight: 600;
}}

/* ===== Text edit ===== */
QTextEdit {{
    background-color: {BG_SECONDARY};
    border: 1px solid {BORDER};
    border-radius: 4px;
    color: {TEXT};
    padding: 4px;
    selection-background-color: {ACCENT};
    selection-color: {BG_PRIMARY};
}}

/* ===== List widgets ===== */
QListWidget {{
    background-color: {BG_SECONDARY};
    border: 1px solid {BORDER};
    border-radius: 4px;
    outline: none;
    padding: 2px;
}}
QListWidget::item {{
    padding: 4px 6px;
    border-radius: 3px;
}}
QListWidget::item:hover {{
    background-color: {SURFACE};
}}
QListWidget::item:selected {{
    background-color: {SELECTION};
    color: {TEXT};
}}

/* ===== Scroll area ===== */
QScrollArea {{
    border: none;
    background-color: {BG_PRIMARY};
}}

/* ===== Labels (form labels inherit from QWidget) ===== */
QLabel {{
    background-color: transparent;
}}

/* ===== Dialog button box ===== */
QDialogButtonBox QPushButton {{
    min-width: 72px;
}}

/* ===== Input dialog ===== */
QInputDialog {{
    background-color: {BG_PRIMARY};
}}

/* ===== Message box ===== */
QMessageBox {{
    background-color: {BG_PRIMARY};
}}

/* ===== Progress bar ===== */
QProgressBar {{
    background-color: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 4px;
    text-align: center;
    color: {TEXT};
    height: 18px;
}}
QProgressBar::chunk {{
    background-color: {ACCENT};
    border-radius: 3px;
}}

/* ===== Frame ===== */
QFrame {{
    background-color: transparent;
}}
"""

# ---------------------------------------------------------------------------
# Background presets for Display Settings
# ---------------------------------------------------------------------------
BACKGROUND_PRESETS = {
    "Dark (Default)": {"mode": "solid", "color": "#1e1e2e"},
    "Light": {"mode": "solid", "color": "#ffffff"},
    "Blue Gradient": {"mode": "gradient", "top": "#1a2a4a", "bottom": "#0a0a1e"},
    "Sky Gradient": {"mode": "gradient", "top": "#87ceeb", "bottom": "#e0f0ff"},
    "Neutral Gray": {"mode": "solid", "color": "#404040"},
    "Charcoal Gradient": {"mode": "gradient", "top": "#2d2d2d", "bottom": "#1a1a1a"},
    "Engineering Blue": {"mode": "gradient", "top": "#1e3a5f", "bottom": "#0a1628"},
    "Black": {"mode": "solid", "color": "#000000"},
}
