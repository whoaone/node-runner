"""professional collapsible section widget.

A QWidget with a clickable header (showing an arrow icon and title)
that hides / shows the body widget on click. Used by the v5.2.0
Post-Processing Toolbox to give the user professional collapsible
Results / Deform / Contour sections.

Layout:

    [▼ Title]            <-- clickable header (QToolButton)
        body widget       <-- any QWidget passed via set_body()

When collapsed the body is hidden and the arrow flips to right (▶).
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QToolButton, QFrame, QSizePolicy,
)


class CollapsibleSection(QWidget):
    """A header + body widget that collapses on header click.

    Signals
    -------
    expanded_changed : bool
        Emitted whenever the section expands (True) or collapses
        (False). Useful for telemetry.
    """

    expanded_changed = Signal(bool)

    def __init__(self, title: str, parent=None, expanded: bool = True):
        super().__init__(parent)
        self._title = title
        self._expanded = bool(expanded)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        # Header is a QToolButton with text + arrow icon. We use
        # setArrowType() rather than a custom icon so it renders
        # consistently across platforms.
        self._header = QToolButton(self)
        self._header.setStyleSheet(
            "QToolButton { border: none; text-align: left; "
            "padding: 4px 4px; font-weight: bold; }"
            "QToolButton:hover { background: rgba(255,255,255,0.05); }"
        )
        self._header.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._header.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._header.setCheckable(True)
        self._header.setChecked(self._expanded)
        self._header.setText(title)
        self._update_arrow()
        self._header.toggled.connect(self._on_toggled)
        layout.addWidget(self._header)

        # Body frame holds the content widget set via set_body().
        self._body_frame = QFrame(self)
        self._body_layout = QVBoxLayout(self._body_frame)
        self._body_layout.setContentsMargins(8, 0, 4, 4)
        self._body_layout.setSpacing(4)
        self._body_frame.setVisible(self._expanded)
        layout.addWidget(self._body_frame)

        self._body_widget: QWidget | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_body(self, body: QWidget):
        """Replace (or set) the content widget shown below the header."""
        if self._body_widget is not None:
            self._body_layout.removeWidget(self._body_widget)
            self._body_widget.setParent(None)
        self._body_widget = body
        if body is not None:
            self._body_layout.addWidget(body)

    def body_layout(self) -> QVBoxLayout:
        """Direct access for callers that prefer to add widgets one-by-one."""
        return self._body_layout

    def set_expanded(self, on: bool):
        """Programmatically collapse / expand the section."""
        on = bool(on)
        if on == self._expanded:
            return
        self._header.setChecked(on)   # triggers _on_toggled

    @property
    def is_expanded(self) -> bool:
        return self._expanded

    @property
    def title(self) -> str:
        return self._title

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _on_toggled(self, on: bool):
        self._expanded = bool(on)
        self._body_frame.setVisible(self._expanded)
        self._update_arrow()
        self.expanded_changed.emit(self._expanded)

    def _update_arrow(self):
        if self._expanded:
            self._header.setArrowType(Qt.DownArrow)
        else:
            self._header.setArrowType(Qt.RightArrow)
