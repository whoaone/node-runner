"""Command Palette: Ctrl+P fuzzy-search over every menu action.

The palette is populated from a flat registry built once after the menu bar
is constructed. Each entry is a tuple of (label_path, keywords, QAction).
Triggering an entry calls ``QAction.trigger()`` so all existing slot wiring
(including command pattern undo/redo, status updates, etc.) flows through
unchanged.
"""

from PySide6.QtCore import Qt, QTimer, QAbstractListModel, QModelIndex, QItemSelectionModel
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLineEdit, QListView, QAbstractItemView,
)


def fuzzy_score(query, candidate):
    """Subsequence fuzzy match score, larger is better.

    Returns None if `query` is not a subsequence (case-insensitive) of
    `candidate`. Otherwise returns a positive score that rewards:
      - matches at word boundaries
      - tightly grouped match positions
      - exact-case matches (small bonus)
    """
    if not query:
        return 0
    q = query.lower()
    c = candidate.lower()
    score = 0
    last_match = -1
    qi = 0
    for i, ch in enumerate(c):
        if qi >= len(q):
            break
        if ch == q[qi]:
            # Word-boundary bonus: first char, or prev char was non-alnum.
            if i == 0 or not c[i - 1].isalnum():
                score += 8
            else:
                score += 2
            # Tightness bonus: closer to last match is better.
            if last_match >= 0:
                gap = i - last_match - 1
                score -= min(gap, 8)
            # Exact-case bonus.
            if candidate[i] == query[qi]:
                score += 1
            last_match = i
            qi += 1
    if qi < len(q):
        return None
    # Shorter candidates win on ties.
    score += max(0, 30 - len(candidate)) // 4
    return score


class _ActionListModel(QAbstractListModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        # List of (label_path, keywords, action) currently shown.
        self._rows = []

    def set_rows(self, rows):
        self.beginResetModel()
        self._rows = rows
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._rows)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or index.row() >= len(self._rows):
            return None
        path, _kw, action = self._rows[index.row()]
        if role == Qt.DisplayRole:
            shortcut = action.shortcut().toString()
            if shortcut:
                return f"{path}     [{shortcut}]"
            return path
        if role == Qt.ToolTipRole:
            return action.toolTip() or path
        if role == Qt.UserRole:
            return action
        return None

    def action_at(self, row):
        if 0 <= row < len(self._rows):
            return self._rows[row][2]
        return None


class CommandPaletteDialog(QDialog):
    """Frameless palette: type to search, Enter to run, Esc to close."""

    MAX_RESULTS = 60

    def __init__(self, registry, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Command Palette")
        self.setWindowFlag(Qt.FramelessWindowHint, True)
        self.setWindowFlag(Qt.Dialog, True)
        self.setModal(True)
        self.resize(600, 360)

        # Center on parent.
        if parent is not None:
            geo = parent.geometry()
            self.move(
                geo.x() + (geo.width() - self.width()) // 2,
                geo.y() + (geo.height() - self.height()) // 3,
            )

        self._registry = [r for r in registry if r[2].isEnabled()]
        self._all_registry = list(registry)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        self.search = QLineEdit(self)
        self.search.setPlaceholderText("Type to search commands...")
        self.search.setStyleSheet(
            "QLineEdit { padding: 8px; font-size: 14px; "
            "border: 1px solid #444; border-radius: 4px; }"
        )
        layout.addWidget(self.search)

        self.list = QListView(self)
        self.list.setUniformItemSizes(True)
        self.list.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.model = _ActionListModel(self)
        self.list.setModel(self.model)
        layout.addWidget(self.list, 1)

        self.search.textChanged.connect(self._refresh)
        self.search.returnPressed.connect(self._run_selected)
        self.list.doubleClicked.connect(lambda _idx: self._run_selected())

        # Initial population (show all enabled actions sorted by path).
        self._refresh("")

    def keyPressEvent(self, ev):
        # Up/Down move the list selection while focus stays in the line edit.
        if ev.key() in (Qt.Key_Down, Qt.Key_Up):
            sm = self.list.selectionModel()
            cur = sm.currentIndex()
            row = cur.row() if cur.isValid() else -1
            if ev.key() == Qt.Key_Down:
                row = min(row + 1, self.model.rowCount() - 1)
            else:
                row = max(row - 1, 0)
            new_idx = self.model.index(max(row, 0), 0)
            sm.setCurrentIndex(new_idx, QItemSelectionModel.SelectCurrent)
            self.list.scrollTo(new_idx)
            ev.accept()
            return
        if ev.key() == Qt.Key_Escape:
            self.reject()
            return
        super().keyPressEvent(ev)

    def _refresh(self, text=""):
        text = (text or "").strip()
        if not text:
            rows = sorted(self._registry, key=lambda r: r[0].lower())[:self.MAX_RESULTS]
        else:
            scored = []
            for path, keywords, action in self._registry:
                hay = f"{path} {keywords}"
                s = fuzzy_score(text, hay)
                if s is not None:
                    scored.append((s, path, keywords, action))
            scored.sort(key=lambda r: (-r[0], r[1].lower()))
            rows = [(p, kw, a) for _s, p, kw, a in scored[:self.MAX_RESULTS]]
        self.model.set_rows(rows)
        # Auto-select the top result.
        if self.model.rowCount() > 0:
            idx = self.model.index(0, 0)
            self.list.selectionModel().setCurrentIndex(
                idx, QItemSelectionModel.SelectCurrent,
            )

    def _run_selected(self):
        idx = self.list.selectionModel().currentIndex()
        if not idx.isValid():
            return
        action = self.model.action_at(idx.row())
        if action is None:
            return
        self.accept()
        # Defer trigger so the dialog finishes closing first; otherwise the
        # triggered action might open another dialog while we are still
        # tearing down.
        QTimer.singleShot(0, action.trigger)


def install_command_palette_shortcut(main_window, registry_provider):
    """Bind Ctrl+P on `main_window` to open the palette.

    `registry_provider` is a no-arg callable that returns the current
    flat action registry. We resolve it at trigger time so newly added
    actions (or enable/disable changes) are reflected.
    """
    shortcut = QShortcut(QKeySequence("Ctrl+P"), main_window)
    shortcut.setContext(Qt.ApplicationShortcut)

    def _open():
        registry = registry_provider() or []
        if not registry:
            return
        dlg = CommandPaletteDialog(registry, parent=main_window)
        dlg.exec()

    shortcut.activated.connect(_open)
    return shortcut
