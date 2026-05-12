import re

from PySide6 import QtCore
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QLineEdit, QTextEdit, QListWidget,
    QDialogButtonBox, QWidget, QComboBox, QButtonGroup, QToolButton,
    QRadioButton, QStackedWidget, QFrame, QSizePolicy,
    QMenu, QPlainTextEdit, QApplication, QGridLayout,
    QAbstractItemView,
)
from PySide6.QtGui import QAction, QActionGroup, QFont, QShortcut, QKeySequence


# ---------------------------------------------------------------------------
# EntitySelectionBar - Femap-style floating selection window (v3.3.0)
# ---------------------------------------------------------------------------
#
#  Layout (3-column):
#
#  +-----------------------------------------------------------------+
#  | (oAdd)(oRem)(oExcl) +-------------+ [SelAll][SelVis][Reset]     |
#  | ID[__] To[__] By[__]|  bucket     | [More][Previous][Delete]    |
#  | Group [________v]   |  list       | [Pick v][Method v][Grow v]  |
#  +-----------------------------------------------------------------+
#  Count: N / Total
#  -----------------------------------------------------------------
#                                           [Hilite]  [OK] [Cancel]
#
#  v3.3.0: each pick action becomes a *bucket entry* with a sign:
#    add range 1..10 step 2 -> "1, 10, 2"
#    remove single 8        -> "-8"
#    exclude range 5..7     -> "x5, x7, 1"
#  Final selection = walk entries in order, expand each range, apply
#  set arithmetic (+ unions, - and x subtract). Delete removes a row
#  from the bucket (literally undoes that pick action).
#


class _RangeEntry:
    """One bucket row.

    Attributes:
        mode:  '+', '-', or 'x'  (add / remove / exclude)
        start, end, step:        the Femap-style range tuple. For a
                                 single ID, start == end and step == 1.
        label:                   display string with sign prefix.
        ids:                     cached frozenset of integer IDs the
                                 entry expands to.
    """

    __slots__ = ('mode', 'start', 'end', 'step', 'label', 'ids')

    def __init__(self, mode, start, end, step, all_entity_ids=None):
        self.mode = mode
        self.start = int(start)
        self.end = int(end)
        self.step = max(1, int(step))
        prefix = {'+': '', '-': '-', 'x': 'x'}.get(mode, '')
        if self.start == self.end and self.step == 1:
            self.label = f"{prefix}{self.start}"
        else:
            self.label = (f"{prefix}{self.start}, "
                          f"{prefix}{self.end}, {self.step}")
        expanded = set(range(self.start, self.end + 1, self.step))
        if all_entity_ids is not None:
            expanded &= all_entity_ids
        self.ids = frozenset(expanded)

    def expand(self):
        return set(self.ids)


class EntitySelectionBar(QDialog):
    """Femap-style entity selection window - floating, movable, non-modal."""

    request_show_selection = QtCore.Signal(str, list)
    request_picking_mode = QtCore.Signal(str, bool)
    request_advanced_selection = QtCore.Signal(str, float)
    request_selection_mode = QtCore.Signal(str)
    request_previous_selection = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.entity_type = 'Node'
        self.all_entity_ids = set()
        self.selected_ids = set()       # computed from _entries (kept in sync)
        self.single_selection_mode = False
        self.select_by_data = {}
        self.groups = {}
        self._current_method = "ID"
        # v3.3.0: each entry is a _RangeEntry with sign + (start,end,step).
        self._entries = []  # type: list[_RangeEntry]

        self.setWindowTitle("Entity Selection")
        self.setWindowFlags(
            QtCore.Qt.Dialog
            | QtCore.Qt.WindowTitleHint
            | QtCore.Qt.CustomizeWindowHint
        )
        self.setModal(False)

        self.setStyleSheet("""
            QLabel { font-size: 11px; }
            QLineEdit { padding: 1px 4px; min-height: 18px; font-size: 11px; }
            QPushButton, QToolButton { font-size: 11px; padding: 2px 6px;
                                       min-height: 20px; }
            QComboBox { font-size: 11px; padding: 1px 4px; min-height: 18px; }
            QRadioButton { font-size: 11px; spacing: 3px; }
            QListWidget { font-family: Consolas, monospace; font-size: 12px;
                          border: 1px solid #aaaaaa; border-radius: 2px;
                          padding: 2px; }
        """)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(5, 5, 5, 5)
        outer.setSpacing(2)

        # Title label
        self._title_label = QLabel("Select Node")
        self._title_label.setStyleSheet("font-weight: bold; font-size: 12px;")
        outer.addWidget(self._title_label)

        body = QHBoxLayout()
        body.setSpacing(4)

        # ----- LEFT COLUMN -----
        left = QVBoxLayout()
        left.setSpacing(2)

        # Action mode radios
        radio_row = QHBoxLayout()
        radio_row.setSpacing(6)
        self._action_group = QButtonGroup(self)
        self._add_radio = QRadioButton("Add")
        self._add_radio.setChecked(True)
        self._remove_radio = QRadioButton("Remove")
        self._exclude_radio = QRadioButton("Exclude")
        self._action_group.addButton(self._add_radio)
        self._action_group.addButton(self._remove_radio)
        self._action_group.addButton(self._exclude_radio)
        for r in (self._add_radio, self._remove_radio, self._exclude_radio):
            radio_row.addWidget(r)
        radio_row.addStretch()
        left.addLayout(radio_row)

        # Method-dependent input (QStackedWidget)
        self._method_stack = QStackedWidget()

        # Page 0: ID / to / by
        id_page = QWidget()
        id_lay = QHBoxLayout(id_page)
        id_lay.setContentsMargins(0, 0, 0, 0)
        id_lay.setSpacing(3)
        id_lay.addWidget(QLabel("ID"))
        self._id_input = QLineEdit()
        self._id_input.setPlaceholderText("ID")
        self._id_input.setMinimumWidth(40)
        id_lay.addWidget(self._id_input, 2)
        id_lay.addWidget(QLabel("To:"))
        self._to_input = QLineEdit()
        self._to_input.setPlaceholderText("to")
        self._to_input.setMinimumWidth(40)
        id_lay.addWidget(self._to_input, 2)
        id_lay.addWidget(QLabel("By"))
        self._by_input = QLineEdit()
        self._by_input.setPlaceholderText("by")
        self._by_input.setMinimumWidth(30)
        id_lay.addWidget(self._by_input, 1)
        self._method_stack.addWidget(id_page)

        # Page 1: value combo (Property / Material / Type / Quality)
        value_page = QWidget()
        val_lay = QHBoxLayout(value_page)
        val_lay.setContentsMargins(0, 0, 0, 0)
        val_lay.addWidget(QLabel("Value:"))
        self._value_combo = QComboBox()
        val_lay.addWidget(self._value_combo, 1)
        self._method_stack.addWidget(value_page)

        # Page 2: stub
        stub_page = QWidget()
        stub_lay = QHBoxLayout(stub_page)
        stub_lay.setContentsMargins(0, 0, 0, 0)
        stub_lbl = QLabel("(Not yet implemented)")
        stub_lbl.setStyleSheet("color: #888; font-style: italic;")
        stub_lay.addWidget(stub_lbl)
        self._method_stack.addWidget(stub_page)

        left.addWidget(self._method_stack)

        # Group combo
        grp_row = QHBoxLayout()
        grp_row.setSpacing(4)
        self._group_label = QLabel("Group")
        grp_row.addWidget(self._group_label)
        self._group_combo = QComboBox()
        self._group_combo.addItem("")
        self._group_combo.setMinimumWidth(80)
        grp_row.addWidget(self._group_combo, 1)
        left.addLayout(grp_row)

        left.addStretch()

        body.addLayout(left, 0)

        # ----- CENTER: bucket list -----
        self._list_widget = QListWidget()
        self._list_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._list_widget.setMinimumWidth(160)
        self._list_widget.setMinimumHeight(110)
        body.addWidget(self._list_widget, 1)

        # ----- RIGHT COLUMN: 3x3 button grid (v3.3.0) -----
        # Row 0: Select All  | Select Visible | Reset
        # Row 1: More        | Previous       | Delete
        # Row 2: Pick v      | Method v       | Grow v
        grid = QGridLayout()
        grid.setSpacing(3)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(2, 1)

        # Row 0
        self._all_btn = QPushButton("Select All")
        self._all_btn.setToolTip("Select all entities (one entry).")
        grid.addWidget(self._all_btn, 0, 0)

        self._visible_btn = QPushButton("Select Visible")
        self._visible_btn.setToolTip(
            "Select every entity currently visible in the viewport "
            "(respects tree on/off state).")
        grid.addWidget(self._visible_btn, 0, 1)

        self._reset_btn = QPushButton("Reset")
        self._reset_btn.setToolTip("Clear entire bucket.")
        grid.addWidget(self._reset_btn, 0, 2)

        # Row 1
        self._more_btn = QPushButton("More")
        self._more_btn.setToolTip(
            "Apply the current ID/To/By input as a new bucket row.")
        grid.addWidget(self._more_btn, 1, 0)

        self._prev_btn = QPushButton("Previous")
        self._prev_btn.setToolTip("Restore previous selection.")
        grid.addWidget(self._prev_btn, 1, 1)

        self._delete_btn = QPushButton("Delete")
        self._delete_btn.setToolTip(
            "Remove the highlighted bucket row(s) - undoes those picks.")
        grid.addWidget(self._delete_btn, 1, 2)

        # Row 2
        self._pick_btn = QToolButton()
        self._pick_btn.setText("Pick ▾")
        self._pick_btn.setToolTip(
            "Pick entities from viewport; dropdown for pick mode "
            "(Box/Circle/Polygon).")
        self._pick_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._pick_btn.setPopupMode(QToolButton.MenuButtonPopup)
        self._pick_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._pick_menu = QMenu(self)
        self._build_pick_menu()
        self._pick_btn.setMenu(self._pick_menu)
        grid.addWidget(self._pick_btn, 2, 0)

        self._method_btn = QToolButton()
        self._method_btn.setText("Method ▾")
        self._method_btn.setToolTip("Choose selection method.")
        self._method_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._method_btn.setPopupMode(QToolButton.InstantPopup)
        self._method_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._method_menu = QMenu(self)
        self._method_btn.setMenu(self._method_menu)
        grid.addWidget(self._method_btn, 2, 1)

        self._grow_menu_btn = QToolButton()
        self._grow_menu_btn.setText("Grow ▾")
        self._grow_menu_btn.setToolTip(
            "Adjacent / Connected / Grow / Shrink ; angle threshold.")
        self._grow_menu_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._grow_menu_btn.setPopupMode(QToolButton.InstantPopup)
        self._grow_menu_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._grow_menu = QMenu(self)
        self._grow_menu.addAction("Adjacent", self._request_select_adjacent)
        self._grow_menu.addAction("Connected", self._request_select_connected)
        self._grow_menu.addAction("Grow by 1 layer", self._request_grow)
        self._grow_menu.addAction("Shrink by 1 layer", self._request_shrink)
        self._grow_menu.addSeparator()
        self._grow_menu.addAction("Set angle threshold...",
                                  self._prompt_angle_threshold)
        self._grow_menu_btn.setMenu(self._grow_menu)
        grid.addWidget(self._grow_menu_btn, 2, 2)

        body.addLayout(grid, 0)

        outer.addLayout(body)

        # Count label
        self._count_label = QLabel("Count: 0 / 0")
        self._count_label.setStyleSheet("font-weight: bold; font-size: 11px;")
        self._count_label.setToolTip("Final selected / Total entities")
        outer.addWidget(self._count_label)

        # Hidden adv frame (kept for back-compat with handlers, never shown)
        self._adv_frame = QFrame()
        adv_lay = QHBoxLayout(self._adv_frame)
        adv_lay.setContentsMargins(0, 0, 0, 0)
        adv_lay.setSpacing(4)
        adv_lay.addWidget(QLabel("Angle:"))
        self._angle_input = QLineEdit("30")
        self._angle_input.setMaximumWidth(35)
        adv_lay.addWidget(self._angle_input)
        adv_lay.addWidget(QLabel("°"))
        self._adj_btn = QPushButton("Adj")
        self._conn_btn = QPushButton("Conn")
        self._grow_btn = QPushButton("Grow")
        self._shrink_btn = QPushButton("Shrink")
        for btn in (self._adj_btn, self._conn_btn, self._grow_btn, self._shrink_btn):
            btn.setMaximumWidth(52)
            adv_lay.addWidget(btn)
        adv_lay.addStretch()
        self._adv_frame.hide()
        self._adv_frame.setVisible(False)

        # Separator + bottom row: Hilite + OK / Cancel (v3.3.0)
        _sep = QFrame()
        _sep.setFrameShape(QFrame.HLine)
        _sep.setFrameShadow(QFrame.Sunken)
        outer.addWidget(_sep)

        _bottom = QHBoxLayout()
        _bottom.setSpacing(3)
        _bottom.addStretch(1)

        self._highlight_btn = QToolButton()
        self._highlight_btn.setText("Hilite")
        self._highlight_btn.setToolTip("Toggle highlighting in viewport.")
        self._highlight_btn.setCheckable(True)
        self._highlight_btn.setChecked(True)
        self._highlight_btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
        self._highlight_btn.setMinimumWidth(70)
        self._highlight_btn.setStyleSheet(
            "QToolButton:checked { font-weight: bold; background: #f9e2af; "
            "color: #1e1e2e; border: 1px solid #cba651; border-radius: 3px; }")
        _bottom.addWidget(self._highlight_btn)

        self._ok_btn = QPushButton("OK")
        self._ok_btn.setStyleSheet("font-weight: bold;")
        self._ok_btn.setMinimumWidth(80)
        _bottom.addWidget(self._ok_btn)

        self._cancel_btn = QPushButton("Cancel")
        self._cancel_btn.setMinimumWidth(80)
        _bottom.addWidget(self._cancel_btn)

        outer.addLayout(_bottom)

        # === Connect internal signals ===
        self._pick_btn.clicked.connect(self._activate_picking)
        self._more_btn.clicked.connect(self._on_more_clicked)
        self._prev_btn.clicked.connect(self._request_previous)
        self._all_btn.clicked.connect(self._select_all)
        self._visible_btn.clicked.connect(self._select_visible)
        self._reset_btn.clicked.connect(self._clear_all)
        self._delete_btn.clicked.connect(self._delete_selected_entries)
        self._highlight_btn.toggled.connect(self._on_highlight_toggled)

        # ID field Enter -> More
        self._id_input.returnPressed.connect(self._on_more_clicked)
        self._to_input.returnPressed.connect(self._on_more_clicked)
        self._by_input.returnPressed.connect(self._on_more_clicked)

        # Advanced selection (legacy hidden buttons)
        self._adj_btn.clicked.connect(self._request_select_adjacent)
        self._conn_btn.clicked.connect(self._request_select_connected)
        self._grow_btn.clicked.connect(self._request_grow)
        self._shrink_btn.clicked.connect(self._request_shrink)

        # OK / Cancel
        self._ok_btn.clicked.connect(self._on_ok)
        self._cancel_btn.clicked.connect(self._on_cancel)

        # v3.4.0 (item 1): Ctrl+V anywhere in the dialog pastes clipboard
        # IDs into the bucket via the same _paste_selection path that
        # the Pick > Paste menu uses. Works regardless of which child
        # widget has focus.
        self._paste_shortcut = QShortcut(
            QKeySequence.Paste, self, activated=self._paste_selection)
        self._paste_shortcut.setContext(QtCore.Qt.WidgetWithChildrenShortcut)

    # ------------------------------------------------------------------
    # Pick menu builder
    # ------------------------------------------------------------------

    def _build_pick_menu(self):
        """Build the Pick button's dropdown menu."""
        menu = self._pick_menu
        menu.clear()

        # Picking modes (checkable, mutually exclusive)
        self._pick_mode_group = QActionGroup(self)
        self._pick_mode_group.setExclusive(True)

        box_act = menu.addAction("Box")
        box_act.setCheckable(True)
        box_act.setChecked(True)
        self._pick_mode_group.addAction(box_act)

        circle_act = menu.addAction("Circle")
        circle_act.setCheckable(True)
        self._pick_mode_group.addAction(circle_act)

        poly_act = menu.addAction("Polygon")
        poly_act.setCheckable(True)
        self._pick_mode_group.addAction(poly_act)

        menu.addSeparator()

        # Advanced picking (stubs - disabled)
        for label in ("Coordinate...", "Around Point...",
                      "Around Vector...", "Around Plane..."):
            act = menu.addAction(label)
            act.setEnabled(False)

        menu.addSeparator()

        # Selection by visual (stubs - disabled)
        for label in ("By Color...", "By Faces..."):
            act = menu.addAction(label)
            act.setEnabled(False)

        menu.addSeparator()

        # Clipboard operations
        copy_act = menu.addAction("Copy")
        copy_list_act = menu.addAction("Copy as List")
        paste_act = menu.addAction("Paste")

        # Connect picking mode actions
        box_act.triggered.connect(lambda: self.request_selection_mode.emit('box'))
        circle_act.triggered.connect(lambda: self.request_selection_mode.emit('circle'))
        poly_act.triggered.connect(lambda: self.request_selection_mode.emit('polygon'))

        # Connect clipboard actions
        copy_act.triggered.connect(self._copy_selection_ranges)
        copy_list_act.triggered.connect(self._copy_selection_list)
        paste_act.triggered.connect(self._paste_selection)

        # Store references
        self._pick_box_action = box_act
        self._pick_circle_action = circle_act
        self._pick_polygon_action = poly_act

    # ------------------------------------------------------------------
    # Method menu builder
    # ------------------------------------------------------------------

    def _build_method_menu(self):
        """Build the Method button's dropdown menu from select_by_data."""
        menu = self._method_menu
        menu.clear()

        self._method_action_group = QActionGroup(self)
        self._method_action_group.setExclusive(True)

        id_act = menu.addAction("ID")
        id_act.setCheckable(True)
        id_act.setChecked(True)
        self._method_action_group.addAction(id_act)
        id_act.triggered.connect(lambda: self._on_method_changed("ID"))

        if self.select_by_data:
            menu.addSeparator()
            for key in self.select_by_data:
                act = menu.addAction(key)
                act.setCheckable(True)
                self._method_action_group.addAction(act)
                act.triggered.connect(
                    lambda checked, k=key: self._on_method_changed(k))

        menu.addSeparator()

        # Stubs (disabled)
        for label in ("Shape...", "Color...", "Using Node..."):
            act = menu.addAction(label)
            act.setEnabled(False)

    # ------------------------------------------------------------------
    # configure - reset for a new selection workflow
    # ------------------------------------------------------------------

    def configure(self, entity_type, all_entity_ids, select_by_data=None,
                  groups=None, single_selection_mode=False):
        """Reset the dialog for a new entity selection workflow."""
        self.entity_type = entity_type
        self.all_entity_ids = set(all_entity_ids)
        self.selected_ids = set()
        self._entries = []
        self.single_selection_mode = single_selection_mode
        self.select_by_data = select_by_data or {}
        self.groups = groups or {}
        self._current_method = "ID"

        self._title_label.setText(f"Select {entity_type}")
        self.setWindowTitle(f"{entity_type} Selection")

        # Reset action mode to Add
        self._add_radio.setChecked(True)

        # Rebuild Method menu
        self._build_method_menu()
        self._method_btn.setText("Method ▾")
        self._method_stack.setCurrentIndex(0)

        # Rebuild Group combo - always visible, shows "GID: Name"
        self._group_combo.blockSignals(True)
        self._group_combo.clear()
        self._group_combo.addItem("")
        for name in sorted(self.groups.keys(),
                           key=lambda n: self.groups[n].get("gid", 0)):
            gid = self.groups[name].get("gid", 0)
            self._group_combo.addItem(f"{gid}: {name}", name)
        self._group_combo.setCurrentIndex(0)
        self._group_combo.blockSignals(False)

        # Adv frame stays hidden; Grow submenu is the path.
        self._adv_frame.setVisible(False)

        self._clear_inputs()
        self._list_widget.clear()
        self._highlight_btn.setChecked(True)
        self._update_count()

        # Position near the bottom-right of the parent window
        self.adjustSize()
        if self.parent():
            pg = self.parent().frameGeometry()
            self.move(
                pg.x() + pg.width() - self.sizeHint().width() - 30,
                pg.y() + pg.height() - self.sizeHint().height() - 80,
            )

    # ------------------------------------------------------------------
    # Action mode
    # ------------------------------------------------------------------

    def get_action_mode(self):
        """Return current action mode: 'add', 'remove', or 'exclude'."""
        if self._remove_radio.isChecked():
            return 'remove'
        if self._exclude_radio.isChecked():
            return 'exclude'
        return 'add'

    def _mode_char(self):
        """Action mode as the single-char prefix used by entries."""
        m = self.get_action_mode()
        return {'add': '+', 'remove': '-', 'exclude': 'x'}[m]

    # ------------------------------------------------------------------
    # Method switching
    # ------------------------------------------------------------------

    def _on_method_changed(self, text):
        self._current_method = text
        if text == "ID":
            self._method_stack.setCurrentIndex(0)
        elif text in self.select_by_data:
            self._method_stack.setCurrentIndex(1)
            self._update_value_combo(text)
        else:
            self._method_stack.setCurrentIndex(2)
        self._method_btn.setText("Method ▾")

    def _update_value_combo(self, method_text):
        """Populate the value combo from select_by_data."""
        self._value_combo.clear()
        data = self.select_by_data.get(method_text)
        if data == '__deferred_quality__':
            parent = self.parent()
            resolver = getattr(parent, '_compute_quality_select_data', None)
            if callable(resolver):
                try:
                    data = resolver()
                    self.select_by_data[method_text] = data
                except Exception:
                    data = {}
            else:
                data = {}
        if isinstance(data, dict):
            for key in sorted(data.keys(),
                              key=lambda k: (isinstance(k, str), k)):
                count = len(data[key])
                self._value_combo.addItem(f"{key} ({count} elems)", key)

    # ------------------------------------------------------------------
    # Bucket entry helpers
    # ------------------------------------------------------------------

    def _append_entry(self, mode_char, start, end, step):
        """Append a new _RangeEntry to the bucket and refresh the UI.

        ``mode_char`` is '+', '-', or 'x'. The entry's expanded ID set
        is intersected against ``self.all_entity_ids`` so dangling IDs
        never leak into the final selection.
        """
        entry = _RangeEntry(mode_char, start, end, step,
                            all_entity_ids=self.all_entity_ids)
        if not entry.ids:
            return  # nothing to add
        self._entries.append(entry)
        self._refresh_list()
        self._show_selection()

    def _append_id_set(self, mode_char, ids):
        """Decompose an arbitrary id set into contiguous ranges and
        append one entry per range. Used by Pick / Select All / Select
        Visible / Paste."""
        if not ids:
            return
        valid = sorted(int(i) for i in self.all_entity_ids.intersection(ids))
        if not valid:
            return
        if self.single_selection_mode and mode_char == '+':
            # Keep only the first id, replace bucket.
            self._entries.clear()
            valid = [valid[0]]
        for start, end, step in self.compress_ids_to_range_tuples(valid):
            entry = _RangeEntry(mode_char, start, end, step,
                                all_entity_ids=self.all_entity_ids)
            if entry.ids:
                self._entries.append(entry)
        self._refresh_list()
        self._show_selection()
        if self.single_selection_mode and self.selected_ids:
            self.accepted.emit()

    # ------------------------------------------------------------------
    # Input collection & application
    # ------------------------------------------------------------------

    def _collect_current_input_range(self):
        """Read current Method + inputs and return either a
        ``(start, end, step)`` tuple OR a set of IDs (for value-based
        methods / groups). Returns None on empty / invalid input.
        """
        # Group dropdown wins if active
        group_name = self._group_combo.currentData()
        if group_name and group_name in self.groups:
            group = self.groups[group_name]
            if self.entity_type == 'Node':
                return set(group.get("nodes", []))
            elif self.entity_type == 'Element':
                return set(group.get("elements", []))
            return None

        method = self._current_method
        if method != "ID":
            value = self._value_combo.currentData()
            if method in self.select_by_data and value in self.select_by_data[method]:
                return set(self.select_by_data[method][value])
            return None

        id_text = self._id_input.text().strip()
        to_text = self._to_input.text().strip()
        by_text = self._by_input.text().strip()
        if not id_text:
            return None
        try:
            id_val = int(id_text)
        except ValueError:
            return None
        if to_text:
            try:
                to_val = int(to_text)
            except ValueError:
                return None
            by_val = 1
            if by_text:
                try:
                    by_val = max(1, int(by_text))
                except ValueError:
                    by_val = 1
            start, end = min(id_val, to_val), max(id_val, to_val)
            return (start, end, by_val)
        return (id_val, id_val, 1)

    def _on_more_clicked(self):
        """Apply current input as a new bucket row (with the mode sign)."""
        result = self._collect_current_input_range()
        if result is None:
            return
        mode = self._mode_char()
        if isinstance(result, tuple):
            start, end, step = result
            self._append_entry(mode, start, end, step)
        else:
            self._append_id_set(mode, result)
        self._clear_inputs()

    def _clear_inputs(self):
        self._id_input.clear()
        self._to_input.clear()
        self._by_input.clear()
        self._group_combo.setCurrentIndex(0)

    # ------------------------------------------------------------------
    # External selection-manipulation hooks (interactor / parent)
    # ------------------------------------------------------------------

    def add_selection(self, ids_to_add):
        """Bulk-add IDs as one or more add-entries (signed '+')."""
        self._append_id_set('+', set(ids_to_add))

    def remove_selection(self, ids_to_remove):
        """Bulk-remove IDs as one or more remove-entries (signed '-')."""
        self._append_id_set('-', set(ids_to_remove))

    def replace_selection(self, new_ids):
        """Clear bucket then add given IDs as one or more add-entries."""
        self._entries.clear()
        self._append_id_set('+', set(new_ids))

    def toggle_selection(self, id_to_toggle):
        """Pick-mode hook for a single-click. Mode comes from the radio:
        Add appends a '+', Remove appends a '-', Exclude appends 'x'.

        Toggle semantics: if Add mode and the id is already in the
        selection, the click is treated as a Remove (so users can pick
        in Add mode and click again to undo). For Remove/Exclude, the
        click always appends an entry of that mode.
        """
        if self.single_selection_mode:
            self._entries.clear()
            if id_to_toggle in self.all_entity_ids:
                self._append_entry('+', id_to_toggle, id_to_toggle, 1)
            return

        mode = self._mode_char()
        if mode == '+':
            if id_to_toggle in self.selected_ids:
                # Convenience toggle: clicking a selected id in Add
                # mode removes it.
                mode = '-'
        if id_to_toggle in self.all_entity_ids:
            self._append_entry(mode, id_to_toggle, id_to_toggle, 1)

    def get_selected_ids(self):
        """Compute final selection by walking the bucket in order.

        + entries union their IDs into the result; - and x entries
        subtract. Order matters: an 'add 1..10' then 'remove 8' yields
        {1..7, 9, 10}. The same two entries in opposite order would
        yield {1..10} (the remove subtracts from an empty result, then
        the add unions in everything).
        """
        result = set()
        for entry in self._entries:
            ids = entry.ids
            if entry.mode == '+':
                result |= ids
            else:
                result -= ids
        return sorted(result)

    # ------------------------------------------------------------------
    # Utility actions
    # ------------------------------------------------------------------

    def _select_all(self):
        """Add an add-entry containing every entity (compressed range)."""
        self._append_id_set('+', set(self.all_entity_ids))

    def _select_visible(self):
        """v3.3.0: add an entry containing every currently-visible entity.

        'Visible' = the parent MainWindow's tree-checked filter. We
        request the visible-id set from the parent via the
        ``visible_entity_ids(entity_type)`` callback. Falls back to
        select-all if the parent doesn't expose that helper.
        """
        parent = self.parent()
        resolver = getattr(parent, 'visible_entity_ids', None)
        if callable(resolver):
            try:
                vis = resolver(self.entity_type)
            except Exception:
                vis = None
        else:
            vis = None
        if not vis:
            # Conservative fallback: every entity is "visible".
            vis = set(self.all_entity_ids)
        self._append_id_set(self._mode_char(), set(vis))

    def _clear_all(self):
        """Empty the bucket (Reset)."""
        self._entries.clear()
        self._refresh_list()
        self.request_show_selection.emit(self.entity_type, [])

    def _delete_selected_entries(self):
        """Remove highlighted bucket rows from _entries."""
        rows = sorted(
            {idx.row() for idx in self._list_widget.selectedIndexes()},
            reverse=True)
        for row in rows:
            if 0 <= row < len(self._entries):
                del self._entries[row]
        self._refresh_list()
        self._show_selection()

    def _show_selection(self):
        """Emit highlight signal only if the highlight toggle is ON."""
        if self._highlight_btn.isChecked():
            self.request_show_selection.emit(
                self.entity_type, self.get_selected_ids())

    def _on_highlight_toggled(self, checked):
        if checked:
            ids = self.get_selected_ids()
            self.request_show_selection.emit(self.entity_type, ids)
            if not ids:
                self._count_label.setText(
                    f"Count: 0 / {len(self.all_entity_ids):,}  "
                    f"(nothing selected to highlight yet)")
        else:
            self.request_show_selection.emit(self.entity_type, [])
            self._count_label.setText(
                f"Count: {len(self.get_selected_ids()):,} / "
                f"{len(self.all_entity_ids):,}  (highlighting OFF)")

    # ------------------------------------------------------------------
    # Picking / Previous
    # ------------------------------------------------------------------

    def _activate_picking(self):
        self.request_picking_mode.emit(self.entity_type, True)

    def _request_previous(self):
        self.request_previous_selection.emit(self.entity_type)

    # ------------------------------------------------------------------
    # Advanced selection (Elements only)
    # ------------------------------------------------------------------

    def _prompt_angle_threshold(self):
        from PySide6.QtWidgets import QInputDialog
        try:
            current = float(self._angle_input.text())
        except (TypeError, ValueError):
            current = 30.0
        value, ok = QInputDialog.getDouble(
            self, "Set angle threshold",
            "Dihedral angle (degrees) for Adjacent / Grow:",
            current, 0.0, 180.0, 1)
        if ok:
            self._angle_input.setText(f"{value:g}")

    def _request_select_adjacent(self):
        try:
            angle = float(self._angle_input.text())
        except (ValueError, AttributeError):
            angle = 30.0
        self.request_advanced_selection.emit('adjacent', angle)

    def _request_select_connected(self):
        self.request_advanced_selection.emit('connected', 180.0)

    def _request_grow(self):
        try:
            angle = float(self._angle_input.text())
        except (ValueError, AttributeError):
            angle = 30.0
        self.request_advanced_selection.emit('grow', angle)

    def _request_shrink(self):
        self.request_advanced_selection.emit('shrink', 0.0)

    # ------------------------------------------------------------------
    # Range compression
    # ------------------------------------------------------------------

    @staticmethod
    def compress_ids_to_range_tuples(sorted_ids):
        """Compress a sorted iterable of ints into list of (start, end, step)
        tuples. Single IDs are returned as (id, id, 1)."""
        ids = list(sorted_ids)
        if not ids:
            return []
        out = []
        i = 0
        n = len(ids)
        while i < n:
            start = ids[i]
            if i + 1 >= n:
                out.append((start, start, 1))
                i += 1
                continue
            step = ids[i + 1] - start
            if step <= 0:
                out.append((start, start, 1))
                i += 1
                continue
            j = i + 1
            while j < n and ids[j] - ids[j - 1] == step:
                j += 1
            end = ids[j - 1]
            out.append((start, end, step))
            i = j
        return out

    @staticmethod
    def compress_ids_to_ranges(sorted_ids):
        """Back-compat: string form for clipboard / legacy callers.

        Returns a list of strings:
          "5"       - single ID
          "1,5,1"   - consecutive IDs 1..5 step 1
          "1,9,2"   - stepped IDs 1, 3, 5, 7, 9
        """
        parts = []
        for s, e, k in EntitySelectionBar.compress_ids_to_range_tuples(
                sorted_ids):
            if s == e:
                parts.append(str(s))
            else:
                parts.append(f"{s},{e},{k}")
        return parts

    # ------------------------------------------------------------------
    # Clipboard operations
    # ------------------------------------------------------------------

    def _copy_selection_ranges(self):
        sorted_ids = self.get_selected_ids()
        ranges = self.compress_ids_to_ranges(sorted_ids)
        QApplication.clipboard().setText('\n'.join(ranges))

    def _copy_selection_list(self):
        sorted_ids = self.get_selected_ids()
        QApplication.clipboard().setText(
            ','.join(str(i) for i in sorted_ids))

    def _paste_selection(self):
        text = QApplication.clipboard().text()
        if not text:
            return
        ids = set()
        for token in re.split(r'[\n;]+', text.strip()):
            token = token.strip()
            if not token:
                continue
            parts = [p.strip() for p in token.split(',')]
            if len(parts) == 3:
                try:
                    start, end, step = (int(parts[0]), int(parts[1]),
                                        int(parts[2]))
                    if step > 0:
                        ids.update(range(start, end + 1, step))
                        continue
                except ValueError:
                    pass
            for p in parts:
                try:
                    ids.add(int(p.strip()))
                except ValueError:
                    pass
        if ids:
            self._append_id_set(self._mode_char(), ids)

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------

    def _refresh_list(self):
        """Rebuild the list widget from entries and recompute selected_ids."""
        self._list_widget.setUpdatesEnabled(False)
        try:
            self._list_widget.clear()
            for entry in self._entries:
                self._list_widget.addItem(entry.label)
        finally:
            self._list_widget.setUpdatesEnabled(True)
        self.selected_ids = set(self.get_selected_ids())
        self._update_count()

    def _update_count(self):
        n_sel = len(self.selected_ids)
        n_total = len(self.all_entity_ids)
        self._count_label.setText(f"Count: {n_sel} / {n_total}")

    def _on_ok(self):
        self.accept()

    def _on_cancel(self):
        self.reject()


# ---------------------------------------------------------------------------
# EntitySelectionDialog - kept for 3 modal Geometry Point usages
# ---------------------------------------------------------------------------

class EntitySelectionDialog(QDialog):
    """Simple modal selection dialog for Geometry Points (exec-based)."""

    request_show_selection = QtCore.Signal(str, list)
    request_picking_mode = QtCore.Signal(str, bool)
    request_advanced_selection = QtCore.Signal(str, float)
    request_selection_mode = QtCore.Signal(str)
    request_previous_selection = QtCore.Signal(str)

    def __init__(self, entity_type, all_entity_ids, parent=None,
                 single_selection_mode=False, select_by_data=None):
        super().__init__(parent)
        self.entity_type = entity_type
        self.setWindowTitle(f"{self.entity_type} Selection")

        self.all_entity_ids, self.selected_ids = set(all_entity_ids), set()
        self.single_selection_mode = single_selection_mode
        self.select_by_data = select_by_data or {}

        self.setMaximumWidth(380)

        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(6)

        title = QLabel(f"Select {entity_type}(s)")
        title.setStyleSheet("font-weight: bold; font-size: 12px;")
        main_layout.addWidget(title)

        self.count_label = QLabel("0 of 0 selected")
        self.count_label.setStyleSheet("font-size: 11px;")
        main_layout.addWidget(self.count_label)

        self.list_widget = QListWidget()
        if self.single_selection_mode:
            self.list_widget.setSelectionMode(QListWidget.SingleSelection)
        self.list_widget.setMaximumHeight(80)
        main_layout.addWidget(self.list_widget)

        action_layout = QHBoxLayout()
        action_layout.setSpacing(2)
        select_all_button = QToolButton(); select_all_button.setText("All")
        invert_button = QToolButton(); invert_button.setText("Invert")
        remove_button = QToolButton(); remove_button.setText("Rm")
        clear_button = QToolButton(); clear_button.setText("Clear")
        for btn in (select_all_button, invert_button, remove_button, clear_button):
            btn.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            action_layout.addWidget(btn)
        action_layout.addStretch()
        main_layout.addLayout(action_layout)

        id_row = QHBoxLayout()
        id_row.addWidget(QLabel("ID:"))
        self.single_id_input = QLineEdit()
        self.single_id_input.setMaximumWidth(80)
        add_single_button = QPushButton("Add")
        add_single_button.setMaximumWidth(50)
        id_row.addWidget(self.single_id_input)
        id_row.addWidget(add_single_button)
        id_row.addWidget(QLabel("  Range:"))
        self.start_id_input = QLineEdit()
        self.start_id_input.setPlaceholderText("Start")
        self.start_id_input.setMaximumWidth(60)
        self.end_id_input = QLineEdit()
        self.end_id_input.setPlaceholderText("End")
        self.end_id_input.setMaximumWidth(60)
        add_range_button = QPushButton("Add Range")
        add_range_button.setMaximumWidth(70)
        id_row.addWidget(self.start_id_input)
        id_row.addWidget(QLabel("-"))
        id_row.addWidget(self.end_id_input)
        id_row.addWidget(add_range_button)
        id_row.addStretch()
        main_layout.addLayout(id_row)

        paste_row = QHBoxLayout()
        self.paste_edit = QTextEdit()
        self.paste_edit.setPlaceholderText("Paste IDs (comma/space/newline)")
        self.paste_edit.setMaximumHeight(40)
        add_list_button = QPushButton("Add List")
        add_list_button.setMaximumWidth(60)
        paste_row.addWidget(self.paste_edit)
        paste_row.addWidget(add_list_button)
        main_layout.addLayout(paste_row)

        if self.single_selection_mode:
            select_all_button.setEnabled(False)
            invert_button.setEnabled(False)
            self.start_id_input.setEnabled(False)
            self.end_id_input.setEnabled(False)
            add_range_button.setEnabled(False)
            self.paste_edit.setEnabled(False)
            add_list_button.setEnabled(False)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        main_layout.addWidget(button_box)

        select_all_button.clicked.connect(self._select_all)
        invert_button.clicked.connect(self._invert_selection)
        remove_button.clicked.connect(self._remove_selected)
        clear_button.clicked.connect(self._clear_all)
        add_single_button.clicked.connect(self._add_single_id)
        add_range_button.clicked.connect(self._add_range_ids)
        add_list_button.clicked.connect(self._add_list_ids)
        self.single_id_input.returnPressed.connect(self._add_single_id)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        self._update_count_label()
        self.adjustSize()

    def add_selection(self, ids_to_add):
        if self.single_selection_mode:
            self.selected_ids.clear()
            if ids_to_add:
                ids_to_add = {list(ids_to_add)[0]}
        valid_ids = self.all_entity_ids.intersection(set(ids_to_add))
        self.selected_ids.update(valid_ids)
        self._refresh_list_widget()
        if self.single_selection_mode and self.selected_ids:
            self.accept()

    def remove_selection(self, ids_to_remove):
        self.selected_ids -= set(ids_to_remove)
        self._refresh_list_widget()

    def replace_selection(self, new_ids):
        self.selected_ids = self.all_entity_ids.intersection(set(new_ids))
        self._refresh_list_widget()

    def toggle_selection(self, id_to_toggle):
        if self.single_selection_mode:
            self.selected_ids.clear()
            if id_to_toggle in self.all_entity_ids:
                self.selected_ids.add(id_to_toggle)
            self._refresh_list_widget()
            if self.selected_ids:
                self.accept()
            return
        if id_to_toggle in self.selected_ids:
            self.selected_ids.remove(id_to_toggle)
        elif id_to_toggle in self.all_entity_ids:
            self.selected_ids.add(id_to_toggle)
        self._refresh_list_widget()

    def get_selected_ids(self):
        return sorted(list(self.selected_ids))

    def _select_all(self):
        self.selected_ids = self.all_entity_ids.copy()
        self._refresh_list_widget()

    def _invert_selection(self):
        self.selected_ids = self.all_entity_ids - self.selected_ids
        self._refresh_list_widget()

    def _refresh_list_widget(self):
        self.list_widget.clear()
        self.list_widget.addItems([str(i) for i in sorted(self.selected_ids)])
        self._update_count_label()

    def _update_count_label(self):
        self.count_label.setText(
            f"{len(self.selected_ids)} of {len(self.all_entity_ids)} selected")

    def _add_single_id(self):
        try:
            self.add_selection({int(self.single_id_input.text())})
            self.single_id_input.clear()
        except ValueError:
            pass

    def _add_range_ids(self):
        try:
            start_id = int(self.start_id_input.text())
            end_id = int(self.end_id_input.text())
            self.add_selection(set(range(min(start_id, end_id), max(start_id, end_id) + 1)))
            self.start_id_input.clear()
            self.end_id_input.clear()
        except ValueError:
            pass

    def _add_list_ids(self):
        try:
            self.add_selection(
                {int(i) for i in re.findall(r'\d+', self.paste_edit.toPlainText())})
            self.paste_edit.clear()
        except ValueError:
            pass

    def _remove_selected(self):
        for item in self.list_widget.selectedItems():
            self.selected_ids.discard(int(item.text()))
        self._refresh_list_widget()

    def _clear_all(self):
        self.selected_ids.clear()
        self._refresh_list_widget()
