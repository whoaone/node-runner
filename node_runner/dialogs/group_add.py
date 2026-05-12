"""v3.4.0 (item 7): Femap-style tabbed Add-to-Group / Remove-from-Group
dialog.

The dialog lets the user push a set of IDs into a group's per-entity
list via a rule-based selector. One tab per entity type (Node /
Element / Property / Material / Coord). Each tab presents a small set
of selectors appropriate for that entity:

    Node      : By ID range, By connected to element, From selection, By group
    Element   : By ID range, By Property, By Material, By Type,
                From selection, By group
    Property  : By ID range, By Material, From selection, By group
    Material  : By ID range, From selection, By group
    Coord     : By ID range, From selection, By group

Selectors return an integer set; the caller merges it into the
target group's entity list (Add mode) or subtracts it (Remove mode).

The dialog is purely a UI shell - all the actual entity lookups go
through callback hooks the parent (MainWindow) provides at construction
time. That keeps the dialog free of pyNastran knowledge and easy to
test.
"""

from __future__ import annotations

from PySide6 import QtCore
from PySide6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QLineEdit, QComboBox, QRadioButton, QButtonGroup, QFrame, QTabWidget,
    QStackedWidget, QSizePolicy,
)


# ----------------------------------------------------------------------
# Selector page (one entity tab is a stack of these, switched by radio)
# ----------------------------------------------------------------------


def _range_widget():
    """ID / To / By inputs returning a (start, end, step) tuple. Returns
    the QWidget and a getter callable."""
    w = QWidget()
    lay = QHBoxLayout(w); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
    lay.addWidget(QLabel("ID"))
    id_in = QLineEdit(); id_in.setMaximumWidth(70); lay.addWidget(id_in)
    lay.addWidget(QLabel("To:"))
    to_in = QLineEdit(); to_in.setMaximumWidth(70); lay.addWidget(to_in)
    lay.addWidget(QLabel("By"))
    by_in = QLineEdit(); by_in.setMaximumWidth(50)
    by_in.setPlaceholderText("1"); lay.addWidget(by_in)
    lay.addStretch(1)

    def parse():
        try:
            s = int(id_in.text().strip())
        except ValueError:
            return set()
        to_text = to_in.text().strip()
        by_text = by_in.text().strip() or "1"
        try:
            step = max(1, int(by_text))
        except ValueError:
            step = 1
        if not to_text:
            return {s}
        try:
            e = int(to_text)
        except ValueError:
            return {s}
        lo, hi = min(s, e), max(s, e)
        return set(range(lo, hi + 1, step))

    return w, parse


def _combo_widget(label_text, items):
    """A labeled combo whose ``.currentData`` returns the picked id-set.

    ``items`` is a list of ``(label, ids_set)`` tuples."""
    w = QWidget()
    lay = QHBoxLayout(w); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
    lay.addWidget(QLabel(f"{label_text}:"))
    combo = QComboBox()
    for label, ids in items:
        combo.addItem(label, ids)
    combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    lay.addWidget(combo, 1)

    def parse():
        data = combo.currentData()
        return set(data) if data else set()

    return w, parse


def _info_widget(text):
    w = QLabel(text)
    w.setStyleSheet("color: #888; font-style: italic;")
    return w


class _EntityTab(QWidget):
    """One tab of the dialog. Owns a radio group of selectors and a
    QStackedWidget of input pages, one per selector."""

    def __init__(self, entity_type: str, selectors: list[dict], parent=None):
        super().__init__(parent)
        self.entity_type = entity_type
        self._parsers = []  # list[callable -> set[int]]

        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # Header
        lay.addWidget(QLabel("Selector:"))
        self._radio_group = QButtonGroup(self)
        radio_row = QVBoxLayout(); radio_row.setSpacing(2)
        self._stack = QStackedWidget()
        for i, sel in enumerate(selectors):
            rb = QRadioButton(sel['label'])
            self._radio_group.addButton(rb, i)
            radio_row.addWidget(rb)
            if i == 0:
                rb.setChecked(True)
            self._stack.addWidget(sel['widget'])
            self._parsers.append(sel['parse'])
        lay.addLayout(radio_row)
        lay.addWidget(self._stack)

        # Live preview label
        self._preview = QLabel("")
        self._preview.setStyleSheet("color: #94e2d5; font-style: italic;")
        lay.addWidget(self._preview)

        lay.addStretch(1)

        # Wire radio to stack + preview refresh
        self._radio_group.buttonToggled.connect(self._on_radio_toggled)

    def _on_radio_toggled(self, btn, checked):
        if checked:
            self._stack.setCurrentIndex(self._radio_group.id(btn))
        self.refresh_preview()

    def refresh_preview(self):
        try:
            n = len(self.current_ids())
        except Exception:
            n = 0
        self._preview.setText(f"({n} {self.entity_type.lower()}s will be applied)")

    def current_ids(self) -> set[int]:
        idx = self._radio_group.checkedId()
        if idx < 0 or idx >= len(self._parsers):
            return set()
        return self._parsers[idx]() or set()


# ----------------------------------------------------------------------
# Public dialog
# ----------------------------------------------------------------------


class GroupAddDialog(QDialog):
    """Tabbed Add (or Remove) entities from a group.

    Construction takes a ``context`` dict supplied by the parent that
    enumerates the model's entities + the current viewport selection +
    the current group registry. Keys (all optional except ``group_name``):

        group_name   : str   (target group)
        mode         : 'add' | 'remove'
        node_ids     : Iterable[int]
        element_ids  : Iterable[int]
        property_ids : Iterable[int]
        material_ids : Iterable[int]
        coord_ids    : Iterable[int]
        eids_by_pid  : dict[int -> set[int]]
        eids_by_mid  : dict[int -> set[int]]
        eids_by_type : dict[str -> set[int]]
        pids_by_mid  : dict[int -> set[int]]
        nids_for_eids: dict[int -> set[int]]   (membership of nodes in elements)
        groups       : dict[name -> group_data]   (for "By group" selectors)
        current_selection_nids : set[int]
        current_selection_eids : set[int]

    The result, retrievable via ``ids_by_entity()``, is a dict mapping
    entity-type keys ('nodes' / 'elements' / 'properties' / 'materials'
    / 'coords') -> set of IDs the user picked across all 5 tabs.
    """

    ENTITY_KEYS = ('nodes', 'elements', 'properties', 'materials', 'coords')

    def __init__(self, context: dict, parent=None):
        super().__init__(parent)
        self.context = context
        mode = context.get('mode', 'add')
        gname = context.get('group_name', '<unnamed>')
        self.setWindowTitle(
            f"{'Add to' if mode == 'add' else 'Remove from'} Group: {gname}")
        self.setMinimumWidth(560)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 8, 8, 8)
        outer.setSpacing(6)

        title = QLabel(
            f"{'Add entities to' if mode == 'add' else 'Remove entities from'} "
            f"group: {gname}")
        title.setStyleSheet(
            "font-weight: bold; font-size: 12px; color: #f9e2af;")
        outer.addWidget(title)

        # Build the 5 tabs.
        self._tabs = QTabWidget()
        self._tab_node = self._build_node_tab()
        self._tab_elem = self._build_element_tab()
        self._tab_prop = self._build_property_tab()
        self._tab_mat = self._build_material_tab()
        self._tab_coord = self._build_coord_tab()
        self._tabs.addTab(self._tab_node, "Node")
        self._tabs.addTab(self._tab_elem, "Element")
        self._tabs.addTab(self._tab_prop, "Property")
        self._tabs.addTab(self._tab_mat, "Material")
        self._tabs.addTab(self._tab_coord, "Coord")
        outer.addWidget(self._tabs)

        # Bottom row
        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        outer.addWidget(sep)
        bot = QHBoxLayout(); bot.addStretch(1)
        self._ok_btn = QPushButton("Add" if mode == 'add' else "Remove")
        self._ok_btn.setStyleSheet("font-weight: bold;")
        self._ok_btn.setMinimumWidth(90)
        self._ok_btn.clicked.connect(self.accept)
        cancel = QPushButton("Cancel"); cancel.setMinimumWidth(80)
        cancel.clicked.connect(self.reject)
        bot.addWidget(self._ok_btn); bot.addWidget(cancel)
        outer.addLayout(bot)

        # Initial preview text on the current tab.
        for t in (self._tab_node, self._tab_elem, self._tab_prop,
                  self._tab_mat, self._tab_coord):
            t.refresh_preview()

    # ------------------------------------------------------------------
    # Tab builders
    # ------------------------------------------------------------------

    def _build_node_tab(self):
        selectors = []
        rw, rp = _range_widget()
        selectors.append({'label': 'By ID range', 'widget': rw, 'parse': rp})
        # By connected to element: use current selection_eids -> nodes
        cw, cp = self._by_connected_widget()
        selectors.append({'label': 'By connected to selected element',
                          'widget': cw, 'parse': cp})
        # From selection
        sw, sp = self._from_selection_widget('nodes')
        selectors.append({'label': 'From current selection',
                          'widget': sw, 'parse': sp})
        # By group
        gw, gp = self._by_group_widget('nodes')
        selectors.append({'label': 'By group', 'widget': gw, 'parse': gp})
        return _EntityTab('Node', selectors, self)

    def _build_element_tab(self):
        selectors = []
        rw, rp = _range_widget()
        selectors.append({'label': 'By ID range', 'widget': rw, 'parse': rp})
        # By Property
        items = [(f"PID {pid}", ids)
                 for pid, ids in sorted(
                     (self.context.get('eids_by_pid') or {}).items())]
        pw, pp = _combo_widget("PID", items) if items else (
            _info_widget("(no properties in model)"), lambda: set())
        selectors.append({'label': 'By Property', 'widget': pw, 'parse': pp})
        # By Material
        items = [(f"MID {mid}", ids)
                 for mid, ids in sorted(
                     (self.context.get('eids_by_mid') or {}).items())]
        mw, mp = _combo_widget("MID", items) if items else (
            _info_widget("(no materials in model)"), lambda: set())
        selectors.append({'label': 'By Material', 'widget': mw, 'parse': mp})
        # By Type
        items = [(f"{etype} ({len(ids)})", ids)
                 for etype, ids in sorted(
                     (self.context.get('eids_by_type') or {}).items())]
        tw, tp = _combo_widget("Type", items) if items else (
            _info_widget("(no elements in model)"), lambda: set())
        selectors.append({'label': 'By Type', 'widget': tw, 'parse': tp})
        # From selection
        sw, sp = self._from_selection_widget('elements')
        selectors.append({'label': 'From current selection',
                          'widget': sw, 'parse': sp})
        # By group
        gw, gp = self._by_group_widget('elements')
        selectors.append({'label': 'By group', 'widget': gw, 'parse': gp})
        return _EntityTab('Element', selectors, self)

    def _build_property_tab(self):
        selectors = []
        rw, rp = _range_widget()
        selectors.append({'label': 'By ID range', 'widget': rw, 'parse': rp})
        items = [(f"MID {mid}", ids)
                 for mid, ids in sorted(
                     (self.context.get('pids_by_mid') or {}).items())]
        mw, mp = _combo_widget("MID", items) if items else (
            _info_widget("(no materials in model)"), lambda: set())
        selectors.append({'label': 'By Material', 'widget': mw, 'parse': mp})
        sw, sp = self._from_selection_widget('properties')
        selectors.append({'label': 'From current selection',
                          'widget': sw, 'parse': sp})
        gw, gp = self._by_group_widget('properties')
        selectors.append({'label': 'By group', 'widget': gw, 'parse': gp})
        return _EntityTab('Property', selectors, self)

    def _build_material_tab(self):
        selectors = []
        rw, rp = _range_widget()
        selectors.append({'label': 'By ID range', 'widget': rw, 'parse': rp})
        sw, sp = self._from_selection_widget('materials')
        selectors.append({'label': 'From current selection',
                          'widget': sw, 'parse': sp})
        gw, gp = self._by_group_widget('materials')
        selectors.append({'label': 'By group', 'widget': gw, 'parse': gp})
        return _EntityTab('Material', selectors, self)

    def _build_coord_tab(self):
        selectors = []
        rw, rp = _range_widget()
        selectors.append({'label': 'By ID range', 'widget': rw, 'parse': rp})
        sw, sp = self._from_selection_widget('coords')
        selectors.append({'label': 'From current selection',
                          'widget': sw, 'parse': sp})
        gw, gp = self._by_group_widget('coords')
        selectors.append({'label': 'By group', 'widget': gw, 'parse': gp})
        return _EntityTab('Coord', selectors, self)

    # ------------------------------------------------------------------
    # Selector widget factories that depend on context
    # ------------------------------------------------------------------

    def _by_connected_widget(self):
        eids = self.context.get('current_selection_eids') or set()
        nids_for_eids = self.context.get('nids_for_eids') or {}
        n_resolved = set()
        for eid in eids:
            n_resolved |= set(nids_for_eids.get(eid, ()))
        info = _info_widget(
            f"(will pull {len(n_resolved)} nodes connected to "
            f"{len(eids)} currently-selected element(s))")

        def parse():
            return set(n_resolved)
        return info, parse

    def _from_selection_widget(self, kind: str):
        """Returns a label + parse function reading the current
        viewport selection. For node/element kinds the parent
        usually has live selection sets; for property/material/coord
        these are derived from whichever entries are highlighted in
        the model tree (caller responsibility)."""
        sel = set()
        if kind == 'nodes':
            sel = set(self.context.get('current_selection_nids') or [])
        elif kind == 'elements':
            sel = set(self.context.get('current_selection_eids') or [])
        elif kind == 'properties':
            sel = set(self.context.get('current_selection_pids') or [])
        elif kind == 'materials':
            sel = set(self.context.get('current_selection_mids') or [])
        elif kind == 'coords':
            sel = set(self.context.get('current_selection_cids') or [])
        info = _info_widget(f"(will pull {len(sel)} {kind} from current selection)")

        def parse():
            return set(sel)
        return info, parse

    def _by_group_widget(self, kind: str):
        groups = self.context.get('groups') or {}
        items = []
        for name, data in sorted(
                groups.items(), key=lambda kv: kv[1].get('gid', 0)):
            ids = set(data.get(kind, []) or [])
            items.append((f"{name} ({len(ids)} {kind})", ids))
        if not items:
            return (_info_widget("(no other groups defined)"),
                    lambda: set())
        return _combo_widget("Group", items)

    # ------------------------------------------------------------------
    # Result
    # ------------------------------------------------------------------

    def ids_by_entity(self) -> dict[str, set[int]]:
        """Return the final id sets per entity type.

        Returns a dict with keys 'nodes', 'elements', 'properties',
        'materials', 'coords'. The set on the currently-active tab
        comes from its current selector; all OTHER tabs return empty
        sets (so the result is unambiguous - only the tab the user
        committed contributes ids)."""
        out = {k: set() for k in self.ENTITY_KEYS}
        idx = self._tabs.currentIndex()
        if idx == 0:
            out['nodes'] = self._tab_node.current_ids()
        elif idx == 1:
            out['elements'] = self._tab_elem.current_ids()
        elif idx == 2:
            out['properties'] = self._tab_prop.current_ids()
        elif idx == 3:
            out['materials'] = self._tab_mat.current_ids()
        elif idx == 4:
            out['coords'] = self._tab_coord.current_ids()
        return out
