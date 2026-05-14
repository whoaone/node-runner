"""v4.0.0 (B1) tree-state snapshot/restore tests.

The QTreeWidgets in the sidebar tabs get rebuilt after every small
edit (add SPC, edit material, etc.). Without state preservation,
the user's checkbox visibility selections and expand state were wiped
on every rebuild. These tests pin the snapshot/restore round-trip.
"""

import pytest


@pytest.fixture
def qapp():
    from PySide6.QtWidgets import QApplication
    import sys
    app = QApplication.instance() or QApplication(sys.argv)
    yield app


def _make_tree_with_items(items):
    """Build a fresh QTreeWidget from [(label, userrole_key, checked, expanded), ...]."""
    from PySide6 import QtCore
    from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem
    tw = QTreeWidget()
    for label, key, checked, expanded in items:
        item = QTreeWidgetItem(tw, [label])
        item.setData(0, QtCore.Qt.UserRole, key)
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(
            0, QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
        item.setExpanded(expanded)
    return tw


def test_snapshot_captures_check_and_expand(qapp):
    from node_runner.mainwindow import MainWindow
    from PySide6 import QtCore
    tw = _make_tree_with_items([
        ("A", ('nodes',), True, True),
        ("B", ('elements',), False, False),
    ])
    snap = MainWindow._snapshot_tree_state(tw)
    assert snap[('nodes',)]['checked'] == QtCore.Qt.Checked
    assert snap[('nodes',)]['expanded'] is True
    assert snap[('elements',)]['checked'] == QtCore.Qt.Unchecked
    assert snap[('elements',)]['expanded'] is False


def test_restore_applies_snapshot(qapp):
    from node_runner.mainwindow import MainWindow
    from PySide6 import QtCore
    # Original tree: A checked, B unchecked.
    tw1 = _make_tree_with_items([
        ("A", ('nodes',), True, True),
        ("B", ('elements',), False, False),
    ])
    snap = MainWindow._snapshot_tree_state(tw1)
    # Rebuilt tree: defaults flipped.
    tw2 = _make_tree_with_items([
        ("A", ('nodes',), False, False),
        ("B", ('elements',), True, True),
    ])
    MainWindow._restore_tree_state(tw2, snap)
    # Walk and assert restored.
    root = tw2.invisibleRootItem()
    by_key = {}
    for i in range(root.childCount()):
        c = root.child(i)
        by_key[c.data(0, QtCore.Qt.UserRole)] = c
    assert by_key[('nodes',)].checkState(0) == QtCore.Qt.Checked
    assert by_key[('nodes',)].isExpanded() is True
    assert by_key[('elements',)].checkState(0) == QtCore.Qt.Unchecked


def test_added_items_keep_defaults(qapp):
    """Items present in the rebuild but absent from snapshot keep the
    state the rebuild gave them (no crash, no overwrite to defaults)."""
    from node_runner.mainwindow import MainWindow
    from PySide6 import QtCore
    tw1 = _make_tree_with_items([("A", ('nodes',), True, True)])
    snap = MainWindow._snapshot_tree_state(tw1)
    tw2 = _make_tree_with_items([
        ("A", ('nodes',), False, False),    # restore will flip back to True
        ("NEW", ('new_section',), True, True),  # absent from snapshot
    ])
    MainWindow._restore_tree_state(tw2, snap)
    root = tw2.invisibleRootItem()
    by_key = {root.child(i).data(0, QtCore.Qt.UserRole): root.child(i)
              for i in range(root.childCount())}
    # 'nodes' restored from snapshot
    assert by_key[('nodes',)].checkState(0) == QtCore.Qt.Checked
    # 'new_section' kept its rebuild default
    assert by_key[('new_section',)].checkState(0) == QtCore.Qt.Checked


def test_removed_items_are_silently_dropped(qapp):
    """Snapshot keys not present in the rebuilt tree are ignored."""
    from node_runner.mainwindow import MainWindow
    tw1 = _make_tree_with_items([
        ("A", ('nodes',), True, True),
        ("GONE", ('removed_section',), True, True),
    ])
    snap = MainWindow._snapshot_tree_state(tw1)
    tw2 = _make_tree_with_items([("A", ('nodes',), False, False)])
    # Must not raise.
    MainWindow._restore_tree_state(tw2, snap)


def test_unkeyed_items_are_skipped(qapp):
    """Items with no UserRole key cannot be matched and are skipped."""
    from PySide6 import QtCore
    from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem
    from node_runner.mainwindow import MainWindow
    tw = QTreeWidget()
    QTreeWidgetItem(tw, ["unkeyed"])  # no setData call
    snap = MainWindow._snapshot_tree_state(tw)
    assert snap == {}


def test_nested_items_round_trip(qapp):
    """Recurse into children."""
    from PySide6 import QtCore
    from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem
    from node_runner.mainwindow import MainWindow

    def make_tree():
        tw = QTreeWidget()
        root = QTreeWidgetItem(tw, ["root"])
        root.setData(0, QtCore.Qt.UserRole, ('root',))
        root.setFlags(root.flags() | QtCore.Qt.ItemIsUserCheckable)
        root.setCheckState(0, QtCore.Qt.Checked)
        child = QTreeWidgetItem(root, ["child"])
        child.setData(0, QtCore.Qt.UserRole, ('child',))
        child.setFlags(child.flags() | QtCore.Qt.ItemIsUserCheckable)
        child.setCheckState(0, QtCore.Qt.Unchecked)
        root.setExpanded(True)
        return tw

    tw1 = make_tree()
    snap = MainWindow._snapshot_tree_state(tw1)
    assert ('root',) in snap
    assert ('child',) in snap
    assert snap[('child',)]['checked'] == QtCore.Qt.Unchecked

    # Rebuild with flipped states and restore.
    tw2 = QTreeWidget()
    r = QTreeWidgetItem(tw2, ["root"])
    r.setData(0, QtCore.Qt.UserRole, ('root',))
    r.setFlags(r.flags() | QtCore.Qt.ItemIsUserCheckable)
    r.setCheckState(0, QtCore.Qt.Unchecked)  # opposite of original
    c = QTreeWidgetItem(r, ["child"])
    c.setData(0, QtCore.Qt.UserRole, ('child',))
    c.setFlags(c.flags() | QtCore.Qt.ItemIsUserCheckable)
    c.setCheckState(0, QtCore.Qt.Checked)  # opposite of original

    MainWindow._restore_tree_state(tw2, snap)
    assert r.checkState(0) == QtCore.Qt.Checked
    assert c.checkState(0) == QtCore.Qt.Unchecked
