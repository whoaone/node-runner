"""Tests for the v3.3.0 Femap-style signed-bucket EntitySelectionBar.

The bucket is a list of (mode, start, end, step) entries. Each entry's
label is signed:
    +  add range/ID  -> "1, 10, 2"   or "12345"
    -  remove        -> "-1, -10, 1" or "-8"
    x  exclude       -> "x5, x7, 1"  or "x100"
Final selection = walk entries in order, expand each range, apply
set arithmetic (+ unions, - / x subtract).
"""

import pytest


@pytest.fixture
def qapp():
    """Create a QApplication once for the whole module."""
    from PySide6.QtWidgets import QApplication
    import sys
    app = QApplication.instance() or QApplication(sys.argv)
    yield app


@pytest.fixture
def bar(qapp):
    """Fresh EntitySelectionBar configured with entity ids 1..100."""
    from node_runner.dialogs.selection import EntitySelectionBar
    b = EntitySelectionBar()
    b.configure('Element', set(range(1, 101)))
    return b


# ----------------------------------------------------------------------
# Label formatting
# ----------------------------------------------------------------------

def test_range_entry_label_single_id():
    from node_runner.dialogs.selection import _RangeEntry
    e_add = _RangeEntry('+', 12345, 12345, 1)
    e_rem = _RangeEntry('-', 12345, 12345, 1)
    e_exc = _RangeEntry('x', 12345, 12345, 1)
    assert e_add.label == "12345"
    assert e_rem.label == "-12345"
    assert e_exc.label == "x12345"


def test_range_entry_label_range_step1():
    from node_runner.dialogs.selection import _RangeEntry
    assert _RangeEntry('+', 1, 10, 1).label == "1, 10, 1"
    assert _RangeEntry('-', 1, 10, 1).label == "-1, -10, 1"
    assert _RangeEntry('x', 1, 10, 1).label == "x1, x10, 1"


def test_range_entry_label_range_step2():
    from node_runner.dialogs.selection import _RangeEntry
    assert _RangeEntry('+', 1, 10, 2).label == "1, 10, 2"
    assert _RangeEntry('-', 1, 10, 2).label == "-1, -10, 2"


def test_range_entry_expand_step2():
    from node_runner.dialogs.selection import _RangeEntry
    e = _RangeEntry('+', 1, 9, 2)
    assert sorted(e.ids) == [1, 3, 5, 7, 9]


# ----------------------------------------------------------------------
# Final-set arithmetic
# ----------------------------------------------------------------------

def test_add_then_remove_subtracts_id(bar):
    """User spec: 'i have ids added 1 - 10 and then remove id 8, my
    set won't include id 8.'"""
    bar._append_entry('+', 1, 10, 1)
    bar._append_entry('-', 8, 8, 1)
    assert bar.get_selected_ids() == [1, 2, 3, 4, 5, 6, 7, 9, 10]


def test_exclude_subtracts_range(bar):
    bar._append_entry('+', 1, 10, 1)
    bar._append_entry('x', 5, 7, 1)
    assert bar.get_selected_ids() == [1, 2, 3, 4, 8, 9, 10]


def test_add_after_remove_re_includes(bar):
    """Order matters: add-after-remove brings back the IDs."""
    bar._append_entry('+', 1, 10, 1)
    bar._append_entry('-', 5, 7, 1)
    bar._append_entry('+', 5, 5, 1)
    assert bar.get_selected_ids() == [1, 2, 3, 4, 5, 8, 9, 10]


def test_delete_row_undoes_pick(bar):
    """Delete the remove-row in the middle and the final set should
    be as if that pick had never happened."""
    bar._append_entry('+', 1, 10, 1)
    bar._append_entry('-', 8, 8, 1)
    assert bar.get_selected_ids() == [1, 2, 3, 4, 5, 6, 7, 9, 10]

    # Highlight just row 1 (the -8 row) and Delete.
    bar._list_widget.clearSelection()
    bar._list_widget.item(1).setSelected(True)
    bar._delete_selected_entries()
    assert bar.get_selected_ids() == list(range(1, 11))


def test_bucket_display_signs(bar):
    """Visible labels in the QListWidget use the Femap-style signed
    triplet (or short single-ID form)."""
    bar._append_entry('+', 1, 10, 2)   # "1, 10, 2"
    bar._append_entry('-', 8, 8, 1)    # "-8"
    bar._append_entry('x', 100, 110, 1)  # "x100, x110, 1"
    bar._append_entry('+', 42, 42, 1)   # "42"
    labels = [bar._list_widget.item(i).text()
              for i in range(bar._list_widget.count())]
    assert labels == ["1, 10, 2", "-8", "x100, x110, 1", "42"]


# ----------------------------------------------------------------------
# Range compression helper (used by Pick / Select Visible / Select All)
# ----------------------------------------------------------------------

def test_compress_ids_to_range_tuples_contiguous():
    from node_runner.dialogs.selection import EntitySelectionBar
    assert (EntitySelectionBar.compress_ids_to_range_tuples([1, 2, 3, 4, 5])
            == [(1, 5, 1)])


def test_compress_ids_to_range_tuples_step2():
    from node_runner.dialogs.selection import EntitySelectionBar
    assert (EntitySelectionBar.compress_ids_to_range_tuples([1, 3, 5, 7, 9])
            == [(1, 9, 2)])


def test_compress_ids_to_range_tuples_split_runs():
    from node_runner.dialogs.selection import EntitySelectionBar
    assert (EntitySelectionBar.compress_ids_to_range_tuples(
        [1, 2, 3, 5, 7, 9, 20, 21, 22])
            == [(1, 3, 1), (5, 9, 2), (20, 22, 1)])


def test_compress_ids_to_range_tuples_single():
    from node_runner.dialogs.selection import EntitySelectionBar
    assert (EntitySelectionBar.compress_ids_to_range_tuples([42])
            == [(42, 42, 1)])


# ----------------------------------------------------------------------
# Bulk operations decompose into per-range entries
# ----------------------------------------------------------------------

def test_select_all_emits_compressed_range(bar):
    """Select All on entities 1..100 produces a single range entry
    '1, 100, 1' rather than 100 individual rows."""
    bar._select_all()
    labels = [bar._list_widget.item(i).text()
              for i in range(bar._list_widget.count())]
    assert labels == ["1, 100, 1"]
    assert bar.get_selected_ids() == list(range(1, 101))


def test_append_id_set_skips_dangling(bar):
    """IDs outside all_entity_ids never leak into the final selection."""
    bar._append_id_set('+', {50, 200, 300})  # 200, 300 are dangling
    assert bar.get_selected_ids() == [50]


# ----------------------------------------------------------------------
# v3.4.0 - Phase A regressions
# ----------------------------------------------------------------------

def test_single_id_more_creates_one_entry(bar):
    """Item 3: typing one ID with To/By blank and clicking More must
    create a single-ID bucket row, not a range or nothing."""
    bar._id_input.setText("12345")
    # Configure has set all_entity_ids = 1..100 by default, so 12345
    # is out of range; widen the entity set for this test.
    bar.all_entity_ids = set(range(1, 100001))
    bar._on_more_clicked()
    labels = [bar._list_widget.item(i).text()
              for i in range(bar._list_widget.count())]
    assert labels == ["12345"]
    assert bar.get_selected_ids() == [12345]


def test_single_id_more_with_stray_whitespace(bar):
    """Whitespace around the ID should not break the single-ID path."""
    bar._id_input.setText("  42  ")
    bar.all_entity_ids = set(range(1, 1000))
    bar._on_more_clicked()
    labels = [bar._list_widget.item(i).text()
              for i in range(bar._list_widget.count())]
    assert labels == ["42"]


def test_paste_via_internal_helper_into_bucket(qapp):
    """Item 1+2: clipboard text in Femap range format gets parsed and
    added to the bucket. This exercises the same code path that the
    Ctrl+V QShortcut and the Pick > Paste menu trigger."""
    from PySide6.QtWidgets import QApplication
    from node_runner.dialogs.selection import EntitySelectionBar

    bar = EntitySelectionBar()
    bar.configure('Element', set(range(1, 1001)))
    QApplication.clipboard().setText("1,10,2\n42\n100,110,1")
    bar._paste_selection()
    # Three rows should appear (one per clipboard line that parsed).
    labels = [bar._list_widget.item(i).text()
              for i in range(bar._list_widget.count())]
    # The paste path collapses everything into a single id-set, then
    # decomposes via compress_ids_to_range_tuples. So the rows may
    # merge or split based on contiguity. We just assert the final
    # selected ids are what we expect.
    expected = set(range(1, 11, 2)) | {42} | set(range(100, 111, 1))
    assert set(bar.get_selected_ids()) == expected
