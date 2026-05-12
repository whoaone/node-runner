"""v3.5.0 (item 2) Load Combination tests.

Covers:
- LoadCombinationDialog modes (create / edit / copy)
- result_payload roundtrip
- EditLoadCombinationCommand undo / redo
- CopyLoadCombinationCommand undo / redo
- SID-rename via EditLoadCombinationCommand
"""

import pytest


@pytest.fixture
def qapp():
    from PySide6.QtWidgets import QApplication
    import sys
    app = QApplication.instance() or QApplication(sys.argv)
    yield app


@pytest.fixture
def small_model_with_combos():
    """Build a tiny BDF with 3 load sets + 2 combinations."""
    from pyNastran.bdf.bdf import BDF
    m = BDF(debug=False)
    m.add_grid(1, [0.0, 0.0, 0.0])
    m.add_grid(2, [1.0, 0.0, 0.0])
    # Three load sets.
    m.add_force(sid=100, node=1, mag=1.0, xyz=[1.0, 0.0, 0.0])
    m.add_force(sid=200, node=2, mag=2.0, xyz=[0.0, 1.0, 0.0])
    m.add_force(sid=300, node=1, mag=3.0, xyz=[0.0, 0.0, 1.0])
    # Two combinations.
    m.load_combinations = {
        1000: {'scale': 1.0, 'scale_factors': [1.0, 0.5],
               'load_ids': [100, 200]},
        2000: {'scale': 2.0, 'scale_factors': [1.0],
               'load_ids': [300]},
    }
    return m


# ----------------------------------------------------------------------
# Dialog modes + payload
# ----------------------------------------------------------------------

def test_edit_mode_populates_initial_state(qapp, small_model_with_combos):
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        small_model_with_combos,
        mode=LoadCombinationDialog.MODE_EDIT,
        combo_sid=1000)
    sid, payload = dlg.result_payload()
    assert sid == 1000
    assert payload['scale'] == 1.0
    assert payload['scale_factors'] == [1.0, 0.5]
    assert payload['load_ids'] == [100, 200]
    # SID input must be read-only in edit mode.
    assert dlg._sid_input.isEnabled() is False


def test_copy_mode_advances_sid(qapp, small_model_with_combos):
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        small_model_with_combos,
        mode=LoadCombinationDialog.MODE_COPY,
        combo_sid=1000)
    sid, payload = dlg.result_payload()
    # New SID must be max(existing) + 1 = max(100,200,300,1000,2000)+1 = 2001
    assert sid == 2001
    # But payload mirrors the source.
    assert payload['scale_factors'] == [1.0, 0.5]
    assert payload['load_ids'] == [100, 200]
    # SID input must be editable in copy mode.
    assert dlg._sid_input.isEnabled() is True


def test_create_mode_starts_empty(qapp, small_model_with_combos):
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        small_model_with_combos,
        mode=LoadCombinationDialog.MODE_CREATE)
    sid, payload = dlg.result_payload()
    assert sid > 0
    assert payload['scale_factors'] == []
    assert payload['load_ids'] == []


def test_add_remove_member_row(qapp, small_model_with_combos):
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        small_model_with_combos,
        mode=LoadCombinationDialog.MODE_EDIT,
        combo_sid=1000)
    assert dlg._table.rowCount() == 2
    dlg._add_row(0.75, 300)
    assert dlg._table.rowCount() == 3
    sid, payload = dlg.result_payload()
    assert payload['scale_factors'] == [1.0, 0.5, 0.75]
    assert payload['load_ids'] == [100, 200, 300]


# ----------------------------------------------------------------------
# Edit / Copy / Rename commands
# ----------------------------------------------------------------------

def test_edit_command_roundtrip(small_model_with_combos):
    from node_runner.commands import EditLoadCombinationCommand
    m = small_model_with_combos
    cmd = EditLoadCombinationCommand(
        1000, 1000,
        {'scale': 3.0, 'scale_factors': [2.0], 'load_ids': [300]})
    cmd.execute(m)
    assert m.load_combinations[1000]['scale'] == 3.0
    assert m.load_combinations[1000]['load_ids'] == [300]
    cmd.undo(m)
    assert m.load_combinations[1000]['scale'] == 1.0
    assert m.load_combinations[1000]['load_ids'] == [100, 200]


def test_edit_command_can_change_sid(small_model_with_combos):
    from node_runner.commands import EditLoadCombinationCommand
    m = small_model_with_combos
    cmd = EditLoadCombinationCommand(
        1000, 1500,
        {'scale': 1.0, 'scale_factors': [1.0, 0.5],
         'load_ids': [100, 200]})
    cmd.execute(m)
    assert 1500 in m.load_combinations
    assert 1000 not in m.load_combinations  # renamed
    cmd.undo(m)
    assert 1500 not in m.load_combinations
    assert 1000 in m.load_combinations  # restored


def test_copy_command_leaves_source(small_model_with_combos):
    from node_runner.commands import CopyLoadCombinationCommand
    m = small_model_with_combos
    cmd = CopyLoadCombinationCommand(
        9000, {'scale': 4.0, 'scale_factors': [1.0],
               'load_ids': [100]})
    cmd.execute(m)
    assert 9000 in m.load_combinations
    assert m.load_combinations[9000]['scale'] == 4.0
    # Source intact.
    assert m.load_combinations[1000]['scale'] == 1.0
    cmd.undo(m)
    assert 9000 not in m.load_combinations
    assert m.load_combinations[1000]['scale'] == 1.0
