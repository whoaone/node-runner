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


# ----------------------------------------------------------------------
# v4.0.0 (A1): read_combo_payload adapter + real-parse fixtures
#
# After read_bdf parses a deck containing LOAD cards, the slot
# `model.load_combinations[sid]` is `list[LOAD]`, not a dict. v3.5.0
# field-tested as crashing the Edit dialog because the code called
# `.get('scale', 1.0)` on a list. These tests pin the adapter's
# behavior across both shapes.
# ----------------------------------------------------------------------

@pytest.fixture
def parsed_model_with_combos(tmp_path):
    """Real `read_bdf` round-trip producing list[LOAD] storage."""
    from pyNastran.bdf.bdf import read_bdf
    bdf_text = (
        "GRID,1,,0.0,0.0,0.0\n"
        "GRID,2,,1.0,0.0,0.0\n"
        "FORCE,100,1,,1.0,1.0,0.0,0.0\n"
        "FORCE,200,2,,2.0,0.0,1.0,0.0\n"
        "FORCE,300,1,,3.0,0.0,0.0,1.0\n"
        "LOAD,1000,1.0,1.0,100,0.5,200\n"
        "LOAD,2000,2.0,1.0,300\n"
    )
    p = tmp_path / "combos.bdf"
    p.write_text(bdf_text)
    m = read_bdf(str(p), xref=False, validate=False, punch=True)
    return m


def test_adapter_handles_list_from_real_parse(parsed_model_with_combos):
    from node_runner.dialogs.load_combination import read_combo_payload
    m = parsed_model_with_combos
    # pyNastran stores parsed LOADs as a list.
    assert isinstance(m.load_combinations[1000], list)
    payload = read_combo_payload(m.load_combinations[1000])
    assert payload['scale'] == 1.0
    assert payload['scale_factors'] == [1.0, 0.5]
    assert payload['load_ids'] == [100, 200]


def test_adapter_handles_dict_from_in_session():
    from node_runner.dialogs.load_combination import read_combo_payload
    payload = read_combo_payload(
        {'scale': 2.5, 'scale_factors': [1.0, 0.5],
         'load_ids': [10, 20], 'title': 'demo'})
    assert payload['scale'] == 2.5
    assert payload['scale_factors'] == [1.0, 0.5]
    assert payload['load_ids'] == [10, 20]
    assert payload['title'] == 'demo'


def test_adapter_handles_empty_fallback():
    from node_runner.dialogs.load_combination import read_combo_payload
    payload = read_combo_payload(None)
    assert payload == {
        'scale': 1.0, 'scale_factors': [], 'load_ids': [], 'title': ''}


def test_edit_dialog_opens_on_parsed_combo(qapp, parsed_model_with_combos):
    """Regression: v3.5.0 dialog raised AttributeError on .get() of list."""
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        parsed_model_with_combos,
        mode=LoadCombinationDialog.MODE_EDIT,
        combo_sid=1000)
    sid, payload = dlg.result_payload()
    assert sid == 1000
    assert payload['scale'] == 1.0
    assert payload['load_ids'] == [100, 200]


def test_copy_dialog_opens_on_parsed_combo(qapp, parsed_model_with_combos):
    """Regression: MODE_COPY also did `dict(list_value)` which raised."""
    from node_runner.dialogs.load_combination import LoadCombinationDialog
    dlg = LoadCombinationDialog(
        parsed_model_with_combos,
        mode=LoadCombinationDialog.MODE_COPY,
        combo_sid=2000)
    sid, payload = dlg.result_payload()
    assert sid > 2000  # new SID assigned
    assert payload['scale_factors'] == [1.0]
    assert payload['load_ids'] == [300]


def test_edit_command_executes_on_parsed_combo(parsed_model_with_combos):
    """Regression: EditLoadCombinationCommand.execute called dict() on
    the old value to snapshot for undo; that raised on list[LOAD]."""
    from node_runner.commands import EditLoadCombinationCommand
    m = parsed_model_with_combos
    cmd = EditLoadCombinationCommand(
        1000, 1000,
        {'scale': 3.0, 'scale_factors': [2.0], 'load_ids': [300]})
    cmd.execute(m)  # was raising TypeError pre-fix
    assert m.load_combinations[1000]['scale'] == 3.0
    cmd.undo(m)
    # Undo restores dict-shape (adapter normalized snapshot).
    restored = m.load_combinations[1000]
    assert restored['scale'] == 1.0
    assert restored['load_ids'] == [100, 200]


def test_mixed_shape_model_displays_consistently(parsed_model_with_combos):
    """Build a model with a parsed combo + an in-session combo and
    verify the adapter normalizes both to identical payload shapes."""
    from node_runner.dialogs.load_combination import read_combo_payload
    m = parsed_model_with_combos
    m.load_combinations[3000] = {
        'scale': 1.0, 'scale_factors': [1.0], 'load_ids': [100],
        'title': 'in-session'}
    parsed = read_combo_payload(m.load_combinations[1000])
    in_session = read_combo_payload(m.load_combinations[3000])
    assert set(parsed) == set(in_session)
    assert all(isinstance(parsed[k], type(in_session[k]))
               for k in parsed)
