"""v5.2.0 regression tests for the Post-Processing Toolbox.

Covers:
* The three collapsible sections exist with the expected titles.
* ``populate_from_results`` immediately drives ``_on_output_set_changed``
  on the first item so ``current_subcase_id`` and the Output Vector
  combo both light up without requiring the user to manually toggle
  the Output Set dropdown (v5.2.0 Round-2 Bug #1).
* Output Vector combo is grouped per kind (Displacement / Stress /
  Eigenvectors) with category headers.
* Switching the Output Set dropdown rebuilds the Output Vector combo
  for the new subcase.
* Animation "Through Modes" radio is enabled only when both Animate
  style AND eigenvector data are present.
* Collapsible sections actually collapse.
"""

from __future__ import annotations

import pytest


BUNDLE_STATIC_ONLY = {
    'subcases': {
        1: {'title': 'Static',
            'displacements': {1: [0.1, 0.0, 0.0, 0, 0, 0]}},
    },
}

BUNDLE_MODES_ONLY = {
    'subcases': {
        2: {'title': 'Modes',
            'eigenvectors': [{1: [0.5, 0, 0, 0, 0, 0]}],
            'frequencies': [12.34]},
    },
}

BUNDLE_BOTH = {
    'subcases': {
        1: {'title': 'Static',
            'displacements': {1: [0.1, 0.0, 0.0, 0, 0, 0]}},
        2: {'title': 'Modes',
            'eigenvectors': [{1: [0.5, 0, 0, 0, 0, 0]}],
            'frequencies': [12.34]},
    },
}


@pytest.fixture
def tab(qtbot):
    from node_runner.dialogs.results_tab import ResultsTab
    t = ResultsTab(None)
    qtbot.addWidget(t)
    return t


def test_three_sections_present(tab):
    titles = [s.title for s in tab.sections]
    assert titles == ['Results', 'Deform', 'Contour']


def test_collapsible_section_toggles(tab):
    sec = tab.sections[0]
    assert sec.is_expanded
    sec.set_expanded(False)
    assert not sec.is_expanded
    sec.set_expanded(True)
    assert sec.is_expanded


def test_populate_triggers_first_subcase(tab):
    """v5.2.0 Round-2 Bug #1: populate_from_results must immediately
    drive _on_output_set_changed so current_subcase_id != None and the
    Output Vector combo is populated."""
    tab._find_bundle_from_parent = lambda: BUNDLE_STATIC_ONLY
    tab.populate_from_results(BUNDLE_STATIC_ONLY, file_label='fake')
    assert tab.current_subcase_id() == 1
    # Output vector combo must have items (8 displacement components +
    # the "── Displacement ──" header).
    assert tab.output_vector_combo.count() > 0
    # First selectable entry should be Displacement - Magnitude.
    kind, comp, mode_idx = tab.current_output_vector()
    assert kind == 'displacement'
    assert comp == 'Magnitude'
    assert mode_idx == -1


def test_populate_with_modes_only(tab):
    """SOL 103-style bundle (only eigenvectors) populates Mode groups."""
    tab._find_bundle_from_parent = lambda: BUNDLE_MODES_ONLY
    tab.populate_from_results(BUNDLE_MODES_ONLY, file_label='fake')
    assert tab.current_subcase_id() == 2
    assert tab._has_eigenvectors is True
    labels = [tab.output_vector_combo.itemText(i)
              for i in range(tab.output_vector_combo.count())]
    # At least one Mode group header + components.
    assert any('Mode 1' in l for l in labels)
    assert any('12.3 Hz' in l for l in labels)


def test_switching_output_set_rebuilds_vectors(tab):
    tab._find_bundle_from_parent = lambda: BUNDLE_BOTH
    tab.populate_from_results(BUNDLE_BOTH, file_label='fake')
    # Subcase 1 is Static -> Displacement components, no Mode groups.
    labels = [tab.output_vector_combo.itemText(i)
              for i in range(tab.output_vector_combo.count())]
    assert not any('Mode' in l for l in labels)
    # Toggle to Subcase 2.
    tab.output_set_combo.setCurrentIndex(1)
    labels = [tab.output_vector_combo.itemText(i)
              for i in range(tab.output_vector_combo.count())]
    assert any('Mode 1' in l for l in labels)


def test_through_modes_radio_enablement(tab):
    """Through Modes is greyed off without eigenvectors, and stays
    greyed when Style != Animate even if eigenvectors are present."""
    # Initial: no bundle -> Through Modes disabled.
    assert not tab._anim_mode_modes.isEnabled()
    # Load static-only bundle, switch to Animate -> still disabled (no eigs).
    tab._find_bundle_from_parent = lambda: BUNDLE_STATIC_ONLY
    tab.populate_from_results(BUNDLE_STATIC_ONLY, file_label='fake')
    tab.deform_style_combo.setCurrentIndex(2)   # 'animate'
    assert not tab._anim_mode_modes.isEnabled()
    # Load a bundle with eigenvectors -> Through Modes enabled
    # regardless of current Deform Style. v5.2 fix: Animation Options
    # spinboxes + radios are now always editable so the user can
    # pre-configure animation parameters before picking Animate.
    tab._find_bundle_from_parent = lambda: BUNDLE_MODES_ONLY
    tab.populate_from_results(BUNDLE_MODES_ONLY, file_label='fake')
    tab._update_animation_options_enabled()
    assert tab._anim_mode_modes.isEnabled()
    # Switching back to Deformed leaves Through Modes available since
    # eigenvectors are still loaded.
    tab.deform_style_combo.setCurrentIndex(1)   # 'deformed'
    assert tab._anim_mode_modes.isEnabled()


def test_clear_results(tab):
    """Calling populate_from_results(None) should leave the combos
    clear / disabled."""
    tab.populate_from_results(None)
    assert tab.current_subcase_id() is None
    assert not tab.output_set_combo.isEnabled()


def test_data_conversion_signal_emits(qtbot, tab):
    """Picking a different conversion radio emits data_conversion_changed."""
    with qtbot.waitSignal(tab.data_conversion_changed, timeout=1000) as blocker:
        tab._conv_buttons['max_node'].setChecked(True)
    assert blocker.args == ['max_node']


def test_populate_via_mainwindow_no_recursion(qtbot):
    """Regression for v5.2.0 Round-2 Bug #4: when populate_from_results
    is called through MainWindow (which auto-engages Results color
    mode), the new toolbox emitting output_vector_changed used to lead
    back into MainWindow._on_results_changed -> populate_from_results
    -> infinite recursion.
    """
    from node_runner.mainwindow import MainWindow
    mw = MainWindow()
    qtbot.addWidget(mw)
    mw.op2_results = {
        'subcases': {
            1: {'title': 'Static',
                'displacements': {1: [0.1, 0, 0, 0, 0, 0]}},
        },
    }
    # If recursion is back this would hit Python's recursion limit
    # and crash the test runner.
    mw.results_tab.populate_from_results(
        mw.op2_results, file_label='regression')
    # Toolbox latched the subcase + selected a default vector.
    assert mw.results_tab.current_subcase_id() == 1
    kind, _comp, _mode = mw.results_tab.current_output_vector()
    assert kind == 'displacement'
    mw.close()
