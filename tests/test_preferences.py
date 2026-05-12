"""v3.5.0 (item 3) Preferences dialog + QSettings roundtrip tests."""

import pytest


@pytest.fixture
def qapp():
    from PySide6.QtWidgets import QApplication
    import sys
    app = QApplication.instance() or QApplication(sys.argv)
    yield app


@pytest.fixture(autouse=True)
def isolated_settings(tmp_path, monkeypatch):
    """Redirect QSettings to a temp folder so tests don't pollute
    the user's actual saved preferences."""
    from PySide6.QtCore import QSettings, QStandardPaths
    QStandardPaths.setTestModeEnabled(True)
    # Wipe any leftover keys at our org/app namespace.
    qs = QSettings("NodeRunner", "NodeRunner")
    qs.remove("colors")
    qs.remove("sizes")
    yield
    qs.remove("colors")
    qs.remove("sizes")
    QStandardPaths.setTestModeEnabled(False)


def test_load_preferences_returns_defaults_when_empty(qapp):
    from node_runner.dialogs.preferences import (
        load_preferences, DEFAULT_TYPE_COLORS, DEFAULT_HIGHLIGHT_COLOR)
    prefs = load_preferences()
    for key, default in DEFAULT_TYPE_COLORS.items():
        assert prefs['colors'][key] == default
    assert prefs['highlight_color'] == DEFAULT_HIGHLIGHT_COLOR
    assert 'mass_glyph_scale_pct' in prefs['sizes']


def test_save_then_load_roundtrip(qapp):
    from node_runner.dialogs.preferences import (
        load_preferences, save_preferences)
    new_payload = {
        'colors': {'RBE2': '#aabbcc', 'Masses': '#112233'},
        'highlight_color': '#deadbe',
        'sizes': {'mass_glyph_scale_pct': 3.5, 'node_size': 7,
                  'highlight_outline_width': 8},
    }
    save_preferences(new_payload)
    loaded = load_preferences()
    assert loaded['colors']['RBE2'] == '#aabbcc'
    assert loaded['colors']['Masses'] == '#112233'
    assert loaded['highlight_color'] == '#deadbe'
    assert abs(loaded['sizes']['mass_glyph_scale_pct'] - 3.5) < 1e-6
    assert loaded['sizes']['node_size'] == 7
    assert loaded['sizes']['highlight_outline_width'] == 8


def test_dialog_result_payload_reflects_initial_state(qapp):
    from node_runner.dialogs.preferences import (
        PreferencesDialog, DEFAULT_TYPE_COLORS, DEFAULT_HIGHLIGHT_COLOR,
        DEFAULT_SIZES)
    initial = {
        'colors': dict(DEFAULT_TYPE_COLORS),
        'highlight_color': DEFAULT_HIGHLIGHT_COLOR,
        'sizes': dict(DEFAULT_SIZES),
    }
    dlg = PreferencesDialog(initial)
    payload = dlg.result_payload()
    # Defaults survive round-trip.
    assert payload['colors']['RBE2'] == DEFAULT_TYPE_COLORS['RBE2']
    assert payload['highlight_color'] == DEFAULT_HIGHLIGHT_COLOR
    assert (payload['sizes']['mass_glyph_scale_pct']
            == DEFAULT_SIZES['mass_glyph_scale_pct'])


def test_dialog_swatch_changes_propagate(qapp):
    from node_runner.dialogs.preferences import (
        PreferencesDialog, DEFAULT_TYPE_COLORS, DEFAULT_HIGHLIGHT_COLOR,
        DEFAULT_SIZES)
    initial = {
        'colors': dict(DEFAULT_TYPE_COLORS),
        'highlight_color': DEFAULT_HIGHLIGHT_COLOR,
        'sizes': dict(DEFAULT_SIZES),
    }
    dlg = PreferencesDialog(initial)
    # Programmatically change the RBE2 swatch and a size spinbox.
    dlg._swatches['RBE2'].set_color('#abcdef')
    dlg._size_spins['node_size'].setValue(9)
    payload = dlg.result_payload()
    assert payload['colors']['RBE2'] == '#abcdef'
    assert payload['sizes']['node_size'] == 9


def test_restore_defaults_clears_modifications(qapp):
    from node_runner.dialogs.preferences import (
        PreferencesDialog, DEFAULT_TYPE_COLORS, DEFAULT_HIGHLIGHT_COLOR,
        DEFAULT_SIZES)
    initial = {
        'colors': {**DEFAULT_TYPE_COLORS, 'RBE2': '#000000'},
        'highlight_color': '#000000',
        'sizes': {**DEFAULT_SIZES, 'node_size': 12},
    }
    dlg = PreferencesDialog(initial)
    dlg._restore_all_defaults()
    payload = dlg.result_payload()
    assert payload['colors']['RBE2'] == DEFAULT_TYPE_COLORS['RBE2']
    assert payload['highlight_color'] == DEFAULT_HIGHLIGHT_COLOR
    assert payload['sizes']['node_size'] == DEFAULT_SIZES['node_size']
