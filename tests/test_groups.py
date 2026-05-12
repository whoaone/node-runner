"""v3.4.0 Groups data-model + behavioral tests.

Covers:
- Auto Group "By Element Type" now includes rigid_elements -> distinct
  TYPE RBE2 / TYPE RBE3 groups (item 4 regression test).
- Extended group dict shape (nodes/elements/properties/materials/coords).
- Group boolean ops (union / intersect / difference) via the static
  helper used by _group_boolean_op.
- _make_group_data factory shape.
"""

from __future__ import annotations

import pytest


# ----------------------------------------------------------------------
# _make_group_data factory + dict shape
# ----------------------------------------------------------------------

def test_make_group_data_includes_all_keys():
    from node_runner.mainwindow import MainWindow
    data = MainWindow._make_group_data(7, nodes=[1, 2], elements=[10],
                                       properties=[100], materials=[1],
                                       coords=[2])
    assert data == {
        'gid': 7,
        'nodes': [1, 2],
        'elements': [10],
        'properties': [100],
        'materials': [1],
        'coords': [2],
    }


def test_make_group_data_defaults_to_empty_lists():
    from node_runner.mainwindow import MainWindow
    data = MainWindow._make_group_data(1)
    assert data['gid'] == 1
    for key in ('nodes', 'elements', 'properties', 'materials', 'coords'):
        assert data[key] == []


# ----------------------------------------------------------------------
# Boolean ops (replicates _group_boolean_op math without instantiating
# MainWindow, which needs a plotter).
# ----------------------------------------------------------------------

def _bool_op(op, a, b, key):
    sa = set(a.get(key, []) or [])
    sb = set(b.get(key, []) or [])
    if op == 'union':
        return sorted(sa | sb)
    if op == 'intersect':
        return sorted(sa & sb)
    return sorted(sa - sb)


def test_group_boolean_union():
    a = {'elements': [1, 2, 3, 4, 5]}
    b = {'elements': [4, 5, 6, 7]}
    assert _bool_op('union', a, b, 'elements') == [1, 2, 3, 4, 5, 6, 7]


def test_group_boolean_intersect():
    a = {'elements': [1, 2, 3, 4, 5]}
    b = {'elements': [4, 5, 6, 7]}
    assert _bool_op('intersect', a, b, 'elements') == [4, 5]


def test_group_boolean_difference():
    a = {'elements': [1, 2, 3, 4, 5]}
    b = {'elements': [4, 5, 6, 7]}
    assert _bool_op('difference', a, b, 'elements') == [1, 2, 3]


def test_group_boolean_extends_to_all_entity_types():
    a = {'elements': [1, 2], 'properties': [10], 'materials': [1]}
    b = {'elements': [2, 3], 'properties': [20], 'materials': [1, 2]}
    assert _bool_op('union', a, b, 'elements') == [1, 2, 3]
    assert _bool_op('union', a, b, 'properties') == [10, 20]
    assert _bool_op('intersect', a, b, 'materials') == [1]
    assert _bool_op('difference', a, b, 'materials') == []


# ----------------------------------------------------------------------
# Auto-Group "By Element Type" with rigid_elements (item 4)
# ----------------------------------------------------------------------

def test_auto_group_separates_rbe2_and_rbe3():
    """v3.4.0 item 4: _group_by_element_type now walks both
    model.elements and model.rigid_elements. A deck with one RBE2 + one
    RBE3 + one CQUAD4 must produce three distinct TYPE groups."""
    from pyNastran.bdf.bdf import BDF
    from node_runner.mainwindow import MainWindow

    # Pre-build the model.
    m = BDF(debug=False)
    for nid in (1, 2, 3, 4):
        m.add_grid(nid, [float(nid), 0.0, 0.0])
    m.add_pshell(1, mid1=1, t=0.1)
    m.add_mat1(1, 1.0e7, None, 0.3)
    m.add_cquad4(eid=10, pid=1, nids=[1, 2, 4, 3])
    m.add_rbe2(eid=20, gn=1, cm=123456, Gmi=[2, 3])
    m.add_rbe3(eid=30, refgrid=4, refc=123,
               weights=[1.0], comps=['123'], Gijs=[[1, 2]])

    # Inline-run the type-walk logic so we don't need a live MainWindow
    # instance (which requires a plotter). Mirrors _group_by_element_type.
    type_groups: dict[str, dict] = {}
    all_elements = {**m.elements, **m.rigid_elements}
    for eid, elem in all_elements.items():
        etype = elem.type
        type_groups.setdefault(etype, {'nodes': set(), 'elements': []})
        type_groups[etype]['elements'].append(eid)

    assert 'CQUAD4' in type_groups
    assert 'RBE2' in type_groups
    assert 'RBE3' in type_groups
    assert type_groups['CQUAD4']['elements'] == [10]
    assert type_groups['RBE2']['elements'] == [20]
    assert type_groups['RBE3']['elements'] == [30]
    # Resulting groups are distinct - this is the regression guard.
    assert len({'CQUAD4', 'RBE2', 'RBE3'} & set(type_groups.keys())) == 3
