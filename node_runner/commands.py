"""Undo/Redo command system using the Command pattern.

Every reversible model mutation is wrapped in a Command subclass. The
CommandManager maintains undo/redo stacks with a configurable history limit.
"""

from __future__ import annotations

import copy
from abc import ABC, abstractmethod
from typing import Any


# ---------------------------------------------------------------------------
# Base classes
# ---------------------------------------------------------------------------

class Command(ABC):
    """Abstract base for all undoable commands."""

    @abstractmethod
    def execute(self, model) -> None:
        """Perform the operation on the pyNastran BDF *model*."""

    @abstractmethod
    def undo(self, model) -> None:
        """Reverse the operation on *model*."""

    @property
    @abstractmethod
    def description(self) -> str:
        """Human-readable label for the status bar (e.g. "Add 3 nodes")."""


class CompoundCommand(Command):
    """Groups multiple commands into a single undo step."""

    def __init__(self, commands: list[Command], desc: str | None = None):
        self._commands = commands
        self._desc = desc

    def execute(self, model) -> None:
        for cmd in self._commands:
            cmd.execute(model)

    def undo(self, model) -> None:
        for cmd in reversed(self._commands):
            cmd.undo(model)

    @property
    def description(self) -> str:
        if self._desc:
            return self._desc
        if self._commands:
            return self._commands[0].description
        return "Compound operation"


class CommandManager:
    """Manages undo/redo stacks with a fixed history depth."""

    def __init__(self, max_history: int = 20):
        self.max_history = max_history
        self.undo_stack: list[Command] = []
        self.redo_stack: list[Command] = []

    def execute(self, cmd: Command, model) -> None:
        """Execute *cmd* and push it onto the undo stack."""
        cmd.execute(model)
        self.undo_stack.append(cmd)
        if len(self.undo_stack) > self.max_history:
            self.undo_stack.pop(0)
        self.redo_stack.clear()

    def undo(self, model) -> bool:
        """Undo the most recent command. Returns True on success."""
        if not self.undo_stack:
            return False
        cmd = self.undo_stack.pop()
        cmd.undo(model)
        self.redo_stack.append(cmd)
        return True

    def redo(self, model) -> bool:
        """Redo the most recently undone command. Returns True on success."""
        if not self.redo_stack:
            return False
        cmd = self.redo_stack.pop()
        cmd.execute(model)
        self.undo_stack.append(cmd)
        return True

    @property
    def can_undo(self) -> bool:
        return bool(self.undo_stack)

    @property
    def can_redo(self) -> bool:
        return bool(self.redo_stack)

    @property
    def undo_description(self) -> str:
        if self.undo_stack:
            return self.undo_stack[-1].description
        return ""

    @property
    def redo_description(self) -> str:
        if self.redo_stack:
            return self.redo_stack[-1].description
        return ""

    def clear(self) -> None:
        """Clear both stacks (called on model replacement)."""
        self.undo_stack.clear()
        self.redo_stack.clear()


# ---------------------------------------------------------------------------
# Node commands
# ---------------------------------------------------------------------------

class AddNodesCommand(Command):
    """Add one or more nodes by (nid, x, y, z)."""

    def __init__(self, nodes: list[tuple[int, float, float, float]]):
        self._nodes = nodes  # [(nid, x, y, z), ...]
        self._created_nids: list[int] = []

    def execute(self, model) -> None:
        self._created_nids.clear()
        for nid, x, y, z in self._nodes:
            actual_nid = model.add_grid(nid, [x, y, z]).nid
            self._created_nids.append(actual_nid)

    def undo(self, model) -> None:
        for nid in self._created_nids:
            model.nodes.pop(nid, None)

    @property
    def description(self) -> str:
        n = len(self._nodes)
        return f"Add {n} node{'s' if n != 1 else ''}"


class DeleteNodesCommand(Command):
    """Delete nodes and cascade-delete connected elements.

    Captures deep copies of all affected entities before deletion so that
    undo can fully restore them.
    """

    def __init__(self, node_ids: list[int]):
        self._node_ids = node_ids
        # Snapshots populated during execute
        self._deleted_nodes: dict[int, Any] = {}
        self._deleted_elements: dict[int, Any] = {}
        self._deleted_rigid_elements: dict[int, Any] = {}
        self._deleted_masses: dict[int, Any] = {}
        self._deleted_plotels: dict[int, Any] = {}

    def execute(self, model) -> None:
        nid_set = set(self._node_ids)
        # Snapshot nodes
        self._deleted_nodes = {nid: copy.deepcopy(model.nodes[nid])
                               for nid in self._node_ids if nid in model.nodes}

        # Find and snapshot connected elements (standard + rigid + masses)
        self._deleted_elements.clear()
        self._deleted_rigid_elements.clear()
        self._deleted_masses.clear()
        for eid, elem in list(model.elements.items()):
            if nid_set.intersection(elem.nodes):
                self._deleted_elements[eid] = copy.deepcopy(elem)
        for eid, elem in list(model.rigid_elements.items()):
            # Rigid elements (RBE2/RBE3) use independent_nodes/dependent_nodes, not .nodes
            try:
                all_nids = set(elem.independent_nodes + elem.dependent_nodes)
            except (AttributeError, TypeError):
                all_nids = set(getattr(elem, 'nodes', []))
            if nid_set.intersection(all_nids):
                self._deleted_rigid_elements[eid] = copy.deepcopy(elem)
        for eid, mass_elem in list(model.masses.items()):
            if mass_elem.nid in nid_set:
                self._deleted_masses[eid] = copy.deepcopy(mass_elem)
        for eid, plotel in list(model.plotels.items()):
            if nid_set.intersection(plotel.nodes):
                self._deleted_plotels[eid] = copy.deepcopy(plotel)

        # Delete elements then nodes
        for eid in self._deleted_elements:
            model.elements.pop(eid, None)
        for eid in self._deleted_rigid_elements:
            model.rigid_elements.pop(eid, None)
        for eid in self._deleted_masses:
            model.masses.pop(eid, None)
        for eid in self._deleted_plotels:
            model.plotels.pop(eid, None)
        for nid in self._deleted_nodes:
            model.nodes.pop(nid, None)

    def undo(self, model) -> None:
        # Restore nodes first, then elements
        for nid, node in self._deleted_nodes.items():
            model.nodes[nid] = node
        for eid, elem in self._deleted_elements.items():
            model.elements[eid] = elem
        for eid, elem in self._deleted_rigid_elements.items():
            model.rigid_elements[eid] = elem
        for eid, mass_elem in self._deleted_masses.items():
            model.masses[eid] = mass_elem
        for eid, plotel in self._deleted_plotels.items():
            model.plotels[eid] = plotel

    @property
    def description(self) -> str:
        n = len(self._node_ids)
        return f"Delete {n} node{'s' if n != 1 else ''}"


class TransformNodesCommand(Command):
    """Apply a spatial transform to a set of nodes.

    Always stores original positions so that any transform type (including
    non-invertible ones like "move to point") can be undone exactly.
    """

    def __init__(self, node_ids: list[int], params: dict):
        self._node_ids = node_ids
        self._params = params
        self._old_positions: dict[int, Any] = {}  # nid -> numpy xyz

    def execute(self, model) -> None:
        import numpy as np
        # Snapshot current positions before transforming
        self._old_positions = {
            nid: model.nodes[nid].get_position().copy()
            for nid in self._node_ids if nid in model.nodes
        }
        # Apply transform using same logic as model.transform_nodes
        _apply_transform(model, self._node_ids, self._params)

    def undo(self, model) -> None:
        for nid, pos in self._old_positions.items():
            if nid in model.nodes:
                model.nodes[nid].set_position(model, pos)

    @property
    def description(self) -> str:
        t = self._params.get('type', 'transform')
        n = len(self._node_ids)
        return f"{t.capitalize()} {n} node{'s' if n != 1 else ''}"


class MergeNodesCommand(Command):
    """Merge coincident node groups.

    This is the most expensive undo command — it must snapshot the full state
    of all affected nodes, elements, SPCs, and loads before merging.
    """

    def __init__(self, groups_to_merge: list[list[int]]):
        self._groups = groups_to_merge
        # Pre-merge snapshots
        self._deleted_nodes: dict[int, Any] = {}
        self._original_elem_nodes: dict[int, list[int]] = {}
        self._original_rigid_elem_nodes: dict[int, list[int]] = {}
        self._original_spc_nodes: list[tuple[int, int, list[int]]] = []
        self._original_load_node_ids: list[tuple[int, int, int | None]] = []

    def execute(self, model) -> None:
        # Build remap dict: slave -> master
        remap: dict[int, int] = {}
        nodes_to_delete: list[int] = []
        for group in self._groups:
            master = min(group)
            for nid in group:
                if nid != master:
                    remap[nid] = master
                    nodes_to_delete.append(nid)

        # Snapshot nodes to delete
        self._deleted_nodes = {nid: copy.deepcopy(model.nodes[nid])
                               for nid in nodes_to_delete if nid in model.nodes}

        # Snapshot and remap standard elements
        self._original_elem_nodes.clear()
        for eid, elem in model.elements.items():
            old_nodes = list(elem.nodes)
            new_nodes = [remap.get(n, n) for n in elem.nodes]
            if old_nodes != new_nodes:
                self._original_elem_nodes[eid] = old_nodes
                elem.nodes = new_nodes

        # Snapshot and remap rigid elements
        self._original_rigid_elem_nodes.clear()
        for eid, elem in model.rigid_elements.items():
            old_nodes = list(elem.nodes)
            new_nodes = [remap.get(n, n) for n in elem.nodes]
            if old_nodes != new_nodes:
                self._original_rigid_elem_nodes[eid] = old_nodes
                elem.nodes = new_nodes

        # Snapshot and remap SPCs
        self._original_spc_nodes.clear()
        for sid, spc_list in model.spcs.items():
            for idx, spc in enumerate(spc_list):
                if hasattr(spc, 'nodes'):
                    old_nodes = list(spc.nodes)
                    new_nodes = list(set(remap.get(n, n) for n in spc.nodes))
                    if old_nodes != new_nodes:
                        self._original_spc_nodes.append((sid, idx, old_nodes))
                        spc.nodes = new_nodes

        # Snapshot and remap loads
        self._original_load_node_ids.clear()
        for sid, load_list in model.loads.items():
            for idx, load in enumerate(load_list):
                if hasattr(load, 'node_id') and load.node_id in remap:
                    self._original_load_node_ids.append((sid, idx, load.node_id))
                    load.node_id = remap[load.node_id]
                elif hasattr(load, 'node') and load.node in remap:
                    self._original_load_node_ids.append((sid, idx, load.node))
                    load.node = remap[load.node]

        # Delete slave nodes
        for nid in nodes_to_delete:
            model.nodes.pop(nid, None)

    def undo(self, model) -> None:
        # Restore deleted nodes
        for nid, node in self._deleted_nodes.items():
            model.nodes[nid] = node

        # Restore original element connectivity
        for eid, old_nodes in self._original_elem_nodes.items():
            if eid in model.elements:
                model.elements[eid].nodes = old_nodes

        for eid, old_nodes in self._original_rigid_elem_nodes.items():
            if eid in model.rigid_elements:
                model.rigid_elements[eid].nodes = old_nodes

        # Restore original SPC nodes
        for sid, idx, old_nodes in self._original_spc_nodes:
            if sid in model.spcs and idx < len(model.spcs[sid]):
                model.spcs[sid][idx].nodes = old_nodes

        # Restore original load node IDs
        for sid, idx, old_node_id in self._original_load_node_ids:
            if sid in model.loads and idx < len(model.loads[sid]):
                load = model.loads[sid][idx]
                if hasattr(load, 'node_id'):
                    load.node_id = old_node_id
                elif hasattr(load, 'node'):
                    load.node = old_node_id

    @property
    def description(self) -> str:
        n = sum(len(g) - 1 for g in self._groups)
        return f"Merge {n} node{'s' if n != 1 else ''}"


# ---------------------------------------------------------------------------
# Element commands
# ---------------------------------------------------------------------------

class AddElementCommand(Command):
    """Add a single element (line, plate, solid, bush, shear, gap, or RBE)."""

    def __init__(self, elem_type: str, params: dict):
        self._elem_type = elem_type  # 'line', 'plate', 'solid', 'bush', 'shear', 'gap', 'rbe'
        self._params = params
        self._created_eid: int | None = None
        self._is_rigid = False

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model

        if self._elem_type == 'line':
            self._created_eid = gen.add_line_element(**self._params)
        elif self._elem_type == 'plate':
            self._created_eid = gen.add_plate_element(**self._params)
        elif self._elem_type == 'solid':
            self._created_eid = gen.add_solid_element(**self._params)
        elif self._elem_type == 'bush':
            self._created_eid = gen.add_cbush(**self._params)
        elif self._elem_type == 'shear':
            self._created_eid = gen.add_cshear(**self._params)
        elif self._elem_type == 'gap':
            self._created_eid = gen.add_cgap(**self._params)
        elif self._elem_type == 'rbe':
            self._created_eid = gen.add_rbe_element(**self._params)
            self._is_rigid = True

    def undo(self, model) -> None:
        if self._created_eid is not None:
            if self._is_rigid:
                model.rigid_elements.pop(self._created_eid, None)
            else:
                model.elements.pop(self._created_eid, None)

    @property
    def description(self) -> str:
        etype = self._params.get('type', self._params.get('elem_type', self._elem_type))
        return f"Add {etype}"


class DeleteElementsCommand(Command):
    """Delete elements by EID, snapshotting them for undo."""

    def __init__(self, eids: list[int]):
        self._eids = eids
        self._deleted_elements: dict[int, Any] = {}
        self._deleted_rigid_elements: dict[int, Any] = {}

    def execute(self, model) -> None:
        self._deleted_elements.clear()
        self._deleted_rigid_elements.clear()
        for eid in self._eids:
            if eid in model.elements:
                self._deleted_elements[eid] = copy.deepcopy(model.elements.pop(eid))
            elif eid in model.rigid_elements:
                self._deleted_rigid_elements[eid] = copy.deepcopy(model.rigid_elements.pop(eid))

    def undo(self, model) -> None:
        for eid, elem in self._deleted_elements.items():
            model.elements[eid] = elem
        for eid, elem in self._deleted_rigid_elements.items():
            model.rigid_elements[eid] = elem

    @property
    def description(self) -> str:
        n = len(self._eids)
        return f"Delete {n} element{'s' if n != 1 else ''}"


class FlipNormalsCommand(Command):
    """Flip shell element normals. Self-inverse operation."""

    def __init__(self, eids: list[int]):
        self._eids = eids

    def execute(self, model) -> None:
        _flip_normals(model, self._eids)

    def undo(self, model) -> None:
        # Self-inverse: flipping again restores original
        _flip_normals(model, self._eids)

    @property
    def description(self) -> str:
        n = len(self._eids)
        return f"Flip normals on {n} element{'s' if n != 1 else ''}"


class AutoOrientNormalsCommand(Command):
    """Auto-orient shell normals for consistency. Self-inverse."""

    def __init__(self, eids_to_flip: list[int]):
        self._eids = eids_to_flip

    def execute(self, model) -> None:
        _flip_normals(model, self._eids)

    def undo(self, model) -> None:
        _flip_normals(model, self._eids)

    @property
    def description(self) -> str:
        n = len(self._eids)
        return f"Auto-orient normals ({n} element{'s' if n != 1 else ''} flipped)"


class EditElementCommand(Command):
    """Edit an element's PID, nodes, and/or orientation."""

    def __init__(self, eid: int, new_pid: int, new_nodes: list[int],
                 new_orientation: dict | None = None):
        self._eid = eid
        self._new_pid = new_pid
        self._new_nodes = new_nodes
        self._new_orientation = new_orientation
        # Old state captured during execute
        self._old_element: Any = None

    def execute(self, model) -> None:
        elem = model.elements[self._eid]
        self._old_element = copy.deepcopy(elem)

        elem.pid = self._new_pid
        elem.nodes = list(self._new_nodes)

        if self._new_orientation is not None:
            elem.g0 = elem.x = None
            if hasattr(elem, 'cid'):
                elem.cid = None
            method = self._new_orientation.get('method')
            if method == 'vector':
                elem.x = self._new_orientation['values']
            elif method == 'node':
                elem.g0 = self._new_orientation['values'][0]
            elif method == 'cid' and hasattr(elem, 'cid'):
                elem.cid = self._new_orientation['values'][0]

    def undo(self, model) -> None:
        if self._old_element is not None:
            model.elements[self._eid] = self._old_element

    @property
    def description(self) -> str:
        return f"Edit element {self._eid}"


# ---------------------------------------------------------------------------
# Material commands
# ---------------------------------------------------------------------------

class AddMaterialCommand(Command):
    """Add a material card from a params dict."""

    def __init__(self, params: dict):
        self._params = params

    def execute(self, model) -> None:
        _add_material(model, self._params)

    def undo(self, model) -> None:
        mid = self._params['mid']
        model.materials.pop(mid, None)

    @property
    def description(self) -> str:
        return f"Add {self._params.get('type', 'material')} {self._params['mid']}"


class EditMaterialCommand(Command):
    """Replace an existing material with new parameters."""

    def __init__(self, mid: int, new_params: dict):
        self._mid = mid
        self._new_params = new_params
        self._old_material: Any = None

    def execute(self, model) -> None:
        if self._mid in model.materials:
            self._old_material = copy.deepcopy(model.materials[self._mid])
        model.materials.pop(self._mid, None)
        _add_material(model, self._new_params)

    def undo(self, model) -> None:
        new_mid = self._new_params['mid']
        model.materials.pop(new_mid, None)
        if self._old_material is not None:
            model.materials[self._mid] = self._old_material

    @property
    def description(self) -> str:
        return f"Edit material {self._mid}"


class DeleteMaterialCommand(Command):
    """Delete a material by MID, snapshotting for undo."""

    def __init__(self, mid: int):
        self._mid = mid
        self._deleted_material: Any = None

    def execute(self, model) -> None:
        if self._mid in model.materials:
            self._deleted_material = copy.deepcopy(model.materials[self._mid])
            del model.materials[self._mid]

    def undo(self, model) -> None:
        if self._deleted_material is not None:
            model.materials[self._mid] = self._deleted_material

    @property
    def description(self) -> str:
        return f"Delete material {self._mid}"


# ---------------------------------------------------------------------------
# Property commands
# ---------------------------------------------------------------------------

class AddPropertyCommand(Command):
    """Add a property card from a params dict."""

    def __init__(self, params: dict):
        self._params = params

    def execute(self, model) -> None:
        _add_property(model, self._params)

    def undo(self, model) -> None:
        pid = self._params['pid']
        model.properties.pop(pid, None)

    @property
    def description(self) -> str:
        return f"Add {self._params.get('type', 'property')} {self._params['pid']}"


class EditPropertyCommand(Command):
    """Replace an existing property with new parameters."""

    def __init__(self, pid: int, new_params: dict):
        self._pid = pid
        self._new_params = new_params
        self._old_property: Any = None

    def execute(self, model) -> None:
        if self._pid in model.properties:
            self._old_property = copy.deepcopy(model.properties[self._pid])
        model.properties.pop(self._pid, None)
        _add_property(model, self._new_params)

    def undo(self, model) -> None:
        new_pid = self._new_params['pid']
        model.properties.pop(new_pid, None)
        if self._old_property is not None:
            model.properties[self._pid] = self._old_property

    @property
    def description(self) -> str:
        return f"Edit property {self._pid}"


class DeletePropertyCommand(Command):
    """Delete a property by PID, snapshotting for undo."""

    def __init__(self, pid: int):
        self._pid = pid
        self._deleted_property: Any = None

    def execute(self, model) -> None:
        if self._pid in model.properties:
            self._deleted_property = copy.deepcopy(model.properties[self._pid])
            del model.properties[self._pid]

    def undo(self, model) -> None:
        if self._deleted_property is not None:
            model.properties[self._pid] = self._deleted_property

    @property
    def description(self) -> str:
        return f"Delete property {self._pid}"


# ---------------------------------------------------------------------------
# Load commands
# ---------------------------------------------------------------------------

class AddLoadCommand(Command):
    """Add load cards (FORCE, MOMENT, PLOAD4, TEMPD, GRAV) for a given SID."""

    def __init__(self, load_type: str, params: dict):
        self._load_type = load_type  # 'nodal', 'pressure', 'temperature', 'gravity'
        self._params = params
        self._sid: int = params['sid']
        self._pre_existing_count: int = 0
        self._pre_existing_tempd = None  # For TEMPD undo (stored in model.tempds)

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model

        if self._load_type == 'temperature':
            # TEMPD stored in model.tempds (single object per SID), not model.loads
            self._pre_existing_tempd = copy.deepcopy(model.tempds.get(self._sid))
            gen.add_default_temperature(**self._params)
        else:
            self._pre_existing_count = len(model.loads.get(self._sid, []))
            if self._load_type == 'nodal':
                gen.add_nodal_load(**self._params)
            elif self._load_type == 'pressure':
                gen.add_pressure_load(**self._params)
            elif self._load_type == 'gravity':
                gen.add_gravity_load(self._params['sid'], self._params['scale'],
                                     self._params['N'], self._params.get('cid', 0))

    def undo(self, model) -> None:
        if self._load_type == 'temperature':
            # Restore or remove TEMPD from model.tempds
            if self._pre_existing_tempd is not None:
                model.tempds[self._sid] = self._pre_existing_tempd
            else:
                model.tempds.pop(self._sid, None)
        else:
            # Remove only the cards we added (slice from pre_existing_count onwards)
            if self._sid in model.loads:
                model.loads[self._sid] = model.loads[self._sid][:self._pre_existing_count]
                if not model.loads[self._sid]:
                    del model.loads[self._sid]

    @property
    def description(self) -> str:
        return f"Add {self._load_type} load SID {self._sid}"


class DeleteLoadCommand(Command):
    """Delete an entire load set by SID (from model.loads and model.tempds)."""

    def __init__(self, sid: int):
        self._sid = sid
        self._deleted_loads: list[Any] = []
        self._deleted_tempd = None

    def execute(self, model) -> None:
        self._deleted_loads = copy.deepcopy(model.loads.get(self._sid, []))
        self._deleted_tempd = copy.deepcopy(model.tempds.get(self._sid))
        model.loads.pop(self._sid, None)
        model.tempds.pop(self._sid, None)

    def undo(self, model) -> None:
        if self._deleted_loads:
            model.loads[self._sid] = self._deleted_loads
        if self._deleted_tempd is not None:
            model.tempds[self._sid] = self._deleted_tempd

    @property
    def description(self) -> str:
        return f"Delete load set {self._sid}"


class DeleteTempdCommand(Command):
    """Delete only the TEMPD entry for a given SID (from model.tempds)."""

    def __init__(self, sid: int):
        self._sid = sid
        self._deleted_tempd = None

    def execute(self, model) -> None:
        self._deleted_tempd = copy.deepcopy(model.tempds.get(self._sid))
        model.tempds.pop(self._sid, None)

    def undo(self, model) -> None:
        if self._deleted_tempd is not None:
            model.tempds[self._sid] = self._deleted_tempd

    @property
    def description(self) -> str:
        return f"Delete TEMPD SID {self._sid}"


# ---------------------------------------------------------------------------
# Constraint commands
# ---------------------------------------------------------------------------

class AddConstraintCommand(Command):
    """Add constraint (SPC1) cards."""

    def __init__(self, params: dict):
        self._params = params
        self._sid: int = params['sid']
        self._pre_existing_count: int = 0

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._pre_existing_count = len(model.spcs.get(self._sid, []))
        gen.add_nodal_constraint(**self._params)

    def undo(self, model) -> None:
        if self._sid in model.spcs:
            model.spcs[self._sid] = model.spcs[self._sid][:self._pre_existing_count]
            if not model.spcs[self._sid]:
                del model.spcs[self._sid]

    @property
    def description(self) -> str:
        return f"Add constraint SID {self._sid}"


class DeleteConstraintCommand(Command):
    """Delete an entire constraint set by SID."""

    def __init__(self, sid: int):
        self._sid = sid
        self._deleted_spcs: list[Any] = []
        self._deleted_spcadds: list[Any] = []

    def execute(self, model) -> None:
        self._deleted_spcs = copy.deepcopy(model.spcs.get(self._sid, []))
        self._deleted_spcadds = copy.deepcopy(model.spcadds.get(self._sid, []))
        model.spcs.pop(self._sid, None)
        model.spcadds.pop(self._sid, None)

    def undo(self, model) -> None:
        if self._deleted_spcs:
            model.spcs[self._sid] = self._deleted_spcs
        if self._deleted_spcadds:
            model.spcadds[self._sid] = self._deleted_spcadds

    @property
    def description(self) -> str:
        return f"Delete constraint set {self._sid}"


# ---------------------------------------------------------------------------
# Coordinate system commands
# ---------------------------------------------------------------------------

class AddCoordCommand(Command):
    """Add a coordinate system (CORD2R, CORD2C, or CORD2S)."""

    # Global/absolute coordinate systems — must never be overwritten.
    _PROTECTED_CIDS = frozenset({0, 1, 2})

    def __init__(self, params: dict):
        self._params = params
        self._cid: int = params['cid']
        self._old_coord: Any = None  # non-None when editing (overwriting)

    def execute(self, model) -> None:
        if self._cid in self._PROTECTED_CIDS:
            return  # refuse to overwrite global coordinate systems
        from node_runner.model import NastranModelGenerator
        # If this CID already exists, snapshot it (edit case)
        if self._cid in model.coords:
            self._old_coord = copy.deepcopy(model.coords[self._cid])
            del model.coords[self._cid]

        gen = NastranModelGenerator()
        gen.model = model

        method = self._params.get('method')
        coord_type = self._params.get('type', 'rectangular')
        cid = self._cid

        comment = self._params.get('comment', '')

        if method == '3 Points':
            origin = self._params['origin']
            z_axis = self._params['z_axis_point']
            xz_plane = self._params['xz_plane_point']
            if coord_type == 'rectangular':
                gen.add_coordinate_system(cid, origin, z_axis, xz_plane, comment=comment)
            elif coord_type == 'cylindrical':
                gen.add_cylindrical_coord_system(cid, origin, z_axis, xz_plane, comment=comment)
            elif coord_type == 'spherical':
                gen.add_spherical_coord_system(cid, origin, z_axis, xz_plane, comment=comment)
        elif method == 'Rotate and Translate':
            gen.add_coord_by_translate_rotate(
                cid, self._params['ref_cid'],
                self._params['translations'], self._params['rotations'],
                coord_type, comment=comment
            )

    def undo(self, model) -> None:
        model.coords.pop(self._cid, None)
        if self._old_coord is not None:
            model.coords[self._cid] = self._old_coord

    @property
    def description(self) -> str:
        action = "Edit" if self._old_coord is not None else "Add"
        return f"{action} coord system {self._cid}"


class DeleteCoordCommand(Command):
    """Delete a coordinate system by CID, snapshotting for undo."""

    # Global/absolute coordinate systems — must never be deleted or modified.
    _PROTECTED_CIDS = frozenset({0, 1, 2})

    def __init__(self, cid: int):
        self._cid = cid
        self._deleted_coord: Any = None

    def execute(self, model) -> None:
        if self._cid in self._PROTECTED_CIDS:
            return  # refuse to delete global coordinate systems
        if self._cid in model.coords:
            self._deleted_coord = copy.deepcopy(model.coords[self._cid])
            del model.coords[self._cid]

    def undo(self, model) -> None:
        if self._deleted_coord is not None:
            model.coords[self._cid] = self._deleted_coord

    @property
    def description(self) -> str:
        return f"Delete coord system {self._cid}"


# ---------------------------------------------------------------------------
# Mass element commands
# ---------------------------------------------------------------------------

class AddMassCommand(Command):
    """Add a CONM2 concentrated mass element."""

    def __init__(self, params: dict):
        self._params = params
        self._created_eid: int | None = None

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._created_eid = gen.add_conm2(**self._params)

    def undo(self, model) -> None:
        if self._created_eid is not None:
            model.masses.pop(self._created_eid, None)

    @property
    def description(self) -> str:
        return f"Add CONM2 at node {self._params.get('nid', '?')}"


class DeleteMassCommand(Command):
    """Delete mass elements by EID, snapshotting them for undo."""

    def __init__(self, eids):
        self._eids = list(eids)
        self._deleted: dict = {}

    def execute(self, model) -> None:
        for eid in self._eids:
            if eid in model.masses:
                self._deleted[eid] = copy.deepcopy(model.masses[eid])
                del model.masses[eid]

    def undo(self, model) -> None:
        for eid, mass_elem in self._deleted.items():
            model.masses[eid] = mass_elem

    @property
    def description(self) -> str:
        return f"Delete {len(self._eids)} mass element(s)"


class AddPlotelCommand(Command):
    """Add a PLOTEL visualization element."""

    def __init__(self, params: dict):
        self._params = params
        self._created_eid: int | None = None

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._created_eid = gen.add_plotel(**self._params)

    def undo(self, model) -> None:
        if self._created_eid is not None:
            model.plotels.pop(self._created_eid, None)

    @property
    def description(self) -> str:
        return f"Add PLOTEL {self._params.get('nodes', '?')}"


class DeletePlotelCommand(Command):
    """Delete PLOTEL elements by EID, snapshotting them for undo."""

    def __init__(self, eids):
        self._eids = list(eids)
        self._deleted: dict = {}

    def execute(self, model) -> None:
        for eid in self._eids:
            if eid in model.plotels:
                self._deleted[eid] = copy.deepcopy(model.plotels[eid])
                del model.plotels[eid]

    def undo(self, model) -> None:
        for eid, elem in self._deleted.items():
            model.plotels[eid] = elem

    @property
    def description(self) -> str:
        return f"Delete {len(self._eids)} PLOTEL(s)"


# ---------------------------------------------------------------------------
# Private helpers (shared logic extracted from MainWindow / model)
# ---------------------------------------------------------------------------

def _add_material(model, params: dict) -> None:
    """Add a material card to *model* from a params dict.

    Mirrors MainWindow._add_material_card_from_params but operates on model
    directly so commands can call it without needing MainWindow.
    """
    p = params.copy()
    mat_type = p.pop('type')
    mid = p.pop('mid')
    comment = p.pop('comment', '')

    if mat_type == 'MAT1':
        model.add_mat1(mid, p['E'], p['G'], p['nu'],
                       rho=p.get('rho', 0.0), a=p.get('a', 0.0),
                       tref=p.get('tref', 0.0), ge=p.get('ge', 0.0),
                       comment=comment)
    elif mat_type == 'MAT8':
        model.add_mat8(mid, p['E1'], p['E2'], p['nu12'],
                       g12=p.get('G12', 0.0), g1z=p.get('G1z', 0.0),
                       g2z=p.get('G2z', 0.0), rho=p.get('rho', 0.0),
                       a1=p.get('a1', 0.0), a2=p.get('a2', 0.0),
                       tref=p.get('tref', 0.0), comment=comment)
    elif mat_type == 'MAT9':
        g_values = {k.lower(): v for k, v in p.items() if k.startswith('G')}
        model.add_mat9(mid, rho=p.get('rho', 0.0), tref=p.get('tref', 0.0),
                       comment=comment, **g_values)


def _add_property(model, params: dict) -> None:
    """Add a property card to *model* from a params dict.

    Mirrors MainWindow._add_property_card_from_params.
    """
    p = params.copy()
    prop_type = p.pop('type')
    pid = p.pop('pid')
    comment = p.pop('comment', '')

    if prop_type == 'PSHELL':
        model.add_pshell(pid, mid1=p['mid1'], t=p['t'],
                         nsm=p['nsm'], comment=comment)
    elif prop_type == 'PCOMP':
        plies = p['plies']
        if plies:
            mids, thicknesses, thetas, souts = zip(*plies)
            model.add_pcomp(pid, list(mids), list(thicknesses), list(thetas),
                            souts=list(souts), nsm=p['nsm'], ft=p['ft'],
                            comment=comment)
        else:
            model.add_pcomp(pid, [], [], [], nsm=p['nsm'], ft=p['ft'],
                            comment=comment)
    elif prop_type == 'PBAR':
        model.add_pbar(pid, p['mid'], p['A'], p['i1'], p['i2'], p['j'],
                       comment=comment)
    elif prop_type == 'PBEAM':
        model.add_pbeam(pid, p['mid'], [0.0], ['C'], [p['A']], [p['i1']],
                        [p['i2']], [0.0], [p['j']], comment=comment)
    elif prop_type == 'PBUSH':
        model.add_pbush(pid, k=p['k'], b=[], ge=[], comment=comment)
    elif prop_type == 'PSOLID':
        model.add_psolid(pid, p['mid'], comment=comment)
    elif prop_type == 'PROD':
        model.add_prod(pid, p['mid'], p['A'], j=p.get('j', 0.0),
                       c=p.get('c', 0.0), nsm=p.get('nsm', 0.0), comment=comment)
    elif prop_type == 'PSHEAR':
        model.add_pshear(pid, p['mid'], p['t'], nsm=p.get('nsm', 0.0),
                         comment=comment)
    elif prop_type == 'PGAP':
        model.add_pgap(pid, u0=p.get('u0', 0.0), f0=p.get('f0', 0.0),
                       ka=p.get('ka', 1.0e8), kb=p.get('kb', None),
                       comment=comment)


def _flip_normals(model, eids: list[int]) -> int:
    """Flip shell element normals in-place. Returns count of flipped elements."""
    count = 0
    for eid in eids:
        elem = model.elements.get(eid)
        if elem is None:
            continue
        if elem.type in ('CQUAD4', 'CMEMBRAN') and len(elem.nodes) == 4:
            n = elem.nodes
            elem.nodes = [n[0], n[3], n[2], n[1]]
            count += 1
        elif elem.type == 'CTRIA3' and len(elem.nodes) == 3:
            n = elem.nodes
            elem.nodes = [n[0], n[2], n[1]]
            count += 1
    return count


def _apply_transform(model, node_ids: list[int], params: dict) -> None:
    """Apply a node transform using the same logic as
    NastranModelGenerator.transform_nodes.
    """
    import numpy as np
    import math

    transform_type = params['type']
    nodes = [model.nodes[nid] for nid in node_ids if nid in model.nodes]
    if not nodes:
        return

    if transform_type == 'translate':
        delta = np.array(params['delta'])
        for node in nodes:
            node.xyz += delta

    elif transform_type == 'move_to_point':
        target = np.array(params['target'])
        positions = np.array([n.get_position() for n in nodes])
        centroid = positions.mean(axis=0)
        delta = target - centroid
        for node in nodes:
            node.xyz += delta

    elif transform_type == 'scale':
        factors = np.array(params['factors'])
        center_type = params.get('center_type', 'origin')
        if center_type == 'centroid':
            positions = np.array([n.get_position() for n in nodes])
            center = positions.mean(axis=0)
        elif center_type == 'custom':
            center = np.array(params['custom_center'])
        else:
            center = np.array([0.0, 0.0, 0.0])
        for node in nodes:
            pos = node.get_position()
            new_pos = center + (pos - center) * factors
            node.set_position(model, new_pos)

    elif transform_type == 'rotate':
        angle_deg = params['angle']
        axis = params['axis']
        angle_rad = math.radians(angle_deg)
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        if axis == 'x':
            R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        elif axis == 'y':
            R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        else:  # z
            R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

        center_type = params.get('center_type', 'origin')
        if center_type == 'centroid':
            positions = np.array([n.get_position() for n in nodes])
            center = positions.mean(axis=0)
        elif center_type == 'custom':
            center = np.array(params['custom_center'])
        else:
            center = np.array([0.0, 0.0, 0.0])

        for node in nodes:
            pos = node.get_position()
            new_pos = center + R @ (pos - center)
            node.set_position(model, new_pos)


# ---------------------------------------------------------------------------
# Find/Replace commands
# ---------------------------------------------------------------------------


class ReassignPropertyCommand(Command):
    """Reassign PID on a set of elements."""

    def __init__(self, eids: list[int], new_pid: int):
        self._eids = eids
        self._new_pid = new_pid
        self._old_pids: dict[int, int] = {}

    def execute(self, model) -> None:
        for eid in self._eids:
            elem = model.elements.get(eid)
            if elem and hasattr(elem, 'pid'):
                self._old_pids[eid] = elem.pid
                elem.pid = self._new_pid

    def undo(self, model) -> None:
        for eid, old_pid in self._old_pids.items():
            if eid in model.elements:
                model.elements[eid].pid = old_pid

    @property
    def description(self) -> str:
        return f"Reassign {len(self._eids)} elements to PID {self._new_pid}"


class ReassignMaterialCommand(Command):
    """Reassign MID on a set of properties."""

    def __init__(self, pids: list[int], new_mid: int):
        self._pids = pids
        self._new_mid = new_mid
        self._old_mids: dict[int, int] = {}

    def execute(self, model) -> None:
        for pid in self._pids:
            prop = model.properties.get(pid)
            if prop:
                if hasattr(prop, 'mid'):
                    self._old_mids[pid] = prop.mid
                    prop.mid = self._new_mid
                elif hasattr(prop, 'mid1'):
                    self._old_mids[pid] = prop.mid1
                    prop.mid1 = self._new_mid

    def undo(self, model) -> None:
        for pid, old_mid in self._old_mids.items():
            prop = model.properties.get(pid)
            if prop:
                if hasattr(prop, 'mid'):
                    prop.mid = old_mid
                elif hasattr(prop, 'mid1'):
                    prop.mid1 = old_mid

    @property
    def description(self) -> str:
        return f"Reassign {len(self._pids)} properties to MID {self._new_mid}"


# ---------------------------------------------------------------------------
# LOAD combination command
# ---------------------------------------------------------------------------


class AddLoadCombinationCommand(Command):
    """Add a LOAD combination card."""

    def __init__(self, params: dict):
        self._params = params
        self._sid: int = params['sid']

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        gen.add_load_combination(
            self._sid, self._params['overall_scale'], self._params['components']
        )

    def undo(self, model) -> None:
        if hasattr(model, 'load_combinations') and self._sid in model.load_combinations:
            del model.load_combinations[self._sid]
        # Also check model.loads in case pyNastran stores it there
        if self._sid in model.loads:
            del model.loads[self._sid]

    @property
    def description(self) -> str:
        return f"Add LOAD combination SID {self._sid}"


# ---------------------------------------------------------------------------
# Renumber command
# ---------------------------------------------------------------------------


class RenumberCommand(Command):
    """Renumber all node and/or element IDs sequentially."""

    def __init__(self, start_nid: int, start_eid: int,
                 do_nodes: bool, do_elements: bool):
        self._start_nid = start_nid
        self._start_eid = start_eid
        self._do_nodes = do_nodes
        self._do_elements = do_elements
        self._old_model_snapshot: Any = None

    def execute(self, model) -> None:
        self._old_model_snapshot = copy.deepcopy(model)
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        gen.renumber_ids(self._start_nid, self._start_eid,
                         self._do_nodes, self._do_elements)

    def undo(self, model) -> None:
        if self._old_model_snapshot is None:
            return
        # Restore all core dicts from the snapshot
        snapshot = self._old_model_snapshot
        model.nodes = snapshot.nodes
        model.elements = snapshot.elements
        model.rigid_elements = snapshot.rigid_elements
        model.properties = snapshot.properties
        model.materials = snapshot.materials
        model.loads = snapshot.loads
        model.spcs = snapshot.spcs
        model.masses = snapshot.masses
        if hasattr(snapshot, 'plotels'):
            model.plotels = snapshot.plotels
        model.coords = snapshot.coords

    @property
    def description(self) -> str:
        return "Renumber model IDs"


# ---------------------------------------------------------------------------
# Group commands  (operate on application-level groups dict, not BDF model)
# ---------------------------------------------------------------------------

class CreateGroupCommand(Command):
    """Create a named group with an optional numeric group ID."""

    def __init__(self, groups_dict: dict, name: str, gid: int = 0):
        self._groups = groups_dict
        self._name = name
        self._gid = gid

    def execute(self, model) -> None:
        self._groups[self._name] = {"nodes": [], "elements": [], "gid": self._gid}

    def undo(self, model) -> None:
        self._groups.pop(self._name, None)

    @property
    def description(self) -> str:
        return f"Create group {self._gid}: '{self._name}'"


class DeleteGroupCommand(Command):
    """Delete a named group, snapshotting for undo."""

    def __init__(self, groups_dict: dict, name: str):
        self._groups = groups_dict
        self._name = name
        self._snapshot: dict | None = None

    def execute(self, model) -> None:
        self._snapshot = self._groups.pop(self._name, None)

    def undo(self, model) -> None:
        if self._snapshot is not None:
            self._groups[self._name] = self._snapshot

    @property
    def description(self) -> str:
        return f"Delete group '{self._name}'"


class RenameGroupCommand(Command):
    """Rename a group."""

    def __init__(self, groups_dict: dict, old_name: str, new_name: str):
        self._groups = groups_dict
        self._old = old_name
        self._new = new_name

    def execute(self, model) -> None:
        if self._old in self._groups:
            self._groups[self._new] = self._groups.pop(self._old)

    def undo(self, model) -> None:
        if self._new in self._groups:
            self._groups[self._old] = self._groups.pop(self._new)

    @property
    def description(self) -> str:
        return f"Rename group '{self._old}' to '{self._new}'"


class ModifyGroupCommand(Command):
    """Add or remove nodes/elements from a group."""

    def __init__(self, groups_dict: dict, name: str,
                 add_nodes: list | None = None,
                 add_elements: list | None = None,
                 remove_nodes: list | None = None,
                 remove_elements: list | None = None):
        self._groups = groups_dict
        self._name = name
        self._add_nodes = add_nodes or []
        self._add_elements = add_elements or []
        self._remove_nodes = remove_nodes or []
        self._remove_elements = remove_elements or []

    def execute(self, model) -> None:
        g = self._groups.get(self._name)
        if not g:
            return
        g["nodes"] = list(
            (set(g["nodes"]) | set(self._add_nodes)) - set(self._remove_nodes)
        )
        g["elements"] = list(
            (set(g["elements"]) | set(self._add_elements)) - set(self._remove_elements)
        )

    def undo(self, model) -> None:
        g = self._groups.get(self._name)
        if not g:
            return
        g["nodes"] = list(
            (set(g["nodes"]) | set(self._remove_nodes)) - set(self._add_nodes)
        )
        g["elements"] = list(
            (set(g["elements"]) | set(self._remove_elements)) - set(self._add_elements)
        )

    @property
    def description(self) -> str:
        return f"Modify group '{self._name}'"


# ---------------------------------------------------------------------------
# Geometry commands
# ---------------------------------------------------------------------------

class AddGeometryPointCommand(Command):
    """Add a geometry point."""

    def __init__(self, geometry_store, point_id: int, xyz: list):
        self._store = geometry_store
        self._point_id = point_id
        self._xyz = xyz

    def execute(self, model) -> None:
        self._store.add_point(self._xyz, point_id=self._point_id)

    def undo(self, model) -> None:
        self._store.points.pop(self._point_id, None)

    @property
    def description(self) -> str:
        return f"Add geometry point P{self._point_id}"


class AddGeometryLineCommand(Command):
    """Add a geometry line between two points."""

    def __init__(self, geometry_store, line_id: int,
                 start_id: int, end_id: int):
        self._store = geometry_store
        self._line_id = line_id
        self._start_id = start_id
        self._end_id = end_id

    def execute(self, model) -> None:
        self._store.add_line(self._start_id, self._end_id,
                             line_id=self._line_id)

    def undo(self, model) -> None:
        self._store.lines.pop(self._line_id, None)

    @property
    def description(self) -> str:
        return f"Add geometry line {self._line_id}"


class AddGeometryArcCommand(Command):
    """Add a geometry arc through three points."""

    def __init__(self, geometry_store, arc_id: int,
                 start_id: int, mid_id: int, end_id: int):
        self._store = geometry_store
        self._arc_id = arc_id
        self._start_id = start_id
        self._mid_id = mid_id
        self._end_id = end_id

    def execute(self, model) -> None:
        self._store.add_arc(self._start_id, self._mid_id, self._end_id,
                            arc_id=self._arc_id)

    def undo(self, model) -> None:
        self._store.arcs.pop(self._arc_id, None)

    @property
    def description(self) -> str:
        return f"Add geometry arc {self._arc_id}"


class AddGeometryCircleCommand(Command):
    """Add a geometry circle."""

    def __init__(self, geometry_store, circle_id: int,
                 center_id: int, radius: float, normal: list):
        self._store = geometry_store
        self._circle_id = circle_id
        self._center_id = center_id
        self._radius = radius
        self._normal = normal

    def execute(self, model) -> None:
        self._store.add_circle(self._center_id, self._radius, self._normal,
                               circle_id=self._circle_id)

    def undo(self, model) -> None:
        self._store.circles.pop(self._circle_id, None)

    @property
    def description(self) -> str:
        return f"Add geometry circle {self._circle_id}"


class AddGeometrySurfaceCommand(Command):
    """Add a geometry surface defined by boundary curves."""

    def __init__(self, geometry_store, surface_id: int, curve_ids: list):
        self._store = geometry_store
        self._surface_id = surface_id
        self._curve_ids = curve_ids

    def execute(self, model) -> None:
        self._store.add_surface(self._curve_ids, surface_id=self._surface_id)

    def undo(self, model) -> None:
        self._store.surfaces.pop(self._surface_id, None)

    @property
    def description(self) -> str:
        return f"Add geometry surface {self._surface_id}"


class MeshCurveCommand(Command):
    """Mesh along geometry curves creating nodes and line elements."""

    def __init__(self, geometry_store, curve_ids: list,
                 n_elements: int, elem_type: str, pid: int):
        self._store = geometry_store
        self._curve_ids = curve_ids
        self._n_elements = n_elements
        self._elem_type = elem_type
        self._pid = pid
        self._created_nids: list[int] = []
        self._created_eids: list[int] = []

    def execute(self, model) -> None:
        import numpy as np
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._created_nids.clear()
        self._created_eids.clear()

        for curve_id in self._curve_ids:
            pts = self._store.evaluate_curve(curve_id, self._n_elements + 1)
            # Create nodes
            nids = []
            for pt in pts:
                nid = gen.add_node(None, float(pt[0]), float(pt[1]), float(pt[2]))
                nids.append(nid)
                self._created_nids.append(nid)
            # Create elements between consecutive nodes
            for i in range(len(nids) - 1):
                eid = gen.add_line_element(
                    None, self._pid, nids[i], nids[i + 1],
                    self._elem_type, None)
                self._created_eids.append(eid)

    def undo(self, model) -> None:
        for eid in self._created_eids:
            model.elements.pop(eid, None)
        for nid in self._created_nids:
            model.nodes.pop(nid, None)

    @property
    def description(self) -> str:
        n = len(self._curve_ids)
        return f"Mesh {n} curve{'s' if n != 1 else ''}"


class NodesAtGeometryPointsCommand(Command):
    """Create FEM nodes at geometry point positions."""

    def __init__(self, geometry_store, point_ids: list):
        self._store = geometry_store
        self._point_ids = point_ids
        self._created_nids: list[int] = []

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._created_nids.clear()
        for pid in self._point_ids:
            pt = self._store.points[pid]
            nid = gen.add_node(None, float(pt.xyz[0]), float(pt.xyz[1]),
                               float(pt.xyz[2]))
            self._created_nids.append(nid)

    def undo(self, model) -> None:
        for nid in self._created_nids:
            model.nodes.pop(nid, None)

    @property
    def description(self) -> str:
        n = len(self._point_ids)
        return f"Create {n} node{'s' if n != 1 else ''} at geometry points"


# ---------------------------------------------------------------------------
# Geometry deletion commands
# ---------------------------------------------------------------------------

class DeleteGeometryPointsCommand(Command):
    """Delete geometry points with cascade deletion of dependent curves and surfaces."""

    def __init__(self, geometry_store, point_ids: list[int]):
        self._store = geometry_store
        self._point_ids = point_ids
        self._deleted_points: dict = {}
        self._deleted_lines: dict = {}
        self._deleted_arcs: dict = {}
        self._deleted_circles: dict = {}
        self._deleted_surfaces: dict = {}

    def execute(self, model) -> None:
        # Find dependent curves
        cascade_curve_ids = set()
        for pid in self._point_ids:
            cascade_curve_ids.update(self._store.curves_referencing_point(pid))

        # Find dependent surfaces
        cascade_surface_ids = set(self._store.surfaces_referencing_curves(cascade_curve_ids))

        # Snapshot everything before deletion
        self._deleted_points = {pid: copy.deepcopy(self._store.points[pid])
                                for pid in self._point_ids if pid in self._store.points}
        self._deleted_lines = {cid: copy.deepcopy(self._store.lines[cid])
                               for cid in cascade_curve_ids if cid in self._store.lines}
        self._deleted_arcs = {cid: copy.deepcopy(self._store.arcs[cid])
                              for cid in cascade_curve_ids if cid in self._store.arcs}
        self._deleted_circles = {cid: copy.deepcopy(self._store.circles[cid])
                                 for cid in cascade_curve_ids if cid in self._store.circles}
        self._deleted_surfaces = {sid: copy.deepcopy(self._store.surfaces[sid])
                                  for sid in cascade_surface_ids if sid in self._store.surfaces}

        # Delete in order: surfaces, curves, points
        for sid in self._deleted_surfaces:
            self._store.surfaces.pop(sid, None)
        for cid in self._deleted_lines:
            self._store.lines.pop(cid, None)
        for cid in self._deleted_arcs:
            self._store.arcs.pop(cid, None)
        for cid in self._deleted_circles:
            self._store.circles.pop(cid, None)
        for pid in self._deleted_points:
            self._store.points.pop(pid, None)

    def undo(self, model) -> None:
        # Restore in reverse: points, curves, surfaces
        for pid, pt in self._deleted_points.items():
            self._store.points[pid] = pt
        for cid, line in self._deleted_lines.items():
            self._store.lines[cid] = line
        for cid, arc in self._deleted_arcs.items():
            self._store.arcs[cid] = arc
        for cid, circle in self._deleted_circles.items():
            self._store.circles[cid] = circle
        for sid, surf in self._deleted_surfaces.items():
            self._store.surfaces[sid] = surf

    @property
    def description(self) -> str:
        n = len(self._point_ids)
        return f"Delete {n} geometry point{'s' if n != 1 else ''}"


class DeleteGeometryCurvesCommand(Command):
    """Delete geometry curves with cascade deletion of dependent surfaces."""

    def __init__(self, geometry_store, curve_ids: list[int]):
        self._store = geometry_store
        self._curve_ids = curve_ids
        self._deleted_lines: dict = {}
        self._deleted_arcs: dict = {}
        self._deleted_circles: dict = {}
        self._deleted_surfaces: dict = {}

    def execute(self, model) -> None:
        # Find dependent surfaces
        cascade_surface_ids = set(self._store.surfaces_referencing_curves(set(self._curve_ids)))

        # Snapshot
        for cid in self._curve_ids:
            if cid in self._store.lines:
                self._deleted_lines[cid] = copy.deepcopy(self._store.lines[cid])
            elif cid in self._store.arcs:
                self._deleted_arcs[cid] = copy.deepcopy(self._store.arcs[cid])
            elif cid in self._store.circles:
                self._deleted_circles[cid] = copy.deepcopy(self._store.circles[cid])
        self._deleted_surfaces = {sid: copy.deepcopy(self._store.surfaces[sid])
                                  for sid in cascade_surface_ids if sid in self._store.surfaces}

        # Delete
        for sid in self._deleted_surfaces:
            self._store.surfaces.pop(sid, None)
        for cid in self._deleted_lines:
            self._store.lines.pop(cid, None)
        for cid in self._deleted_arcs:
            self._store.arcs.pop(cid, None)
        for cid in self._deleted_circles:
            self._store.circles.pop(cid, None)

    def undo(self, model) -> None:
        for cid, line in self._deleted_lines.items():
            self._store.lines[cid] = line
        for cid, arc in self._deleted_arcs.items():
            self._store.arcs[cid] = arc
        for cid, circle in self._deleted_circles.items():
            self._store.circles[cid] = circle
        for sid, surf in self._deleted_surfaces.items():
            self._store.surfaces[sid] = surf

    @property
    def description(self) -> str:
        n = len(self._curve_ids)
        return f"Delete {n} geometry curve{'s' if n != 1 else ''}"


class DeleteGeometrySurfacesCommand(Command):
    """Delete geometry surfaces."""

    def __init__(self, geometry_store, surface_ids: list[int]):
        self._store = geometry_store
        self._surface_ids = surface_ids
        self._deleted_surfaces: dict = {}

    def execute(self, model) -> None:
        self._deleted_surfaces = {sid: copy.deepcopy(self._store.surfaces[sid])
                                  for sid in self._surface_ids if sid in self._store.surfaces}
        for sid in self._deleted_surfaces:
            self._store.surfaces.pop(sid, None)

    def undo(self, model) -> None:
        for sid, surf in self._deleted_surfaces.items():
            self._store.surfaces[sid] = surf

    @property
    def description(self) -> str:
        n = len(self._surface_ids)
        return f"Delete {n} geometry surface{'s' if n != 1 else ''}"


class MeshSurfaceCommand(Command):
    """Mesh a geometry surface using gmsh, creating CQUAD4/CTRIA3 elements."""

    def __init__(self, geometry_store, surface_id: int,
                 n_per_edge: int, pid: int, elem_preference: str = 'quad'):
        self._store = geometry_store
        self._surface_id = surface_id
        self._n_per_edge = n_per_edge
        self._pid = pid
        self._elem_preference = elem_preference
        self._created_nids: list[int] = []
        self._created_eids: list[int] = []

    def execute(self, model) -> None:
        import numpy as np
        try:
            import gmsh
        except ImportError:
            raise RuntimeError("gmsh is not installed. Run: pip install gmsh")

        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model
        self._created_nids.clear()
        self._created_eids.clear()

        surface = self._store.surfaces[self._surface_id]
        curve_ids = surface.boundary_curve_ids

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("surface_mesh")

        try:
            # Build boundary from geometry curves
            gmsh_points = {}
            gmsh_lines = []
            tag_counter = 1

            all_boundary_pts = []
            for cid in curve_ids:
                pts = self._store.evaluate_curve(cid, self._n_per_edge + 1)
                all_boundary_pts.append(pts)

            # Estimate target mesh size from boundary
            total_length = 0
            for pts in all_boundary_pts:
                for i in range(len(pts) - 1):
                    total_length += np.linalg.norm(pts[i + 1] - pts[i])
            n_total_edges = self._n_per_edge * len(curve_ids)
            target_size = total_length / max(n_total_edges, 1)

            # Add boundary points and lines to gmsh
            ordered_point_tags = []
            for pts in all_boundary_pts:
                for pt in pts[:-1]:  # skip last (duplicate of next curve start)
                    key = tuple(np.round(pt, 10))
                    if key not in gmsh_points:
                        gmsh.model.occ.addPoint(
                            float(pt[0]), float(pt[1]), float(pt[2]),
                            meshSize=target_size, tag=tag_counter)
                        gmsh_points[key] = tag_counter
                        tag_counter += 1
                    ordered_point_tags.append(gmsh_points[key])

            # Create line segments between consecutive boundary points
            line_tags = []
            n_pts = len(ordered_point_tags)
            for i in range(n_pts):
                p1 = ordered_point_tags[i]
                p2 = ordered_point_tags[(i + 1) % n_pts]
                if p1 != p2:
                    lt = gmsh.model.occ.addLine(p1, p2)
                    line_tags.append(lt)

            # Create curve loop and plane surface
            cl = gmsh.model.occ.addCurveLoop(line_tags)
            ps = gmsh.model.occ.addPlaneSurface([cl])
            gmsh.model.occ.synchronize()

            # Mesh sizing
            gmsh.option.setNumber("Mesh.MeshSizeMin", target_size * 0.5)
            gmsh.option.setNumber("Mesh.MeshSizeMax", target_size * 1.5)
            if self._elem_preference == 'quad':
                gmsh.option.setNumber("Mesh.RecombineAll", 1)
                gmsh.option.setNumber("Mesh.Algorithm", 8)

            gmsh.model.mesh.generate(2)

            # Extract results
            node_tags, coords, _ = gmsh.model.mesh.getNodes()
            coords = coords.reshape(-1, 3)

            # Map gmsh node tags to FEM node IDs
            gmsh_to_nid = {}
            for i, gtag in enumerate(node_tags):
                nid = gen.add_node(
                    None, float(coords[i, 0]), float(coords[i, 1]),
                    float(coords[i, 2]))
                gmsh_to_nid[int(gtag)] = nid
                self._created_nids.append(nid)

            # Extract elements
            elem_types, elem_tags_list, elem_node_tags_list = \
                gmsh.model.mesh.getElements(dim=2)

            for et, e_tags, en_tags in zip(
                    elem_types, elem_tags_list, elem_node_tags_list):
                n_per_elem = len(en_tags) // len(e_tags)
                for j in range(len(e_tags)):
                    gnodes = [int(en_tags[j * n_per_elem + k])
                              for k in range(n_per_elem)]
                    fem_nodes = [gmsh_to_nid[g] for g in gnodes]
                    if n_per_elem == 4:
                        eid = gen.add_plate_element(
                            None, self._pid, fem_nodes, 'CQUAD4')
                    elif n_per_elem == 3:
                        eid = gen.add_plate_element(
                            None, self._pid, fem_nodes, 'CTRIA3')
                    else:
                        continue
                    self._created_eids.append(eid)

        finally:
            gmsh.finalize()

    def undo(self, model) -> None:
        for eid in self._created_eids:
            model.elements.pop(eid, None)
        for nid in self._created_nids:
            model.nodes.pop(nid, None)

    @property
    def description(self) -> str:
        return f"Mesh surface {self._surface_id}"


class ImportCADCommand(Command):
    """Import a STEP/IGES/STL file, creating GRID + shell elements."""

    def __init__(self, filepath: str, mesh_size: float, elem_preference: str,
                 pid: int):
        self._filepath = filepath
        self._mesh_size = mesh_size
        self._elem_preference = elem_preference
        self._pid = pid
        self._created_nids: list[int] = []
        self._created_eids: list[int] = []

    def execute(self, model) -> None:
        import os
        from node_runner.model import NastranModelGenerator

        gen = NastranModelGenerator()
        gen.model = model

        self._created_nids.clear()
        self._created_eids.clear()

        ext = os.path.splitext(self._filepath)[1].lower()
        if ext == '.stl':
            nids, eids = gen.import_stl_file(self._filepath, self._pid)
        elif ext in ('.step', '.stp', '.iges', '.igs'):
            nids, eids = gen.import_cad_file(
                self._filepath, self._mesh_size,
                self._elem_preference, self._pid)
        else:
            raise ValueError(f"Unsupported file type: {ext}")

        self._created_nids = list(nids)
        self._created_eids = list(eids)

    def undo(self, model) -> None:
        for eid in self._created_eids:
            model.elements.pop(eid, None)
        for nid in self._created_nids:
            model.nodes.pop(nid, None)

    @property
    def description(self) -> str:
        import os
        return f"Import CAD: {os.path.basename(self._filepath)}"


class AddSpiderCommand(Command):
    """Create a spider connection: center node + RBE2/RBE3 to ring nodes."""

    def __init__(self, ring_nids: list[int], rbe_type: str, dof: str):
        self._ring_nids = ring_nids
        self._rbe_type = rbe_type  # 'RBE2' or 'RBE3'
        self._dof = dof
        self._center_nid: int | None = None
        self._rbe_eid: int | None = None

    def execute(self, model) -> None:
        import numpy as np
        from node_runner.model import NastranModelGenerator

        gen = NastranModelGenerator()
        gen.model = model

        # Compute centroid of ring nodes
        coords = [model.nodes[nid].get_position() for nid in self._ring_nids
                  if nid in model.nodes]
        if not coords:
            return
        centroid = np.mean(coords, axis=0)

        # Create center node
        self._center_nid = gen.add_node(
            None, float(centroid[0]), float(centroid[1]), float(centroid[2]))

        # Create RBE element
        eid = gen.get_next_available_id('element')
        if self._rbe_type == 'RBE2':
            model.add_rbe2(eid, self._center_nid, self._dof, self._ring_nids)
        elif self._rbe_type == 'RBE3':
            model.add_rbe3(eid, self._center_nid, self._dof, self._ring_nids)
        self._rbe_eid = eid

    def undo(self, model) -> None:
        if self._rbe_eid is not None:
            model.rigid_elements.pop(self._rbe_eid, None)
        if self._center_nid is not None:
            model.nodes.pop(self._center_nid, None)

    @property
    def description(self) -> str:
        return (f"Spider {self._rbe_type}: center + "
                f"{len(self._ring_nids)} ring nodes")


class AddWeldCommand(Command):
    """Create a CWELD fastener connection between two surfaces."""

    def __init__(self, nid_a: int, nid_b: int, diameter: float, pid: int):
        self._nid_a = nid_a
        self._nid_b = nid_b
        self._diameter = diameter
        self._pid = pid
        self._created_pid: int | None = None
        self._created_eid: int | None = None

    def execute(self, model) -> None:
        from node_runner.model import NastranModelGenerator
        gen = NastranModelGenerator()
        gen.model = model

        # Create PWELD property if it doesn't exist
        if self._pid not in model.properties:
            mid = min(model.materials.keys()) if model.materials else 1
            if mid not in model.materials:
                model.add_mat1(mid, 70000.0, 0.33, rho=2.7e-9,
                               comment='Auto-created for CWELD')
            # Use PBUSH as a lightweight connector property (pyNastran compatible)
            model.add_pbush(self._pid,
                            [1e6, 1e6, 1e6, 1e4, 1e4, 1e4],
                            comment=f'CWELD D={self._diameter}')
            self._created_pid = self._pid

        # Create CBUSH element as the weld connector
        eid = gen.get_next_available_id('element')
        model.add_cbush(eid, self._pid, [self._nid_a, self._nid_b],
                        comment=f'CWELD D={self._diameter}')
        self._created_eid = eid

    def undo(self, model) -> None:
        if self._created_eid is not None:
            model.elements.pop(self._created_eid, None)
        if self._created_pid is not None:
            model.properties.pop(self._created_pid, None)

    @property
    def description(self) -> str:
        return f"Weld: N{self._nid_a}-N{self._nid_b} D={self._diameter}"


# ---------------------------------------------------------------------------
# Geometry modify commands
# ---------------------------------------------------------------------------

class EditGeometryPointCommand(Command):
    """Move a geometry point to a new position, updating referencing curves."""

    def __init__(self, geometry_store, point_id: int, new_xyz):
        self._store = geometry_store
        self._point_id = point_id
        self._new_xyz = new_xyz
        self._old_xyz = None
        self._old_curves: dict = {}

    def execute(self, model) -> None:
        self._old_xyz = self._store.points[self._point_id].xyz.copy()
        for cid in self._store.curves_referencing_point(self._point_id):
            curve = self._store.get_curve(cid)
            self._old_curves[cid] = copy.deepcopy(curve)
        self._store.update_point(self._point_id, self._new_xyz)

    def undo(self, model) -> None:
        self._store.points[self._point_id].xyz = self._old_xyz
        for cid, curve in self._old_curves.items():
            if cid in self._store.lines:
                self._store.lines[cid] = curve
            elif cid in self._store.arcs:
                self._store.arcs[cid] = curve
            elif cid in self._store.circles:
                self._store.circles[cid] = curve

    @property
    def description(self) -> str:
        return f"Edit geometry point P{self._point_id}"


class EditGeometryCurveCommand(Command):
    """Edit a curve's point references or parameters."""

    def __init__(self, geometry_store, curve_id: int, params: dict):
        self._store = geometry_store
        self._curve_id = curve_id
        self._params = params
        self._old_curve = None

    def execute(self, model) -> None:
        curve = self._store.get_curve(self._curve_id)
        self._old_curve = copy.deepcopy(curve)

        if self._curve_id in self._store.lines:
            self._store.update_line_refs(
                self._curve_id,
                self._params['start_point_id'],
                self._params['end_point_id'])
        elif self._curve_id in self._store.arcs:
            self._store.update_arc_refs(
                self._curve_id,
                self._params['start_point_id'],
                self._params['mid_point_id'],
                self._params['end_point_id'])
        elif self._curve_id in self._store.circles:
            self._store.update_circle(
                self._curve_id,
                center_id=self._params.get('center_point_id'),
                radius=self._params.get('radius'),
                normal=self._params.get('normal'))

    def undo(self, model) -> None:
        if self._old_curve is not None:
            if self._curve_id in self._store.lines:
                self._store.lines[self._curve_id] = self._old_curve
            elif self._curve_id in self._store.arcs:
                self._store.arcs[self._curve_id] = self._old_curve
            elif self._curve_id in self._store.circles:
                self._store.circles[self._curve_id] = self._old_curve

    @property
    def description(self) -> str:
        return f"Edit geometry curve {self._curve_id}"


class EditGeometrySurfaceBoundariesCommand(Command):
    """Change the boundary curve list of a geometry surface."""

    def __init__(self, geometry_store, surface_id: int, new_curve_ids: list):
        self._store = geometry_store
        self._surface_id = surface_id
        self._new_curve_ids = new_curve_ids
        self._old_curve_ids = None

    def execute(self, model) -> None:
        self._old_curve_ids = list(
            self._store.surfaces[self._surface_id].boundary_curve_ids)
        self._store.update_surface_boundaries(self._surface_id,
                                              self._new_curve_ids)

    def undo(self, model) -> None:
        self._store.update_surface_boundaries(self._surface_id,
                                              self._old_curve_ids)

    @property
    def description(self) -> str:
        return f"Edit surface {self._surface_id} boundaries"


# ---------------------------------------------------------------------------
# Geometry transform commands
# ---------------------------------------------------------------------------

def _apply_geometry_transform(store, point_ids, params):
    """Apply a spatial transform to geometry points."""
    import numpy as np
    import math as _math

    transform_type = params['type']
    points = [store.points[pid] for pid in point_ids if pid in store.points]
    if not points:
        return

    if transform_type == 'translate':
        delta = np.array(params['delta'], dtype=float)
        for pt in points:
            pt.xyz = pt.xyz + delta

    elif transform_type == 'rotate':
        angle_rad = _math.radians(params['angle'])
        axis = params['axis']
        c, s = _math.cos(angle_rad), _math.sin(angle_rad)
        if axis == 'x':
            R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        elif axis == 'y':
            R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        else:
            R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        center = np.array(params.get('center', [0, 0, 0]), dtype=float)
        for pt in points:
            pt.xyz = center + R @ (pt.xyz - center)

    elif transform_type == 'mirror':
        plane_point = np.array(params['plane_point'], dtype=float)
        plane_normal = np.array(params['plane_normal'], dtype=float)
        plane_normal = plane_normal / np.linalg.norm(plane_normal)
        for pt in points:
            d = np.dot(pt.xyz - plane_point, plane_normal)
            pt.xyz = pt.xyz - 2 * d * plane_normal

    elif transform_type == 'scale':
        factor = params['factor']
        center = np.array(params.get('center', [0, 0, 0]), dtype=float)
        for pt in points:
            pt.xyz = center + (pt.xyz - center) * factor


class TransformGeometryPointsCommand(Command):
    """Apply translate/rotate/mirror/scale to geometry points, updating curves."""

    def __init__(self, geometry_store, point_ids: list, params: dict):
        self._store = geometry_store
        self._point_ids = list(point_ids)
        self._params = params
        self._old_positions: dict = {}
        self._old_curves: dict = {}

    def execute(self, model) -> None:
        import numpy as np
        self._old_positions = {
            pid: self._store.points[pid].xyz.copy()
            for pid in self._point_ids if pid in self._store.points
        }
        affected_cids = set()
        for pid in self._point_ids:
            affected_cids.update(self._store.curves_referencing_point(pid))
        for cid in affected_cids:
            self._old_curves[cid] = copy.deepcopy(self._store.get_curve(cid))

        _apply_geometry_transform(self._store, self._point_ids, self._params)

        refreshed = set()
        for pid in self._point_ids:
            for cid in self._store.curves_referencing_point(pid):
                if cid not in refreshed:
                    refreshed.add(cid)
        for pid in self._point_ids:
            self._store._refresh_curves_for_point(pid)

    def undo(self, model) -> None:
        for pid, pos in self._old_positions.items():
            self._store.points[pid].xyz = pos
        for cid, curve in self._old_curves.items():
            if cid in self._store.lines:
                self._store.lines[cid] = curve
            elif cid in self._store.arcs:
                self._store.arcs[cid] = curve
            elif cid in self._store.circles:
                self._store.circles[cid] = curve

    @property
    def description(self) -> str:
        t = self._params.get('type', 'transform')
        n = len(self._point_ids)
        return f"{t.capitalize()} {n} geometry point{'s' if n != 1 else ''}"


class CopyGeometryCommand(Command):
    """Copy geometry points and their dependent curves with a translation offset."""

    def __init__(self, geometry_store, point_ids: list, delta):
        import numpy as np
        self._store = geometry_store
        self._point_ids = list(point_ids)
        self._delta = np.array(delta, dtype=float)
        self._created_point_ids: list = []
        self._created_curve_ids: list = []

    def execute(self, model) -> None:
        import numpy as np
        self._created_point_ids.clear()
        self._created_curve_ids.clear()

        point_map = {}
        for pid in self._point_ids:
            if pid not in self._store.points:
                continue
            new_xyz = self._store.points[pid].xyz + self._delta
            new_pid = self._store.next_point_id()
            self._store.add_point(new_xyz, point_id=new_pid)
            self._created_point_ids.append(new_pid)
            point_map[pid] = new_pid

        selected = set(self._point_ids)
        for lid, line in list(self._store.lines.items()):
            if (line.start_point_id in selected and
                    line.end_point_id in selected):
                cid = self._store.next_curve_id()
                self._store.add_line(
                    point_map[line.start_point_id],
                    point_map[line.end_point_id],
                    line_id=cid)
                self._created_curve_ids.append(cid)

        for aid, arc in list(self._store.arcs.items()):
            if (arc.start_point_id in selected and
                    arc.mid_point_id in selected and
                    arc.end_point_id in selected):
                cid = self._store.next_curve_id()
                self._store.add_arc(
                    point_map[arc.start_point_id],
                    point_map[arc.mid_point_id],
                    point_map[arc.end_point_id],
                    arc_id=cid)
                self._created_curve_ids.append(cid)

        for ccid, circle in list(self._store.circles.items()):
            if circle.center_point_id in selected:
                cid = self._store.next_curve_id()
                self._store.add_circle(
                    point_map[circle.center_point_id],
                    circle.radius,
                    circle.normal.copy(),
                    circle_id=cid)
                self._created_curve_ids.append(cid)

    def undo(self, model) -> None:
        for cid in self._created_curve_ids:
            self._store.lines.pop(cid, None)
            self._store.arcs.pop(cid, None)
            self._store.circles.pop(cid, None)
        for pid in self._created_point_ids:
            self._store.points.pop(pid, None)

    @property
    def description(self) -> str:
        n = len(self._point_ids)
        return f"Copy {n} geometry point{'s' if n != 1 else ''}"


# ---------------------------------------------------------------------------
# CAD-like geometry operations
# ---------------------------------------------------------------------------

class SplitCurveCommand(Command):
    """Split a curve into N segments, creating new points and sub-curves."""

    def __init__(self, geometry_store, curve_id: int, n_segments: int = 2):
        self._store = geometry_store
        self._curve_id = curve_id
        self._n_segments = n_segments
        self._created_point_ids: list = []
        self._created_curve_ids: list = []
        self._old_curve = None
        self._old_curve_type: str = None
        self._old_surface_boundaries: dict = {}

    def execute(self, model) -> None:
        curve = self._store.get_curve(self._curve_id)
        self._old_curve = copy.deepcopy(curve)

        if self._curve_id in self._store.lines:
            self._old_curve_type = 'line'
        elif self._curve_id in self._store.arcs:
            self._old_curve_type = 'arc'
        elif self._curve_id in self._store.circles:
            raise ValueError("Cannot split a full circle. Convert to arc first.")

        ts = [i / self._n_segments for i in range(1, self._n_segments)]
        split_xyzs = [curve.evaluate(t) for t in ts]

        for xyz in split_xyzs:
            pid = self._store.next_point_id()
            self._store.add_point(xyz, point_id=pid)
            self._created_point_ids.append(pid)

        if self._old_curve_type == 'line':
            all_pids = ([curve.start_point_id] +
                        self._created_point_ids +
                        [curve.end_point_id])
            for i in range(len(all_pids) - 1):
                cid = self._store.next_curve_id()
                self._store.add_line(all_pids[i], all_pids[i + 1],
                                     line_id=cid)
                self._created_curve_ids.append(cid)

        elif self._old_curve_type == 'arc':
            segment_ts = [0.0] + ts + [1.0]
            split_point_idx = 0
            for i in range(len(segment_ts) - 1):
                t_start = segment_ts[i]
                t_end = segment_ts[i + 1]
                t_mid = (t_start + t_end) / 2.0

                if i == 0:
                    start_pid = curve.start_point_id
                else:
                    start_pid = self._created_point_ids[split_point_idx - 1]

                if i == len(segment_ts) - 2:
                    end_pid = curve.end_point_id
                else:
                    end_pid = self._created_point_ids[split_point_idx]
                    split_point_idx += 1

                mid_xyz = curve.evaluate(t_mid)
                mid_pid = self._store.next_point_id()
                self._store.add_point(mid_xyz, point_id=mid_pid)
                self._created_point_ids.append(mid_pid)

                cid = self._store.next_curve_id()
                self._store.add_arc(start_pid, mid_pid, end_pid,
                                    arc_id=cid)
                self._created_curve_ids.append(cid)

        for sid in self._store.surfaces_referencing_curve(self._curve_id):
            surf = self._store.surfaces[sid]
            self._old_surface_boundaries[sid] = list(
                surf.boundary_curve_ids)
            idx = surf.boundary_curve_ids.index(self._curve_id)
            surf.boundary_curve_ids[idx:idx + 1] = self._created_curve_ids

        if self._old_curve_type == 'line':
            self._store.lines.pop(self._curve_id, None)
        elif self._old_curve_type == 'arc':
            self._store.arcs.pop(self._curve_id, None)

    def undo(self, model) -> None:
        for cid in self._created_curve_ids:
            self._store.lines.pop(cid, None)
            self._store.arcs.pop(cid, None)
        for pid in self._created_point_ids:
            self._store.points.pop(pid, None)

        if self._old_curve_type == 'line':
            self._store.lines[self._curve_id] = self._old_curve
        elif self._old_curve_type == 'arc':
            self._store.arcs[self._curve_id] = self._old_curve

        for sid, old_bounds in self._old_surface_boundaries.items():
            self._store.surfaces[sid].boundary_curve_ids = old_bounds

    @property
    def description(self) -> str:
        return f"Split curve {self._curve_id} into {self._n_segments} segments"


class OffsetCurveCommand(Command):
    """Create a parallel curve at a given distance."""

    def __init__(self, geometry_store, curve_id: int, distance: float,
                 direction=None):
        self._store = geometry_store
        self._curve_id = curve_id
        self._distance = distance
        self._direction = direction
        self._created_point_ids: list = []
        self._created_curve_id = None

    def execute(self, model) -> None:
        import numpy as np

        self._created_point_ids.clear()
        self._created_curve_id = None

        if self._curve_id in self._store.lines:
            line = self._store.lines[self._curve_id]
            direction = line._end_xyz - line._start_xyz
            if self._direction is not None:
                ref = np.array(self._direction, dtype=float)
            else:
                ref = np.array([0, 0, 1], dtype=float)
            perp = np.cross(direction, ref)
            norm = np.linalg.norm(perp)
            if norm < 1e-12:
                ref = np.array([0, 1, 0], dtype=float)
                perp = np.cross(direction, ref)
                norm = np.linalg.norm(perp)
            perp = perp / norm * self._distance

            p1 = self._store.next_point_id()
            self._store.add_point(line._start_xyz + perp, point_id=p1)
            self._created_point_ids.append(p1)

            p2 = self._store.next_point_id()
            self._store.add_point(line._end_xyz + perp, point_id=p2)
            self._created_point_ids.append(p2)

            cid = self._store.next_curve_id()
            self._store.add_line(p1, p2, line_id=cid)
            self._created_curve_id = cid

        elif self._curve_id in self._store.arcs:
            arc = self._store.arcs[self._curve_id]
            new_radius = arc._radius + self._distance
            if new_radius <= 0:
                raise ValueError("Offset would create negative radius")
            for t in [0.0, 0.5, 1.0]:
                pt_on_arc = arc.evaluate(t)
                radial = pt_on_arc - arc._center
                radial_unit = radial / np.linalg.norm(radial)
                offset_pt = pt_on_arc + self._distance * radial_unit
                pid = self._store.next_point_id()
                self._store.add_point(offset_pt, point_id=pid)
                self._created_point_ids.append(pid)

            cid = self._store.next_curve_id()
            self._store.add_arc(
                self._created_point_ids[0],
                self._created_point_ids[1],
                self._created_point_ids[2],
                arc_id=cid)
            self._created_curve_id = cid

        elif self._curve_id in self._store.circles:
            circle = self._store.circles[self._curve_id]
            new_radius = circle.radius + self._distance
            if new_radius <= 0:
                raise ValueError("Offset would create negative radius")
            cid = self._store.next_curve_id()
            self._store.add_circle(
                circle.center_point_id, new_radius,
                circle.normal.copy(), circle_id=cid)
            self._created_curve_id = cid

    def undo(self, model) -> None:
        if self._created_curve_id is not None:
            self._store.lines.pop(self._created_curve_id, None)
            self._store.arcs.pop(self._created_curve_id, None)
            self._store.circles.pop(self._created_curve_id, None)
        for pid in self._created_point_ids:
            self._store.points.pop(pid, None)

    @property
    def description(self) -> str:
        return f"Offset curve {self._curve_id} by {self._distance}"


class ProjectPointToCurveCommand(Command):
    """Find closest point on a curve and create a geometry point there."""

    def __init__(self, geometry_store, curve_id: int, source_xyz):
        import numpy as np
        self._store = geometry_store
        self._curve_id = curve_id
        self._source_xyz = np.array(source_xyz, dtype=float)
        self._created_point_id = None

    def execute(self, model) -> None:
        import numpy as np
        curve = self._store.get_curve(self._curve_id)
        pts = curve.evaluate_array(200)
        dists = np.linalg.norm(pts - self._source_xyz, axis=1)
        closest_xyz = pts[int(np.argmin(dists))]

        pid = self._store.next_point_id()
        self._store.add_point(closest_xyz, point_id=pid)
        self._created_point_id = pid

    def undo(self, model) -> None:
        if self._created_point_id is not None:
            self._store.points.pop(self._created_point_id, None)

    @property
    def description(self) -> str:
        return f"Project point onto curve {self._curve_id}"


class FilletCommand(Command):
    """Create a fillet arc between two lines sharing an endpoint."""

    def __init__(self, geometry_store, line_id_1: int, line_id_2: int,
                 radius: float):
        self._store = geometry_store
        self._line_id_1 = line_id_1
        self._line_id_2 = line_id_2
        self._radius = radius
        self._created_point_ids: list = []
        self._created_arc_id = None
        self._old_line_1 = None
        self._old_line_2 = None

    def execute(self, model) -> None:
        import numpy as np
        import math as _math

        line1 = self._store.lines[self._line_id_1]
        line2 = self._store.lines[self._line_id_2]
        self._old_line_1 = copy.deepcopy(line1)
        self._old_line_2 = copy.deepcopy(line2)

        l1_pids = {line1.start_point_id, line1.end_point_id}
        l2_pids = {line2.start_point_id, line2.end_point_id}
        common = l1_pids & l2_pids
        if not common:
            raise ValueError("Lines do not share an endpoint")
        shared_pid = common.pop()

        corner = self._store.points[shared_pid].xyz
        other1_pid = (line1.start_point_id if line1.end_point_id == shared_pid
                      else line1.end_point_id)
        other2_pid = (line2.start_point_id if line2.end_point_id == shared_pid
                      else line2.end_point_id)

        dir1 = self._store.points[other1_pid].xyz - corner
        dir1_len = np.linalg.norm(dir1)
        if dir1_len < 1e-12:
            raise ValueError("Line 1 has zero length")
        dir1 = dir1 / dir1_len

        dir2 = self._store.points[other2_pid].xyz - corner
        dir2_len = np.linalg.norm(dir2)
        if dir2_len < 1e-12:
            raise ValueError("Line 2 has zero length")
        dir2 = dir2 / dir2_len

        cos_angle = np.clip(np.dot(dir1, dir2), -1.0, 1.0)
        full_angle = _math.acos(cos_angle)
        half_angle = full_angle / 2.0
        if half_angle < 1e-6 or abs(half_angle - _math.pi / 2) < 1e-6:
            if half_angle < 1e-6:
                raise ValueError("Lines are parallel, cannot fillet")

        tangent_length = self._radius / _math.tan(half_angle)
        if tangent_length > dir1_len or tangent_length > dir2_len:
            raise ValueError(
                "Fillet radius too large for the line lengths")

        tp1 = corner + tangent_length * dir1
        tp2 = corner + tangent_length * dir2

        bisector = dir1 + dir2
        bis_len = np.linalg.norm(bisector)
        if bis_len < 1e-12:
            raise ValueError("Lines are anti-parallel, cannot fillet")
        bisector = bisector / bis_len
        mid_dist = self._radius / _math.sin(half_angle)
        arc_mid_pt = corner + (mid_dist - self._radius) * bisector

        tp1_id = self._store.next_point_id()
        self._store.add_point(tp1, point_id=tp1_id)
        self._created_point_ids.append(tp1_id)

        mid_id = self._store.next_point_id()
        self._store.add_point(arc_mid_pt, point_id=mid_id)
        self._created_point_ids.append(mid_id)

        tp2_id = self._store.next_point_id()
        self._store.add_point(tp2, point_id=tp2_id)
        self._created_point_ids.append(tp2_id)

        arc_cid = self._store.next_curve_id()
        self._store.add_arc(tp1_id, mid_id, tp2_id, arc_id=arc_cid)
        self._created_arc_id = arc_cid

        if line1.end_point_id == shared_pid:
            self._store.update_line_refs(
                self._line_id_1, line1.start_point_id, tp1_id)
        else:
            self._store.update_line_refs(
                self._line_id_1, tp1_id, line1.end_point_id)

        if line2.end_point_id == shared_pid:
            self._store.update_line_refs(
                self._line_id_2, line2.start_point_id, tp2_id)
        else:
            self._store.update_line_refs(
                self._line_id_2, tp2_id, line2.end_point_id)

    def undo(self, model) -> None:
        if self._created_arc_id is not None:
            self._store.arcs.pop(self._created_arc_id, None)
        for pid in self._created_point_ids:
            self._store.points.pop(pid, None)
        if self._old_line_1 is not None:
            self._store.lines[self._line_id_1] = self._old_line_1
        if self._old_line_2 is not None:
            self._store.lines[self._line_id_2] = self._old_line_2

    @property
    def description(self) -> str:
        return (f"Fillet R={self._radius} between curves "
                f"{self._line_id_1}, {self._line_id_2}")


# ---------------------------------------------------------------------------
# Load/Constraint set replacement commands
# ---------------------------------------------------------------------------

class ReplaceLoadSetCommand(Command):
    """Atomically replace an entire load set with a new list of entries."""

    def __init__(self, sid: int, new_entries: list):
        self._sid = sid
        self._new_entries = new_entries
        self._old_entries = None

    def execute(self, model) -> None:
        self._old_entries = copy.deepcopy(model.loads.get(self._sid, []))
        if self._new_entries:
            model.loads[self._sid] = self._new_entries
        else:
            model.loads.pop(self._sid, None)

    def undo(self, model) -> None:
        if self._old_entries:
            model.loads[self._sid] = self._old_entries
        else:
            model.loads.pop(self._sid, None)

    @property
    def description(self) -> str:
        return f"Modify load set {self._sid}"


class ReplaceConstraintSetCommand(Command):
    """Atomically replace an entire constraint set with a new list of entries."""

    def __init__(self, sid: int, new_entries: list):
        self._sid = sid
        self._new_entries = new_entries
        self._old_entries = None

    def execute(self, model) -> None:
        self._old_entries = copy.deepcopy(model.spcs.get(self._sid, []))
        if self._new_entries:
            model.spcs[self._sid] = self._new_entries
        else:
            model.spcs.pop(self._sid, None)

    def undo(self, model) -> None:
        if self._old_entries:
            model.spcs[self._sid] = self._old_entries
        else:
            model.spcs.pop(self._sid, None)

    @property
    def description(self) -> str:
        return f"Modify constraint set {self._sid}"
