"""Cross-section / clipping plane controller and plane definitions.

Femap-style "chop the model and look inside" feature, implemented at the
VTK mapper level so it costs next to nothing at render time and works
uniformly across every actor in the scene (skins, beams, nodes, ghost
overlays, selection highlights).

Mechanics:
  * One ``vtkPlane`` (origin + normal) defines a half-space.
  * When enabled, the plane is added to every actor's mapper via
    ``mapper.AddClippingPlane(plane)``. Anything in the negative half-space
    is fragment-discarded by the GL pipeline - no geometry rebuild, no
    extra dataset copies.
  * Live updates to the origin and normal trigger a debounced render.
  * After a scene rebuild (e.g. _update_plot_visibility), call
    ``reapply()`` to push the plane into the new actors' mappers.

The plane itself can be specified several ways via ``PlaneDefinition``:
global axis, 3 picked nodes, normal+point, coordinate value, two nodes
(line normal). The controller stays method-agnostic; the dock and dialog
do the picking and resolve to (origin, normal) before pushing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Tuple

import vtk


# Cardinal axes for the dropdown. Each entry is (label, normal-vector).
AXIS_PRESETS = (
    ("+X", (1.0, 0.0, 0.0)),
    ("-X", (-1.0, 0.0, 0.0)),
    ("+Y", (0.0, 1.0, 0.0)),
    ("-Y", (0.0, -1.0, 0.0)),
    ("+Z", (0.0, 0.0, 1.0)),
    ("-Z", (0.0, 0.0, -1.0)),
)


# ---------------------------------------------------------------------------
# Plane definition value object
# ---------------------------------------------------------------------------

# Methods, in priority order from the picker dialog. The label is what the
# dock summary shows; the key is the persisted/serialized identifier.
METHOD_AXIS = "axis"
METHOD_THREE_POINT = "three_point"
METHOD_NORMAL_POINT = "normal_point"
METHOD_COORD_VALUE = "coord_value"
METHOD_TWO_NODES = "two_nodes"

METHOD_LABELS = {
    METHOD_AXIS: "Global Axis",
    METHOD_THREE_POINT: "3 Points",
    METHOD_NORMAL_POINT: "Normal + Point",
    METHOD_COORD_VALUE: "Coordinate Value",
    METHOD_TWO_NODES: "Two Nodes (Line Normal)",
}


@dataclass
class PlaneDefinition:
    """Value object describing how a cross-section plane is defined.

    The dock keeps one of these and the controller resolves it to
    (origin, normal) before pushing into vtkPlane.

    Fields are method-specific:
      * axis: ``axis_label`` is one of '+X', '-X', '+Y', '-Y', '+Z', '-Z'.
      * three_point: ``node_ids`` is a 3-tuple of grid IDs.
      * normal_point: ``normal`` and ``point`` (each a 3-tuple of floats).
      * coord_value: ``axis_label`` (which axis, '+X' family) and ``value``
        (where along that axis to cut).
      * two_nodes: ``node_ids`` is a 2-tuple; plane is perpendicular to the
        line between them, through their midpoint.

    All variants resolve to (origin, normal) via ``resolve()``.
    """

    method: str = METHOD_AXIS
    axis_label: str = "+X"
    value: float = 0.0
    point: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    normal: Tuple[float, float, float] = (1.0, 0.0, 0.0)
    node_ids: Tuple[int, ...] = field(default_factory=tuple)

    # Cached resolved geometry (filled by resolve()) - the dock reads
    # these so the slider can scrub along the resolved normal.
    _resolved_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    _resolved_normal: Tuple[float, float, float] = (1.0, 0.0, 0.0)

    def summary(self) -> str:
        """Human-readable single-line summary for the dock label."""
        if self.method == METHOD_AXIS:
            return f"Axis {self.axis_label}"
        if self.method == METHOD_COORD_VALUE:
            return f"Axis {self.axis_label} at {self.value:g}"
        if self.method == METHOD_THREE_POINT:
            ids = ", ".join(str(n) for n in self.node_ids[:3])
            return f"3 Points: nodes {ids}"
        if self.method == METHOD_NORMAL_POINT:
            n = self.normal
            p = self.point
            return (f"Normal ({n[0]:g}, {n[1]:g}, {n[2]:g}) "
                    f"at ({p[0]:g}, {p[1]:g}, {p[2]:g})")
        if self.method == METHOD_TWO_NODES:
            ids = ", ".join(str(n) for n in self.node_ids[:2])
            return f"Two Nodes: {ids}"
        return self.method

    def resolve(self, model=None) -> Tuple[Tuple[float, float, float],
                                           Tuple[float, float, float]]:
        """Return (origin, normal) for this definition.

        ``model`` is a pyNastran BDF model; required for any method that
        references node IDs. Methods that don't need it accept None.
        Raises ``PlaneDefinitionError`` if the definition can't be resolved.
        """
        if self.method == METHOD_AXIS:
            origin = (0.0, 0.0, 0.0)
            normal = _axis_label_to_vec(self.axis_label)
        elif self.method == METHOD_COORD_VALUE:
            normal = _axis_label_to_vec(self.axis_label)
            # Origin sits on the chosen axis at value, with the other two
            # components zero. This is the simplest and most intuitive
            # choice: "X = 1500" => plane at x=1500 normal to X.
            v = float(self.value)
            origin = (v if abs(normal[0]) > 0.5 else 0.0,
                      v if abs(normal[1]) > 0.5 else 0.0,
                      v if abs(normal[2]) > 0.5 else 0.0)
        elif self.method == METHOD_NORMAL_POINT:
            normal = _normalize(self.normal)
            origin = tuple(float(c) for c in self.point)
        elif self.method == METHOD_THREE_POINT:
            if len(self.node_ids) != 3 or model is None:
                raise PlaneDefinitionError("3-point method needs three node IDs and a model")
            p1 = _node_xyz(model, self.node_ids[0])
            p2 = _node_xyz(model, self.node_ids[1])
            p3 = _node_xyz(model, self.node_ids[2])
            v1 = _sub(p2, p1)
            v2 = _sub(p3, p1)
            n = _cross(v1, v2)
            if _norm(n) < 1e-12:
                raise PlaneDefinitionError("3 points are colinear; pick three non-colinear nodes")
            normal = _normalize(n)
            origin = p1
        elif self.method == METHOD_TWO_NODES:
            if len(self.node_ids) != 2 or model is None:
                raise PlaneDefinitionError("Two-nodes method needs two node IDs and a model")
            p1 = _node_xyz(model, self.node_ids[0])
            p2 = _node_xyz(model, self.node_ids[1])
            v = _sub(p2, p1)
            if _norm(v) < 1e-12:
                raise PlaneDefinitionError("The two nodes are coincident")
            normal = _normalize(v)
            origin = ((p1[0] + p2[0]) * 0.5,
                      (p1[1] + p2[1]) * 0.5,
                      (p1[2] + p2[2]) * 0.5)
        else:
            raise PlaneDefinitionError(f"Unknown plane method: {self.method}")

        self._resolved_origin = origin
        self._resolved_normal = normal
        return origin, normal


class PlaneDefinitionError(ValueError):
    """Raised when a PlaneDefinition can't be resolved (bad nodes, etc.)."""


# ---------------------------------------------------------------------------
# Small geometry helpers (kept here so the dock doesn't pull numpy)
# ---------------------------------------------------------------------------

def _axis_label_to_vec(label: str) -> Tuple[float, float, float]:
    for lab, vec in AXIS_PRESETS:
        if lab == label:
            return vec
    return (1.0, 0.0, 0.0)


def _node_xyz(model, nid: int) -> Tuple[float, float, float]:
    node = model.nodes.get(int(nid)) if model is not None else None
    if node is None:
        raise PlaneDefinitionError(f"Node {nid} not found in model")
    try:
        xyz = node.get_position()
    except Exception:
        xyz = node.xyz
    return (float(xyz[0]), float(xyz[1]), float(xyz[2]))


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _normalize(v):
    n = _norm(v)
    if n < 1e-12:
        return (1.0, 0.0, 0.0)
    return (v[0] / n, v[1] / n, v[2] / n)


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

class CrossSectionController:
    """Owns one shared clipping plane and bookkeeping state."""

    def __init__(self, plotter):
        self.plotter = plotter
        self._plane = vtk.vtkPlane()
        self._plane.SetOrigin(0.0, 0.0, 0.0)
        self._plane.SetNormal(1.0, 0.0, 0.0)
        self._enabled = False
        # Cached current normal label and slider fraction (0..1) so the
        # dock can rebuild the slider after a model load.
        self._normal_label = "+X"
        self._slider_fraction = 0.5
        # The plane definition currently driving (origin, normal). Set by
        # apply_definition(); the dock keeps a copy too for UI display.
        self._definition: Optional[PlaneDefinition] = None
        # Optional outline actor that visualizes the plane as a translucent
        # rectangle in the scene. Created on demand.
        self._outline_actor = None
        self._outline_visible = False
        self._flipped = False  # toggled by Flip Normal in the dock

    # ----- state queries -----

    @property
    def is_enabled(self) -> bool:
        return self._enabled

    @property
    def plane(self):
        return self._plane

    @property
    def normal_label(self) -> str:
        return self._normal_label

    @property
    def slider_fraction(self) -> float:
        return self._slider_fraction

    @property
    def definition(self) -> Optional[PlaneDefinition]:
        return self._definition

    @property
    def is_flipped(self) -> bool:
        return self._flipped

    def set_slider_fraction(self, frac: float) -> None:
        self._slider_fraction = max(0.0, min(1.0, float(frac)))

    # ----- plane mutators -----

    def set_origin(self, x: float, y: float, z: float) -> None:
        self._plane.SetOrigin(float(x), float(y), float(z))
        self._update_outline()
        if self._enabled:
            self.plotter.render()

    def set_normal(self, nx: float, ny: float, nz: float, label: str = "") -> None:
        self._plane.SetNormal(float(nx), float(ny), float(nz))
        if label:
            self._normal_label = label
        self._update_outline()
        if self._enabled:
            self.plotter.render()

    def apply_definition(self, definition: PlaneDefinition, model=None) -> None:
        """Resolve a PlaneDefinition into (origin, normal) and push it.

        The slider in the dock then scrubs along the resolved normal,
        offsetting the origin to follow the slider fraction.
        """
        origin, normal = definition.resolve(model)
        self._definition = definition
        self._normal_label = (definition.axis_label
                              if definition.method in (METHOD_AXIS, METHOD_COORD_VALUE)
                              else f"({normal[0]:.2g}, {normal[1]:.2g}, {normal[2]:.2g})")
        if self._flipped:
            normal = (-normal[0], -normal[1], -normal[2])
        self._plane.SetNormal(*normal)
        self._plane.SetOrigin(*origin)
        self._update_outline()
        if self._enabled:
            self.plotter.render()

    def flip(self) -> None:
        """Invert the side that gets clipped."""
        self._flipped = not self._flipped
        nx, ny, nz = self._plane.GetNormal()
        self._plane.SetNormal(-nx, -ny, -nz)
        self._update_outline()
        if self._enabled:
            self.plotter.render()

    # ----- enable / disable -----

    def enable(self) -> None:
        if self._enabled:
            return
        self._enabled = True
        self._add_plane_to_all()
        self.plotter.render()

    def disable(self) -> None:
        if not self._enabled:
            return
        self._enabled = False
        self._remove_plane_from_all()
        self._set_outline_visible(False)
        self.plotter.render()

    def reapply(self) -> None:
        """Re-attach the clipping plane to every current actor.

        Called by MainWindow after the scene is rebuilt while cross-section
        is active. Idempotent: removes the plane from any actors that
        already have it before re-adding, so the plane never appears
        twice in a single mapper's clip-plane collection.
        """
        if not self._enabled:
            return
        self._remove_plane_from_all()
        self._add_plane_to_all()

    # ----- outline actor (visual feedback for where the plane sits) -----

    def set_outline_visible(self, visible: bool) -> None:
        self._outline_visible = bool(visible)
        self._set_outline_visible(self._outline_visible and self._enabled)

    def _set_outline_visible(self, on: bool):
        if self._outline_actor is None:
            return
        try:
            self._outline_actor.SetVisibility(bool(on))
            self.plotter.render()
        except Exception:
            pass

    def _ensure_outline_actor(self):
        """Create the translucent plane rectangle once, sized to the model."""
        if self._outline_actor is not None:
            return
        try:
            grid_bounds = self._scene_bounds()
            if grid_bounds is None:
                return
            xmin, xmax, ymin, ymax, zmin, zmax = grid_bounds
            diag = math.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2 + (zmax - zmin) ** 2)
            half = diag * 0.5 if diag > 0 else 1.0
            src = vtk.vtkPlaneSource()
            src.SetOrigin(-half, -half, 0.0)
            src.SetPoint1(half, -half, 0.0)
            src.SetPoint2(-half, half, 0.0)
            src.Update()
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(src.GetOutputPort())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(1.0, 0.85, 0.2)
            actor.GetProperty().SetOpacity(0.18)
            actor.GetProperty().SetEdgeVisibility(True)
            actor.GetProperty().SetEdgeColor(1.0, 0.65, 0.0)
            actor.SetVisibility(False)
            self.plotter.add_actor(actor, name="_cross_section_outline", reset_camera=False)
            self._outline_actor = actor
            self._outline_size = half
        except Exception:
            self._outline_actor = None

    def _scene_bounds(self):
        try:
            return self.plotter.renderer.ComputeVisiblePropBounds()
        except Exception:
            return None

    def _update_outline(self):
        if not self._outline_visible or self._outline_actor is None:
            return
        try:
            ox, oy, oz = self._plane.GetOrigin()
            nx, ny, nz = self._plane.GetNormal()
            # Build a transform that takes the canonical XY square at z=0
            # and rotates+translates it onto the plane.
            up = (0.0, 0.0, 1.0)
            if abs(nx) < 1e-6 and abs(ny) < 1e-6:
                up = (1.0, 0.0, 0.0)  # avoid degenerate case
            # Rodrigues rotation from (0,0,1) to (nx,ny,nz)
            mat = _rotation_matrix_from_z((nx, ny, nz))
            t = vtk.vtkTransform()
            t.PostMultiply()
            m = vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    m.SetElement(i, j, mat[i][j])
            t.SetMatrix(m)
            t.Translate(ox, oy, oz)
            self._outline_actor.SetUserTransform(t)
        except Exception:
            pass

    # ----- internals -----

    def _all_mappers(self):
        for name, actor in list(self.plotter.actors.items()):
            if name == "_cross_section_outline":
                continue
            mapper = actor.GetMapper() if actor is not None else None
            if mapper is None:
                continue
            yield mapper

    def _add_plane_to_all(self) -> None:
        for mapper in self._all_mappers():
            try:
                mapper.AddClippingPlane(self._plane)
            except Exception:
                # Some mappers (e.g., 2D text overlays) don't support
                # clipping planes. Silently skip them.
                pass
        self._ensure_outline_actor()
        self._set_outline_visible(self._outline_visible)

    def _remove_plane_from_all(self) -> None:
        for mapper in self._all_mappers():
            try:
                mapper.RemoveClippingPlane(self._plane)
            except Exception:
                pass


def _rotation_matrix_from_z(target):
    """Rotation matrix that takes (0,0,1) onto unit-vector ``target``."""
    tx, ty, tz = _normalize(target)
    # Rotation axis = z x target; angle = acos(z . target)
    axis = (-ty, tx, 0.0)  # cross((0,0,1), (tx,ty,tz))
    s = math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2])
    c = tz
    if s < 1e-9:
        if c > 0:
            return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        return [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    ux, uy, uz = axis[0] / s, axis[1] / s, axis[2] / s
    angle = math.atan2(s, c)
    cs = math.cos(angle)
    sn = math.sin(angle)
    one_c = 1.0 - cs
    return [
        [cs + ux * ux * one_c,       ux * uy * one_c - uz * sn, ux * uz * one_c + uy * sn],
        [uy * ux * one_c + uz * sn,  cs + uy * uy * one_c,      uy * uz * one_c - ux * sn],
        [uz * ux * one_c - uy * sn,  uz * uy * one_c + ux * sn, cs + uz * uz * one_c],
    ]
