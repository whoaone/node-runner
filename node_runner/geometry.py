"""Geometry primitives for the preprocessor.

Geometry entities are independent of FEM nodes/elements. They provide
construction curves and surfaces that can later be meshed into FEM entities.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class GeometryPoint:
    id: int
    xyz: np.ndarray  # shape (3,)


@dataclass
class GeometryLine:
    id: int
    start_point_id: int
    end_point_id: int
    _start_xyz: np.ndarray = field(default=None, repr=False)
    _end_xyz: np.ndarray = field(default=None, repr=False)

    def evaluate(self, t: float) -> np.ndarray:
        return self._start_xyz + t * (self._end_xyz - self._start_xyz)

    def evaluate_array(self, n_points: int) -> np.ndarray:
        ts = np.linspace(0, 1, n_points)
        return np.outer(1 - ts, self._start_xyz) + np.outer(ts, self._end_xyz)

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self._end_xyz - self._start_xyz))


@dataclass
class GeometryArc:
    id: int
    start_point_id: int
    mid_point_id: int
    end_point_id: int
    _center: np.ndarray = field(default=None, repr=False)
    _radius: float = field(default=0.0, repr=False)
    _normal: np.ndarray = field(default=None, repr=False)
    _start_angle: float = field(default=0.0, repr=False)
    _sweep_angle: float = field(default=0.0, repr=False)
    _basis_u: np.ndarray = field(default=None, repr=False)
    _basis_v: np.ndarray = field(default=None, repr=False)

    def compute_arc_params(self, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray):
        """Compute arc centre, radius, and angular span from three 3-D points."""
        # Plane normal
        v12 = p2 - p1
        v13 = p3 - p1
        normal = np.cross(v12, v13)
        norm_len = np.linalg.norm(normal)
        if norm_len < 1e-12:
            raise ValueError("Arc points are collinear.")
        normal = normal / norm_len
        self._normal = normal

        # Circumscribed circle via perpendicular bisectors in the plane
        m1 = (p1 + p2) / 2.0
        m2 = (p2 + p3) / 2.0
        d1 = np.cross(v12, normal)
        d2 = np.cross(p3 - p2, normal)

        # Solve m1 + s*d1 = m2 + t*d2 (in least-squares sense)
        A = np.column_stack([d1, -d2])
        b = m2 - m1
        # Use pseudoinverse for robustness
        params, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        center = m1 + params[0] * d1
        self._center = center
        self._radius = float(np.linalg.norm(p1 - center))

        # Local basis in the arc plane
        u = (p1 - center)
        u_len = np.linalg.norm(u)
        if u_len < 1e-12:
            raise ValueError("Start point coincides with arc centre.")
        self._basis_u = u / u_len
        self._basis_v = np.cross(normal, self._basis_u)

        # Angles
        self._start_angle = 0.0
        angle_mid = math.atan2(
            np.dot(p2 - center, self._basis_v),
            np.dot(p2 - center, self._basis_u))
        angle_end = math.atan2(
            np.dot(p3 - center, self._basis_v),
            np.dot(p3 - center, self._basis_u))

        # Ensure sweep goes through the midpoint
        # Normalise angles to [0, 2*pi)
        angle_mid = angle_mid % (2 * math.pi)
        angle_end = angle_end % (2 * math.pi)

        # Check both possible sweep directions
        if angle_mid <= angle_end:
            self._sweep_angle = angle_end
        else:
            self._sweep_angle = angle_end - 2 * math.pi  # negative sweep

        # Verify midpoint is within the sweep
        if self._sweep_angle > 0:
            if not (0 <= angle_mid <= self._sweep_angle):
                self._sweep_angle = angle_end - 2 * math.pi
        else:
            if not (self._sweep_angle <= angle_mid <= 0):
                self._sweep_angle = angle_end

    def evaluate(self, t: float) -> np.ndarray:
        angle = self._start_angle + t * self._sweep_angle
        return (self._center +
                self._radius * math.cos(angle) * self._basis_u +
                self._radius * math.sin(angle) * self._basis_v)

    def evaluate_array(self, n_points: int) -> np.ndarray:
        ts = np.linspace(0, 1, n_points)
        angles = self._start_angle + ts * self._sweep_angle
        return (self._center[np.newaxis, :] +
                self._radius * np.cos(angles)[:, np.newaxis] * self._basis_u[np.newaxis, :] +
                self._radius * np.sin(angles)[:, np.newaxis] * self._basis_v[np.newaxis, :])


@dataclass
class GeometryCircle:
    id: int
    center_point_id: int
    radius: float
    normal: np.ndarray  # plane normal, shape (3,)
    _center_xyz: np.ndarray = field(default=None, repr=False)
    _basis_u: np.ndarray = field(default=None, repr=False)
    _basis_v: np.ndarray = field(default=None, repr=False)

    def compute_basis(self):
        """Compute u, v basis vectors in the circle plane from the normal."""
        n = self.normal / np.linalg.norm(self.normal)
        # Pick a non-parallel reference vector
        ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(n, ref)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        self._basis_u = np.cross(n, ref)
        self._basis_u /= np.linalg.norm(self._basis_u)
        self._basis_v = np.cross(n, self._basis_u)

    def evaluate(self, t: float) -> np.ndarray:
        angle = 2 * math.pi * t
        return (self._center_xyz +
                self.radius * math.cos(angle) * self._basis_u +
                self.radius * math.sin(angle) * self._basis_v)

    def evaluate_array(self, n_points: int) -> np.ndarray:
        ts = np.linspace(0, 1, n_points, endpoint=False)
        angles = 2 * math.pi * ts
        return (self._center_xyz[np.newaxis, :] +
                self.radius * np.cos(angles)[:, np.newaxis] * self._basis_u[np.newaxis, :] +
                self.radius * np.sin(angles)[:, np.newaxis] * self._basis_v[np.newaxis, :])


@dataclass
class GeometrySurface:
    id: int
    boundary_curve_ids: list  # ordered curve IDs forming the boundary


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

class GeometryStore:
    """Container for all geometry entities with ID management."""

    def __init__(self):
        self.points: dict[int, GeometryPoint] = {}
        self.lines: dict[int, GeometryLine] = {}
        self.arcs: dict[int, GeometryArc] = {}
        self.circles: dict[int, GeometryCircle] = {}
        self.surfaces: dict[int, GeometrySurface] = {}
        self._next_point_id = 1
        self._next_curve_id = 1
        self._next_surface_id = 1

    # --- ID management ---

    def next_point_id(self) -> int:
        while self._next_point_id in self.points:
            self._next_point_id += 1
        return self._next_point_id

    def next_curve_id(self) -> int:
        all_ids = set(self.lines) | set(self.arcs) | set(self.circles)
        while self._next_curve_id in all_ids:
            self._next_curve_id += 1
        return self._next_curve_id

    def next_surface_id(self) -> int:
        while self._next_surface_id in self.surfaces:
            self._next_surface_id += 1
        return self._next_surface_id

    # --- Add methods ---

    def add_point(self, xyz, point_id=None) -> GeometryPoint:
        if point_id is None:
            point_id = self.next_point_id()
        pt = GeometryPoint(id=point_id, xyz=np.array(xyz, dtype=float))
        self.points[point_id] = pt
        return pt

    def add_line(self, start_id: int, end_id: int, line_id=None) -> GeometryLine:
        if line_id is None:
            line_id = self.next_curve_id()
        line = GeometryLine(id=line_id, start_point_id=start_id,
                            end_point_id=end_id)
        line._start_xyz = self.points[start_id].xyz.copy()
        line._end_xyz = self.points[end_id].xyz.copy()
        self.lines[line_id] = line
        return line

    def add_arc(self, start_id: int, mid_id: int, end_id: int,
                arc_id=None) -> GeometryArc:
        if arc_id is None:
            arc_id = self.next_curve_id()
        arc = GeometryArc(id=arc_id, start_point_id=start_id,
                          mid_point_id=mid_id, end_point_id=end_id)
        arc.compute_arc_params(
            self.points[start_id].xyz,
            self.points[mid_id].xyz,
            self.points[end_id].xyz)
        self.arcs[arc_id] = arc
        return arc

    def add_circle(self, center_id: int, radius: float, normal,
                   circle_id=None) -> GeometryCircle:
        if circle_id is None:
            circle_id = self.next_curve_id()
        circle = GeometryCircle(id=circle_id, center_point_id=center_id,
                                radius=radius,
                                normal=np.array(normal, dtype=float))
        circle._center_xyz = self.points[center_id].xyz.copy()
        circle.compute_basis()
        self.circles[circle_id] = circle
        return circle

    def add_surface(self, curve_ids, surface_id=None) -> GeometrySurface:
        if surface_id is None:
            surface_id = self.next_surface_id()
        surf = GeometrySurface(id=surface_id,
                               boundary_curve_ids=list(curve_ids))
        self.surfaces[surface_id] = surf
        return surf

    # --- Lookup ---

    def get_curve(self, curve_id: int):
        """Look up a curve by ID across all curve types."""
        if curve_id in self.lines:
            return self.lines[curve_id]
        if curve_id in self.arcs:
            return self.arcs[curve_id]
        if curve_id in self.circles:
            return self.circles[curve_id]
        raise KeyError(f"Curve {curve_id} not found")

    def evaluate_curve(self, curve_id: int, n_points: int) -> np.ndarray:
        """Return *n_points* equally spaced XYZ positions along a curve."""
        return self.get_curve(curve_id).evaluate_array(n_points)

    def all_curve_ids(self) -> list[int]:
        """Return sorted list of all curve IDs across types."""
        return sorted(set(self.lines) | set(self.arcs) | set(self.circles))

    # --- Dependency lookup (for cascade deletion) ---

    def curves_referencing_point(self, point_id: int) -> list[int]:
        """Return curve IDs that reference the given point ID."""
        result = []
        for lid, line in self.lines.items():
            if point_id in (line.start_point_id, line.end_point_id):
                result.append(lid)
        for aid, arc in self.arcs.items():
            if point_id in (arc.start_point_id, arc.mid_point_id, arc.end_point_id):
                result.append(aid)
        for cid, circle in self.circles.items():
            if point_id == circle.center_point_id:
                result.append(cid)
        return result

    def surfaces_referencing_curve(self, curve_id: int) -> list[int]:
        """Return surface IDs whose boundary references the given curve ID."""
        return [sid for sid, surf in self.surfaces.items()
                if curve_id in surf.boundary_curve_ids]

    def surfaces_referencing_curves(self, curve_ids: set) -> list[int]:
        """Return surface IDs referencing any of the given curve IDs."""
        return [sid for sid, surf in self.surfaces.items()
                if curve_ids.intersection(surf.boundary_curve_ids)]

    # --- Update methods (for geometry modify operations) ---

    def update_point(self, point_id: int, new_xyz) -> None:
        """Move a point to new coordinates and refresh all referencing curves."""
        self.points[point_id].xyz = np.array(new_xyz, dtype=float)
        self._refresh_curves_for_point(point_id)

    def _refresh_curves_for_point(self, point_id: int) -> None:
        """Recompute cached data in all curves referencing this point."""
        for cid in self.curves_referencing_point(point_id):
            if cid in self.lines:
                line = self.lines[cid]
                line._start_xyz = self.points[line.start_point_id].xyz.copy()
                line._end_xyz = self.points[line.end_point_id].xyz.copy()
            elif cid in self.arcs:
                arc = self.arcs[cid]
                arc.compute_arc_params(
                    self.points[arc.start_point_id].xyz,
                    self.points[arc.mid_point_id].xyz,
                    self.points[arc.end_point_id].xyz)
            elif cid in self.circles:
                circle = self.circles[cid]
                circle._center_xyz = self.points[circle.center_point_id].xyz.copy()
                circle.compute_basis()

    def update_line_refs(self, line_id: int, start_id: int, end_id: int) -> None:
        """Change which points a line references and refresh caches."""
        line = self.lines[line_id]
        line.start_point_id = start_id
        line.end_point_id = end_id
        line._start_xyz = self.points[start_id].xyz.copy()
        line._end_xyz = self.points[end_id].xyz.copy()

    def update_arc_refs(self, arc_id: int, start_id: int, mid_id: int,
                        end_id: int) -> None:
        """Change which points an arc references and recompute arc params."""
        arc = self.arcs[arc_id]
        arc.start_point_id = start_id
        arc.mid_point_id = mid_id
        arc.end_point_id = end_id
        arc.compute_arc_params(
            self.points[start_id].xyz,
            self.points[mid_id].xyz,
            self.points[end_id].xyz)

    def update_circle(self, circle_id: int, center_id: int = None,
                      radius: float = None, normal=None) -> None:
        """Update circle parameters and recompute basis."""
        circle = self.circles[circle_id]
        if center_id is not None:
            circle.center_point_id = center_id
            circle._center_xyz = self.points[center_id].xyz.copy()
        if radius is not None:
            circle.radius = radius
        if normal is not None:
            circle.normal = np.array(normal, dtype=float)
        circle.compute_basis()

    def update_surface_boundaries(self, surface_id: int, curve_ids: list) -> None:
        """Replace the boundary curve list of a surface."""
        self.surfaces[surface_id].boundary_curve_ids = list(curve_ids)

    @property
    def is_empty(self) -> bool:
        return not (self.points or self.lines or self.arcs or
                    self.circles or self.surfaces)

    def clear(self):
        self.points.clear()
        self.lines.clear()
        self.arcs.clear()
        self.circles.clear()
        self.surfaces.clear()
        self._next_point_id = 1
        self._next_curve_id = 1
        self._next_surface_id = 1
