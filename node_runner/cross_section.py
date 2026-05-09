"""Cross-section / clipping plane controller.

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
"""

from __future__ import annotations

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


class CrossSectionController:
    """Owns one shared clipping plane and bookkeeping state."""

    def __init__(self, plotter):
        self.plotter = plotter
        self._plane = vtk.vtkPlane()
        self._plane.SetOrigin(0.0, 0.0, 0.0)
        self._plane.SetNormal(1.0, 0.0, 0.0)
        self._enabled = False
        # Cached current normal label and slider fraction (0..1) so the
        # dialog can rebuild the slider after a model load.
        self._normal_label = "+X"
        self._slider_fraction = 0.5

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

    def set_slider_fraction(self, frac: float) -> None:
        self._slider_fraction = max(0.0, min(1.0, float(frac)))

    # ----- plane mutators -----

    def set_origin(self, x: float, y: float, z: float) -> None:
        self._plane.SetOrigin(float(x), float(y), float(z))
        if self._enabled:
            self.plotter.render()

    def set_normal(self, nx: float, ny: float, nz: float, label: str = "") -> None:
        self._plane.SetNormal(float(nx), float(ny), float(nz))
        if label:
            self._normal_label = label
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

    # ----- internals -----

    def _all_mappers(self):
        for actor in list(self.plotter.actors.values()):
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

    def _remove_plane_from_all(self) -> None:
        for mapper in self._all_mappers():
            try:
                mapper.RemoveClippingPlane(self._plane)
            except Exception:
                pass
