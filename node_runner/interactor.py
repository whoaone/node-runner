from __future__ import annotations

import time
from typing import TYPE_CHECKING

import numpy as np
import vtk

if TYPE_CHECKING:
    from node_runner.mainwindow import MainWindow


class ClickAndDragInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, main_window: MainWindow):
        super().__init__()
        self.main_window = main_window
        self.mouse_press_time = 0
        self.mouse_press_pos = (0, 0)
        self.click_duration_threshold = 0.2
        self.drag_distance_threshold = 5
        self._area_select_active = False

        # Selection mode: 'box', 'polygon', 'circle'
        self._selection_mode = 'box'

        # Polygon mode state
        self._polygon_points = []  # list of (x, y) screen coords
        self._polygon_active = False

        # Circle mode state
        self._circle_center = None  # (x, y) screen coords

        # Middle-button navigation state (Femap-style)
        self._nav_mode = None          # 'rotate' or 'pan' when middle button active
        self._middle_btn_active = False

        self.AddObserver(vtk.vtkCommand.LeftButtonPressEvent, self.left_button_press_event)
        self.AddObserver(vtk.vtkCommand.LeftButtonReleaseEvent, self.left_button_release_event)
        self.AddObserver(vtk.vtkCommand.MiddleButtonPressEvent, self.middle_button_press_event)
        self.AddObserver(vtk.vtkCommand.MiddleButtonReleaseEvent, self.middle_button_release_event)
        self.AddObserver(vtk.vtkCommand.MouseMoveEvent, self.mouse_move_event)
        self.AddObserver(vtk.vtkCommand.KeyPressEvent, self.key_press_event)
        self.AddObserver(vtk.vtkCommand.LeftButtonDoubleClickEvent, self.double_click_event)

    def set_selection_mode(self, mode):
        """Set the selection mode: 'box', 'polygon', 'circle'."""
        self._selection_mode = mode
        self.main_window._pick_selection_mode = mode  # Sync back to mainwindow
        # Reset any in-progress polygon
        self._polygon_points = []
        self._polygon_active = False
        self._circle_center = None
        self.main_window._hide_rubber_band()
        self.main_window._update_status(f"Pick mode: {mode.capitalize()}")

    def key_press_event(self, obj, event):
        key = self.GetInteractor().GetKeySym()
        if key == "Return":
            if self._polygon_active and len(self._polygon_points) >= 3:
                # Finalize polygon selection
                self._finalize_polygon()
            else:
                self.main_window._request_accept_picking_mode()
        elif key == "Escape":
            if self._polygon_active:
                self._polygon_points = []
                self._polygon_active = False
                self.main_window._hide_rubber_band()
                self.main_window._update_status("Polygon selection cancelled.")
            else:
                self.main_window._request_disable_picking_mode()
        elif key in ('b', 'B'):
            self.set_selection_mode('box')
        elif key in ('p', 'P'):
            self.set_selection_mode('polygon')
        elif key in ('c', 'C'):
            self.set_selection_mode('circle')
        else:
            self.OnKeyPress()

    def double_click_event(self, obj, event):
        """Handle double-click: finalize polygon, or clear highlights and refresh view."""
        if self._polygon_active and len(self._polygon_points) >= 3:
            self._finalize_polygon()
        else:
            # Clear any selection highlights and refresh view (no auto-zoom)
            if hasattr(self.main_window, '_clear_highlight'):
                self.main_window._clear_highlight()

    def left_button_press_event(self, obj, event):
        self.mouse_press_time = time.time()
        self.mouse_press_pos = self.GetInteractor().GetEventPosition()

        # In picking mode with a selection dialog, intercept for area selection
        if (self.main_window.active_selection_dialog
                and not self.main_window.picking_target_callback):

            if self._selection_mode == 'polygon' and self._polygon_active:
                # Adding another polygon vertex - handled in release
                self._area_select_active = True
                return

            if self._selection_mode == 'circle':
                # Start circle: store center
                self._circle_center = self.mouse_press_pos

            self._area_select_active = True
            return
        self.OnLeftButtonDown()

    def middle_button_press_event(self, obj, event):
        """Femap-style: Middle=Rotate, Ctrl+Middle=Pan. Always works."""
        self._middle_btn_active = True
        x, y = self.GetInteractor().GetEventPosition()
        self.FindPokedRenderer(x, y)
        ctrl_key = self.GetInteractor().GetControlKey()
        if ctrl_key:
            self._nav_mode = 'pan'
            self.OnMiddleButtonDown()
        else:
            self._nav_mode = 'rotate'
            self.OnLeftButtonDown()

    def middle_button_release_event(self, obj, event):
        if self._nav_mode == 'rotate':
            self.OnLeftButtonUp()
        elif self._nav_mode == 'pan':
            self.OnMiddleButtonUp()
        self._nav_mode = None
        self._middle_btn_active = False

    def mouse_move_event(self, obj, event):
        # Middle-button navigation always has priority
        if self._nav_mode:
            self.OnMouseMove()
            return

        if self._area_select_active:
            current_pos = self.GetInteractor().GetEventPosition()
            dist = np.linalg.norm(
                np.array(current_pos) - np.array(self.mouse_press_pos))
            if dist > self.drag_distance_threshold:
                if self._selection_mode == 'box':
                    self.main_window._update_rubber_band(
                        self.mouse_press_pos, current_pos)
                elif self._selection_mode == 'circle' and self._circle_center:
                    radius = dist
                    self.main_window._update_circle_overlay(
                        self._circle_center, radius)
            return

        # Polygon live preview: show edges + cursor line while building polygon
        if self._polygon_active and self._polygon_points:
            current_pos = self.GetInteractor().GetEventPosition()
            self.main_window._update_polygon_overlay(
                self._polygon_points, current_pos)
            return

        self.OnMouseMove()

    def left_button_release_event(self, obj, event):
        current_pos = self.GetInteractor().GetEventPosition()

        if self._area_select_active:
            self._area_select_active = False
            self.main_window._hide_rubber_band()

            is_click = (
                time.time() - self.mouse_press_time
                < self.click_duration_threshold
                and np.linalg.norm(
                    np.array(current_pos) - np.array(self.mouse_press_pos))
                < self.drag_distance_threshold
            )

            if self._selection_mode == 'polygon':
                if is_click:
                    # Add point to polygon
                    self._polygon_points.append(current_pos)
                    self._polygon_active = True
                    n = len(self._polygon_points)
                    self.main_window._update_polygon_overlay(
                        self._polygon_points)
                    self.main_window._update_status(
                        f"Polygon: {n} point{'s' if n != 1 else ''} - "
                        f"double-click or Enter to finish, Esc to cancel")
                    return
                else:
                    # Drag in polygon mode - ignore
                    return
            elif self._selection_mode == 'circle':
                if is_click:
                    # Single click in circle mode - do a normal pick
                    self.perform_pick(current_pos)
                else:
                    # Drag completed - do circle area pick
                    self.perform_circle_pick(self._circle_center, current_pos)
                self._circle_center = None
                return
            else:
                # Box mode
                if is_click:
                    self.perform_pick(current_pos)
                else:
                    self.perform_area_pick(self.mouse_press_pos, current_pos)
                return

        is_click = (
            time.time() - self.mouse_press_time < self.click_duration_threshold
            and np.linalg.norm(
                np.array(current_pos) - np.array(self.mouse_press_pos))
            < self.drag_distance_threshold
        )

        self.OnLeftButtonUp()

        if is_click:
            self.perform_pick(self.GetInteractor().GetEventPosition())

    def _finalize_polygon(self):
        """Complete polygon selection and pick entities inside."""
        if len(self._polygon_points) < 3:
            self._polygon_points = []
            self._polygon_active = False
            self.main_window._hide_rubber_band()
            return

        polygon = np.array(self._polygon_points, dtype=float)
        self._perform_polygon_pick(polygon)
        self._polygon_points = []
        self._polygon_active = False
        self.main_window._hide_rubber_band()

    def _perform_polygon_pick(self, polygon):
        """Select entities whose screen projections fall inside the polygon."""
        renderer = self.main_window.plotter.renderer
        if not renderer:
            return

        dialog = self.main_window.active_selection_dialog
        if not dialog:
            return

        entity_type = dialog.entity_type
        if not self.main_window.current_generator:
            return

        shift_key = self.GetInteractor().GetShiftKey()
        ctrl_key = self.GetInteractor().GetControlKey()

        selected_ids = set()

        try:
            if entity_type == 'Node':
                nids = self.main_window.current_node_ids_sorted
                if not nids:
                    return
                model = self.main_window.current_generator.model
                coords_3d = np.array(
                    [model.nodes[nid].get_position() for nid in nids])
                screen = self._project_to_screen(coords_3d, renderer)
                mask = self._points_in_polygon(screen, polygon)
                selected_ids = {nids[i] for i in np.where(mask)[0]}

            elif entity_type == 'Element':
                grid = self.main_window.current_grid
                if grid is None or 'EID' not in grid.cell_data:
                    return
                eids = grid.cell_data['EID']
                centers = grid.cell_centers().points
                screen = self._project_to_screen(centers, renderer)
                mask = self._points_in_polygon(screen, polygon)
                selected_ids = {int(eids[i]) for i in np.where(mask)[0]}

        except Exception as e:
            self.main_window._update_status(
                f"Polygon selection error: {e}", is_error=True)
            return

        self._apply_selection(dialog, selected_ids, entity_type,
                              shift_key, ctrl_key, "polygon")

    def perform_circle_pick(self, center_pos, edge_pos):
        """Select entities within a circle defined by center and edge point."""
        renderer = self.main_window.plotter.renderer
        if not renderer:
            return

        dialog = self.main_window.active_selection_dialog
        if not dialog:
            return

        entity_type = dialog.entity_type
        if not self.main_window.current_generator:
            return

        radius = np.linalg.norm(
            np.array(edge_pos, dtype=float) - np.array(center_pos, dtype=float))
        if radius < 3:
            return

        center = np.array(center_pos, dtype=float)
        shift_key = self.GetInteractor().GetShiftKey()
        ctrl_key = self.GetInteractor().GetControlKey()

        selected_ids = set()

        try:
            if entity_type == 'Node':
                nids = self.main_window.current_node_ids_sorted
                if not nids:
                    return
                model = self.main_window.current_generator.model
                coords_3d = np.array(
                    [model.nodes[nid].get_position() for nid in nids])
                screen = self._project_to_screen(coords_3d, renderer)
                dists = np.linalg.norm(screen - center, axis=1)
                mask = dists <= radius
                selected_ids = {nids[i] for i in np.where(mask)[0]}

            elif entity_type == 'Element':
                grid = self.main_window.current_grid
                if grid is None or 'EID' not in grid.cell_data:
                    return
                eids = grid.cell_data['EID']
                centers = grid.cell_centers().points
                screen = self._project_to_screen(centers, renderer)
                dists = np.linalg.norm(screen - center, axis=1)
                mask = dists <= radius
                selected_ids = {int(eids[i]) for i in np.where(mask)[0]}

        except Exception as e:
            self.main_window._update_status(
                f"Circle selection error: {e}", is_error=True)
            return

        self._apply_selection(dialog, selected_ids, entity_type,
                              shift_key, ctrl_key, "circle")

    def _apply_selection(self, dialog, selected_ids, entity_type,
                         shift_key, ctrl_key, mode_name):
        """Apply selection results with modifier key logic.

        Shift/Ctrl keyboard modifiers override the dialog's radio state.
        When no modifier is held, the dialog's Add/Remove/Exclude radio
        determines the behavior.
        """
        if not selected_ids:
            self.main_window._update_status(
                f"No entities found in {mode_name} selection.")
            return

        if ctrl_key:
            dialog.remove_selection(selected_ids)
            verb = "Removed"
        elif shift_key:
            dialog.add_selection(selected_ids)
            verb = "Added"
        else:
            # Respect the dialog's action mode radio (Add/Remove/Exclude)
            action_mode = getattr(dialog, 'get_action_mode', lambda: 'add')()
            if action_mode == 'remove':
                dialog.remove_selection(selected_ids)
                verb = "Removed"
            elif action_mode == 'exclude':
                if hasattr(dialog, '_excluded_ids'):
                    valid = dialog.all_entity_ids.intersection(set(selected_ids))
                    dialog._excluded_ids.update(valid)
                    dialog._refresh_list()
                    dialog._show_selection()
                verb = "Excluded"
            else:  # 'add' (default)
                dialog.selected_ids.clear()
                dialog.add_selection(selected_ids)
                verb = "Selected"

        count = len(selected_ids)
        self.main_window._update_status(
            f"{verb} {count} {entity_type.lower()}"
            f"{'s' if count != 1 else ''}")
        self.main_window._highlight_entities(
            dialog.entity_type, dialog.get_selected_ids())

    @staticmethod
    def _points_in_polygon(points_2d, polygon):
        """Ray-casting point-in-polygon test (vectorized with numpy).

        Args:
            points_2d: (N, 2) array of screen coordinates to test
            polygon: (M, 2) array of polygon vertex coordinates

        Returns:
            (N,) boolean array - True for points inside the polygon
        """
        n_poly = len(polygon)
        inside = np.zeros(len(points_2d), dtype=bool)
        j = n_poly - 1
        for i in range(n_poly):
            yi, yj = polygon[i, 1], polygon[j, 1]
            xi, xj = polygon[i, 0], polygon[j, 0]

            cond1 = (yi > points_2d[:, 1]) != (yj > points_2d[:, 1])
            # Avoid division by zero
            dy = yj - yi
            if abs(dy) < 1e-12:
                j = i
                continue
            slope = (xj - xi) * (points_2d[:, 1] - yi) / dy + xi
            cond2 = points_2d[:, 0] < slope
            inside ^= (cond1 & cond2)
            j = i
        return inside

    def perform_pick(self, click_pos):
        if not (renderer := self.main_window.plotter.renderer): return

        entity_type = self.main_window.current_selection_type
        if not entity_type:
            entity_type = "Node"

        try:
            entity_id = -1
            if entity_type == 'Node':
                picker = vtk.vtkPointPicker()
                picker.Pick(click_pos[0], click_pos[1], 0, renderer)

                point_id_local = picker.GetPointId()
                picked_actor = picker.GetActor()

                if picked_actor and point_id_local != -1:
                    dataset = picked_actor.GetMapper().GetInput()
                    if dataset and 'vtkOriginalPointIds' in dataset.point_data:
                        original_point_id = int(dataset.point_data['vtkOriginalPointIds'][point_id_local])
                        if 0 <= original_point_id < len(self.main_window.current_node_ids_sorted):
                            entity_id = self.main_window.current_node_ids_sorted[original_point_id]

            elif entity_type == 'Element':
                picker = vtk.vtkCellPicker()
                picker.Pick(click_pos[0], click_pos[1], 0, renderer)
                cell_id_local = picker.GetCellId()
                picked_actor = picker.GetActor()

                if picked_actor and cell_id_local != -1:
                    dataset = picked_actor.GetMapper().GetInput()
                    if dataset and 'EID' in dataset.cell_data:
                        entity_id = int(dataset.cell_data["EID"][cell_id_local])

            if entity_id == -1:
                 self.main_window._update_status(f"No {entity_type.lower()} found at this location.")
                 return

            if self.main_window.picking_target_callback:
                self.main_window.picking_target_callback(str(entity_id))
                self.main_window._set_picking_mode("Node", False)
                self.main_window.picking_target_callback = None
                return

            if self.main_window.active_selection_dialog:
                dialog = self.main_window.active_selection_dialog
                dialog.toggle_selection(entity_id)
                status_verb = 'Selected' if entity_id in dialog.selected_ids else 'Deselected'
                self.main_window._update_status(f"{status_verb} {entity_type}: {entity_id}")
                self.main_window._highlight_entities(dialog.entity_type, dialog.get_selected_ids())

        except (IndexError, TypeError, KeyError) as e:
            self.main_window._update_status(f"Could not select {entity_type.lower()}: {e}", is_error=True)

    def perform_area_pick(self, start_pos, end_pos):
        """Select entities whose screen projections fall within the drag rectangle."""
        renderer = self.main_window.plotter.renderer
        if not renderer:
            return

        dialog = self.main_window.active_selection_dialog
        if not dialog:
            return

        entity_type = dialog.entity_type
        if not self.main_window.current_generator:
            return

        # Rectangle bounds in VTK display coords
        x0, x1 = min(start_pos[0], end_pos[0]), max(start_pos[0], end_pos[0])
        y0, y1 = min(start_pos[1], end_pos[1]), max(start_pos[1], end_pos[1])

        shift_key = self.GetInteractor().GetShiftKey()
        ctrl_key = self.GetInteractor().GetControlKey()

        selected_ids = set()

        try:
            if entity_type == 'Node':
                nids = self.main_window.current_node_ids_sorted
                if not nids:
                    return
                model = self.main_window.current_generator.model
                coords_3d = np.array(
                    [model.nodes[nid].get_position() for nid in nids])
                screen = self._project_to_screen(coords_3d, renderer)
                mask = ((screen[:, 0] >= x0) & (screen[:, 0] <= x1)
                        & (screen[:, 1] >= y0) & (screen[:, 1] <= y1))
                selected_ids = {nids[i] for i in np.where(mask)[0]}

            elif entity_type == 'Element':
                grid = self.main_window.current_grid
                if grid is None or 'EID' not in grid.cell_data:
                    return
                eids = grid.cell_data['EID']
                centers = grid.cell_centers().points
                screen = self._project_to_screen(centers, renderer)
                mask = ((screen[:, 0] >= x0) & (screen[:, 0] <= x1)
                        & (screen[:, 1] >= y0) & (screen[:, 1] <= y1))
                selected_ids = {int(eids[i]) for i in np.where(mask)[0]}

        except Exception as e:
            self.main_window._update_status(
                f"Area selection error: {e}", is_error=True)
            return

        self._apply_selection(dialog, selected_ids, entity_type,
                              shift_key, ctrl_key, "box")

    def _project_to_screen(self, coords_3d, renderer):
        """Project 3D world coordinates to 2D display coordinates (vectorized)."""
        camera = renderer.GetActiveCamera()
        aspect = renderer.GetTiledAspectRatio()
        vtk_mat = camera.GetCompositeProjectionTransformMatrix(aspect, -1, 1)
        m = np.array([[vtk_mat.GetElement(i, j) for j in range(4)]
                       for i in range(4)])

        n = coords_3d.shape[0]
        pts = np.hstack([coords_3d, np.ones((n, 1))])
        clip = pts @ m.T

        w = clip[:, 3:4].copy()
        w[w == 0] = 1e-10
        ndc = clip[:, :2] / w

        size = renderer.GetRenderWindow().GetSize()
        return np.column_stack([
            (ndc[:, 0] + 1) * 0.5 * size[0],
            (ndc[:, 1] + 1) * 0.5 * size[1],
        ])
