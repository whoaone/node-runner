# Node Runner v2.0.0

A lightweight pre-processor for creating, editing, and visualizing Nastran models. Built with Python, PySide6, and PyVista.

## Changelog for v2.0.0

This is a major release with coordinate system protection, enhanced display options, improved view controls, geometry transforms, and numerous quality-of-life improvements.

### ✨ New Features

- **Coordinate System Protection:**
  - Global Rectangular (CID 0), Cylindrical (CID 1), and Spherical (CID 2) coordinate systems are now immutable — they cannot be deleted or overwritten by any tool, import, or generator.
  - Protection enforced at the command level (`AddCoordCommand`, `DeleteCoordCommand`) and guaranteed after every fresh model/generator creation.

- **CSys Display Options (Display Tab):**
  - Master toggle to show/hide all coordinate system actors.
  - Separate toggle for global coordinate systems (CID 0, 1, 2).
  - Label visibility toggle for coordinate system names.
  - Axis size slider for scaling coordinate triad arrows.
  - Centralized `_refresh_coord_actors()` pipeline handles full actor lifecycle.

- **Isometric View Cube Icons:**
  - Standard view buttons (Top, Bottom, Front, Back, Left, Right, Isometric) now display 3D isometric cube icons generated with QPainter.
  - Blue highlight for positive directions, orange for negative.
  - Single-row layout with separators between paired views.

- **Geometry Transform & Operations:**
  - Translate, Rotate, Mirror, Scale, and Copy geometry entities.
  - Split, Offset, Project, and Fillet curve operations.
  - Full geometry modify support for points, curves, and surfaces.

- **Load System Enhancements:**
  - Load Set Manager dialog for viewing/removing individual load entries.
  - `CreateLoadDialog` now always appends (never overwrites entire SID on edit).
  - Tree right-click context menu: "Edit Load Set" opens manager, "Add to Load Set" opens create dialog with pre-filled SID.

- **Full CBUSH Element Support:**
  - Create `CBUSH` elements with a dedicated dialog.
  - Orient bushes using a vector, third node, or coordinate system.
  - Edit existing `CBUSH` elements.
  - Visualize `CBUSH` orientation with a 3D triad glyph.

- **Enhanced Beam Visualization:**
  - `CBEAM` and `CBAR` elements now display 3D orientation glyphs.
  - Element Editor supports modifying beam orientation.

- **File Import/Append:**
  - "Import" function appends cards from a `.bdf` file to the current model, with logic to overwrite conflicting entities.

### 🐛 Bug Fixes & Refinements
- Fixed Group "List Contents" table showing empty cells (classic Qt `setSortingEnabled` bug).
- Rewrote the file parser to robustly handle both full and partial Nastran decks.
- Fixed numerous bugs related to UI dialogs, rendering, and entity creation.
- Corrected the display logic to ensure free-floating nodes are always visible.
- Made title-parsing from comments more intelligent for a cleaner model tree.
- Fixed global CSys toggles not working (actors now properly created on demand).
- Eliminated duplicate render calls in `_update_plot_visibility`.

---

## LICENSE

This project is licensed under the Apache-2.0 License.
