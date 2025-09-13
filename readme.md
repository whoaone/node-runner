# Node Runner v1.1.0

A lightweight pre-processor for creating, editing, and visualizing Nastran models. Built with Python, PySide6, and PyVista.

## Changelog for v1.1.0

This is a significant feature release that adds full support for `CBUSH` elements, enhanced visualization tools, a powerful coordinate system editor, and major quality-of-life improvements.

### ✨ New Features
- **Full CBUSH Element Support:**
  - Create `CBUSH` elements with a new dedicated dialog.
  - Orient bushes using a vector, a third node, or a coordinate system.
  - Edit existing `CBUSH` elements.
  - Visualize `CBUSH` orientation with a 3D triad glyph.
- **Enhanced Beam Visualization:**
  - `CBEAM` and `CBAR` elements now also have 3D orientation glyphs.
  - The Element Editor now supports modifying beam orientation.
- **Default Coordinate Systems:**
  - The application now guarantees the presence of three undeletable global coordinate systems (Rectangular, Cylindrical, and Spherical) in every model.
- **File Import/Append:**
  - A new "Import" function allows appending cards from a `.bdf` file to the current model, with logic to overwrite conflicting entities.

### 🐛 Bug Fixes & Refinements
- Rewrote the file parser to robustly handle both full and partial Nastran decks.
- Fixed numerous bugs related to UI dialogs, rendering, and entity creation.
- Corrected the display logic to ensure free-floating nodes are always visible.
- Made title-parsing from comments more intelligent for a cleaner model tree.

---

## LICENSE

This project is licensed under the Apache-2.0 License.