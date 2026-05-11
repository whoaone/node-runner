# Node Runner v3.1.1

A lightweight pre-processor for creating, editing, and visualizing Nastran models. Built with Python, PySide6, and PyVista.

## Running

The simplest path on Windows is to grab `Node.Runner.exe` from the latest release and run it - no Python install needed.

If you'd rather run from source, use the bundled venv to avoid Qt platform-plugin errors that come from accidentally launching with a different Python:

```
run.bat            (or)
.\run.ps1
```

`run.py` works too, as long as you've activated the project venv first.

## Changelog for v3.1.1

Patch release for big multi-INCLUDE decks where pyNastran's strict parse fails and the importer used to grind through the slow lenient card-by-card path with a frozen progress bar.

### Chunked-strict parser
- New `_read_bdf_chunked` splits the bulk-data section into ~20k-card chunks and runs pyNastran's `read_bdf(punch=True, validate=False)` on each chunk independently. Each chunk is one fast `read_bdf` call instead of N slow per-card `add_card` calls - typically 10-50x faster than whole-deck lenient on a real aerospace deck. Chunks that fail (e.g. one bad card type) fall back to lenient on just their ~20k lines so the rest of the deck still gets the fast path.
- When pyNastran's strict whole-deck parse fails, the importer now tries chunked-strict next, before falling through to whole-deck lenient as a last resort.
- Imported chunk BDFs are merged into a single master model by dict-level update of nodes, elements, properties, materials, coords, loads, SPCs, MPCs, DMIGs, sets, params, and ~20 other entity dicts. First definition wins on duplicate IDs (permissive).

### Faster + visible lenient path
- The per-card `sys.stdout` redirect in `_read_bdf_lenient` is gone - that swap was happening on every card and added meaningful overhead on million-card decks. Now stdout is redirected once at the top of the loop and restored once at the bottom.
- The lenient and chunked paths both emit per-card progress callbacks every ~500 cards. The import dialog's progress bar now advances smoothly from 26% through 98% during the fallback parse phases. Previously the bar sat at a static 50% the entire time lenient ran, looking hung.
- Status text shows current position: "Chunked parsing: 240,000 / 1,250,000 cards" or "Lenient parsing card 240,000 of 1,250,000".

## Changelog for v3.1.0

Minor feature release: full support for real-world multi-file aerospace Nastran decks. v3.1.0 supersedes the briefly-tagged v3.0.1 and v3.0.2 (their work is consolidated here under proper semver - features land in a minor bump, not a patch).

### INCLUDE-aware import
- Decks with INCLUDE statements (in executive, case, or bulk sections) now import reliably. Relative paths like `..\BULK\foo.bdf` resolve against the original deck's directory even when the importer falls back to temp-file strategies.
- New `_inline_includes` pre-expands every INCLUDE into a single flat buffer as a fallback path. If pyNastran's strict parse fails on the original, the importer tries again on the flattened deck before resorting to lenient card-by-card parsing.
- Lenient parsing now sees cards from included files too (previously the lenient path scanned only the top-level file and silently dropped everything that lived behind an INCLUDE).
- File-open / import dialogs now accept `.bdf`, `.dat`, `.nas`, `.pch`, and `.asm` (the latter two are common as INCLUDE targets in real decks).

### Multi-line INCLUDE statements
- New `_parse_include_statement` handles every continuation form seen in real decks: single-line, unclosed-quote spanning lines (the most common multi-line form, used to keep paths under MSC Nastran's 72-column field limit), continuation lines with leading `+` or `*` Nastran markers, and indented continuation.
- The previous regex required the close quote on the same line, silently dropping files referenced by 2-line INCLUDEs. The new character-by-character parser accumulates the path across lines until the matching close quote is found.

### Streaming parser with per-file progress
- New `_read_bdf_streaming` walks the INCLUDE chain first to enumerate every reachable file, then inlines them one at a time, then hands the flat buffer to pyNastran. The import dialog now shows a determinate progress bar (0..100%) that advances file-by-file during the inline pass.
- Per-file status: "Reading BULK\\foo.bdf (3 of 12)". The previous import dialog showed only "Parsing BDF..." with no indication of progress on multi-file decks.
- The progress bar is deliberately determinate (a static rectangle that fills), not the marching-ants spinner that caused the GIL-starvation bug in v3.0.0. The window may still freeze briefly during pyNastran's parse phase, but you can see how far the inline pass got before that.
- New `_gather_include_chain` walks every reachable INCLUDE without inlining, so the streaming parser can enumerate files (and total bytes) up front.

### Performance: no cross-reference on import
- `BDF.read_bdf` is now called with `xref=False` by default. Node Runner is a pre-processor and doesn't need pyNastran's cross-reference web at open time, where it can easily double or triple parse time on large decks. Callers that need it can run `model.cross_reference()` on demand.

### Tests
- 12 new pytest cases across `TestIncludeStatements`, `TestMultiLineInclude`, and `TestStreamingParser` cover INCLUDE resolution, every multi-line continuation form, end-to-end streaming, progress-callback emission, and the include-chain enumerator. 118 tests pass.

## Changelog for v3.0.0

A multi-phase upgrade across the export engine, UX layer, rendering pipeline, performance, mesh editing, loads/BCs, result interaction, and the cross-section UI.

### Engine: flexible export and format detection
- Choose Short (8-char), Long (16-char), or Free (comma-separated) field format per save via the new Export Options dialog.
- BDF field format is auto-detected on open (`detect_bdf_field_format`) and pre-fills the next save dialog.
- New Settings menu: Export Defaults..., Units... Both persist via QSettings.
- Persistent status bar shows `Model | Nodes | Export | Units` at all times.
- Group export respects the chosen field format.

### UX: command palette, threaded import, render debouncer
- Ctrl+P opens a fuzzy-search Command Palette across all menu actions (124 auto-harvested from the menu bar).
- File Open / Import run on a QThread with a modal progress dialog and Cancel button. UI stays responsive during heavy parses.
- 16ms render debouncer (`_request_render`) coalesces render bursts.

### Rendering: kill the green blob
- Node cloud renders through `vtkPointGaussianMapper` for crisp, anti-aliased circular splats instead of GL_POINTS rectangles.
- ScaleFactor sized to the model diagonal; cloud reads at any zoom level.
- Default node size 5 -> 3, default node color neon green -> theme accent (Catppuccin blue).
- Selection highlight uses Catppuccin peach instead of yellow.

### Performance: adaptive LOD, ghost mode, parallel parsing
- Adaptive node decimation: above 50,000 points the displayed cloud strides down to ~10,000 for snappy rendering. Toggleable in Settings; full picking-resolution preserved on the underlying mesh.
- Ghost mode (View > Ghost Hidden Groups): hidden groups render as translucent wireframe overlays so you keep spatial context.
- Experimental parallel BDF parsing (Settings > Parallel BDF Parsing): chunks the bulk-data section across worker threads and dict-merges the results. Falls back to single-threaded on small files, INCLUDE directives, or any chunked-parse failure. Default OFF.

### Mesh editing (Theme A)
- Split CQUAD4 into 2 CTRIA3, refine quads/trias 1-into-4 with shared midside nodes, combine two trias into a quad.
- Laplacian smoothing of node sets, mirror elements with reversed connectivity, copy elements with translate/rotate.
- Insert node on edge with neighbor splitting; element quality metrics extended to CHEXA / CTETRA / CPENTA.

### Loads and BCs (Theme C)
- Creation dialogs for PLOAD1, PLOAD2, SPCD, MPC, RBAR, RBE1, RSPLINE, BOLT.
- SPCD wires through to the analysis-set case-control LOAD entry.

### Result interaction (Theme B)
- Probe tool: hover for node/element values; floating tooltip with displacement and stress components.
- Result Browser dock with sortable nodal/element tables, expression bar (sandboxed numpy eval), 3D and table cursor sync.
- Vector overlays (displacement, principal stress, reaction force) and strain contour from OP2.
- Animation timeline with play / pause / scrub / speed; export-GIF retained alongside.

### Cross-Section dock + Define Plane methods
- Replaced the always-on-top Cross Section dialog with a dockable widget that tabifies alongside Results / Animation / Vectors.
- New Define Plane picker with five methods: Global Axis, Coordinate Value, 3 Points, Normal + Point, Two Nodes (line normal).
- Translucent plane outline actor for visual feedback at the cut. Plane definitions are serializable so the dock can summarize and persist the last-used method.

### Big-model robustness
- Switched the BDF importer back to a synchronous parse with a busy cursor and lightweight status dialog. The previous QProgressDialog spinner starved the worker thread of CPU on Windows; pyNastran's bulk parse is mostly Python and the GIL serialized worker progress with the spinner, so the UI hung on 60k+ files. Sync parse is briefly unresponsive (~6 s for 60k nodes) but always finishes.
- 1M-node b747 streaming generator (`examples/generate_b747.py`) writes BDFs without ever building a giant pyNastran object in memory.

### Launcher
- `run.bat` and `run.ps1` pin the project venv Python so PySide6's Qt platform plugins resolve correctly. Avoids the "no Qt platform plugin could be initialized" dialog you get when `python run.py` picks up the wrong interpreter.

## Changelog for v2.0.0

This is a major release with coordinate system protection, enhanced display options, improved view controls, geometry transforms, and numerous quality-of-life improvements.

### ✨ New Features

- **Coordinate System Protection:**
  - Global Rectangular (CID 0), Cylindrical (CID 1), and Spherical (CID 2) coordinate systems are now immutable - they cannot be deleted or overwritten by any tool, import, or generator.
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
