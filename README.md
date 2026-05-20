# Node Runner v5.3

A lightweight pre, solve, and post environment for Nastran models. Built with Python, PySide6, PyVista, and the open-source MYSTRAN solver.

## Changelog for v5.3.0

UX consolidation + entity editor + group polish.

**v5.3.0 same-day patch (2026-05-17):**

* **RBE2 / RBE3 spider follows the deformed mesh.** Previously the rigid-element lines stayed pinned to undeformed node positions while the slave grids moved with the plate. The actor is now rebuilt with the deformed point coordinates on each results render. New `results.rbe_redrawn` profile-log catcher.
* **Group editor populates every row.** The 5 tables in the Edit Group dialog had `setSortingEnabled(True)` set before the populate loop, causing Qt to scramble row / column correspondence mid-insert. Sorting is now disabled during populate and re-enabled after — matching the canonical pattern used by `NodeElementInfoDialog`.
* **Default deform scale = `Actual` at `1.0`.** Previously `% of model` at `10.0` auto-fit the deformation to ~10% of the bbox, hiding the engineering magnitude. The new default shows the raw values; switch to `% of model` for auto-scaling on small-displacement decks.
* **Output Vector dropdown shows CQUAD4 stress.** The F06 fallback parser now extracts top/bottom-fiber stress from CQUAD4 elements (Normal-X / Normal-Y / Shear-XY / Max / Min Principal / von Mises). CTRIA3 / CBEAM stress and strain / force parsing remain deferred to v5.3.1.
* **STRESS / STRAIN / FORCE Output Requests reach the deck.** `_apply_case_control_from_set` now falls back to the analysis-set-level Output Requests when a per-subcase entry isn't explicitly set, so MYSTRAN actually computes what the user requested.

**v5.3.0 patch 2 (2026-05-17, same day):**

Results-reading expansion. Closes the F06-fallback gap that left users seeing only Displacement in the Output Vector dropdown:

* **F06 CTRIA3 stress parser.** Mirrors the v5.3.0 CQUAD4 parser; handles the "Anywhere"/"in elem" row markers MYSTRAN emits. Combined with the existing CQUAD4 parser this covers shell stress for mixed-element decks.
* **F06 CBAR stress + CBAR engineering forces parsers.** Two-station fiber stresses (SA1-SA4 + SB1-SB4 + axial + max/min) plus the one-row force table (M1/M2 at each end + V1/V2 + axial + torque). Beam-element decks now show stress and force contours.
* **F06 SPC Force vector parsing.** The section was previously detected but its rows weren't stored; the `spc_forces` slot now populates. SPC Force - T1/T2/T3/Magnitude appear in the Output Vector dropdown.
* **Bundle schema gains `strains` + `forces` slots.** F06 path matches the OP2 path so the dropdown's Strain/Force category gates work uniformly. Strain stays empty on F06-only paths (MYSTRAN doesn't emit strain to F06; OP2 fills it when pyNastran's MYSTRAN dialect is available).
* **Output Vector dropdown gains Strain and Force categories.** Auto-detected from the bundle slots; categories only appear when their dict is non-empty.
* **Scalar-bar `(geom × N)` annotation in `% of model` scale mode.** Makes the auto-fit multiplier explicit so users don't read 43.45 on the bar and wonder why the visible geometry is 1.4 inches. No annotation in `Actual` mode (the v5.3.0 default).
* **Femap-style Vector Info dialog.** New `View → Vector Info...` opens a modeless readout of the active Output Vector's min/max with their Node or Element IDs. Refreshes live as you flip Output Set / Output Vector. Copy button writes a tab-separated row to the clipboard.

Still deferred to v5.3.1+: CBEAM/CROD stress (no fixtures yet), broader OP2 parser expansion, strain F06 (MYSTRAN-side limitation).

* **Model menu regrouped** into themed submenus (Geometry / Elements / Sections / Loads / Constraints / Modify). The 19-item flat list is gone; mesh-editing operations now live inside the Modify submenu. Duplicate `Orient Normals` removed from Tools.
* **Settings menu removed.** Export Defaults, Units, Adaptive Node LOD, Element LOD Threshold, and the retired Parallel BDF Parsing toggle moved into `Edit -> Preferences` as new General / Export / Performance tabs. One place for all configuration.
* **Solution Type submenu removed from Analysis menu.** SOL type lives only on the active AnalysisSet now (set via `Analysis -> Analysis Set Manager`), removing a second source of truth.
* **Group editor expanded.** The Edit Group dialog gains Materials, Properties, and Coords tabs to match the Femap-faithful group schema. Stale node/element IDs now show `(not in current model)` instead of leaving cells blank, and group IDs are coerced to int so JSON-stringified group payloads resolve cleanly. `CreateGroupCommand` now produces the full 6-field schema; `_group_by_material` and `_group_by_element_type` auto-collect dependencies (properties / materials / coord systems).
* **Probe Mode** gains a tooltip, status-bar help, and the `Ctrl+Shift+P` shortcut.
* **Edit -> Select Entities** workflow reworked. The old "what type?" `QInputDialog` is replaced with a richer picker that offers Node / Element / Property / Material / Coord. Selecting a single entity opens a read-only details dialog with an `Edit...` button that routes to the existing editor for that type.
* **View cube duplicate-actor bug fixed.** Toggling between light and dark mode used to accumulate orphaned orientation-marker widgets, producing the cluster of overlapping cubes the user reported. `_toggle_theme` now tears down the old widget and sweeps stale axes actors before adding the fresh marker. New `view_cube.rebuilt` profile-log event tracks the actor count.

## What it does

Node Runner reads a Nastran bulk-data file (`.bdf` or `.dat`), renders the mesh interactively, and lets you inspect or edit it before running an analysis. After a solve, it loads the results (OP2 or F06), animates deformations, and paints contour gradients on the elements.

Core capabilities:

* Imports a BDF / DAT deck (with INCLUDE chains) and renders nodes, elements, properties, materials, constraints, and loads in a 3D viewport.
* Pre-flight scanner inspects a deck for cards MYSTRAN does not support before launching the solver.
* Drives MYSTRAN as a subprocess with live stage feedback, then auto-loads the resulting OP2 or F06 back into a Results Toolbox.
* Post-Processing Toolbox with three collapsible sections (Results, Deform, Contour) for picking an output vector and toggling deformation or contour rendering independently.
* Data Table dialog (under Tools) for tabular inspection of nodal or element results.
* Group-based scoping for runs, exports, and the BDF "Save BDF" path that round-trips Node Runner metadata via `$ NR-META` comment lines.

## Running

The easiest path on Windows is to grab `Node Runner.exe` from the latest release and run it. No Python install is required.

To run from source, use the bundled venv:

```
run.bat
```

(or `.\run.ps1` from PowerShell, or `run.py` after activating the project venv).

## Running analyses with MYSTRAN

Node Runner supports the solution types MYSTRAN runs:

* SOL 1 / 101 (linear static)
* SOL 3 / 103 (normal modes)
* SOL 5 / 105 (linear buckling)

Quickstart:

1. Download a MYSTRAN binary from the [upstream project](https://www.mystran.com/). Node Runner does **not** bundle it.
2. Open **Edit → Preferences → MYSTRAN** and set the executable path (or click *Detect* to search `PATH` and common install locations).
3. As a first-run sanity check, import [`examples/cantilever_beam.bdf`](examples/cantilever_beam.bdf) and run **Analysis → Run Analysis (MYSTRAN)** (Ctrl+R). The tip displacement at node 11 should land around `T3 = -0.0407 in`, matching the Euler-Bernoulli analytical solution to within one percent.
4. For your own work, load a BDF, then **Analysis → Run Analysis (MYSTRAN)**.
5. The pre-flight scan runs first. It flags any cards MYSTRAN cannot handle (aero, nonlinear, contact, optimization, etc.). Resolve blocking issues, accept warnings.
6. The solver runs in a non-blocking dialog with stage tracking pulled from MYSTRAN's stdout.
7. On finish, results auto-load into the Post-Processing Toolbox. The status bar updates with the SOL number and subcase count.

## Card and PARAM compatibility

The translator at `node_runner/solve/mystran_export.py` is the single source of truth for what Node Runner passes through, rewrites, or drops when exporting to MYSTRAN. A summary:

* **Pass through:** `PARAM,WTMASS`, `PARAM,GRDPNT`, `PARAM,EPSIL`, `PCOMP`, `MAT8`, all linear element cards.
* **Dropped with a warning:** MSC-only PARAMs (`COUPMASS`, `AUTOSPC`, `BAILOUT`, `MAXRATIO`, etc.).
* **Blocking (pre-flight refuses to run):** aero cards (`CAERO*`, `AESTAT`, `TRIM`, `FLUTTER`), nonlinear (`NLPARM`), contact (`BCTABLE`, `BCONTACT`), optimization (`DESVAR`, `DCONSTR`, `DRESP*`), and unsupported element types (`CQUADR`, `CTRIAR`, `CWELD`, `CSEAM`).
* **Injected on write:** `PARAM,SOLLIB` (default `SPARSE`), `PARAM,QUAD4TYP` (default `MIN4T`), plus any user-configured overrides from Preferences.

## Post-Processing Toolbox

The Results tab in the left sidebar opens after a solve (or after `File → Load Results`). It has three collapsible sections:

* **Results** holds two dropdowns: the **Output Set** (a subcase from the bundle) and the **Output Vector** (displacement components, stress components, SPC forces, eigenvector modes). The Output Vector you pick here drives the contour gradient and is the canonical "what you are looking at" selection.
* **Deform** controls whether the model deforms. Style is Undeformed (default), Deformed, or Animate. Scale mode is `% of model` (auto-fits the deflection to a percentage of the bounding box diagonal) or `Actual`. Animation options (Sine wave or Through Modes, frame count, delay) are configurable in advance.
* **Contour** controls whether a gradient is painted. Style is No Contour (default), Filled, Filled + Edges, or Discrete Bands. Data Conversion picks one of Average (nodal), No Averaging (centroid), Max at Node, or Min at Node. Levels and palette control the colormap. Auto Range vs Manual gates the min/max spin boxes. Annotations add Min/Max markers and per-element value labels.

Deform and Contour are independent. You can deform a stress contour, deform without a contour, contour without deforming, or animate either or both.

## Data Table (Tools menu)

A dockable dialog that shows tabular nodal or element data for whatever entities you have selected in the viewport, with one column per output vector. Selection updates the table live. The current view exports to clipboard, CSV, or XLSX (requires `openpyxl` for the latter).

## Tests

```
pytest tests/ -q
```

Roughly 285 tests pass on a clean tree. One end-to-end MYSTRAN test is gated behind `pytest -m mystran_e2e` and skipped by default.

## License

See `LICENSE`.
