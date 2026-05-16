# Node Runner v5.2

A lightweight pre, solve, and post environment for Nastran models. Built with Python, PySide6, PyVista, and the open-source MYSTRAN solver.

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
