# Node Runner

> A lightweight, open-source pre, solve, and post environment for Nastran FEA models.

![python](https://img.shields.io/badge/python-3.9%2B-blue) ![license](https://img.shields.io/badge/license-Apache--2.0-green) ![platform](https://img.shields.io/badge/platform-Windows-lightgrey) ![tests](https://img.shields.io/badge/tests-321%20passing-brightgreen) ![version](https://img.shields.io/badge/version-5.3.0-informational)

Node Runner reads a Nastran bulk-data deck, renders the mesh interactively, drives the open-source MYSTRAN solver as a subprocess, and loads the results back into a Femap-style post-processing toolbox. It is built for engineers who want a real FEA workflow without paying a commercial-license-per-seat-per-year.

---

## Why Node Runner

Commercial FEA pre / post-processors are powerful but cost-locked. Open-source alternatives (PrePoMax, FreeCAD FEM, Gmsh) target their own solvers and don't play cleanly with Nastran-dialect decks. Node Runner sits in the gap:

| | Commercial (Femap / Patran / HyperMesh) | Node Runner | Other open-source |
|---|---|---|---|
| Nastran BDF / DAT import (with INCLUDE chains) | yes | yes | partial |
| Open-source solver bundled | no | MYSTRAN (linked) | varies |
| OP2 + F06 results parsing | yes | yes | rare |
| Group-based scoping (Femap-style) | yes | yes | no |
| Interactive 3D viewport (pan / zoom / pick) | yes | yes (PyVista) | yes |
| Deform + contour independently | yes | yes | partial |
| Per-element / per-property colouring | yes | yes | no |
| Tabular results browser | yes | yes | partial |
| License | proprietary | Apache-2.0 | varies |
| Per-seat annual cost | thousands | zero | zero |

The trade-off is breadth: commercial tools cover dozens of solvers, hundreds of element formulations, decades of edge cases. Node Runner targets a focused slice (linear static, normal modes, linear buckling) and does that slice well.

---

## Quick start

```powershell
# 1. Grab the latest release.
#    https://github.com/whoaone/node-runner/releases/latest
#    Download "Node Runner.exe" (~242 MB, Windows).
#
# 2. Run it. No Python or dependency install required.
#
# 3. File > Open > examples/cantilever_beam.bdf
#
# 4. (Optional, for solving) Edit > Preferences > MYSTRAN tab,
#    point at the MYSTRAN executable from https://www.mystran.com/.
#    Click Detect to auto-find it on PATH.
#
# 5. Analysis > Run Analysis (MYSTRAN), or Ctrl+R.
#    The post-processing toolbox auto-opens on success.
```

That should take under 60 seconds end-to-end once the exe is on disk.

---

## Feature highlights

### Pre-processing

* **BDF / DAT import** with full INCLUDE-chain resolution and lenient-mode fallback for decks with non-standard cards.
* **6 themed Model menus** (Geometry / Elements / Sections / Loads / Constraints / Modify) so you find what you need without scrolling.
* **Interactive 3D viewport** built on PyVista. Pan / zoom / pick / clipping plane / shading toggle / theme switch all live.
* **Per-property and per-material colouring**, with per-cell RGB rendering on every shell and beam.
* **Femap-style Groups** that house nodes, elements, properties, materials, AND coordinate systems. Auto-group by Property, Material, or Element Type. Group-scoped exports and runs.
* **AnalysisSet** as the unit of work: SOL type, subcases, EIGRL parameters, output requests, params, scope (groups / load SIDs / SPC SIDs). One set, one run.
* **Pre-flight scanner** that walks the deck before launching the solver, flagging cards MYSTRAN cannot handle (aero, nonlinear, contact, optimization) with a clean modal report.

### Solving (MYSTRAN integration)

* **SOL 1 / 101** linear static
* **SOL 3 / 103** normal modes (real eigenvalue)
* **SOL 5 / 105** linear buckling
* **Non-blocking QProcess driver** with live stage tracking (Reading BDF, Assembling K, Decomposing, Solving, Writing output).
* **Auto-load results** on success: OP2 first (via pyNastran), F06 fallback (in-house parser) for decks where pyNastran's MYSTRAN dialect coverage hasn't caught up.
* **Run history** under `~/.node_runner/mystran_runs/runs/` with stamped folder names and a `run_meta.json` summary. Browse previous runs from the Analysis menu.
* **Case-control resolver** that bridges per-AnalysisSet Output Requests into per-subcase CASE CONTROL lines (STRESS / STRAIN / FORCE / DISPLACEMENT).

### Post-processing toolbox

* **Three collapsible sections** (Results / Deform / Contour) you toggle independently. Want a stress contour on a non-deformed mesh? Easy. Animate a displacement vector without contouring? Also easy.
* **Output Vector dropdown** auto-detects what's in the bundle: Displacement, SPC Force, Stress (CQUAD4, CTRIA3, CBAR fiber-level), Strain (OP2 only), Element Forces (CBAR moments / shears / axial / torque), Eigenvector modes.
* **Vector Info dialog** (Femap-style) showing Min / Max for the active vector with their Node or Element IDs and a copy-to-clipboard button.
* **Scalar bar with "(geom x N)" tag** in `% of model` scale mode so the engineering value on the bar and the visible geometry never feel disconnected.
* **Discrete bands**, levels, palette selection (jet / viridis / plasma / coolwarm / grayscale), and auto vs manual range gating.
* **Min / Max marker callouts** and per-element value labels, both opt-in.
* **Animation** in two modes: Sine wave (a single deformed mode oscillating) and Through Modes (steps through every mode in a normal-modes bundle). Configurable frame count and delay.

### Tools and tabular access

* **Data Table** under Tools, a dockable dialog showing nodal or element data for whatever entities you have selected, with one column per output vector. Live-syncs with selection. Exports to clipboard, CSV, or XLSX (XLSX needs `openpyxl`).
* **Find Entities** with a Ctrl+F shortcut, a rich entity-type picker (Node / Element / Property / Material / Coord), and a read-only details dialog on single-id selection with an Edit... button that routes to the existing editor.
* **Probe Mode** (Ctrl+Shift+P) for hover-tooltip readouts on the 3D viewport.

---

## Card and PARAM compatibility

The translator at `node_runner/solve/mystran_export.py` is the single source of truth for what Node Runner passes through, rewrites, or drops when exporting to MYSTRAN.

| Category | Treatment |
|---|---|
| Standard linear element cards (CQUAD4, CTRIA3, CBAR, CBEAM, CROD, CBUSH, CHEXA, CTETRA, CPENTA, CONROD, RBE2, RBE3, RBAR) | Pass through |
| Material cards (MAT1, MAT2, MAT8, MAT9) | Pass through |
| Property cards (PSHELL, PCOMP, PBAR, PBEAM, PBARL, PROD, PSOLID, PBUSH) | Pass through |
| `PARAM,WTMASS`, `PARAM,GRDPNT`, `PARAM,EPSIL`, `PARAM,POST` | Pass through |
| MSC-only PARAMs (`COUPMASS`, `AUTOSPC`, `BAILOUT`, `MAXRATIO`, `K6ROT`) | Dropped with a warning |
| Aero (`CAERO*`, `AESTAT`, `TRIM`, `FLUTTER`, etc.) | Blocking, pre-flight refuses to run |
| Nonlinear (`NLPARM`, `NLPCI`) | Blocking |
| Contact (`BCTABLE`, `BCONTACT`, `BCSET`) | Blocking |
| Optimization (`DESVAR`, `DCONSTR`, `DRESP*`, `DOPTPRM`) | Blocking |
| Unsupported elements (`CQUADR`, `CTRIAR`, `CWELD`, `CSEAM`) | Blocking |
| `PARAM,SOLLIB` (default `SPARSE`) | Injected if not present |
| `PARAM,QUAD4TYP` (default `MIN4T`) | Injected if not present |

The pre-flight dialog explains every dropped or blocked card before the solver fires.

---

## Under the hood

```
┌──────────────────────┐
│   Node Runner.exe    │
│   (PyInstaller, ~242 │
│    MB single-file)   │
└──────────┬───────────┘
           │
           ├── PySide6           ← Qt UI (menus, dialogs, sidebar tabs)
           ├── PyVista           ← 3D viewport + actor management
           ├── VTK               ← low-level rendering
           ├── numpy             ← per-cell scalar arrays
           ├── pyNastran         ← BDF parse / write + OP2 results
           ├── openpyxl          ← Data Table xlsx export
           └── mystran.exe       ← driven as a QProcess subprocess
                                   (not bundled; user-configured path)
```

The codebase is organized into a few major modules:

* `node_runner/mainwindow.py` (~20k lines) houses the `MainWindow` class, the rendering pipeline, every menu wiring, and all the cross-cutting handlers.
* `node_runner/dialogs/` is one file per dialog (creators, editors, group editor, results tab, preferences, etc.).
* `node_runner/solve/` carries the MYSTRAN export translator, the QProcess runner, the F06 fallback parser, and the run-history helpers.
* `node_runner/commands.py` provides the Command pattern (undo / redo) for every model mutation.
* `node_runner/scene_build.py` is the vectorized cell-array builder for the 3D viewport.
* `tests/` carries 321 pytest cases, run headless via pytest-qt.

---

## Building from source

```powershell
# Clone
git clone https://github.com/whoaone/node-runner.git
cd node-runner

# Install dependencies into a venv
python -m venv venv
.\venv\Scripts\python.exe -m pip install -r requirements.txt

# Run from source
.\venv\Scripts\python.exe -m node_runner

# Or run the tests
.\venv\Scripts\python.exe -m pytest tests/ -q

# Build the standalone exe
.\venv\Scripts\pyinstaller.exe --noconfirm "Node Runner.spec"
# Output lands in dist/Node Runner.exe
```

Python 3.9+ required. The PyInstaller build needs a working Qt runtime in the venv (pulled in automatically by `pip install pyside6`).

---

## Tests

```powershell
.\venv\Scripts\python.exe -m pytest tests/ -q
```

321 pytest cases pass on a clean tree. One end-to-end MYSTRAN test is gated behind `pytest -m mystran_e2e` and skipped by default. Headless rendering uses pytest-qt's offscreen Qt platform.

---

## Roadmap

Plans for the next several releases (subject to change):

### v5.3.1 (next)
* CTRIA3 + CBEAM + CROD F06 stress parsers.
* F06 strain block expansion once MYSTRAN starts emitting strains.
* Broader OP2 parser coverage (beyond CQUAD4 + CTRIA3 + CBAR).

### v5.4
* Element-local coordinate-system contours.
* Free-body diagrams on demand for any node set (Femap section 5.11.1.1).
* Element quality histogram with jacobian / skew / aspect-ratio thresholds.
* CSV results export with full statistics per output vector.

### v6.0
* External solver hooks (NX Nastran, MSC Nastran) via subprocess driver. Same workflow as MYSTRAN.
* Cloud / remote solver execution.
* Composite ply-by-ply post-processing UI.
* PDF report generation from saved results.

---

## Contributing

Issues and pull requests welcome. Before opening a PR:

1. Run `pytest tests/ -q`. The full suite must pass.
2. Add a regression test for any bug fix or new feature.
3. Match the existing code style (4-space indent, no trailing whitespace, no em-dashes).
4. Update the changelog block at the top of this README in the section for the next release.

The codebase favors instrumentation over guesswork. If you touch a hot path, leave a `perf_event` catcher behind so the next developer can see what changed.

---

## Acknowledgements

Built on the shoulders of these open-source projects:

* [MYSTRAN](https://www.mystran.com/) by Bill Case and contributors. The actual solver.
* [pyNastran](https://github.com/SteveDoyle2/pyNastran) by Steve Doyle. BDF parse / write + OP2 reader.
* [PyVista](https://pyvista.org/) and the underlying [VTK](https://vtk.org/). The 3D viewport.
* [PySide6](https://doc.qt.io/qtforpython-6/) (Qt for Python). The UI.
* [openpyxl](https://openpyxl.readthedocs.io/) for the Data Table xlsx exporter.

---

## License

Apache License 2.0. See `LICENSE`.

---

## Changelog

### v5.3.0 (2026-05-17)

UX consolidation + entity editor + group polish + results-reading expansion.

**Major UX consolidation:**

* **Model menu regrouped** into themed submenus (Geometry / Elements / Sections / Loads / Constraints / Modify). The 19-item flat list is gone; mesh-editing operations now live inside the Modify submenu. Duplicate `Orient Normals` removed from Tools.
* **Settings menu removed.** Export Defaults, Units, Adaptive Node LOD, Element LOD Threshold, and the retired Parallel BDF Parsing toggle moved into `Edit > Preferences` as new General / Export / Performance tabs. One place for all configuration.
* **Solution Type submenu removed from Analysis menu.** SOL type lives only on the active AnalysisSet now (set via `Analysis > Analysis Set Manager`), removing a second source of truth.
* **Edit > Select Entities** workflow reworked. The old "what type?" `QInputDialog` is replaced with a richer picker that offers Node / Element / Property / Material / Coord. Selecting a single entity opens a read-only details dialog with an `Edit...` button that routes to the existing editor for that type.
* **Probe Mode** gains a tooltip, status-bar help, and the `Ctrl+Shift+P` shortcut.
* **Group editor expanded.** The Edit Group dialog gains Materials, Properties, and Coords tabs to match the Femap-faithful group schema. Stale node/element IDs now show `(not in current model)` instead of leaving cells blank, and group IDs are coerced to int so JSON-stringified group payloads resolve cleanly. `CreateGroupCommand` now produces the full 6-field schema; `_group_by_material` and `_group_by_element_type` auto-collect dependencies (properties / materials / coord systems).

**Bug fixes:**

* **RBE2 / RBE3 spider follows the deformed mesh.** Previously the rigid-element lines stayed pinned to undeformed node positions while the slave grids moved with the plate. The actor is now rebuilt with the deformed point coordinates on each results render.
* **Group editor populates every row.** The 5 tables in the Edit Group dialog had `setSortingEnabled(True)` set before the populate loop, causing Qt to scramble row / column correspondence mid-insert. Sorting is now disabled during populate and re-enabled after, matching the canonical pattern used by `NodeElementInfoDialog`.
* **View cube duplicate-actor bug.** Toggling between light and dark mode used to accumulate orphaned orientation-marker widgets, producing a cluster of overlapping cubes. `_toggle_theme` now tears down the old widget and sweeps stale axes actors before adding the fresh marker.
* **Default deform scale = `Actual` at `1.0`.** Previously `% of model` at `10.0` auto-fit the deformation to ~10% of the bbox, hiding the engineering magnitude. The new default shows the raw values; switch to `% of model` for auto-scaling on small-displacement decks.

**Results-reading expansion:**

* **CQUAD4 + CTRIA3 + CBAR F06 stress parsers** (top/bottom fiber, max-of-fiber von Mises for shells; SA1-SA4 + axial + max/min for bars).
* **CBAR engineering forces parser** (bending moments at A and B, shear, axial, torque).
* **SPC Force vector parsing** that populates the bundle. Previously the section was detected but rows were dropped.
* **Bundle schema gains `strains` + `forces` slots.** F06 path matches the OP2 path so the Strain / Force dropdown categories work uniformly.
* **Output Vector dropdown gains Strain and Force categories.** Auto-detected from the bundle slots; categories only appear when their dict is non-empty.
* **STRESS / STRAIN / FORCE Output Requests reach the deck.** `_apply_case_control_from_set` now falls back to the analysis-set-level Output Requests when a per-subcase entry isn't explicitly set, so MYSTRAN actually computes what the user requested.
* **Scalar-bar `(geom x N)` annotation** in `% of model` scale mode. Makes the auto-fit multiplier explicit so users don't read 43.45 on the bar and wonder why the visible geometry is 1.4 inches.
* **Vector Info dialog.** New `View > Vector Info...` opens a modeless readout of the active Output Vector's min/max with their Node or Element IDs. Refreshes live as you flip Output Set / Output Vector. Copy button writes a tab-separated row to the clipboard.

**Still deferred to v5.3.1+:** CBEAM / CROD F06 stress (no test fixtures yet), F06 strain (MYSTRAN-side limitation), broader OP2 parser expansion.

---

Older releases: see the [Releases page](https://github.com/whoaone/node-runner/releases) for v1.0.0 through v5.2.0 changelogs.
