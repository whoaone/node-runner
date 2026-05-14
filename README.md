# Node Runner v4.1.0

A lightweight pre-processor for creating, editing, and visualizing Nastran models. Built with Python, PySide6, and PyVista.

## Running

The simplest path on Windows is to grab `Node Runner.exe` from the latest release and run it — no Python install needed.

If you'd rather run from source, use the bundled venv to avoid Qt platform-plugin errors that come from accidentally launching with a different Python:

```
run.bat            (or)
.\run.ps1
```

`run.py` works too, as long as you've activated the project venv first.

## Changelog for v4.1.0

Consolidated release of the v4.0.x perf push. Builds on v3.5.0 with
a complete rework of the visibility hot path, the moment-glyph
convention, the camera-reset behaviour, the property-coloring
pipeline, and a per-Nastran-type actor architecture that takes
category toggles on the MEGA `MEGA-XWFF_04A_DUL_FUC-Entry.dat` deck
(2.4 M elements / 1,962 load SIDs / 25 INCLUDEs) from ~5 s per
click to sub-100 ms.

### Premium per-click feel on MEGA

| User action | Pre-v4.1.0 | v4.1.0 |
|---|---|---|
| Toggle a category (Plates, Beams, Solids, Bars, Rods, Bushes, Shear, Gap) | ~5 s | **sub-100 ms** |
| Toggle Rigid / Masses / Plotels | full visibility cycle | sub-100 ms |
| Toggle coord system | full visibility cycle | sub-100 ms |
| Toggle a single load_set / constraint_set | full visibility cycle | sub-100 ms (or fall-through ≤ 1 s) |
| Click any tree item with the visibility-load loop in the path | up to 3.7 s (1,962 SIDs walked) | ~40 ms (bounded to visible SIDs only) |

### Visibility correctness

- **Bushes category** (CBUSH) now masks correctly. Previously the
  Bushes checkbox ran the full cycle but never modified the
  visibility mask, so CBUSH elements stayed visible regardless.
- **Camera stays put** on tree toggles. v3.5.0 would
  snap the view to fit the whole model on every glyph add or
  legacy mesh rebuild. Pre-built actors and every glyph add path
  now pass `reset_camera=False`; the camera only resets on
  explicit view buttons (Iso, Top, Front, etc.) and the first
  model load.
- **Nodes are independent of category toggles.** The fast-path
  visibility router used to hide `nodes_actor` alongside legacy
  mesh leftovers; it now leaves the user-controlled Nodes display
  alone.

### Visual fidelity

- **Shading toggle** produces a visible lit / flat delta again.
  `plotter.clear()` inside `_update_viewer` was wiping the 4-light
  rig installed once in `__init__`. The rig now re-installs on
  every viewer rebuild.
- **Moment vectors** use the physics-textbook double-tipped arrow
  (▶▶) instead of the symmetric ←→ double-arrow. The shaft
  direction is the right-hand-rule thumb; the double tip
  distinguishes a moment vector from a force vector (single ▶).
- **PID coloring is value-keyed.** v4.1.0 writes per-cell RGB
  directly to `cell_data['cell_rgb']` and renders via
  `add_mesh(scalars='cell_rgb', rgb=True)`. The same PID renders
  the same color regardless of which actor / which rebuild / which
  mode-transition. Previously PyVista's `categories=True` +
  list-cmap pattern produced different colors for the same PID
  across actors with different unique-PID sets.

### Architecture

- **Pre-built per-Nastran-type actors.** At viewer-build time, the
  full `current_grid` is partitioned by `cell_data['type']` and
  one named PyVista actor is added per Nastran type (`cat_CQUAD4`,
  `cat_CBEAM`, `cat_CBUSH`, `cat_CHEXA`, etc.). Category toggles
  become `actor.SetVisibility(bool)` — sub-100 ms on MEGA.
- **Fast-path router** in `_handle_tree_item_changed` dispatches
  load_set / constraint_set / Rigid / Masses / Plotels / coord
  toggles to targeted handlers that bypass the full
  `_update_plot_visibility` cycle.
- **Filter-mode coexistence.** When PID / MID / isolate / hidden-
  groups state is active, the category fast path returns False
  and the call falls through to the legacy `_rebuild_plot` path
  for cell-level filtering. Legacy and pre-built actors never
  double-render (each path hides the other side's actors).
- **Color-mode rebuild.** Switching `View > Color by >
  Property ID | Element Type | Quality | Results` invalidates and
  rebuilds the pre-built actors so the new mode's coloring is
  applied. Quality / Results modes fall through to legacy
  (their scalars depend on external data flow).
- **Module-level constants.** `_NTYPE_TO_TYPE_CATEGORY`,
  `_NTYPE_TO_SHAPE_CATEGORY`,
  `_PRE_BUILT_SUPPORTED_COLOR_MODES` are the single source of
  truth for category-to-Nastran-type mapping and which color modes
  the pre-built path handles.

### Perf optimizations under the hood

- **`shells_compute_normals` cache.** The phong-shading
  `compute_normals` call on 2.36 M LOD'd shell cells costs ~1.1 s
  per cycle. Cache key:
  `(id(current_grid), n_cells, EID-hash, elem_shrink)`. Saves
  ~1.1 s on every cache hit. Invalidates automatically on
  `_update_viewer` (new `current_grid`) or PID/Type filter
  changes (different EID set).
- **Bounded SID-loop iteration.** A maintained
  `self._visible_load_sids` set replaces walking
  `model.loads.keys()` on every visibility cycle. On MEGA the
  load_sid_loop drops from 3.7 s to ~40 ms (97× faster).
- **Per-SID O(1) fast path.** `load_set` and `constraint_set`
  toggles build / remove glyphs only for the affected SID — no
  full visibility cycle.
- **Vectorized RBE actor build.** `_create_rbe_actors` runs in
  ~400 ms on MEGA's 46,782 rigid elements (was ~144 s of the
  legacy 434 s viewer build via Python-level cell iteration).
- **Parallel BDF parse.** Decks with 4+ INCLUDEs and ≥ 50 MB use
  a `ProcessPoolExecutor` per-INCLUDE parser
  (`node_runner/parallel_parse.py`); falls back to serial on
  smaller decks or any DMIG presence.
- **Status-bar element count** alongside the node count
  (`Elements: 2,452,222` next to `Nodes: 2,484,333` on MEGA).

### Diagnostics

`NR_PROFILE=1` produces a structured profile log at
`~/.node_runner/profile_<timestamp>.log`. ~30 named catchers
across visibility / actor / glyph / shading paths emit timing and
state. Examples:

```
[perf] viewer.apply_display_lod        wall=0.377s n_input_cells=2405441
[perf] viewer.lod_result               active=True displayed=248261 total=2405441
[perf] actors.build_per_category       wall=3.4s n_categories=7
[perf] actors.category_added           ntype=CQUAD4 n_cells=2347347 is_shell=True mode=property
[perf] fast_path.category_toggle       group_name=Plates n_actors_changed=2 wall<100ms
[perf] load_sid_loop.iter_bounded      n_iter=1 n_total=1962
[perf] shells_normals.cache_hit        n_cells=235765
[perf] glyph.moment_double_tip_built   n_points=68 n_cells=50
[perf] mem.rss_mb                      where=viewer.total mb=...
[perf] lights.rig_installed            n_lights=4
```

The `~30` catchers are intended to stay in forever — they catch
regressions silently on every subsequent build.

### Tests

The 202-test pytest suite continues to pass. New v4.1.0 tests
cover the parallel parse path and the tree snapshot / restore
round-trip.

---

## Changelog for v3.5.0

Feature release covering four items reported on
`MEGA-XWFF_04A_DUL_FUC-Entry.dat` (512,203 load entries, 3,911 load
combinations).

### Entity selection respects tree visibility (item 1)
Hiding a category in the Model tree (e.g. unchecking "Plates") now
filters them out of the selection pool automatically. The selection
bar gets a new `[ ] Include hidden` checkbox so the user can opt back
into the full pool when needed. Internally the bar stores the
*unfiltered* range on each entry and intersects against the current
view pool at `get_selected_ids()` time, so toggling Include-hidden
expands the result set without losing the original authored intent.

### Loads tab populates instantly even on huge decks (item 4)
`_populate_loads_tab` previously created one `QTreeWidgetItem` per
load entry at populate time - on a deck with 512k entries that was a
5-30 s hang. v3.5.0 makes Load Set and Constraint Set children
**lazy**: the headers get a placeholder child like `(N entries -
expand to load)`, and the real entries are built on first expand via
`_on_load_tab_item_expanded` (batched via `addChildren`). Tree
populate is now `<100ms` regardless of total entry count.

### Femap-style Load Combination editing (item 2)
New `LoadCombinationDialog` (modes: create / edit / copy) with:
- Combo SID + Title fields.
- Overall scale factor.
- (Scale, Load Set SID, Summary) member table with `+ Add` / `- Remove`.
- Live "Resolves to N member sets, totalling M load entries" preview.
- Restore / OK / Cancel.

Right-click on a load combination in the Loads tab now offers
**Edit... / Copy... / Rename... / Delete** instead of only Delete.
Edit and Copy are full undo/redo via new
`EditLoadCombinationCommand` and `CopyLoadCombinationCommand`.

### Preferences dialog for entity colors and sizes (item 3)
New `Edit > Preferences...` dialog (Ctrl+,) with three tabs:

- **Entity Colors** - per-type color picker for Shells / Beams / Bars
  / Rods / Bushes / Plates / RBE2 / RBE3 / Masses / Solids / Shear /
  Gap / Plotel / Free Nodes. Per-row reset to default.
- **Entity Sizes** - mass glyph scale (as % of model length, replaces
  the hardcoded `0.015`), node point size, beam line width, edge
  width, RBE line width, free-edge line width.
- **Highlight** - selection highlight color (was hardcoded peach;
  now user pick) + selection outline width.

All values persist via `QSettings("NodeRunner", "NodeRunner")` and
are read on app startup so changes survive restart. Restore Defaults
on the dialog resets every row to its baked-in default.

### MEGA-XWFF deck profile (headless)
On `MEGA-XWFF_04A_DUL_FUC-Entry.dat` (2.48M nodes, 2.40M elements,
1,962 load SIDs, 3,912 load combos, 512k load cards, 269k PLOAD4
face refs):

| Stage | v3.5.0 wall-clock |
|---|---|
| pyNastran parse (chunked-strict) | 253.7 s |
| _finalize_for_viewer | 1.3 s |
| build_node_coords_vectorized | 2.6 s |
| build_element_arrays_vectorized | 19.4 s |
| grid build + EID cell_data | 0.2 s |
| compute_grid_centers_and_normals | 0.8 s |
| Loads tab populate (lazy) | <0.1 s |
| **end-to-end first display** | **~278 s (≈ 4.6 min)** |

The parse phase is now the dominant cost (pyNastran's text parse on
a multi-GB BDF). Once parse completes, scene build + tree + load
actors are all under 30s combined. Future work: investigate
parallel/streaming parse to bring the parse phase down further.

### Tests
185 pass (170 prior + 7 load combo + 5 prefs + 3 visibility).

## Changelog for v3.4.2

Hotfix for an import hang on decks with thousands of load SIDs (e.g. `FUSE-XWFF_04A_DUL_FUL-Entry.dat` has 1,962 SIDs).

### Symptom
Even after v3.4.1's PLOAD4 fix, import would still creep through "Stage 6/6: Building 3D scene -> Building load actors (840/1962 SIDs)..." over minutes.

### Root cause
`_create_all_load_actors` created three `add_mesh(..., opacity=0)` actors per SID up-front. On 1,962 SIDs that's ~5,886 Qt/VTK actor allocations even though almost all of them were invisible (every SID's checkbox is off by default).

### Fix
v3.4.2 splits the function:
- **`_create_all_load_actors`** now only does the scaling pass (FORCE/MOMENT/PLOAD4 max magnitudes for arrow scaling) and clears caches. Zero per-SID actor work.
- **`_ensure_load_actor_for_sid(sid)`** (new) builds the three actors for ONE SID. Idempotent. Centers/normals/EID-lookup are cached on first call and reused.
- **`_update_plot_visibility`** calls `_ensure_load_actor_for_sid(sid)` lazily when the user actually toggles a SID checkbox on.

Net effect on the user's deck:
- Import time: cut from minutes-hung to seconds.
- First SID shown: small build cost (~50ms for typical SID).
- All subsequent toggles: instantaneous.

Workflows that visualize 1-5 SIDs at a time (the common case) now pay only for those SIDs. Workflows that want every SID built can still get there by checking everything in the Loads tab - the cost is just deferred to when the user actually requests it.

170 tests still pass.

## Changelog for v3.4.1

Hotfix for an import hang on decks with hundreds of thousands of PLOAD4 pressure-load faces (e.g. `FUSE-XWFF_04A_DUL_FDC-Boost-2.dat`).

### Symptom
Import would reach "Stage 6/6: Building 3D scene" and stick on "Building load actors (FORCE/MOMENT/PLOAD4)..." with Windows marking the app "Not Responding" indefinitely.

### Root cause
`_create_all_load_actors` iterates every PLOAD4 face and called `self.current_grid.get_cell(cell_idx)` to fetch the cell's center + points for the pressure-arrow normal. Each call wraps a VTK cell in a `pyvista.Cell` object - on a deck with 600k+ pressure face references that's 600k Python objects allocated, plus a Python-level cross-product per cell, plus per-cell normalize. With the GIL held the main thread blocks long enough that Windows raises the "Not Responding" overlay.

### Fix
v3.4.1 precomputes per-cell centers + shell normals once per grid in vectorized numpy:

- `cell_centers = grid.cell_centers().points` — one bulk VTK call.
- Shell normals (QUAD + TRIANGLE) computed via a single batched `np.cross(p1-p0, p2-p0)` over all shell cell indices simultaneously, using `grid.cell_connectivity` + `grid.offset` to gather the first three vertex indices per cell.
- Non-shell cells get a fallback `+Z` normal (PLOAD4-on-solid arrows still get *a* direction; magnitude/coloring is correct).

The per-PLOAD4 loop in `_create_all_load_actors` now does dict lookups + a single `np.asarray` + a batched index into the precomputed arrays. On the user's deck this drops the load-actor build from minutes-to-hung to under a second.

Per-SID progress is also emitted into the import dialog so monster decks with many load sets still show forward motion.

170 tests still pass.

## Changelog for v3.4.0

Feature release covering nine items from end-user testing on
`Fuse_BFEM_Complete.bdf` (1.28M elements).

### Entity selection
- **Ctrl+V into the bucket** - the EntitySelectionBar now has a
  Qt-standard paste shortcut. Pastes Femap-style range tuples
  (`1,10,2`), single IDs (`12345`), and multi-line lists; the IDs
  land in the bucket as signed range rows respecting the current
  Add/Remove/Exclude radio.
- **Pick > Paste** menu confirmed to work end-to-end; smoke test added.
- **Single-ID + More** - typing one ID and clicking More with To/By
  blank now reliably produces a `12345` bucket row (whitespace-tolerant
  parser + pytest regression).

### Model tree
- **Materials / Properties / Coordinate Systems branches always shown**
  even when empty - the count `(0)` makes absence visible at a glance
  instead of the branch silently disappearing.
- **Model summary in the import status** - new line
  `"N nodes, N elements, N properties, N materials, N coord systems"`
  so missing entity classes (typically in a missing INCLUDE) are
  obvious.

### List > Element Information
- **Type-dispatched node extraction** - RBE2 rows show
  `gn -> [Gmi]`, RBE3 rows show `[indep nodes] -> refgrid`, CONM2
  rows show `nid`. Before v3.4.0 the dialog blindly read `elem.nodes`,
  which doesn't exist on rigid elements - so RBE2/RBE3 rows had an
  empty Nodes column. PID/MID columns now show `-` for rigid/mass
  elements that have no PID.

### RBE integrity
- `_create_rbe_actors` tracks counts of RBEs whose center or any leg
  isn't in `model.nodes`. The import summary now appends a
  `"rigid integrity: N missing center, N with partial legs, ..."`
  line, so the "spider goes to random places" case can be
  attributed to dangling references in the deck vs the renderer.

### Auto Group fix (RBE2 vs RBE3 distinction)
- `_group_by_element_type`, `_group_by_property`, `_group_by_material`
  now walk both `model.elements` AND `model.rigid_elements`. Auto
  Group "By Element Type" produces distinct `TYPE RBE2` and
  `TYPE RBE3` groups - previously rigids were skipped entirely
  because the loops only saw `model.elements`.

### Shading
- Added a Shading checkbox to the Display tab (next to Style / Color
  By). Toggles in sync with the existing View menu action.
- The lit/unlit difference is now obviously visible: `flat`
  interpolation + full ambient (off) vs `phong` interpolation with
  high diffuse (on). Status line shows `"Shading: flat (unlit)"` /
  `"Shading: phong (lit)"` so the user knows it took effect.

### Groups - Femap-style redesign
- Bottom button grid replaced by a **top toolbar**:
  ```
  [New] [Delete] [Rename] | [Add ▾] [Remove ▾] [Clear]
  [Show All] [Isolate] [Highlight] | [Auto ▾] [Boolean ▾]
  ```
- Groups now carry **5 entity types**: `nodes`, `elements`,
  `properties`, `materials`, `coords`. Group rows in the list now
  show `(N elems, N nodes, N props, N mats)` when non-empty.
- New **`GroupAddDialog`** (Femap-style tabbed dialog) with one tab
  per entity type and rule-based selectors:
  - Node: By ID range / By connected to selected element /
    From selection / By group
  - Element: By ID range / By Property / By Material / By Type /
    From selection / By group
  - Property: By ID range / By Material / From selection / By group
  - Material: By ID range / From selection / By group
  - Coord: By ID range / From selection / By group
  Live count refreshes as the user changes selectors.
- New **Boolean** menu: Union / Intersect / Difference between two
  existing groups, applied across all 5 entity-type lists. Result
  goes into a new auto-named group `A + B`, `A & B`, `A - B`.
- **Isolate / hide** filters now honor `properties` and `materials`
  too - a group of `materials=[1]` automatically isolates every
  element whose property's material is 1.

### Sketch-first workflow followed
- The new Groups tab and the Add-to-Group dialog were prototyped
  in standalone PySide6 sketch scripts before the production code
  landed. Real widgets were re-rendered after edits and compared
  against the sketches.

### Tests
- 156 pre-existing + 4 new RBE listing tests + 3 new bucket
  regression tests + 7 new Groups tests = 170 total, all green.

## Changelog for v3.3.1

Hotfix for a v3.3.0 regression on decks with RBE elements.

### `Invalid array shape. Array 'EID' has length (L) but a length of (3L) was expected`
Symptom: importing a deck with RBE2/RBE3 elements (e.g. `Fuse_BFEM_Complete.bdf`) failed with the popup above. Pre-v3.3.0, RBE actors were `pickable=False` so this code path didn't exist.

Root cause: `_create_rbe_actors` built the PolyData as `pv.PolyData(np.array(points))` then assigned `.lines = ...`. The first call auto-generates one VERTEX cell per point, so a polydata with 2L endpoints + L lines ends up with 3L cells — and the L-entry `cell_data['EID']` array length didn't match. PyVista's array-shape validator raised, the exception bubbled to the post-parse popup.

Fix: construct PolyData with the lines on the ctor (`pv.PolyData(points, lines=flat_lines)`) so no auto-vertices are generated. RBE picking still works (the line cells carry the EID array as before), and the cell count finally matches the EID array length.

## Changelog for v3.3.0

Minor release: real-deck stability fixes, faster model-tree population, picking now works on every element type, and a Femap-style selection bucket.

### CORD2R safety
On the user's `FUSE-XWFF_04A_DUL_FDC-Boost-1.dat` deck, import failed with `Post-parse handling failed: Local unit vectors haven't been set. Type='CORD2R' cid=920000 rid=200000`. The deck had a CORD2R whose `rid` referenced a coord that wasn't defined in the bulk (likely in a missing INCLUDE).

`_finalize_for_viewer` is rewritten to topologically sort coords by `rid` dependency (Kahn's algorithm), assign `rid_ref` in parent-first order, and call `cs.setup()` explicitly. Dangling parents fall back to basic (cid=0) instead of leaving `rid_ref = None`, and the fallback is reported on `model._coord_resolution_warnings`. The import-summary status line now mentions `"N coord system(s) with dangling parent - treated as basic"` when this happens.

### Faster model-tree population
On a 600k-element deck, `_populate_tree` was running two full `Counter(... for e in all_elements.values())` walks just to build the `By Type` / `By Shape` group labels - ~1.2M Python attribute lookups. Those counts now piggy-back on `scene_build.build_element_arrays_vectorized`'s existing single pass over all elements: the helper returns `by_type_counts` / `by_shape_counts` dicts and `_populate_tree` reads them in constant time. The tree build is also wrapped with `setUpdatesEnabled(False)`, batched insertion via `addChildren()`, and explicit per-group expansion (no more `expandAll()` + collapse-just-coords).

### Picking works on every element type
RBE2/RBE3 spider lines and CONM2 mass glyphs were on dedicated actors with `pickable=False` and no `EID` cell_data. Both actors are now `pickable=True`, each cell carries the parent element's eid, and the area-pick paths (box / circle / polygon) also project RBE/mass center positions to screen so they're included in the selection rectangle.

### Femap-style signed bucket
The Entity Selection bucket now records each pick as a signed range entry (`mode`, `start`, `end`, `step`) and shows them as Femap-style triplets:

| Action | Bucket row |
|---|---|
| Add range 1..10 step 2 | `1, 10, 2` |
| Add single 12345 | `12345` |
| Remove range 1..10 step 1 | `-1, -10, 1` |
| Remove single 8 | `-8` |
| Exclude range 5..7 | `x5, x7, 1` |

Final selection on OK = walk entries in order, expand each range, unioning `+` entries and subtracting `-` / `x` entries. Add 1..10 then remove 8 yields {1..7, 9, 10}. Delete on a bucket row literally undoes that pick.

### Layout changes
- New `Select Visible` button in the action grid - adds every currently-visible entity (respecting tree on/off state) to the bucket as a compressed range or set of ranges.
- Action grid reorganized:
  - row 0: `Select All` | `Select Visible` | `Reset`
  - row 1: `More` | `Previous` | `Delete`
  - row 2: `Pick ▾` | `Method ▾` | `Grow ▾`
- `Hilite` button moves out of the 3x3 grid into the bottom row, left of OK / Cancel.

### Tests
- 156 tests pass (138 prior + 3 finalize/coord + 15 bucket).

## Changelog for v3.2.7

Patch release: Entity Selection window redesigned (take 2).

### Layout
The dialog now follows standard Qt dialog convention with OK and Cancel separated at the bottom-right. Action buttons are uniform-sized in a 3x3 grid:

```
+---------------------+ +-------------------+ +----------------------+
| Add / Remove / Exc  | |  entries list     | | [Sel All][Reset][Prev]|
| ID / To / By        | |                   | | [Pick v][Method v][More]|
| Group               | |                   | | [Delete][Hilite][Grow v]|
+---------------------+ +-------------------+ +----------------------+
Count: 0 / N
-----------------------------------------------------------------
                                              [OK]      [Cancel]
```

Changes vs v3.2.6:
- **Uniform button sizes** in the right column via equal column-stretch on the 3x3 grid. No `setFixedSize` anywhere. The H button (now `Hilite`) sizes itself to match the others.
- **OK/Cancel moved out of the action grid** to their own row at the bottom of the dialog, separated by a horizontal divider. Standard Qt convention.
- **Grow submenu** at row 2 col 2 replaces the inline `Angle | Adj | Conn | Grow | Shrink` row. Click `Grow ▾` for Adjacent / Connected / Grow / Shrink, plus `Set angle threshold...` to change the dihedral tolerance. Femap also tucks these behind menus.
- **`Hilite` button** paints amber only when checked, blending with the other buttons in its default state.
- **Count label kept** (it's useful) - sits between the body and the bottom OK/Cancel row.

### v3.2.5 mistake addressed
v3.2.5 created the new layout but never attached it to the body, leaving the right column blank. v3.2.7 was screenshot-tested before shipping: a standalone PySide6 sketch script builds the proposed dialog and saves a PNG. The real `EntitySelectionBar` was also rendered to PNG before the release was tagged.

### Tests
- 138 tests still pass.

## Changelog for v3.2.6

Revert release. v3.2.5's Entity Selection redesign had a bug where the new `right_col` layout was created but never added to the dialog body, so the entire right side (Select All / Reset / Previous / Pick / Method / More / Delete / Hilite / OK / Cancel) was missing on screen.

v3.2.6 reverts `node_runner/dialogs/selection.py` to the v3.2.4 state. The original Femap-style layout (3x3 grid: Select All / Reset / Pick+H | Previous / Delete / OK | More / Method / Cancel) is back. Everything else from v3.2.4 stays: the deferred quality-compute fix, the Hilite-count-label feedback, the tighter spacing, the model-tree search bar, the Stage 6/6 scene-build progress.

The redesign idea (uniform button sizes, OK/Cancel separated) will be redone properly in a future release once we have screen-tested the layout end-to-end before shipping.

## Changelog for v3.2.5 (withdrawn)

This release had a layout bug that hid the right-column buttons. It has been superseded by v3.2.6.

## Changelog for v3.2.4

Patch release: Entity Selection window redesign.

### Uniform-size buttons, OK/Cancel separated at bottom-right
- Every action button (`Select All`, `Reset`, `Previous`, `Pick`, `Method`, `More`, `Delete`, `Hilite`) is now the same size, driven by equal column-stretch on a 3-column grid. No more H button at 26x26 px next to wide buttons.
- OK and Cancel sit together on their own row at the bottom-right of the right column, separated from the action grid by a vertical stretch.
- Layout:
  ```
  [Select All]  [Reset]     [Previous]
  [Pick ▾]      [Method ▾]  [More]
  [Delete]      [Hilite]
  ---------- (spacer) ---------------
                              [OK]  [Cancel]
  ```
- The Hilite button (renamed from "H") only paints amber when checked, so it visually matches the other buttons in its default state instead of standing out.

### Result
- The dialog's default height matches the manually-resized version from v3.2.4 — no more excess vertical whitespace.
- Buttons are logically grouped: bulk-actions on top, input-methods in the middle, list-ops on the bottom; decision (OK/Cancel) clearly separated.

### Tests
- 138 tests still pass (no semantic changes; layout-only).

## Changelog for v3.2.4

Patch release covering six user-reported issues from the v3.2.3 review.

### `List > Element Information` hang (~5 min)
- Root cause: `_build_select_by_data` eagerly computed element-quality metrics (aspect, skew, warp, Jacobian, taper) for *every* element on every selection-bar open. On a 700k-element deck that's many minutes of pure-Python work for a "By Quality" category most users never click.
- Fix: defer the quality compute. The Method menu now lists `By Quality` as a placeholder; the slow calc runs only when the user actually picks that method, with a wait cursor + status message. Result is cached for the session.
- The selection bar's `Highlight` path also got a performance guard: when the user's selection covers more than 50% of the grid, we overlay the *full* `current_grid` wireframe instead of doing an `extract_cells` deep-copy. Visually equivalent, ~100x faster on whole-mesh selections.

### Highlight (H) button in the Entity Selection window
- The button was wired correctly but produced no visible feedback when nothing was selected yet, which read as "broken". Now the count label updates with a hint: `Count: 0 / 700,000  (nothing selected to highlight yet)` when H is toggled on with no entries, and `(highlighting OFF)` when toggled off.

### Entity Selection window vertical spacing
- Tightened the per-row padding (button `min-height` 24 → 20, `QLineEdit` 20 → 18, layout spacing 4 → 2) so the dialog isn't dominated by whitespace between rows.

### Coordinate Systems tree group: actually collapsed
- v3.2.3 marked the coords group `setExpanded(False)` but `_populate_tree` calls `tree.expandAll()` afterward, so the group re-expanded. v3.2.4 explicitly collapses the coords group *after* the tree-wide expand so it stays closed.

### Model-tree search bar
- New `QLineEdit` above the Model tree with **Filter** / **Find** radio buttons:
  - **Filter mode** (default): hides tree items whose text doesn't match the query (case-insensitive substring; matches GRIDs / CIDs / element types / property labels). Parents of matching children stay visible and auto-expand.
  - **Find mode**: Enter scrolls the next matching item into view + selects it. Wraps to top after the last match.
- Use it to find `CID 10`, `CQUAD4`, `Property 42`, or `Material steel` instantly without scrolling.

### "Still rendering" feedback after import completes
- The import dialog used to close at `Imported N nodes` and the user then saw the window freeze for a few seconds while `_update_viewer` built per-element actors (loads, constraints, RBEs, masses, PLOTELs). On a 268k-PLOAD4 deck that gap was 5-15 s of silent waiting.
- v3.2.4 keeps the import dialog open through the post-parse scene-build phase. New `Stage 6/6: Building 3D scene` with detail lines for each actor-creation step (loads, constraints, rigid elements, masses, PLOTELs, final render). Dialog closes only after the scene is on-screen.

### Tests
- 138 tests still pass (changes are wiring/UX; no semantic-level test changes).

## Changelog for v3.2.3

Patch release. Two small but visible fixes.

### Post-parse `'Cell' object has no attribute 'normal'` fix
- Symptom: import completes (vector geometry comes through), then a popup says `Post-parse handling failed: 'Cell' object has no attribute 'normal'`.
- Root cause: `_create_all_load_actors` builds pressure-load arrows for every PLOAD4 by reading `cell.normal` off a `pyvista.Cell` returned from `grid.get_cell()`. PyVista 0.46 doesn't expose a `.normal` attribute on `Cell`. The user's deck has 268,000 PLOAD4 cards so this fires on every real aerospace deck with shell-pressure loads.
- Fix: compute the shell normal manually from the cell's points (`np.cross(p1-p0, p2-p0)` then normalize). One-time cost per PLOAD4, no API dependency on PyVista internals.

### Coordinate Systems group: collapsed and off by default
- Big aerospace decks routinely have dozens of CORD2R systems (one per superelement boundary, one per skin panel). Showing all of them on first open clutters the tree and the 3D view.
- v3.2.3: the "Coordinate Systems" tree group now defaults to **unchecked** (all coord actors hidden) and **collapsed** (children not expanded). The user can expand and check individual systems as needed.

## Changelog for v3.2.2

Patch release. Two fixes; the headline is that **decks dominated by analysis-only cards (TEMP, TLOAD, tables, etc.) now import in seconds instead of "Not Responding forever"**.

### Strip analysis-only cards on import

The user's `Fuse_BFEM_Light.bdf` was 76% TEMP cards (4.6 million of 6.08 million total) — nodal temperature loads for downstream thermal analysis. Node Runner has no UI to render, edit, or otherwise use TEMP cards (or TLOAD, RLOAD, TABLED*, DCONSTR, DRESP, etc.). But v3.2.1 still fed all 4.6M to pyNastran's parser, which is what made the dialog say "Not Responding".

v3.2.2 strips a hardcoded set of analysis-only cards from the inlined buffer **before** pyNastran sees them. The raw lines are stashed verbatim and appended back on export, so the round-trip preserves every card. The skip list is intentionally not user-toggleable — a toggle would just be a footgun (flip it off, get a hang, no idea why).

The skip set covers:
- Thermal loads: `TEMP`, `TEMPD`, `TEMPF`, `TEMPB`, `TEMPRB`, `TEMPAX`, `TEMPP1..4`
- Dynamic loads: `TLOAD1/2`, `RLOAD1/2`, `DAREA`, `DPHASE`, `DELAY`
- Function tables: `TABLED1..4`, `TABLEM1..4`, `TABLES1`, `TABRND1`, `TABRNDG`, `TABDMP1`
- Optimization: `DCONSTR`, `DESVAR`, `DLINK`, `DRESP1..3`, `DVPREL1/2`, `DVCREL1/2`, `DVMREL1/2`, `DOPTPRM`
- Output sets / eigenvalue control: `SET1/2/3`, `BOUTPUT`, `EIGR`, `EIGRL`, `EIGB`, `EIGC`, `EIGP`

The post-import status line reports the count:
```
Opened Fuse_BFEM_Light.bdf - includes 4,419 rigid elements (RBE/RBAR),
  4,638,690 analysis-only cards stashed (TEMP: 4,638,690) - written back on export
```

When a future version of Node Runner adds UI for any of these card types, that type comes off the skip list and starts going through the full parser.

### "Last file read" regression fixed
- v3.2.1's dialog line "Last file read: 6,079,465 cards in ~8k-card chunks. Bad chunks fall back to lenient ..." was a regression: my `_emit` adapter was putting the message tail into the `source_file` field for legacy emits. v3.2.2 keeps `source_file` strictly for actual include-file basenames; a new `detail` field on `ImportProgress` carries the descriptive text. The detail line under the bar uses `record.detail` directly.

### Tests
- 12 new pytest cases in `tests/test_skip_cards.py` covering the skip set composition, card-name detection for all three formats (fixed, comma, 16-char), strip semantics, continuation handling, and an end-to-end round-trip that verifies TEMP cards survive `import → save → re-read`.
- **138 tests pass** (was 126).

## Changelog for v3.2.1

Patch release that fixes three v3.2.0 papercuts.

### Threaded BDF import (the real "Not Responding" fix)
- v3.2.0 vectorized the *scene build* but the *parse* was still synchronous - pyNastran held the GIL through long calls and Windows declared the app "Not Responding" on big multi-INCLUDE decks.
- v3.2.1 actually runs the import in a `QThread`. The dialog stays responsive throughout because the determinate progress bar updates only on signal emission (no marching-ants spinner; that was what drove us sync back in v3.0.0).
- Cancel button now wires to `thread.cancel()`. The worker's progress emit checks the flag and raises a `_Cancelled` exception that the thread catches cleanly. The viewer keeps its previous state.

### Persistent "Last file read" line
- The dialog now has a dedicated `Last file read: <filename>` line under the header. It's only ever updated, never cleared, so the user can always see which include the parser got furthest with - across stage transitions, errors, and silent stretches.
- Detail line under the bar no longer repeats the source file (it was duplicating the persistent line). Detail line now carries only the in-flight info: card type, line number, error reason, counter, rate, ETA.

### Honest card count when there's no BEGIN BULK
- v3.2.0's `_count_card_starts_in_text` returned 0 if the inlined buffer had no `BEGIN BULK` marker (which happens when a user opens a partial deck like `Fuse_BFEM_Light.bdf` whose entire content is INCLUDE statements). The dialog then read "from 0 cards in ~8k-card chunks", which is wrong and confusing.
- v3.2.1 falls back to counting all bulk-data card-start lines when no `BEGIN BULK` is present, so the dialog shows the real card count instead of zero.

## Changelog for v3.2.0

Minor feature release: the import dialog now actually tells you what's going on, and the viewer handles 600k+ node/element decks without going "Not Responding".

### Import dialog rebuilt
- **Header**: persistent "Importing: <filename>" so you always see which file the operation is working on.
- **Stage label**: high-level step (Stage 1/5 ... 5/5).
- **Progress bar**: determinate, 0..100%.
- **Detail line**: a real status — current source include file, card type, line number (when known from lenient), most recent pyNastran error, cards/sec rate, ETA. No more "Parsing flattened deck..." in both labels.
- **Cancel button**: actually cancels (sets a flag the parser checks between chunks). The parser exits cleanly, the viewer keeps its previous state.
- **Watchdog**: if no progress events arrive for 30 s, the detail line stamps "(no parser activity for N s)" so you know the parser is in a single big pyNastran call. The stamp updates live.

### Structured progress pipeline
New `ImportProgress` dataclass replaces the v3.1.x `(stage, message, fraction)` callback signature. Every emit site - inline, chunked, lenient - now publishes `source_file`, `card_type`, `line_number`, `error`, `counter_done`, `counter_total`, `rate`, `eta_seconds` as separate fields. The dialog composes the display lines from those fields, instead of trying to split a hand-formatted string on `" - "`.

When a chunk falls into lenient fallback, the lenient parser now forwards the latest skipped card's name and error message up to the dialog detail line in real time. So you see e.g. `from AllLoads.bdf | CBUSH | line 12345 | field 5 must be integer` while the slow fallback runs, instead of a frozen-looking bar.

### Vectorized scene build (the "Not Responding" fix)
The viewer's `_update_viewer` was doing two giant Python loops on every import:
1. `node.get_position()` per node (with internal matmul for each non-basic-coord node)
2. Per-element type dispatch with `[node_map[nid] for nid in elem.nodes]` per element

On a 600k-element deck this was ~30+ seconds of pure-Python work, with Qt only getting CPU between yield points. Windows would declare the app "Not Responding" before scene-build finished, and the user had to kill it from Task Manager.

New `node_runner/scene_build.py` does both in bulk numpy:
- **Node coords**: pulls every `node.xyz` into one flat array, then groups by `cp` and does one batched matmul per unique coord system. For decks where all nodes use the basic coord (almost all decks), it's just `np.array([...])`. Measured at 0.05 s on a 60k-node deck.
- **Element arrays**: sorts elements by kind once, builds per-kind `(n_elem, n_nodes)` arrays, then uses `np.searchsorted` to convert raw node IDs into grid indices in one vectorized op. Measured at 0.25 s on a 100k-element deck.

Combined: ~50-100x faster than the old loops. A 600k-element deck now builds in ~1-3 s instead of "never".

### Femap-style element LOD
New `apply_display_lod` in `scene_build.py`. Above a threshold (default 500k cells), the displayed mesh:
- Stride-samples shells (CQUAD4 / CTRIA3 / CSHEAR) so the rendered count stays around 200k.
- Extracts the outer surface of solids (CHEXA / CTETRA / CPENTA) so a 600k-cell solid mesh becomes ~150-250k surface triangles.
- Keeps beams / rods / bushes as-is (typically small counts).

The underlying `current_grid` still has every cell, so cell-pick, queries, mass props, and export work unchanged. Only the displayed grid (`current_display_grid`) is decimated.

Threshold is configurable via **Settings > Element LOD Threshold...**; persists in QSettings. Status bar shows "LOD active: showing 200k of 600k cells (full mesh kept for picking)" when active.

### Tests
- Three new pytest cases: `test_vectorized_node_coords_match_get_position` (vectorized coords match per-node `get_position` to 1e-9), `test_vectorized_element_arrays` (right cell counts + EID order), `test_apply_display_lod_passes_through_small_grids` (LOD is no-op below threshold).
- 126 tests pass.

## Changelog for v3.1.5

Patch release: fixes the post-parse crash on decks that use non-basic coord systems, and cleans up the import dialog so the top label and detail line aren't showing the same text.

### Post-parse `'NoneType' object has no attribute 'transform_node_to_global'`
- Symptom: a deck imports successfully (progress bar reaches 100%, no errors during parse), then the application throws a `'NoneType'.transform_node_to_global` error in a popup. Nothing displays.
- Root cause: Node Runner reads with `xref=False` for speed (full pyNastran cross_reference takes minutes on big decks). But `pyNastran.node.get_position()` requires `node.cp_ref` to be set, and `xref=False` leaves it as `None`. The scene-build path calls `get_position()` and crashes.
- Fix: new `_finalize_for_viewer(model)` helper does the minimum cross-reference work needed for visualization: walks every CORD2x's RID chain, then wires every GRID's `cp_ref` and `cd_ref` to the matching COORD object. Unknown coord IDs fall back to basic (CID 0) so we never crash, just degrade gracefully.
- The finalize step shows as **Stage 5/5: Resolving coordinate references** in the import dialog.

### Cleaner import dialog
- Before: the top label and the line under the progress bar showed the same text ("Parsing flattened deck..." in both places).
- After: the top label is the high-level stage ("Stage 2/4: Reading include files"), the line under the bar is the concrete detail ("foo.bdf (5 / 7 files, 237.1 / 246.6 MB)"). The split happens on the first `" - "` in the model-layer status string, so the model layer controls both fields explicitly.
- All Stage 1-5 messages now follow the `"Stage X/Y: <stage> - <detail>"` convention.

### Tests
- New `test_finalize_for_viewer_resolves_cp_ref` builds a deck with both basic-coord and non-basic-coord GRIDs, plus a GRID referencing a non-existent CID, and confirms `get_position()` works on all of them after the finalize step.
- 123 tests pass.

## Changelog for v3.1.4

Patch release: external-superelement support. Real aerospace decks with `.pch` punch files containing DMIG stiffness/mass matrices now import without silently losing entries.

### What was broken
A single DMIG matrix in a `.pch` file can have thousands of column-entry cards. With the chunked-strict path, those entries could span multiple chunks and the chunk-merge logic (first-write-wins by matrix name) would drop entries from chunks 2..N. Even worse, if the chunk fell into lenient fallback, pyNastran's `add_card` would queue the DMIG entries into `_dmig_temp` but `fill_dmigs()` was never called to materialize them - so the matrix ended up with zero entries.

### Fixes

#### DMIG-family card pre-extract
- Before chunking, the importer now scans the inlined buffer for every DMIG / DMI / DMIJ / DMIJI / DMIK card (plus their continuation lines) and pulls them out into a dedicated punch sub-deck.
- The DMIG sub-deck is parsed in one shot with `pyNastran.read_bdf(punch=True, validate=False)` - which DOES call `fill_dmigs()` automatically. Every DMIG entry survives.
- The non-DMIG bulk data goes through chunked-strict as before.

#### DMIG detection in both formats
- The DMIG-family detector handles both fixed-format (`DMIG    K25 ...`) and comma-delimited free-field (`DMIG,K25,...`) cards, as well as the 16-char `DMIG*` form used in real `.pch` files.

#### fill_dmigs in lenient path
- `_read_bdf_lenient` now calls `fill_dmigs()` at the end so any DMIG entries that arrived through the lenient path are materialized into actual matrix entries.

#### List-extending merge for SID-keyed dicts
- `loads`, `spcs`, `mpcs`, `tables_sdamping`, `random_tables` are now merged by extending the list when the same key (SID) appears across chunks, instead of dropping subsequent chunks' contributions.

#### Defensive DMIG matrix merge
- Same-name DMIG matrices across chunks (should be rare with pre-extract, but kept as a safety net) now combine their GCi/GCj/Real/Complex/GCs arrays instead of first-write-wins.

#### Import summary now reports DMIG / RBE / MPC / SE counts
- The status bar message after a successful import now reads e.g.: "Opened deck.dat - includes 4 DMIG matrix(es), 1,247 rigid elements (RBE/RBAR), 89 MPC equation(s), 3 superelement card(s)". So you can verify the SE data actually came through.

### Tests
- New `test_dmig_matrix_survives_chunked_path` in `TestChunkedCrossRefs` builds a 12,000-entry DMIG matrix (well over the 8k chunk size), injects a duplicate GRID to force the chunked path, and asserts all 12,000 entries are present in the final master. 122 tests pass.

## Changelog for v3.1.3

Patch release that fixes the v3.1.2 "stuck at 26% (Not Responding) with no real status" symptom on big multi-INCLUDE decks.

### Root cause
At 26%, the importer was attempting a whole-deck strict parse (one blocking `pyNastran.read_bdf()` call on the entire flattened buffer). On a 100k+ card deck this can take many minutes, the call holds the GIL the whole time so Qt can't repaint, the dialog title goes to "(Not Responding)", and our status text is stuck on whatever it said last. Worst of all, whole-deck strict almost always fails on real industrial aerospace decks (one bad card type kills the whole call), so we paid the freeze cost only to fall through to chunked-strict anyway.

### Skip the freeze entirely on big decks
- After the inline phase, we count card-starts in the flattened buffer. If the deck has >= 100k cards or >= 50 MB of flat text, we **skip whole-deck strict entirely** and go straight to chunked-strict. The 26% freeze disappears.
- Smaller decks still try whole-deck strict (it's the fastest happy path) and now show a clear pre-warning: "Stage 3/4: Parsing flat deck strictly (12,000 cards, 4.2 MB). Window will appear frozen for this step. Last file read: foo.bdf."

### When strict fails, show WHY
- We now capture pyNastran's actual exception message ("eid=12345 must be positive; elem=...") and put it in the status text: "Stage 4/4: Whole-deck strict failed (eid=12345 must be positive). Switching to chunked-strict...".
- Per-chunk failures inside chunked-strict are also captured. The first failure per source file becomes a `SkippedCard` entry visible in the post-import report, so you know which file (e.g. `Wing_SE_Complete.pch`) contains the problematic card and what pyNastran said about it.

### Stage labels + novice-friendly status
- All status messages now lead with a stage label:
  - Stage 1/4: Detecting field format
  - Stage 2/4: Reading N include files (M MB total)
  - Stage 3/4: Parsing flat deck strictly OR Stage 3/4: Big deck detected, skipping whole-deck strict
  - Stage 4/4: Chunked-strict parse (N cards in ~8k-card chunks)
- The chunked phase shows: "Parsing bulk data - Card 6,360,000 / 6,375,641 | from Wing_SE_Complete.pch | 22,500 cards/s | ETA 0.6 s".
- Whatever the model layer emits is now passed through to the dialog verbatim. v3.1.2 was overwriting the rich message with a generic "Parsing flattened deck (pyNastran)..." string in `workers.py`; that hardcoded override is gone.

## Changelog for v3.1.2

Patch release that fixes the v3.1.1 "stuck at 97%" symptom: when a chunk falls into lenient fallback, the bar now advances within the chunk instead of waiting for the whole chunk to finish.

### Real status during slow imports
- The import dialog now shows the current source file alongside the card counter, taken from the `$ -- begin INCLUDE` markers we already write into the inlined buffer. So instead of "Parsing flattened deck..." you see "Parsing - 6,360,000 / 6,375,641 cards | file: Wing_SE_Complete.pch | 22,500 cards/s | ETA 0.6 s".
- Cards-per-second rate is computed from a running EMA. ETA is derived from that rate and the cards remaining.
- Progress updates are time-rate-limited to ~4/sec so the dialog doesn't spend its budget in Qt repaints.

### Lenient fallback inside chunked emits per-card progress
- The v3.1.1 chunked parser only emitted progress between chunks. If chunk K ran into lenient (e.g. a slow card type), the bar sat at the chunk-(K-1) value for the entire lenient pass over those ~20k cards. That's why your screenshot froze at 97% for 20 minutes - one chunk near the end was grinding through lenient silently.
- v3.1.2 passes the progress callback INTO the chunk's lenient fallback so the bar interpolates within that chunk. You now see motion every ~500 cards even during slow phases.

### Smaller default chunk size
- Default chunk size dropped 20,000 -> 8,000. Worst-case stall on a single bad chunk shrinks from ~60 s lenient to ~25 s, and the bar advances roughly 2.5x more often during a healthy run.

### Inline-phase status shows file size + tree depth
- During the inline phase (the first 25% of the bar), status now reads "Reading   Fuse_Base_0.3.3.bdf (12 / 47 files, 92.3 / 184.0 MB)" with two-space indentation per INCLUDE depth so you can see the nesting.

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
