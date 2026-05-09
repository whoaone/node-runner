# Node Runner v3.0.0 - Manual test plan

Things the automated suite can't verify (mostly visual / interactive
behavior). Work through this top-to-bottom on a fresh launch. Each
checklist item should pass.

Pre-flight:
- Activate venv: `.\venv\Scripts\Activate.ps1`
- Run: `python run.py`
- Window title says "Node Runner v3.0.0"
- Status bar shows `Model: Untitled | Nodes: 0 | Export: Short` (no Units field by default)

---

## 1. File I/O - small model

### 1.1 Open `examples/b747_60k.bdf`
- [ ] File > Open > select b747_60k.bdf
- [ ] Progress dialog appears with "Parsing BDF..." label
- [ ] Dialog disappears on its own when parse finishes (~3-5s)
- [ ] Status bar updates to `Model: b747_60k.bdf | Nodes: 58,968 | Export: Short`
- [ ] Model renders: fuselage + wings + tails + floor visible
- [ ] Adaptive LOD message flashes: `showing ~10k of 58,968 nodes`
- [ ] Ctrl+Z (Undo) does nothing visible (no model edits yet)

### 1.2 Save round-trip
- [ ] File > Save > pick a temp filename
- [ ] Export Options dialog appears with three radios pre-set to "Short"
- [ ] Hit Export
- [ ] Status bar message confirms save with field-format note

### 1.3 Long-field export
- [ ] File > Save > new filename
- [ ] Pick "Long Field (16-character)"
- [ ] Check "Remember this choice for future exports"
- [ ] Open the saved file in a text editor: every GRID line ends with `*` continuation
- [ ] File > Open the long-field BDF: status bar shows `Export: Long`

### 1.4 Free-field export
- [ ] File > Save > new filename
- [ ] Pick "Free Field (comma-separated)"
- [ ] Open in text editor: GRIDs are comma-separated
- [ ] Re-import: status bar shows `Export: Free`

---

## 2. Mesh editing (Theme A)

Use the b747_60k model.

### 2.1 Split QUAD into TRIA
- [ ] Model > Modify > Split QUAD into TRIA...
- [ ] Selection bar appears at the top
- [ ] Box-select ~50 elements on the wing
- [ ] Hit Accept
- [ ] Status: "Split N quads into trias"
- [ ] Visually: those elements now show as triangles
- [ ] Edit > Undo: quads return

### 2.2 Refine 1-into-4
- [ ] Model > Modify > Refine Elements (1-into-4)...
- [ ] Pick a small patch of elements
- [ ] Status: "Refined N element(s) (1 -> 4 each)"
- [ ] Free-edge report (Tools > Free Edge Report) shows no new free edges along the refined patch's interior

### 2.3 Smooth Nodes
- [ ] Model > Modify > Smooth Nodes...
- [ ] Set iterations=10, factor=0.5, "Pin nodes on free edges" checked
- [ ] Pick a patch of interior nodes
- [ ] Mass properties (Tools > Model Check > Mass Properties...) before vs after: total mass shifts by < 0.1%

### 2.4 Mirror Elements
- [ ] Generate a half-model: pick all elements with y > 0
- [ ] Model > Modify > Mirror Elements...
- [ ] Plane = Y, offset = 0
- [ ] Status: "Mirrored N element(s) across Y=0"
- [ ] Visually: full symmetric model present
- [ ] Tools > Coincident Nodes (equivalence): zero duplicates after merge

### 2.5 Copy Elements with rotation
- [ ] Pick a small group
- [ ] Model > Modify > Copy Elements...
- [ ] Set translate (5, 0, 0), rotation Z = 90 deg, center (0, 0, 0)
- [ ] Confirm copies appear shifted and rotated as expected

### 2.6 Insert Node on Edge
- [ ] Model > Modify > Insert Node on Edge...
- [ ] Type two existing connected node IDs (use List > Nodes if needed)
- [ ] Status: "Inserted node N on edge (a-b); split M element(s)"
- [ ] Affected quads now show as triangles around the new node

### 2.7 Solid quality (uses any solid model)
- [ ] Open or generate a model with CHEXA / CTETRA elements
- [ ] Tools > Model Check > Element Quality
- [ ] Volume / aspect / jacobian columns appear for solids

---

## 3. Loads & BCs (Theme C)

### 3.1 PLOAD1
- [ ] Model > Distributed Load (PLOAD1)...
- [ ] Pick a CBAR/CBEAM, fill in load type FZ, scale LE, X1=0/P1=100, X2=1/P2=100
- [ ] Save BDF, reopen, confirm PLOAD1 card preserved

### 3.2 PLOAD2
- [ ] Model > Element Pressure (PLOAD2)...
- [ ] Pick shell EIDs, set pressure 500
- [ ] Round-trip via Save > Open

### 3.3 SPCD
- [ ] Model > Enforced Displacement (SPCD)...
- [ ] Set DOFs=13, value=0.001 on a couple of nodes
- [ ] Round-trip; SPCD card present

### 3.4 MPC
- [ ] Model > Multi-Point Constraint (MPC)...
- [ ] Add three terms: (1, 1, 1.0), (2, 1, -1.0), (3, 1, 0.5)
- [ ] Status: "Added MPC SID 1 with 3 terms"
- [ ] Round-trip preserves it

### 3.5 RBAR / RBE1 / RSPLINE
- [ ] Model > Connections > Rigid Bar (RBAR)...: EID 999, GA=1, GB=2, CMA empty, CMB=123456
- [ ] Model > Connections > General Rigid (RBE1)...: a few independent + dependent nodes
- [ ] Model > Connections > Rigid Spline (RSPLINE)...: a chain of 4+ nodes
- [ ] All three appear in the model tree under Rigid Elements
- [ ] Round-trip via Save / Open

### 3.6 BOLT
- [ ] Model > Bolt Preload (BOLT)...
- [ ] Set BID, preload, EIDs
- [ ] Status message acknowledges (may note pyNastran version limitation)

---

## 4. Convert Units

### 4.1 m -> mm preset
- [ ] Open a small model with known dimensions
- [ ] Note PSHELL.t, MAT1.E, a node coordinate
- [ ] Tools > Convert Units...
- [ ] Pick "m -> mm  (length only)" preset
- [ ] Apply
- [ ] Coordinates ×1000, thickness ×1000, E unchanged (force factor 1 -> stress factor 1e-6)
- [ ] Edit > Undo: every value restores

### 4.2 m,kg,N -> mm,t,N preset
- [ ] Verify on a known steel material:
  - E from 2.1e11 Pa to 2.1e5 MPa
  - rho from 7800 kg/m^3 to 7.8e-9 t/mm^3

### 4.3 Unit hint
- [ ] Settings > Units...
- [ ] Type "mm / N / t", save
- [ ] Status bar shows `Units: mm / N / t`
- [ ] Settings > Units... again, blank the field
- [ ] Status bar Units field disappears entirely

---

## 5. Cross Section

### 5.1 Basic clipping
- [ ] Load b747_60k
- [ ] View > Cross Section... (or Ctrl+Shift+X)
- [ ] Enable cross section
- [ ] Pick normal +X
- [ ] Slide the slider: model "chops" through, internal frames/stringers/stanchions visible
- [ ] Pick normal +Y: now slicing side-to-side
- [ ] Pick normal +Z: horizontal slice; floor reveals floor beams + stanchions
- [ ] Close dialog: model returns to whole

### 5.2 With ghost mode
- [ ] Create a group named "Wings" with all wing elements
- [ ] Hide the Wings group via the Groups panel
- [ ] View > Ghost Hidden Groups: ON
- [ ] Cross Section enabled
- [ ] Hidden wing renders as translucent wireframe inside the section view

---

## 6. Result browser + probe + animation (Theme B)

You need an OP2 file for these. If you don't have one, skip to section 7.

### 6.1 Load OP2
- [ ] File > Load Results (OP2)...
- [ ] Three docks appear on the right: Results, Animation, Vectors

### 6.2 Result browser
- [ ] Results dock: Nodes tab populates with NID, |U|, Ux, Uy, Uz columns
- [ ] Click "Ux" header: rows sort by Ux
- [ ] Double-click the top row: 3D view zooms to that node, highlight peach
- [ ] Click an element in 3D view: corresponding row highlights in the table

### 6.3 Probe tool
- [ ] Tools > Probe Mode (toggle on)
- [ ] Hover a node: tooltip with NID, xyz, displacement
- [ ] Hover an element: tooltip with EID, type, PID, von Mises

### 6.4 Custom contour expression
- [ ] Color By > Results
- [ ] In the Results dock expression bar, type: `1 - stress_vm / 450e6`
- [ ] Press Enter
- [ ] Status bar: "Contour expression: 1 - stress_vm / 450e6"
- [ ] 3D view recolors with the safety-factor scalar
- [ ] Clear the expression: original contour returns

### 6.5 Vector overlays
- [ ] Vectors dock: tick "Disp"
- [ ] Blue arrows appear scaled to displacement magnitudes
- [ ] Tick "Principal": peach arrows at element centroids in real principal directions
- [ ] Adjust length scale: arrows resize live

### 6.6 Animation timeline
- [ ] Pick a modal subcase
- [ ] Animation dock: hit Play
- [ ] Mode shape oscillates smoothly
- [ ] Pause, drag scrubber: model freezes at that phase
- [ ] Switch speed combo (0.25x..4x): playback speed changes

---

## 7. Performance / 1M model

### 7.1 1M load
- [ ] File > Open > examples/b747.bdf
- [ ] Progress dialog runs for 1-3 minutes (parse) plus 30-90s (scene build)
- [ ] Window stays responsive throughout (drag, resize work)
- [ ] Status bar updates with progress messages during scene build
- [ ] Final state: 1,025,008 nodes, LOD active
- [ ] Cross-section still scrubs smoothly at 1M

### 7.2 Settings persistence
- [ ] Toggle: Settings > Adaptive Node LOD off
- [ ] Close app, reopen
- [ ] Setting persisted (Adaptive Node LOD still off)

---

## 8. Command Palette (Ctrl+P)

- [ ] Press Ctrl+P
- [ ] Type "split": Split QUAD into TRIA tops the list
- [ ] Type "convert": Convert Units... appears
- [ ] Type "ghost": Ghost Hidden Groups appears
- [ ] Press Enter on the highlighted entry: action fires
- [ ] Esc closes the palette without acting

---

## 9. Undo/redo across themes

- [ ] Run several Theme A operations in a row (split, refine, smooth)
- [ ] Run several Theme C operations (PLOAD1, MPC, RBAR)
- [ ] Run Convert Units
- [ ] Press Ctrl+Z several times
- [ ] Each undo step rewinds one command
- [ ] Press Ctrl+Y same number of times: model returns to most recent state
- [ ] No dialog gets stuck open between undos

---

## 10. Edge cases / failure modes

- [ ] Open a non-BDF file (e.g. a .py): error dialog, no crash
- [ ] Cancel a parse mid-way: progress dialog closes, status updates, model not loaded
- [ ] Smooth Nodes with nothing selected: status says "no nodes selected", no crash
- [ ] Mirror with no elements selected: ditto
- [ ] Convert Units with factor 0: dialog rejects ("factors must be positive")
- [ ] Custom expression with `__import__`: status bar shows "Expression error: Disallowed syntax"

---

## When you're done

If everything above passes, the v3.0.0 release is solid. If any item fails,
note which check and the steps to reproduce - that's the punch list.
