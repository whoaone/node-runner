"""Tooltip strings for MYSTRAN-touching UI surfaces.

Central registry so the wording stays consistent between Preferences,
the Run dialog, the status bar, and the menu items. Every string is
written for an FEA user who knows MSC Nastran but has never
touched MYSTRAN -- the tooltips explicitly call out where MYSTRAN
differs from the MSC/NX dialect.

All strings are HTML-compatible (Qt tooltips support a subset of
HTML 4 / CSS 2). Keep them short enough that the tooltip box doesn't
swallow the screen; ~6 lines max.
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Solutions
# ---------------------------------------------------------------------------

SOL_101 = (
    "<b>SOL 101. Linear Static</b><br>"
    "Apply loads, get nodal displacements + element stresses. "
    "Same meaning as <i>MSC Nastran SOL 101</i>. "
    "Equivalent to the <i>Static</i> analysis type in common "
    "pre/post tools. "
    "Use for sizing, stress checks, deflection limits."
)

SOL_103 = (
    "<b>SOL 103. Normal Modes</b><br>"
    "Eigenvalue extraction: natural frequencies + mode shapes "
    "(no loads). Same as <i>MSC SOL 103</i>. "
    "Requires an EIGRL card (handled automatically when you pick "
    "this SOL on a deck with no eigenvalue extractor)."
)

SOL_105 = (
    "<b>SOL 105. Linear Buckling</b><br>"
    "Critical-load eigenvalue extraction. Same as <i>MSC SOL 105</i>. "
    "Needs <b>both</b> a static load case <b>and</b> a buckling "
    "extractor (EIGRL). The eigenvalue λ multiplied by the applied "
    "load gives the buckling load P<sub>cr</sub>."
)

SOLS_NOT_SUPPORTED = (
    "<b>Not supported by MYSTRAN:</b><br>"
    "SOL 106 (nonlinear static), SOL 107-112 (complex eigenvalue / "
    "transient / frequency response), SOL 144 (static aeroelasticity), "
    "SOL 200 (optimization). The pre-flight scanner will block runs "
    "containing these solutions or their cards."
)


# ---------------------------------------------------------------------------
# Tuning PARAMs
# ---------------------------------------------------------------------------

SOLLIB = (
    "<b>SOLLIB</b>. solver library for the K·x = F system.<br>"
    "• <b>SPARSE</b>: scales to large DOF, low memory. Recommended.<br>"
    "• <b>BANDED</b>: faster for tiny tightly-banded decks; blows up "
    "on production-size models.<br>"
    "<i>NOT</i> the same as MSC's IntMKL/LAPACK selector. MYSTRAN's "
    "only two options are SPARSE and BANDED."
)

QUAD4TYP = (
    "<b>QUAD4TYP</b>. CQUAD4 shell element formulation.<br>"
    "• <b>MIN4T</b>: Mindlin–Reissner <i>with</i> transverse shear "
    "(recommended for most structural work).<br>"
    "• <b>MIN4</b>: Mindlin without transverse shear (thin shells; "
    "faster, less accurate near supports).<br>"
    "MSC Nastran's equivalent is the <code>QUAD4TYP</code> param too, "
    "but its default is MIN4."
)

WTMASS = (
    "<b>WTMASS</b>. weight-to-mass conversion factor.<br>"
    "If your model density is in <b>weight</b> units (lbm/in³, "
    "kg/mm³), set WTMASS = 1/g so MYSTRAN converts to mass units. "
    "Common values:<br>"
    "• in / lbm: <code>1 / 386.4 = 0.002588</code><br>"
    "• mm / kg : <code>1.0e-9</code> (or <code>1.0</code> if already mass)<br>"
    "Same meaning as <i>MSC PARAM,WTMASS</i>. Default 1.0 assumes "
    "you've already entered mass densities."
)

GRDPNT = (
    "<b>GRDPNT</b>. grid point id for mass-properties summary.<br>"
    "MYSTRAN computes rigid-body inertia about this grid and writes "
    "it to the <code>G R I D   P O I N T   W E I G H T   G E N E R A T O R</code> "
    "section of the F06. Use 0 to summarize about the origin. "
    "Same as <i>MSC PARAM,GRDPNT</i>."
)

EPSIL = (
    "<b>EPSIL</b>. numerical singularity tolerance (advanced).<br>"
    "Node Runner does <b>not</b> inject this automatically: MYSTRAN's "
    "syntax is two-field (<code>PARAM,EPSIL,&lt;eqn_set_id&gt;,&lt;tol&gt;</code>) "
    "and pyNastran's writer doesn't reproduce that shape. If you "
    "need to tune it, put a hand-written EPSIL card in your deck. "
    "MYSTRAN's default is fine for most work."
)


# ---------------------------------------------------------------------------
# Run-options dialog fields
# ---------------------------------------------------------------------------

ANALYSIS_SET = (
    "<b>AnalysisSet</b> to use for this run.<br>"
    "AnalysisSets bundle SOL, subcases, output requests, and PARAMs. "
    "Pick <i>(Default)</i> to use the SOL number on the deck itself. "
    "Edit AnalysisSets via <i>Analysis → Analysis Set Manager…</i>."
)

OUTPUT_DISP = (
    "<b>DISPLACEMENT</b> output request.<br>"
    "Per-node 6-DOF displacement vector. Always on for SOL 101/103/105 "
    ". this is the data that drives Color-By-Results contours and "
    "the Animation timeline."
)

OUTPUT_STRESS = (
    "<b>STRESS</b> output request.<br>"
    "Per-element stresses (CQUAD4, CTRIA3, CBAR, CBEAM, CHEXA, "
    "CTETRA, CPENTA, CSHEAR, CROD, CBUSH). Required for stress "
    "contours in the Result Browser. Increases F06/OP2 file size."
)

OUTPUT_STRAIN = (
    "<b>STRAIN</b> output request.<br>"
    "Per-element strains (same element coverage as STRESS). "
    "Needed only if you plan to plot strain contours."
)

OUTPUT_FORCE = (
    "<b>FORCE</b> output request.<br>"
    "Per-element internal forces / moments (in element coordinates). "
    "Useful for hand-checks, free-body diagrams, and beam-margin "
    "calculations."
)

WORKING_DIR = (
    "<b>Working directory</b><br>"
    "Where the translated BDF, F06, OP2, and run_meta.json land. "
    "Defaults to a timestamped folder under "
    "<i>Preferences → MYSTRAN → Scratch directory</i>. "
    "Your source decks are never modified."
)


# ---------------------------------------------------------------------------
# Preferences fields (extras the run dialog doesn't show)
# ---------------------------------------------------------------------------

SCRATCH_ROOT = (
    "<b>Scratch directory</b><br>"
    "Root folder where MYSTRAN run subfolders are created. "
    "Each run goes into <code>&lt;scratch_root&gt;/runs/&lt;timestamp&gt;_&lt;deck&gt;/</code>. "
    "Default: <code>~/.node_runner/mystran_runs</code>. "
    "Source models are never written to."
)

KEEP_INTERMEDIATES = (
    "<b>Keep intermediate files</b><br>"
    "When ON: MYSTRAN's working files (.L1*, .OP2, .ERR) stay in the "
    "run folder. When OFF: only BDF + F06 + log are kept. "
    "Recommend ON while you're learning MYSTRAN. the .ERR file is "
    "essential for diagnosing failed runs."
)

AUTO_LOAD_RESULTS = (
    "<b>Auto-load results</b><br>"
    "When a run finishes, automatically parse the OP2/F06 and push "
    "results into the Result Browser. Provides "
    "<i>Read Results on Analyze</i>."
)

AUTO_OPEN_BROWSER = (
    "<b>Auto-open Result Browser</b><br>"
    "Show the Result Browser dock automatically when results load. "
    "Turn off if you keep the browser docked already and don't need "
    "the focus change."
)

CLEANUP_AFTER_DAYS = (
    "<b>Auto-cleanup runs older than N days</b><br>"
    "Run folders in the scratch directory older than this are "
    "deleted on app start. Set <b>0</b> to disable cleanup. "
    "Doesn't touch the run_meta.json index used by Analysis History."
)

EXE_PATH = (
    "<b>MYSTRAN executable path</b><br>"
    "Full path to your MYSTRAN binary. Click <i>Detect</i> to search "
    "PATH and common install locations "
    "(<code>C:\\MYSTRAN\\mystran.exe</code>, "
    "<code>/usr/local/bin/mystran</code>, etc.). "
    "Node Runner does not bundle MYSTRAN. download it from "
    "the upstream project."
)


# ---------------------------------------------------------------------------
# Menu actions / status bar
# ---------------------------------------------------------------------------

ANALYSIS_RUN_MENU = (
    "<b>Run Analysis (MYSTRAN)…</b>   (Ctrl+R)<br>"
    "Pre-flight scan, then translate + solve via the configured "
    "MYSTRAN binary, then auto-load results. The MYSTRAN equivalent "
    "of MSC Nastran's analysis runner."
)

ANALYSIS_HISTORY_MENU = (
    "<b>Analysis History…</b><br>"
    "Browse prior MYSTRAN runs in the scratch dir. Double-click a "
    "row to reload its OP2/F06 results into the current model "
    "(matches by EID/GID hash; warns on mismatch)."
)

ANALYSIS_CONFIGURE_MENU = (
    "<b>Configure MYSTRAN…</b><br>"
    "Shortcut to <i>Preferences → MYSTRAN</i>. Set the executable "
    "path, scratch directory, SOLLIB / QUAD4TYP defaults, and "
    "auto-cleanup behavior."
)

STATUS_RESULTS = (
    "<b>Loaded results</b><br>"
    "When MYSTRAN (or an external OP2) feeds results into Node "
    "Runner, the active solution, source (OP2 vs F06 fallback), "
    "subcase count, and total displacement-vector entries land "
    "here. Empty when no results loaded."
)

PREFLIGHT_BANNER = (
    "<b>Pre-flight scanner</b><br>"
    "Walks your deck looking for cards MYSTRAN can't run (aero, "
    "nonlinear, contact, optimization, unsupported elements). "
    "<b>Blocking</b> issues stop the run; <b>warnings</b> mean the "
    "translator will drop something on export. Resolve blocking "
    "issues before clicking Run."
)
