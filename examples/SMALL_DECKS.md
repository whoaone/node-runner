# Small example decks (v5.0.0)

Five hand-built decks for quickly exercising Node Runner's
**pre + solve + post** workflow on something small enough to run in a
few seconds. All in US customary units (psi / in / lbm) so they pair
with the built-in MMPDS material library.

| File | Solver | Purpose | Analytical check |
|---|---|---|---|
| [`cantilever_beam.bdf`](cantilever_beam.bdf) | SOL 101 | 11-node / 10-CBAR cantilever, tip force | δ_tip ≈ **0.0404 in** (Euler-Bernoulli, PL³ / 3EI) |
| [`cantilever_modes.bdf`](cantilever_modes.bdf) | SOL 103 | Same geometry, normal modes via EIGRL | f₁ ≈ **81 Hz** for Al 6061-T6 (1.875² / 2πL² · √(EI / ρA)) |
| [`plate_uniform_pressure.bdf`](plate_uniform_pressure.bdf) | SOL 101 | 5×5 CQUAD4 simply-supported plate under 1 psi | centre δ ≈ **0.0438 in** (Kirchhoff α=0.00406) |
| [`column_buckling.bdf`](column_buckling.bdf) | SOL 105 | 20-element pinned-pinned column under 100 lbf | λ₁ ≈ **12.7** → P_cr ≈ **1271 lbf** (π²EI / L²) |
| [`plate_with_rbe2_spider.bdf`](plate_with_rbe2_spider.bdf) | SOL 101 | Cantilever plate with mixed CQUAD4 + CTRIA3 mesh + RBE2 spider | Smoke test for v5.0.0 shell-bucket split + RBE handling |

## Running through MYSTRAN

1. Configure MYSTRAN once: **Edit → Preferences… → MYSTRAN** → set the
   executable path → *Detect* to verify.
2. **File → Import…** one of these decks.
3. **Analysis → Run Analysis (MYSTRAN)…** (Ctrl+R).
4. Accept the pre-flight report (these decks pass with zero blocking
   issues and zero warnings).
5. Watch the progress dialog. Results auto-load when the run finishes.
6. Status bar lands at e.g. `Results: SOL 1 (OP2) | 1 subcase | 15 disp`.

## Material

Every deck uses **Al 6061-T6** (MAT1: E = 9.9e6 psi, ν = 0.33,
ρ = 0.098 lbm/in³). To regenerate from the MMPDS preset library:
File → New → New Material → **Load…** → "Al 6061-T6 sheet/plate".

## Why these specific decks?

- **Cantilever** is the textbook validation case — every FE textbook
  has the closed-form for it. If MYSTRAN's tip deflection doesn't land
  near 0.0404 in on `cantilever_beam.bdf` something is wrong with the
  integration, not your real model.
- **Plate with pressure** validates the QUAD shell formulation (MIN4T
  vs MIN4, configurable in Preferences → MYSTRAN). Centre node 22 / 23
  carries the maximum displacement.
- **Column buckling** validates the geometric-stiffness pipeline; the
  eigenvalue λ₁ should equal P_cr / P_applied = 1271 / 100 ≈ 12.7.
- **Modes deck** validates the Lanczos extractor (EIGRL). First five
  cycles should be sub-1000 Hz for this geometry.
- **RBE2 spider** exercises the mixed-shell mesh path that v5.0.0
  split into `shell_q4` and `shell_tri3` buckets — both element types
  appear in `current_grid.cell_data['type']` after import and both go
  through their own uniform-width vectorized compose path.
