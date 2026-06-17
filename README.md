# DIC Analysis Pipeline — FSR Tensile Coupons

Three-stage pipeline that turns raw VIC-3D DIC exports and MTS load-frame
data into ASTM D638 mechanical properties and plots. Run in order:
`DIC_Level1.py` → `DIC_Level2.py` → `DIC_Level3.py`.

Coupon IDs follow the pattern `P01-T<EXPOSURE><DIRECTION>-<REPLICATE>`,
e.g. `P01-TCL00-01` (Print 01, Control exposure, 0° direction, replicate 1).
Exposures: `CL` (Control), `SW` (Seawater), `UV`, `IS` (In-Situ).

**Pipeline split**: Level 1 extracts the raw DIC fields *and* pairs them
with the MTS load frame into a per-coupon stress/strain-ready record
(nothing truncated, nothing scaled yet). Level 2 owns all of the tunable
signal processing — scaling, failure truncation, smoothing — and computes
the D638 mechanical properties. Level 3 only plots; it reads everything it
needs from Level 2's outputs and never recomputes anything. The scalar
mechanical properties (E, yield, UTS, Poisson's ratio, ...) are written
**once**, into `FSR-SpecimenTesting.xlsx`, and that workbook is the single
source of truth every downstream script (Level 3, `tensile_group_plots.py`,
`printStatsAll.py`, `matlab/tensile_plots.m`) reads from — there is no
separate per-coupon summary CSV to keep in sync.

---

## Level 1 — `DIC_Level1.py`

**Step A** — converts raw VIC-3D `.out` files into per-frame CSVs, written
next to each `.out` file on the raw data drive (unchanged location —
that's where VIC-3D's project files already live, and where the per-frame
CSVs belong).

**Step B** — pairs those per-frame CSVs with the MTS force/displacement
record, computes virtual axial + transverse extensometers, and writes one
full (untruncated) result CSV per coupon. No failure truncation or load
scaling happens here — that's Level 2's job, so it can be tuned without
re-running the slow per-`.out` export or the extensometer pass.

**What it does**
- Step A: walks each selected coupon's data directory, finds every `.out`
  file, loads it with `vicpyx`, filters out invalid points (`sigma < 0`),
  and exports the requested full-field variables (coords, displacements,
  strains, pixel coords) to a CSV.
- Step B: reads the per-frame CSVs and the VIC sync CSV (analog channels
  captured at the DIC frame rate); reads the MTS `.txt` raw file to get
  peak force and looks up gauge cross-sectional area from the specimen
  spreadsheet (thickness × width); builds two point extensometers from the
  reference-frame AOI centroid (axial gauge length 4.36 in / 110.7 mm,
  transverse 1.0 in / 25.4 mm, per ASTM D638 §5.2.1 / Annex A3.5.2) and
  computes engineering strain from marker displacement each frame; saves a
  raw MTS force-displacement sanity-check plot.

**Switches** (independent, so the slow step never reruns just to rebuild
the fast one):
- `DO_EXPORT_FRAMES` / `OVERWRITE_FRAMES` — Step A. `OVERWRITE_FRAMES`
  defaults `False`: skip `.out` files whose CSV already exists.
- `DO_BUILD_L1` / `OVERWRITE_L1` — Step B. `OVERWRITE_L1` defaults `False`:
  skip coupons whose `_L1.csv` already exists. Flip to `True` to rebuild
  just the consolidated CSV (e.g. after changing a gauge-length constant)
  without re-exporting any `.out` files.

**Inputs**
- `<coupon_dir>/*.out` — VIC-3D full-field export, one per DIC frame.
- `<coupon_dir>/<coupon_id>.csv` — VIC sync CSV (analog channels @ DIC frame rate).
- `<MTS_DIR>/<coupon_id>*.txt` — MTS raw file: `disp_mm, force_N, output_V, time_s`.
- `FSR-SpecimenTesting.xlsx` — gauge thickness × width → cross-sectional area.

**Outputs**
- `<coupon_dir>/<out_filename>.csv` — one CSV per `.out` file, written next
  to it, with columns `X, Y, Z, U, V, W, exx, eyy, exy, e1, e2, gamma, x, y,
  u, v, q, r, q_ref, r_ref` (sigma is used only to filter rows, not exported).
- `<DIC_DIR>/<coupon_id>_L1.csv` — full per-frame record: `step, time_s,
  load_raw, disp_mm, strain_axial, strain_transverse, mts_peak_N, area_mm2`.
  (`DIC_DIR` is the `DIC/` folder next to `MTS/`, not inside the coupon's
  raw data folder.)
- `<FIGS_ROOT>/<coupon_id>/MTS_force_disp.png` — raw MTS curve sanity check.

---

## Level 2 — `DIC_Level2.py`

Reads Level-1's per-coupon CSVs, scales raw load to force/stress, applies
failure truncation and a light rolling-median smoothing pass, computes
ASTM D638 mechanical properties, and writes everything downstream needs —
no plotting here (see Level 3).

**What it does**
- Scales `load_raw` to `force_N` using the per-coupon scale factor
  (`mts_peak_N / max(|load_raw|)`, both already in `_L1.csv` — no separate
  calibration pass needed), falling back to a combined `SCALE_N_PER_UNIT`
  if `mts_peak_N` is missing, and divides by `area_mm2` to get `stress_MPa`.
- Truncates each record: drops pre-load slack (load < 2% of peak) and
  post-fracture rebound (first post-UTS frame where load < 50% of peak).
- Smooths `force_N`/`stress_MPa`, `strain_axial`, and `strain_transverse`
  with a rolling median (`MEDIAN_WINDOW` = 11 frames), applied only to this
  truncated window. A median was chosen over a linear filter (e.g.
  Butterworth) because it doesn't ring or systematically undershoot a sharp
  peak the way averaging-based filters do — though some undershoot at UTS
  is still possible if the window straddles into the retained post-fracture
  decline; raise/lower `MEDIAN_WINDOW` if the plotted curve looks over- or
  under-smoothed.
- Computes, from this truncated, smoothed signal:
  - **Modulus E** (D638 §11.4) — slope of the linear region (0.05–0.3% strain).
  - **Toe compensation** (D638 Annex A1) — shifts strain origin using the
    modulus line's x-intercept.
  - **UTS** (D638 §11.2) — max stress on original area.
  - **0.2% offset yield** (D638 §A2.6).
  - **Poisson's ratio** (D638 §A3.10.1.3) — chord at εₐ = 0.002 over
    0.0005–0.0025, plus a least-squares slope for reference.
  - **Group stats** (D638 §11.7/§12.1) — mean/std per exposure × direction.
- Always recomputes/overwrites on every run — this step is cheap (pure
  pandas/numpy over an already-built CSV), so there's no overwrite switch;
  re-run freely while tuning truncation/smoothing/property settings.

**Inputs**
- `<DIC_DIR>/<coupon_id>_L1.csv` — Level-1 output (one per coupon).

**Outputs**
- `<DIC_DIR>/<coupon_id>_L2.csv` — per-frame curve data only: `step, eps,
  sig, eps_t, i_uts, eps_raw, sig_raw`. `eps`/`sig`/`eps_t` are
  toe-corrected and smoothed; `eps_raw`/`sig_raw` are the same truncated
  window *before* smoothing, kept only so Level-3 can draw its raw-vs-
  smoothed diagnostic overlay without recomputing anything. Scalar
  properties are *not* repeated here — they live once per coupon in
  `FSR-SpecimenTesting.xlsx` (see below).
- `FSR-SpecimenTesting.xlsx` (`SPECIMEN_SHEET`) — each coupon's scalar
  properties (E, toe strain, yield stress/strain, UTS, strain at UTS,
  Poisson's ratio) are written into new columns on that coupon's existing
  row, matched by Specimen ID. Only those columns are touched — other
  rows, formulas, and formatting in the workbook are left alone. If the
  file is open elsewhere when Level 2 runs, this step is skipped with a
  warning rather than failing the whole run. **This is the only place
  scalar properties are stored** — Level 3 and the group-plot scripts read
  them back out of here rather than recomputing.
- `<DIC_DIR>/level2_group_stats.csv` — D638 §11.7 mean/std/count per
  (exposure, direction) group.

---

## Level 3 — `DIC_Level3.py`

Plot-only. Reads Level-2's per-frame curve CSVs and the scalar properties
Level-2 already wrote into `FSR-SpecimenTesting.xlsx`, and saves per-coupon
plots. **Writes no CSVs.**

**Inputs**
- `<DIC_DIR>/<coupon_id>_L2.csv` — Level-2 per-frame curves.
- `FSR-SpecimenTesting.xlsx` — Level-2's scalar properties, read back by
  column header (e.g. `"E (GPa)"`, `"UTS (MPa)"`) and matched on Specimen ID.

**Outputs**
- `<FIGS_ROOT>/<coupon_id>/stress_strain_DIC.png` — per-coupon σ-ε curve
  (toe-corrected, truncated at UTS) with modulus/yield/UTS markers, plus
  the pre-smoothing raw curve in light gray for comparison.
- `<FIGS_ROOT>/<coupon_id>/poisson_DIC.png` — per-coupon −ε_xx vs ε_yy plot.

Group overlay plots (`tensile_curves_DIC.png`, `tensile_summary_DIC.png`)
are produced separately by `tensile_group_plots.py`, which (like
`matlab/tensile_plots.m`) also reads scalar properties from
`FSR-SpecimenTesting.xlsx` rather than from any per-coupon CSV.

---

## Common configuration

All three scripts share a `SWITCHES` block at the top (`PRINTS`,
`EXPOSURES`, `DIRECTIONS`, `REPLICATES`) used to select which coupons to
process, and a `PATHS` block — trimmed to only what each level actually
touches (e.g. Level 2 and Level 3 don't need the raw-data `DATA_ROOTS` or
`MTS_DIR` at all, since by the time they run everything they need is
already in `DIC_DIR` or the specimen sheet).

## File naming history

If you're looking at old notes or file listings: the consolidated
per-coupon CSV used to be called `_L2.csv` (written by a separate
`DIC_Level2.py` that only did the MTS/extensometer pairing), and the
truncated per-frame property CSV used to be `_L3.csv` with a
`level3_summary.csv`/`level3_group_stats.csv` pair holding scalar
properties. That pairing step is now part of Level 1, the truncation/
property step is now Level 2, `_L2.csv` → `_L1.csv`, `_L3.csv` → `_L2.csv`,
and `level3_summary.csv` was retired in favor of `FSR-SpecimenTesting.xlsx`
as the single scalar-property store.
