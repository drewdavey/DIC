# DIC Analysis Pipeline — FSR Tensile Coupons

Three-stage pipeline that turns raw VIC-3D DIC exports and MTS load-frame
data into ASTM D638 mechanical properties and plots. Run in order:
`DIC_Level1.py` → `DIC_Level2.py` → `DIC_Level3.py`.

Coupon IDs follow the pattern `P01-T<EXPOSURE><DIRECTION>-<REPLICATE>`,
e.g. `P01-TCL00-01` (Print 01, Control exposure, 0° direction, replicate 1).
Exposures: `CL` (Control), `SW` (Seawater), `UV`, `IS` (In-Situ).

---

## Level 1 — `DIC_Level1.py`

Converts raw VIC-3D `.out` files into per-frame CSVs.

**What it does**
- Walks each selected coupon's data directory and finds every `.out` file.
- Loads each one with `vicpyx`, filters out invalid points (`sigma < 0`),
  and exports the requested full-field variables (coords, displacements,
  strains, pixel coords) to a CSV.

**Inputs**
- `<coupon_dir>/*.out` — VIC-3D full-field export, one per DIC frame.

**Outputs**
- `<coupon_dir>/<out_filename>.csv` — one CSV per `.out` file, written next
  to it, with columns `X, Y, Z, U, V, W, exx, eyy, exy, e1, e2, gamma, x, y,
  u, v, q, r, q_ref, r_ref` (sigma is used only to filter rows, not exported).

---

## Level 2 — `DIC_Level2.py`

Pairs the Level-1 per-frame DIC data with the MTS force/displacement
record, computes virtual axial + transverse extensometers, and writes one
full (untruncated) result CSV per coupon. No failure truncation or load
scaling happens here — that's deferred to Level 3 so it can be tuned without
re-running this slower step.

**What it does**
- Reads the per-frame Level-1 CSVs and the VIC sync CSV (analog channels
  captured at the DIC frame rate).
- Reads the MTS `.txt` raw file to get peak force and looks up gauge
  cross-sectional area from the specimen spreadsheet (thickness × width).
- Builds two point extensometers from the reference-frame AOI centroid
  (axial gauge length 4.36 in / 110.7 mm, transverse 1.0 in / 25.4 mm,
  per ASTM D638 §5.2.1 / Annex A3.5.2) and computes engineering strain from
  marker displacement each frame.
- Saves a raw MTS force-displacement sanity-check plot.

**Inputs**
- `<coupon_dir>/<coupon_id>.csv` — VIC sync CSV (analog channels @ DIC frame rate).
- `<coupon_dir>/<coupon_id>-*.csv` — per-frame full-field DIC CSVs (Level-1 output).
- `<MTS_DIR>/<coupon_id>*.txt` — MTS raw file: `disp_mm, force_N, output_V, time_s`.
- `FSR-SpecimenTesting.xlsx` — gauge thickness × width → cross-sectional area.

**Outputs**
- `<DIC_DIR>/<coupon_id>_L2.csv` — full per-frame record: `step, time_s,
  load_raw, disp_mm, strain_axial, strain_transverse, mts_peak_N, area_mm2`.
  (`DIC_DIR` is the `DIC/` folder next to `MTS/`, not inside the coupon's
  raw data folder.)
- `<FIGS_ROOT>/<coupon_id>/MTS_force_disp.png` — raw MTS curve sanity check.

---

## Level 3 — `DIC_Level3.py`

Reads the Level-2 result CSVs, applies failure truncation and a light
rolling-median smoothing pass, computes ASTM D638 mechanical properties,
and produces per-coupon and group plots/summaries.

**What it does**
- Scales raw `load_raw` to `force_N` — per-coupon `mts_peak_N / raw sync
  peak`, falling back to a combined `SCALE_N_PER_UNIT` — and divides by
  `area_mm2` to get `stress_MPa`.
- Truncates each record: drops pre-load slack (load < 2% of peak) and
  post-fracture rebound (first post-UTS frame where load < 50% of peak).
- Smooths `force_N`/`stress_MPa`, `strain_axial`, and `strain_transverse`
  with a rolling median (`MEDIAN_WINDOW` = 51 frames), applied only to this
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

**Inputs**
- `<DIC_DIR>/<coupon_id>_L2.csv` — Level-2 output (one per coupon).

**Outputs**
- `<FIGS_ROOT>/<coupon_id>/stress_strain_DIC.png` — per-coupon σ-ε curve
  (toe-corrected, truncated at UTS) with modulus/yield/UTS markers.
- `<FIGS_ROOT>/<coupon_id>/poisson_DIC.png` — per-coupon −ε_xx vs ε_yy plot.
- `<DIC_DIR>/<coupon_id>_L3.csv` — per-frame signals only (`step, eps, sig,
  eps_t, i_uts`), for downstream MATLAB group plots. Scalar properties are
  *not* repeated here — they live once per coupon in `level3_summary.csv`
  and in `SPECIMEN_SHEET` (see below).
- `<FIGS_ROOT>/level3_summary.csv` (also written to `DIC_DIR`) — one row
  per coupon with all computed properties.
- `<FIGS_ROOT>/level3_group_stats.csv` (also written to `DIC_DIR`) — mean/
  std/count per (exposure, direction) group.
- `FSR-SpecimenTesting.xlsx` (`SPECIMEN_SHEET`) — each coupon's scalar
  properties (E, toe strain, yield stress/strain, UTS, strain at UTS,
  Poisson's ratio) are written into new columns on that coupon's existing
  row, matched by Specimen ID. Only those columns are touched — other
  rows, formulas, and formatting in the workbook are left alone. If the
  file is open elsewhere when Level 3 runs, this step is skipped with a
  warning rather than failing the whole run.

Group overlay plots (`tensile_curves_DIC.png`, `tensile_summary_DIC.png`)
are produced separately by `tensile_group_plots.py`, which (like
`matlab/tensile_plots.m`) now reads scalar properties from
`level3_summary.csv` rather than from the per-frame `_L3.csv` files.

---

## Common configuration

All three scripts share a `SWITCHES` block at the top (`PRINTS`,
`EXPOSURES`, `DIRECTIONS`, `REPLICATES`) used to select which coupons to
process, and a `PATHS` block pointing at the raw data roots, MTS/DIC/figs
directories, and the specimen spreadsheet.
