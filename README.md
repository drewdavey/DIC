# DIC Analysis Pipeline — FSR Tensile Coupons

Three-stage pipeline that turns raw VIC-3D DIC exports and MTS load-frame
data into ASTM D638 mechanical properties and plots. Run in order:
`DIC_Level1.py` → `DIC_Level2.py` → `DIC_Level3.py`.

Coupon IDs follow the pattern `P01-T<EXPOSURE><DIRECTION>-<REPLICATE>`,
e.g. `P01-TCL00-01` (Print 01, Control exposure, 0° direction, replicate 1).
Exposures: `CL` (Control), `SW` (Seawater), `UV`, `IS` (In-Situ).

**Pipeline split**: Level 1 extracts the raw DIC fields *and* pairs them
with the MTS load frame into a per-coupon stress/strain-ready record
(nothing truncated, nothing scaled yet), plus a raw-signal sanity check.
Level 2 owns all of the tunable signal processing — scaling, failure
truncation, optional smoothing — and computes the D638 mechanical
properties. Level 3 only plots and summarizes; it reads everything it
needs from Level 2's outputs (plus raw MTS files for bearing, which has
no DIC step of its own) and never recomputes any tensile property. The
scalar mechanical properties (E, yield, UTS, Poisson's ratio, ...) are
written **once**, into `FSR-SpecimenTesting.xlsx`, and that workbook is
the single source of truth every downstream part of Level 3 (per-coupon
plots, group plots, and the printed stat tables) reads from — there is no
separate per-coupon summary CSV to keep in sync.

**One CSV per coupon, columns appended as it moves through the levels**:
Level 1 writes `<DIC_DIR>/<coupon_id>.csv`, one row per DIC frame,
untruncated. Level 2 reads that same file and appends its own columns to
it in place — there's no separate `_L1.csv`/`_L2.csv` pair, no duplicate
`step`/`time_s` columns, and no reindexing to track between levels.
Failure truncation doesn't drop rows; it sets a boolean `kept` column
(`True` inside Level-2's analysis window) and leaves every Level-2-derived
column `NaN` outside it, so a coupon's full untruncated record and its
truncated analysis window live side by side in the same file. Per-coupon
scalars that used to be repeated down every row (`mts_peak_N`, `area_mm2`)
now live once each in `<DIC_DIR>/coupon_scalars.csv`; `i_uts` (the index of
UTS) was dropped entirely since it's just `argmax(stress_MPa)` within the
kept rows — recomputing it is cheaper than storing and keeping it in sync.

**Merge history**: `tensile_group_plots.py`, `printStatsAll.py`, and
`force_displacement_signal_plots.py` used to be separate scripts. Their
functionality now lives inside `DIC_Level3.py` (group plots + stat tables)
and `DIC_Level1.py` (Step C, raw-signal plots), respectively — those
standalone files have been removed. `matlab/tensile_plots.m` and
`matlab/bearing_plots.m` are an independent MATLAB re-implementation of
part of Level 3's plotting, kept for reference; they are not run as part
of this pipeline and aren't guaranteed to stay in sync with it.

---

## Level 1 — `DIC_Level1.py`

**Step A** — converts raw VIC-3D `.out` files into per-frame CSVs, written
next to each `.out` file in its coupon's raw-data folder.

**Step B** — pairs those per-frame CSVs with the MTS force/displacement
record, computes virtual axial + transverse extensometers, and writes one
full (untruncated) result CSV per coupon. No failure truncation or load
scaling happens here — that's Level 2's job, so it can be tuned without
re-running the slow per-`.out` export or the extensometer pass.

**Step C** — signal-inspection plots and a cross-check report comparing
the raw MTS file and the raw DIC-sync CSV's force/displacement channels
directly against each other. Independent of Steps A/B: it's a sanity
check on the *inputs*, not on anything Level 1 computes, and always
reruns when enabled (cheap — no vicpyx calls).

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
- Step C: for each coupon, plots the raw MTS force/displacement vs. time,
  and a twin-axis plot of the DIC sync CSV's raw (volts) vs. scaled
  (engineering-unit) force and displacement channels side by side; then
  writes one aggregate report row per coupon cross-checking the DIC-sync
  and MTS peak forces, the raw→scaled linear calibration for each DIC
  channel, and the dominant FFT frequency of each raw channel (useful for
  spotting sensor noise/aliasing).

**Switches** (independent, so the slow step never reruns just to rebuild
the fast one):
- `DO_EXPORT_FRAMES` / `OVERWRITE_FRAMES` — Step A. `OVERWRITE_FRAMES`
  defaults `False`: skip `.out` files whose CSV already exists.
- `DO_BUILD_L1` / `OVERWRITE_L1` — Step B. `OVERWRITE_L1` defaults `False`:
  skip coupons whose per-coupon CSV already exists. Flip to `True` to
  rebuild just the consolidated CSV (e.g. after changing a gauge-length
  constant) without re-exporting any `.out` files.
- `DO_SIGNAL_PLOTS` — Step C. `ZERO_DISPLACEMENT` (default `True`) shifts
  the displacement traces in Step C's plots so the first finite value is
  zero; doesn't affect the per-coupon CSV (its `disp_mm` is zeroed
  separately in Step B).

**Inputs**
- `<coupon_dir>/*.out` — VIC-3D full-field export, one per DIC frame.
- `<coupon_dir>/<coupon_id>.csv` — VIC sync CSV (analog channels @ DIC frame rate).
- `<MTS_DIR>/<coupon_id>*.txt` — MTS raw file: `disp_mm, force_N, output_V, time_s`.
- `FSR-SpecimenTesting.xlsx` — gauge thickness × width → cross-sectional area.
- `<coupon_dir>` lives under `DATA_ROOTS[exposure]`, which points at
  `DIC_DIR/raw/<project-folder>` — a local mirror of the raw VIC-3D project
  data. (Older notes may reference `G:\DrewDavey\...`; that external-drive
  location is stale — the raw data now lives under `DIC_DIR/raw` instead.)

**Outputs**
- `<coupon_dir>/<out_filename>.csv` — one CSV per `.out` file, written next
  to it, with columns `X, Y, Z, U, V, W, exx, eyy, exy, e1, e2, gamma, x, y,
  u, v, q, r, q_ref, r_ref` (sigma is used only to filter rows, not exported).
- `<DIC_DIR>/<coupon_id>.csv` — full per-frame record: `step, time_s,
  disp_mm, load_raw, strain_axial_raw, strain_transverse_raw`. Level 2
  reads this same file and appends its own columns to it (see below).
  (`DIC_DIR` is the `DIC/` folder next to `MTS/`, not inside the coupon's
  raw data folder.)
- `<DIC_DIR>/coupon_scalars.csv` — one row per coupon: `coupon, mts_peak_N,
  area_mm2`. The only place these two scalars are stored — Level 2 reads
  them from here instead of finding them repeated down every row of the
  per-frame CSV.
- `<FIGS_ROOT>/<coupon_id>/MTS_force_disp.png` — raw MTS force-vs-displacement
  curve, Step B's sanity check.
- `<FIGS_ROOT>/<coupon_id>/MTS_force_displacement_signals.png` — Step C:
  raw MTS force and displacement vs. time.
- `<FIGS_ROOT>/<coupon_id>/DIC_sync_force_displacement_signals.png` — Step C:
  raw vs. scaled DIC-sync force/displacement channels vs. time.
- `<DIC_DIR>/raw_dic_force_displacement_signal_report.csv` — Step C's
  per-coupon cross-check table (one row per coupon).

---

## Level 2 — `DIC_Level2.py`

Reads Level-1's per-coupon CSVs, scales raw load to force/stress, applies
failure truncation, and computes ASTM D638 mechanical properties. No
plotting here (see Level 3).

**What it does**
- Scales `load_raw` to `force_N` using the per-coupon scale factor
  (`mts_peak_N / max(|load_raw|)`, read from `coupon_scalars.csv` — no
  separate calibration pass needed), falling back to a combined
  `SCALE_N_PER_UNIT` if `mts_peak_N` is missing, and divides by `area_mm2`
  (also from `coupon_scalars.csv`) to get `stress_MPa`.
- Truncates each record: marks pre-load slack (load < 2% of peak) and
  post-fracture rebound (first post-UTS frame where load < 50% of peak) as
  outside the analysis window (`kept = False`) rather than dropping rows.
- **Smoothing is optional and off by default** (`APPLY_SMOOTHING = False`)
  — ASTM D638 doesn't call for filtering the stress-strain record, and the
  modulus/UTS/yield windows sit well clear of where the raw signal is
  noisiest. If a batch's raw signal genuinely needs it, flip
  `APPLY_SMOOTHING` on; `FILTER_METHOD` selects the algorithm:
  - `"median"` (recommended when enabled) — rolling median
    (`MEDIAN_WINDOW` = 31 frames). Doesn't ring or systematically
    undershoot a sharp peak the way an averaging-based filter does — though
    some undershoot at UTS is still possible if the window straddles into
    the retained post-fracture decline; raise/lower `MEDIAN_WINDOW` if the
    plotted curve looks over- or under-smoothed.
  - `"butterworth"` — zero-phase low-pass filter (`BUTTER_ORDER`,
    `BUTTER_CUTOFF`, fraction of Nyquist — must stay `< 1`), clipped to the
    raw data's range to suppress `filtfilt` ringing at the truncation
    edges. Because clipping only chops the amplitude of any overshoot
    rather than removing it, it can still leave a flat "shelf" artifact
    right at a sharp peak — median is the safer default if you turn
    smoothing on.
- Computes, from this truncated (and, if enabled, smoothed) signal:
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
- `<DIC_DIR>/<coupon_id>.csv` — Level-1 output (one per coupon).
- `<DIC_DIR>/coupon_scalars.csv` — per-coupon `mts_peak_N`, `area_mm2`.

**Outputs**
- `<DIC_DIR>/<coupon_id>.csv` — the same file Level 1 wrote, with new
  columns appended: `kept` (bool, inside the truncated analysis window),
  `force_N`, `stress_MPa`, `strain_axial`, `strain_transverse`
  (toe-corrected, and smoothed if `APPLY_SMOOTHING` is on), and
  `stress_MPa_unsmoothed`/`strain_axial_unsmoothed` (the same window
  *before* smoothing — identical to the smoothed columns when smoothing is
  off), kept only so Level-3 can draw its raw-vs-smoothed diagnostic
  overlay without recomputing anything. All of these are `NaN` where
  `kept` is `False`. Scalar properties are *not* repeated here — they live
  once per coupon in `FSR-SpecimenTesting.xlsx` (see below).
- `FSR-SpecimenTesting.xlsx` (`SPECIMEN_SHEET`) — each coupon's scalar
  properties (E, toe strain, yield stress/strain, UTS, strain at UTS,
  Poisson's ratio) are written into new columns on that coupon's existing
  row, matched by Specimen ID. Only those columns are touched — other
  rows, formulas, and formatting in the workbook are left alone. If the
  file is open elsewhere when Level 2 runs, this step is skipped with a
  warning rather than failing the whole run. **This is the only place
  scalar properties are stored** — Level 3 reads them back out of here
  rather than recomputing.
- `<DIC_DIR>/level2_group_stats.csv` — D638 §11.7 mean/std/count per
  (exposure, direction) group.

### Variant — `DIC_Level2_tmp.py`

A one-off alternate Level-2 run, used when the DIC sync CSV's own "Load"
channel for a batch is only trustworthy as a *shape* and its peak-only
rescale (Level 2's normal method) isn't good enough. Instead of scaling
`load_raw`, it maps the *entire* raw MTS force/displacement time series
onto each DIC frame: the DIC-sync peak-force row and the raw-MTS
peak-force row are treated as the same physical instant (the two files
don't share a clock), and each DIC frame's time offset from that anchor
is used to interpolate MTS force/displacement at the matching MTS-clock
time. Truncation, optional smoothing, and property calculations are
otherwise identical to `DIC_Level2.py`, and it writes to the exact same
files (the per-coupon CSV, `FSR-SpecimenTesting.xlsx`,
`level2_group_stats.csv`) — running it after `DIC_Level2.py` overwrites
Level 2's results with this method's. It also adds one column of its own,
`disp_mm_mts` (the peak-anchored MTS displacement it derives), without
touching Level-1's own `disp_mm` column. `DIC_Level2.py` itself is never
modified by running this; it remains the standard pipeline for future
batches. Needs the same raw DIC sync CSVs as Level 1 (`DATA_ROOTS`,
pointing at `DIC_DIR/raw/...`) plus `MTS_DIR`, in addition to Level-2's
usual per-coupon CSV input.

---

## Level 3 — `DIC_Level3.py`

Plot- and stats-only. Reads Level-2's per-frame curve CSVs and the scalar
properties Level-2 already wrote into `FSR-SpecimenTesting.xlsx`, and
produces per-coupon plots, group overlay/summary plots, and mean ± std
(CV%) property tables — for both tensile (D638) and pin-bearing (D953)
coupons. **Writes no CSVs**; the only files it writes besides plots are
`P01_MechanicalStats.xlsx` (the stat tables). Nothing tensile is
recomputed — bearing statistics are the one exception, computed here
directly from the raw MTS `.txt` files since bearing has no DIC/Level-2
step of its own.

**Switches**
- `DO_PER_COUPON_PLOTS` — per-coupon σ-ε and Poisson plots.
- `DO_GROUP_PLOTS` — the raw-MTS force-displacement group plot and all
  three DIC-derived group plots (curve overlay, property scatter,
  peak-strength).
- `DO_PRINT_STATS` — the D638 tensile + D953 bearing mean±std(CV%) tables
  (stdout) and the `P01_MechanicalStats.xlsx` export.
- `DIC_EXCLUDE` — coupon IDs to drop from the group plots and stat tables
  only (per-coupon plots still include them). Currently `P01-TCL45-01`,
  excluded for an anomalously large toe correction (see the `TODO` next to
  it) pending backup DIC data.

**What it does**
- Per-coupon: draws each coupon's toe-corrected σ-ε curve (truncated at
  UTS) with modulus/0.2%-offset/yield/UTS markers and an Airtech
  reference line, plus the pre-smoothing raw curve for comparison; and a
  −ε_xx vs ε_yy Poisson plot. Both read straight from the per-coupon CSV's
  `kept == True` rows — no spline fit, resampling, or other alteration of
  Level-2's data.
- Group (MTS): reads every `P01-T*.txt` directly and plots raw
  force-vs-displacement curves for all coupons, colored by exposure,
  styled by direction, with peak-force markers.
- Group (DIC): reads every coupon's `kept == True` rows (truncated at UTS,
  again unaltered) for the σ-ε overlay; reads `UTS_MPa`/`E_GPa` from the
  specimen sheet for the property-scatter and peak-strength plots.
- Print stats: builds tensile property rows from the specimen sheet
  (E, yield stress, UTS, strain at UTS, Poisson's chord ratio) plus
  peak load read directly from the raw MTS tensile files; builds bearing
  rows by toe-correcting each raw MTS bearing curve (D953 Appendix X1
  tangent method) and computing pin-bearing strength at 4% hole
  deformation (D953-19 §13.3) and at failure (§3.2.5). Both tables are
  aggregated by (exposure, direction) and by direction alone (all
  exposures pooled), printed to stdout, and exported to Excel.

**Inputs**
- `<DIC_DIR>/<coupon_id>.csv` — Level-2's per-frame curves (`kept == True` rows).
- `FSR-SpecimenTesting.xlsx` — Level-2's scalar properties, read back by
  column header (e.g. `"E (GPa)"`, `"UTS (MPa)"`) and matched on Specimen ID.
- `<MTS_DIR>/P01-T*.txt` — raw MTS tensile files (group force-displacement
  plot + peak-load stat column).
- `<MTS_DIR>/P01-B*.txt` — raw MTS pin-bearing files (bearing stats only;
  no DIC/Level-2 counterpart).

**Outputs**
- `<FIGS_ROOT>/<coupon_id>/stress_strain_DIC.png` — per-coupon σ-ε curve
  (toe-corrected, truncated at UTS) with modulus/yield/UTS markers, plus
  the pre-smoothing raw curve in light gray for comparison.
- `<FIGS_ROOT>/<coupon_id>/poisson_DIC.png` — per-coupon −ε_xx vs ε_yy plot.
- `<FIGS_ROOT>/tensile_mts_FD.png` — group raw-MTS force vs. displacement.
- `<FIGS_ROOT>/tensile_curves_DIC.png` — group σ-ε overlay.
- `<FIGS_ROOT>/tensile_summary_DIC.png` — group UTS/E property scatter.
- `<FIGS_ROOT>/tensile_peak_strength_DIC.png` — group UTS by exposure.
- stdout — tensile (D638) + bearing (D953) mean±std(CV%) tables.
- `P01_MechanicalStats.xlsx` — same stats, `Tensile` + `Bearing` sheets.

---

## Common configuration

All scripts share a `SWITCHES` block at the top (`PRINTS`, `EXPOSURES`,
`DIRECTIONS`, `REPLICATES`) used to select which coupons to process, and
a `PATHS` block — trimmed to only what each script actually touches.
`DIC_Level2.py` and `DIC_Level3.py` don't need the raw-data `DATA_ROOTS`
at all, since by the time they run everything they need is already in
`DIC_DIR` or the specimen sheet; `DIC_Level2_tmp.py` is the exception —
it needs both the raw DIC sync CSVs (`DATA_ROOTS`) and `MTS_DIR`, same as
Level 1, because it re-derives force from the raw MTS record.

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

More recently, `_L1.csv`/`_L2.csv` were merged into one file per coupon:
`<coupon_id>.csv`, written by Level 1 and appended to by Level 2 (see
"One CSV per coupon" above). If you have old `_L1.csv`/`_L2.csv` pairs on
disk from before this change, they're stale — re-run Level 1 then Level 2
to regenerate `<coupon_id>.csv` and `coupon_scalars.csv` (which replaces
the `mts_peak_N`/`area_mm2` columns that used to be repeated down every
row of `_L1.csv`, and the `i_uts` column that used to be repeated down
every row of `_L2.csv`).
