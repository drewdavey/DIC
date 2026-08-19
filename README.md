# DIC Analysis Pipeline — FSR Test Coupons

Turns raw VIC-3D DIC exports and MTS load-frame data into ASTM mechanical
properties and plots. **Two pipelines**, one per test type, sharing the same
level structure, the same per-coupon-CSV convention, and the same output
directory:

| | scripts | standard | status |
|---|---|---|---|
| **Tensile** | `TensileDIC_Level1.py` → `TensileDIC_Level2.py` → `TensileDIC_Level3.py` | ASTM D638 | complete |
| **Flexural** | `FlexuralDIC_Level1.py` → `FlexuralDIC_Level2.py` | ASTM D790 | Levels 1–2; **no Level 3 yet** |

`tensile_modulus_sensitivity.py` sits beside the tensile pipeline rather than
in it: it is read-only, measures how much of the scatter in E comes from
processing choices rather than the material, and is documented under
[Result — how much of the tensile modulus scatter is processing?](#result--how-much-of-the-tensile-modulus-scatter-is-processing)
below. Changes made to the tensile pipeline on 2026-08-29, and the reasoning
behind each, are in `CHANGELOG_tensile_2026-08-29.md`.

Run each pipeline's levels in order. The two never read or write each other's
files — coupon IDs keep them apart — so a flexural run cannot disturb tensile
results and vice versa.

Pin-bearing (ASTM D953) coupons have **no DIC step at all**; their statistics
are computed straight from the raw MTS files inside `TensileDIC_Level3.py`.

Coupon IDs follow the pattern `P01-<TYPE><EXPOSURE><DIRECTION>-<REPLICATE>`,
e.g. `P01-TCL00-01` (Print 01, **T**ensile, Control exposure, 0° direction,
replicate 1). Type: `T` tensile, `F` flexural, `B` bearing. Exposures: `CL`
(Control), `SW` (Seawater), `UV`, `IS` (In-Situ). Only `CL` and `IS` were
bend-tested, and only at `00`/`90` — 12 flexural coupons in all.

---

## ⚠ The P01 DIC sync CSVs carry bad analog data — all three test types

Every P01 coupon's per-coupon VIC sync CSV (`<coupon_dir>/<folder>.csv`) has
unreliable analog channels: there was a **connection issue** during the batch.
This is the single most important thing to know before reading either pipeline,
because it is why the two derive force so differently.

The sync CSV's one *good* column is its clock. `Time_0_0` is a real epoch
timestamp stamped on each frame trigger, and it is the only thing that puts the
DIC frames on a time axis at all. Everything else in that file is suspect.

⚠ **That clock did not survive the trip into the per-coupon CSVs, and it has
now been fixed.** Tensile Level 1 wrote `Time_0_0` verbatim through
`to_csv(float_format="%.6g")`; six significant figures rounds a ~1.7756e9 epoch
to the nearest 10 000 s, so **every row of every tensile per-coupon CSV written
before 2026-08-29 carries the same `time_s`**. Level 1 now stores *elapsed
seconds from the first frame* (spans are 40–70 s, so `%.6g` keeps
sub-millisecond resolution) and keeps the absolute start as `t0_epoch_s` in
`coupon_scalars.csv`. Level 2 detects a degenerate `time_s` and re-reads
`Time_0_0` from the raw sync CSV, so it works on the old files too. The
flexural pipeline always stored elapsed time and was never affected.

**Tensile** — the load channel is trustworthy in *shape* but not in scale, so
it is rescaled. `TensileDIC_Level2.py` now does this by **regressing** the raw
MTS force against `load_raw` over the rising ramp (see Level 2 below); it used
to do it on peak alone (`mts_peak_N / max(|load_raw|)`), which is still
reachable via `LOAD_SCALE_MODE = "peak"`. `TensileDIC_Level2_tmp.py` goes
further and maps the entire raw MTS series onto the frames by anchoring the two
peak-force rows and back-tracking to the start of test.

⚠ **Three tensile coupons have DIC records that stop before the specimen
fails** — `P01-TCL45-01` (spans 90.6 % of the MTS peak force), `P01-TSW00-01`
(96.3 %) and `P01-TSW00-02` (97.2 %); every other coupon is above 99.8 %. For
those three, `max(|load_raw|)` is not the peak, so the peak-anchored scale is
inflated by ~10 %, ~4 % and ~3 %, and their reported **UTS and strain-at-UTS
are not the specimen's under any scale**. Level 2 detects and reports this
(`DIC_COVERAGE_MIN`); nothing can recover it from these recordings. This is why
`P01-TCL45-01` looked anomalous — not its toe correction.

**Flexural** — worse: there is **no load channel at all**. `Dev1/ai2`, the
input the DAQ config labels LOAD, is a scaled copy of `Dev1/ai1` (the
displacement input) at r > 0.99 — nothing was wired to it, so what appears on
it is multiplexer crosstalk from the input beside it. It shows no touchdown
knee where the MTS sits flat on its tare, and it ramps at 81–103 % of its
loaded rate during the approach travel. A load cell cannot do that. Peak
anchoring is therefore not available either, so the flexural pipeline aligns on
the **last correlated frame** instead — see Flexural Level 1 below. Force comes
entirely from the MTS load cell.

**Bearing** — no DIC step exists, so nothing depends on it.

*Worth fixing before the next batch*: wire the load-cell analog output to an
input, and either slow the DAQ scan rate or leave a grounded channel between
the two live inputs to kill the ghosting. The tensile batch has the same
crosstalk; it is only harmless there because the ghost lands on an unused
channel.

**Frame rate.** All flexural DIC recordings should be 5 Hz. Verified: all 12
coupons measure 5.000 Hz median on `Time_0_0`. Four of them (`FCL9001`,
`FIS9001`, `FIS9002`, `FIS9003`) have a handful of short intervals — dropped or
early triggers — but the median is unaffected. Flexural Level 1 checks the rate
per coupon against `EXPECTED_FRAME_RATE_HZ` and warns rather than assuming;
either way it works from each coupon's own measured clock, so a coupon that
drifts off 5 Hz still processes correctly.

---

## Tensile pipeline — ASTM D638

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
functionality now lives inside `TensileDIC_Level3.py` (group plots + stat tables)
and `TensileDIC_Level1.py` (Step C, raw-signal plots), respectively — those
standalone files have been removed. `matlab/tensile_plots.m` and
`matlab/bearing_plots.m` are an independent MATLAB re-implementation of
part of Level 3's plotting, kept for reference; they are not run as part
of this pipeline and aren't guaranteed to stay in sync with it.

---

### Level 1 — `TensileDIC_Level1.py`

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

  **On the 4.36 in axial gauge**: these are *not* D638 Type I dogbones. The
  gauge section is ~34 mm wide × ~14 mm thick (≈480 mm²; a Type I is ≈41 mm²)
  and the correlated ROI runs 137–149 mm along the loading axis, so 110.7 mm
  sits comfortably inside it on every P01 coupon. The 2.00 in gauge in D638
  Fig. 1 is specified for the Type I geometry and does not carry over. (A
  comment in the source used to claim a 50 mm gauge, contradicting the
  constant beside it; that comment was the error, not the number.)

  The extensometer carries two guards that **do not change any P01 result** —
  both failure modes were checked against the frame CSVs and neither fires —
  but protect a batch with worse correlation: endpoint candidates are filtered
  on `sigma` (inert on P01, whose exports carry no `sigma` column despite
  `EXPORT_VARS` requesting it — it takes effect only after a Step-A re-export),
  and the two endpoints are located once in the reference frame and thereafter
  tracked by their reference-frame coordinates rather than re-searched every
  frame, so a lost point cannot silently substitute a neighbouring subset and
  put a discrete step in the strain record. Both report counts if they trigger.
  `ext_endpoints` likewise now warns when a gauge is clipped to the ROI (no P01
  coupon clips) and the length actually used is recorded per coupon.
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
- `FSR-SpecimenTesting.csv` — gauge thickness × width → cross-sectional area.
  The **CSV**, not the `.xlsx` — see the geometry warning near the end of this
  README. Falls back to the workbook only if the CSV can't be read, and refuses
  outright rather than proceeding with NaN geometry.
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
  raw data folder.) `time_s` is **elapsed seconds from the first frame**, not
  the raw epoch — see the sync-CSV warning above.
- `<DIC_DIR>/coupon_scalars.csv` — one row per coupon: `coupon, mts_peak_N,
  area_mm2, t0_epoch_s` (the absolute start `time_s` is measured from), and
  `axial_gauge_mm` / `trans_gauge_mm` (the gauge lengths *actually* used, which
  differ from the requested ones only if `ext_endpoints` clipped them to the
  ROI). The only place these scalars are stored — Level 2 reads them from here
  instead of finding them repeated down every row of the per-frame CSV. Written
  with `float_format="%.12g"`, wider than the `%.6g` used elsewhere, because
  `t0_epoch_s` is an epoch and would otherwise be rounded to the nearest
  10 000 s — the same rounding that flattened `time_s`.
- `<FIGS_ROOT>/<coupon_id>/MTS_force_disp.png` — raw MTS force-vs-displacement
  curve, Step B's sanity check.
- `<FIGS_ROOT>/<coupon_id>/MTS_force_displacement_signals.png` — Step C:
  raw MTS force and displacement vs. time.
- `<FIGS_ROOT>/<coupon_id>/DIC_sync_force_displacement_signals.png` — Step C:
  raw vs. scaled DIC-sync force/displacement channels vs. time.
- `<DIC_DIR>/raw_dic_force_displacement_signal_report.csv` — Step C's
  per-coupon cross-check table (one row per coupon).

---

### Level 2 — `TensileDIC_Level2.py`

Reads Level-1's per-coupon CSVs, scales raw load to force/stress, applies
failure truncation, and computes ASTM D638 mechanical properties. No
plotting here (see Level 3).

**What it does**
- Scales `load_raw` to `force_N` by **regressing the raw MTS force against
  `load_raw` over the rising ramp** (10–85 % of peak), using every point in it
  — `LOAD_SCALE_MODE = "regress"`. It divides by `area_mm2` (from
  `coupon_scalars.csv`) to get `stress_MPa`.

  **Why this replaced the peak ratio.** The old rule
  (`mts_peak_N / max(|load_raw|)`, still available as
  `LOAD_SCALE_MODE = "peak"`) set the entire stress axis — and therefore E, UTS
  and yield, all linearly — from the ratio of two single samples: the largest
  of a few thousand noisy MTS samples over the largest of a few thousand noisy
  sync samples, from two records that don't share a clock. Both maxima are
  biased high by their own noise floors and there is no reason the biases
  match. The load cell and DAQ gain are hardware constants, so every coupon in
  a batch must return the same scale. Measured over the 35 P01 tensile coupons:

  | | mean | CV | range |
  |---|---|---|---|
  | peak ratio | 555.3 | **2.55 %** | 540.0 – 613.8 |
  | regressed | 559.4 | **0.71 %** | 551.0 – 566.6 |

  The fitted **intercept is deliberately not applied** — a DC offset on the
  stress axis cannot change a slope; it is absorbed by the toe correction, and
  is printed only as a check on the sync channel's zero. Forcing the fit
  through the origin instead (`regress0` in the sensitivity sweep) is both
  biased high and noisier, because `load_raw` carries a tare offset.

  The two clocks are anchored peak-to-peak and the residual lag refined by
  cross-correlation. **`LAG_SEARCH_FRAC` must stay generous (0.25).** At ±5 %
  of test duration the search rails against its own limit on 3 of 35 coupons
  (true lags −12.3, −6.0, −5.5 s) and returns a scale 5–17 % wrong while R²
  stays at 0.996–0.997 — so an R²-only guard passes it. A railed optimum is
  detected explicitly and falls back to the peak ratio rather than being used.
- **Checks DIC coverage** (`DIC_COVERAGE_MIN = 0.98`): the largest MTS force
  the DIC record actually spans, as a fraction of the MTS peak. Below that, the
  sync CSV stopped before the specimen did — see the warning at the top of this
  README for the three P01 coupons this catches. Recorded per coupon in the
  workbook.
- Truncates each record: marks pre-load slack and post-fracture rebound (first
  post-UTS frame where load < 50% of peak) as outside the analysis window
  (`kept = False`) rather than dropping rows. The window **starts on the rising
  edge** — the last sample *below* 2 % of peak before the peak — not the first
  sample above it, which triggers on pre-touchdown baseline noise. On P01 this
  moves the start from frame 0–12 to frame 41–84, i.e. the old rule was keeping
  40–80 frames of pre-touchdown noise, some of it landing inside the modulus
  fit window. `FlexuralDIC_Level2.truncation_mask` already worked this way;
  the two now agree, as does `TensileDIC_Level2_tmp.py`'s copy.
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
  - **Modulus E** (D638 §11.4) — slope of the linear region (0.05–0.3% strain),
    and **toe compensation** (D638 Annex A1) in the same step, because the two
    are not separable: the fit window is defined in *toe-corrected* strain but
    the correction is what the fit produces, so `fit_modulus` iterates (fit,
    shift the origin, re-select the window, refit). This is
    `FlexuralDIC_Level2.fit_modulus` — the two pipelines now compute a modulus
    the same way.

    The old code selected the window on *raw* strain and subtracted the toe
    afterwards, which contradicted the comment on `MODULUS_STRAIN_RANGE`
    ("without including the toe") and made the toe correction a mathematical
    no-op for E — subtracting a constant from the abscissa after the window is
    fixed shifts the origin without changing the slope. Confirmed empirically:
    an uncorrected E and the old "toe-corrected" E agree to the last digit on
    all 35 coupons. **The numerical effect of fixing it is small**: measured
    toe offsets are 1e-5 to 1e-4 strain, an order of magnitude below the
    window's 5e-4 floor, so iterating moves mean E by 0.16 %. The exception is
    `P01-TCL45-01` (toe 1.1e-3). Fixed for correctness, not for the number.
  - **Chord modulus** — secant across `CHORD_STRAIN_RANGE` (0.10–0.30 % strain),
    reported alongside the tangent value the way `FlexuralDIC_Level2` reports
    `Ef_*_chord_GPa`. Not a D638 requirement; it is here because the sensitivity
    sweep (below) shows the tangent modulus is not a stable quantity on this
    material. Its window is deliberately *not* `MODULUS_STRAIN_RANGE` — see the
    `LOAD_START_FRAC` warning in that section. It returns NaN rather than
    extrapolating when the record doesn't reach both endpoints, and Level 2
    prints the strain span it did reach.
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
  properties (E tangent and chord, toe strain, modulus fit points, yield
  stress/strain, UTS, strain at UTS, Poisson's ratio, plus the load-scale
  provenance — `Load Scale (N/unit)`, `Load Scale R2`, `DIC Load Coverage`, so
  a suspect stress axis is visible in the workbook and not only in a console
  log) are written into new columns on that coupon's existing
  row, matched by Specimen ID. Only those columns are touched — other
  rows, formulas, and formatting in the workbook are left alone. If the
  file is open elsewhere when Level 2 runs, this step is skipped with a
  warning rather than failing the whole run. **This is the only place
  scalar properties are stored** — Level 3 reads them back out of here
  rather than recomputing.
- `<DIC_DIR>/level2_group_stats.csv` — D638 §11.7 mean/std/count per
  (exposure, direction) group, under a `test` index level of `"tensile"`.
  **Shared with `FlexuralDIC_Level2.py`.** That index level is load-bearing:
  both test types have CL and IS exposures at 00 and 90, so without it the two
  would overwrite each other. Each script replaces only its own test's rows.
  Files written before the `test` level existed are read as tensile-only and
  upgraded in place.

#### Variant — `TensileDIC_Level2_tmp.py`

A one-off alternate Level-2 run, used when the DIC sync CSV's own "Load"
channel for a batch is only trustworthy as a *shape* and its peak-only
rescale (Level 2's normal method) isn't good enough. Instead of scaling
`load_raw`, it maps the *entire* raw MTS force/displacement time series
onto each DIC frame: the DIC-sync peak-force row and the raw-MTS
peak-force row are treated as the same physical instant (the two files
don't share a clock), and each DIC frame's time offset from that anchor
is used to interpolate MTS force/displacement at the matching MTS-clock
time. Truncation, optional smoothing, and property calculations are
otherwise identical to `TensileDIC_Level2.py`, and it writes to the exact same
files (the per-coupon CSV, `FSR-SpecimenTesting.xlsx`,
`level2_group_stats.csv`) — running it after `TensileDIC_Level2.py` overwrites
Level 2's results with this method's. It also adds one column of its own,
`disp_mm_mts` (the peak-anchored MTS displacement it derives), without
touching Level-1's own `disp_mm` column. `TensileDIC_Level2.py` itself is never
modified by running this; it remains the standard pipeline for future
batches. Needs the same raw DIC sync CSVs as Level 1 (`DATA_ROOTS`,
pointing at `DIC_DIR/raw/...`) plus `MTS_DIR`, in addition to Level-2's
usual per-coupon CSV input.

---

### Level 3 — `TensileDIC_Level3.py`

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
  excluded for what the `TODO` beside it calls an anomalously large toe
  correction, pending backup DIC data.

  **That diagnosis was wrong, and the real one is worse.** The 2026-08-29
  sensitivity run found that this coupon's DIC record covers only 90.6 % of the
  MTS peak force — the sync CSV stopped before the specimen failed. Its
  peak-anchored stress axis was inflated ~10 %, which is where the large toe
  came from. Level 2's regressed scale fixes its *modulus* (−10.8 %), but its
  **UTS and strain-at-UTS are still not the specimen's**, and no amount of
  reprocessing recovers them from this recording. `P01-TSW00-01` (96.3 %) and
  `P01-TSW00-02` (97.2 %) have the same problem more mildly and are **not**
  currently excluded.

  Recommendation: keep `P01-TCL45-01` out of UTS and strain-at-UTS statistics
  and consider excluding the other two from those as well; its modulus is now
  usable. Backup DIC data would only help if it is a *longer* recording — a
  re-export of the same frames will not change the coverage. `DIC_EXCLUDE` is
  all-or-nothing per coupon, so splitting modulus from UTS means either
  splitting that switch or noting the caveat in the write-up.

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

## Flexural pipeline — ASTM D790

Two levels so far, mirroring the tensile pipeline's split and conventions:
`FlexuralDIC_Level1.py` → `FlexuralDIC_Level2.py`. **Bare bones on purpose** —
the structure is in place and the numbers are validated, but the plotting and
reporting layer is deliberately absent so results can be added one at a time.

**One CSV per coupon, columns appended as it moves through the levels** —
exactly the tensile convention. Level 1 writes `<DIC_DIR>/<coupon_id>.csv`, one
row per DIC frame, untruncated; Level 2 reads that same file and appends its
own columns in place. Truncation sets a boolean `kept` column rather than
dropping rows. Per-coupon scalars live once each in
`<DIC_DIR>/coupon_scalars.csv` (Level 1) and in `FSR-SpecimenTesting.xlsx`
(Level 2) — the same two places the tensile pipeline puts them. Flexural and
tensile per-coupon CSVs share the same `DIC_DIR`, and the flexural rows share
`coupon_scalars.csv` and the workbook with the tensile ones, without colliding,
because the coupon IDs differ (`P01-F…` vs `P01-T…`). There are no
`flexural_*`-prefixed files any more; the two test types use one set of
files throughout.

**Why this isn't the tensile pipeline with different constants.** The tensile
pipeline reduces each frame to two virtual point extensometers, which is the
wrong reduction for a bend test. The flexural ROI is the *side profile* — the
0.50 in specimen depth seen edge-on over the full span — so each frame carries
the entire bending kinematic field. Each frame is reduced instead to midspan
deflection (referenced to the chord through the two supports, which removes
rigid settling and fixture tilt), curvature and neutral-axis height from a
straight-line fit of `exx` against `Y` at midspan, and the extreme-fibre strains
that fit extrapolates to the two faces. The ROI is inset from both faces by the
correlation subset radius — 7.68 of 12.50 mm on `P01-FCL00-01`, 61 % — so the
surface strain D790 asks for is never directly measured.

**Three strain channels are carried side by side, and they are meant to
disagree**: `eps_curvature` (κ·d/2, DIC only, no span in it), `eps_deflection`
(6Dd/L², D790 Eq. 4, needs the span), `eps_crosshead` (same from MTS travel, so
it also carries machine compliance and support indentation). They share one
force and one specimen, so every difference between the curves is measurement
method, not material. `eps_curvature` is the reference; both gaps are printed.
On `P01-FCL00-01` the deflection modulus runs 10.9 % below the curvature
modulus and the crosshead one 12.8 % below.

### Span — 8.00 in, confirmed

`FLEX_SPAN_MM = 203.2` (8.00 in, D790's 16:1 for the nominal 0.50 in depth).
**Confirmed against the fixture.** This was an open assumption for a long time —
no record of the setting exists in the repo or the workbook — so older notes,
and the removed scratch scripts, hedge about it. They are out of date.

Two related numbers, not the same thing, both correct:

| | | |
|---|---|---|
| **203.2 mm** | the fixture setting — where the support rollers are | D790 defines flexural stress and strain on this nominal span, so it is what every formula here uses |
| **~200 mm** | where the bending moment actually crosses zero at low load | measured off the DIC curvature diagram; shrinks further under load |

The second is measurable because a roller cannot transmit moment, so each limb
of the curvature-vs-position diagram runs straight down to zero exactly at the
support, and the two zero crossings locate them. On `P01-FCL00-01` that gave
200.0 ± 1.2 mm extrapolated to zero load, falling to 190.9 mm by peak as the
contact points ride inboard over the rollers (−5.96 mm/kN).

With the fixture now confirmed at 203.2 mm, **that 1.6 % gap is contact
geometry, not a fixture error** — the specimen touches each roller slightly
inboard of its centre, and rolls further inboard as it bends. Two consequences
worth knowing:

- It is *not* an error to correct. D790 defines σ and ε on the nominal span,
  which is what the pipeline uses.
- It does mean the `eps_deflection` − `eps_curvature` gap (−7.6 to −10.9 %
  across the batch) can **no longer be partly attributed to a wrong span**, the
  way the older notes did. The effective span being ~1.6 % short accounts for
  some of it; the rest is D790's small-deflection kinematics, shear, and
  indentation under the loading nose. `eps_curvature` has no span in it at all
  and stays the reference regardless.

The measurement itself is not yet ported into this pipeline — see the list
below. The constant is defined separately in both levels and the two must match.

---

### Level 1 — `FlexuralDIC_Level1.py`

**Step A** — dumps each `.out` to a CSV next to it, the way tensile Step A
does, with the identical `EXPORT_VARS` column list, so a flexural per-frame CSV
and a tensile one are interchangeable. **On** (`DO_EXPORT_FRAMES = True`), to
match tensile. Step B does not read these — it loads each `.out` itself — but
`dic_heatmap_animation.py` does, and reading a CSV needs no vicpyx. Budget
~1 GB per coupon.

**Step B** — reduces every `.out` frame to the bending kinematics, pairs them
with the MTS record, and writes the per-coupon CSV. The slow step (~50–60 s per
coupon, one vicpyx load per frame).

**Step C** — reports which variables each coupon's `.out` files actually
contain. Cheap, vicpyx-only, and the reason it exists: **run it after
reprocessing a project in VIC-3D** to confirm new inspector items are reaching
the export before spending an hour on Step B. It prints the full variable list,
flags any `REQUIRED_VARS` that are absent, and lists variables that are present
but not being read — anything in that last group can be added to
`OPTIONAL_VARS` and it will be carried straight through Step B.

**Variable handling, built for the reprocess.** `REQUIRED_VARS`
(`sigma, X, Y, V, exx`) must be present or a frame is unusable. `OPTIONAL_VARS`
are read when they exist and silently skipped when they don't, so adding a new
inspector item there costs nothing on coupons that predate it. As of now every
flexural `.out` carries the standard 21-variable VIC-3D set and no inspector
items.

**Alignment — the last correlated frame.** With no usable sync load channel
(see the warning at the top), the two clocks cannot be aligned by peak
anchoring the way `TensileDIC_Level2_tmp.py` does, and they cannot be aligned
by scanning either: force is very nearly linear in time in both records, and a
time shift between two straight lines is absorbed exactly by a regression's
intercept, so the residual has no minimum to find. Displacement is worse — a
constant-rate ramp is straighter still.

Fracture is the one sharp event in the test, and on a bend specimen it can be
read straight off the images: correlation is lost the instant the specimen
cracks — *every* point in the frame goes invalid at once, not gradually — so the
last correlated frame **is** the fracture instant, measured with no analog
channel in it. That frame is aligned to the MTS force peak. The crosshead-travel
cross-check is reported alongside: the sync CSV's displacement input (`ai1`, the
one channel genuinely connected) and the MTS crosshead record measure the same
motion by separate paths, so after alignment they should agree — 0.07–0.28 mm
RMS over 9–14 mm of travel on the coupons run so far. That does not
*independently* confirm the offset (a constant-rate ramp cannot), but a large
residual would be a loud signal that the two files are not the same test.

**Load-cell tare.** The flexural records open on a constant ~860 N held through
the approach travel — the loading nose hanging on an un-tared cell. It shows on
all 12 flexural specimens in both orientations, and left in it inflates flexural
strength ~1.6×. `find_force_baseline()` detects and removes it (ported from
`mts_plots.py`). Tensile and bearing records don't have it.

**ROI orientation check — three coupons are blocked on it.** Every reduction in
this script assumes the VIC-3D world frame is oriented the same way on every
coupon: `X` along the specimen (the span), `Y` through the depth. That is not
automatic — it comes from the calibration and the alignment plane chosen when
the project was built — and **three coupons in this batch were reconstructed on
a different frame**:

| coupon | ROI X extent | ROI Y extent | vs. depth `d` |
|---|---|---|---|
| `P01-FCL90-03` | 89 mm | **196 mm** | 15.2× — span and depth swapped |
| `P01-FIS90-01` | 157 mm | **309 mm** | 24.0× — span and depth swapped |
| `P01-FIS00-01` | 222 mm | **40 mm** | 3.2× — X is right, Y is not the depth |

The other nine sit at 7.7 mm of a ~12.5 mm depth, as they should.

Nothing downstream can detect this on its own, which is why it is caught here.
Left unguarded, the curvature fit happily regresses `exx` against a coordinate
that isn't depth and returns R² ≈ 0 with a neutral axis tens of mm off the
specimen (+56, +140 mm on two of these), and Level 2 turns that into a
plausible-looking row of numbers — `P01-FIS00-01` produced a 10.02 GPa modulus
against the batch's 7.9–8.2 GPa before the check existed. So Level 1 tests the
ROI's Y extent against the specimen depth (`ROI_DEPTH_FRAC_RANGE`) and its
aspect ratio (`ROI_MIN_ASPECT`) and **skips the coupon, writing nothing** — a
missing coupon is recoverable, a silently wrong one is not. The check runs
before anything expensive, so a bad coupon costs ~1 s.

**Fix at the source**: re-align those three VIC-3D projects' coordinate systems
and re-export, then re-run Level 1. Worth doing as part of the inspector-item
reprocess.

**Fixture location.** Midspan is found as the interior extremum of the deflected
shape — the one place along a 3-point beam where `dV/dx = 0`, by symmetry — then
refined twice as the largest departure from the support chord. This needs no
prior knowledge of the sign of `V` or of where the specimen sits in the world
frame. Set `MIDSPAN_X_MM` to pin it manually instead; leaving it `None` (the
default) is recommended, since the located position is itself a fixture check.
The deflection sampling windows sit a few mm inboard of the supports whenever
the ROI doesn't reach that far, where the true deflection isn't zero, so the
midspan value is scaled back up using the ideal 3-point shape
(`APPLY_SUPPORT_OFFSET_CORRECTION`, ~1.6 % here).

**Switches**
- `DO_LIST_VARS` — Step C. Cheap; leave on.
- `DO_EXPORT_FRAMES` / `OVERWRITE_FRAMES` — Step A. On by default, as tensile.
  `OVERWRITE_FRAMES` defaults `False`: skip `.out` files whose `.csv` exists.
- `DO_BUILD_L1` / `OVERWRITE_L1` — Step B. `OVERWRITE_L1` defaults `False`:
  skip coupons whose per-coupon CSV already exists.
- `PRINTS` / `EXPOSURES` / `DIRECTIONS` / `REPLICATES` — coupon selection.
  `EXPOSURES` is `{CL, IS}` and `DIRECTIONS` is `{00, 90}`; the other tensile
  exposures and the 45° direction weren't bend-tested.

**Inputs**
- `<coupon_dir>/*.out` — VIC-3D full-field export, one per DIC frame.
- `<coupon_dir>/<folder>.csv` — VIC sync CSV. **Read for `Time_0_0` only** (plus
  the displacement channel for the alignment cross-check). None of its analog
  channels reach the per-coupon CSV.
- `<MTS_DIR>/<coupon_id>.txt` — MTS raw: `disp_mm, force_N, output_V, time_s`.
- `FSR-SpecimenTesting.xlsx` — depth `d` and width `b`. Read only, never written.
- The raw folders **drop the print prefix and the dashes**: `P01-FCL00-01` lives
  in `DIC/raw/2026_FSR_Flexural_FCL_FIS/FCL0001`, so unlike tensile the coupon
  ID cannot be used to find the directory directly. `raw_folder()` maps it.

**Outputs**
- `<DIC_DIR>/<coupon_id>.csv` — full per-frame record: `step, time_s, force_N,
  disp_mts_mm, n_pts, defl_mm, kappa_1pmm, na_Y_mm, profile_r2, eps_bot,
  eps_top, eps_membrane`. Level 2 appends to this same file.
- `<DIC_DIR>/coupon_scalars.csv` — one row per coupon: `b_mm`, `d_mm`, the
  located fixture (`x_mid_mm`, `x_left_mm`, `x_right_mm`, `defl_sign`), the ROI
  and assumed face positions, `tare_mts_N`, `peak_net_N`, `frame_rate_hz`,
  `mts_offset_s`, `break_frame`, `disp_check_rmse_mm`. The only place these
  per-coupon scalars are stored. **The same file `TensileDIC_Level1.py`
  writes** — the merge is keyed on coupon ID and the two test types never share
  one, so each script leaves the other's rows untouched and pandas fills the
  columns a given test type doesn't have. Upserted, so coupons skipped on a run
  keep their existing row.
- `<coupon_dir>/<out_filename>.csv` — Step A, one per `.out`.

---

### Level 2 — `FlexuralDIC_Level2.py`

Reads Level-1's per-coupon CSV and `coupon_scalars.csv`, applies failure
truncation, and computes the D790 properties. Cheap (pure pandas/numpy over an
already-built CSV), so it always recomputes — re-run freely while tuning
truncation and modulus windows.

**What it does**
- Truncates each record: marks pre-load slack and post-peak decay outside the
  analysis window (`kept = False`) rather than dropping rows. Two differences
  from the tensile version: the window additionally requires **DIC validity**,
  since a bend specimen loses correlation at fracture several frames before the
  load channel finishes falling; and the start is taken on the *rising edge*
  (last sample below threshold before the peak) rather than the first sample
  over it, because with the tare removed the pre-touchdown baseline sits on zero
  with noise either side and a plain threshold can trigger tens of seconds early.
- **No smoothing pass.** D790 doesn't call for filtering the stress-strain
  record, the MTS force is the primary channel and is already clean (~0.7 % of
  span), and the DIC-derived channels are the cleanest signals in the test.
- Computes:
  - **Flexural stress** (D790 §12.2 Eq. 3) — `σ = 3PL / 2bd²`, on the
    tare-removed force.
  - **Three flexural strains** — `eps_curvature`, `eps_deflection`,
    `eps_crosshead`, each toe-compensated on its own.
  - **Toe compensation** (D790 §12.1 / D638 Annex A1) — fitted **twice**,
    because the fit window is specified in corrected strain but the correction
    is what the fit produces. Without the second pass the tangent and chord
    moduli are quietly measured over two different strain ranges, offset by the
    toe, and are not comparable.
  - **Tangent modulus** (D790 §12.4) and **chord modulus** (§12.5) on each of
    the three channels, over 0.05–0.25 % strain.
  - **Curvature modulus** — `E = (M/I)/κ`, fitted in moment-curvature space.
    Not in D790; the DIC-native equivalent of `E_B`, and a check on the
    arithmetic of `eps_curvature` (the two agree to 0.6 %).
  - **Flexural strength** (D790 §3.2.7) — max flexural stress, with a warning if
    the specimen did *not* break before 5 % strain (in which case D790 §12.2
    asks for the stress at 5 % instead), and another if `D/L > 0.10` at peak,
    where §12.3's large-deflection correction is called for.
  - **Bending-quality diagnostics** — R² of the through-depth `exx(Y)` fit
    (summarised only above 25 % load, since at low load the strain range is
    below the DIC noise floor), the neutral-axis offset, and the membrane
    strain.
  - **Group stats** — mean/std/count per exposure × direction.

**The neutral-axis offset has two causes, and they're reported separately.** A
mis-centred ROI is a constant, present from the very first frame: an ROI offset
by `δ` reads as a neutral axis `δ` off mid-depth and a membrane strain of
exactly `κδ`, indistinguishable from a genuinely asymmetric specimen.
Asymmetric material nonlinearity only develops under load. So the low-load
baseline (`na_offset_at_low_load_mm`) and the growth above it (`na_drift_mm`)
are stored separately — **only the drift is a statement about the material.**
Neither touches `κ`, and therefore neither touches `eps_curvature` or any
modulus.

**Inputs**
- `<DIC_DIR>/<coupon_id>.csv` — Level-1 output.
- `<DIC_DIR>/coupon_scalars.csv` — Level-1's per-coupon scalars.

**Outputs**
- `<DIC_DIR>/<coupon_id>.csv` — Level-1's columns, unchanged, plus `kept`,
  `stress_MPa`, `M_over_I_MPa_per_mm`, `eps_curvature`, `eps_deflection`,
  `eps_crosshead` (all toe-corrected). All `NaN` where `kept` is `False`.
- `FSR-SpecimenTesting.xlsx` — the D790 scalars written into each coupon's row:
  the four moduli (tangent + chord), toe offsets, `sigma_fM_MPa`,
  `sigma_fB_MPa`, strain at max on each channel, and the bending-quality
  diagnostics. Every header is prefixed `Flex ` so nothing is confusable with
  the tensile columns in the same sheet; see `SPECIMEN_SHEET_COLUMNS`. The
  workbook is the single source of truth for per-coupon scalars, exactly as it
  is for tensile. Level-1 measurements (`b_mm`, `d_mm`, `tare_mts_N`,
  `break_frame`, …) are deliberately *not* repeated here — they are in
  `coupon_scalars.csv`, and a second copy could disagree with the first.
- `<DIC_DIR>/level2_group_stats.csv` — mean/std/count per (exposure, direction),
  under a `test` index level of `"flexural"`. **Shared with
  `TensileDIC_Level2.py`.** That index level is load-bearing: both test types
  have CL and IS exposures at 00 and 90, so without it a flexural `(CL, 00)` row
  would overwrite the tensile one. Each script replaces only its own test's rows
  and reads the other's back unchanged.

### Level 3 — not written yet

There is no `FlexuralDIC_Level3.py`. Plot from `<DIC_DIR>/<coupon_id>.csv`
(`kept == True` rows) and the `Flex …` columns of `FSR-SpecimenTesting.xlsx`
until there is — the same two sources `TensileDIC_Level3.py` reads.

### Validation

Levels 1–2 reproduce, exactly, the numbers from the earlier standalone
`Flexural_tmp.py` worked example on `P01-FCL00-01`:

| | |
|---|---|
| tare removed / net peak | 860 N / 1526 N |
| break frame (correlation lost) | 847 of 881 |
| MTS clock leads frames by | 25.81 s |
| analysis window | frames 141–847, 707 of 881 kept |
| flexural strength `σ_fM` | 117.2 MPa at 2.61 % strain |
| `E_f` curvature / deflection / crosshead / (M/I)-vs-κ | 7.94 / 7.08 / 6.93 / 7.99 GPa |
| through-depth fit R² | 0.9992 median above 25 % load |

**All 12 coupons have been run.** Nine processed cleanly; the three listed under
the ROI orientation check above are skipped pending a VIC-3D re-alignment. The
nine group tightly by orientation, which is the real check that the method
travels across coupons — mean ± std:

| | n | `σ_fM` (MPa) | `E_f` curvature (GPa) | `E_f` deflection | `E_f` crosshead | strain at max |
|---|---|---|---|---|---|---|
| **0°** | 5 | 115.3 ± 3.7 | 8.08 ± 0.11 | 7.35 ± 0.16 | 7.01 ± 0.12 | 2.6 % |
| **90°** | 4 | 52.7 ± 3.9 | 4.01 ± 0.09 | 3.67 ± 0.07 | 3.47 ± 0.05 | 1.5 % |

The curvature modulus scatters by ~1.3 % within an orientation across two
different exposures, and the channel gaps stay consistent coupon to coupon
(deflection −7.6 to −10.9 %, crosshead −11.7 to −15.1 %) — they are systematic,
as the argument above says they should be, not noise. Through-depth fit R² is
0.998–0.9996 on all nine. The 90° modulus lines up with the 3.43–3.57 GPa the
MTS ramp slopes and the tensile workbook give for that orientation.

### Still to do

- **Re-align the three blocked VIC-3D projects' coordinate systems**
  (`FCL9003`, `FIS0001`, `FIS9001`) so `X` is the span and `Y` the depth, then
  re-run Level 1 for them. Worth folding into the inspector-item reprocess.
- **Reprocess the VIC-3D projects with more inspector items**, then re-check
  Step C's variable listing and add whatever appears to `OPTIONAL_VARS`.
- **Port the span measurement** off the curvature diagram (see the span section
  above). No longer needed to validate the fixture — that is confirmed — but it
  measures where the moment actually crosses zero, which is a real diagnostic of
  contact behaviour and feeds the interpretation of the channel gaps.
- **Decide where the D790 scalars live** — their own columns in
  `FSR-SpecimenTesting.xlsx` alongside the tensile ones, or the properties CSV
  as now.
- **Write Level 3** (per-coupon and group plots, stat tables), one result at a
  time.

---

## Result — how much of the tensile modulus scatter is processing?

`tensile_modulus_sensitivity.py` turns every knob in the Level-2 reduction, one
at a time and then all together, and recomputes E for each setting: 810
settings × 35 coupons. It is **read-only on the pipeline** — it writes only to
`<ROOT>/figs/sensitivity/` and never touches `FSR-SpecimenTesting.xlsx` or any
Level-1/2 output. Run it whenever a processing choice is in dispute; it turns
the argument into a number.

E is a *slope*, and that alone decides which factors can matter. A scale error
on the stress axis multiplies E directly. A DC offset on the stress axis cannot
touch it — that is absorbed by the toe correction. Smoothing a nearly-straight
segment is nearly the identity. The sweep confirms all three.

### Ranked: how far E moves when each factor alone is flipped

| factor | max \|ΔE\| | median | p90 |
|---|---|---|---|
| **window** — modulus fit window | 12.28 % | 2.73 % | 6.52 % |
| **scale** — how `load_raw` becomes N | 9.37 % | 2.86 % | 7.14 % |
| **toe** — Annex A1 handling | 7.07 % | 0.00 % | 2.58 % |
| **trunc** — where the window starts | 5.99 % | 0.00 % | 2.63 % |
| **smooth** — none / median / butter | 5.66 % | 0.77 % | 2.56 % |
| **area_fac** — ±2 % on A₀ | 2.04 % | 2.00 % | 2.04 % |

`area_fac` returning exactly ∓2 % for a ±2 % perturbation is the harness
checking itself: E ∝ 1/A₀ exactly, so anything else would mean a bug in the
sweep. Set `AREA_PERTURBATIONS` from your own caliper repeat spread and that row
becomes your measurement uncertainty, on the same figure and at the same scale
as the processing choices.

### Four things it settled

**1. The load scale was a real first-order error, and it is fixed.** CV across
coupons 2.55 % → 0.71 % (details in Level 2 above). Worth doing on its own
terms: a hardware constant should not vary 13 % coupon to coupon.

**2. Three coupons have truncated DIC records.** The finding with the largest
consequences, and the one nothing in the pipeline previously looked for — see
the warning at the top of this README. `P01-TCL45-01` is anomalous because its
DIC record covers 90.6 % of the test, **not** because of its toe correction.
The TODO beside `DIC_EXCLUDE` in `TensileDIC_Level3.py` now has an answer, and
it is not the one it expected: the regressed scale fixes that coupon's modulus
(−10.8 %), but its UTS and strain-at-UTS remain unrecoverable from this
recording, so it should stay out of any group statistic involving them.

**3. The toe fix was a correctness fix, not a numerical one.** 0.16 % on mean
E. See Level 2 above.

**4. The recommended settings do NOT reduce the group CV.** Mean within-group
CV% of E goes **3.83 → 3.87** — it does not improve. This is the opposite of
what such an exercise is usually set up to show, and it is reported here
because it is what the data says. The processing choices being argued over are
not what drives the coupon-to-coupon scatter. They are still the right choices
— a per-coupon scale error of up to 10 % is wrong whether or not it happens to
cancel in a CV over three replicates — but **they are not a scatter fix and
should not be presented as one.**

### What actually drives the scatter: there is no linear region

| | |
|---|---|
| E spread across the five-window grid | median **7.7 %** per coupon, max 21.8 % |
| modulus fit R² | median **0.879**, min 0.808 (worst on 45°/90°) |
| tangent (0.05–0.30 %) vs chord (0.10–0.30 %) | median **+6.1 %**, max 33.6 % |
| within-group CV of E, for comparison | **3.9 %** |

A quantity that moves more when you change the fit window than when you change
which coupon you measured is not a tangent modulus. This is why Level 2 now
writes `E_chord_GPa` next to `E_GPa` — but **the chord is a cross-check, not a
replacement**, and it would be wrong to report it as the better number. A chord
is read from two interpolated points, so where the tangent averages strain
noise over ~60 points the chord takes it at face value. Measured, by direction:

| | E tangent (GPa) | E chord (GPa) | UTS (MPa) |
|---|---|---|---|
| 00° | 6.931 ± 0.124 | 6.321 ± **0.477** | 76.30 ± 2.25 |
| 45° | 3.187 ± 0.153 | 3.081 ± **0.338** | 48.18 ± 2.10 |
| 90° | 3.437 ± 0.207 | 3.301 ± **0.299** | 39.65 ± 1.83 |

The chord's standard deviation is 1.4–3.8× the tangent's. It is more
*transparent* — "the average stiffness between these two strains" survives
being asked what window you used — but it is a noisier estimate.

**Recommendation**, in order:

1. **Whatever you report, state the window.** This is the whole finding.
2. **Quote the window sensitivity as the error bar on E** — ~8 % median, not
   the 3.9 % within-group CV, which understates it by a factor of two. That is
   the honest uncertainty on a modulus from this data.
3. **Use the tangent/chord pair as a straightness test**, the way
   `FlexuralDIC_Level2` does. Agreement means the segment really is straight;
   the 5 % median (14 % max) disagreement here means it is not.
4. **Don't switch windows hoping to fix it.** The code's
   `MODULUS_STRAIN_RANGE = (0.0005, 0.003)` and `solid_mechanics_core.pdf`
   §8.1's 0.05–0.25 % disagree, and that has been left alone deliberately: the
   sweep shows *every* window in the grid is equally arbitrary on this
   material. Switching improves nothing.
5. **The real fix is less strain noise** — see the next section.

⚠ **`LOAD_START_FRAC` and `MODULUS_STRAIN_RANGE` are inconsistent, and this is
an open decision.** The analysis window starts at the last frame *below* 2 % of
peak load, and the next frame — the first one kept — can already be well past
it: at 10 Hz on a ~50 s ramp one frame is ~2 % of the ramp, so the first kept
frame lands anywhere between 2 % and 8.5 % of peak depending on phase. On **18
of the 35 P01 coupons** that puts the lowest available corrected strain between
5.6e-4 and 9.5e-4 — *above* the modulus window's 5e-4 floor. The tangent fit
tolerates it (it uses whatever points fall inside the window, which is itself a
per-coupon window inconsistency of exactly the kind the toe fix removed); a
chord cannot, so `CHORD_STRAIN_RANGE` is `(0.001, 0.003)` — 1e-3 being the
lowest round floor every coupon actually reaches. Either lower
`LOAD_START_FRAC` or raise the modulus window's floor to make the two agree;
`chord_modulus` refuses rather than extrapolating, so this fails loudly instead
of silently.

### The next thing to measure, not yet done

A low R² on a σ-ε line means noise on the **strain** axis — which also biases a
least-squares slope *downward* (errors-in-variables attenuation), so it is not
only scatter. `DO_GAUGE_SWEEP = True` re-derives axial strain from the
per-frame full-field CSVs at several gauge lengths by both the two-point method
Level 1 uses and a least-squares fit of `V` against `Y` over the whole ROI. The
field fit is ~√N quieter and its R² says whether strain was uniform across the
gauge — which for the 45° coupons is the measurement, not a diagnostic.
(`FlexuralDIC_Level1` already works this way, fitting `exx` against Y for κ.)

It is off by default because it re-reads every per-frame CSV (17–20 k rows
each, on a network share). It has not been run. Replacing the extensometer
changes what the strain channel *is*, not how well it is computed, so it wants
the noise-floor numbers in hand before it is made.

### Outputs — `<ROOT>/figs/sensitivity/`

| file | what |
|---|---|
| `modulus_sensitivity_full.csv` | every combination, one row each |
| `modulus_sensitivity_oat.csv` | one-factor-at-a-time deltas (the tornado data) |
| `modulus_sensitivity_groups.csv` | CV% of E by exposure × direction, old vs new |
| `scale_check.csv` | peak vs regressed scale per coupon, with R², intercept, lag, DIC coverage |
| `window_check.csv` | E vs fit window, tangent vs chord, fit R² |
| `tornado_<coupon>.png` | per-coupon tornado |
| `window_sensitivity.png` | E vs fit window, all coupons |
| `gauge_sensitivity.png` | E vs gauge length (only if `DO_GAUGE_SWEEP`) |

Full per-change detail — including what in the original review was checked and
*not* applied, and why — is in `CHANGELOG_tensile_2026-08-29.md`.

---

## Outside the pipelines

`mts_plots.py` (renamed from `mts_quick_plots.py`) is a standalone first-look
script — MTS channels only, no DIC, no property extraction. It plots
force–displacement and stress–strain for all three test types (one figure each,
colour by exposure, linestyle by orientation) straight from the raw `.txt` files
plus the geometry columns in the specimen sheet. Use it to eyeball a batch; use
the levelled pipelines for anything quotable.

- **Geometry input** — it reads `FSR-SpecimenTesting.csv`, the CSV export kept
  alongside `FSR-SpecimenTesting.xlsx` under the same stem, and falls back to
  the `.xlsx` if the CSV is missing or unreadable. The CSV needs no Excel engine
  and is not locked while the workbook is open. It is a Windows export, so it is
  cp1252, not UTF-8 — the reader tries `utf-8-sig`, `cp1252`, `latin-1` in turn.
- **ASTM reductions** — every formula it applies is written out in full, with
  its source, in the `ASTM EQUATIONS` block at the top of the file: D638 §11.2
  σ = P/A₀ and Annex A1 toe compensation; D790 §12.2 Eq.3 σ_f = 3PL/2bd²,
  §12.3 Eq.4 ε_f = 6Dd/L², and the §12.3 large-support-span stress correction
  (`--large-deflection`, off by default — this fixture is 16:1 and D/L peaks at
  0.067, below the 0.10 threshold, and it prints max D/L on every run);
  D953 §13.2 σ_b = P/(D·t) on the projected area with the §13.3 4 % hole
  deformation ordinate drawn on the panel.
- **Toe compensation** follows D638 Annex A1.3 literally — the *steepest*
  straight segment of the load–deflection record, searched between 5 % and 60 %
  of peak load, extrapolated to zero load. That agrees with the fixed 10–40 %
  window it used before to within 0.02–0.12 mm on every P01 coupon.
- **The tensile right-hand axis is not a D638 strain** and is labelled as such.
  D638 §3.2.5 nominal strain is referred to the original *grip separation*,
  which was never recorded for this batch, so crosshead travel is divided by the
  DIC axial gauge length (4.36 in) instead, purely so the axis is on the same
  scale as the DIC record. Fill in `TENSILE_GRIP_SEPARATION_MM` if it is ever
  measured and the panel becomes a true D638 nominal strain.
- It is where the flexural load-cell tare was first documented, and it detects
  and removes the tare itself (`find_force_baseline`, with `--keep-tare` to
  disable). Tare removal runs *before* the toe fit, since a constant offset
  moves the Hookean line's zero-load intercept.
- Outputs are `figs/mts_tensile.png`, `figs/mts_flexural.png`,
  `figs/mts_bearing.png` (was `mts_quick_*.png`), plus a per-specimen stdout
  table carrying the ASTM validity flags: D790 §12.2's 5 % strain rule, D790
  §12.3's D/L threshold, and D953 §13.3's 4 % deformation.

`matlab/tensile_plots.m` and `matlab/bearing_plots.m` are an independent MATLAB
re-implementation of part of tensile Level 3's plotting, kept for reference;
they are not run as part of either pipeline and aren't guaranteed to stay in
sync.

---
## ⚠ Specimen geometry comes from the CSV, not the workbook

`Width / Dia. (in)` and `Computed Area (in²)` in `FSR-SpecimenTesting.xlsx` are
**formula** columns:

```text
=IF(C2="Tensile",1.5,IF(C2="Bearing",0.5625,IF(C2="Flexural",1,"")))
```

openpyxl does not evaluate formulas. Every time a Level 2 writes its scalars
back into the workbook it saves the formula and **drops the cached value**.
`write_specimen_sheet` sets `wb.calculation.fullCalcOnLoad = True` to ask for a
recalculation — but only *Excel* honours that flag. `pandas.read_excel` and
openpyxl both read the cache, so those columns come back **blank** until
somebody opens the workbook by hand and saves it.

That is a read/write cycle between the two levels quietly poisoning the
geometry Level 1 depends on. When it bites, the failure is silent and
misleading: Level 1 computes `area = thickness × NaN`, writes `area_mm2 = NaN`
to `coupon_scalars.csv`, and Level 2 reports **`insufficient data` for every
coupon** — because the whole stress axis is NaN and `compute_properties` never
gets 10 finite points.

**Both Level 1 scripts therefore read `FSR-SpecimenTesting.csv`** — the export
kept alongside the workbook under the same stem — and fall back to the `.xlsx`
only if the CSV is missing or unreadable. The CSV holds evaluated values, needs
no Excel engine, and is not locked while the workbook is open. It is a Windows
export, so it is cp1252, not UTF-8; all three readers try `utf-8-sig`,
`cp1252`, `latin-1` in turn. `mts_plots.py` already worked this way; the same
reader now lives in `TensileDIC_Level1.py` and `FlexuralDIC_Level1.py`,
duplicated rather than imported so each stays runnable alone (the convention
`FLEX_SPAN_MM` already follows).

Both now **refuse** rather than proceeding with NaN geometry, naming the cause
and the fix, and Level 1 prints which source it read.

**If you edit geometry in Excel, re-export the CSV.** It is the source of truth
for the pipeline now; the workbook is the source of truth for you. Any Level-2
run re-strips the workbook's formula caches, so the `.xlsx` fallback should not
be relied on.

Unaffected: `TensileDIC_Level3.py` reads only plain-value columns (thickness,
and the scalars Level 2 wrote), and `TensileDIC_Level2_tmp.py` takes area from
`coupon_scalars.csv`.

---

## Common configuration

All scripts share a `SWITCHES` block at the top (`PRINTS`, `EXPOSURES`,
`DIRECTIONS`, `REPLICATES`) used to select which coupons to process, and
a `PATHS` block — trimmed to only what each script actually touches.
`TensileDIC_Level3.py` doesn't need the raw-data `DATA_ROOTS` at all, since by
the time it runs everything it needs is already in `DIC_DIR` or the specimen
sheet. `TensileDIC_Level2.py` **now does** need `MTS_DIR` (it regresses the load
scale against the raw MTS record) and `DATA_ROOTS` (as a fallback clock source,
for per-coupon CSVs written before the `time_s` fix — see the sync-CSV warning
near the top). `TensileDIC_Level2_tmp.py` and `tensile_modulus_sensitivity.py`
need both for the same reasons.

The flexural scripts follow the same shape. Their `EXPOSURES` is `{CL, IS}` and
`DIRECTIONS` is `{00, 90}` — the other exposures and the 45° direction weren't
bend-tested. `FlexuralDIC_Level2.py` needs only `DIC_DIR`; Level 1 needs
`RAW_ROOT`, `MTS_DIR` and the specimen sheet.

**`FLEX_SPAN_MM` is defined in both flexural levels and the two must match.**
It is duplicated rather than shared because neither script imports the other —
keeping each runnable on its own, the way the tensile levels are. If you change
it, change it in both.

## File naming history

**The tensile scripts were renamed** `DIC_Level1/2/3.py` →
`TensileDIC_Level1/2/3.py` (and `DIC_Level2_tmp.py` →
`TensileDIC_Level2_tmp.py`) when the flexural pipeline was added, so the two
pipelines' file names say which test type they handle. Nothing about their
behaviour changed. Old notes and commit messages still use the short names.

The standalone flexural scratch scripts — `Flexural_tmp.py` (a self-contained
worked example of the whole bend-test chain on one coupon),
`Flexural_figs_tmp.py` (its figures, signal report and HTML method doc) and
`export_fields_for_matlab.py` (a one-off `.out` → CSV export for MATLAB) — are
**no longer in this directory**, and neither are the `flexural_walkthrough`
`.md`/`.pdf`/`.m` notes that went with them. `FlexuralDIC_Level1/2.py` supersede
the analysis and reproduce its numbers exactly (see Validation above). Anything
left under `<DIC_DIR>/tmp_flexural/` is output from those scripts and is stale;
the live flexural outputs are in `<DIC_DIR>` proper. If you still have the
walkthrough notes elsewhere they remain the best long-form derivation of the
method, but they describe `Flexural_tmp.py`'s file layout, not this one's.

If you're looking at old notes or file listings: the consolidated
per-coupon CSV used to be called `_L2.csv` (written by a separate
`TensileDIC_Level2.py` that only did the MTS/extensometer pairing), and the
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
