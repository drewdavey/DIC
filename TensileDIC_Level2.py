#!/usr/bin/env python3
"""
DIC_Level2.py  —  FSR Tensile Coupons
======================================
Reads Level-1's per-coupon frame CSV, scales raw load to force/stress, applies
failure truncation, and computes ASTM D638 mechanical properties. An
optional smoothing pass (median or Butterworth, see APPLY_SMOOTHING /
FILTER_METHOD) is available but off by default — ASTM D638 doesn't call
for filtering the stress-strain record; only enable it if a batch's raw
signal is genuinely too noisy. No plotting here — see DIC_Level3.py.

See also DIC_Level2_tmp.py: a one-off variant of this script that derives
force/stress from peak-force-anchored raw MTS data instead of scaling the
DIC sync CSV's own load channel, used for batches where that channel is
unreliable. It writes the same outputs as this script but leaves this
file untouched.

INPUT NOTE — the load scale
  load_raw becomes force_N by regressing the raw MTS force against it over the
  rising ramp (LOAD_SCALE_MODE = "regress"), not by the ratio of the two peak
  samples the earlier versions of this script used. See LOAD_SCALE_MODE for the
  measurement that motivated the change, and tensile_modulus_sensitivity.py for
  the full sweep. The old behaviour is still reachable with
  LOAD_SCALE_MODE = "peak".

Standards compliance — what each calculation cites
  Toe compensation     : D638 Annex A1 (mandatory unless toe is real material response)
  Modulus              : D638 §11.4   (slope of initial linear region of σ-ε)
  Chord modulus        : not a D638 requirement — reported alongside the tangent
                         value because the tangent is window-sensitive on this
                         material (see chord_modulus)
  0.2 % offset yield    : D638 §A2.6 / Fig. A2.1 (offset from toe-corrected origin)
  Tensile strength UTS  : D638 §11.2   (max stress / original area)
  Poisson's ratio       : D638 Annex A3.10.1.3 (chord at ε_a=0.002 over 0.0005-0.0025)
  Group statistics      : D638 §11.7 / §12.1   (mean, std per series)

PROCESSING NOTE
  Level-1 writes the full, untruncated per-frame record (see its docstring).
  Level-2 reads that SAME file and appends its own columns to it in place —
  no separate per-coupon output file. Failure truncation happens here via
  truncation_mask(): pre-load slack and post-fracture rebound (past 50%
  post-UTS load drop) are marked out of the analysis window rather than
  dropped from the file — every L2-derived column is NaN outside that
  window, and the boolean 'kept' column marks which rows are inside it.
  Properties are computed from this window (smoothed, if APPLY_SMOOTHING is
  on); pre-smoothing values (still restricted to the window) are kept
  alongside either way for Level-3's diagnostic overlay — they're identical
  to the smoothed columns when smoothing is off.

INPUT per coupon
  <DIC_DIR>/<coupon_id>.csv      Level-1's per-frame record (see its docstring)
  <DIC_DIR>/coupon_scalars.csv   per-coupon mts_peak_N, area_mm2 (Level-1 output)

OUTPUT per coupon
  <DIC_DIR>/<coupon_id>.csv   Level-1's columns, unchanged, plus:
                               kept, force_N, stress_MPa, strain_axial,
                               strain_transverse (toe-corrected, and smoothed
                               if APPLY_SMOOTHING is on), stress_MPa_unsmoothed,
                               strain_axial_unsmoothed (same window, pre-smoothing).
                               All L2 columns are NaN where kept is False.
  FSR-SpecimenTesting.csv     scalar properties written into each coupon's
                               row (E tangent and chord, toe strain, modulus
                               fit points, yield stress/strain, UTS, strain at
                               UTS, Poisson's ratio, plus the load-scale
                               provenance: scale, its R², and the DIC coverage
                               fraction) — the single source of truth for
                               per-coupon scalars, read back out by Level-3.
  <DIC_DIR>/level2_group_stats.csv   D638 §11.7 mean/std/count by exposure×direction
"""

from __future__ import annotations
import os
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.signal import butter, filtfilt

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# PATHS
# =============================================================================
DIC_DIR = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\DIC"
)
# Specimen sheet. THE CSV *IS* THE SHEET — there is no .xlsx any more.
#
# "Width / Dia. (in)" and "Computed Area (in²)" used to be FORMULA cells in
# FSR-SpecimenTesting.xlsx, and openpyxl does not evaluate formulas: every time
# a Level 2 saved its scalars back into the workbook it wrote the formula and
# dropped the cached value, so Level 1 read those columns back as blank, wrote
# area_mm2 = NaN, and Level 2 then reported "insufficient data" for every
# coupon. The workbook is retired. FSR-SpecimenTesting.csv holds evaluated
# values, needs no Excel engine, is not locked while something else has it
# open, and is what every script in this folder now both reads and writes.
SPECIMEN_CSV = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.csv"
)

# The CSV started life as a Windows Excel export, so its superscript characters
# ("Computed Area (in²)") may still be cp1252 rather than UTF-8. Try in this
# order; write_specimen_sheet re-writes the file as utf-8-sig.
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "UV": True, "SW": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

APPLY_SMOOTHING = False  # ASTM D638 does not call for filtering the stress-strain record;
                          # toggle on only if a batch's raw signal is genuinely too noisy.
FILTER_METHOD   = "butterworth"  # "median" or "butterworth" — see SMOOTHING section below

# =============================================================================
# FAILURE TRUNCATION  — restricts Level-1 data to the valid test window
# before property extraction. Trims pre-load slack and post-fracture
# rebound; frames outside this window get NaN in every L2-derived column.
# =============================================================================
LOAD_START_FRAC = 0.02     # pre-load: drop frames before load exceeds this × peak
LOAD_END_FRAC   = 0.50     # post-fracture: cut first post-UTS frame where load < this × peak

# Scale factor applied to load_raw to produce force_N (N per sync-CSV unit).
#
# LOAD_SCALE_MODE selects how it is found:
#   "regress" (default) — least squares of the raw MTS force against load_raw
#                over the rising ramp, using every point in it.
#   "peak"     — the original rule: mts_peak_N / max(|load_raw|).
#   Either falls back to the other, then to SCALE_N_PER_UNIT, rather than
#   failing closed on a coupon whose MTS file is missing or unreadable.
#
# WHY THE DEFAULT CHANGED. The peak ratio sets the entire stress axis — and
# therefore E, UTS and yield, all linearly — from the ratio of two single
# samples: the largest of a few thousand noisy MTS samples over the largest of
# a few thousand noisy sync samples, taken from two records that do not share a
# clock. Both maxima are biased high by their own noise floors and there is no
# reason the two biases match. It also assumes the DIC record contains the peak
# at all, and on three P01 coupons it does not (see DIC_COVERAGE_MIN).
#
# The load cell and DAQ gain are hardware constants, so every coupon in a batch
# should return the same scale. Measured across the 35 P01 tensile coupons:
#
#     peak ratio    mean 555.3  CV 2.55 %   range 540.0 – 613.8
#     regressed     mean 559.4  CV 0.71 %   range 551.0 – 566.6
#
# — a 3.6x reduction in the scatter of a factor that multiplies E directly.
# The regression also returns an R² as a free diagnostic on sync quality, which
# a ratio of two numbers cannot.
#
# The fitted intercept is deliberately NOT applied. A DC offset on the stress
# axis does not change a slope; it is absorbed by the toe correction. It is
# printed only as a check on the sync channel's zero. (It does bias UTS and
# yield, but those are anchored to the MTS peak by other means.)
LOAD_SCALE_MODE = "regress"      # "regress" | "peak"
SCALE_N_PER_UNIT: float = 555.5928

MTS_DIR = DIC_DIR.parent / "MTS"   # raw MTS .txt, for the regressed load scale
MTS_HEADERS = 8                    # header rows before row 1 of data

# Raw VIC-3D project dirs, same mapping Level 1 and Level2_tmp use. Needed only
# as a FALLBACK clock source: per-coupon CSVs written by a Level 1 older than
# 2026-08-29 have a flattened time_s (an absolute epoch rounded to 6
# significant figures — every row identical), which cannot anchor anything.
# When that is detected, the full-precision Time_0_0 is re-read from the sync
# CSV here, so Level 2 works on existing files without a full Level-1 re-run.
DATA_ROOTS = {
    "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
    "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
}

# Band of the ramp used for the regressed scale, as a fraction of peak force.
# Bounded below to stay off the noise floor and above to stay out of the
# roll-over near peak, where any residual clock drift shows up as curvature.
SCALE_FIT_BAND = (0.10, 0.85)

# The two records don't share a clock, so they are anchored peak-to-peak and
# the residual lag is refined by cross-correlation over ±LAG_SEARCH_FRAC of the
# test duration.
#
# LAG_SEARCH_FRAC MUST BE GENEROUS. At 0.05 the search rails against its own
# limit on 3 of 35 P01 coupons (TCL45-01, TSW00-01, TSW00-02 — true lags of
# −12.3, −6.0 and −5.5 s) and returns a scale 5–17 % wrong while R² stays above
# 0.996, so the R² guard does not catch it. At 0.25 every P01 coupon converges
# with the optimum in the interior. A railed optimum is reported and rejected.
LAG_SEARCH_FRAC = 0.25
LAG_SEARCH_STEPS = 1001

SCALE_R2_MIN = 0.98        # below this, warn: sync/MTS agreement is poor

# Fraction of the MTS peak force the DIC record must actually reach. The sync
# CSV can stop before the specimen fails, and then max(|load_raw|) is not the
# peak at all and the peak-ratio scale is inflated by however much of the ramp
# was missed. On P01: TCL45-01 reaches 90.6 %, TSW00-01 96.3 %, TSW00-02
# 97.2 %; every other coupon is above 99.8 %. TCL45-01 is the coupon Level 3
# excludes for an "anomalous toe" — its stress axis is simply 10 % too high.
DIC_COVERAGE_MIN = 0.98

# =============================================================================
# PROPERTY SETTINGS
# =============================================================================
# Modulus fit window (axial strain, dimensionless).
# D638 §11.4: "initial linear portion of the load-extension curve".
# A window of 0.05–0.3% covers the typical linear region for stiff polymers
# and composites without including the toe. Adjust if the fit line on the
# generated plot doesn't sit on the linear segment.
MODULUS_STRAIN_RANGE = (0.0005, 0.003)

# Window for the chord modulus. DELIBERATELY NOT MODULUS_STRAIN_RANGE.
#
# The truncated record does not always reach down to 5e-4 corrected strain. The
# analysis window starts at the last frame below LOAD_START_FRAC x peak, and the
# NEXT frame — the first one kept — can already be well past that threshold: at
# 10 Hz on a ~50 s ramp one frame is ~2 % of the ramp, so the first kept frame
# lands anywhere between 2 % and 8.5 % of peak load depending on phase. On 18 of
# the 35 P01 coupons that puts the lowest available corrected strain between
# 5.6e-4 and 9.5e-4, i.e. ABOVE the modulus window's floor.
#
# The tangent fit tolerates this (it uses whatever points fall inside the
# window) but a chord cannot: it is defined by its two endpoints, so an
# unavailable endpoint has to be either extrapolated or refused. chord_modulus
# refuses — np.interp would otherwise clamp at the end of the data and return a
# confident-looking number computed from the wrong point.
#
# 1e-3 is the lowest round floor every P01 coupon actually reaches (worst case
# 9.47e-4). Widen MODULUS_STRAIN_RANGE's floor or lower LOAD_START_FRAC and
# this can come back down; see the README for that open decision.
CHORD_STRAIN_RANGE = (0.001, 0.003)

# D638 §A2.6 — 0.2% offset yield strength
YIELD_OFFSET = 0.002

# D638 §A3.10.1.3 — Poisson chord method window (when no clear proportionality)
# Chord computed at ε_a = 0.002 over the range 0.0005 to 0.0025 strain.
POISSON_RANGE = (0.0005, 0.0025)
POISSON_CHORD_AT = 0.002

# Scalar property columns written into SPECIMEN_CSV, keyed by coupon
# ("Specimen ID") — maps the property dict key to the sheet column header.
# Level-3 reads these same headers back out (per-coupon plots, group plots,
# and the printed/exported stat tables all read from this one file).
SPECIMEN_SHEET_COLUMNS = {
    "E_GPa":         "E (GPa)",
    "E_chord_GPa":   "E chord (GPa)",
    "eps_toe":       "Toe Strain",
    "n_fit":         "Modulus Fit Points",
    "scale_N_per_unit": "Load Scale (N/unit)",
    "scale_r2":      "Load Scale R2",
    "dic_coverage":  "DIC Load Coverage",
    "sigma_y_MPa":   "Yield Stress (MPa)",
    "eps_y":         "Yield Strain",
    "UTS_MPa":       "UTS (MPa)",
    "eps_at_UTS":    "Strain at UTS",
    "poisson_chord": "Poisson's Ratio (chord)",
    "poisson_slope": "Poisson's Ratio (slope)",
}

# =============================================================================
# SMOOTHING  — FILTER_METHOD selects which of these is used
#   "median"     : rolling median. Window must be odd. Raise MEDIAN_WINDOW
#                  for noisier data. Preferred default — doesn't ring or
#                  systematically undershoot a sharp peak the way an
#                  averaging-based filter like Butterworth can.
#   "butterworth": zero-phase low-pass (filtfilt). Lower BUTTER_CUTOFF for
#                  heavier smoothing; raise BUTTER_ORDER for a sharper
#                  rolloff. Output is clipped to the raw data's range to
#                  suppress filtfilt ringing at the truncation edges.
# =============================================================================
MEDIAN_WINDOW = 31  # frames

BUTTER_ORDER  = 3     # filter order
BUTTER_CUTOFF = 0.1   # cutoff frequency, fraction of Nyquist (0-1) — must be < 1

# =============================================================================
# HELPERS
# =============================================================================
def coupon_id(p, e, d, r): return f"{p}-T{e}{d}-{r}"

def selected_coupons():
    return [coupon_id(p, e, d, r)
            for p in PRINTS
            for e, on in EXPOSURES.items() if on
            for d, on2 in DIRECTIONS.items() if on2
            for r in REPLICATES]

def parse_id(cid):
    """Return (exposure, direction_str) e.g. ('CL', '00')."""
    part = cid.split("-")[1]
    return part[1:-2], part[-2:]

def find_frames_csv(cid):
    p = DIC_DIR / f"{cid}.csv"
    return p if p.exists() else None

def load_coupon_scalars() -> pd.DataFrame:
    """coupon_scalars.csv, indexed by coupon — mts_peak_N and area_mm2,
    written once per coupon by Level-1 (never repeated per row)."""
    fp = DIC_DIR / "coupon_scalars.csv"
    if not fp.exists():
        return pd.DataFrame(columns=["mts_peak_N", "area_mm2"])
    return pd.read_csv(fp).set_index("coupon")

def frame_clock(cid: str, df: pd.DataFrame) -> np.ndarray | None:
    """Elapsed seconds per DIC frame, or None if no usable clock exists.

    Prefers the per-coupon CSV's own time_s. Falls back to re-reading
    Time_0_0 from the raw sync CSV when that column is degenerate — which it is
    on every file written by a Level 1 older than 2026-08-29, where an absolute
    Unix epoch was pushed through float_format="%.6g" and came out as one
    repeated constant. A clock with fewer than two distinct values cannot place
    a frame in time, so it is rejected outright rather than silently producing a
    flat interpolation.
    """
    if "time_s" in df:
        t = pd.to_numeric(df["time_s"], errors="coerce").to_numpy(dtype=float)
        finite = t[np.isfinite(t)]
        if finite.size and np.unique(finite).size > 1:
            return t - finite[0]

    exposure, _ = parse_id(cid)
    root = DATA_ROOTS.get(exposure)
    if root is None:
        return None
    sync_fp = root / cid / f"{cid}.csv"
    if not sync_fp.exists():
        return None
    try:
        sync = pd.read_csv(sync_fp)
    except Exception:
        return None
    col = next((c for c in sync.columns if "time" in c.lower()), None)
    if col is None:
        return None
    t = pd.to_numeric(sync[col], errors="coerce").to_numpy(dtype=float)[:len(df)]
    finite = t[np.isfinite(t)]
    if not finite.size or np.unique(finite).size <= 1:
        return None
    t = t - finite[0]
    if len(t) < len(df):                       # sync shorter than the record
        t = np.concatenate([t, np.full(len(df) - len(t), np.nan)])
    return t


def read_mts(cid: str) -> pd.DataFrame | None:
    """Raw MTS record for one coupon, or None. Columns are already mm/N/V/s.

    Tab-separated with an 8-line header, same as Level 1's load_mts_txt — the
    separator matters: r"\\s+" also splits on the decimal-aligned padding some
    rows carry and silently shifts columns.
    """
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    if not hits:
        return None
    try:
        return (pd.read_csv(hits[0], sep="\t", skiprows=MTS_HEADERS, header=None,
                            names=["disp_mm", "force_N", "output_V", "time_s"],
                            encoding="utf-8-sig", on_bad_lines="skip")
                  .apply(pd.to_numeric, errors="coerce")
                  .dropna(subset=["force_N"]))
    except Exception:
        return None


def load_scale_for(cid, load_raw, df, mts_peak, raw_peak):
    """Newtons per sync-CSV unit. Returns (scale, r2, coverage, note).

    See LOAD_SCALE_MODE above for why the regression is preferred over the peak
    ratio, and what it measured on P01.

    `coverage` is the largest MTS force the DIC record actually spans, as a
    fraction of the MTS peak. Below DIC_COVERAGE_MIN the sync CSV stopped
    before the specimen did, max(|load_raw|) is not the peak, and any
    peak-anchored scale for that coupon is inflated. It is NaN when the
    regression could not run.

    Falls back to the peak ratio, then to SCALE_N_PER_UNIT, so this can never
    fail closed on a coupon whose MTS file is missing.
    """
    fallback = (mts_peak / raw_peak
                if mts_peak is not None and raw_peak > 0 else SCALE_N_PER_UNIT)
    if LOAD_SCALE_MODE == "peak":
        return fallback, np.nan, np.nan, "peak ratio (LOAD_SCALE_MODE)"

    mts = read_mts(cid)
    if mts is None:
        return fallback, np.nan, np.nan, "peak ratio (no MTS file)"
    time_s = frame_clock(cid, df)
    if time_s is None:
        return fallback, np.nan, np.nan, "peak ratio (no usable frame clock)"

    f, t = mts["force_N"].to_numpy(), mts["time_s"].to_numpy()
    if not (np.isfinite(f).any() and np.isfinite(t).any()):
        return fallback, np.nan, np.nan, "peak ratio (no usable MTS clock)"
    t = t - t[np.isfinite(t)][0]

    # Anchor the two clocks peak-to-peak, then refine the lag by correlation.
    # The argmax is itself a single sample and the force curve is flat near
    # peak, so on a noisy record it can land well off the true peak; it would be
    # poor form to replace a single-sample scale and keep a single-sample time
    # anchor. And when the DIC record is truncated (coverage below), the DIC
    # "peak" is not the peak at all and the anchor is out by many seconds —
    # which is exactly the case a narrow lag search cannot recover from.
    base = (time_s - float(time_s[int(np.nanargmax(np.abs(load_raw)))])
            + float(t[int(np.nanargmax(np.abs(f)))]))
    duration = float(np.nanmax(time_s) - np.nanmin(time_s))
    span = LAG_SEARCH_FRAC * duration if np.isfinite(duration) and duration > 0 else 0.0
    best_lag, best_r = 0.0, -np.inf
    for lag in (np.linspace(-span, span, LAG_SEARCH_STEPS) if span > 0 else [0.0]):
        y = np.interp(base + lag, t, f, left=np.nan, right=np.nan)
        m = np.isfinite(y) & np.isfinite(load_raw)
        if m.sum() < 50:
            continue
        a, b = load_raw[m], y[m]
        if a.std() == 0 or b.std() == 0:
            continue
        r = float(np.mean((a - a.mean()) * (b - b.mean())) / (a.std() * b.std()))
        if r > best_r:
            best_r, best_lag = r, float(lag)
    if span > 0 and abs(abs(best_lag) - span) < span / (LAG_SEARCH_STEPS - 1):
        # Optimum sits on the edge of the search: the true lag is outside it and
        # the fit is against a misaligned ramp. Do not use it.
        return (fallback, np.nan, np.nan,
                f"peak ratio (lag search railed at {best_lag:+.2f} s — "
                f"raise LAG_SEARCH_FRAC)")

    f_mts = np.interp(base + best_lag, t, f, left=np.nan, right=np.nan)
    mts_peak_abs = float(np.nanmax(np.abs(f)))
    coverage = (float(np.nanmax(np.abs(f_mts))) / mts_peak_abs
                if mts_peak_abs > 0 and np.isfinite(f_mts).any() else np.nan)

    pk = float(np.nanmax(np.abs(f_mts))) if np.isfinite(f_mts).any() else 0.0
    lo_f, hi_f = SCALE_FIT_BAND
    band = (np.isfinite(f_mts) & np.isfinite(load_raw)
            & (np.abs(f_mts) > lo_f * pk) & (np.abs(f_mts) < hi_f * pk))
    if band.sum() < 20:
        return fallback, np.nan, coverage, "peak ratio (too few in-band points)"

    x, y = load_raw[band], f_mts[band]
    slope, icept = (float(v) for v in np.polyfit(x, y, 1))
    resid = y - (slope * x + icept)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else np.nan
    return slope, r2, coverage, (
        f"regressed on {int(band.sum())} pts, lag {best_lag:+.3f} s, "
        f"offset {icept:+.1f} N (not applied)")


def smooth_signal(x):
    """Dispatches to the filter selected by FILTER_METHOD. No-op when
    APPLY_SMOOTHING is False."""
    x = np.asarray(x, dtype=float)
    if not APPLY_SMOOTHING:
        return x.copy()
    if FILTER_METHOD == "butterworth":
        return _smooth_butterworth(x)
    return _smooth_median(x)

def _smooth_median(x):
    """Rolling median. mode='nearest' avoids the zero-padding edge artifacts
    scipy.signal.medfilt has."""
    win = MEDIAN_WINDOW
    if win % 2 == 0:
        win -= 1
    if win < 1 or len(x) < win:
        return x.copy()
    return median_filter(x, size=win, mode="constant", cval=np.nan)

def _smooth_butterworth(x):
    """Zero-phase low-pass (filtfilt). NaN gaps are bridged by linear
    interpolation before filtering, then the original NaN positions are
    restored. Output is clipped to the raw data's range — filtfilt can ring
    past it, especially near the truncation edges."""
    n = len(x)
    padlen = 3 * BUTTER_ORDER
    nan_mask = ~np.isfinite(x)
    if nan_mask.all() or n <= padlen:
        return x.copy()
    xi = x.copy()
    if nan_mask.any():
        idx = np.arange(n)
        xi[nan_mask] = np.interp(idx[nan_mask], idx[~nan_mask], x[~nan_mask])
    b, a = butter(BUTTER_ORDER, BUTTER_CUTOFF, btype="low", analog=False)
    out = filtfilt(b, a, xi)
    out = np.clip(out, np.nanmin(xi), np.nanmax(xi))
    out[nan_mask] = np.nan
    return out

def read_specimen_csv() -> pd.DataFrame:
    """The specimen sheet as raw text — every cell a string, nothing coerced.

    dtype=str with keep_default_na=False is what makes the file safe to write
    back: cells this script does not touch round-trip character for character,
    so a Level-2 run cannot reformat a hand-entered geometry value or turn an
    all-integer column into floats on its way through pandas.
    """
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc,
                               dtype=str, keep_default_na=False)
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f"{SPECIMEN_CSV.name}: not decodable as "
                       f"{'/'.join(CSV_ENCODINGS)}")

def _cell(v) -> str:
    """One specimen-sheet cell as text. Missing or non-finite writes blank, so
    a property that could not be computed clears the cell instead of leaving
    the previous run's value standing."""
    if v is None:
        return ""
    v = float(v)
    return "" if not np.isfinite(v) else f"{v:.12g}"

def write_specimen_sheet(rows: list[dict]) -> None:
    """Write each coupon's scalar properties into its row of SPECIMEN_CSV,
    matched by Specimen ID. Adds any missing property columns at the end;
    every other cell is written back exactly as it was read. Skipped with a
    warning if the file can't be read or replaced — e.g. open in Excel.
    """
    try:
        df = read_specimen_csv()
    except FileNotFoundError:
        print(f"[!] {SPECIMEN_CSV} not found — skipping specimen sheet update")
        return
    except Exception as exc:
        print(f"[!] {SPECIMEN_CSV.name}: {exc} — skipping specimen sheet update")
        return

    if "Specimen ID" not in df.columns:
        print("[!] 'Specimen ID' column not found in specimen sheet — skipping update")
        return

    for label in SPECIMEN_SHEET_COLUMNS.values():
        if label not in df.columns:
            df[label] = ""
    col_pos  = {c: j for j, c in enumerate(df.columns)}
    row_by_id = {cid: i for i, cid in enumerate(df["Specimen ID"])}

    n_written = 0
    for row in rows:
        i = row_by_id.get(row["coupon"])
        if i is None:
            print(f"[!] {row['coupon']} has no row in the specimen sheet — "
                  f"its properties were not written")
            continue
        for key, label in SPECIMEN_SHEET_COLUMNS.items():
            df.iat[i, col_pos[label]] = _cell(row.get(key))
        n_written += 1

    # Write alongside the target and rename over it. The specimen sheet is now
    # the only copy of the hand-entered geometry, so a half-written file would
    # be real data loss rather than just a failed run.
    tmp = SPECIMEN_CSV.with_name(SPECIMEN_CSV.name + ".tmp")
    try:
        df.to_csv(tmp, index=False, encoding="utf-8-sig")
        os.replace(tmp, SPECIMEN_CSV)
        print(f"Specimen sheet: {n_written} coupon(s) → {SPECIMEN_CSV.name}")
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        print(f"[!] could not update {SPECIMEN_CSV.name} ({exc}) — "
              f"properties were not saved")

def read_existing_group_stats(fp: Path) -> pd.DataFrame | None:
    """Read level2_group_stats.csv, upgrading the pre-'test' two-level layout.

    The index depth has to be READ, not assumed. pandas writes a named
    MultiIndex as a third header line holding the level names, so that line says
    how many index columns the file on disk really has — and reading a
    two-level file with index_col=[0,1,2] does not raise. It quietly swallows
    the first data column into the index and shifts every value one place left,
    which looks like a successful merge and silently corrupts the other test
    type's rows.

    Returns None if the file can't be understood, so the caller can say so
    rather than write a half-merged table.
    """
    try:
        with open(fp, encoding="utf-8") as fh:
            head = [fh.readline() for _ in range(3)]
        names = [c.strip() for c in head[2].rstrip(chr(13) + chr(10)).split(",") if c.strip()]
        if not names:
            return None
        old = pd.read_csv(fp, header=[0, 1], index_col=list(range(len(names))))
        if "test" not in names:
            # Written before the 'test' level existed, when only the tensile
            # pipeline wrote this file. Everything in it is tensile.
            old = pd.concat({"tensile": old}, names=["test"])
        # A direction of "00" round-trips through read_csv as the integer 0, so
        # normalise every level back to the string form both scripts write.
        old.index = pd.MultiIndex.from_tuples(
            [(str(t), str(e), f"{int(d):02d}" if str(d).strip().isdigit() else str(d))
             for t, e, d in old.index],
            names=["test", "exposure", "direction"])
        return old
    except Exception:
        return None


def write_group_stats(rows: list[dict]) -> Path:
    """Upsert this run's group statistics into level2_group_stats.csv.

    The file is shared with FlexuralDIC_Level2 and indexed by
    (test, exposure, direction). Both test types have CL and IS exposures at 00
    and 90, so without the 'test' level the flexural rows and the tensile ones
    would land on top of each other. Only rows whose test == "tensile" are
    replaced; anything else already in the file is read back and written out
    unchanged.
    """
    fp = DIC_DIR / "level2_group_stats.csv"
    df_sum = pd.DataFrame(rows)
    df_sum["test"] = "tensile"
    df_sum["exposure"]  = df_sum["coupon"].map(lambda c: parse_id(c)[0])
    df_sum["direction"] = df_sum["coupon"].map(lambda c: parse_id(c)[1])
    agg_cols = ["E_GPa", "E_chord_GPa", "sigma_y_MPa", "UTS_MPa", "eps_at_UTS",
                "poisson_chord"]
    group = (df_sum.groupby(["test", "exposure", "direction"])[agg_cols]
                   .agg(["mean", "std", "count"]))

    if fp.exists():
        try:
            old = pd.read_csv(fp, header=[0, 1], index_col=[0, 1, 2])
            old = old.drop(index="tensile", level=0, errors="ignore")
            group = pd.concat([old, group]).sort_index()
        except Exception as ex:
            # An older two-level file, written before the 'test' level existed,
            # can't be merged column-wise. Say so rather than silently dropping
            # it; re-running FlexuralDIC_Level2 restores the flexural rows.
            print(f"[!] could not merge existing {fp.name} ({ex}) — it is being "
                  f"replaced with tensile rows only.\n    Re-run "
                  f"FlexuralDIC_Level2.py to put the flexural rows back.")
    group.to_csv(fp)
    return fp


def truncation_mask(force_N: np.ndarray) -> np.ndarray:
    """Boolean mask, same length as force_N, selecting the valid test
    window: pre-load slack and post-fracture rebound are False."""
    n = len(force_N)
    peak = float(np.nanmax(np.abs(force_N)))
    if peak <= 0:
        return np.ones(n, dtype=bool)
    i_uts = int(np.nanargmax(np.abs(force_N)))
    # Start on the *rising edge*, not the first sample over threshold: the
    # pre-touchdown baseline sits near zero with noise either side of it, so
    # "first sample above 2 % of peak" triggers on noise long before the
    # specimen is loaded. Take the last sample below threshold before the peak.
    # Same fix FlexuralDIC_Level2.truncation_mask already carries.
    #
    # On P01 this moves the start from frame 0–12 to frame 41–84 — i.e. the old
    # rule was keeping 40–80 frames of pre-touchdown noise, and some of those
    # frames carry near-zero stress with strain scattered into the modulus fit
    # window, where they flatten the fit.
    f_abs = np.abs(np.nan_to_num(force_N, nan=0.0))
    below = np.flatnonzero(f_abs[:i_uts + 1] < LOAD_START_FRAC * peak)
    i0 = int(below[-1]) + 1 if below.size else 0
    post = np.where(np.abs(force_N[i_uts:]) < LOAD_END_FRAC * peak)[0]
    # stop one frame BEFORE the first post-UTS drop-off so the failure point
    # itself isn't kept (it corrupts the smoothing pass)
    i1 = int(i_uts + post[0]) - 1 if len(post) else n - 1
    mask = np.zeros(n, dtype=bool)
    # max(i0, i1) so a pathological record with i1 < i0 yields an empty window
    # rather than a reversed slice that selects nothing while looking like it
    # worked.
    mask[i0:max(i0, i1) + 1] = True
    return mask


def fit_modulus(eps: np.ndarray, sig: np.ndarray,
                window: tuple[float, float]) -> tuple[float, float, int]:
    """Least-squares slope of sigma against eps inside `window`.

    Returns (slope MPa, toe offset in strain, n points). The toe offset is the
    fit line's x-intercept — D638 Annex A1: extend the straight line back to
    zero stress and take that as the true strain origin.

    Fitted twice, and this is the same routine FlexuralDIC_Level2.fit_modulus
    uses, so the two pipelines now compute a modulus the same way.
    MODULUS_STRAIN_RANGE is a window in *corrected* strain, but the correction
    is what the first fit produces, so the second pass re-selects the window
    after shifting the origin.

    Without the second pass, E is measured over corrected strain
    [lo - toe, hi - toe] rather than [lo, hi], which contradicts the comment on
    MODULUS_STRAIN_RANGE ("without including the toe"). It also makes the toe
    correction a mathematical no-op for E: subtracting a constant from the
    abscissa after the window is fixed shifts the origin without changing the
    slope, so a single-pass "toe-corrected" E and an uncorrected one are equal
    to the last digit.

    SCALE OF THE EFFECT ON P01: small. Measured toe offsets are 1e-5 to 1e-4
    strain, an order of magnitude below the window's 5e-4 lower bound, so the
    iterated fit moves mean E by 0.16 %. The exception is P01-TCL45-01
    (toe 1.1e-3, twice the window floor) — the coupon Level 3 excludes. This is
    a correctness fix, not a large numerical one; see
    tensile_modulus_sensitivity.py for the measurement.
    """
    toe = 0.0
    slope = np.nan
    n = 0
    for _ in range(2):
        e = eps - toe
        m = (e >= window[0]) & (e <= window[1]) & np.isfinite(e) & np.isfinite(sig)
        n = int(m.sum())
        if n < 3:
            return np.nan, np.nan, n
        slope, icept = np.polyfit(e[m], sig[m], 1)
        if slope == 0:
            return np.nan, np.nan, n
        toe = toe + (-icept / slope)      # x-intercept, accumulated
    return float(slope), float(toe), n


def chord_modulus(eps: np.ndarray, sig: np.ndarray,
                  window: tuple[float, float]) -> float:
    """Secant between the two ends of `window`, in already-toe-corrected strain.

    Reported alongside the tangent modulus for the reason FlexuralDIC_Level2
    reports Ef_*_chord_GPa alongside its tangent value: on this material the
    tangent modulus is not a stable quantity. Sweeping the fit window over five
    plausible ranges moves E by a median 7.7 % and up to 21.8 % per coupon —
    larger than the 3.9 % within-group scatter it is being used to compare —
    and the tangent fit's own R² is a median 0.879 (min 0.808), worst on the
    45°/90° coupons.

    A CHORD IS NOT A BETTER ESTIMATOR, and this one is not offered as the
    number to report instead. It is read from two interpolated points, so where
    the tangent averages strain noise over ~60 points the chord takes it at
    face value: on P01 the chord itself moves 15–20 % between two neighbouring
    windows on some coupons. What it gives you is transparency — a chord is
    exactly "the average stiffness between these two strains", which is a claim
    that survives being asked what window you used.

    The pair is the useful output: tangent and chord agreeing means the segment
    really is straight, and disagreeing (median +6.1 %, up to 33.6 % here)
    means it isn't. Both go in the specimen sheet. See the README for what to do
    about that.

    Returns NaN rather than extrapolating when the data does not reach both
    endpoints — np.interp would otherwise clamp at the end of the record and
    return a confident-looking number computed from the wrong point.
    """
    m = np.isfinite(eps) & np.isfinite(sig)
    if m.sum() < 3:
        return np.nan
    order = np.argsort(eps[m])
    e, s = eps[m][order], sig[m][order]
    lo, hi = window
    if not (e[0] <= lo and hi <= e[-1]):
        return np.nan
    return float((np.interp(hi, e, s) - np.interp(lo, e, s)) / (hi - lo))


# =============================================================================
# COMPUTE PROPERTIES
# =============================================================================
def compute_properties(eps_axial, sig, eps_transverse,
                       eps_axial_unsmoothed=None, sig_unsmoothed=None):
    """
    D638-compliant property extraction from one coupon's truncated (and, if
    enabled, smoothed) axial-strain / stress / transverse-strain arrays.

    Returns a dict with E_GPa, eps_toe (toe-correction offset applied),
    sigma_y_MPa, eps_y, UTS_MPa, eps_at_UTS, poisson_chord, poisson_slope,
    the toe-corrected _eps/_sig/_eps_t arrays, a _valid boolean mask (which
    input rows survived the finite-value filter — needed by the caller to
    scatter these arrays back into the full per-frame file), and, when
    eps_axial_unsmoothed/sig_unsmoothed were given, their toe-corrected
    pre-smoothing counterparts _eps_raw/_sig_raw for Level-3's diagnostic
    overlay (identical to _eps/_sig when smoothing is off).
    """
    eps_axial = np.asarray(eps_axial, dtype=float)
    sig = np.asarray(sig, dtype=float)
    eps_transverse = np.asarray(eps_transverse, dtype=float)

    valid = np.isfinite(eps_axial) & np.isfinite(sig)
    if valid.sum() < 10:
        return None

    eps_raw   = eps_axial[valid]
    sig       = sig[valid]
    eps_t_raw = eps_transverse[valid]
    eps_unsmoothed = (np.asarray(eps_axial_unsmoothed, dtype=float)[valid]
                       if eps_axial_unsmoothed is not None else None)
    sig_unsmoothed = (np.asarray(sig_unsmoothed, dtype=float)[valid]
                       if sig_unsmoothed is not None else None)

    # ---- 1+2. Modulus (D638 §11.4) and toe compensation (Annex A1) ----------
    # One call, because they are not separable: the fit window is defined in
    # toe-corrected strain and the correction comes out of the fit, so
    # fit_modulus iterates. See its docstring.
    E_MPa, eps_offset, n_fit = fit_modulus(eps_raw, sig, MODULUS_STRAIN_RANGE)
    if not np.isfinite(E_MPa):
        return None
    eps   = eps_raw   - eps_offset
    eps_unsmoothed_corr = (eps_unsmoothed - eps_offset) if eps_unsmoothed is not None else None
    # Transverse strain: subtract its value at the corrected zero of axial strain.
    # Find the index where corrected axial ≈ 0 and subtract that ε_t.
    if np.any(np.isfinite(eps_t_raw)):
        i0 = int(np.nanargmin(np.abs(eps)))
        eps_t = eps_t_raw - eps_t_raw[i0]
    else:
        eps_t = eps_t_raw.copy()

    # ---- 3. UTS (D638 §11.2 — max stress) -----------------------------------
    i_uts   = int(np.nanargmax(sig))
    uts     = float(sig[i_uts])
    eps_ult = float(eps[i_uts])

    # ---- 4. 0.2% offset yield (D638 §A2.6, Fig. A2.1) -----------------------
    # First crossing of σ-ε curve with the line σ = E·(ε − YIELD_OFFSET).
    sigma_y, eps_y = np.nan, np.nan
    diff  = sig - E_MPa * (eps - YIELD_OFFSET)
    valid_y = np.where(eps > YIELD_OFFSET)[0]
    if len(valid_y) > 1:
        d = diff[valid_y]
        crossings = np.where(np.diff(np.sign(d)) < 0)[0]
        if len(crossings):
            k = valid_y[crossings[0]]
            denom = diff[k] - diff[k+1]
            f = diff[k] / denom if denom != 0 else 0.0
            eps_y   = float(eps[k] + f * (eps[k+1] - eps[k]))
            sigma_y = float(sig[k] + f * (sig[k+1] - sig[k]))

    # ---- 5. Poisson's ratio (D638 §A3.10.1.3, chord at ε_a = 0.002) --------
    #   ν = − ε_t(at ε_a = 0.002) / 0.002
    # Range 0.0005 – 0.0025 strain. Also report least-squares slope (§A3.10.1.1)
    # for transparency when proportionality holds.
    nu_chord = nu_slope = np.nan
    pm = ((eps >= POISSON_RANGE[0]) & (eps <= POISSON_RANGE[1]) &
          np.isfinite(eps) & np.isfinite(eps_t))
    if pm.sum() >= 3:
        nu_slope = float(-np.polyfit(eps[pm], eps_t[pm], 1)[0])
        order = np.argsort(eps[pm])
        ea, et = eps[pm][order], eps_t[pm][order]
        if ea[0] <= POISSON_CHORD_AT <= ea[-1]:
            nu_chord = float(-np.interp(POISSON_CHORD_AT, ea, et) / POISSON_CHORD_AT)

    # ---- 6. Chord modulus ---------------------------------------------------
    # Not a D638 requirement; reported as a cross-check on the tangent value,
    # which is window-sensitive on this material. Over CHORD_STRAIN_RANGE, not
    # MODULUS_STRAIN_RANGE — see that constant for why the two differ.
    E_chord_MPa = chord_modulus(eps, sig, CHORD_STRAIN_RANGE)
    # The record not spanning even the chord window means the truncation start
    # has eaten into the strain range the properties are defined over, which is
    # worth saying rather than emitting a bare NaN.
    if not np.isfinite(E_chord_MPa):
        finite_eps = eps[np.isfinite(eps)]
        if finite_eps.size:
            print(f"    [!] record spans corrected strain "
                  f"[{finite_eps.min():.2e}, {finite_eps.max():.2e}] — does not "
                  f"cover CHORD_STRAIN_RANGE {CHORD_STRAIN_RANGE}; chord modulus "
                  f"not computed")

    return {
        "E_GPa":         E_MPa / 1000.0,
        "E_chord_GPa":   E_chord_MPa / 1000.0 if np.isfinite(E_chord_MPa) else np.nan,
        "eps_toe":       eps_offset,
        "n_fit":         n_fit,
        "sigma_y_MPa":   sigma_y,
        "eps_y":         eps_y,
        "UTS_MPa":       uts,
        "eps_at_UTS":    eps_ult,
        "poisson_chord": nu_chord,
        "poisson_slope": nu_slope,
        "_eps":          eps,        # toe-corrected strain, used in property calcs
        "_sig":          sig,
        "_eps_t":        eps_t,
        "_eps_raw":      eps_unsmoothed_corr,  # diagnostic overlay only
        "_sig_raw":      sig_unsmoothed,
        "_valid":        valid,      # which input rows survived the finite-value filter
    }


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("DIC_Level2 — scale, truncate, smooth, compute D638 properties")
    print("=" * 70)
    rows = []
    scalars = load_coupon_scalars()

    for cid in selected_coupons():
        frames_fp = find_frames_csv(cid)
        if frames_fp is None:
            print(f"[{cid}] no per-coupon CSV — run Level 1 first")
            continue
        df = pd.read_csv(frames_fp)
        n = len(df)
        load_raw = df["load_raw"].to_numpy()

        # ---- scale load_raw -> force_N --------------------------------------
        raw_peak = float(np.nanmax(np.abs(load_raw)))
        mts_peak = (float(scalars.loc[cid, "mts_peak_N"])
                    if cid in scalars.index and pd.notna(scalars.loc[cid, "mts_peak_N"])
                    else None)
        area = (float(scalars.loc[cid, "area_mm2"])
                if cid in scalars.index and pd.notna(scalars.loc[cid, "area_mm2"])
                else np.nan)
        scale, scale_r2, coverage, scale_note = load_scale_for(
            cid, load_raw, df, mts_peak, raw_peak)
        print(f"[{cid}] scale {scale:.4f} N/unit  ({scale_note}"
              + (f", R²={scale_r2:.5f}" if np.isfinite(scale_r2) else "") + ")")
        if np.isfinite(scale_r2) and scale_r2 < SCALE_R2_MIN:
            print(f"    [!] sync/MTS agreement is poor (R² < {SCALE_R2_MIN}) — "
                  f"check the alignment before trusting this coupon's stresses")
        if np.isfinite(coverage) and coverage < DIC_COVERAGE_MIN:
            print(f"    [!] the DIC record only spans {coverage*100:.1f} % of the "
                  f"MTS peak force — it stopped before the specimen did.\n"
                  f"        UTS and strain-at-UTS for this coupon are NOT the "
                  f"specimen's; a peak-anchored scale would be inflated ~"
                  f"{(1/coverage - 1)*100:.0f} %.")
        force_N_all = load_raw * scale

        # ---- failure truncation ---------------------------------------------
        kept = truncation_mask(force_N_all)
        kept_idx = np.flatnonzero(kept)

        force_kept    = force_N_all[kept]
        stress_kept   = force_kept / area if np.isfinite(area) else np.full(kept.sum(), np.nan)
        eps_a_kept    = df["strain_axial_raw"].to_numpy()[kept]
        eps_t_kept    = df["strain_transverse_raw"].to_numpy()[kept]

        # Smoothing pass (no-op unless APPLY_SMOOTHING is on) — kept
        # pre-smoothing copies alongside for Level-3's diagnostic overlay.
        force_smoothed  = smooth_signal(force_kept)
        stress_smoothed = force_smoothed / area if np.isfinite(area) else np.full(kept.sum(), np.nan)
        eps_a_smoothed  = smooth_signal(eps_a_kept)
        eps_t_smoothed  = smooth_signal(eps_t_kept)

        p = compute_properties(eps_a_smoothed, stress_smoothed, eps_t_smoothed,
                               eps_axial_unsmoothed=eps_a_kept, sig_unsmoothed=stress_kept)
        if not p:
            print(f"[{cid}]  insufficient data")
            continue

        print(f"[{cid}]  E={p['E_GPa']:.2f} GPa  "
              f"E_chord={p['E_chord_GPa']:.2f} GPa  "
              f"σ_y={p['sigma_y_MPa']:.1f} MPa  "
              f"UTS={p['UTS_MPa']:.1f} MPa  "
              f"ε_UTS={p['eps_at_UTS']*100:.2f}%  "
              f"ν_chord={p['poisson_chord']:.3f}  "
              f"toe={p['eps_toe']*100:.3f}%")

        # ---- scatter L2 columns back onto the full per-frame file -----------
        final_idx = kept_idx[p["_valid"]]

        def scatter(values):
            col = np.full(n, np.nan)
            col[final_idx] = values
            return col

        df["kept"]                     = kept
        df["force_N"]                  = scatter(force_smoothed[p["_valid"]])
        df["stress_MPa"]               = scatter(p["_sig"])
        df["strain_axial"]             = scatter(p["_eps"])
        df["strain_transverse"]        = scatter(p["_eps_t"])
        df["stress_MPa_unsmoothed"]    = scatter(p["_sig_raw"])
        df["strain_axial_unsmoothed"]  = scatter(p["_eps_raw"])
        df.to_csv(frames_fp, index=False, float_format="%.6g")

        rows.append({"coupon": cid,
                     **{k: v for k, v in p.items() if not k.startswith("_")},
                     # Provenance for the stress axis, recorded per coupon so a
                     # suspect scale is visible in the specimen sheet rather than only
                     # in this run's console output.
                     "scale_N_per_unit": scale,
                     "scale_r2": scale_r2,
                     "dic_coverage": coverage})

    if rows:
        write_specimen_sheet(rows)

        # ---- D638 §11.7 / §12.1: mean & std per (exposure, direction) -------
        write_group_stats(rows)

        print(f"\n{len(rows)} coupon(s) → DIC/*.csv, {SPECIMEN_CSV.name}, "
              f"DIC/level2_group_stats.csv")

    print(f"\nDone. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
