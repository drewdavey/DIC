#!/usr/bin/env python3
"""
TensileDIC_Level2.py  —  FSR Tensile Coupons (ASTM D638)
=========================================================
Reads Level-1's per-frame CSV, puts the raw MTS force onto the DIC frames,
truncates to the valid test window, and computes D638 properties.
No plotting — that is Level 3's job.

FORCE COMES FROM THE MTS RECORD, NOT THE SYNC LOAD CHANNEL. The two files do
not share a clock, so the sync CSV's peak-force row and the raw MTS peak-force
row are treated as the same physical instant, and MTS force and displacement
are interpolated at each DIC frame's offset from that anchor. The sync load
channel is used only to locate that peak — its values are never carried into a
stress. See map_mts_force_disp and README.md.

Cited: toe compensation D638 Annex A1, modulus §11.4, 0.2 % offset yield
§A2.6, tensile strength §11.2, Poisson chord §A3.10.1.3, group statistics
§11.7. The chord modulus is not a D638 requirement; it is reported as a
cross-check because the tangent value is window-sensitive on this material.

INPUT
  <DIC_DIR>/<coupon_id>.csv      Level-1's per-frame record
  <DIC_DIR>/coupon_scalars.csv   per-coupon area_mm2
  <coupon_dir>/<coupon_id>.csv   raw sync CSV — read only for its clock and its
                                 peak-force row
  <MTS_DIR>/<coupon_id>*.txt     raw MTS force/displacement

OUTPUT
  <DIC_DIR>/<coupon_id>.csv      Level-1's columns plus: kept, force_N,
                                 stress_MPa, strain_axial, strain_transverse,
                                 stress_MPa_unsmoothed, strain_axial_unsmoothed,
                                 disp_mm_mts (the peak-anchored MTS
                                 displacement; Level-1's own disp_mm is left
                                 untouched). All NaN where kept is False.
  FSR-SpecimenTesting.csv        the D638 scalars, under the
                                 SPECIMEN_SHEET_COLUMNS headers below.
  <DIC_DIR>/level2_group_stats.csv   mean/std/count per exposure x direction,
                                 under a 'test' index level of "tensile".
                                 Shared with FlexuralDIC_Level2; each script
                                 replaces only its own test's rows.
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
MTS_DIR = DIC_DIR.parent / "MTS"
# Raw VIC-3D project dirs — the sync CSV's clock and peak-force row live here.
DATA_ROOTS = {
    "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
    "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
}
# The CSV is the specimen sheet — there is no .xlsx any more. See README.md.
SPECIMEN_CSV = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.csv"
)
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")

MTS_HEADERS = 8
DIC_FORCE_SCALED_COL = "LOAD_[kip]_|_CH07_/ai2_scaled"   # sync peak locator only

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "UV": True, "SW": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

APPLY_SMOOTHING = False           # D638 does not require filtering
FILTER_METHOD   = "butterworth"  # "median" or "butterworth"

# The peak anchor assumes the DIC record actually contains the peak. When the
# sync CSV stops before the specimen fails, its "peak" is only the largest force
# it happened to catch, and the anchor is wrong.
#
# That is measurable without trusting the sync channel's calibration. The load
# cell and DAQ gain are hardware constants, so sync_peak / mts_peak should be
# the same number for every coupon in a batch; a coupon that reads low did not
# record its own peak. add_coverage() takes the batch median as the nominal
# gain and reports each coupon against it. On P01 this flags TCL45-01 (0.90),
# TSW00-01 (0.94) and TSW00-02 (0.96) and nothing else.
DIC_COVERAGE_MIN = 0.95

# The peak anchor is refined by correlating the two force ramps over +/- this
# fraction of the test duration. It must stay generous: on P01-TSW00-01 and
# TSW00-02 the true correction is ~6 s, and a narrow search would not reach it.
LAG_SEARCH_FRAC  = 0.25
LAG_SEARCH_STEPS = 1001

# =============================================================================
# ANALYSIS  — must match TensileDIC_Level3 (and FlexuralDIC_Level2's windows)
# =============================================================================
# Analysis window, as a fraction of peak load.
LOAD_START_FRAC = 0.02     # outside the window until the load rises past this
LOAD_END_FRAC   = 0.50     # cut at the first post-UTS frame below this

# Modulus fit windows, in toe-corrected strain. D638 §11.4 "initial linear
# portion"; 0.05-0.3 % covers it for stiff polymers without including the toe.
MODULUS_STRAIN_RANGE = (0.0005, 0.003)

# Chord window. The floor is 1e-3, not 5e-4: the truncated record does not
# always reach down that far (worst P01 case 9.47e-4), and a chord cannot
# tolerate a missing endpoint the way the tangent fit can.
CHORD_STRAIN_RANGE = (0.001, 0.003)

YIELD_OFFSET     = 0.002             # D638 §A2.6
POISSON_RANGE    = (0.0005, 0.0025)  # D638 §A3.10.1.3
POISSON_CHORD_AT = 0.002

# Smoothing settings (used only when APPLY_SMOOTHING is True).
MEDIAN_WINDOW = 31    # frames, forced odd
BUTTER_ORDER  = 3
BUTTER_CUTOFF = 0.1   # fraction of Nyquist, must be < 1

# Scalar columns written into SPECIMEN_CSV, keyed by "Specimen ID".
# Level 3 reads these same headers back out.
SPECIMEN_SHEET_COLUMNS = {
    "E_GPa":         "E (GPa)",
    "E_chord_GPa":   "E chord (GPa)",
    "eps_toe":       "Toe Strain",
    "n_fit":         "Modulus Fit Points",
    "mts_peak_N":    "MTS Peak Load (N)",
    "dic_coverage":  "DIC Load Coverage",
    "sigma_y_MPa":   "Yield Stress (MPa)",
    "eps_y":         "Yield Strain",
    "UTS_MPa":       "UTS (MPa)",
    "eps_at_UTS":    "Strain at UTS",
    "poisson_chord": "Poisson's Ratio (chord)",
    "poisson_slope": "Poisson's Ratio (slope)",
}


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
    """Return (exposure, direction_str), e.g. ('CL', '00')."""
    part = cid.split("-")[1]
    return part[1:-2], part[-2:]

def coupon_dir(cid):
    return DATA_ROOTS[parse_id(cid)[0]] / cid

def load_coupon_scalars():
    """coupon_scalars.csv indexed by coupon — only area_mm2 is needed here,
    since force now comes straight from the mapped MTS record."""
    fp = DIC_DIR / "coupon_scalars.csv"
    if not fp.exists():
        return pd.DataFrame(columns=["area_mm2"])
    return pd.read_csv(fp).set_index("coupon")

def find_mts_txt(cid):
    for fp in (MTS_DIR / f"{cid}.txt", MTS_DIR / f"{cid}-TEST.txt"):
        if fp.exists():
            return fp
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    return hits[0] if hits else None

def find_sync_csv(cid):
    cdir = coupon_dir(cid)
    direct = cdir / f"{cid}.csv"
    if direct.exists():
        return direct
    if not cdir.is_dir():
        return None
    hits = sorted(cdir.rglob(f"{cid}.csv"))
    return hits[0] if hits else None

def pick_col(df, hint):
    for col in df.columns:
        if hint.lower() in str(col).lower():
            return col
    return None

def numeric(values):
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)

def relative_time(values, fallback_len):
    arr = numeric(values)
    finite = np.flatnonzero(np.isfinite(arr))
    if finite.size:
        return arr - arr[finite[0]]
    return np.arange(fallback_len, dtype=float)

def sample_rate_hz(time_s):
    dt = np.diff(time_s[np.isfinite(time_s)])
    dt = dt[np.isfinite(dt) & (dt > 0)]
    return float(1.0 / np.median(dt)) if dt.size else np.nan

def load_mts_txt(fp):
    """MTS .txt: 8-line header, tab separated, cols disp_mm force_N output_V
    time_s. The tab separator matters: r"\\s+" also splits on decimal-aligned
    padding and silently shifts columns."""
    return (pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                        names=["disp_mm", "force_N", "output_V", "time_s"],
                        encoding="utf-8-sig", on_bad_lines="skip")
              .apply(pd.to_numeric, errors="coerce")
              .dropna(how="all"))


# =============================================================================
# PEAK-ANCHORED MTS MAPPING
# =============================================================================
def refine_lag(base_time, sync_force, mts_time, mts_force):
    """Seconds to add to the peak anchor for the best correlation with the MTS
    force. Returns (lag, r).

    The peak anchor alone trusts a single sample of the sync load channel, and
    on this batch that channel is unreliable: on P01-TSW00-01 and TSW00-02 its
    maximum sits ~6 s away from the true peak instant, which anchored the whole
    record 6 s out and paired an unloaded specimen with 8 kN of force. Matching
    the two ramps as shapes, over every sample, is what catches that.
    """
    span = LAG_SEARCH_FRAC * float(np.nanmax(base_time) - np.nanmin(base_time))
    if not np.isfinite(span) or span <= 0:
        return 0.0, np.nan
    lag_best, r_best = 0.0, -np.inf
    for lag in np.linspace(-span, span, LAG_SEARCH_STEPS):
        y = np.interp(base_time + lag, mts_time, mts_force, left=np.nan, right=np.nan)
        m = np.isfinite(y) & np.isfinite(sync_force)
        if m.sum() < 50:
            continue
        a, b = sync_force[m], y[m]
        if a.std() == 0 or b.std() == 0:
            continue
        r = float(np.mean((a - a.mean()) * (b - b.mean())) / (a.std() * b.std()))
        if r > r_best:
            r_best, lag_best = r, float(lag)
    return lag_best, r_best


def map_mts_force_disp(cid, n_l1):
    """Map raw MTS force and displacement onto the DIC frame times.

    The sync CSV and the raw MTS file do not share a clock, so their two
    peak-force rows are treated as the same physical instant, and the residual
    lag is then refined by correlating the two ramps (see refine_lag — the peak
    alone is one sample and gets it badly wrong on some coupons). Each DIC
    frame's offset from that anchor is used to interpolate the MTS record,
    which puts the trustworthy MTS channel behind the WHOLE curve, not just its
    peak. Displacement is zeroed to its first finite mapped value.

    Returns (force_N, disp_mm, info), each array length n_l1 and NaN-padded
    past the sync CSV's length if that is shorter, or None.
    """
    sync_fp = find_sync_csv(cid)
    mts_fp = find_mts_txt(cid)
    if sync_fp is None or mts_fp is None:
        print(f"[{cid}] skip: sync={sync_fp is not None} MTS={mts_fp is not None}")
        return None

    sync = pd.read_csv(sync_fp)
    time_col = pick_col(sync, "time")
    if DIC_FORCE_SCALED_COL not in sync.columns or time_col is None:
        print(f"[{cid}] skip: sync CSV missing force or time column")
        return None

    n = min(n_l1, len(sync))
    sync = sync.iloc[:n].reset_index(drop=True)
    sync_time = relative_time(sync[time_col], len(sync))
    sync_force_kip = numeric(sync[DIC_FORCE_SCALED_COL])
    sync_peak_i = int(np.nanargmax(np.abs(sync_force_kip)))

    mts = load_mts_txt(mts_fp)
    mts_time = relative_time(mts["time_s"], len(mts))
    mts_force_N = numeric(mts["force_N"])
    mts_disp_mm = numeric(mts["disp_mm"])
    mts_peak_i = int(np.nanargmax(np.abs(mts_force_N)))

    anchored = mts_time[mts_peak_i] + (sync_time - sync_time[sync_peak_i])
    lag, lag_r = refine_lag(anchored, sync_force_kip, mts_time, mts_force_N)
    mapped_time = anchored + lag

    force_N = np.interp(mapped_time, mts_time, mts_force_N, left=np.nan, right=np.nan)
    disp_mm = np.interp(mapped_time, mts_time, mts_disp_mm, left=np.nan, right=np.nan)
    finite_disp = np.flatnonzero(np.isfinite(disp_mm))
    if finite_disp.size:
        disp_mm = disp_mm - disp_mm[finite_disp[0]]

    info = {"mts_peak_N": float(abs(mts_force_N[mts_peak_i])),
            "sync_peak_raw": float(abs(sync_force_kip[sync_peak_i])),
            "lag_s": lag, "lag_r": lag_r}

    print(f"[{cid}] sync_fs={sample_rate_hz(sync_time):.3f} Hz  "
          f"mts_fs={sample_rate_hz(mts_time):.2f} Hz  "
          f"peak_i sync/MTS={sync_peak_i}/{mts_peak_i}  "
          f"MTS_peak={info['mts_peak_N']:.1f} N  "
          f"lag {lag:+.2f} s (r={lag_r:.4f})")
    if abs(lag) > 1.0:
        print(f"    [i] the peak anchor was {abs(lag):.1f} s out on this coupon — "
              f"the sync load channel's maximum is not at the true peak instant.")

    # The specimen must be unloaded when the DIC reference frame is taken, or
    # the strain axis starts from an already-strained state while the stress
    # axis does not. A record that opens above the analysis-window threshold has
    # not been registered correctly.
    first = force_N[np.isfinite(force_N)]
    peak = float(np.nanmax(np.abs(mts_force_N)))
    if first.size and peak > 0 and abs(first[0]) > LOAD_START_FRAC * peak:
        print(f"    [!] the first mapped frame already carries {abs(first[0]):.0f} N "
              f"({100 * abs(first[0]) / peak:.1f} % of peak) — the DIC record does "
              f"not contain the\n        unloaded specimen, so E and the strain "
              f"axis for this coupon are not trustworthy.")

    if n < n_l1:
        pad = np.full(n_l1 - n, np.nan)
        force_N = np.concatenate([force_N, pad])
        disp_mm = np.concatenate([disp_mm, pad])
    return force_N, disp_mm, info


def add_coverage(rows):
    """Fill in each row's dic_coverage: its sync-to-MTS peak ratio against the
    batch median, which is the nominal DAQ gain. Below DIC_COVERAGE_MIN the DIC
    recording stopped before the specimen did, so that coupon's peak anchor —
    and therefore its UTS and strain-at-UTS — are not the specimen's."""
    ratios = [r["sync_peak_raw"] / r["mts_peak_N"] for r in rows
              if r["mts_peak_N"] > 0]
    gain = float(np.median(ratios)) if ratios else np.nan
    for r in rows:
        r["dic_coverage"] = ((r["sync_peak_raw"] / r["mts_peak_N"]) / gain
                             if gain > 0 and r["mts_peak_N"] > 0 else np.nan)
    low = [r for r in rows if np.isfinite(r["dic_coverage"])
           and r["dic_coverage"] < DIC_COVERAGE_MIN]
    for r in sorted(low, key=lambda r: r["dic_coverage"]):
        print(f"[!] {r['coupon']}: the sync load channel peaked at "
              f"{r['dic_coverage'] * 100:.1f} % of the gain the rest of the batch")
        print( "    shows — either the DIC record stopped before the specimen "
               "failed, or that channel")
        print( "    is distorted on this coupon. Check the reported lag and this "
               "coupon's UTS.")


# =============================================================================
# SMOOTHING
# =============================================================================
def smooth_signal(x):
    """No-op unless APPLY_SMOOTHING is True."""
    x = np.asarray(x, dtype=float)
    if not APPLY_SMOOTHING:
        return x.copy()
    if FILTER_METHOD == "butterworth":
        return smooth_butterworth(x)
    return smooth_median(x)

def smooth_median(x):
    """Rolling median. Doesn't ring or undershoot a sharp peak the way an
    averaging filter can."""
    win = MEDIAN_WINDOW - 1 if MEDIAN_WINDOW % 2 == 0 else MEDIAN_WINDOW
    if win < 1 or len(x) < win:
        return x.copy()
    return median_filter(x, size=win, mode="constant", cval=np.nan)

def smooth_butterworth(x):
    """Zero-phase low-pass. NaN gaps are bridged before filtering and restored
    after; the output is clipped to the raw range because filtfilt rings past it
    at the truncation edges."""
    n = len(x)
    nan_mask = ~np.isfinite(x)
    if nan_mask.all() or n <= 3 * BUTTER_ORDER:
        return x.copy()
    xi = x.copy()
    if nan_mask.any():
        idx = np.arange(n)
        xi[nan_mask] = np.interp(idx[nan_mask], idx[~nan_mask], x[~nan_mask])
    b, a = butter(BUTTER_ORDER, BUTTER_CUTOFF, btype="low", analog=False)
    out = np.clip(filtfilt(b, a, xi), np.nanmin(xi), np.nanmax(xi))
    out[nan_mask] = np.nan
    return out


# =============================================================================
# SPECIMEN SHEET AND GROUP STATS
# =============================================================================
def read_specimen_csv():
    """The specimen sheet as raw text. dtype=str with keep_default_na=False is
    what makes it safe to write back: cells this script does not touch
    round-trip character for character."""
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc,
                               dtype=str, keep_default_na=False)
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f"{SPECIMEN_CSV.name}: not decodable as "
                       f"{'/'.join(CSV_ENCODINGS)}")

def cell_text(v):
    """One specimen-sheet cell. Missing or non-finite writes blank, so a
    property that could not be computed clears the cell instead of leaving the
    previous run's value standing."""
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"
    if v is None:
        return ""
    v = float(v)
    return "" if not np.isfinite(v) else f"{v:.12g}"

def write_specimen_sheet(rows):
    """Write each coupon's scalars into its row of SPECIMEN_CSV, matched by
    Specimen ID. Missing columns are added at the end; every other cell is
    written back as it was read. Skipped with a warning if the file can't be
    read or replaced — e.g. open in Excel."""
    try:
        df = read_specimen_csv()
    except FileNotFoundError:
        print(f"[!] {SPECIMEN_CSV} not found — skipping specimen sheet update")
        return
    except Exception as exc:
        print(f"[!] {SPECIMEN_CSV.name}: {exc} — skipping specimen sheet update")
        return

    if "Specimen ID" not in df.columns:
        print("[!] no 'Specimen ID' column in the specimen sheet — skipping update")
        return

    for label in SPECIMEN_SHEET_COLUMNS.values():
        if label not in df.columns:
            df[label] = ""
    col_pos = {c: j for j, c in enumerate(df.columns)}
    row_by_id = {cid: i for i, cid in enumerate(df["Specimen ID"])}

    n_written = 0
    for row in rows:
        i = row_by_id.get(row["coupon"])
        if i is None:
            print(f"[!] {row['coupon']} has no row in the specimen sheet — "
                  f"its properties were not written")
            continue
        for key, label in SPECIMEN_SHEET_COLUMNS.items():
            df.iat[i, col_pos[label]] = cell_text(row.get(key))
        n_written += 1

    # Write alongside and rename over: the sheet is the only copy of the
    # hand-entered geometry, so a half-written file would be real data loss.
    tmp = SPECIMEN_CSV.with_name(SPECIMEN_CSV.name + ".tmp")
    try:
        df.to_csv(tmp, index=False, encoding="utf-8-sig")
        os.replace(tmp, SPECIMEN_CSV)
        print(f"Specimen sheet: {n_written} coupon(s) → {SPECIMEN_CSV.name}")
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        print(f"[!] could not update {SPECIMEN_CSV.name} ({exc}) — "
              f"properties were not saved")


def read_existing_group_stats(fp):
    """Read level2_group_stats.csv, upgrading the pre-'test' two-level layout.

    The index depth is READ, not assumed: pandas writes a named MultiIndex as a
    third header line, and reading a two-level file with index_col=[0,1,2] does
    not raise — it swallows the first data column into the index and shifts
    every value one place left. Returns None if the file can't be understood.
    Same routine as FlexuralDIC_Level2.read_existing_group_stats.
    """
    try:
        with open(fp, encoding="utf-8") as fh:
            head = [fh.readline() for _ in range(3)]
        names = [c.strip() for c in head[2].strip("\r\n").split(",") if c.strip()]
        if not names:
            return None
        old = pd.read_csv(fp, header=[0, 1], index_col=list(range(len(names))))
        if "test" not in names:
            old = pd.concat({"tensile": old}, names=["test"])   # pre-'test' file
        old.index = pd.MultiIndex.from_tuples(
            [(str(t), str(e), f"{int(d):02d}" if str(d).strip().isdigit() else str(d))
             for t, e, d in old.index],
            names=["test", "exposure", "direction"])
        return old
    except Exception:
        return None


def write_group_stats(rows):
    """Upsert this run's group statistics into level2_group_stats.csv.
    Only rows whose test == "tensile" are replaced."""
    fp = DIC_DIR / "level2_group_stats.csv"
    df_sum = pd.DataFrame(rows)
    df_sum["test"] = "tensile"
    df_sum["exposure"] = df_sum["coupon"].map(lambda c: parse_id(c)[0])
    df_sum["direction"] = df_sum["coupon"].map(lambda c: parse_id(c)[1])
    agg_cols = ["E_GPa", "E_chord_GPa", "sigma_y_MPa", "UTS_MPa", "eps_at_UTS",
                "poisson_chord"]
    group = (df_sum.groupby(["test", "exposure", "direction"])[agg_cols]
                   .agg(["mean", "std", "count"]))

    if fp.exists():
        old = read_existing_group_stats(fp)
        if old is None:
            print(f"[!] could not read existing {fp.name} — it is being replaced "
                  f"with tensile rows only.\n    Re-run FlexuralDIC_Level2.py to "
                  f"put the flexural rows back.")
        else:
            old = old.drop(index="tensile", level=0, errors="ignore")
            group = pd.concat([old, group]).sort_index()
    group.to_csv(fp)
    return fp


# =============================================================================
# TRUNCATION AND MODULUS
# =============================================================================
def truncation_mask(force_N):
    """Boolean mask selecting the valid test window: pre-load slack and
    post-fracture rebound are False."""
    n = len(force_N)
    peak = float(np.nanmax(np.abs(force_N)))
    if peak <= 0:
        return np.ones(n, dtype=bool)
    i_uts = int(np.nanargmax(np.abs(force_N)))
    # Start on the rising edge, not the first sample over threshold: the
    # pre-touchdown baseline sits near zero with noise either side, so a plain
    # threshold triggers long before the specimen is loaded. On P01 this moves
    # the start from frame 0-12 to frame 41-84.
    f_abs = np.abs(np.nan_to_num(force_N, nan=0.0))
    below = np.flatnonzero(f_abs[:i_uts + 1] < LOAD_START_FRAC * peak)
    i0 = int(below[-1]) + 1 if below.size else 0
    post = np.where(np.abs(force_N[i_uts:]) < LOAD_END_FRAC * peak)[0]
    # Stop one frame BEFORE the first post-UTS drop-off so the failure point
    # itself isn't kept.
    i1 = int(i_uts + post[0]) - 1 if len(post) else n - 1
    mask = np.zeros(n, dtype=bool)
    # max(i0, i1) so a pathological record gives an empty window rather than a
    # reversed slice that selects nothing while looking like it worked.
    mask[i0:max(i0, i1) + 1] = True
    return mask


def fit_modulus(eps, sig, window):
    """Least-squares slope of sig against eps inside `window`.
    Returns (slope MPa, toe offset in strain, n points).

    The toe offset is the fit line's x-intercept — D638 Annex A1. Fitted twice
    because the window is specified in CORRECTED strain but the correction is
    what the first fit produces; without the second pass E is measured over
    [lo - toe, hi - toe] and the toe correction becomes a no-op for E.
    Same routine as FlexuralDIC_Level2.fit_modulus.
    """
    toe = 0.0
    best = (np.nan, np.nan, 0)
    for step in range(2):
        e = eps - toe
        m = (e >= window[0]) & (e <= window[1]) & np.isfinite(e) & np.isfinite(sig)
        n = int(m.sum())
        if n < 3 or not np.any(m):
            # The second pass is a refinement, so losing it must not throw away
            # the first pass's answer. It only happens on a coupon whose toe is
            # comparable to the window itself (P01-TCL45-01).
            if step == 1:
                print(f"    [!] modulus re-fit after toe correction left {n} "
                      f"point(s) in {window} — keeping the first-pass fit")
            break
        slope, icept = np.polyfit(e[m], sig[m], 1)
        if slope == 0:
            break
        toe = toe + (-icept / slope)      # x-intercept, accumulated
        best = (float(slope), float(toe), n)
    return best


def chord_modulus(eps, sig, window):
    """Secant between the two ends of `window`, in already-toe-corrected strain.
    Returns NaN rather than extrapolating when the record does not reach both
    endpoints — np.interp would otherwise clamp and return a confident-looking
    number computed from the wrong point."""
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
    """D638 property extraction from one coupon's truncated arrays.

    Returns a dict of the scalars plus the toe-corrected _eps/_sig/_eps_t
    arrays, a _valid mask saying which input rows survived the finite filter,
    and the pre-smoothing _eps_raw/_sig_raw for Level 3's overlay. None if
    there is too little data.
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
    eps_unsm = (np.asarray(eps_axial_unsmoothed, dtype=float)[valid]
                if eps_axial_unsmoothed is not None else None)
    sig_unsm = (np.asarray(sig_unsmoothed, dtype=float)[valid]
                if sig_unsmoothed is not None else None)

    # Modulus (§11.4) and toe compensation (Annex A1) — one call, because the
    # fit window is defined in toe-corrected strain and the correction comes out
    # of the fit. See fit_modulus.
    E_MPa, eps_offset, n_fit = fit_modulus(eps_raw, sig, MODULUS_STRAIN_RANGE)
    if not np.isfinite(E_MPa):
        return None
    eps = eps_raw - eps_offset
    eps_unsm_corr = (eps_unsm - eps_offset) if eps_unsm is not None else None

    # Transverse strain: subtract its value at the corrected zero of axial strain.
    if np.any(np.isfinite(eps_t_raw)):
        eps_t = eps_t_raw - eps_t_raw[int(np.nanargmin(np.abs(eps)))]
    else:
        eps_t = eps_t_raw.copy()

    # UTS (§11.2 — max stress).
    i_uts = int(np.nanargmax(sig))
    uts, eps_ult = float(sig[i_uts]), float(eps[i_uts])

    # 0.2 % offset yield (§A2.6): first crossing with sigma = E (eps - offset).
    sigma_y = eps_y = np.nan
    diff = sig - E_MPa * (eps - YIELD_OFFSET)
    valid_y = np.where(eps > YIELD_OFFSET)[0]
    if len(valid_y) > 1:
        crossings = np.where(np.diff(np.sign(diff[valid_y])) < 0)[0]
        if len(crossings):
            k = valid_y[crossings[0]]
            denom = diff[k] - diff[k + 1]
            f = diff[k] / denom if denom != 0 else 0.0
            eps_y   = float(eps[k] + f * (eps[k + 1] - eps[k]))
            sigma_y = float(sig[k] + f * (sig[k + 1] - sig[k]))

    # Poisson's ratio (§A3.10.1.3): chord at eps_a = 0.002 over POISSON_RANGE.
    # The least-squares slope (§A3.10.1.1) is reported alongside it.
    nu_chord = nu_slope = np.nan
    pm = ((eps >= POISSON_RANGE[0]) & (eps <= POISSON_RANGE[1])
          & np.isfinite(eps) & np.isfinite(eps_t))
    if pm.sum() >= 3:
        nu_slope = float(-np.polyfit(eps[pm], eps_t[pm], 1)[0])
        order = np.argsort(eps[pm])
        ea, et = eps[pm][order], eps_t[pm][order]
        if ea[0] <= POISSON_CHORD_AT <= ea[-1]:
            nu_chord = float(-np.interp(POISSON_CHORD_AT, ea, et) / POISSON_CHORD_AT)

    # Chord modulus — a cross-check on the tangent value, not a D638 requirement.
    E_chord_MPa = chord_modulus(eps, sig, CHORD_STRAIN_RANGE)
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
        "_eps":     eps,
        "_sig":     sig,
        "_eps_t":   eps_t,
        "_eps_raw": eps_unsm_corr,     # diagnostic overlay only
        "_sig_raw": sig_unsm,
        "_valid":   valid,
    }


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("TensileDIC_Level2 — peak-anchored MTS force, truncate, D638 properties")
    print("=" * 70)
    rows = []
    scalars = load_coupon_scalars()

    for cid in selected_coupons():
        frames_fp = DIC_DIR / f"{cid}.csv"
        if not frames_fp.exists():
            print(f"[{cid}] no per-coupon CSV — run Level 1 first")
            continue
        df = pd.read_csv(frames_fp)
        n = len(df)

        mapped = map_mts_force_disp(cid, n)
        if mapped is None:
            continue
        force_all, disp_mm_mts, info = mapped

        area = (float(scalars.loc[cid, "area_mm2"])
                if cid in scalars.index and pd.notna(scalars.loc[cid, "area_mm2"])
                else np.nan)

        # ---- failure truncation ----
        kept = truncation_mask(force_all)
        kept_idx = np.flatnonzero(kept)
        n_kept = int(kept.sum())

        force_kept  = force_all[kept]
        stress_kept = force_kept / area if np.isfinite(area) else np.full(n_kept, np.nan)
        eps_a_kept  = df["strain_axial_raw"].to_numpy()[kept]
        eps_t_kept  = df["strain_transverse_raw"].to_numpy()[kept]

        # Smoothing pass, with pre-smoothing copies kept alongside for Level 3.
        force_sm  = smooth_signal(force_kept)
        stress_sm = force_sm / area if np.isfinite(area) else np.full(n_kept, np.nan)
        eps_a_sm  = smooth_signal(eps_a_kept)
        eps_t_sm  = smooth_signal(eps_t_kept)

        p = compute_properties(eps_a_sm, stress_sm, eps_t_sm,
                               eps_axial_unsmoothed=eps_a_kept,
                               sig_unsmoothed=stress_kept)
        if not p:
            print(f"[{cid}]  insufficient data")
            continue

        print(f"[{cid}]  E={p['E_GPa']:.2f} GPa  "
              f"E_chord={p['E_chord_GPa']:.2f} GPa  "
              f"σ_y={p['sigma_y_MPa']:.1f} MPa  "
              f"UTS={p['UTS_MPa']:.1f} MPa  "
              f"ε_UTS={p['eps_at_UTS'] * 100:.2f}%  "
              f"ν_chord={p['poisson_chord']:.3f}  "
              f"toe={p['eps_toe'] * 100:.3f}%")

        # ---- scatter the L2 columns back onto the full per-frame file ----
        final_idx = kept_idx[p["_valid"]]

        def scatter(values):
            col = np.full(n, np.nan)
            col[final_idx] = values
            return col

        df["kept"]                    = kept
        df["force_N"]                 = scatter(force_sm[p["_valid"]])
        df["stress_MPa"]              = scatter(p["_sig"])
        df["strain_axial"]            = scatter(p["_eps"])
        df["strain_transverse"]       = scatter(p["_eps_t"])
        df["stress_MPa_unsmoothed"]   = scatter(p["_sig_raw"])
        df["strain_axial_unsmoothed"] = scatter(p["_eps_raw"])
        df["disp_mm_mts"]             = disp_mm_mts
        df.to_csv(frames_fp, index=False, float_format="%.6g")

        rows.append({"coupon": cid,
                     **{k: v for k, v in p.items() if not k.startswith("_")},
                     **info})

    if rows:
        add_coverage(rows)
        write_specimen_sheet(rows)
        write_group_stats(rows)          # D638 §11.7 mean & std per series
        print(f"\n{len(rows)} coupon(s) → DIC/*.csv, {SPECIMEN_CSV.name}, "
              f"DIC/level2_group_stats.csv")

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
