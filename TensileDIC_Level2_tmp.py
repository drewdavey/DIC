#!/usr/bin/env python3
"""
DIC_Level2_tmp.py  —  FSR Tensile Coupons (one-off batch run)
===============================================================
Production Level-2 run for THIS batch only, using peak-force-anchored raw
MTS force/displacement instead of DIC_Level2.py's load_raw*scale approach.

Why this exists: the DIC sync CSV's own "Load" channel for this batch is
only usable as a *shape*, rescaled to the true MTS peak (DIC_Level2.py's
method). Here we instead map the full raw MTS force/displacement time
series onto each DIC frame by treating the DIC-sync peak-force row and the
raw-MTS peak-force row as the same physical instant (they don't share a
clock), then interpolating MTS force/displacement at each DIC frame's time
offset from that anchor. This uses the trustworthy MTS channels for the
whole curve, not just its peak.

This script writes the same production outputs as DIC_Level2.py so
DIC_Level3.py can consume them unchanged. DIC_Level2.py itself is left
untouched — it remains the standard pipeline for future batches.

OUTPUT per coupon (identical contract to DIC_Level2.py, plus one extra
column of its own)
  <DIC_DIR>/<coupon_id>.csv   Level-1's columns, unchanged, plus:
                               kept, force_N, stress_MPa, strain_axial,
                               strain_transverse, stress_MPa_unsmoothed,
                               strain_axial_unsmoothed, disp_mm_mts (the
                               peak-anchored MTS displacement this script
                               derives — kept separate from Level-1's own
                               disp_mm, which stays untouched as the raw
                               sync-CSV reference).
  FSR-SpecimenTesting.csv    scalar properties written into each coupon's row
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
MTS_DIR = DIC_DIR.parent / "MTS"

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

# Raw DIC-sync CSVs (unreliable displacement channel, but their Time column
# and force *shape* are what we use to locate the peak-force row per coupon).
DATA_ROOTS = {
    "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
    "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
}

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "UV": True, "SW": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

APPLY_SMOOTHING = True  # ASTM D638 does not call for filtering the stress-strain record;
                          # toggle on only if a batch's raw signal is genuinely too noisy.
FILTER_METHOD   = "butterworth"  # "median" or "butterworth" — see SMOOTHING section below

# =============================================================================
# CONSTANTS / COLUMNS — DIC-sync raw file
# =============================================================================
HEADERS = 8
IN2MM   = 25.4

DIC_FORCE_SCALED_COL = "LOAD_[kip]_|_CH07_/ai2_scaled"

# =============================================================================
# FAILURE TRUNCATION  — restricts Level-1 data to the valid test window
# before property extraction (mirrors DIC_Level2.py).
# =============================================================================
LOAD_START_FRAC = 0.02
LOAD_END_FRAC   = 0.50

# =============================================================================
# PROPERTY SETTINGS  (mirrors DIC_Level2.py)
# =============================================================================
MODULUS_STRAIN_RANGE = (0.0005, 0.003)
YIELD_OFFSET = 0.002
POISSON_RANGE = (0.0005, 0.0025)
POISSON_CHORD_AT = 0.002

SPECIMEN_SHEET_COLUMNS = {
    "E_GPa":         "E (GPa)",
    "eps_toe":       "Toe Strain",
    "sigma_y_MPa":   "Yield Stress (MPa)",
    "eps_y":         "Yield Strain",
    "UTS_MPa":       "UTS (MPa)",
    "eps_at_UTS":    "Strain at UTS",
    "poisson_chord": "Poisson's Ratio (chord)",
    "poisson_slope": "Poisson's Ratio (slope)",
}

# =============================================================================
# SMOOTHING  (mirrors DIC_Level2.py)
# =============================================================================
MEDIAN_WINDOW = 31

BUTTER_ORDER  = 3
BUTTER_CUTOFF = 0.1   # fraction of Nyquist (0-1) — must be < 1

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
    part = cid.split("-")[1]
    return part[1:-2], part[-2:]

def coupon_dir(cid):
    exposure, _ = parse_id(cid)
    return DATA_ROOTS[exposure] / cid

def find_frames_csv(cid):
    p = DIC_DIR / f"{cid}.csv"
    return p if p.exists() else None

def load_coupon_scalars() -> pd.DataFrame:
    """coupon_scalars.csv, indexed by coupon — only area_mm2 is needed here
    (mts_peak_N isn't: force comes straight from the mapped MTS record)."""
    fp = DIC_DIR / "coupon_scalars.csv"
    if not fp.exists():
        return pd.DataFrame(columns=["mts_peak_N", "area_mm2"])
    return pd.read_csv(fp).set_index("coupon")

def find_mts_txt(cid):
    exact = [MTS_DIR / f"{cid}.txt", MTS_DIR / f"{cid}-TEST.txt"]
    for fp in exact:
        if fp.exists():
            return fp
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    return hits[0] if hits else None

def find_dic_sync_csv(cid):
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
    if not dt.size:
        return np.nan
    return float(1.0 / np.median(dt))

def load_mts_txt(fp):
    return (
        pd.read_csv(
            fp,
            sep="\t",
            skiprows=HEADERS,
            header=None,
            names=["disp_mm", "force_N", "output_V", "time_s"],
            encoding="utf-8-sig",
            on_bad_lines="skip",
        )
        .apply(pd.to_numeric, errors="coerce")
        .dropna(how="all")
    )

def smooth_signal(x):
    x = np.asarray(x, dtype=float)
    if not APPLY_SMOOTHING:
        return x.copy()
    if FILTER_METHOD == "butterworth":
        return _smooth_butterworth(x)
    return _smooth_median(x)

def _smooth_median(x):
    win = MEDIAN_WINDOW
    if win % 2 == 0:
        win -= 1
    if win < 1 or len(x) < win:
        return x.copy()
    return median_filter(x, size=win, mode="constant", cval=np.nan)

def _smooth_butterworth(x):
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
    unchanged. Same routine as TensileDIC_Level2.write_group_stats — this
    variant writes the same tensile rows, from a different alignment path.
    """
    fp = DIC_DIR / "level2_group_stats.csv"
    df_sum = pd.DataFrame(rows)
    df_sum["test"] = "tensile"
    df_sum["exposure"]  = df_sum["coupon"].map(lambda c: parse_id(c)[0])
    df_sum["direction"] = df_sum["coupon"].map(lambda c: parse_id(c)[1])
    agg_cols = ["E_GPa", "sigma_y_MPa", "UTS_MPa", "eps_at_UTS", "poisson_chord"]
    group = (df_sum.groupby(["test", "exposure", "direction"])[agg_cols]
                   .agg(["mean", "std", "count"]))

    if fp.exists():
        try:
            old = pd.read_csv(fp, header=[0, 1], index_col=[0, 1, 2])
            old = old.drop(index="tensile", level=0, errors="ignore")
            group = pd.concat([old, group]).sort_index()
        except Exception as ex:
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
    # Rising edge, not the first sample over threshold — kept byte-identical in
    # behaviour to TensileDIC_Level2.truncation_mask on purpose, since this
    # variant writes to the same outputs. See that function for the reasoning.
    f_abs = np.abs(np.nan_to_num(force_N, nan=0.0))
    below = np.flatnonzero(f_abs[:i_uts + 1] < LOAD_START_FRAC * peak)
    i0 = int(below[-1]) + 1 if below.size else 0
    post = np.where(np.abs(force_N[i_uts:]) < LOAD_END_FRAC * peak)[0]
    i1 = int(i_uts + post[0]) - 1 if len(post) else n - 1
    mask = np.zeros(n, dtype=bool)
    mask[i0:max(i0, i1) + 1] = True
    return mask


# =============================================================================
# PEAK-FORCE-ANCHORED MTS MAPPING
# =============================================================================
def map_mts_force_disp(cid, n_l1):
    """Map raw MTS force/displacement onto DIC-frame times.

    The DIC sync CSV and the raw MTS file don't share a clock, so the
    DIC-sync peak-force row and the raw-MTS peak-force row are treated as
    the same physical instant. Each DIC frame's time offset from that
    anchor is used to interpolate MTS force/displacement at the
    corresponding MTS-clock time. Returns (force_N, disp_mm), each length
    n_l1 (NaN-padded past the sync CSV's length if it's shorter), or None.
    """
    sync_fp = find_dic_sync_csv(cid)
    mts_fp = find_mts_txt(cid)
    if sync_fp is None or mts_fp is None:
        print(f"[{cid}] skip: missing sync={bool(sync_fp)} MTS={bool(mts_fp)}")
        return None

    sync = pd.read_csv(sync_fp)
    time_col = pick_col(sync, "time")
    if DIC_FORCE_SCALED_COL not in sync.columns or time_col is None:
        print(f"[{cid}] skip: sync CSV missing force or time column")
        return None

    n = min(n_l1, len(sync))
    sync = sync.iloc[:n].reset_index(drop=True)
    sync_time = relative_time(sync[time_col], len(sync))
    sync_fs = sample_rate_hz(sync_time)
    sync_force_kip = numeric(sync[DIC_FORCE_SCALED_COL])
    sync_peak_i = int(np.nanargmax(np.abs(sync_force_kip)))

    mts = load_mts_txt(mts_fp)
    mts_time = relative_time(mts["time_s"], len(mts))
    mts_fs = sample_rate_hz(mts_time)
    mts_force_N = numeric(mts["force_N"])
    mts_disp_mm = numeric(mts["disp_mm"])
    mts_peak_i = int(np.nanargmax(np.abs(mts_force_N)))

    dic_dt_from_peak = sync_time - sync_time[sync_peak_i]
    mts_mapped_time = mts_time[mts_peak_i] + dic_dt_from_peak

    mapped_force_N = np.interp(mts_mapped_time, mts_time, mts_force_N, left=np.nan, right=np.nan)
    mapped_disp_mm = np.interp(mts_mapped_time, mts_time, mts_disp_mm, left=np.nan, right=np.nan)
    finite_disp = np.flatnonzero(np.isfinite(mapped_disp_mm))
    if finite_disp.size:
        mapped_disp_mm = mapped_disp_mm - mapped_disp_mm[finite_disp[0]]

    print(f"[{cid}] sync_fs={sync_fs:.3f} Hz  mts_fs={mts_fs:.2f} Hz  "
          f"peak_i sync/MTS={sync_peak_i}/{mts_peak_i}  "
          f"MTS_peak={abs(mts_force_N[mts_peak_i]):.1f} N")

    if n < n_l1:
        pad = np.full(n_l1 - n, np.nan)
        mapped_force_N = np.concatenate([mapped_force_N, pad])
        mapped_disp_mm = np.concatenate([mapped_disp_mm, pad])
    return mapped_force_N, mapped_disp_mm


# =============================================================================
# COMPUTE PROPERTIES  (mirrors DIC_Level2.py)
# =============================================================================
def compute_properties(eps_axial, sig, eps_transverse,
                       eps_axial_unsmoothed=None, sig_unsmoothed=None):
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

    lo, hi = MODULUS_STRAIN_RANGE
    mfit = (eps_raw >= lo) & (eps_raw <= hi) & np.isfinite(eps_raw) & np.isfinite(sig)
    if mfit.sum() < 3:
        return None
    slope, intercept = np.polyfit(eps_raw[mfit], sig[mfit], 1)
    E_MPa = float(slope)

    eps_offset = -intercept / E_MPa if E_MPa != 0 else 0.0
    eps   = eps_raw   - eps_offset
    eps_unsmoothed_corr = (eps_unsmoothed - eps_offset) if eps_unsmoothed is not None else None
    if np.any(np.isfinite(eps_t_raw)):
        i0 = int(np.nanargmin(np.abs(eps)))
        eps_t = eps_t_raw - eps_t_raw[i0]
    else:
        eps_t = eps_t_raw.copy()

    i_uts   = int(np.nanargmax(sig))
    uts     = float(sig[i_uts])
    eps_ult = float(eps[i_uts])

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

    nu_chord = nu_slope = np.nan
    pm = ((eps >= POISSON_RANGE[0]) & (eps <= POISSON_RANGE[1]) &
          np.isfinite(eps) & np.isfinite(eps_t))
    if pm.sum() >= 3:
        nu_slope = float(-np.polyfit(eps[pm], eps_t[pm], 1)[0])
        order = np.argsort(eps[pm])
        ea, et = eps[pm][order], eps_t[pm][order]
        if ea[0] <= POISSON_CHORD_AT <= ea[-1]:
            nu_chord = float(-np.interp(POISSON_CHORD_AT, ea, et) / POISSON_CHORD_AT)

    return {
        "E_GPa":         E_MPa / 1000.0,
        "eps_toe":       eps_offset,
        "sigma_y_MPa":   sigma_y,
        "eps_y":         eps_y,
        "UTS_MPa":       uts,
        "eps_at_UTS":    eps_ult,
        "poisson_chord": nu_chord,
        "poisson_slope": nu_slope,
        "_eps":          eps,
        "_sig":          sig,
        "_eps_t":        eps_t,
        "_eps_raw":      eps_unsmoothed_corr,
        "_sig_raw":      sig_unsmoothed,
        "_valid":        valid,
    }


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("DIC_Level2_tmp — batch Level-2 run using peak-anchored MTS force mapping")
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

        mapped = map_mts_force_disp(cid, n)
        if mapped is None:
            continue
        force_N_all, disp_mm_mts = mapped

        area = (float(scalars.loc[cid, "area_mm2"])
                if cid in scalars.index and pd.notna(scalars.loc[cid, "area_mm2"])
                else np.nan)

        # ---- failure truncation ---------------------------------------------
        kept = truncation_mask(force_N_all)
        kept_idx = np.flatnonzero(kept)

        force_kept  = force_N_all[kept]
        stress_kept = force_kept / area if np.isfinite(area) else np.full(kept.sum(), np.nan)
        eps_a_kept  = df["strain_axial_raw"].to_numpy()[kept]
        eps_t_kept  = df["strain_transverse_raw"].to_numpy()[kept]

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
        df["disp_mm_mts"]              = disp_mm_mts
        df.to_csv(frames_fp, index=False, float_format="%.6g")

        rows.append({"coupon": cid, **{k: v for k, v in p.items() if not k.startswith("_")}})

    if rows:
        write_specimen_sheet(rows)

        write_group_stats(rows)

        print(f"\n{len(rows)} coupon(s) → DIC/*.csv, {SPECIMEN_CSV.name}, "
              f"DIC/level2_group_stats.csv")

    print(f"\nDone. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
