#!/usr/bin/env python3
"""
DIC_Level1.py  —  FSR Tensile Coupons
======================================
Step A — exports each coupon's raw VIC-3D .out files to per-frame CSVs
(next to the .out files, on the raw data drive).
Step B — pairs those per-frame CSVs with the MTS force/displacement record,
inserts virtual axial + transverse extensometers, and writes one compact
per-coupon CSV (full, untruncated record) to DIC_DIR.
Step C — signal-inspection plots + a cross-check report comparing the raw
MTS and raw DIC-sync force/displacement channels directly (no Level-1
computation involved; a sanity check on the inputs themselves).

No truncation or load scaling happens here — that's Level-2's job, so it
can be tuned without re-running Step A (slow, one .out load per frame) or
Step B (needs every per-frame CSV, slower than Level-2/3).

USAGE
-----
1. Toggle the SWITCHES section below to pick which coupons to run.
2. Make sure vicpyx is installed in the active Python environment:
       pip install vicpyx
3. Run:
       python DIC_Level1.py

INPUT per coupon
  <coupon_dir>/*.out                    VIC-3D full-field export, one per DIC frame
  <coupon_dir>/<coupon_id>.csv           VIC sync CSV (analog channels @ DIC frame rate)
  <MTS_DIR>/<coupon_id>*.txt             MTS raw file: cols disp_mm, force_N, output_V, time_s
  FSR-SpecimenTesting.xlsx               gauge thickness × width  →  area

OUTPUTS
-------
- <coupon_dir>/<out_filename>.csv        one CSV per .out file, written next to it
- <DIC_DIR>/<coupon_id>_L1.csv           step, time_s, load_raw, disp_mm,
                                         strain_axial, strain_transverse,
                                         mts_peak_N, area_mm2
- <FIGS_ROOT>/<coupon_id>/MTS_force_disp.png   raw MTS curve (sanity check)
- <FIGS_ROOT>/<coupon_id>/MTS_force_displacement_signals.png       raw MTS force/disp vs. time
- <FIGS_ROOT>/<coupon_id>/DIC_sync_force_displacement_signals.png  raw vs. scaled DIC-sync channels
- <DIC_DIR>/raw_dic_force_displacement_signal_report.csv           per-coupon raw-signal cross-check

UNIT NOTES
  MTS .txt is already in mm / N / V / sec (verified against tensile_analysis.py).
  VIC sync CSV "Load" column units are device-dependent; rather than guess,
  Level-2 derives a per-coupon scale factor from mts_peak_N / max(load_raw).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vicpyx import VicDataSet as VICDataSet

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# PATHS
# =============================================================================
MTS_DIR = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\MTS"
)
DIC_DIR = MTS_DIR.parent / "DIC"   # consolidated _L1.csv files land here
# Raw VIC-3D project data (.out/.tif per frame + per-coupon sync CSV), mirrored
# locally under DIC_DIR/raw — G:\DrewDavey\... was the old external-drive location.
DATA_ROOTS = {
    "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
    "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
}
SPECIMEN_SHEET = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.xlsx"
)
FIGS_ROOT = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\figs"
)

# =============================================================================
# SWITCHES — toggle which coupons to process
# =============================================================================
PRINTS = ["P01"]
EXPOSURES = {"CL": True, "SW": True, "UV": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

# Pipeline steps. Steps A and B are each independently skippable/overwritable
# so the slow .out export never has to re-run just to rebuild the
# consolidated CSV. Step C is cheap (plots + a report CSV, no vicpyx calls)
# and always regenerates its outputs when enabled — no overwrite switch needed.
DO_EXPORT_FRAMES = True     # Step A: export each .out to a CSV next to it
OVERWRITE_FRAMES = False    # if False, skip .out files whose .csv already exists
DO_BUILD_L1      = True     # Step B: pair frames + MTS, build extensometers
OVERWRITE_L1     = False    # if False, skip coupons whose _L1.csv already exists
DO_SIGNAL_PLOTS  = True     # Step C: raw MTS + DIC-sync signal-inspection plots + report CSV

# If True, displacement signals in the Step C plots are shifted so the first
# finite value is zero. Doesn't affect _L1.csv (disp_mm there is zeroed separately).
ZERO_DISPLACEMENT = True

# Variables to export from each .out (Step A)
# Standard full-field DIC variables; sigma is needed to filter invalid points.
EXPORT_VARS = [
    "sigma",                              # correlation confidence (filter on this)
    "X", "Y", "Z",                        # world coords (mm)
    "U", "V", "W",                        # displacements (mm)
    "exx", "eyy", "exy",                  # in-plane strains
    "e1", "e2", "gamma",                  # principal & max-shear strains
    "x", "y", "u", "v",                   # pixel coords / pixel disps
    "q", "r", "q_ref", "r_ref",           # subset coords
]

# =============================================================================
# VIRTUAL EXTENSOMETER  — ASTM D638 §5.2.1 (Class B-2 equivalent for modulus)
# Gauge length 50 mm (2 in) per D638 Type I Fig. 1.
# =============================================================================
AXIAL_GAUGE_IN = 4.36       # axial (Y, loading) gauge length, inches  [D638 G = 2.00 in]
TRANS_GAUGE_IN = 1.0        # transverse (X) gauge length, inches      [Annex A3.5.2]

# =============================================================================
# CONSTANTS
# =============================================================================
IN2MM   = 25.4
HEADERS = 8               # MTS .txt has an 8-line header (verified)

# DIC sync CSV column names (Step C signal-inspection plots + report only —
# Step B above reads the same load/drift columns generically via pick_col()).
DIC_FORCE_RAW_COL    = "Dev1/ai2"
DIC_FORCE_SCALED_COL = "LOAD_[kip]_|_CH07_/ai2_scaled"
DIC_DISP_RAW_COL     = "Dev1/ai1"
DIC_DISP_SCALED_COL  = "DRIFT_[in]|_CH06/ai1_scaled"


# =============================================================================
# HELPERS — shared
# =============================================================================
def coupon_id(p, e, d, r): return f"{p}-T{e}{d}-{r}"

def selected_coupons():
    return [coupon_id(p, e, d, r)
            for p in PRINTS
            for e, on in EXPOSURES.items() if on
            for d, on2 in DIRECTIONS.items() if on2
            for r in REPLICATES]

def coupon_dir(cid):
    exp = cid.split("-")[1][1:-2]
    return DATA_ROOTS[exp] / cid

def find_out_files(cdir):
    """All .out files inside a coupon's project directory (searched recursively
    since VIC-3D's exact output location can vary by project structure)."""
    return sorted(cdir.rglob("*.out"))

def find_first(paths):
    for p in paths:
        if p.exists():
            return p
    return None

def find_mts_txt(cid):
    p = find_first([MTS_DIR / f"{cid}.txt", MTS_DIR / f"{cid}-TEST.txt"])
    if p is None:
        hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
        p = hits[0] if hits else None
    return p

def find_sync_csv(cdir, cid):
    return find_first([cdir / f"{cid}.csv"])

def find_frame_csvs(cdir, cid):
    return sorted(cdir.glob(f"{cid}-????????_0.csv"))

def pick_col(df, hint):
    for c in df.columns:
        if hint.lower() in c.lower():
            return c
    return None

def load_mts_txt(fp):
    """MTS .txt: 8-line header, then tab-separated cols disp_mm, force_N, output_V, time_s."""
    return (pd.read_csv(fp, sep="\t", skiprows=HEADERS, header=None,
                        names=["disp_mm", "force_N", "output_V", "time_s"],
                        encoding="utf-8-sig", on_bad_lines="skip")
              .apply(pd.to_numeric, errors="coerce")
              .dropna(subset=["force_N"]))

def load_specimen_sheet():
    df = pd.read_excel(SPECIMEN_SHEET)
    t_col = next((c for c in df.columns if "thickness" in c.lower()), None)
    w_col = next((c for c in df.columns if "width" in c.lower() and "dia" in c.lower()), None)
    if t_col is None or w_col is None:
        raise RuntimeError(f"could not find thickness/width cols in {SPECIMEN_SHEET}")
    df = df.rename(columns={t_col: "t_in", w_col: "w_in"})
    return df.set_index("Specimen ID")

def get_area_mm2(spec, cid):
    """ASTM D638 §11.2: stress uses *original* cross-sectional area."""
    if cid not in spec.index:
        return None
    row = spec.loc[cid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    return float(row["t_in"]) * float(row["w_in"]) * IN2MM * IN2MM


# =============================================================================
# STEP A — .out → per-frame CSV  (vicpyx export)
# =============================================================================
def export_out_to_csv(out_path: Path, csv_path: Path,
                      var_names: Iterable[str]) -> bool:
    """Convert a single .out file to CSV using vicpyx.

    Writes one row per valid AOI data point with columns named after
    var_names (sigma column dropped from the output, only used to
    filter invalid points). Returns True on success, False otherwise.
    """
    ds = VICDataSet()
    try:
        ds.load(str(out_path))
    except Exception as ex:
        print(f"    [!] could not load {out_path.name}: {ex}")
        return False

    try:
        available = list(ds.variables())
    except Exception:
        available = []
    wanted = [v for v in var_names if (not available) or (v in available)]
    if not wanted:
        print(f"    [!] none of the requested variables found in {out_path.name}")
        return False
    try:
        values = ds.get_values(wanted)  # numpy structured array
    except Exception as ex:
        print(f"    [!] get_values failed on {out_path.name}: {ex}")
        return False

    mask = np.ones(len(values), dtype=bool)
    if "sigma" in values.dtype.names:
        mask &= values["sigma"] >= 0

    export_cols = [v for v in wanted if v != "sigma"]
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write(",".join(export_cols) + "\n")
        arrs = [np.asarray(values[v])[mask] for v in export_cols]
        for row in zip(*arrs):
            fh.write(",".join(f"{x}" for x in row) + "\n")
    return True


def export_frames(cid, cdir) -> dict:
    """Step A for one coupon: export every .out to a sibling CSV."""
    stats = {"n_out": 0, "n_csv_written": 0, "n_csv_skipped": 0}
    out_files = find_out_files(cdir)
    stats["n_out"] = len(out_files)
    if not out_files:
        print(f"  [skip] no .out files in {cdir}")
        return stats

    for out_path in out_files:
        csv_path = out_path.with_suffix(".csv")
        if csv_path.exists() and not OVERWRITE_FRAMES:
            stats["n_csv_skipped"] += 1
            continue
        if export_out_to_csv(out_path, csv_path, EXPORT_VARS):
            stats["n_csv_written"] += 1

    print(f"  Step A: {stats['n_out']} .out files, "
          f"wrote {stats['n_csv_written']}, skipped {stats['n_csv_skipped']} (already existed)")
    return stats


# =============================================================================
# POINT EXTENSOMETER  — mirrors VIC-3D InspectorItemSet.add_extensometer()
# Two fixed endpoints; nearest DIC point to each (mirrors at_global_xy);
# =============================================================================
def ext_endpoints(frame_csv0, axial_mm, trans_mm):
    """
    Compute the four endpoint world-coordinate positions from the AOI centroid
    of the reference frame.  Returns (Xc, Y_bot, Y_top, X_lft, X_rgt, Yc).
    """
    ref  = pd.read_csv(frame_csv0).dropna(subset=["X", "Y"])
    Yc   = float(ref["Y"].median())
    Xc   = float(ref["X"].median())
    Ymin, Ymax = float(ref["Y"].min()), float(ref["Y"].max())
    Xmin, Xmax = float(ref["X"].min()), float(ref["X"].max())
    Y_top = min(Yc + axial_mm / 2, Ymax)
    Y_bot = max(Yc - axial_mm / 2, Ymin)
    X_rgt = min(Xc + trans_mm / 2, Xmax)
    X_lft = max(Xc - trans_mm / 2, Xmin)
    return Xc, Y_bot, Y_top, X_lft, X_rgt, Yc


def point_extensometer(frame_csvs, x0, y0, x1, y1):
    """
    VIC-3D style extensometer: two fixed endpoint markers at (x0,y0) and (x1,y1).
    For every frame find the nearest DIC point to each endpoint (mirrors
    VICDataSet.at_global_xy), read its (U, V) displacement, and compute
    engineering strain from the change in distance between the displaced markers.

        ε = (L_deformed − L₀) / L₀
        L₀ = √((x1−x0)² + (y1−y0)²)
    """
    L0 = float(np.sqrt((x1 - x0)**2 + (y1 - y0)**2))
    if L0 == 0:
        return np.full(len(frame_csvs), np.nan)

    eps = []
    for fp in frame_csvs:
        try:
            df = pd.read_csv(fp).dropna(subset=["X", "Y", "U", "V"])
        except Exception:
            eps.append(np.nan); continue
        if len(df) < 2:
            eps.append(np.nan); continue

        # Nearest DIC point to each endpoint (at_global_xy equivalent)
        r0 = df.loc[((df["X"] - x0)**2 + (df["Y"] - y0)**2).idxmin()]
        r1 = df.loc[((df["X"] - x1)**2 + (df["Y"] - y1)**2).idxmin()]

        dx = (x1 + float(r1["U"])) - (x0 + float(r0["U"]))
        dy = (y1 + float(r1["V"])) - (y0 + float(r0["V"]))
        eps.append((np.sqrt(dx**2 + dy**2) - L0) / L0)

    return np.array(eps)


# =============================================================================
# RAW MTS PLOT (sanity check; displacement zeroed to start of test)
# =============================================================================
def plot_mts(cid, mts):
    fig_dir = FIGS_ROOT / cid
    fig_dir.mkdir(parents=True, exist_ok=True)
    disp_rel = mts["disp_mm"].to_numpy() - float(mts["disp_mm"].iloc[0])  # relative
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(disp_rel, mts["force_N"]/1000.0, lw=1.0)
    ax.set_xlabel("Crosshead Displacement (mm, relative)")
    ax.set_ylabel("Force (kN)")
    ax.set_title(f"{cid}  —  raw MTS")
    ax.grid(alpha=0.3, linestyle="--")
    fig.tight_layout()
    out = fig_dir / "MTS_force_disp.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# =============================================================================
# STEP B — pair frames + MTS, build extensometers, write _L1.csv
# =============================================================================
def build_l1(cid, cdir, spec) -> bool:
    """Step B for one coupon. Returns True if _L1.csv was written."""
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    out_fp = DIC_DIR / f"{cid}_L1.csv"
    if out_fp.exists() and not OVERWRITE_L1:
        print(f"  Step B: [skip] _L1.csv exists")
        return False

    sync_fp   = find_sync_csv(cdir, cid)
    frame_fps = find_frame_csvs(cdir, cid)
    mts_fp    = find_mts_txt(cid)
    area_mm2  = get_area_mm2(spec, cid)

    if sync_fp is None:    print(f"  Step B: [skip] no sync CSV {cid}.csv");           return False
    if not frame_fps:      print(f"  Step B: [skip] no per-frame CSVs (run Step A first)"); return False
    if mts_fp is None:     print(f"  Step B: [skip] no MTS .txt — needed for mts_peak_N"); return False
    if area_mm2 is None:   print(f"  Step B: [warn] no area for {cid} — stress NaN")

    print(f"  Step B: {len(frame_fps)} frames | area = "
          + (f"{area_mm2:.2f} mm²" if area_mm2 else "N/A"))

    # ---- raw MTS load ----
    mts = load_mts_txt(mts_fp)
    plot_mts(cid, mts)
    peak_force_N = float(mts["force_N"].abs().max())
    print(f"  MTS peak force: {peak_force_N/1000:.2f} kN")

    # ---- sync CSV: read raw load signal (unscaled) ----
    # Scaling is deferred to Level-2: smooth(load_raw) → scale → force_N.
    sync = pd.read_csv(sync_fp)
    load_col = pick_col(sync, "load")
    disp_col = pick_col(sync, "drift")
    time_col = pick_col(sync, "time")
    if load_col is None:
        print(f"  Step B: [skip] no LOAD column. Cols: {list(sync.columns)}"); return False

    load_raw = pd.to_numeric(sync[load_col], errors="coerce").to_numpy()
    raw_peak = float(np.nanmax(np.abs(load_raw)))
    if raw_peak <= 0 or not np.isfinite(raw_peak):
        print(f"  Step B: [skip] sync load peak invalid"); return False
    print(f"  sync raw peak: {raw_peak:.4f} units  |  MTS peak: {peak_force_N:.0f} N  "
          f"(implied scale {peak_force_N/raw_peak:.4f})")

    # ---- displacement: from sync CSV "Drift" column, zeroed to start ----
    if disp_col:
        disp_raw = pd.to_numeric(sync[disp_col], errors="coerce").to_numpy()
        disp_mm = disp_raw * IN2MM
        disp_mm = disp_mm - disp_mm[0]
    else:
        disp_mm = np.full_like(load_raw, np.nan)

    time_s = (pd.to_numeric(sync[time_col], errors="coerce").to_numpy()
              if time_col else np.arange(len(load_raw), dtype=float))

    # ---- align to DIC frame count ----
    n = min(len(load_raw), len(frame_fps))
    load_raw, disp_mm, time_s = load_raw[:n], disp_mm[:n], time_s[:n]
    frame_fps_used = frame_fps[:n]

    # ---- point extensometers: E0 axial, E1 transverse (D638 §5.2 / Annex A3) ----
    axial_mm = AXIAL_GAUGE_IN * IN2MM
    trans_mm = TRANS_GAUGE_IN * IN2MM
    Xc, Y_bot, Y_top, X_lft, X_rgt, Yc = ext_endpoints(
        frame_fps_used[0], axial_mm, trans_mm)
    print(f"    E0 axial     (Xc={Xc:.1f})  Y: {Y_bot:.1f} → {Y_top:.1f}  "
          f"L={Y_top-Y_bot:.1f} mm")
    print(f"    E1 transverse (Yc={Yc:.1f})  X: {X_lft:.1f} → {X_rgt:.1f}  "
          f"L={X_rgt-X_lft:.1f} mm")
    eps_a = point_extensometer(frame_fps_used, Xc, Y_bot, Xc, Y_top)
    eps_t = point_extensometer(frame_fps_used, X_lft, Yc, X_rgt, Yc)

    # ---- write ----
    # force_N and stress_MPa are NOT saved — Level-2 computes them after
    # smoothing and per-coupon scaling. load_raw, mts_peak_N, area_mm2 provide
    # everything Level-2 needs without re-running this slower step.
    # strain stays "raw" — toe compensation per D638 Annex A1 done in Level-2.
    pd.DataFrame({
        "step":              np.arange(n),
        "time_s":            time_s,
        "load_raw":          load_raw,
        "disp_mm":           disp_mm,
        "strain_axial":      eps_a,
        "strain_transverse": eps_t,
        "mts_peak_N":        peak_force_N,
        "area_mm2":          area_mm2 if area_mm2 else np.nan,
    }).to_csv(out_fp, index=False, float_format="%.6g")
    print(f"  → DIC/{out_fp.name} ({n} rows)")
    return True


# =============================================================================
# STEP C — raw signal-inspection plots + cross-check report
# Independent of Steps A/B (reads only the raw MTS .txt and raw DIC sync
# CSV) — a sanity check on the input signals themselves, not on anything
# Level-1 computes.
# =============================================================================
def first_finite_zeroed(values) -> np.ndarray:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not ZERO_DISPLACEMENT:
        return arr
    finite = np.flatnonzero(np.isfinite(arr))
    if finite.size:
        arr = arr - arr[finite[0]]
    return arr

def relative_time(values, fallback_len) -> np.ndarray:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = np.flatnonzero(np.isfinite(arr))
    if finite.size:
        return arr - arr[finite[0]]
    return np.arange(fallback_len, dtype=float)

def numeric_col(df, col) -> np.ndarray:
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)

def linear_scale(raw, scaled) -> tuple[float, float]:
    """Return slope and intercept for scaled ~= slope * raw + intercept."""
    raw = np.asarray(raw, dtype=float)
    scaled = np.asarray(scaled, dtype=float)
    mask = np.isfinite(raw) & np.isfinite(scaled)
    if mask.sum() < 2:
        return np.nan, np.nan
    slope, intercept = np.polyfit(raw[mask], scaled[mask], 1)
    return float(slope), float(intercept)

def dominant_frequency_hz(values, time_s) -> float:
    """Dominant nonzero FFT frequency for a raw signal, in Hz."""
    values = np.asarray(values, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    mask = np.isfinite(values) & np.isfinite(time_s)
    if mask.sum() < 4:
        return np.nan

    y = values[mask]
    t = time_s[mask]
    dt = np.nanmedian(np.diff(t))
    if not np.isfinite(dt) or dt <= 0:
        return np.nan

    y = y - np.nanmean(y)
    if not np.any(np.isfinite(y)) or np.nanmax(np.abs(y)) == 0:
        return np.nan

    amp = np.abs(np.fft.rfft(y))
    freq = np.fft.rfftfreq(len(y), d=dt)
    if len(freq) <= 1:
        return np.nan
    peak_i = int(np.nanargmax(amp[1:]) + 1)
    return float(freq[peak_i])

def mts_peak_disp_zeroed(cid) -> tuple[float, float, float]:
    """Return MTS peak force, zeroed displacement at peak force, and zero offset."""
    fp = find_mts_txt(cid)
    if fp is None:
        return np.nan, np.nan, np.nan
    df = load_mts_txt(fp)
    if df.empty or "force_N" not in df or "disp_mm" not in df:
        return np.nan, np.nan, np.nan

    force_N = pd.to_numeric(df["force_N"], errors="coerce").to_numpy(dtype=float)
    disp_mm = pd.to_numeric(df["disp_mm"], errors="coerce").to_numpy(dtype=float)
    finite_force = np.flatnonzero(np.isfinite(force_N))
    finite_disp = np.flatnonzero(np.isfinite(disp_mm))
    if not finite_force.size or not finite_disp.size:
        return np.nan, np.nan, np.nan

    peak_i = int(finite_force[int(np.nanargmax(np.abs(force_N[finite_force])))])
    zero_factor = float(disp_mm[finite_disp[0]])
    peak_force = float(abs(force_N[peak_i]))
    peak_disp_zeroed = (float(disp_mm[peak_i] - zero_factor)
                        if np.isfinite(disp_mm[peak_i]) else np.nan)
    return peak_force, peak_disp_zeroed, zero_factor

def require_columns(df, cols, fp) -> bool:
    missing = [col for col in cols if col not in df.columns]
    if missing:
        print(f"  [skip] {fp.name} missing columns: {missing}")
        return False
    return True

def save_signal_plot(cid, source_label, time_s, force_or_load, force_label,
                     disp_mm, out_fp) -> None:
    fig, (ax_force, ax_disp) = plt.subplots(2, 1, figsize=(8, 5.8), sharex=True)

    ax_force.plot(time_s, force_or_load, lw=0.9, color="#2457a6")
    ax_force.set_ylabel(force_label)
    ax_force.grid(alpha=0.25, linestyle="--")

    if disp_mm is None or np.all(~np.isfinite(disp_mm)):
        ax_disp.text(0.5, 0.5, "No displacement column found",
                     transform=ax_disp.transAxes, ha="center", va="center", color="0.35")
    else:
        ax_disp.plot(time_s, disp_mm, lw=0.9, color="#b3412c")
    ax_disp.set_ylabel("Displacement (mm)")
    ax_disp.set_xlabel("Time from first sample (s)")
    ax_disp.grid(alpha=0.25, linestyle="--")

    fig.suptitle(f"{cid} - {source_label}")
    fig.tight_layout()
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

def save_dic_twin_axis_plot(cid, time_s, force_raw_v, force_scaled_kip,
                            disp_raw_v, disp_scaled_in, out_fp) -> None:
    fig, (ax_force_raw, ax_disp_raw) = plt.subplots(2, 1, figsize=(8.8, 6.2), sharex=True)

    ax_force_scaled = ax_force_raw.twinx()
    line_f_raw, = ax_force_raw.plot(time_s, force_raw_v, lw=0.9, color="#2457a6",
                                    label=DIC_FORCE_RAW_COL)
    line_f_scaled, = ax_force_scaled.plot(time_s, force_scaled_kip, lw=1.0, linestyle="--",
                                          color="#b3412c", alpha=0.8, label=DIC_FORCE_SCALED_COL)
    ax_force_raw.set_ylabel("Dev1/ai2 (V)", color="#2457a6")
    ax_force_scaled.set_ylabel("Load (kip)", color="#b3412c")
    ax_force_raw.tick_params(axis="y", labelcolor="#2457a6")
    ax_force_scaled.tick_params(axis="y", labelcolor="#b3412c")
    ax_force_raw.grid(alpha=0.25, linestyle="--")
    ax_force_raw.legend([line_f_raw, line_f_scaled],
                        [line_f_raw.get_label(), line_f_scaled.get_label()],
                        fontsize=8, loc="best")

    ax_disp_scaled = ax_disp_raw.twinx()
    line_d_raw, = ax_disp_raw.plot(time_s, disp_raw_v, lw=0.9, color="#2457a6",
                                   label=DIC_DISP_RAW_COL)
    line_d_scaled, = ax_disp_scaled.plot(time_s, disp_scaled_in, lw=1.0, linestyle="--",
                                         color="#b3412c", alpha=0.8, label=DIC_DISP_SCALED_COL)
    ax_disp_raw.set_ylabel("Dev1/ai1 (V)", color="#2457a6")
    ax_disp_scaled.set_ylabel("Drift (in)", color="#b3412c")
    ax_disp_raw.tick_params(axis="y", labelcolor="#2457a6")
    ax_disp_scaled.tick_params(axis="y", labelcolor="#b3412c")
    ax_disp_raw.set_xlabel("Time from first sample (s)")
    ax_disp_raw.grid(alpha=0.25, linestyle="--")
    ax_disp_raw.legend([line_d_raw, line_d_scaled],
                       [line_d_raw.get_label(), line_d_scaled.get_label()],
                       fontsize=8, loc="best")

    fig.suptitle(f"{cid} - DIC raw sync CSV")
    fig.tight_layout()
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_mts_signals(cid) -> Path | None:
    fp = find_mts_txt(cid)
    if fp is None:
        print(f"[{cid}] MTS: skip, no .txt found")
        return None

    df = load_mts_txt(fp)
    time_s = relative_time(df["time_s"], len(df))
    force_N = pd.to_numeric(df["force_N"], errors="coerce").to_numpy(dtype=float)
    disp_mm = first_finite_zeroed(df["disp_mm"])

    out_fp = FIGS_ROOT / cid / "MTS_force_displacement_signals.png"
    save_signal_plot(cid, "MTS raw file", time_s, force_N, "Force (N)", disp_mm, out_fp)
    print(f"[{cid}] MTS: {out_fp}")
    return out_fp

def plot_dic_sync_signals(cid, cdir) -> Path | None:
    fp = find_sync_csv(cdir, cid)
    if fp is None:
        print(f"[{cid}] DIC sync: skip, no raw per-coupon CSV found")
        return None

    df = pd.read_csv(fp)
    time_col = pick_col(df, "time")
    required = [DIC_FORCE_RAW_COL, DIC_FORCE_SCALED_COL, DIC_DISP_RAW_COL, DIC_DISP_SCALED_COL]
    if not require_columns(df, required, fp):
        return None

    time_s = relative_time(df[time_col], len(df)) if time_col else np.arange(len(df), dtype=float)

    out_fp = FIGS_ROOT / cid / "DIC_sync_force_displacement_signals.png"
    save_dic_twin_axis_plot(cid, time_s,
                            numeric_col(df, DIC_FORCE_RAW_COL),
                            numeric_col(df, DIC_FORCE_SCALED_COL),
                            numeric_col(df, DIC_DISP_RAW_COL),
                            numeric_col(df, DIC_DISP_SCALED_COL),
                            out_fp)
    print(f"[{cid}] DIC sync: {out_fp}")
    return out_fp

def dic_report_row(cid, cdir) -> dict:
    row = {
        "coupon": cid,  # Coupon/specimen ID.
        "dic_disp_raw2scaled": np.nan,  # Slope from Dev1/ai1 volts to scaled drift inches.
        "dic_disp_scaled2mm": IN2MM,  # Converts the scaled DIC displacement from inches to mm.
        "dic_disp_zeroed": np.nan,  # DIC displacement at peak-force/break index after zeroing, in mm.
        "mts_disp_zeroed": np.nan,  # MTS displacement at MTS peak-force index after zeroing, in mm.
        "mts_disp_zeroFactor": np.nan,  # First MTS displacement value in mm; subtract this to zero the trace.
        "dic_disp_zeroFactor": np.nan,  # First scaled displacement value in mm; subtract this to zero the trace.
        "dic_force_raw2scaled": np.nan,  # Slope from Dev1/ai2 volts to scaled load kips.
        "dic_force_scaled2mts": np.nan,  # DIC peak force after scaling to match the MTS peak force, in N.
        "mts_raw_peak": np.nan,  # Peak absolute force from the matching raw MTS .txt file, in N.
        "scale_dic2mts_peaks": np.nan,  # Per-coupon scale: MTS peak force (N) / DIC scaled peak load (kip).
        "dic_disp_raw_peakFreq": np.nan,  # Dominant nonzero FFT frequency of raw Dev1/ai1 displacement signal.
        "dic_force_raw_peakFreq": np.nan,  # Dominant nonzero FFT frequency of raw Dev1/ai2 force signal.
    }

    fp = find_sync_csv(cdir, cid)
    if fp is None:
        return row

    df = pd.read_csv(fp)
    time_col = pick_col(df, "time")
    required = [DIC_FORCE_RAW_COL, DIC_FORCE_SCALED_COL, DIC_DISP_RAW_COL, DIC_DISP_SCALED_COL]
    if any(col not in df.columns for col in required):
        return row

    force_raw_v = numeric_col(df, DIC_FORCE_RAW_COL)
    force_scaled_kip = numeric_col(df, DIC_FORCE_SCALED_COL)
    disp_raw_v = numeric_col(df, DIC_DISP_RAW_COL)
    disp_scaled_in = numeric_col(df, DIC_DISP_SCALED_COL)
    time_s = relative_time(df[time_col], len(df)) if time_col else np.arange(len(df), dtype=float)

    disp_scale, _ = linear_scale(disp_raw_v, disp_scaled_in)
    force_scale, _ = linear_scale(force_raw_v, force_scaled_kip)
    row["dic_disp_raw2scaled"] = disp_scale
    row["dic_force_raw2scaled"] = force_scale
    row["dic_disp_raw_peakFreq"] = dominant_frequency_hz(disp_raw_v, time_s)
    row["dic_force_raw_peakFreq"] = dominant_frequency_hz(force_raw_v, time_s)
    peak_N, mts_disp_zeroed, mts_disp_zeroFactor = mts_peak_disp_zeroed(cid)
    row["mts_raw_peak"] = peak_N
    row["mts_disp_zeroed"] = mts_disp_zeroed
    row["mts_disp_zeroFactor"] = mts_disp_zeroFactor

    finite_disp = np.flatnonzero(np.isfinite(disp_scaled_in))
    disp_zero_mm = np.nan
    if finite_disp.size:
        disp_zero_mm = float(disp_scaled_in[finite_disp[0]] * IN2MM)
        row["dic_disp_zeroFactor"] = disp_zero_mm

    finite_force = np.flatnonzero(np.isfinite(force_scaled_kip))
    if finite_force.size:
        local_i = int(np.nanargmax(np.abs(force_scaled_kip[finite_force])))
        peak_i = int(finite_force[local_i])
        dic_peak_kip = float(abs(force_scaled_kip[peak_i]))
        if np.isfinite(disp_scaled_in[peak_i]) and np.isfinite(disp_zero_mm):
            row["dic_disp_zeroed"] = float(disp_scaled_in[peak_i] * IN2MM - disp_zero_mm)

        if np.isfinite(peak_N) and dic_peak_kip > 0:
            scale_dic2mts_peaks = peak_N / dic_peak_kip
            row["dic_force_scaled2mts"] = dic_peak_kip * scale_dic2mts_peaks
            row["scale_dic2mts_peaks"] = scale_dic2mts_peaks

    return row

def write_dic_report(coupons) -> Path:
    rows = [dic_report_row(cid, coupon_dir(cid)) for cid in coupons]
    out_fp = DIC_DIR / "raw_dic_force_displacement_signal_report.csv"
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_fp, index=False, float_format="%.8g")
    print(f"DIC report: {out_fp}")
    return out_fp


# =============================================================================
# MAIN
# =============================================================================
def process_coupon(cid, spec) -> dict:
    print(f"[{cid}]")
    cdir = coupon_dir(cid)
    if not cdir.is_dir():
        print(f"  [skip] directory not found: {cdir}")
        return {"coupon": cid}

    stats = {"coupon": cid}
    if DO_EXPORT_FRAMES:
        stats.update(export_frames(cid, cdir))
    if DO_BUILD_L1:
        build_l1(cid, cdir, spec)
    if DO_SIGNAL_PLOTS:
        plot_mts_signals(cid)
        plot_dic_sync_signals(cid, cdir)
    return stats


def main():
    t0 = time.time()
    print("=" * 70)
    print("DIC_Level1 — export .out frames + pair with MTS, virtual extensometers")
    print("=" * 70)
    for exp, root in DATA_ROOTS.items():
        print(f"Data root ({exp}): {root}")
    print(f"DIC dir   : {DIC_DIR}")
    print(f"Figs root : {FIGS_ROOT}")
    print()

    spec = load_specimen_sheet()
    coupons = selected_coupons()
    print(f"Processing {len(coupons)} coupon(s)\n")

    summary = []
    for cid in coupons:
        try:
            summary.append(process_coupon(cid, spec))
        except Exception as ex:
            print(f"[{cid}] [error] {ex}")
        print()

    if DO_SIGNAL_PLOTS:
        write_dic_report(coupons)

    print(f"Done. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
