#!/usr/bin/env python3
"""
TensileDIC_Level1.py  —  FSR Tensile Coupons (ASTM D638)
=========================================================
Step A — export each coupon's VIC-3D .out files to per-frame CSVs next to the
         .out files.
Step B — pair those per-frame CSVs with the MTS record, insert virtual axial
         and transverse extensometers, and write one per-coupon CSV to DIC_DIR.
Step C — signal-inspection plots and a cross-check report comparing the raw MTS
         and raw DIC-sync channels directly. A sanity check on the inputs
         themselves; nothing Level 1 computes goes into it.

No truncation and no load scaling here — that is Level 2's job, so it can be
tuned without re-running Step A or Step B.

USAGE
  1. Toggle the SWITCHES section below to pick which coupons to run.
  2. pip install vicpyx
  3. python TensileDIC_Level1.py

INPUT per coupon
  <coupon_dir>/*.out             VIC-3D full-field export, one per DIC frame
  <coupon_dir>/<coupon_id>.csv   VIC sync CSV (analog channels at the frame rate)
  <MTS_DIR>/<coupon_id>*.txt     MTS raw: disp_mm, force_N, output_V, time_s
  FSR-SpecimenTesting.csv        gauge thickness x width -> area

OUTPUT
  <coupon_dir>/<out_name>.csv    Step A: one CSV per .out, EXPORT_VARS columns
  <DIC_DIR>/<coupon_id>.csv      full per-frame record: step, time_s, disp_mm,
                                 load_raw, strain_axial_raw, strain_transverse_raw.
                                 time_s is ELAPSED seconds from the first frame,
                                 not the raw epoch — see build_l1. Level 2
                                 appends to this same file.
  <DIC_DIR>/coupon_scalars.csv   one row per coupon: mts_peak_N, area_mm2,
                                 t0_epoch_s, and the gauge lengths actually used.
                                 Shared with FlexuralDIC_Level1, keyed by coupon.
  <FIGS_ROOT>/<coupon_id>/MTS_force_disp.png                     raw MTS curve
  <FIGS_ROOT>/<coupon_id>/MTS_force_displacement_signals.png     MTS vs time
  <FIGS_ROOT>/<coupon_id>/DIC_sync_force_displacement_signals.png  sync raw vs scaled
  <DIC_DIR>/raw_dic_force_displacement_signal_report.csv         Step C report

UNITS
  MTS .txt is already mm / N / V / s. The sync CSV's Load column is in
  device-dependent units; Level 2 derives the scale per coupon.
"""

from __future__ import annotations
import sys
import time
from pathlib import Path
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
DIC_DIR = MTS_DIR.parent / "DIC"
FIGS_ROOT = MTS_DIR.parent / "figs"
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

# =============================================================================
# SWITCHES — which coupons, and which steps
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "SW": True, "UV": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

DO_EXPORT_FRAMES = True     # Step A
OVERWRITE_FRAMES = False    # if False, skip .out files whose .csv already exists
DO_BUILD_L1      = True     # Step B
OVERWRITE_L1     = True     # if False, skip coupons whose per-coupon CSV exists
DO_SIGNAL_PLOTS  = True     # Step C

# Zero the Step C displacement traces to their first finite value. Does not
# affect the per-coupon CSV, whose disp_mm is zeroed separately.
ZERO_DISPLACEMENT = True

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
# VIRTUAL EXTENSOMETER  — D638 §5.2.1
#
# These are FSR reef coupons, not D638 Type I dogbones: the gauge section is
# ~34 mm wide x ~14 mm thick and the correlated ROI runs 137-149 mm along the
# loading axis, so D638 Fig. 1's 2.00 in gauge does not carry over.
# AXIAL_GAUGE_IN is sized to this specimen and sits inside the ROI on every P01
# coupon; ext_endpoints warns if that ever stops being true.
# =============================================================================
AXIAL_GAUGE_IN = 4.36       # axial (Y, loading) gauge length, in  (110.7 mm)
TRANS_GAUGE_IN = 1.0        # transverse (X) gauge length, in      [Annex A3.5.2]

# Correlation-confidence cutoff for extensometer endpoints. The P01 exports on
# disk carry no sigma column, so this is inert on that batch and takes effect
# only after a Step-A re-export that includes it. None disables it.
SIGMA_MAX = 0.10

IN2MM   = 25.4
HEADERS = 8               # MTS .txt header lines

# Sync CSV columns, Step C only (Step B picks its columns generically).
DIC_FORCE_RAW_COL    = "Dev1/ai2"
DIC_FORCE_SCALED_COL = "LOAD_[kip]_|_CH07_/ai2_scaled"
DIC_DISP_RAW_COL     = "Dev1/ai1"
DIC_DISP_SCALED_COL  = "DRIFT_[in]|_CH06/ai1_scaled"


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

def coupon_dir(cid):
    return DATA_ROOTS[cid.split("-")[1][1:-2]] / cid

def find_out_files(cdir):
    """Searched recursively — VIC-3D's exact output location varies by project."""
    return sorted(cdir.rglob("*.out"))

def find_mts_txt(cid):
    for p in (MTS_DIR / f"{cid}.txt", MTS_DIR / f"{cid}-TEST.txt"):
        if p.exists():
            return p
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    return hits[0] if hits else None

def find_sync_csv(cdir, cid):
    p = cdir / f"{cid}.csv"
    return p if p.exists() else None

def find_frame_csvs(cdir, cid):
    return sorted(cdir.glob(f"{cid}-????????_0.csv"))

def pick_col(df, hint):
    for c in df.columns:
        if hint.lower() in c.lower():
            return c
    return None

def numeric_col(df, col):
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)

def load_mts_txt(fp):
    """MTS .txt: 8-line header, tab separated, cols disp_mm force_N output_V time_s."""
    return (pd.read_csv(fp, sep="\t", skiprows=HEADERS, header=None,
                        names=["disp_mm", "force_N", "output_V", "time_s"],
                        encoding="utf-8-sig", on_bad_lines="skip")
              .apply(pd.to_numeric, errors="coerce")
              .dropna(subset=["force_N"]))

def read_specimen_table():
    """The specimen sheet. Duplicated across scripts rather than imported so
    each stays runnable on its own."""
    if not SPECIMEN_CSV.exists():
        raise RuntimeError(f"Could not read the specimen sheet: "
                           f"{SPECIMEN_CSV.name} not found")
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception as exc:              # malformed CSV, locked file
            raise RuntimeError(f"Could not read {SPECIMEN_CSV.name}: {exc}")
    raise RuntimeError(f"{SPECIMEN_CSV.name}: not decodable as "
                       f"{'/'.join(CSV_ENCODINGS)}")

def load_specimen_sheet():
    df = read_specimen_table()
    t_col = next((c for c in df.columns if "thickness" in c.lower()), None)
    w_col = next((c for c in df.columns
                  if "width" in c.lower() and "dia" in c.lower()), None)
    if t_col is None or w_col is None:
        raise RuntimeError(f"could not find thickness/width cols in {SPECIMEN_CSV.name}")
    df = df.rename(columns={t_col: "t_in", w_col: "w_in"}).set_index("Specimen ID")

    # Fail loudly instead of quietly producing NaN areas: continuing writes
    # area_mm2 = NaN for every coupon and Level 2 then reports "insufficient
    # data" with no clue why.
    w = pd.to_numeric(df["w_in"], errors="coerce")
    if not w.notna().any():
        raise RuntimeError(
            f"'{w_col}' is empty for every specimen in {SPECIMEN_CSV.name}.\n"
            f"    Fix: fill it in. That file is a plain CSV — edit it directly.")
    print(f"Specimen geometry: {SPECIMEN_CSV.name} "
          f"({int(w.notna().sum())} rows with a width)")
    return df

def get_area_mm2(spec, cid):
    """D638 §11.2: stress uses the ORIGINAL cross-sectional area. Returns None —
    not NaN — when either dimension is missing, because build_l1 tests for None."""
    if cid not in spec.index:
        return None
    row = spec.loc[cid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    t = pd.to_numeric(row["t_in"], errors="coerce")
    w = pd.to_numeric(row["w_in"], errors="coerce")
    if not (np.isfinite(t) and np.isfinite(w)):
        return None
    return float(t) * float(w) * IN2MM * IN2MM


# =============================================================================
# STEP A — .out -> per-frame CSV
# =============================================================================
def export_out_to_csv(out_path, csv_path, var_names):
    """Convert a single .out to CSV, one row per valid AOI point. sigma filters
    invalid points and is then dropped. Returns True on success."""
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
        values = ds.get_values(wanted)          # numpy structured array
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


def export_frames(cid, cdir):
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
    print(f"  Step A: {stats['n_out']} .out files, wrote {stats['n_csv_written']}, "
          f"skipped {stats['n_csv_skipped']} (already existed)")
    return stats


# =============================================================================
# POINT EXTENSOMETER  — mirrors VIC-3D's add_extensometer()
# =============================================================================
def ext_endpoints(frame_csv0, axial_mm, trans_mm):
    """The four endpoint world coordinates, from the reference frame's AOI
    centroid. Returns (Xc, Y_bot, Y_top, X_lft, X_rgt, Yc)."""
    ref = pd.read_csv(frame_csv0).dropna(subset=["X", "Y"])
    Yc, Xc = float(ref["Y"].median()), float(ref["X"].median())
    Ymin, Ymax = float(ref["Y"].min()), float(ref["Y"].max())
    Xmin, Xmax = float(ref["X"].min()), float(ref["X"].max())
    Y_top = min(Yc + axial_mm / 2, Ymax)
    Y_bot = max(Yc - axial_mm / 2, Ymin)
    X_rgt = min(Xc + trans_mm / 2, Xmax)
    X_lft = max(Xc - trans_mm / 2, Xmin)
    # A clipped gauge is whatever this coupon's correlated field happened to
    # span, not AXIAL_GAUGE_IN, and it varies coupon to coupon. Say so, and let
    # coupon_scalars.csv carry the length actually used. No P01 coupon clips.
    if Yc + axial_mm / 2 > Ymax or Yc - axial_mm / 2 < Ymin:
        print(f"    [!] axial gauge clipped to ROI: requested {axial_mm:.1f} mm, "
              f"got {Y_top - Y_bot:.1f} mm")
    if Xc + trans_mm / 2 > Xmax or Xc - trans_mm / 2 < Xmin:
        print(f"    [!] transverse gauge clipped to ROI: requested {trans_mm:.1f} mm, "
              f"got {X_rgt - X_lft:.1f} mm")
    return Xc, Y_bot, Y_top, X_lft, X_rgt, Yc


def point_extensometer(frame_csvs, x0, y0, x1, y1, sigma_max=SIGMA_MAX):
    """Two fixed endpoint markers; engineering strain from the change in
    distance between them.  eps = (L - L0) / L0.

    Two hardening changes over a literal VIC-3D transcription, neither of which
    alters the P01 numbers (both failure modes were checked and neither fires):
      1. Points are filtered on sigma before the nearest search, so a
         low-confidence nearest point loses to the next-nearest good one.
      2. The endpoints are located ONCE in the reference frame and tracked by
         those coordinates. Re-running the search each frame can silently
         substitute a neighbouring subset, which puts a discrete step in the
         strain record that survives smoothing. n_switched counts substitutions.

    Returns (eps, diag).
    """
    L0 = float(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
    diag = {"L0_mm": L0, "n_switched": 0, "n_dropped": 0,
            "p0_xy": None, "p1_xy": None}
    if L0 == 0:
        return np.full(len(frame_csvs), np.nan), diag

    def clean(fp):
        try:
            df = pd.read_csv(fp).dropna(subset=["X", "Y", "U", "V"])
        except Exception:
            return None
        if sigma_max is not None and "sigma" in df.columns:
            df = df[df["sigma"].between(0, sigma_max)]
        return df if len(df) >= 2 else None

    def nearest(df, ax, ay):
        d2 = (df["X"] - ax) ** 2 + (df["Y"] - ay) ** 2
        i = d2.idxmin()
        return i, float(np.sqrt(d2.loc[i]))

    ref = clean(frame_csvs[0])
    if ref is None:
        return np.full(len(frame_csvs), np.nan), diag
    i0r, _ = nearest(ref, x0, y0)
    i1r, _ = nearest(ref, x1, y1)
    a0 = (float(ref.loc[i0r, "X"]), float(ref.loc[i0r, "Y"]))
    a1 = (float(ref.loc[i1r, "X"]), float(ref.loc[i1r, "Y"]))
    diag["p0_xy"], diag["p1_xy"] = a0, a1

    # Subset spacing, to tell "same point" from "neighbour".
    ys = np.sort(ref["Y"].unique())
    step = float(np.median(np.diff(ys))) if ys.size > 1 else 0.0
    tol = 0.25 * step if step > 0 else 0.0

    eps = []
    for fp in frame_csvs:
        df = clean(fp)
        if df is None:
            eps.append(np.nan)
            diag["n_dropped"] += 1
            continue
        j0, d0 = nearest(df, *a0)
        j1, d1 = nearest(df, *a1)
        if tol > 0 and (d0 > tol or d1 > tol):
            diag["n_switched"] += 1
        # Requested endpoint plus the matched point's displacement — the match
        # is at most half a subset step from the request.
        dx = (x1 + float(df.loc[j1, "U"])) - (x0 + float(df.loc[j0, "U"]))
        dy = (y1 + float(df.loc[j1, "V"])) - (y0 + float(df.loc[j0, "V"]))
        eps.append((np.sqrt(dx ** 2 + dy ** 2) - L0) / L0)

    return np.array(eps), diag


# =============================================================================
# STEP B — pair frames + MTS, build extensometers, write the per-coupon CSV
# =============================================================================
def plot_mts(cid, mts):
    """Raw MTS force-displacement, displacement zeroed to the start of test."""
    fig_dir = FIGS_ROOT / cid
    fig_dir.mkdir(parents=True, exist_ok=True)
    disp_rel = mts["disp_mm"].to_numpy() - float(mts["disp_mm"].iloc[0])
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(disp_rel, mts["force_N"] / 1000.0, lw=1.0)
    ax.set_xlabel("Crosshead Displacement (mm, relative)")
    ax.set_ylabel("Force (kN)")
    ax.set_title(f"{cid}  —  raw MTS")
    ax.grid(alpha=0.3, linestyle="--")
    fig.tight_layout()
    out = fig_dir / "MTS_force_disp.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def build_l1(cid, cdir, spec):
    """Step B for one coupon. Returns its coupon_scalars row, or None if skipped."""
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    out_fp = DIC_DIR / f"{cid}.csv"
    if out_fp.exists() and not OVERWRITE_L1:
        print(f"  Step B: [skip] {out_fp.name} exists")
        return None

    sync_fp   = find_sync_csv(cdir, cid)
    frame_fps = find_frame_csvs(cdir, cid)
    mts_fp    = find_mts_txt(cid)
    area_mm2  = get_area_mm2(spec, cid)

    if sync_fp is None:
        print(f"  Step B: [skip] no sync CSV {cid}.csv");  return None
    if not frame_fps:
        print(f"  Step B: [skip] no per-frame CSVs (run Step A first)");  return None
    if mts_fp is None:
        print(f"  Step B: [skip] no MTS .txt — needed for mts_peak_N");   return None
    if area_mm2 is None:
        print(f"  Step B: [warn] no area for {cid} — stress NaN")

    print(f"  Step B: {len(frame_fps)} frames | area = "
          + (f"{area_mm2:.2f} mm²" if area_mm2 else "N/A"))

    mts = load_mts_txt(mts_fp)
    plot_mts(cid, mts)
    peak_force_N = float(mts["force_N"].abs().max())
    print(f"  MTS peak force: {peak_force_N / 1000:.2f} kN")

    # Sync CSV: raw, unscaled. Level 2 does the scaling.
    sync = pd.read_csv(sync_fp)
    load_col = pick_col(sync, "load")
    disp_col = pick_col(sync, "drift")
    time_col = pick_col(sync, "time")
    if load_col is None:
        print(f"  Step B: [skip] no LOAD column. Cols: {list(sync.columns)}")
        return None

    load_raw = numeric_col(sync, load_col)
    raw_peak = float(np.nanmax(np.abs(load_raw)))
    if raw_peak <= 0 or not np.isfinite(raw_peak):
        print(f"  Step B: [skip] sync load peak invalid")
        return None
    print(f"  sync raw peak: {raw_peak:.4f} units  |  MTS peak: {peak_force_N:.0f} N  "
          f"(implied scale {peak_force_N / raw_peak:.4f})")

    if disp_col:
        disp_mm = numeric_col(sync, disp_col) * IN2MM
        disp_mm = disp_mm - disp_mm[0]
    else:
        disp_mm = np.full_like(load_raw, np.nan)

    # Time_0_0 is an absolute Unix epoch (~1.7756e9 s). Written verbatim through
    # float_format="%.6g" it quantises to steps of 10 000 s and every row comes
    # out with the same time_s. Store ELAPSED seconds instead — ~50 s spans, so
    # %.6g keeps sub-millisecond resolution — and put the absolute start in
    # coupon_scalars.csv as t0_epoch_s so nothing is lost.
    time_abs = (numeric_col(sync, time_col) if time_col
                else np.arange(len(load_raw), dtype=float))
    finite_t = np.flatnonzero(np.isfinite(time_abs))
    t0_epoch = float(time_abs[finite_t[0]]) if finite_t.size else np.nan
    time_s = time_abs - t0_epoch if np.isfinite(t0_epoch) else time_abs

    n = min(len(load_raw), len(frame_fps))
    load_raw, disp_mm, time_s = load_raw[:n], disp_mm[:n], time_s[:n]
    frame_fps_used = frame_fps[:n]

    # Point extensometers: E0 axial, E1 transverse (D638 §5.2 / Annex A3).
    axial_mm = AXIAL_GAUGE_IN * IN2MM
    trans_mm = TRANS_GAUGE_IN * IN2MM
    Xc, Y_bot, Y_top, X_lft, X_rgt, Yc = ext_endpoints(
        frame_fps_used[0], axial_mm, trans_mm)
    print(f"    E0 axial      (Xc={Xc:.1f})  Y: {Y_bot:.1f} → {Y_top:.1f}  "
          f"L={Y_top - Y_bot:.1f} mm")
    print(f"    E1 transverse (Yc={Yc:.1f})  X: {X_lft:.1f} → {X_rgt:.1f}  "
          f"L={X_rgt - X_lft:.1f} mm")
    eps_a, diag_a = point_extensometer(frame_fps_used, Xc, Y_bot, Xc, Y_top)
    eps_t, diag_t = point_extensometer(frame_fps_used, X_lft, Yc, X_rgt, Yc)
    for nm, dg in (("axial", diag_a), ("transverse", diag_t)):
        if dg["n_switched"]:
            print(f"    [!] {nm}: endpoint substituted on {dg['n_switched']} "
                  f"frames — ROI is losing correlation at the gauge ends")
        if dg["n_dropped"]:
            print(f"    [!] {nm}: {dg['n_dropped']} frames had no usable field")

    # force_N and stress_MPa are NOT saved — Level 2 computes them after scaling
    # and appends them here. mts_peak_N and area_mm2 are per-coupon scalars, so
    # they go in coupon_scalars.csv rather than down every row. Strain stays
    # "raw": toe compensation is Level 2's job.
    pd.DataFrame({
        "step":                  np.arange(n),
        "time_s":                time_s,
        "disp_mm":               disp_mm,
        "load_raw":              load_raw,
        "strain_axial_raw":      eps_a,
        "strain_transverse_raw": eps_t,
    }).to_csv(out_fp, index=False, float_format="%.6g")
    print(f"  → DIC/{out_fp.name} ({n} rows)")

    return {"coupon": cid, "mts_peak_N": peak_force_N,
            "area_mm2": area_mm2 if area_mm2 else np.nan,
            "t0_epoch_s": t0_epoch,
            "axial_gauge_mm": Y_top - Y_bot,
            "trans_gauge_mm": X_rgt - X_lft}


# =============================================================================
# STEP C — raw signal-inspection plots + cross-check report
# Reads only the raw MTS .txt and raw sync CSV, so it checks the inputs
# themselves rather than anything Level 1 computes.
# =============================================================================
def first_finite_zeroed(values):
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not ZERO_DISPLACEMENT:
        return arr
    finite = np.flatnonzero(np.isfinite(arr))
    return arr - arr[finite[0]] if finite.size else arr

def relative_time(values, fallback_len):
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = np.flatnonzero(np.isfinite(arr))
    if finite.size:
        return arr - arr[finite[0]]
    return np.arange(fallback_len, dtype=float)

def linear_scale(raw, scaled):
    """Slope and intercept for scaled ~= slope * raw + intercept."""
    raw = np.asarray(raw, dtype=float)
    scaled = np.asarray(scaled, dtype=float)
    mask = np.isfinite(raw) & np.isfinite(scaled)
    if mask.sum() < 2:
        return np.nan, np.nan
    slope, intercept = np.polyfit(raw[mask], scaled[mask], 1)
    return float(slope), float(intercept)

def dominant_frequency_hz(values, time_s):
    """Dominant nonzero FFT frequency of a raw signal, in Hz."""
    values = np.asarray(values, dtype=float)
    time_s = np.asarray(time_s, dtype=float)
    mask = np.isfinite(values) & np.isfinite(time_s)
    if mask.sum() < 4:
        return np.nan
    y, t = values[mask], time_s[mask]
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
    return float(freq[int(np.nanargmax(amp[1:]) + 1)])

def mts_peak_disp_zeroed(cid):
    """(peak force N, zeroed displacement at peak mm, the zero offset mm)."""
    fp = find_mts_txt(cid)
    if fp is None:
        return np.nan, np.nan, np.nan
    df = load_mts_txt(fp)
    if df.empty:
        return np.nan, np.nan, np.nan
    force_N = numeric_col(df, "force_N")
    disp_mm = numeric_col(df, "disp_mm")
    finite_force = np.flatnonzero(np.isfinite(force_N))
    finite_disp = np.flatnonzero(np.isfinite(disp_mm))
    if not finite_force.size or not finite_disp.size:
        return np.nan, np.nan, np.nan
    peak_i = int(finite_force[int(np.nanargmax(np.abs(force_N[finite_force])))])
    zero_factor = float(disp_mm[finite_disp[0]])
    peak_disp = (float(disp_mm[peak_i] - zero_factor)
                 if np.isfinite(disp_mm[peak_i]) else np.nan)
    return float(abs(force_N[peak_i])), peak_disp, zero_factor


def save_signal_plot(cid, source_label, time_s, force, force_label, disp_mm, out_fp):
    fig, (ax_force, ax_disp) = plt.subplots(2, 1, figsize=(8, 5.8), sharex=True)
    ax_force.plot(time_s, force, lw=0.9, color="#2457a6")
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
                            disp_raw_v, disp_scaled_in, out_fp):
    """Raw volts and the sync CSV's own scaled channel on twin axes, so a
    mis-scaled channel is visible as a shape difference rather than an offset."""
    fig, (ax_f, ax_d) = plt.subplots(2, 1, figsize=(8.8, 6.2), sharex=True)
    panels = [
        (ax_f, force_raw_v, force_scaled_kip, DIC_FORCE_RAW_COL,
         DIC_FORCE_SCALED_COL, "Dev1/ai2 (V)", "Load (kip)", None),
        (ax_d, disp_raw_v, disp_scaled_in, DIC_DISP_RAW_COL,
         DIC_DISP_SCALED_COL, "Dev1/ai1 (V)", "Drift (in)",
         "Time from first sample (s)"),
    ]
    for ax, raw, scaled, raw_lbl, scaled_lbl, y_raw, y_scaled, xlabel in panels:
        ax2 = ax.twinx()
        l1, = ax.plot(time_s, raw, lw=0.9, color="#2457a6", label=raw_lbl)
        l2, = ax2.plot(time_s, scaled, lw=1.0, linestyle="--", color="#b3412c",
                       alpha=0.8, label=scaled_lbl)
        ax.set_ylabel(y_raw, color="#2457a6")
        ax2.set_ylabel(y_scaled, color="#b3412c")
        ax.tick_params(axis="y", labelcolor="#2457a6")
        ax2.tick_params(axis="y", labelcolor="#b3412c")
        ax.grid(alpha=0.25, linestyle="--")
        ax.legend([l1, l2], [raw_lbl, scaled_lbl], fontsize=8, loc="best")
        if xlabel:
            ax.set_xlabel(xlabel)
    fig.suptitle(f"{cid} - DIC raw sync CSV")
    fig.tight_layout()
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_mts_signals(cid):
    fp = find_mts_txt(cid)
    if fp is None:
        print(f"[{cid}] MTS: skip, no .txt found")
        return None
    df = load_mts_txt(fp)
    out_fp = FIGS_ROOT / cid / "MTS_force_displacement_signals.png"
    save_signal_plot(cid, "MTS raw file", relative_time(df["time_s"], len(df)),
                     numeric_col(df, "force_N"), "Force (N)",
                     first_finite_zeroed(df["disp_mm"]), out_fp)
    print(f"[{cid}] MTS: {out_fp}")
    return out_fp


def sync_signal_columns(cdir, cid):
    """(dataframe, time array) for a coupon's sync CSV, or (None, None) if it is
    missing or lacks the four Step-C channels."""
    fp = find_sync_csv(cdir, cid)
    if fp is None:
        return None, None
    df = pd.read_csv(fp)
    required = [DIC_FORCE_RAW_COL, DIC_FORCE_SCALED_COL,
                DIC_DISP_RAW_COL, DIC_DISP_SCALED_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"  [skip] {fp.name} missing columns: {missing}")
        return None, None
    time_col = pick_col(df, "time")
    t = (relative_time(df[time_col], len(df)) if time_col
         else np.arange(len(df), dtype=float))
    return df, t


def plot_dic_sync_signals(cid, cdir):
    df, time_s = sync_signal_columns(cdir, cid)
    if df is None:
        return None
    out_fp = FIGS_ROOT / cid / "DIC_sync_force_displacement_signals.png"
    save_dic_twin_axis_plot(cid, time_s,
                            numeric_col(df, DIC_FORCE_RAW_COL),
                            numeric_col(df, DIC_FORCE_SCALED_COL),
                            numeric_col(df, DIC_DISP_RAW_COL),
                            numeric_col(df, DIC_DISP_SCALED_COL), out_fp)
    print(f"[{cid}] DIC sync: {out_fp}")
    return out_fp


def dic_report_row(cid, cdir):
    """One row of the Step C report. Columns:
      dic_disp_raw2scaled / dic_force_raw2scaled  volts -> the sync CSV's own
                                                  scaled channel (in, kip)
      dic_disp_scaled2mm                          in -> mm, a constant
      dic_disp_zeroed / mts_disp_zeroed           displacement at peak force, mm
      dic_disp_zeroFactor / mts_disp_zeroFactor   what was subtracted to zero it
      dic_force_scaled2mts                        DIC peak scaled onto the MTS peak, N
      mts_raw_peak                                MTS peak force, N
      scale_dic2mts_peaks                         MTS peak (N) / DIC peak (kip)
      *_peakFreq                                  dominant FFT frequency, Hz
    """
    row = {"coupon": cid, "dic_disp_raw2scaled": np.nan,
           "dic_disp_scaled2mm": IN2MM, "dic_disp_zeroed": np.nan,
           "mts_disp_zeroed": np.nan, "mts_disp_zeroFactor": np.nan,
           "dic_disp_zeroFactor": np.nan, "dic_force_raw2scaled": np.nan,
           "dic_force_scaled2mts": np.nan, "mts_raw_peak": np.nan,
           "scale_dic2mts_peaks": np.nan, "dic_disp_raw_peakFreq": np.nan,
           "dic_force_raw_peakFreq": np.nan}

    df, time_s = sync_signal_columns(cdir, cid)
    if df is None:
        return row

    force_raw_v      = numeric_col(df, DIC_FORCE_RAW_COL)
    force_scaled_kip = numeric_col(df, DIC_FORCE_SCALED_COL)
    disp_raw_v       = numeric_col(df, DIC_DISP_RAW_COL)
    disp_scaled_in   = numeric_col(df, DIC_DISP_SCALED_COL)

    row["dic_disp_raw2scaled"]    = linear_scale(disp_raw_v, disp_scaled_in)[0]
    row["dic_force_raw2scaled"]   = linear_scale(force_raw_v, force_scaled_kip)[0]
    row["dic_disp_raw_peakFreq"]  = dominant_frequency_hz(disp_raw_v, time_s)
    row["dic_force_raw_peakFreq"] = dominant_frequency_hz(force_raw_v, time_s)

    peak_N, mts_disp_zeroed, mts_zero = mts_peak_disp_zeroed(cid)
    row["mts_raw_peak"] = peak_N
    row["mts_disp_zeroed"] = mts_disp_zeroed
    row["mts_disp_zeroFactor"] = mts_zero

    finite_disp = np.flatnonzero(np.isfinite(disp_scaled_in))
    disp_zero_mm = np.nan
    if finite_disp.size:
        disp_zero_mm = float(disp_scaled_in[finite_disp[0]] * IN2MM)
        row["dic_disp_zeroFactor"] = disp_zero_mm

    finite_force = np.flatnonzero(np.isfinite(force_scaled_kip))
    if finite_force.size:
        peak_i = int(finite_force[int(np.nanargmax(np.abs(force_scaled_kip[finite_force])))])
        dic_peak_kip = float(abs(force_scaled_kip[peak_i]))
        if np.isfinite(disp_scaled_in[peak_i]) and np.isfinite(disp_zero_mm):
            row["dic_disp_zeroed"] = float(disp_scaled_in[peak_i] * IN2MM - disp_zero_mm)
        if np.isfinite(peak_N) and dic_peak_kip > 0:
            scale = peak_N / dic_peak_kip
            row["dic_force_scaled2mts"] = dic_peak_kip * scale
            row["scale_dic2mts_peaks"] = scale
    return row


def write_dic_report(coupons):
    rows = [dic_report_row(cid, coupon_dir(cid)) for cid in coupons]
    out_fp = DIC_DIR / "raw_dic_force_displacement_signal_report.csv"
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_fp, index=False, float_format="%.8g")
    print(f"DIC report: {out_fp}")
    return out_fp


def write_coupon_scalars(rows):
    """Upsert this run's per-coupon scalars into coupon_scalars.csv — the single
    place these are stored. Coupons skipped this run keep their existing row.
    Shared with FlexuralDIC_Level1, keyed by coupon ID.

    float_format is %.12g, not the %.6g used elsewhere: t0_epoch_s is a Unix
    epoch near 1.78e9, and 6 significant figures rounds it to the nearest
    10 000 s — the same rounding that flattened time_s. See build_l1.
    """
    out_fp = DIC_DIR / "coupon_scalars.csv"
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    merged = {}
    if out_fp.exists():
        for row in pd.read_csv(out_fp).to_dict("records"):
            merged[row["coupon"]] = row
    for row in rows:
        merged[row["coupon"]] = row
    pd.DataFrame(sorted(merged.values(), key=lambda r: r["coupon"])
                 ).to_csv(out_fp, index=False, float_format="%.12g")
    print(f"Coupon scalars: {out_fp}")
    return out_fp


# =============================================================================
# MAIN
# =============================================================================
def process_coupon(cid, spec):
    print(f"[{cid}]")
    cdir = coupon_dir(cid)
    if not cdir.is_dir():
        print(f"  [skip] directory not found: {cdir}")
        return {"coupon": cid}, None

    stats = {"coupon": cid}
    scalar_row = None
    if DO_EXPORT_FRAMES:
        stats.update(export_frames(cid, cdir))
    if DO_BUILD_L1:
        scalar_row = build_l1(cid, cdir, spec)
    if DO_SIGNAL_PLOTS:
        plot_mts_signals(cid)
        plot_dic_sync_signals(cid, cdir)
    return stats, scalar_row


def main():
    t0 = time.time()
    print("=" * 70)
    print("TensileDIC_Level1 — export .out frames + pair with MTS, virtual extensometers")
    print("=" * 70)
    for exp, root in DATA_ROOTS.items():
        print(f"Data root ({exp}): {root}")
    print(f"DIC dir   : {DIC_DIR}")
    print(f"Figs root : {FIGS_ROOT}")
    print()

    spec = load_specimen_sheet()
    coupons = selected_coupons()
    print(f"Processing {len(coupons)} coupon(s)\n")

    scalar_rows = []
    for i, cid in enumerate(coupons, start=1):
        t_c = time.time()
        try:
            _, scalar_row = process_coupon(cid, spec)
            if scalar_row is not None:
                scalar_rows.append(scalar_row)
        except Exception as ex:
            print(f"[{cid}] [error] {ex}")
        elapsed = time.time() - t0
        print(f"  [{i}/{len(coupons)}] {time.time() - t_c:.1f} s this coupon | "
              f"elapsed {elapsed:.1f} s | ETA {elapsed / i * (len(coupons) - i):.1f} s\n")

    if DO_BUILD_L1:
        write_coupon_scalars(scalar_rows)
    if DO_SIGNAL_PLOTS:
        write_dic_report(coupons)

    print(f"Done. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
