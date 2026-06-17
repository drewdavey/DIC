#!/usr/bin/env python3
"""
DIC_Level1.py  —  FSR Tensile Coupons
======================================
Step A — exports each coupon's raw VIC-3D .out files to per-frame CSVs
(next to the .out files, on the raw data drive).
Step B — pairs those per-frame CSVs with the MTS force/displacement record,
inserts virtual axial + transverse extensometers, and writes one compact
per-coupon CSV (full, untruncated record) to DIC_DIR.

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
DATA_ROOTS = {
    "CL": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TCL"),
    "SW": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
    "UV": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
    "IS": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
}
MTS_DIR = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\MTS"
)
DIC_DIR = MTS_DIR.parent / "DIC"   # consolidated _L1.csv files land here
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

# Pipeline steps — each independently skippable/overwritable so the slow
# .out export never has to re-run just to rebuild the consolidated CSV.
DO_EXPORT_FRAMES = True     # Step A: export each .out to a CSV next to it
OVERWRITE_FRAMES = False    # if False, skip .out files whose .csv already exists
DO_BUILD_L1      = True     # Step B: pair frames + MTS, build extensometers
OVERWRITE_L1     = False    # if False, skip coupons whose _L1.csv already exists

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

    print(f"Done. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
