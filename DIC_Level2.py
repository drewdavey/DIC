#!/usr/bin/env python3
"""
DIC_Level2.py  —  FSR Tensile Coupons
======================================
Pairs per-frame DIC data (Level-1 CSVs) with force/displacement,
inserts virtual axial + transverse extensometers, and writes one compact
result CSV per coupon (full untruncated record). Also saves a raw
MTS force-displacement plot (relative displacement) to figs.

No truncation is applied here — the full DIC record is preserved so that
Level-3 can apply (and tune) failure truncation independently.

INPUT per coupon
  <coupon_dir>/<coupon_id>.csv          VIC sync CSV (analog channels @ DIC frame rate)
  <coupon_dir>/<coupon_id>-*.csv        per-frame full-field DIC CSVs (Level 1 output)
  <MTS_DIR>/<coupon_id>*.txt            MTS raw file: cols disp_mm, force_N, output_V, time_s
  FSR-SpecimenTesting.xlsx              gauge thickness × width  →  area

OUTPUT per coupon
  <DIC_DIR>/<coupon_id>_L2.csv          step, time_s, force_N, disp_mm,
                                        stress_MPa, strain_axial, strain_transverse
  <FIGS_ROOT>/<coupon_id>/MTS_force_disp.png    raw MTS curve (sanity check)

UNIT NOTES
  MTS .txt is already in mm / N / V / sec (verified against tensile_analysis.py).
  VIC sync CSV "Load" column units are device-dependent; rather than guess,
  we use a combined scale factor (SCALE_N_PER_UNIT) derived from all coupons.
"""

from __future__ import annotations
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
DIC_DIR = MTS_DIR.parent / "DIC"   # consolidated _L2.csv files land here
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
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "SW": True, "UV": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]
OVERWRITE  = True

# =============================================================================
# VIRTUAL EXTENSOMETER  — ASTM D638 §5.2.1 (Class B-2 equivalent for modulus)
# Gauge length 50 mm (2 in) per D638 Type I Fig. 1.
# =============================================================================
AXIAL_GAUGE_IN = 4.36       # axial (Y, loading) gauge length, inches  [D638 G = 2.00 in]
TRANS_GAUGE_IN = 1.0       # transverse (X) gauge length, inches      [Annex A3.5.2]

# =============================================================================
# CONSTANTS
# =============================================================================
IN2MM    = 25.4
HEADERS  = 8               # MTS .txt has an 8-line header (verified)

# Scaling is intentionally NOT applied in Level-2. The raw sync voltage
# (load_raw) is saved to the L2 CSV alongside mts_peak_N and area_mm2 so
# that Level-3 can: smooth → scale → truncate, without re-running Level-2.
# The calibration pass below reports per-coupon scale factors for reference;
# set SCALE_N_PER_UNIT in Level-3 once you have a stable value.

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
    exp = cid.split("-")[1][1:-2]
    return DATA_ROOTS[exp] / cid

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

def compute_scale_factor(cid) -> float | None:
    """Return N-per-sync-unit scale factor for cid, or None on any failure."""
    cdir = coupon_dir(cid)
    sync_fp = find_sync_csv(cdir, cid)
    mts_fp  = find_mts_txt(cid)
    if sync_fp is None or mts_fp is None:
        return None
    try:
        mts = load_mts_txt(mts_fp)
        peak_N = float(mts["force_N"].abs().max())
        sync = pd.read_csv(sync_fp)
        load_col = pick_col(sync, "load")
        if load_col is None:
            return None
        load_raw = pd.to_numeric(sync[load_col], errors="coerce").to_numpy()
        raw_peak = float(np.nanmax(np.abs(load_raw)))
        if raw_peak <= 0 or not np.isfinite(raw_peak):
            return None
        return peak_N / raw_peak
    except Exception:
        return None


def get_area_mm2(spec, cid):
    """ASTM D638 §11.2: stress uses *original* cross-sectional area."""
    if cid not in spec.index:
        return None
    row = spec.loc[cid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    return float(row["t_in"]) * float(row["w_in"]) * IN2MM * IN2MM

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
# CORE
# =============================================================================
def process_coupon(cid, spec):
    print(f"[{cid}]")
    cdir = coupon_dir(cid)
    if not cdir.is_dir():
        print(f"  [skip] {cdir} not found"); return

    DIC_DIR.mkdir(parents=True, exist_ok=True)
    out_fp = DIC_DIR / f"{cid}_L2.csv"
    if out_fp.exists() and not OVERWRITE:
        print(f"  [skip] _L2.csv exists"); return

    sync_fp   = find_sync_csv(cdir, cid)
    frame_fps = find_frame_csvs(cdir, cid)
    mts_fp    = find_mts_txt(cid)
    area_mm2  = get_area_mm2(spec, cid)

    if sync_fp is None:    print(f"  [skip] no sync CSV {cid}.csv");           return
    if not frame_fps:      print(f"  [skip] no Level-1 per-frame CSVs");       return
    if mts_fp is None:     print(f"  [skip] no MTS .txt — needed for mts_peak_N"); return
    if area_mm2 is None:   print(f"  [warn] no area for {cid} — stress NaN")

    print(f"  {len(frame_fps)} frames | area = "
          + (f"{area_mm2:.2f} mm²" if area_mm2 else "N/A"))

    # ---- raw MTS load ----
    mts = load_mts_txt(mts_fp)
    plot_mts(cid, mts)
    peak_force_N = float(mts["force_N"].abs().max())
    print(f"  MTS peak force: {peak_force_N/1000:.2f} kN")

    # ---- sync CSV: read raw load signal (unscaled) ----
    # Scaling is deferred to Level-3: smooth(load_raw) → scale → force_N.
    sync = pd.read_csv(sync_fp)
    load_col = pick_col(sync, "load")
    disp_col = pick_col(sync, "drift")
    time_col = pick_col(sync, "time")
    if load_col is None:
        print(f"  [skip] no LOAD column. Cols: {list(sync.columns)}"); return

    load_raw = pd.to_numeric(sync[load_col], errors="coerce").to_numpy()
    raw_peak = float(np.nanmax(np.abs(load_raw)))
    if raw_peak <= 0 or not np.isfinite(raw_peak):
        print(f"  [skip] sync load peak invalid"); return
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
    # force_N and stress_MPa are NOT saved — Level-3 computes them after
    # smoothing and per-coupon scaling. load_raw, mts_peak_N, area_mm2 provide
    # everything Level-3 needs without re-running this slow script.
    # strain stays "raw" — toe compensation per D638 Annex A1 done in Level 3.
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


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("DIC_Level2 — pair MTS + DIC, virtual extensometers  (no truncation)")
    print("=" * 70)
    spec = load_specimen_sheet()
    coupons = selected_coupons()

    # ---- scale factor calibration pass (informational — scaling done in Level-3) ----
    sf_map = {cid: compute_scale_factor(cid) for cid in coupons}
    valid  = np.array([v for v in sf_map.values() if v is not None])
    if len(valid):
        mean_sf = float(np.mean(valid))
        std_sf  = float(np.std(valid, ddof=1) if len(valid) > 1 else 0.0)
        print(f"Scale factors across {len(valid)} coupon(s) [set in Level-3]:")
        print(f"  mean = {mean_sf:.4f} N/unit,  std = {std_sf:.4f}  "
              f"({100*std_sf/mean_sf:.2f}% CV)")
        for cid in coupons:
            v = sf_map[cid]
            print(f"    {cid}: {v:.4f}" if v is not None else f"    {cid}: N/A")
        print(f"\n  → Set SCALE_N_PER_UNIT = {mean_sf:.4f} in Level-3\n")

    print(f"Processing {len(coupons)} coupon(s)\n")
    for cid in coupons:
        try:
            process_coupon(cid, spec)
        except Exception as ex:
            print(f"[{cid}] [error] {ex}")
        print()
    print(f"Done. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
