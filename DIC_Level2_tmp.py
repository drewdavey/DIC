#!/usr/bin/env python3
"""
DIC_Level2_tmp.py
=================
Temporary, non-destructive alignment test for replacing bad DIC-sync
displacement with raw MTS displacement.

This script does not overwrite Level-2 outputs, does not update Excel, and does
not save any CSV data. By default it only prints diagnostics. Set SAVE_PLOTS =
True to write comparison figures into figs/_DIC_Level2_tmp_alignment.

Alignment idea
--------------
The raw DIC sync CSV has a trustworthy force channel and one row per DIC frame,
but its displacement channel is unreliable. The raw MTS .txt file has reliable
force and displacement, but it was started independently and does not share a
timestamp with the DIC sync CSV.

So, for each coupon:
  1. Find the peak force row in the DIC sync CSV.
  2. Find the peak force row in the raw MTS file.
  3. Estimate the DIC sync sampling rate from its Time column.
  4. Treat those two peak-force rows as the same physical instant.
  5. For each DIC frame, compute its time offset from DIC peak force.
  6. Interpolate raw MTS force/displacement at MTS_peak_time + that offset.
  7. Combine mapped MTS force with DIC Level-1 strain and area to inspect
     stress-strain values without writing production outputs.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
FIGS_ROOT = DIC_DIR.parent / "figs"
TMP_FIG_DIR = FIGS_ROOT / "_DIC_Level2_tmp_alignment"

DATA_ROOTS = {
    "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
    "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
}

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS = ["P01"]
EXPOSURES = {"CL": True, "SW": True, "UV": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

SAVE_PLOTS = True
PLOT_COUPONS = {"P01-TCL00-01"}  # used only when SAVE_PLOTS is True

# =============================================================================
# CONSTANTS / COLUMNS
# =============================================================================
HEADERS = 8
IN2MM = 25.4
KIP2N = 4448.2216152605

DIC_FORCE_SCALED_COL = "LOAD_[kip]_|_CH07_/ai2_scaled"
DIC_DISP_SCALED_COL = "DRIFT_[in]|_CH06/ai1_scaled"

LOAD_START_FRAC = 0.02
LOAD_END_FRAC = 0.50
MODULUS_STRAIN_RANGE = (0.0005, 0.003)


def coupon_id(p: str, e: str, d: str, r: str) -> str:
    return f"{p}-T{e}{d}-{r}"


def selected_coupons() -> list[str]:
    return [
        coupon_id(p, e, d, r)
        for p in PRINTS
        for e, on in EXPOSURES.items() if on
        for d, on2 in DIRECTIONS.items() if on2
        for r in REPLICATES
    ]


def parse_id(cid: str) -> tuple[str, str]:
    part = cid.split("-")[1]
    return part[1:-2], part[-2:]


def coupon_dir(cid: str) -> Path:
    exposure, _ = parse_id(cid)
    return DATA_ROOTS[exposure] / cid


def find_l1(cid: str) -> Path | None:
    fp = DIC_DIR / f"{cid}_L1.csv"
    return fp if fp.exists() else None


def find_mts_txt(cid: str) -> Path | None:
    exact = [MTS_DIR / f"{cid}.txt", MTS_DIR / f"{cid}-TEST.txt"]
    for fp in exact:
        if fp.exists():
            return fp
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    return hits[0] if hits else None


def find_dic_sync_csv(cid: str) -> Path | None:
    cdir = coupon_dir(cid)
    direct = cdir / f"{cid}.csv"
    if direct.exists():
        return direct
    if not cdir.is_dir():
        return None
    hits = sorted(cdir.rglob(f"{cid}.csv"))
    return hits[0] if hits else None


def pick_col(df: pd.DataFrame, hint: str) -> str | None:
    for col in df.columns:
        if hint.lower() in str(col).lower():
            return col
    return None


def numeric(values) -> np.ndarray:
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def relative_time(values, fallback_len: int) -> np.ndarray:
    arr = numeric(values)
    finite = np.flatnonzero(np.isfinite(arr))
    if finite.size:
        return arr - arr[finite[0]]
    return np.arange(fallback_len, dtype=float)


def sample_rate_hz(time_s: np.ndarray) -> float:
    dt = np.diff(time_s[np.isfinite(time_s)])
    dt = dt[np.isfinite(dt) & (dt > 0)]
    if not dt.size:
        return np.nan
    return float(1.0 / np.median(dt))


def load_mts_txt(fp: Path) -> pd.DataFrame:
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


def truncate_window(force_N: np.ndarray) -> tuple[int, int]:
    peak = float(np.nanmax(np.abs(force_N)))
    if not np.isfinite(peak) or peak <= 0:
        return 0, len(force_N) - 1
    i_peak = int(np.nanargmax(np.abs(force_N)))
    starts = np.where(np.abs(force_N) > LOAD_START_FRAC * peak)[0]
    i0 = int(starts[0]) if len(starts) else 0
    post = np.where(np.abs(force_N[i_peak:]) < LOAD_END_FRAC * peak)[0]
    # stop one frame BEFORE the first post-UTS drop-off (exclude the failure point)
    i1 = int(i_peak + post[0]) - 1 if len(post) else len(force_N) - 1
    return i0, i1


def modulus_gpa(strain: np.ndarray, stress_mpa: np.ndarray) -> tuple[float, float]:
    lo, hi = MODULUS_STRAIN_RANGE
    ok = (
        (strain >= lo)
        & (strain <= hi)
        & np.isfinite(strain)
        & np.isfinite(stress_mpa)
    )
    if ok.sum() < 3:
        return np.nan, np.nan
    slope, intercept = np.polyfit(strain[ok], stress_mpa[ok], 1)
    toe = -intercept / slope if slope != 0 else np.nan
    return float(slope / 1000.0), float(toe)


def analyze_coupon(cid: str) -> dict | None:
    l1_fp = find_l1(cid)
    sync_fp = find_dic_sync_csv(cid)
    mts_fp = find_mts_txt(cid)
    if l1_fp is None or sync_fp is None or mts_fp is None:
        print(f"[{cid}] skip missing file(s): L1={bool(l1_fp)} sync={bool(sync_fp)} MTS={bool(mts_fp)}")
        return None

    l1 = pd.read_csv(l1_fp)
    sync = pd.read_csv(sync_fp)
    mts = load_mts_txt(mts_fp)

    time_col = pick_col(sync, "time")
    if DIC_FORCE_SCALED_COL not in sync.columns or time_col is None:
        print(f"[{cid}] skip: sync CSV missing force or time column")
        return None

    n = min(len(l1), len(sync))
    l1 = l1.iloc[:n].reset_index(drop=True)
    sync = sync.iloc[:n].reset_index(drop=True)

    sync_time = relative_time(sync[time_col], len(sync))
    sync_fs = sample_rate_hz(sync_time)
    sync_force_kip = numeric(sync[DIC_FORCE_SCALED_COL])
    sync_force_N_nominal = sync_force_kip * KIP2N
    sync_peak_i = int(np.nanargmax(np.abs(sync_force_kip)))

    mts_time = relative_time(mts["time_s"], len(mts))
    mts_fs = sample_rate_hz(mts_time)
    mts_force_N = numeric(mts["force_N"])
    mts_disp_mm = numeric(mts["disp_mm"])
    mts_peak_i = int(np.nanargmax(np.abs(mts_force_N)))

    # Peak-force anchor map: DIC frame i maps to this MTS relative time.
    dic_dt_from_peak = sync_time - sync_time[sync_peak_i]
    mts_mapped_time = mts_time[mts_peak_i] + dic_dt_from_peak

    mapped_force_N = np.interp(mts_mapped_time, mts_time, mts_force_N, left=np.nan, right=np.nan)
    mapped_disp_mm = np.interp(mts_mapped_time, mts_time, mts_disp_mm, left=np.nan, right=np.nan)
    first_valid_disp = mapped_disp_mm[np.flatnonzero(np.isfinite(mapped_disp_mm))[0]]
    mapped_disp_zeroed = mapped_disp_mm - first_valid_disp

    area = float(l1["area_mm2"].iloc[0]) if "area_mm2" in l1.columns else np.nan
    strain = numeric(l1["strain_axial"])
    stress_mpa = mapped_force_N / area if np.isfinite(area) else np.full(n, np.nan)

    i0, i1 = truncate_window(mapped_force_N)
    sl = slice(i0, i1 + 1)
    e_gpa, toe = modulus_gpa(strain[sl], stress_mpa[sl])
    uts_i_local = int(np.nanargmax(stress_mpa[sl]))
    uts_i = i0 + uts_i_local
    uts_mpa = float(stress_mpa[uts_i])
    eps_uts_raw = float(strain[uts_i])
    eps_uts_toe = float(strain[uts_i] - toe) if np.isfinite(toe) else np.nan

    peak_force_error_pct = (
        100.0 * (mapped_force_N[sync_peak_i] - mts_force_N[mts_peak_i]) / mts_force_N[mts_peak_i]
        if mts_force_N[mts_peak_i] != 0
        else np.nan
    )

    out = {
        "coupon": cid,
        "n_frames": n,
        "sync_fs_Hz": sync_fs,
        "mts_fs_Hz": mts_fs,
        "sync_peak_i": sync_peak_i,
        "mts_peak_i": mts_peak_i,
        "sync_peak_time_s": float(sync_time[sync_peak_i]),
        "mts_peak_time_s": float(mts_time[mts_peak_i]),
        "sync_peak_nominal_N": float(abs(sync_force_N_nominal[sync_peak_i])),
        "mts_peak_N": float(abs(mts_force_N[mts_peak_i])),
        "mapped_peak_N": float(abs(mapped_force_N[sync_peak_i])),
        "peak_force_error_pct": float(peak_force_error_pct),
        "mapped_disp_at_peak_mm": float(mapped_disp_zeroed[sync_peak_i]),
        "area_mm2": area,
        "UTS_MPa": uts_mpa,
        "eps_UTS_raw_pct": eps_uts_raw * 100.0,
        "eps_UTS_toe_pct": eps_uts_toe * 100.0,
        "E_GPa": e_gpa,
        "toe_strain_pct": toe * 100.0 if np.isfinite(toe) else np.nan,
        "truncate_i0": i0,
        "truncate_i1": i1,
    }

    print(
        f"[{cid}] sync_fs={sync_fs:.3f} Hz  mts_fs={mts_fs:.2f} Hz  "
        f"peak_i sync/MTS={sync_peak_i}/{mts_peak_i}  "
        f"MTS_peak={out['mts_peak_N']:.1f} N  mapped_peak={out['mapped_peak_N']:.1f} N  "
        f"disp@peak={out['mapped_disp_at_peak_mm']:.3f} mm  "
        f"UTS={uts_mpa:.2f} MPa  eps_UTS={out['eps_UTS_toe_pct']:.3f}%  E={e_gpa:.3f} GPa"
    )

    if SAVE_PLOTS and cid in PLOT_COUPONS:
        plot_alignment(
            cid, sync_time, sync_force_N_nominal, mts_time, mts_force_N,
            mts_mapped_time, mapped_force_N, mapped_disp_zeroed, strain, stress_mpa,
            sync_peak_i, mts_peak_i, out
        )

    return out


def plot_alignment(
    cid: str,
    sync_time: np.ndarray,
    sync_force_N: np.ndarray,
    mts_time: np.ndarray,
    mts_force_N: np.ndarray,
    mts_mapped_time: np.ndarray,
    mapped_force_N: np.ndarray,
    mapped_disp_mm: np.ndarray,
    strain: np.ndarray,
    stress_mpa: np.ndarray,
    sync_peak_i: int,
    mts_peak_i: int,
    summary: dict,
) -> None:
    TMP_FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))

    ax00 = ax[0, 0]
    ax00b = ax00.twinx()
    l00a, = ax00.plot(sync_time - sync_time[sync_peak_i], sync_force_N,
                      color="C0", label="DIC sync force nominal N")
    l00b, = ax00b.plot(mts_time - mts_time[mts_peak_i], mts_force_N,
                       color="C1", alpha=0.8, label="raw MTS force N")
    ax00.axvline(0, color="k", lw=0.8, linestyle="--")
    ax00.set_xlabel("Time from peak force (s)")
    ax00.set_ylabel("DIC sync force (N)", color="C0")
    ax00b.set_ylabel("raw MTS force (N)", color="C1")
    ax00.tick_params(axis="y", labelcolor="C0")
    ax00b.tick_params(axis="y", labelcolor="C1")
    ax00.legend(handles=[l00a, l00b], fontsize=8)
    ax00.grid(alpha=0.25)

    ax01 = ax[0, 1]
    ax01b = ax01.twinx()
    l01a, = ax01.plot(sync_time, mapped_force_N,
                      color="C0", label="MTS force mapped to DIC frames")
    l01b, = ax01b.plot(sync_time, mapped_disp_mm,
                       color="C1", label="MTS disp mapped/zeroed (mm)")
    ax01.set_xlabel("DIC sync time from first frame (s)")
    ax01.set_ylabel("Force (N)", color="C0")
    ax01b.set_ylabel("Displacement (mm)", color="C1")
    ax01.tick_params(axis="y", labelcolor="C0")
    ax01b.tick_params(axis="y", labelcolor="C1")
    ax01.legend(handles=[l01a, l01b], fontsize=8)
    ax01.grid(alpha=0.25)

    ax[1, 0].plot(strain * 100, stress_mpa, lw=1)
    ax[1, 0].set_xlabel("DIC axial strain (%)")
    ax[1, 0].set_ylabel("Mapped MTS stress (MPa)")
    ax[1, 0].grid(alpha=0.25)

    ax[1, 1].axis("off")
    txt = "\n".join(f"{k}: {v:.6g}" if isinstance(v, (int, float, np.floating)) else f"{k}: {v}"
                    for k, v in summary.items())
    ax[1, 1].text(0, 1, txt, va="top", family="monospace", fontsize=8)

    fig.suptitle(f"{cid} peak-force anchored MTS -> DIC-frame mapping")
    fig.tight_layout()
    fig.savefig(TMP_FIG_DIR / f"{cid}_alignment_tmp.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print("=" * 80)
    print("DIC_Level2_tmp: peak-force alignment test only; no production outputs written")
    print("=" * 80)
    rows = [r for cid in selected_coupons() if (r := analyze_coupon(cid)) is not None]

    if rows:
        df = pd.DataFrame(rows)
        print("\nSummary ranges:")
        for col in ["sync_fs_Hz", "mts_fs_Hz", "peak_force_error_pct", "mapped_disp_at_peak_mm", "UTS_MPa", "eps_UTS_toe_pct", "E_GPa"]:
            vals = pd.to_numeric(df[col], errors="coerce")
            print(f"  {col}: min={vals.min():.6g}  median={vals.median():.6g}  max={vals.max():.6g}")
    if SAVE_PLOTS:
        print(f"\nPlots written to {TMP_FIG_DIR}")
    else:
        print("\nSAVE_PLOTS=False, so no plots or data files were written.")


if __name__ == "__main__":
    main()
