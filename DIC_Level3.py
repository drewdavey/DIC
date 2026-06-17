#!/usr/bin/env python3
"""
DIC_Level3.py  —  FSR Tensile Coupons
======================================
Plot-only. Reads Level-2's per-frame curve CSVs and the scalar mechanical
properties Level-2 already wrote into FSR-SpecimenTesting.xlsx (the single
source of truth for scalars — nothing is recomputed here), and saves
per-coupon σ-ε / Poisson plots. No CSV output.

Group overlay plots (tensile_curves_DIC.png, tensile_summary_DIC.png) are
produced separately by tensile_group_plots.py (like matlab/tensile_plots.m),
which also reads scalars from FSR-SpecimenTesting.xlsx.

INPUT per coupon
  <DIC_DIR>/<coupon_id>_L2.csv   step, eps, sig, eps_t, i_uts, eps_raw, sig_raw
                                 (written by DIC_Level2.py)
  FSR-SpecimenTesting.xlsx       scalar properties (E, toe strain, yield
                                 stress/strain, UTS, strain at UTS, Poisson's
                                 ratio), one row per coupon, matched by
                                 "Specimen ID" — written by DIC_Level2.py

OUTPUT
  {FIGS_ROOT}/{coupon_id}/stress_strain_DIC.png   per-coupon σ–ε (toe-corrected, to UTS)
  {FIGS_ROOT}/{coupon_id}/poisson_DIC.png         per-coupon −ε_xx vs ε_yy (to UTS)
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

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# PATHS
# =============================================================================
FIGS_ROOT = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\figs"
)
DIC_DIR = FIGS_ROOT.parent / "DIC"   # per-frame _L2.csv files written by Level 2
SPECIMEN_SHEET = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.xlsx"
)

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "UV": True, "SW": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

# =============================================================================
# PLOT-ANNOTATION RANGES — keep in sync with DIC_Level2.py
# (used only to draw the modulus/offset/Poisson reference lines; the actual
# scalar values come from the Excel sheet, not recomputed here)
# =============================================================================
MODULUS_STRAIN_RANGE = (0.0005, 0.003)
YIELD_OFFSET = 0.002
POISSON_RANGE = (0.0005, 0.0025)
POISSON_CHORD_AT = 0.002

# Excel column headers for each scalar property — keep in sync with
# DIC_Level2.py's SPECIMEN_SHEET_COLUMNS (this is the inverse mapping).
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

# Airtech reference values (printed material spec — comparison lines)
AIRTECH_UTS = {0: 79.3, 45: None, 90: 25.9}    # MPa
AIRTECH_E   = {0: 6.6,  45: None, 90: 3.7}     # GPa

# =============================================================================
# DISPLAY
# =============================================================================
EXPOSURE_COLORS = {"CL": "#1f77b4", "SW": "#17becf", "UV": "#ff7f0e", "IS": "#2ca02c"}
EXPOSURE_LABELS = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
DIRECTION_MARKERS = {"00": "o", "45": "s", "90": "^"}

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

def find_l2(cid):
    p = DIC_DIR / f"{cid}_L2.csv"
    return p if p.exists() else None

def load_specimen_scalars() -> pd.DataFrame:
    """Read scalar properties back out of SPECIMEN_SHEET, indexed by coupon,
    with columns renamed from Excel headers back to property keys."""
    df = pd.read_excel(SPECIMEN_SHEET).set_index("Specimen ID")
    inverse = {label: key for key, label in SPECIMEN_SHEET_COLUMNS.items()}
    cols = {c: inverse[c] for c in df.columns if c in inverse}
    return df.rename(columns=cols)[list(cols.values())]

def load_props(cid, scalars: pd.DataFrame):
    """Combine the per-frame curve CSV with this coupon's scalar row into
    the same props dict shape DIC_Level2.py's compute_properties() used to
    hand to the plotter."""
    l2 = find_l2(cid)
    if l2 is None or cid not in scalars.index:
        return None
    curve = pd.read_csv(l2)
    row = scalars.loc[cid]
    props = {k: float(row[k]) if pd.notna(row[k]) else np.nan
              for k in SPECIMEN_SHEET_COLUMNS}
    props["_eps"]     = curve["eps"].to_numpy()
    props["_sig"]     = curve["sig"].to_numpy()
    props["_eps_t"]   = curve["eps_t"].to_numpy()
    props["_i_uts"]   = int(curve["i_uts"].iloc[0])
    props["_eps_raw"] = curve["eps_raw"].to_numpy() if "eps_raw" in curve.columns else None
    props["_sig_raw"] = curve["sig_raw"].to_numpy() if "sig_raw" in curve.columns else None
    return props


# =============================================================================
# PER-COUPON PLOTS
# =============================================================================
def plot_stress_strain(cid, props, fig_dir):
    """σ-ε curve, toe-corrected, truncated at UTS. The raw (pre-smoothing,
    but still truncated) signal is drawn behind it in light gray, with its
    own peak marked, for comparison."""
    exp, d_str = parse_id(cid)
    direction  = int(d_str)

    sl  = slice(0, props["_i_uts"] + 1)
    e_p = props["_eps"][sl] * 100   # % strain
    s_p = props["_sig"][sl]

    fig, ax = plt.subplots(figsize=(7, 4.8))

    eps_r, sig_r = props.get("_eps_raw"), props.get("_sig_raw")
    if eps_r is not None and sig_r is not None and np.any(np.isfinite(sig_r)):
        ax.plot(eps_r * 100, sig_r, lw=0.8, color="0.8", zorder=1,
                label="raw (unsmoothed)")
        i_raw = int(np.nanargmax(sig_r))
        ax.plot(eps_r[i_raw] * 100, sig_r[i_raw], "^", color="0.6", ms=8,
                zorder=2, label=f"raw UTS = {sig_r[i_raw]:.1f} MPa")

    ax.plot(e_p, s_p, lw=1.4, color=EXPOSURE_COLORS.get(exp, "#333"),
            label=cid, zorder=3)

    E_MPa = props["E_GPa"] * 1000.0
    if np.isfinite(E_MPa):
        # Elastic line through toe-corrected origin (D638 Annex A1)
        x_e = np.array([0.0, MODULUS_STRAIN_RANGE[1] * 1.5])
        ax.plot(x_e * 100, E_MPa * x_e, "k--", lw=0.8, alpha=0.7,
                label=f"E = {props['E_GPa']:.1f} GPa")
        # 0.2% offset line — start at YIELD_OFFSET on the toe-corrected axis
        x_o_end = max(YIELD_OFFSET + 0.005,
                      props["eps_y"] if np.isfinite(props["eps_y"]) else YIELD_OFFSET + 0.005)
        x_o = np.linspace(YIELD_OFFSET, x_o_end, 50)
        ax.plot(x_o * 100, E_MPa * (x_o - YIELD_OFFSET), "k:", lw=0.8, alpha=0.6,
                label="0.2% offset")
    if np.isfinite(props["sigma_y_MPa"]):
        ax.plot(props["eps_y"] * 100, props["sigma_y_MPa"], "o",
                color="orange", ms=7, zorder=5,
                label=f"σ_y = {props['sigma_y_MPa']:.1f} MPa")
    ax.plot(props["eps_at_UTS"] * 100, props["UTS_MPa"], "^",
            color="red", ms=8, zorder=5,
            label=f"UTS = {props['UTS_MPa']:.1f} MPa")

    ref_uts = AIRTECH_UTS.get(direction)
    if ref_uts is not None:
        ax.axhline(ref_uts, color="grey", linestyle=":", lw=1.0, alpha=0.7,
                   label=f"Airtech UTS = {ref_uts} MPa")

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Axial Strain (%)")
    ax.set_ylabel("Engineering Stress (MPa)")
    ax.set_title(f"{cid}  —  {EXPOSURE_LABELS.get(exp, exp)}, {direction}°")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    fig.tight_layout()
    out = fig_dir / "stress_strain_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out

def plot_poisson(cid, props, fig_dir):
    """−ε_xx vs ε_yy, truncated at UTS."""
    eps   = props["_eps"]
    eps_t = props["_eps_t"]
    if not np.any(np.isfinite(eps_t)):
        return None
    i_uts = props["_i_uts"]
    sl = slice(0, i_uts + 1)
    e_p, et_p = eps[sl], eps_t[sl]

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(e_p * 100, -et_p * 100, lw=1.2, label="data")
    nu_c = props["poisson_chord"]
    nu_s = props["poisson_slope"]
    if np.isfinite(nu_s):
        # show fit line over Poisson range
        x = np.linspace(POISSON_RANGE[0], POISSON_RANGE[1], 20)
        ax.plot(x * 100, nu_s * x * 100, "k--", lw=0.8, alpha=0.7,
                label=f"slope ν = {nu_s:.3f}")
    if np.isfinite(nu_c):
        ax.axvline(POISSON_CHORD_AT * 100, color="orange", ls=":", lw=0.8, alpha=0.6)
    else:
        ax.set_title(cid)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Axial strain ε_yy (%)")
    ax.set_ylabel("−Transverse strain  −ε_xx (%)")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = fig_dir / "poisson_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out

# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    FIGS_ROOT.mkdir(parents=True, exist_ok=True)
    scalars = load_specimen_scalars()
    n_plotted = 0

    for cid in selected_coupons():
        props = load_props(cid, scalars)
        if props is None:
            print(f"[{cid}] no _L2.csv or no specimen-sheet row — run Level 2 first")
            continue

        print(f"[{cid}]  E={props['E_GPa']:.2f} GPa  "
              f"σ_y={props['sigma_y_MPa']:.1f} MPa  "
              f"UTS={props['UTS_MPa']:.1f} MPa  "
              f"ε_UTS={props['eps_at_UTS']*100:.2f}%  "
              f"ν_chord={props['poisson_chord']:.3f}")

        fig_dir = FIGS_ROOT / cid
        fig_dir.mkdir(parents=True, exist_ok=True)
        plot_stress_strain(cid, props, fig_dir)
        plot_poisson(cid, props, fig_dir)
        n_plotted += 1

    print(f"\n{n_plotted} coupon(s) plotted -> {FIGS_ROOT}  "
          f"(run tensile_group_plots.py for group figures)")
    print(f"Done. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
