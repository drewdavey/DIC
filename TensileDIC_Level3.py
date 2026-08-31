#!/usr/bin/env python3
"""
TensileDIC_Level3.py  —  FSR Tensile Coupons (ASTM D638)
=========================================================
Plots and statistics only. Reads Level-2's per-frame CSVs and the scalar
properties Level 2 wrote into the specimen sheet; nothing tensile is
recomputed here. Pin-bearing statistics (ASTM D953) ARE computed here, from
the raw MTS files — bearing has no DIC/Level-2 step of its own.
See README.md for the method and for the P01-specific caveats.

INPUT
  <DIC_DIR>/<coupon_id>.csv    per-frame record, Level 1 + Level 2 columns
  FSR-SpecimenTesting.csv      D638 scalars and coupon thickness, one row per coupon
  <MTS_DIR>/P01-T*.txt         raw MTS tensile force/displacement
  <MTS_DIR>/P01-B*.txt         raw MTS bearing force/displacement

OUTPUT
  <FIGS_ROOT>/<coupon_id>/stress_strain_DIC.png  per-coupon stress-strain
  <FIGS_ROOT>/<coupon_id>/poisson_DIC.png        per-coupon −ε_xx vs ε_yy
  <FIGS_ROOT>/tensile_mts_FD.png                 group force vs displacement
  <FIGS_ROOT>/tensile_curves_DIC.png             group stress-strain overlay
  <FIGS_ROOT>/tensile_summary_DIC.png            group property scatter
  <FIGS_ROOT>/tensile_peak_strength_DIC.png      group UTS by exposure
  stdout                                         D638 + D953 stat tables
  P01_MechanicalStats.xlsx                       same stats, Tensile + Bearing sheets
"""

from __future__ import annotations
import re
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# PATHS
# =============================================================================
FIGS_ROOT = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\figs"
)
DIC_DIR = FIGS_ROOT.parent / "DIC"
MTS_DIR = FIGS_ROOT.parent / "MTS"
SPECIMEN_CSV = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.csv"
)
OUT_STATS_XLSX = FIGS_ROOT.parent / "P01_MechanicalStats.xlsx"
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "UV": True, "SW": True, "IS": True}
DIRECTIONS = {"00": True, "45": True, "90": True}
REPLICATES = ["01", "02", "03"]

DO_PER_COUPON_PLOTS = True
DO_GROUP_PLOTS      = True
DO_PRINT_STATS      = True

# P01-TCL45-01's DIC record stops at 90.6 % of the MTS peak force, so its UTS
# and strain-at-UTS are not the specimen's. See README.md.
DIC_EXCLUDE = {"P01-TCL45-01"}

# =============================================================================
# ANALYSIS  — must match TensileDIC_Level2
# =============================================================================
MODULUS_STRAIN_RANGE = (0.0005, 0.003)   # used only to draw the tangent line
YIELD_OFFSET         = 0.002             # D638 A2.6
POISSON_RANGE        = (0.0005, 0.0025)  # D638 A3.10.1.3
POISSON_CHORD_AT     = 0.002
MTS_HEADERS = 8

# Sheet headers written by TensileDIC_Level2, inverted back to property keys.
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

# Airtech printed-material spec, drawn as comparison lines.
AIRTECH_UTS = {0: 79.3, 45: None, 90: 25.9}    # MPa
AIRTECH_E   = {0: 6.6,  45: None, 90: 3.7}     # GPa

# =============================================================================
# BEARING  — ASTM D953-19 Procedure A
# =============================================================================
T_DIRECTIONS = [0, 45, 90]
B_DIRECTIONS = [0, 90]

HOLE_D_MM   = 0.5625 * 25.4       # 14.29 mm reamed hole
DEF_4PCT_MM = 0.04 * HOLE_D_MM    # 0.572 mm, the 4 % hole-deformation point
FIT_LO, FIT_HI = 0.10, 0.40       # toe-correction fit window, fraction of F_max

# Fixture malfunction on BCL00-01; P01-BGM00 is the replacement Control 0 deg run.
SUBST_SRC   = "P01-BGM00"
SUBST_LABEL = "P01-BCL00-01"
SKIP_BEAR   = {"P01-BCL00-01", "P01-BCL00-01-TEST"}

# =============================================================================
# DISPLAY  — same palette as mts_plots.py and FlexuralDIC_Level3
# =============================================================================
EXPOSURE_ORDER    = ["CL", "UV", "SW", "IS"]
EXPOSURE_LABELS   = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
EXPOSURE_COLORS   = {"CL": "#2a78d6", "UV": "#eb6834", "SW": "#1baf7a", "IS": "#4a3aa7"}
DIRECTION_MARKERS = {"00": "o", "45": "s", "90": "^"}
DIR_STYLES        = {0: "-", 45: (0, (6.5, 3.5)), 90: (0, (1.5, 4.0))}


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

def read_specimen_csv():
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc)
        except UnicodeDecodeError:
            continue
    raise SystemExit(f"{SPECIMEN_CSV.name}: not decodable as {'/'.join(CSV_ENCODINGS)}")

def load_specimen_scalars():
    """Specimen sheet indexed by coupon, columns renamed back to property keys."""
    df = read_specimen_csv().set_index("Specimen ID")
    inverse = {label: key for key, label in SPECIMEN_SHEET_COLUMNS.items()}
    cols = {c: inverse[c] for c in df.columns if c in inverse}
    return df.rename(columns=cols)[list(cols.values())]

def load_kept_curve(cid):
    """Rows inside Level-2's analysis window, or None."""
    fp = DIC_DIR / f"{cid}.csv"
    if not fp.exists():
        return None
    df = pd.read_csv(fp)
    if "kept" not in df.columns:
        return None
    return df.loc[df["kept"].astype(bool)].reset_index(drop=True)

def load_props(cid, scalars):
    """Per-frame curve plus this coupon's scalar row, or None."""
    curve = load_kept_curve(cid)
    if curve is None or curve.empty or cid not in scalars.index:
        return None
    row = scalars.loc[cid]
    props = {k: (float(row[k]) if pd.notna(row[k]) else np.nan)
             for k in SPECIMEN_SHEET_COLUMNS}
    props["_eps"]   = curve["strain_axial"].to_numpy()
    props["_sig"]   = curve["stress_MPa"].to_numpy()
    props["_eps_t"] = curve["strain_transverse"].to_numpy()
    props["_i_uts"] = int(np.nanargmax(props["_sig"]))
    props["_eps_raw"] = (curve["strain_axial_unsmoothed"].to_numpy()
                         if "strain_axial_unsmoothed" in curve.columns else None)
    props["_sig_raw"] = (curve["stress_MPa_unsmoothed"].to_numpy()
                         if "stress_MPa_unsmoothed" in curve.columns else None)
    return props

def read_mts_txt(fp):
    """MTS .txt: 8-line header, tab separated, cols disp_mm force_N output_V time_s."""
    raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                      names=["disp_mm", "force_N", "output_V", "time_s"],
                      encoding="utf-8-sig", on_bad_lines="skip")
    return raw.apply(pd.to_numeric, errors="coerce")


# =============================================================================
# PER-COUPON PLOTS
# =============================================================================
def plot_stress_strain(cid, props, fig_dir):
    """Stress-strain, toe-corrected and truncated at UTS. The pre-smoothing
    signal is drawn behind it in grey when Level 2 kept one."""
    exp, d_str = parse_id(cid)
    direction = int(d_str)
    sl = slice(0, props["_i_uts"] + 1)
    eps = props["_eps"][sl] * 100
    sig = props["_sig"][sl]

    fig, ax = plt.subplots(figsize=(7, 4.8))

    eps_r, sig_r = props["_eps_raw"], props["_sig_raw"]
    if eps_r is not None and sig_r is not None and np.any(np.isfinite(sig_r)):
        ax.plot(eps_r * 100, sig_r, lw=0.8, color="0.8", zorder=1,
                label="raw (unsmoothed)")
        i_raw = int(np.nanargmax(sig_r))
        ax.plot(eps_r[i_raw] * 100, sig_r[i_raw], "^", color="0.6", ms=8,
                zorder=2, label=f"raw UTS = {sig_r[i_raw]:.1f} MPa")

    ax.plot(eps, sig, lw=1.4, color=EXPOSURE_COLORS.get(exp, "#333"),
            label=cid, zorder=3)

    E_MPa = props["E_GPa"] * 1000.0
    if np.isfinite(E_MPa):
        x_e = np.array([0.0, MODULUS_STRAIN_RANGE[1] * 1.5])
        ax.plot(x_e * 100, E_MPa * x_e, "k--", lw=0.8, alpha=0.7,
                label=f"E = {props['E_GPa']:.1f} GPa")
        end = props["eps_y"] if np.isfinite(props["eps_y"]) else YIELD_OFFSET + 0.005
        x_o = np.linspace(YIELD_OFFSET, max(YIELD_OFFSET + 0.005, end), 50)
        ax.plot(x_o * 100, E_MPa * (x_o - YIELD_OFFSET), "k:", lw=0.8, alpha=0.6,
                label="0.2% offset")
    if np.isfinite(props["sigma_y_MPa"]):
        ax.plot(props["eps_y"] * 100, props["sigma_y_MPa"], "o", color="orange",
                ms=7, zorder=5, label=f"σ_y = {props['sigma_y_MPa']:.1f} MPa")
    ax.plot(props["eps_at_UTS"] * 100, props["UTS_MPa"], "^", color="red",
            ms=8, zorder=5, label=f"UTS = {props['UTS_MPa']:.1f} MPa")

    ref = AIRTECH_UTS.get(direction)
    if ref is not None:
        ax.axhline(ref, color="grey", linestyle=":", lw=1.0, alpha=0.7,
                   label=f"Airtech UTS = {ref} MPa")

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Axial Strain (%)")
    ax.set_ylabel("Engineering Stress (MPa)")
    ax.set_title(f"{cid}  —  {EXPOSURE_LABELS.get(exp, exp)}, {direction} deg")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    fig.tight_layout()
    out = fig_dir / "stress_strain_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_poisson(cid, props, fig_dir):
    """−ε_xx vs ε_yy, truncated at UTS."""
    eps_t = props["_eps_t"]
    if not np.any(np.isfinite(eps_t)):
        return None
    sl = slice(0, props["_i_uts"] + 1)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(props["_eps"][sl] * 100, -eps_t[sl] * 100, lw=1.2, label="data")
    nu_s = props["poisson_slope"]
    if np.isfinite(nu_s):
        x = np.linspace(POISSON_RANGE[0], POISSON_RANGE[1], 20)
        ax.plot(x * 100, nu_s * x * 100, "k--", lw=0.8, alpha=0.7,
                label=f"slope ν = {nu_s:.3f}")
    if np.isfinite(props["poisson_chord"]):
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
# GROUP PLOTS — raw MTS force/displacement
# =============================================================================
def group_load_mts_coupons():
    coupons = []
    for fp in sorted(MTS_DIR.glob("P01-T*.txt")):
        stem = re.sub(r"-TEST$", "", fp.stem, flags=re.IGNORECASE)
        parts = stem.split("-")
        if len(parts) < 2 or len(parts[1]) < 5:
            continue
        key = parts[1].upper()
        exp_code = key[1:-2]
        try:
            d_int = int(key[-2:])
        except ValueError:
            continue
        if exp_code not in EXPOSURE_COLORS or d_int not in T_DIRECTIONS:
            continue

        raw = read_mts_txt(fp).dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue
        d = raw["disp_mm"].to_numpy() - raw["disp_mm"].iloc[0]
        f = raw["force_N"].to_numpy() / 1000
        i_peak = int(np.argmax(f))
        coupons.append({"exp": exp_code, "dir": d_int, "d": d, "f": f, "i_peak": i_peak})
        print(f"[{stem}]  peak {f[i_peak]:.2f} kN")

    print(f"\n{len(coupons)} MTS coupons loaded")
    return coupons


def group_plot_mts_force_displacement(coupons):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        ax.plot(c["d"], c["f"], color=EXPOSURE_COLORS[c["exp"]],
                ls=DIR_STYLES[c["dir"]], lw=0.9, alpha=0.75)
        ax.scatter(c["d"][c["i_peak"]], c["f"][c["i_peak"]],
                   color=EXPOSURE_COLORS[c["exp"]], s=30, zorder=5)

    ax.set_xlabel("Displacement (mm)", fontsize=12)
    ax.set_ylabel("Force (kN)", fontsize=12)
    ax.set_title("P01 — Tensile: Force vs. Displacement", fontsize=13)
    ax.grid(alpha=0.25, ls="--")
    ax.set_xlim(left=0)

    exp_active = [e for e in EXPOSURE_ORDER if any(c["exp"] == e for c in coupons)]
    handles = [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
               for e in exp_active]
    for d_int in sorted({c["dir"] for c in coupons}):
        handles.append(mlines.Line2D([], [], color="k", ls=DIR_STYLES[d_int],
                                     lw=1.2, label=f"{d_int}°"))
    handles.append(plt.scatter([], [], color="k", s=30, label="Peak force"))
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)

    fig.tight_layout()
    out = FIGS_ROOT / "tensile_mts_FD.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# =============================================================================
# GROUP PLOTS — DIC-derived
# =============================================================================
def group_load_dic_coupons(scalars):
    coupons = []
    for cid in selected_coupons():
        if cid in DIC_EXCLUDE:
            continue
        df = load_kept_curve(cid)
        if df is None or df.empty or cid not in scalars.index:
            continue
        exp, d_str = parse_id(cid)
        sig = df["stress_MPa"].to_numpy()
        sl = slice(0, int(np.nanargmax(sig)) + 1)
        row = scalars.loc[cid]
        coupons.append({
            "cid": cid, "exp": exp, "d_str": d_str, "d_int": int(d_str),
            "eps_plot": df["strain_axial"].to_numpy()[sl],
            "sig_plot": sig[sl],
            "UTS_MPa": float(row["UTS_MPa"]),
            "E_GPa": float(row["E_GPa"]),
        })
        print(f"[{cid}]  UTS={coupons[-1]['UTS_MPa']:.1f} MPa  "
              f"E={coupons[-1]['E_GPa']:.2f} GPa")
    print(f"\n{len(coupons)} DIC coupons loaded (group plots)\n")
    return coupons


def group_plot_dic_curves(coupons):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        col = EXPOSURE_COLORS[c["exp"]]
        ax.plot(c["eps_plot"] * 100, c["sig_plot"], color=col,
                ls=DIR_STYLES[c["d_int"]], lw=0.9, alpha=0.75)
        ax.scatter(c["eps_plot"][-1] * 100, c["sig_plot"][-1], color=col,
                   marker=DIRECTION_MARKERS[c["d_str"]], s=35, zorder=5)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Axial Strain (%)", fontsize=12)
    ax.set_ylabel("Stress (MPa)", fontsize=12)
    ax.set_title("P01 — Tensile: Stress vs. Strain", fontsize=13)
    ax.grid(alpha=0.25, ls="--")

    exp_active = [e for e in EXPOSURE_ORDER if EXPOSURES.get(e)]
    handles = [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
               for e in exp_active]
    for d_int in sorted({c["d_int"] for c in coupons}):
        handles.append(mlines.Line2D([], [], color="k", ls="None",
                                     marker=DIRECTION_MARKERS[f"{d_int:02d}"],
                                     markersize=6, label=f"{d_int}°"))
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)

    fig.tight_layout()
    out = FIGS_ROOT / "tensile_curves_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def group_plot_dic_property_scatter(coupons):
    dirs_present = sorted({c["d_int"] for c in coupons})
    exp_active = [e for e in EXPOSURE_ORDER if EXPOSURES.get(e)]
    panels = [("UTS_MPa", "UTS (MPa)", "Ultimate Tensile Strength", AIRTECH_UTS),
              ("E_GPa", "E (GPa)", "Young's Modulus", AIRTECH_E)]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, panel in zip(axes, panels):
        key, ylabel, title, refs = panel
        for c in coupons:
            ei = exp_active.index(c["exp"])
            x = c["d_int"] + (ei - (len(exp_active) - 1) / 2) * 1.5
            ax.scatter(x, c[key], color=EXPOSURE_COLORS[c["exp"]],
                       marker=DIRECTION_MARKERS.get(c["d_str"], "o"), s=55, zorder=5)
        for ref in refs.values():
            if ref is not None:
                ax.axhline(ref, color="grey", ls="--", lw=0.8, alpha=0.6)
        ax.set_xticks(dirs_present)
        ax.set_xticklabels([f"{d}°" for d in dirs_present])
        ax.set_xlabel("Print Direction")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25, ls="--", axis="y")

    handles = [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
               for e in exp_active]
    axes[-1].legend(handles=handles, fontsize=8, loc="best")
    fig.tight_layout()
    out = FIGS_ROOT / "tensile_summary_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def group_plot_dic_peak_strength(coupons):
    x_pos = np.arange(len(EXPOSURE_ORDER))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axvspan(x_pos[0] - 0.5, x_pos[1] + 0.5, color="#f7f5e8", alpha=0.35)
    ax.axvspan(x_pos[2] - 0.5, x_pos[3] + 0.5, color="#e8f7f5", alpha=0.35)

    for c in coupons:
        x = x_pos[EXPOSURE_ORDER.index(c["exp"])]
        ax.scatter(x, c["UTS_MPa"], color=EXPOSURE_COLORS[c["exp"]],
                   marker=DIRECTION_MARKERS.get(c["d_str"], "o"),
                   s=70, edgecolor="black", linewidth=0.4, zorder=5)

    for d_str in DIRECTIONS:
        vals = [c["UTS_MPa"] for c in coupons if c["d_str"] == d_str]
        if vals:
            ax.axhline(np.mean(vals), color="k", ls="--", lw=0.7, alpha=0.6, zorder=3)

    ax.set_xlim(-0.5, len(EXPOSURE_ORDER) - 0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([EXPOSURE_LABELS[e] for e in EXPOSURE_ORDER])
    ax.set_xlabel("Exposure Condition")
    ax.set_ylabel("UTS (MPa)")
    ax.set_title("P01 — Tensile: Max Stress")
    ax.grid(alpha=0.25, ls="--", axis="y")

    handles = [mlines.Line2D([], [], color="black", marker=DIRECTION_MARKERS[d],
                             linestyle="None", markersize=7, label=f"{d}°")
               for d in DIRECTIONS if DIRECTIONS.get(d)]
    ax.legend(handles=handles, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 1), borderaxespad=0)

    fig.tight_layout()
    out = FIGS_ROOT / "tensile_peak_strength_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# =============================================================================
# STATS — mean +/- std (CV%) tables
# =============================================================================
def stat(vals):
    """(mean, std, CV%) of the finite values; NaN for an empty group."""
    arr = np.array([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan, np.nan
    m = float(np.mean(arr))
    s = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
    cv = 100.0 * s / abs(m) if abs(m) > 1e-12 else 0.0
    return m, s, cv

def stat_cell(m, s, cv, w_m=6, w_s=5, dec_m=2, dec_s=2):
    if not np.isfinite(m):
        return f"{'NaN':>{w_m}} ± {'NaN':>{w_s}} (  N/A)"
    return f"{m:{w_m}.{dec_m}f} ± {s:{w_s}.{dec_s}f} ({cv:4.1f}%)"

def stat_groups(rows_data, exp_order, directions, exp_labels):
    """(label, direction, subset) per exposure x direction, then All exps per direction."""
    groups = []
    for exp in exp_order:
        for d in directions:
            groups.append((exp_labels.get(exp, exp), d,
                           [r for r in rows_data if r["exp"] == exp and r["dir"] == d]))
    for d in directions:
        groups.append(("All exps", d, [r for r in rows_data if r["dir"] == d]))
    return groups

def print_table(title, rows_data, exp_order, directions, props, exp_labels):
    """props : list of (key, header, dec_m, dec_s, w_m, w_s)"""
    col_w, id_w = 22, 24
    bar = "=" * (id_w + col_w * len(props))
    sep = "-" * (id_w + col_w * len(props))

    print(f"\n{title}")
    print(bar)
    print(f"{'Exposure':<12s} {'Dir':>6s}  {'n':3s}   "
          + "".join(f"{p[1]:^{col_w}}" for p in props))
    print(f"{'':12s} {'':6s}  {'':3s}   "
          + "".join(f"{'mean ± std (CV%)':^{col_w}}" for _ in props))
    print(sep)

    groups = stat_groups(rows_data, exp_order, directions, exp_labels)
    for i, (label, d, subset) in enumerate(groups):
        if i == len(groups) - len(directions):
            print(sep)
        if not subset:
            continue
        cells = []
        for key, _, dec_m, dec_s, w_m, w_s in props:
            m, s, cv = stat([r[key] for r in subset])
            cells.append(stat_cell(m, s, cv, w_m, w_s, dec_m, dec_s).ljust(col_w - 3))
        d_label = f"{d}°"
        print(f"{label:<12s} {d_label:>6s}  {len(subset):3d}   " + "   ".join(cells))
    print(bar)

def build_stats_df(rows_data, exp_order, directions, props, exp_labels):
    """mean/std/CV columns per group, plus an 'All exps' row per direction."""
    records = []
    for label, d, subset in stat_groups(rows_data, exp_order, directions, exp_labels):
        if not subset:
            continue
        rec = {"Exposure": label, "Dir (deg)": d, "n": len(subset)}
        for key, header, *_ in props:
            m, s, cv = stat([r[key] for r in subset])
            rec[header + " mean"]   = round(m, 4) if np.isfinite(m) else None
            rec[header + " std"]    = round(s, 4) if np.isfinite(s) else None
            rec[header + " CV (%)"] = round(cv, 2) if np.isfinite(cv) else None
        records.append(rec)
    return pd.DataFrame(records)

def write_stats_sheet(path, sheet, df):
    """Write df to one sheet of the workbook, leaving the other sheets alone.
    FlexuralDIC_Level3 writes a Flexural sheet into this same file."""
    if path.exists():
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="a",
                                if_sheet_exists="replace")
    else:
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="w")
    with writer:
        df.to_excel(writer, sheet_name=sheet, index=False)


# =============================================================================
# STATS — tensile rows (D638) and bearing rows (D953)
# =============================================================================
def load_spec_sheet_raw():
    """Specimen sheet with its original headers, indexed by Specimen ID —
    needed for coupon thickness (bearing area) as well as the tensile scalars."""
    df = read_specimen_csv()
    if "Print ID" in df.columns:
        df = df[df["Print ID"] == PRINTS[0]]
    return df.set_index("Specimen ID")

def load_tensile_stats_rows(spec_df):
    rows = []
    for cid in selected_coupons():
        if cid in DIC_EXCLUDE or cid not in spec_df.index:
            continue
        spec_row = spec_df.loc[cid]
        cols = {}
        for key, label in SPECIMEN_SHEET_COLUMNS.items():
            v = spec_row.get(label)
            cols[key] = float(v) if pd.notna(v) else np.nan
        # Individual properties can legitimately be NaN (a brittle 90 deg coupon
        # has no offset-yield crossing); stat() drops those per property. Only
        # E_GPa missing means the coupon has not been through Level 2 at all.
        if np.isnan(cols["E_GPa"]):
            continue
        exp, d_str = parse_id(cid)
        rows.append({
            "exp": exp, "dir": int(d_str), "cid": cid,
            "E_GPa":          cols["E_GPa"],
            "sigma_y_MPa":    cols["sigma_y_MPa"],
            "UTS_MPa":        cols["UTS_MPa"],
            "eps_at_UTS_pct": cols["eps_at_UTS"] * 100.0,
            "poisson_chord":  cols["poisson_chord"],
        })
    return rows

def load_tensile_fmax():
    """Peak load per tensile coupon, straight from the raw MTS files."""
    fmax = {}
    for fp in sorted(MTS_DIR.glob(f"{PRINTS[0]}-T*.txt")):
        stem = re.sub(r"-TEST$", "", fp.stem, flags=re.IGNORECASE)
        parts = stem.split("-")
        if len(parts) < 3 or len(parts[1]) < 5 or not parts[2].isdigit():
            continue
        key = parts[1].upper()
        exp_code = key[1:-2]
        try:
            d_int = int(key[-2:])
        except ValueError:
            continue
        if exp_code not in EXPOSURE_LABELS or d_int not in T_DIRECTIONS:
            continue
        raw = read_mts_txt(fp).dropna(subset=["force_N"])
        if len(raw) < 10:
            continue
        cid = f"{PRINTS[0]}-T{exp_code}{d_int:02d}-{int(parts[2]):02d}"
        fmax[cid] = float(np.max(raw["force_N"].to_numpy()))
    return fmax

def get_t_mm(spec_df, t_col, sid):
    if t_col is None or sid not in spec_df.index:
        return None
    v = spec_df.loc[sid, t_col]
    if isinstance(v, pd.Series):
        v = v.iloc[0]
    return float(v) * 25.4 if pd.notna(v) else None

def parse_bearing_stem(stem):
    """Return (exp_code, d_int, spec_id) or (None, None, None)."""
    canon = re.sub(r"-TEST$", "", stem, flags=re.IGNORECASE)
    if canon.upper() == SUBST_SRC:
        return "CL", 0, SUBST_LABEL
    parts = stem.split("-")
    if len(parts) < 2 or not parts[1].upper().startswith("B") or len(parts[1]) < 5:
        return None, None, None
    key = parts[1].upper()
    exp_code = key[1:-2]
    try:
        d_int = int(key[-2:])
    except ValueError:
        return None, None, None
    if exp_code not in EXPOSURE_LABELS or d_int not in B_DIRECTIONS:
        return None, None, None
    rep = ""
    if len(parts) >= 3:
        p3 = re.sub(r"^TEST", "", parts[2], flags=re.IGNORECASE)
        if p3.isdigit():
            rep = p3
    base = f"{PRINTS[0]}-B{exp_code}{d_int:02d}"
    return exp_code, d_int, (f"{base}-{int(rep):02d}" if rep else base)

def load_bearing_rows(spec_df, t_col):
    rows = []
    for fp in sorted(MTS_DIR.glob(f"{PRINTS[0]}-B*.txt")):
        stem = fp.stem
        if stem.upper() in {s.upper() for s in SKIP_BEAR}:
            continue
        exp_code, d_int, spec_id = parse_bearing_stem(stem)
        if exp_code is None:
            continue
        t_mm = get_t_mm(spec_df, t_col, spec_id)
        if t_mm is None:
            continue

        raw = read_mts_txt(fp).dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue
        d = raw["disp_mm"].to_numpy() - raw["disp_mm"].iloc[0]
        f = raw["force_N"].to_numpy()
        F_max = float(np.max(f))

        # Toe correction (D953 Appendix X1): fit the 10-40 % F_max region and
        # project back to F = 0. That x-intercept is the machine take-up.
        fit_m = (f >= FIT_LO * F_max) & (f <= FIT_HI * F_max)
        if fit_m.sum() < 3:
            continue
        slope, intercept = np.polyfit(d[fit_m], f[fit_m], 1)
        d_corr = d - (-intercept / slope)

        # P at 4 % hole deformation (D953-19 3.2.2): pre-peak only, made
        # monotone so np.interp is valid.
        i_max = int(np.argmax(f))
        dc, fc = d_corr[:i_max + 1], f[:i_max + 1]
        keep = np.concatenate(([True], np.diff(np.maximum.accumulate(dc)) > 0))
        dc_m, fc_m = dc[keep], fc[keep]
        failed = float(np.max(dc_m)) < DEF_4PCT_MM
        P_4pct = F_max if failed else float(np.interp(DEF_4PCT_MM, dc_m, fc_m))

        area = t_mm * HOLE_D_MM       # D953-19 13.3
        rows.append({"exp": exp_code, "dir": d_int, "cid": spec_id,
                     "F_max_N": F_max,
                     "S_b":    P_4pct / area,    # D953-19 13.3 Eq. 1
                     "S_max":  F_max / area,     # D953-19 3.2.5
                     "failed": failed})
    return rows


def run_print_stats():
    spec_df = load_spec_sheet_raw()
    t_col = next((c for c in spec_df.columns if "thickness" in c.lower()), None)

    tensile_rows = load_tensile_stats_rows(spec_df)
    fmax = load_tensile_fmax()
    for r in tensile_rows:
        r["F_max_N"] = fmax.get(r["cid"])
    bearing_rows = load_bearing_rows(spec_df, t_col)

    #             (key,            header,        dec_m, dec_s, w_m, w_s)
    tensile_props = [
        ("E_GPa",          "E (GPa)",     2, 2, 5, 4),
        ("sigma_y_MPa",    "σ_y (MPa)", 1, 1, 5, 4),
        ("UTS_MPa",        "UTS (MPa)",   1, 1, 5, 4),
        ("eps_at_UTS_pct", "ε_UTS (%)", 2, 2, 4, 4),
        ("poisson_chord",  "ν_chord",    3, 3, 5, 5),
        ("F_max_N",        "F_max (N)",   1, 1, 6, 5),
    ]
    print_table("TENSILE MECHANICAL PROPERTIES — ASTM D638-14",
                tensile_rows, EXPOSURE_ORDER, T_DIRECTIONS,
                tensile_props, EXPOSURE_LABELS)

    bearing_props = [
        ("F_max_N", "F_max (N)",   1, 1, 6, 5),
        ("S_b",     "S_b (MPa)",   1, 1, 5, 4),
        ("S_max",   "S_max (MPa)", 1, 1, 5, 4),
    ]
    print_table("\nPIN-BEARING PROPERTIES — ASTM D953-19 Procedure A",
                bearing_rows, EXPOSURE_ORDER, B_DIRECTIONS,
                bearing_props, EXPOSURE_LABELS)

    n_failed = sum(r["failed"] for r in bearing_rows)
    if n_failed:
        print(f"\n  Note: {n_failed} coupon(s) failed before reaching 4% hole deformation.")
        print("        For those, P_4pct = F_max (conservative upper bound for S_b).")

    write_stats_sheet(OUT_STATS_XLSX, "Tensile",
                      build_stats_df(tensile_rows, EXPOSURE_ORDER, T_DIRECTIONS,
                                     tensile_props, EXPOSURE_LABELS))
    write_stats_sheet(OUT_STATS_XLSX, "Bearing",
                      build_stats_df(bearing_rows, EXPOSURE_ORDER, B_DIRECTIONS,
                                     bearing_props, EXPOSURE_LABELS))
    print(f"\nExported: {OUT_STATS_XLSX} [Tensile, Bearing]")


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    FIGS_ROOT.mkdir(parents=True, exist_ok=True)
    scalars = load_specimen_scalars()

    if DO_PER_COUPON_PLOTS:
        n_plotted = 0
        for cid in selected_coupons():
            props = load_props(cid, scalars)
            if props is None:
                print(f"[{cid}] no per-coupon CSV, no kept rows, or no specimen-sheet "
                      f"row — run Level 2 first")
                continue
            print(f"[{cid}]  E={props['E_GPa']:.2f} GPa  "
                  f"σ_y={props['sigma_y_MPa']:.1f} MPa  "
                  f"UTS={props['UTS_MPa']:.1f} MPa  "
                  f"ε_UTS={props['eps_at_UTS'] * 100:.2f}%  "
                  f"ν_chord={props['poisson_chord']:.3f}")
            fig_dir = FIGS_ROOT / cid
            fig_dir.mkdir(parents=True, exist_ok=True)
            plot_stress_strain(cid, props, fig_dir)
            plot_poisson(cid, props, fig_dir)
            n_plotted += 1
        print(f"\n{n_plotted} coupon(s) plotted → {FIGS_ROOT}")

    if DO_GROUP_PLOTS:
        mts_coupons = group_load_mts_coupons()
        if mts_coupons:
            group_plot_mts_force_displacement(mts_coupons)
        else:
            print("No MTS coupons found; skipping MTS group figure.")

        dic_coupons = group_load_dic_coupons(scalars)
        if dic_coupons:
            group_plot_dic_curves(dic_coupons)
            group_plot_dic_property_scatter(dic_coupons)
            group_plot_dic_peak_strength(dic_coupons)
        else:
            print("No DIC coupons found; skipping DIC group figures.")

    if DO_PRINT_STATS:
        run_print_stats()

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
