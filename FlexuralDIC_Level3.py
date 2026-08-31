#!/usr/bin/env python3
"""
FlexuralDIC_Level3.py  —  FSR Flexural Coupons (ASTM D790, 3-point bend)
========================================================================
Plots and statistics only. Reads Level-2's per-frame CSVs and the "Flex ..."
scalar columns Level 2 wrote into the specimen sheet; nothing is recomputed
here. See README.md for the method and for what the three strain channels mean.

INPUT
  <DIC_DIR>/<coupon_id>.csv    per-frame record, Level 1 + Level 2 columns
  FSR-SpecimenTesting.csv      D790 scalars, one row per coupon
  <MTS_DIR>/<coupon_id>.txt    raw MTS force/displacement (group MTS plot only)

OUTPUT
  <FIGS_ROOT>/<coupon_id>/flexural_stress_strain_DIC.png    stress vs strain
  <FIGS_ROOT>/<coupon_id>/flexural_strain_channels_DIC.png  the three channels
  <FIGS_ROOT>/flexural_mts_FD.png             group force vs displacement
  <FIGS_ROOT>/flexural_curves_DIC.png         group stress-strain overlay
  <FIGS_ROOT>/flexural_summary_DIC.png        group property scatter
  <FIGS_ROOT>/flexural_peak_strength_DIC.png  group strength by exposure
  stdout                                      D790 stat table
  P01_MechanicalStats.xlsx                    same stats, "Flexural" sheet
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
EXPOSURES  = {"CL": True, "IS": True}
DIRECTIONS = {"00": True, "90": True}
REPLICATES = ["01", "02", "03"]

DO_PER_COUPON_PLOTS = True
DO_GROUP_PLOTS      = True
DO_PRINT_STATS      = True

DIC_EXCLUDE: set[str] = set()

# =============================================================================
# ANALYSIS  — must match FlexuralDIC_Level1 and FlexuralDIC_Level2
# =============================================================================
MODULUS_STRAIN_RANGE = (0.0005, 0.003)   # used only to draw the tangent line
MTS_HEADERS = 8

# Sheet headers written by FlexuralDIC_Level2, inverted back to property keys.
SPECIMEN_SHEET_COLUMNS = {
    "sigma_fM_MPa":         "Flex Strength (MPa)",
    "P_max_N":              "Flex Peak Load (N)",
    "Ef_curvature_GPa":     "Flex Ef curvature (GPa)",
    "Ef_deflection_GPa":    "Flex Ef deflection (GPa)",
    "Ef_crosshead_GPa":     "Flex Ef crosshead (GPa)",
    "Ef_MI_kappa_GPa":      "Flex Ef M/I-kappa (GPa)",
    "eps_at_max_curvature": "Flex Strain at Max (curvature)",
    "profile_r2_median":    "Flex Profile R2 (median)",
}

# =============================================================================
# DISPLAY  — same palette as mts_plots.py and TensileDIC_Level3
# =============================================================================
EXPOSURE_ORDER    = ["CL", "IS"]
EXPOSURE_LABELS   = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
EXPOSURE_COLORS   = {"CL": "#2a78d6", "UV": "#eb6834", "SW": "#1baf7a", "IS": "#4a3aa7"}
DIRECTION_MARKERS = {"00": "o", "45": "s", "90": "^"}
DIR_STYLES        = {0: "-", 45: (0, (6.5, 3.5)), 90: (0, (1.5, 4.0))}
CHANNEL_COLORS    = {"curvature": "#2a78d6", "deflection": "#eb6834",
                     "crosshead": "#1baf7a"}


# =============================================================================
# HELPERS
# =============================================================================
def coupon_id(p, e, d, r): return f"{p}-F{e}{d}-{r}"

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
    props = {}
    for key in SPECIMEN_SHEET_COLUMNS:
        v = row[key] if key in row.index else np.nan
        props[key] = float(v) if pd.notna(v) else np.nan
    props["_sig"] = curve["stress_MPa"].to_numpy()
    for name in ("curvature", "deflection", "crosshead"):
        col = f"eps_{name}"
        props["_eps_" + name] = (curve[col].to_numpy() if col in curve.columns
                                 else np.full(len(curve), np.nan))
    props["_i_max"] = int(np.nanargmax(props["_sig"]))
    return props

def find_force_baseline(d, f):
    """Level of the leading flat run in f, or 0.0 if there isn't one.

    The flexural records open on ~860 N of loading nose hanging on an un-tared
    cell. Same detector as FlexuralDIC_Level1 and mts_plots. Expects f already
    truncated at peak and sign-flipped positive.
    """
    n = len(f)
    if n < 200:
        return 0.0
    f_max = f[-1]
    rise = f_max - f[0]
    band = (f >= f[0] + 0.40 * rise) & (f <= f[0] + 0.70 * rise)
    if band.sum() < 10:
        return 0.0
    load_slope = float(np.polyfit(d[band], f[band], 1)[0])
    if load_slope <= 0:
        return 0.0
    w = max(25, n // 100)
    nb = n // w
    if nb < 8:
        return 0.0
    fs = f[:nb * w].reshape(nb, w).mean(axis=1)
    ds = d[:nb * w].reshape(nb, w).mean(axis=1)
    above = np.flatnonzero(np.gradient(fs, ds) > 0.40 * load_slope)
    if above.size == 0:
        return 0.0
    kb = int(above[0])
    k = kb * w
    if k < 50 or (d[k] - d[0]) < 0.2 or kb < 2:
        return 0.0
    baseline = float(np.median(f[:k]))
    if np.ptp(fs[:kb]) > 0.10 * (f_max - baseline):
        return 0.0
    if abs(baseline) < 0.05 * (f_max - baseline):
        return 0.0
    return baseline


# =============================================================================
# PER-COUPON PLOTS
# =============================================================================
def plot_stress_strain(cid, props, fig_dir):
    """Flexural stress vs curvature strain, truncated at peak stress."""
    exp, d_str = parse_id(cid)
    sl = slice(0, props["_i_max"] + 1)
    eps = props["_eps_curvature"][sl] * 100
    sig = props["_sig"][sl]

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(eps, sig, lw=1.4, color=EXPOSURE_COLORS.get(exp, "#333"), label=cid)

    E_MPa = props["Ef_curvature_GPa"] * 1000.0
    if np.isfinite(E_MPa):
        x = np.array([0.0, MODULUS_STRAIN_RANGE[1] * 1.5])
        ax.plot(x * 100, E_MPa * x, "k--", lw=0.8, alpha=0.7,
                label=f"E_f = {props['Ef_curvature_GPa']:.2f} GPa")
    ax.plot(props["eps_at_max_curvature"] * 100, props["sigma_fM_MPa"], "^",
            color="red", ms=8, zorder=5,
            label=f"σ_fM = {props['sigma_fM_MPa']:.1f} MPa")

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Flexural Strain (%)")
    ax.set_ylabel("Flexural Stress (MPa)")
    ax.set_title(f"{cid}  —  {EXPOSURE_LABELS.get(exp, exp)}, {int(d_str)}°")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    fig.tight_layout()
    out = fig_dir / "flexural_stress_strain_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_strain_channels(cid, props, fig_dir):
    """The three strain channels on one stress axis. Curvature is the reference;
    the gaps to deflection and crosshead are D790's small-deflection kinematics
    and the machine compliance."""
    sl = slice(0, props["_i_max"] + 1)
    sig = props["_sig"][sl]

    fig, ax = plt.subplots(figsize=(7, 4.8))
    for name in ("curvature", "deflection", "crosshead"):
        eps = props["_eps_" + name][sl]
        if not np.any(np.isfinite(eps)):
            continue
        E = props["Ef_" + name + "_GPa"]
        label = f"{name}  E_f = {E:.2f} GPa" if np.isfinite(E) else name
        ax.plot(eps * 100, sig, lw=1.2, color=CHANNEL_COLORS[name], label=label)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Flexural Strain (%)")
    ax.set_ylabel("Flexural Stress (MPa)")
    ax.set_title(f"{cid}  —  strain channels")
    ax.grid(alpha=0.25, linestyle="--")
    ax.legend(fontsize=8, framealpha=0.85, loc="best")
    fig.tight_layout()
    out = fig_dir / "flexural_strain_channels_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# GROUP PLOTS — raw MTS force/displacement
# =============================================================================
def group_load_mts_coupons():
    """Raw MTS records, flipped to the first quadrant with the tare removed —
    the same reduction FlexuralDIC_Level1 does."""
    coupons = []
    for cid in selected_coupons():
        fp = MTS_DIR / f"{cid}.txt"
        if not fp.exists():
            continue
        raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                          names=["disp_mm", "force_N", "output_V", "time_s"],
                          encoding="utf-8-sig", on_bad_lines="skip")
        raw = raw.apply(pd.to_numeric, errors="coerce").dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue
        d = -raw["disp_mm"].to_numpy()
        f = -raw["force_N"].to_numpy()
        d = d - d[0]
        i_peak = int(np.argmax(f))
        tare = find_force_baseline(d[:i_peak + 1], f[:i_peak + 1])
        f = (f - tare) / 1000.0
        exp, d_str = parse_id(cid)
        coupons.append({"cid": cid, "exp": exp, "dir": int(d_str),
                        "d": d, "f": f, "i_peak": i_peak})
        print(f"[{cid}]  peak {f[i_peak]:.2f} kN  (tare {tare:.0f} N removed)")
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
    ax.set_title("P01 — Flexural: Force vs. Displacement", fontsize=13)
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
    out = FIGS_ROOT / "flexural_mts_FD.png"
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
            "eps_plot": df["eps_curvature"].to_numpy()[sl],
            "sig_plot": sig[sl],
            "sigma_fM_MPa": float(row["sigma_fM_MPa"]),
            "Ef_curvature_GPa": float(row["Ef_curvature_GPa"]),
        })
        print(f"[{cid}]  σ_fM={coupons[-1]['sigma_fM_MPa']:.1f} MPa  "
              f"E_f={coupons[-1]['Ef_curvature_GPa']:.2f} GPa")
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
    ax.set_xlabel("Flexural Strain (%)", fontsize=12)
    ax.set_ylabel("Flexural Stress (MPa)", fontsize=12)
    ax.set_title("P01 — Flexural: Stress vs. Strain", fontsize=13)
    ax.grid(alpha=0.25, ls="--")

    exp_active = [e for e in EXPOSURE_ORDER if any(c["exp"] == e for c in coupons)]
    handles = [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
               for e in exp_active]
    for d_int in sorted({c["d_int"] for c in coupons}):
        handles.append(mlines.Line2D([], [], color="k", ls="None",
                                     marker=DIRECTION_MARKERS[f"{d_int:02d}"],
                                     markersize=6, label=f"{d_int}°"))
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)

    fig.tight_layout()
    out = FIGS_ROOT / "flexural_curves_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def group_plot_dic_property_scatter(coupons):
    dirs_present = sorted({c["d_int"] for c in coupons})
    exp_active = [e for e in EXPOSURE_ORDER if EXPOSURES.get(e)]
    panels = [("sigma_fM_MPa", "Flexural Strength (MPa)", "Flexural Strength"),
              ("Ef_curvature_GPa", "E_f (GPa)", "Flexural Modulus (curvature)")]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, panel in zip(axes, panels):
        key, ylabel, title = panel
        for c in coupons:
            ei = exp_active.index(c["exp"])
            x = c["d_int"] + (ei - (len(exp_active) - 1) / 2) * 1.5
            ax.scatter(x, c[key], color=EXPOSURE_COLORS[c["exp"]],
                       marker=DIRECTION_MARKERS.get(c["d_str"], "o"), s=55, zorder=5)
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
    out = FIGS_ROOT / "flexural_summary_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def group_plot_dic_peak_strength(coupons):
    x_pos = np.arange(len(EXPOSURE_ORDER))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        x = x_pos[EXPOSURE_ORDER.index(c["exp"])]
        ax.scatter(x, c["sigma_fM_MPa"], color=EXPOSURE_COLORS[c["exp"]],
                   marker=DIRECTION_MARKERS.get(c["d_str"], "o"),
                   s=70, edgecolor="black", linewidth=0.4, zorder=5)

    for d_str in DIRECTIONS:
        vals = [c["sigma_fM_MPa"] for c in coupons if c["d_str"] == d_str]
        if vals:
            ax.axhline(np.mean(vals), color="k", ls="--", lw=0.7, alpha=0.6, zorder=3)

    ax.set_xlim(-0.5, len(EXPOSURE_ORDER) - 0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([EXPOSURE_LABELS[e] for e in EXPOSURE_ORDER])
    ax.set_xlabel("Exposure Condition")
    ax.set_ylabel("Flexural Strength (MPa)")
    ax.set_title("P01 — Flexural: Max Stress")
    ax.grid(alpha=0.25, ls="--", axis="y")

    handles = [mlines.Line2D([], [], color="black", marker=DIRECTION_MARKERS[d],
                             linestyle="None", markersize=7, label=f"{d}°")
               for d in DIRECTIONS if DIRECTIONS.get(d)]
    ax.legend(handles=handles, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 1), borderaxespad=0)

    fig.tight_layout()
    out = FIGS_ROOT / "flexural_peak_strength_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# =============================================================================
# PRINT STATS — mean +/- std (CV%), ASTM D790-17
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
    """(label, subset) for each exposure x direction, then All exps per direction."""
    groups = []
    for exp in exp_order:
        for d in directions:
            subset = [r for r in rows_data if r["exp"] == exp and r["dir"] == d]
            groups.append((exp_labels.get(exp, exp), d, subset))
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
    TensileDIC_Level3 writes Tensile and Bearing into this same file."""
    if path.exists():
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="a",
                                if_sheet_exists="replace")
    else:
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="w")
    with writer:
        df.to_excel(writer, sheet_name=sheet, index=False)

def load_stats_rows(scalars):
    rows = []
    for cid in selected_coupons():
        if cid in DIC_EXCLUDE or cid not in scalars.index:
            continue
        row = scalars.loc[cid]
        if pd.isna(row.get("Ef_curvature_GPa")):
            continue
        exp, d_str = parse_id(cid)
        rec = {"exp": exp, "dir": int(d_str), "cid": cid}
        for key in ("sigma_fM_MPa", "P_max_N", "Ef_curvature_GPa",
                    "Ef_deflection_GPa", "Ef_crosshead_GPa"):
            v = row.get(key)
            rec[key] = float(v) if pd.notna(v) else np.nan
        v = row.get("eps_at_max_curvature")
        rec["eps_at_max_pct"] = float(v) * 100.0 if pd.notna(v) else np.nan
        rows.append(rec)
    return rows

def run_print_stats(scalars):
    rows = load_stats_rows(scalars)
    if not rows:
        print("No flexural scalars in the specimen sheet — run Level 2 first.")
        return

    #     (key,                header,           dec_m, dec_s, w_m, w_s)
    props = [
        ("Ef_curvature_GPa",  "E_f curv (GPa)",  2, 2, 5, 4),
        ("Ef_deflection_GPa", "E_f defl (GPa)",  2, 2, 5, 4),
        ("Ef_crosshead_GPa",  "E_f xhead (GPa)", 2, 2, 5, 4),
        ("sigma_fM_MPa",      "σ_fM (MPa)",  1, 1, 5, 4),
        ("eps_at_max_pct",    "ε_max (%)",     2, 2, 4, 4),
        ("P_max_N",           "P_max (N)",       1, 1, 6, 5),
    ]
    directions = sorted({r["dir"] for r in rows})
    print_table("FLEXURAL PROPERTIES — ASTM D790-17 Procedure A, 3-point",
                rows, EXPOSURE_ORDER, directions, props, EXPOSURE_LABELS)

    df = build_stats_df(rows, EXPOSURE_ORDER, directions, props, EXPOSURE_LABELS)
    write_stats_sheet(OUT_STATS_XLSX, "Flexural", df)
    print(f"\nExported: {OUT_STATS_XLSX} [Flexural]")


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
            print(f"[{cid}]  E_f={props['Ef_curvature_GPa']:.2f} GPa  "
                  f"σ_fM={props['sigma_fM_MPa']:.1f} MPa  "
                  f"ε_max={props['eps_at_max_curvature'] * 100:.2f}%  "
                  f"R²={props['profile_r2_median']:.4f}")
            fig_dir = FIGS_ROOT / cid
            fig_dir.mkdir(parents=True, exist_ok=True)
            plot_stress_strain(cid, props, fig_dir)
            plot_strain_channels(cid, props, fig_dir)
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
        run_print_stats(scalars)

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
