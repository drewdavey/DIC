#!/usr/bin/env python3
"""
DIC_Level3.py  —  FSR Tensile Coupons
======================================
Plot- and stats-only. Reads Level-2's per-frame curve CSVs and the scalar
mechanical properties Level-2 already wrote into FSR-SpecimenTesting.xlsx
(the single source of truth for tensile scalars — nothing tensile is
recomputed here), and produces per-coupon plots, group overlay/summary
plots, and mean ± std (CV%) property tables. Pin-bearing statistics are
computed here directly from the raw MTS files (bearing has no DIC/Level-2
step of its own).

INPUT per coupon
  <DIC_DIR>/<coupon_id>_L2.csv   step, eps, sig, eps_t, i_uts, eps_raw, sig_raw
                                 (written by DIC_Level2.py)
  FSR-SpecimenTesting.xlsx       scalar properties (E, toe strain, yield
                                 stress/strain, UTS, strain at UTS, Poisson's
                                 ratio), one row per coupon, matched by
                                 "Specimen ID" — written by DIC_Level2.py
  <MTS_DIR>/P01-T*.txt           raw MTS tensile force/displacement (group MTS plot + F_max stats)
  <MTS_DIR>/P01-B*.txt           raw MTS bearing force/displacement (bearing stats only)

OUTPUT
  {FIGS_ROOT}/{coupon_id}/stress_strain_DIC.png   per-coupon σ–ε (toe-corrected, to UTS)
  {FIGS_ROOT}/{coupon_id}/poisson_DIC.png         per-coupon −ε_xx vs ε_yy (to UTS)
  {FIGS_ROOT}/tensile_mts_FD.png                  group force vs. displacement (raw MTS)
  {FIGS_ROOT}/tensile_curves_DIC.png              group σ–ε overlay
  {FIGS_ROOT}/tensile_summary_DIC.png             group property scatter (UTS, E)
  {FIGS_ROOT}/tensile_peak_strength_DIC.png       group UTS by exposure
  stdout                                          tensile (D638) + bearing (D953) stat tables
  P01_MechanicalStats.xlsx                        same stats, Tensile + Bearing sheets
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
DIC_DIR = FIGS_ROOT.parent / "DIC"   # per-frame _L2.csv files written by Level 2
MTS_DIR = FIGS_ROOT.parent / "MTS"   # raw MTS .txt files (group MTS plot only)
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

DO_PER_COUPON_PLOTS = True
DO_GROUP_PLOTS      = True   # MTS force-disp + DIC curve/summary/peak-strength group figures
DO_PRINT_STATS      = True   # mean ± std (CV%) tables (tensile + bearing) + P01_MechanicalStats.xlsx

# TODO: restore P01-TCL45-01 once backup DIC data is loaded — excluded because
#       its toe correction is anomalously large, distorting group curves/stats.
#       Per-coupon plots above still include it; only the group plots below don't.
DIC_EXCLUDE = {"P01-TCL45-01"}

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
EXPOSURE_ORDER    = ["CL", "UV", "SW", "IS"]
EXPOSURE_COLORS   = {"CL": "#1f77b4", "SW": "#17becf", "UV": "#ff7f0e", "IS": "#2ca02c"}
EXPOSURE_LABELS   = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
DIRECTION_MARKERS = {"00": "o", "45": "s", "90": "^"}
DIR_STYLES        = {0: "-", 45: "--", 90: ":"}   # MTS group plot line styles, int-keyed

MTS_HEADERS = 8

# =============================================================================
# PRINT-STATS SETTINGS  (ASTM D638-14 tensile + ASTM D953-19 Procedure A bearing)
# =============================================================================
T_DIRECTIONS = [0, 45, 90]   # tensile
B_DIRECTIONS = [0, 90]       # bearing (no 45° coupons)

# D953 geometry
HOLE_D_MM   = 0.5625 * 25.4    # 14.29 mm reamed hole diameter
DEF_4PCT_MM = 0.04 * HOLE_D_MM # 0.572 mm — 4% hole deformation threshold
FIT_LO, FIT_HI = 0.10, 0.40    # toe-correction linear fit window (fraction of F_max)

# BGM00 substitute (fixture malfunction on BCL00-01)
SUBST_SRC   = "P01-BGM00"
SUBST_LABEL = "P01-BCL00-01"
SKIP_BEAR   = {"P01-BCL00-01", "P01-BCL00-01-TEST"}

OUT_STATS_XLSX = FIGS_ROOT.parent / "P01_MechanicalStats.xlsx"

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
        if exp_code not in EXPOSURE_COLORS or d_int not in {0, 45, 90}:
            continue

        raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                          names=["disp_mm", "force_N", "output_V", "time_s"],
                          encoding="utf-8-sig", on_bad_lines="skip")
        raw = raw.apply(pd.to_numeric, errors="coerce").dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue

        d = raw["disp_mm"].to_numpy() - raw["disp_mm"].iloc[0]
        f = raw["force_N"].to_numpy() / 1000
        i_peak = int(np.argmax(f))

        coupons.append({"exp": exp_code, "dir": d_int,
                        "d": d, "f": f, "i_peak": i_peak})
        print(f"[{stem}]  peak {f[i_peak]:.2f} kN")

    print(f"\n{len(coupons)} MTS coupons loaded")
    return coupons

def group_plot_mts_force_displacement(coupons):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        ax.plot(c["d"], c["f"],
                color=EXPOSURE_COLORS[c["exp"]],
                ls=DIR_STYLES[c["dir"]],
                lw=0.9, alpha=0.75)
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
    for d_int, ls in DIR_STYLES.items():
        handles.append(mlines.Line2D([], [], color="k", ls=ls, lw=1.2,
                                     label=f"{d_int}°"))
    handles.append(plt.scatter([], [], color="k", s=30, label="Peak force"))

    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)
    fig.tight_layout()

    out = FIGS_ROOT / "tensile_mts_FD.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# =============================================================================
# GROUP PLOTS — DIC-derived σ-ε and scalar properties
# =============================================================================
def group_load_dic_coupons(scalars: pd.DataFrame):
    coupons = []
    for cid in selected_coupons():
        if cid in DIC_EXCLUDE:
            continue
        l2 = find_l2(cid)
        if l2 is None or cid not in scalars.index:
            continue
        exp, d_str = parse_id(cid)
        df = pd.read_csv(l2)
        i_uts = int(df["i_uts"].iloc[0])
        sl = slice(0, i_uts + 1)
        ef = df["eps"].to_numpy()[sl]
        sf = df["sig"].to_numpy()[sl]

        row = scalars.loc[cid]
        coupons.append({
            "cid": cid,
            "exp": exp,
            "d_str": d_str,
            "d_int": int(d_str),
            "eps_plot": ef,
            "sig_plot": sf,
            "UTS_MPa": float(row["UTS_MPa"]),
            "E_GPa": float(row["E_GPa"]),
        })
        print(f"[{cid}]  UTS={coupons[-1]['UTS_MPa']:.1f} MPa  "
              f"E={coupons[-1]['E_GPa']:.2f} GPa")

    print(f"\n{len(coupons)} DIC coupons loaded (group plots)\n")
    return coupons

def group_plot_dic_curves(coupons):
    exp_active = [e for e in EXPOSURE_COLORS if EXPOSURES.get(e)]
    dirs_present = sorted({c["d_int"] for c in coupons})
    if not dirs_present:
        print("No DIC coupons found; skipping group curves.")
        return

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        col = EXPOSURE_COLORS[c["exp"]]
        ax.plot(c["eps_plot"] * 100, c["sig_plot"],
                color=col, ls="-", lw=0.9, alpha=0.75)
        ax.scatter(c["eps_plot"][-1] * 100, c["sig_plot"][-1],
                   color=col, marker=DIRECTION_MARKERS[c["d_str"]],
                   s=35, zorder=5)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Axial Strain (%)", fontsize=12)
    ax.set_ylabel("Stress (MPa)", fontsize=12)
    ax.set_title("P01 — Tensile: Stress vs. Strain", fontsize=13)
    ax.grid(alpha=0.25, ls="--")

    handles = [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
               for e in exp_active]
    for d_int in dirs_present:
        handles.append(mlines.Line2D([], [], color="k", ls="None",
                                     marker=DIRECTION_MARKERS[f"{d_int:02d}"], markersize=6,
                                     label=f"{d_int}°"))
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)

    fig.tight_layout()
    out = FIGS_ROOT / "tensile_curves_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

def group_plot_dic_property_scatter(coupons):
    dirs_present = sorted({c["d_int"] for c in coupons})
    if not dirs_present:
        print("No DIC coupons found; skipping property scatter.")
        return

    exp_active = [e for e in EXPOSURE_COLORS if EXPOSURES.get(e)]
    prop_keys = ["UTS_MPa", "E_GPa"]
    ylabels = ["UTS (MPa)", "E (GPa)"]
    titles = ["Ultimate Tensile Strength", "Young's Modulus"]
    airtech = [AIRTECH_UTS, AIRTECH_E]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, key, ylabel, title, refs in zip(axes, prop_keys, ylabels, titles, airtech):
        for c in coupons:
            ei = exp_active.index(c["exp"])
            x = c["d_int"] + (ei - (len(exp_active) - 1) / 2) * 1.5
            ax.scatter(x, c[key],
                       color=EXPOSURE_COLORS[c["exp"]],
                       marker=DIRECTION_MARKERS.get(c["d_str"], "o"),
                       s=55, zorder=5)
        for d_int, ref in refs.items():
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
    exp_active = [e for e in EXPOSURE_ORDER if EXPOSURES.get(e)]
    if not exp_active:
        print("No DIC exposures enabled; skipping peak strength plot.")
        return

    x_pos = np.arange(len(EXPOSURE_ORDER))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axvspan(x_pos[0] - 0.5, x_pos[1] + 0.5, color="#f7f5e8", alpha=0.35)
    ax.axvspan(x_pos[2] - 0.5, x_pos[3] + 0.5, color="#e8f7f5", alpha=0.35)

    for c in coupons:
        x = x_pos[EXPOSURE_ORDER.index(c["exp"])]
        ax.scatter(x, c["UTS_MPa"],
                   color=EXPOSURE_COLORS[c["exp"]],
                   marker=DIRECTION_MARKERS.get(c["d_str"], "o"),
                   s=70, edgecolor="black", linewidth=0.4, zorder=5)

    for d_str in ["00", "45", "90"]:
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

    dir_handles = [mlines.Line2D([], [], color="black", marker=DIRECTION_MARKERS[d], linestyle="None",
                                 markersize=7, label=f"{d}°")
                   for d in ["00", "45", "90"] if DIRECTIONS.get(d)]
    ax.legend(handles=dir_handles, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)

    fig.tight_layout()
    out = FIGS_ROOT / "tensile_peak_strength_DIC.png"
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


# =============================================================================
# PRINT STATS — mean ± std (CV%) tables, tensile (D638) + bearing (D953)
# =============================================================================
def load_spec_sheet_raw() -> pd.DataFrame:
    """Full specimen sheet with original Excel headers (not renamed), indexed
    by Specimen ID — needed here for coupon thickness (bearing area) as well
    as the tensile scalar columns."""
    df = pd.read_excel(SPECIMEN_SHEET)
    if "Print ID" in df.columns:
        df = df[df["Print ID"] == PRINTS[0]]
    if "Specimen ID" in df.columns:
        df = df.set_index("Specimen ID")
    return df

def load_tensile_stats_rows(spec_df: pd.DataFrame) -> list[dict]:
    rows = []
    for exp in EXPOSURE_ORDER:
        for d in T_DIRECTIONS:
            for rep in REPLICATES:
                cid = f"{PRINTS[0]}-T{exp}{d:02d}-{rep}"
                if cid in DIC_EXCLUDE or cid not in spec_df.index:
                    continue
                spec_row = spec_df.loc[cid]
                # Individual properties (e.g. sigma_y for a brittle 90° coupon
                # with no offset-yield crossing) may legitimately be NaN —
                # stat() already drops NaN per-property, so only skip the
                # whole coupon if E_GPa itself (needed for every coupon) is missing.
                cols = {k: spec_row.get(v) for k, v in SPECIMEN_SHEET_COLUMNS.items()}
                cols = {k: (float(v) if pd.notna(v) else np.nan) for k, v in cols.items()}
                if np.isnan(cols["E_GPa"]):
                    continue
                rows.append({
                    "exp": exp,
                    "dir": d,
                    "cid": cid,
                    "E_GPa":          cols["E_GPa"],
                    "sigma_y_MPa":    cols["sigma_y_MPa"],
                    "UTS_MPa":        cols["UTS_MPa"],
                    "eps_at_UTS_pct": cols["eps_at_UTS"] * 100.0,
                    "poisson_chord":  cols["poisson_chord"],
                })
    return rows

def load_tensile_fmax() -> dict[str, float]:
    """Peak load per tensile coupon, read directly from the raw MTS .txt
    files (independent of the DIC pipeline)."""
    fmax: dict[str, float] = {}
    for fp in sorted(MTS_DIR.glob(f"{PRINTS[0]}-T*.txt")):
        stem = re.sub(r"\.txt$", "", fp.name, flags=re.IGNORECASE)
        stem_clean = re.sub(r"-TEST$", "", stem, flags=re.IGNORECASE)
        parts = stem_clean.split("-")
        if len(parts) < 3 or not parts[1].upper().startswith("T") or len(parts[1]) < 5:
            continue
        key = parts[1].upper()
        exp_code = key[1:-2]
        try:
            d_int = int(key[-2:])
        except ValueError:
            continue
        if exp_code not in EXPOSURE_LABELS or d_int not in T_DIRECTIONS:
            continue
        rep_str = parts[2]
        if not rep_str.isdigit():
            continue
        cid = f"{PRINTS[0]}-T{exp_code}{d_int:02d}-{int(rep_str):02d}"

        raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                          names=["disp_mm", "force_N", "output_V", "time_s"],
                          encoding="utf-8-sig", on_bad_lines="skip")
        raw = raw.apply(pd.to_numeric, errors="coerce").dropna(subset=["force_N"])
        if len(raw) < 10:
            continue
        fmax[cid] = float(np.max(raw["force_N"].to_numpy()))
    return fmax

def get_t_mm(spec_df: pd.DataFrame, t_col: str | None, sid: str) -> float | None:
    if t_col is None or sid not in spec_df.index:
        return None
    v = spec_df.loc[sid, t_col]
    if isinstance(v, pd.Series):
        v = v.iloc[0]
    return float(v) * 25.4 if pd.notna(v) else None

def parse_bearing_stem(stem: str):
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
    base = f"{PRINTS[0]}-B{exp_code}{d_int:02d}"
    rep_part = ""
    if len(parts) >= 3:
        p3 = re.sub(r"^TEST", "", parts[2], flags=re.IGNORECASE)
        if p3.isdigit():
            rep_part = p3
    spec_id = f"{base}-{int(rep_part):02d}" if rep_part else base
    return exp_code, d_int, spec_id

def load_bearing_rows(spec_df: pd.DataFrame, t_col: str | None) -> list[dict]:
    rows = []
    for fp in sorted(MTS_DIR.glob(f"{PRINTS[0]}-B*.txt")):
        stem = re.sub(r"\.txt$", "", fp.name, flags=re.IGNORECASE)
        if stem.upper() in {s.upper() for s in SKIP_BEAR}:
            continue
        exp_code, d_int, spec_id = parse_bearing_stem(stem)
        if exp_code is None:
            continue
        t_mm = get_t_mm(spec_df, t_col, spec_id)
        if t_mm is None:
            continue

        raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                          names=["disp_mm", "force_N", "output_V", "time_s"],
                          encoding="utf-8-sig", on_bad_lines="skip")
        raw = raw.apply(pd.to_numeric, errors="coerce").dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue

        d = raw["disp_mm"].to_numpy() - raw["disp_mm"].iloc[0]
        f = raw["force_N"].to_numpy()
        F_max = float(np.max(f))

        # Toe correction (D953 Appendix X1 tangent-method equivalent):
        # Fit a line to the 10-40% F_max region and project back to F = 0.
        # That x-intercept is the toe offset (machine take-up).
        fit_m = (f >= FIT_LO * F_max) & (f <= FIT_HI * F_max)
        if fit_m.sum() < 3:
            continue
        slope, intercept = np.polyfit(d[fit_m], f[fit_m], 1)
        d_corr = d - (-intercept / slope)

        # P at 4% hole deformation (D953-19 §3.2.2): only pre-peak data,
        # enforced monotone for interp.
        i_max = int(np.argmax(f))
        dc, fc = d_corr[:i_max + 1], f[:i_max + 1]
        cum_m  = np.maximum.accumulate(dc)
        keep   = np.concatenate(([True], np.diff(cum_m) > 0))
        dc_m, fc_m = dc[keep], fc[keep]
        failed = float(np.max(dc_m)) < DEF_4PCT_MM
        P_4pct = F_max if failed else float(np.interp(DEF_4PCT_MM, dc_m, fc_m))

        area = t_mm * HOLE_D_MM   # D953-19 §13.3
        rows.append({
            "exp":     exp_code,
            "dir":     d_int,
            "cid":     spec_id,
            "F_max_N": F_max,
            "S_b":     P_4pct / area,   # D953-19 §13.3 Eq. 1
            "S_max":   F_max  / area,   # D953-19 §3.2.5
            "failed":  failed,
        })
    return rows

def stat(vals: list) -> tuple[float, float, float]:
    """Return (mean, std, cv%) from a list; NaN for invalid/empty groups."""
    arr = np.array([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan, np.nan
    m = float(np.mean(arr))
    s = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
    cv = 100.0 * s / abs(m) if abs(m) > 1e-12 else 0.0
    return m, s, cv

def stat_cell(m, s, cv, w_m=6, w_s=5, dec_m=2, dec_s=2) -> str:
    """Format one statistic cell: 'mean ± std (CV%)'."""
    if not np.isfinite(m):
        return f"{'NaN':>{w_m}} ± {'NaN':>{w_s}} (  N/A)"
    return f"{m:{w_m}.{dec_m}f} ± {s:{w_s}.{dec_s}f} ({cv:4.1f}%)"

def stat_row(label_exp, label_dir, n, cells: list[str]) -> str:
    return f"{label_exp:<12s} {label_dir:>6s}  {n:3d}   " + "   ".join(cells)

def print_table(title: str, rows_data: list[dict], exp_order: list[str],
                directions: list[int], props: list[tuple],
                exp_labels: dict[str, str]) -> None:
    """props : list of (key, header, dec_m, dec_s, w_m, w_s)"""
    col_w = 20   # width of each stat column (including padding)
    id_w  = 24   # width of 'Exposure Dir n' prefix

    bar  = "=" * (id_w + col_w * len(props))
    sep  = "-" * (id_w + col_w * len(props))
    hdrs = "".join(f"{p[1]:^{col_w}}" for p in props)
    cv_h = "".join(f"{'mean ± std (CV%)':^{col_w}}" for p in props)

    print(f"\n{title}")
    print(bar)
    print(f"{'Exposure':<12s} {'Dir':>6s}  {'n':3s}   " + hdrs)
    print(f"{'':12s} {'':6s}  {'':3s}   " + cv_h)
    print(sep)

    def print_group(label_exp, label_dir, subset):
        n = len(subset)
        if n == 0:
            return
        cells = []
        for key, _, dec_m, dec_s, w_m, w_s in props:
            vals = [r[key] for r in subset]
            m, s, cv = stat(vals)
            cells.append(stat_cell(m, s, cv, w_m, w_s, dec_m, dec_s).ljust(col_w - 3))
        print(stat_row(label_exp, label_dir, n, cells))

    for exp in exp_order:
        for d in directions:
            subset = [r for r in rows_data if r["exp"] == exp and r["dir"] == d]
            print_group(exp_labels.get(exp, exp), f"{d}°", subset)

    print(sep)
    for d in directions:
        subset = [r for r in rows_data if r["dir"] == d]
        print_group("All exps", f"{d}°", subset)

    print(bar)

def build_stats_df(rows_data: list[dict], exp_order: list[str],
                   directions: list[int], props: list[tuple],
                   exp_labels: dict[str, str]) -> pd.DataFrame:
    """Build a DataFrame with mean/std/CV columns for each property."""
    records = []
    for exp in exp_order:
        for d in directions:
            subset = [r for r in rows_data if r["exp"] == exp and r["dir"] == d]
            if not subset:
                continue
            rec = {"Exposure": exp_labels.get(exp, exp), "Dir (deg)": d, "n": len(subset)}
            for key, header, *_ in props:
                vals = [r[key] for r in subset]
                m, s, cv = stat(vals)
                rec[f"{header} mean"]   = round(m, 4) if np.isfinite(m) else None
                rec[f"{header} std"]    = round(s, 4) if np.isfinite(s) else None
                rec[f"{header} CV (%)"] = round(cv, 2) if np.isfinite(cv) else None
            records.append(rec)
    for d in directions:
        subset = [r for r in rows_data if r["dir"] == d]
        if not subset:
            continue
        rec = {"Exposure": "All exps", "Dir (deg)": d, "n": len(subset)}
        for key, header, *_ in props:
            vals = [r[key] for r in subset]
            m, s, cv = stat(vals)
            rec[f"{header} mean"]   = round(m, 4) if np.isfinite(m) else None
            rec[f"{header} std"]    = round(s, 4) if np.isfinite(s) else None
            rec[f"{header} CV (%)"] = round(cv, 2) if np.isfinite(cv) else None
        records.append(rec)
    return pd.DataFrame(records)

def run_print_stats():
    spec_df = load_spec_sheet_raw()
    t_col = next((c for c in spec_df.columns if "thickness" in c.lower()), None)

    tensile_rows = load_tensile_stats_rows(spec_df)
    fmax = load_tensile_fmax()
    for r in tensile_rows:
        r["F_max_N"] = fmax.get(r["cid"])

    bearing_rows = load_bearing_rows(spec_df, t_col)

    #  (key,            header,        dec_m, dec_s, w_m, w_s)
    tensile_props = [
        ("E_GPa",          "E (GPa)",   2, 2, 5, 4),
        ("sigma_y_MPa",    "σ_y (MPa)", 1, 1, 5, 4),
        ("UTS_MPa",        "UTS (MPa)", 1, 1, 5, 4),
        ("eps_at_UTS_pct", "ε_UTS (%)", 2, 2, 4, 4),
        ("poisson_chord",  "ν_chord",   3, 3, 5, 5),
        ("F_max_N",        "F_max (N)", 1, 1, 6, 5),
    ]
    print_table(
        title="TENSILE MECHANICAL PROPERTIES — ASTM D638-14",
        rows_data=tensile_rows, exp_order=EXPOSURE_ORDER, directions=T_DIRECTIONS,
        props=tensile_props, exp_labels=EXPOSURE_LABELS,
    )

    bearing_props = [
        ("F_max_N", "F_max (N)",   1, 1, 6, 5),
        ("S_b",     "S_b (MPa)",   1, 1, 5, 4),
        ("S_max",   "S_max (MPa)", 1, 1, 5, 4),
    ]
    print_table(
        title="\nPIN-BEARING PROPERTIES — ASTM D953-19 Procedure A",
        rows_data=bearing_rows, exp_order=EXPOSURE_ORDER, directions=B_DIRECTIONS,
        props=bearing_props, exp_labels=EXPOSURE_LABELS,
    )

    n_failed = sum(r["failed"] for r in bearing_rows)
    if n_failed:
        print(f"\n  Note: {n_failed} coupon(s) failed before reaching 4% hole deformation.")
        print( "        For those, P_4pct = F_max (conservative upper bound for S_b).")

    df_tensile = build_stats_df(tensile_rows, EXPOSURE_ORDER, T_DIRECTIONS,
                                 tensile_props, EXPOSURE_LABELS)
    df_bearing = build_stats_df(bearing_rows, EXPOSURE_ORDER, B_DIRECTIONS,
                                 bearing_props, EXPOSURE_LABELS)

    with pd.ExcelWriter(OUT_STATS_XLSX, engine="openpyxl") as writer:
        df_tensile.to_excel(writer, sheet_name="Tensile", index=False)
        df_bearing.to_excel(writer, sheet_name="Bearing", index=False)
    print(f"\nExported: {OUT_STATS_XLSX}")


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
        print(f"\n{n_plotted} coupon(s) plotted -> {FIGS_ROOT}")

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

    print(f"\nDone. {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
