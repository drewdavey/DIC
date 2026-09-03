#!/usr/bin/env python3
"""
Bearing_Level3.py  —  FSR Pin-Bearing Coupons (ASTM D953-19 Procedure A)
========================================================================
Figures and statistics for the pin-bearing coupons. Bearing has no DIC step, so
there is no Level 1 or Level 2 to read from: this script does the whole job,
reducing the raw MTS files itself and then producing the same kind of output
the tensile and flexural Level 3s do.

Per coupon it toe-corrects the load-deflection record, reads the load at 4 %
hole deformation, and divides by the projected bearing area:

    area   = t * D_hole                     D953-19 §13.3
    S_b    = P(4 % deformation) / area      D953-19 §13.3 Eq. 1
    S_max  = F_max / area                   D953-19 §3.2.5

A coupon that breaks before reaching 4 % deformation has no P_4pct; F_max is
used instead, which is a conservative upper bound, and it is flagged in the
table and drawn with a different marker.

INPUT
  <MTS_DIR>/P01-B*.txt         raw MTS bearing force/displacement
  FSR-SpecimenTesting.csv      coupon thickness, for the bearing area

OUTPUT
  <FIGS_ROOT>/bearing_curves.png         load-deflection, one panel per direction
  <FIGS_ROOT>/bearing_curves_single.png  stress-deflection, all coupons together
  <FIGS_ROOT>/bearing_summary.png        S_b and S_max scatter
  <FIGS_ROOT>/bearing_max_stress.png     S_max by exposure
  stdout                                 D953 stat table
  P01_MechanicalStats.xlsx               same stats, "Bearing" sheet
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
DIRECTIONS = [0, 90]         

DO_GROUP_PLOTS = True
DO_PRINT_STATS = True

# Fixture malfunction on BCL00-01 (coupon expanded inside the hole, no
# bushings); P01-BGM00 is the replacement Control 0 deg run and is reported
# under that coupon's ID.
SUBST_SRC   = "P01-BGM00"
SUBST_LABEL = "P01-BCL00-01"
SKIP_BEAR   = {"P01-BCL00-01", "P01-BCL00-01-TEST"}

# =============================================================================
# ANALYSIS  — ASTM D953-19 Procedure A
# =============================================================================
IN2MM       = 25.4
HOLE_D_MM   = 0.5625 * IN2MM      # 14.29 mm reamed hole
DEF_4PCT_MM = 0.04 * HOLE_D_MM    # 0.572 mm, the 4 % hole-deformation point
FIT_LO, FIT_HI = 0.10, 0.40       # toe-correction fit window, fraction of F_max
MTS_HEADERS = 8

# =============================================================================
# DISPLAY  — same palette as mts_plots.py and the other Level 3s
# =============================================================================
EXPOSURE_ORDER    = ["CL", "UV", "SW", "IS"]
EXPOSURE_LABELS   = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
EXPOSURE_COLORS   = {"CL": "#2a78d6", "UV": "#eb6834", "SW": "#1baf7a", "IS": "#4a3aa7"}
DIRECTION_MARKERS = {0: "o", 90: "^"}


# =============================================================================
# HELPERS
# =============================================================================
def read_specimen_csv():
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc)
        except UnicodeDecodeError:
            continue
    raise SystemExit(f"{SPECIMEN_CSV.name}: not decodable as {'/'.join(CSV_ENCODINGS)}")


def load_spec_sheet():
    """Specimen sheet indexed by Specimen ID, for coupon thickness."""
    df = read_specimen_csv()
    if "Print ID" in df.columns:
        df = df[df["Print ID"] == PRINTS[0]]
    return df.set_index("Specimen ID")


def get_t_mm(spec_df, t_col, sid):
    """Coupon thickness in mm, or None."""
    if t_col is None or sid not in spec_df.index:
        return None
    v = spec_df.loc[sid, t_col]
    if isinstance(v, pd.Series):
        v = v.iloc[0]
    return float(v) * IN2MM if pd.notna(v) else None


def read_mts_txt(fp):
    """MTS .txt: 8-line header, tab separated, cols disp_mm force_N output_V time_s."""
    raw = pd.read_csv(fp, sep="\t", skiprows=MTS_HEADERS, header=None,
                      names=["disp_mm", "force_N", "output_V", "time_s"],
                      encoding="utf-8-sig", on_bad_lines="skip")
    return raw.apply(pd.to_numeric, errors="coerce")


def parse_bearing_stem(stem):
    """Return (exp_code, d_int, spec_id) from a filename stem, or (None, None, None)."""
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
    if exp_code not in EXPOSURE_LABELS or d_int not in DIRECTIONS:
        return None, None, None
    rep = ""
    if len(parts) >= 3:
        p3 = re.sub(r"^TEST", "", parts[2], flags=re.IGNORECASE)
        if p3.isdigit():
            rep = p3
    base = f"{PRINTS[0]}-B{exp_code}{d_int:02d}"
    return exp_code, d_int, (f"{base}-{int(rep):02d}" if rep else base)


# =============================================================================
# REDUCTION
# =============================================================================
def load_bearing_coupons(spec_df, t_col):
    """Reduce every bearing MTS file to its D953 properties and curve."""
    coupons = []
    for fp in sorted(MTS_DIR.glob(f"{PRINTS[0]}-B*.txt")):
        stem = fp.stem
        if stem.upper() in {s.upper() for s in SKIP_BEAR}:
            continue
        exp_code, d_int, spec_id = parse_bearing_stem(stem)
        if exp_code is None:
            continue
        t_mm = get_t_mm(spec_df, t_col, spec_id)
        if t_mm is None:
            print(f"[skip] {spec_id} — no thickness in {SPECIMEN_CSV.name}")
            continue

        raw = read_mts_txt(fp).dropna(subset=["disp_mm", "force_N"])
        if len(raw) < 10:
            continue
        d = raw["disp_mm"].to_numpy() - raw["disp_mm"].iloc[0]
        f = raw["force_N"].to_numpy()
        F_max = float(np.max(f))

        # Toe correction (D638 Annex A1.3, which D790 §12.1 and this reduction
        # both borrow): fit the linear region and project it back to F = 0.
        # That x-intercept is the machine take-up, and all deformation is
        # measured from it. The linear region is the 10-40 % F_max band — D953
        # defines no modulus and so no modulus window, so a load band is what
        # is left. mts_plots.toe_shift applies the identical band, so the two
        # scripts correct a bearing coupon the same way.
        #
        # LOADING BRANCH ONLY. A1.3 constructs the continuation of the loading
        # curve; the record runs well past peak, and those unloading points sit
        # at large displacement with the force back inside the 10-40 % band, so
        # including them tips the fitted line over and projects it to a NEGATIVE
        # intercept. That is what used to happen on the three BIS00 coupons —
        # they were being clamped to a zero toe and reading S_b = 17-28 MPa
        # against an S_max of ~150. Fitted on the loading branch alone their toe
        # is 0.28-0.39 mm, in family with every other coupon, and S_b lands at
        # ~50 MPa. The clamp below is kept as a guard, but it should no longer
        # have anything to catch.
        i_max = int(np.argmax(f))
        pre_peak = np.arange(len(f)) <= i_max
        fit_m = pre_peak & (f >= FIT_LO * F_max) & (f <= FIT_HI * F_max)
        if fit_m.sum() < 3:
            continue
        slope, intercept = np.polyfit(d[fit_m], f[fit_m], 1)
        d_zero = max(0.0, -intercept / slope)
        d_corr = d - d_zero

        # P at 4 % hole deformation (D953-19 §3.2.2): pre-peak only, made
        # monotone so np.interp is valid.
        dc, fc = d_corr[:i_max + 1], f[:i_max + 1]
        keep = np.concatenate(([True], np.diff(np.maximum.accumulate(dc)) > 0))
        dc_m, fc_m = dc[keep], fc[keep]
        failed = float(np.max(dc_m)) < DEF_4PCT_MM
        P_4pct = F_max if failed else float(np.interp(DEF_4PCT_MM, dc_m, fc_m))

        area = t_mm * HOLE_D_MM       # D953-19 §13.3
        coupons.append({
            "cid": spec_id, "exp": exp_code, "dir": d_int,
            "d_corr": d_corr, "force_N": f, "i_peak": i_max,
            "area": area, "d_zero": d_zero, "failed": failed,
            "F_max_N": F_max,
            "S_b":   P_4pct / area,   # D953-19 §13.3 Eq. 1
            "S_max": F_max / area,    # D953-19 §3.2.5
        })
        note = "  (failed <4%)" if failed else ""
        print(f"[{spec_id}]  S_b={P_4pct / area:.1f} MPa  "
              f"S_max={F_max / area:.1f} MPa  toe_shift={d_zero:.3f} mm{note}")

    print(f"\n{len(coupons)} bearing coupon(s) loaded\n")
    return coupons


# =============================================================================
# GROUP PLOTS
# =============================================================================
def active_exposures(coupons):
    return [e for e in EXPOSURE_ORDER if any(c["exp"] == e for c in coupons)]


def exposure_handles(coupons):
    return [mpatches.Patch(color=EXPOSURE_COLORS[e], label=EXPOSURE_LABELS[e])
            for e in active_exposures(coupons)]


def save(fig, name):
    out = FIGS_ROOT / name
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def plot_curves(coupons):
    """Load-deflection overlay, one panel per print direction."""
    fig, axes = plt.subplots(1, len(DIRECTIONS),
                             figsize=(5.5 * len(DIRECTIONS), 4.8), squeeze=False)
    for ax, d_int in zip(axes[0], DIRECTIONS):
        for c in coupons:
            if c["dir"] != d_int:
                continue
            col = EXPOSURE_COLORS[c["exp"]]
            ax.plot(c["d_corr"], c["force_N"] / 1000, color=col, lw=1.1, alpha=0.85)
            ax.scatter(c["d_corr"][c["i_peak"]], c["force_N"][c["i_peak"]] / 1000,
                       color=col, s=35, zorder=5)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Hole deformation (mm)")
        ax.set_ylabel("Bearing load (kN)")
        ax.set_title(f"Direction {d_int}°")
        ax.grid(alpha=0.25, ls="--")

    axes[0][-1].legend(handles=exposure_handles(coupons), fontsize=8, loc="best")
    fig.tight_layout()
    save(fig, "bearing_curves.png")


def plot_curves_single(coupons):
    """Bearing stress against deflection, every coupon on one axis."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for c in coupons:
        col = EXPOSURE_COLORS[c["exp"]]
        sl = slice(0, c["i_peak"] + 1)
        stress = c["force_N"][sl] / c["area"]
        ax.plot(c["d_corr"][sl], stress, color=col, lw=0.9, alpha=0.75)
        ax.scatter(c["d_corr"][c["i_peak"]], stress[-1], color=col,
                   marker=DIRECTION_MARKERS[c["dir"]], s=35, zorder=5)

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Displacement (mm)", fontsize=12)
    ax.set_ylabel("Bearing Stress (MPa)", fontsize=12)
    ax.set_title("P01 — Bearing: Stress vs. Displacement", fontsize=13)
    ax.grid(alpha=0.25, ls="--")

    handles = exposure_handles(coupons)
    for d_int in DIRECTIONS:
        handles.append(mlines.Line2D([], [], color="k", ls="None",
                                     marker=DIRECTION_MARKERS[d_int], markersize=6,
                                     label=f"{d_int}°"))
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.85)
    fig.tight_layout()
    save(fig, "bearing_curves_single.png")


def plot_summary(coupons):
    """S_b and S_max scatter, exposure x direction."""
    exp_active = active_exposures(coupons)
    panels = [("S_b", "S_b (MPa)", "Pin-Bearing Strength (4% def)"),
              ("S_max", "S_max (MPa)", "Max Bearing Stress")]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, panel in zip(axes, panels):
        key, ylabel, title = panel
        for c in coupons:
            ei = exp_active.index(c["exp"])
            x = c["dir"] + (ei - (len(exp_active) - 1) / 2) * 1.5
            ax.scatter(x, c[key], color=EXPOSURE_COLORS[c["exp"]],
                       marker="v" if c["failed"] else "o", s=55, zorder=5)
        x_half = (len(exp_active) - 1) / 2 * 1.5 + 0.8
        for d_int in DIRECTIONS:
            vals = [c[key] for c in coupons if c["dir"] == d_int]
            if vals:
                ax.hlines(np.mean(vals), d_int - x_half, d_int + x_half,
                          colors="k", linewidths=0.8, linestyles="--", zorder=4)
        ax.set_xticks(DIRECTIONS)
        ax.set_xticklabels([f"{d}°" for d in DIRECTIONS])
        ax.set_xlabel("Print Direction")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25, ls="--", axis="y")

    handles = exposure_handles(coupons)
    handles.append(plt.scatter([], [], marker="v", color="k", label="failed <4%"))
    axes[-1].legend(handles=handles, fontsize=8, loc="best")
    fig.suptitle("P01 Pin-Bearing — Property Summary (D953-19)", fontsize=12)
    fig.tight_layout()
    save(fig, "bearing_summary.png")


def plot_max_stress(coupons):
    """S_max by exposure, shaded dry / wet."""
    exp_active = active_exposures(coupons)
    x_pos = np.arange(len(exp_active))

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, e in enumerate(exp_active):
        # CL and UV are the dry conditions, SW and IS the wet ones.
        ax.axvspan(x_pos[i] - 0.5, x_pos[i] + 0.5,
                   color="#f7f5e8" if e in ("CL", "UV") else "#e8f7f5", alpha=0.35)

    for c in coupons:
        ax.scatter(x_pos[exp_active.index(c["exp"])], c["S_max"],
                   color=EXPOSURE_COLORS[c["exp"]],
                   marker=DIRECTION_MARKERS[c["dir"]],
                   s=70, edgecolor="black", linewidth=0.4, zorder=5)

    for d_int in DIRECTIONS:
        vals = [c["S_max"] for c in coupons if c["dir"] == d_int]
        if vals:
            ax.axhline(np.mean(vals), color="k", ls="--", lw=0.7, alpha=0.6, zorder=3)

    ax.set_xlim(-0.5, len(exp_active) - 0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([EXPOSURE_LABELS[e] for e in exp_active])
    ax.set_xlabel("Exposure Condition", fontsize=12)
    ax.set_ylabel("Max Bearing Stress (MPa)", fontsize=12)
    ax.set_title("P01 — Bearing: Max Stress", fontsize=13)
    ax.grid(alpha=0.25, ls="--", axis="y")

    handles = [mlines.Line2D([], [], color="black", marker=DIRECTION_MARKERS[d],
                             linestyle="None", markersize=7, label=f"{d}°")
               for d in DIRECTIONS]
    ax.legend(handles=handles, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.tight_layout()
    save(fig, "bearing_max_stress.png")


# =============================================================================
# STATS — mean +/- std (CV%)
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
    The Level 3s write Tensile and Flexural into this same file."""
    if path.exists():
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="a",
                                if_sheet_exists="replace")
    else:
        writer = pd.ExcelWriter(path, engine="openpyxl", mode="w")
    with writer:
        df.to_excel(writer, sheet_name=sheet, index=False)


def run_print_stats(coupons):
    #     (key,       header,        dec_m, dec_s, w_m, w_s)
    props = [
        ("F_max_N", "F_max (N)",   1, 1, 6, 5),
        ("S_b",     "S_b (MPa)",   1, 1, 5, 4),
        ("S_max",   "S_max (MPa)", 1, 1, 5, 4),
    ]
    print_table("PIN-BEARING PROPERTIES — ASTM D953-19 Procedure A",
                coupons, EXPOSURE_ORDER, DIRECTIONS, props, EXPOSURE_LABELS)

    n_failed = sum(c["failed"] for c in coupons)
    if n_failed:
        print(f"\n  Note: {n_failed} coupon(s) failed before reaching 4% hole deformation.")
        print("        For those, P_4pct = F_max (conservative upper bound for S_b).")

    write_stats_sheet(OUT_STATS_XLSX, "Bearing",
                      build_stats_df(coupons, EXPOSURE_ORDER, DIRECTIONS,
                                     props, EXPOSURE_LABELS))
    print(f"\nExported: {OUT_STATS_XLSX} [Bearing]")


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 70)
    print("Bearing_Level3 — D953 reduction, figures and statistics")
    print("=" * 70)
    FIGS_ROOT.mkdir(parents=True, exist_ok=True)

    spec_df = load_spec_sheet()
    t_col = next((c for c in spec_df.columns if "thickness" in c.lower()), None)
    coupons = load_bearing_coupons(spec_df, t_col)
    if not coupons:
        print("No bearing coupons found — nothing to do.")
        return

    if DO_GROUP_PLOTS:
        plot_curves(coupons)
        plot_curves_single(coupons)
        plot_summary(coupons)
        plot_max_stress(coupons)

    if DO_PRINT_STATS:
        run_print_stats(coupons)

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
