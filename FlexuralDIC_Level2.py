#!/usr/bin/env python3
"""
FlexuralDIC_Level2.py  —  FSR Flexural Coupons (ASTM D790, 3-point bend)
=========================================================================
Reads Level-1's per-frame CSV, truncates it to the valid test window, and
computes ASTM D790 properties. No plotting — that is Level 3's job.

Three strain channels are computed for the same outer fibre, so any gap
between them is measurement method, not material:
  eps_curvature   kappa * d/2         DIC only, no span. THE REFERENCE.
  eps_deflection  6 D d / L^2         D790 Eq.4 from the DIC midspan deflection.
  eps_crosshead   6 D_mts d / L^2     the same, from crosshead travel.
See README.md for what the gaps mean and for the ROI-centring caveat.

Cited: stress D790 §12.2 Eq.3, strain §12.3 Eq.4, tangent modulus §12.4 Eq.5,
chord modulus §12.5, toe compensation §12.1, strength §3.2.7. The curvature
modulus E = (M/I)/kappa is Euler-Bernoulli, not D790.

INPUT
  <DIC_DIR>/<coupon_id>.csv       Level-1's per-frame record
  <DIC_DIR>/coupon_scalars.csv    per-coupon b, d, fixture, tare, alignment

OUTPUT
  <DIC_DIR>/<coupon_id>.csv          Level-1's columns plus: kept, stress_MPa,
                                     M_over_I_MPa_per_mm, eps_curvature,
                                     eps_deflection, eps_crosshead. All NaN
                                     where kept is False.
  FSR-SpecimenTesting.csv            the D790 scalars, under the "Flex ..."
                                     headers in SPECIMEN_SHEET_COLUMNS below.
  <DIC_DIR>/level2_group_stats.csv   mean/std/count per exposure x direction,
                                     under a 'test' index level of "flexural".
                                     Shared with TensileDIC_Level2; each script
                                     replaces only its own test's rows.
"""

from __future__ import annotations
import os
import sys
import time
from pathlib import Path
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
# The CSV is the specimen sheet — there is no .xlsx any more. See README.md.
SPECIMEN_CSV = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting.csv"
)
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "IS": True}
DIRECTIONS = {"00": True, "90": True}
REPLICATES = ["01", "02", "03"]

# =============================================================================
# ANALYSIS  — must match FlexuralDIC_Level1 and FlexuralDIC_Level3
# =============================================================================
IN2MM        = 25.4
FLEX_SPAN_MM = 8.00 * IN2MM     # L, support span — confirmed against the fixture

# Analysis window, as a fraction of peak load. Same values as TensileDIC_Level2.
LOAD_START_FRAC = 0.02   # outside the window until the load rises past this
LOAD_END_FRAC   = 0.50   # cut at the first post-peak frame below this

# Modulus fit windows, in toe-corrected strain. D790 §12.4 "initial
# straight-line portion" / §12.5 chord endpoints. Same values as
# TensileDIC_Level2 so the two pipelines measure a modulus the same way.
MODULUS_STRAIN_RANGE = (0.0005, 0.003)
CHORD_STRAIN_RANGE   = (0.001, 0.003)

# =============================================================================
# Scalar columns written into SPECIMEN_CSV, keyed by "Specimen ID".
# Every header is prefixed "Flex " so nothing here reads as a tensile column.
# Level-1 measurements (b, d, span, tare, alignment) are NOT repeated here —
# they live in coupon_scalars.csv.
# =============================================================================
SPECIMEN_SHEET_COLUMNS = {
    # headline D790 results
    "sigma_fM_MPa":             "Flex Strength (MPa)",
    "sigma_fB_MPa":             "Flex Stress at Break (MPa)",
    "broke_before_5pct":        "Flex Broke Before 5% Strain",
    "P_max_N":                  "Flex Peak Load (N)",
    # modulus, per strain channel (curvature is the reference)
    "Ef_curvature_GPa":         "Flex Ef curvature (GPa)",
    "Ef_curvature_chord_GPa":   "Flex Ef curvature chord (GPa)",
    "eps_at_max_curvature":     "Flex Strain at Max (curvature)",
    "toe_curvature":            "Flex Toe Strain (curvature)",
    "Ef_deflection_GPa":        "Flex Ef deflection (GPa)",
    "Ef_deflection_chord_GPa":  "Flex Ef deflection chord (GPa)",
    "eps_at_max_deflection":    "Flex Strain at Max (deflection)",
    "toe_deflection":           "Flex Toe Strain (deflection)",
    "Ef_crosshead_GPa":         "Flex Ef crosshead (GPa)",
    "Ef_crosshead_chord_GPa":   "Flex Ef crosshead chord (GPa)",
    "eps_at_max_crosshead":     "Flex Strain at Max (crosshead)",
    "toe_crosshead":            "Flex Toe Strain (crosshead)",
    "Ef_MI_kappa_GPa":          "Flex Ef M/I-kappa (GPa)",
    # quality diagnostics: whether the bending assumptions held
    "profile_r2_median":        "Flex Profile R2 (median)",
    "profile_r2_at_max":        "Flex Profile R2 (at max)",
    "na_offset_at_max_mm":      "Flex NA Offset at Max (mm)",
    "na_offset_at_low_load_mm": "Flex NA Offset at Low Load (mm)",
    "na_drift_mm":              "Flex NA Drift (mm)",
    "eps_membrane_at_max":      "Flex Membrane Strain at Max",
    "defl_over_span_at_max":    "Flex Deflection/Span at Max",
    "n_kept":                   "Flex Frames Kept",
    "n_fit_curvature":          "Flex Modulus Fit Points (curvature)",
    "n_fit_deflection":         "Flex Modulus Fit Points (deflection)",
    "n_fit_crosshead":          "Flex Modulus Fit Points (crosshead)",
    "n_fit_MI_kappa":           "Flex Modulus Fit Points (M/I-kappa)",
}


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

def load_coupon_scalars():
    """coupon_scalars.csv indexed by coupon — b, d, fixture, tare, alignment.
    Shared with the tensile pipeline; only flexural IDs are looked up here."""
    fp = DIC_DIR / "coupon_scalars.csv"
    if not fp.exists():
        return pd.DataFrame()
    return pd.read_csv(fp).set_index("coupon")

def read_specimen_csv():
    """The specimen sheet as raw text. dtype=str with keep_default_na=False is
    what makes it safe to write back: cells this script does not touch
    round-trip character for character."""
    for enc in CSV_ENCODINGS:
        try:
            return pd.read_csv(SPECIMEN_CSV, encoding=enc,
                               dtype=str, keep_default_na=False)
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f"{SPECIMEN_CSV.name}: not decodable as "
                       f"{'/'.join(CSV_ENCODINGS)}")

def cell_text(v):
    """One specimen-sheet cell. Missing or non-finite writes blank, so a
    property that could not be computed clears the cell instead of leaving the
    previous run's value standing."""
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"      # what Excel reads back as boolean
    if v is None:
        return ""
    v = float(v)
    return "" if not np.isfinite(v) else f"{v:.12g}"

def write_specimen_sheet(rows):
    """Write each coupon's scalars into its row of SPECIMEN_CSV, matched by
    Specimen ID. Missing columns are added at the end; every other cell is
    written back as it was read. Skipped with a warning if the file can't be
    read or replaced — e.g. open in Excel."""
    try:
        df = read_specimen_csv()
    except FileNotFoundError:
        print(f"[!] {SPECIMEN_CSV} not found — skipping specimen sheet update")
        return
    except Exception as exc:
        print(f"[!] {SPECIMEN_CSV.name}: {exc} — skipping specimen sheet update")
        return

    if "Specimen ID" not in df.columns:
        print("[!] no 'Specimen ID' column in the specimen sheet — skipping update")
        return

    for label in SPECIMEN_SHEET_COLUMNS.values():
        if label not in df.columns:
            df[label] = ""
    col_pos = {c: j for j, c in enumerate(df.columns)}
    row_by_id = {cid: i for i, cid in enumerate(df["Specimen ID"])}

    n_written = 0
    for row in rows:
        i = row_by_id.get(row["coupon"])
        if i is None:
            print(f"[!] {row['coupon']} has no row in the specimen sheet — "
                  f"its properties were not written")
            continue
        for key, label in SPECIMEN_SHEET_COLUMNS.items():
            df.iat[i, col_pos[label]] = cell_text(row.get(key))
        n_written += 1

    # Write alongside and rename over: the sheet is the only copy of the
    # hand-entered geometry, so a half-written file would be real data loss.
    tmp = SPECIMEN_CSV.with_name(SPECIMEN_CSV.name + ".tmp")
    try:
        df.to_csv(tmp, index=False, encoding="utf-8-sig")
        os.replace(tmp, SPECIMEN_CSV)
        print(f"Specimen sheet: {n_written} coupon(s) → {SPECIMEN_CSV.name}")
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        print(f"[!] could not update {SPECIMEN_CSV.name} ({exc}) — "
              f"properties were not saved")


def read_existing_group_stats(fp):
    """Read level2_group_stats.csv, upgrading the pre-'test' two-level layout.

    The index depth is READ, not assumed: pandas writes a named MultiIndex as a
    third header line, and reading a two-level file with index_col=[0,1,2] does
    not raise — it swallows the first data column into the index and shifts
    every value one place left. Returns None if the file can't be understood.
    """
    try:
        with open(fp, encoding="utf-8") as fh:
            head = [fh.readline() for _ in range(3)]
        names = [c.strip() for c in head[2].strip("\r\n").split(",") if c.strip()]
        if not names:
            return None
        old = pd.read_csv(fp, header=[0, 1], index_col=list(range(len(names))))
        if "test" not in names:
            old = pd.concat({"tensile": old}, names=["test"])   # pre-'test' file
        # "00" round-trips through read_csv as the integer 0, so normalise.
        old.index = pd.MultiIndex.from_tuples(
            [(str(t), str(e), f"{int(d):02d}" if str(d).strip().isdigit() else str(d))
             for t, e, d in old.index],
            names=["test", "exposure", "direction"])
        return old
    except Exception:
        return None


def write_group_stats(rows):
    """Upsert this run's group statistics into level2_group_stats.csv.
    Only rows whose test == "flexural" are replaced."""
    fp = DIC_DIR / "level2_group_stats.csv"
    df_sum = pd.DataFrame(rows)
    df_sum["test"] = "flexural"
    df_sum["exposure"] = df_sum["coupon"].map(lambda c: parse_id(c)[0])
    df_sum["direction"] = df_sum["coupon"].map(lambda c: parse_id(c)[1])
    agg_cols = ["sigma_fM_MPa", "Ef_curvature_GPa", "Ef_deflection_GPa",
                "Ef_crosshead_GPa", "eps_at_max_curvature", "P_max_N"]
    group = (df_sum.groupby(["test", "exposure", "direction"])[agg_cols]
                   .agg(["mean", "std", "count"]))

    if fp.exists():
        old = read_existing_group_stats(fp)
        if old is None:
            print(f"[!] could not read existing {fp.name} — it is being replaced "
                  f"with flexural rows only.\n    Re-run TensileDIC_Level2.py to "
                  f"put the tensile rows back.")
        else:
            old = old.drop(index="flexural", level=0, errors="ignore")
            group = pd.concat([old, group]).sort_index()
    group.to_csv(fp)
    return fp


def truncation_mask(force, has_dic):
    """Valid analysis window: after the load picks up, up to the peak's decay,
    and only where DIC actually correlated. The DIC term is the difference from
    the tensile version — a bend specimen decorrelates at fracture, several
    frames before the load channel finishes falling."""
    n = len(force)
    f = np.where(has_dic, force, np.nan)
    peak = float(np.nanmax(np.abs(f))) if np.isfinite(f).any() else 0.0
    if peak <= 0:
        return np.zeros(n, dtype=bool)
    i_pk = int(np.nanargmax(np.abs(f)))
    # Start on the rising edge, not the first sample over threshold: with the
    # tare removed the pre-touchdown baseline sits on zero with noise either
    # side, so a plain threshold triggers tens of seconds early.
    below = np.flatnonzero(np.abs(np.nan_to_num(f[:i_pk + 1])) < LOAD_START_FRAC * peak)
    i0 = int(below[-1]) + 1 if below.size else 0
    post = np.flatnonzero(np.abs(np.nan_to_num(f[i_pk:], nan=0.0)) < LOAD_END_FRAC * peak)
    i1 = int(i_pk + post[0]) - 1 if post.size else n - 1
    mask = np.zeros(n, dtype=bool)
    mask[i0:max(i0, i1) + 1] = True
    return mask & has_dic & np.isfinite(force)


def fit_modulus(eps, sig, window):
    """Least-squares slope of sig against eps inside `window`.
    Returns (slope MPa, toe offset in strain, n points).

    The toe offset is the fit line's x-intercept — D790 §12.1 / D638 Annex A1.
    Fitted twice because the window is specified in CORRECTED strain but the
    correction is what the first fit produces; without the second pass the
    tangent and the chord are measured over two different strain ranges.
    Same routine as TensileDIC_Level2.fit_modulus.
    """
    toe = 0.0
    best = (np.nan, np.nan, 0)
    for step in range(2):
        e = eps - toe
        m = (e >= window[0]) & (e <= window[1]) & np.isfinite(e) & np.isfinite(sig)
        n = int(m.sum())
        if n < 3 or not np.any(m):
            # The second pass is a refinement, so losing it must not throw away
            # the first pass's answer. It only happens on a coupon whose toe is
            # comparable to the window itself (P01-TCL45-01).
            if step == 1:
                print(f"    [!] modulus re-fit after toe correction left {n} "
                      f"point(s) in {window} — keeping the first-pass fit")
            break
        slope, icept = np.polyfit(e[m], sig[m], 1)
        if slope == 0:
            break
        toe = toe + (-icept / slope)      # x-intercept, accumulated
        best = (float(slope), float(toe), n)
    return best


def chord_modulus(eps, sig, window):
    """D790 §12.5 chord modulus: the secant between the two ends of `window`,
    in already-toe-corrected strain. Returns NaN rather than extrapolating when
    the record does not reach both endpoints."""
    m = np.isfinite(eps) & np.isfinite(sig)
    if m.sum() < 3:
        return np.nan
    order = np.argsort(eps[m])
    e, s = eps[m][order], sig[m][order]
    lo, hi = window
    if not (e[0] <= lo and hi <= e[-1]):
        return np.nan
    return float((np.interp(hi, e, s) - np.interp(lo, e, s)) / (hi - lo))


# =============================================================================
# COMPUTE PROPERTIES
# =============================================================================
def compute_properties(cid, frames, geom):
    """D790 property extraction for one coupon. Appends the Level-2 columns to
    `frames` in place and returns (frames, props), or None if too few frames
    survive truncation."""
    b, d, L = float(geom["b_mm"]), float(geom["d_mm"]), FLEX_SPAN_MM
    I_mm4 = b * d ** 3 / 12.0

    force = frames["force_N"].to_numpy(dtype=float)
    has_dic = ((frames["n_pts"].to_numpy() > 0)
               & np.isfinite(frames["kappa_1pmm"].to_numpy()))
    kept = truncation_mask(force, has_dic)
    frames["kept"] = kept
    ki = np.flatnonzero(kept)
    if ki.size < 10:
        print(f"[{cid}]  truncation left too few frames ({ki.size})")
        return None

    n_rows = len(frames)
    sigma = np.full(n_rows, np.nan)
    sigma[ki] = 3.0 * force[ki] * L / (2.0 * b * d ** 2)        # §12.2 Eq.3
    m_over_i = np.full(n_rows, np.nan)
    m_over_i[ki] = (force[ki] * L / 4.0) / I_mm4                # for E = (M/I)/kappa

    eps_curv = np.full(n_rows, np.nan)
    eps_curv[ki] = frames["kappa_1pmm"].to_numpy()[ki] * (d / 2.0)
    eps_defl = np.full(n_rows, np.nan)
    eps_defl[ki] = 6.0 * frames["defl_mm"].to_numpy()[ki] * d / L ** 2   # §12.3 Eq.4
    eps_xhead = np.full(n_rows, np.nan)
    dx = frames["disp_mts_mm"].to_numpy()[ki]
    dx = dx - np.nanmin(dx) if np.isfinite(dx).any() else dx
    eps_xhead[ki] = 6.0 * dx * d / L ** 2

    props = {"coupon": cid, "b_mm": b, "d_mm": d, "span_mm": L,
             "I_mm4": I_mm4, "span_to_depth": L / d,
             "n_frames": n_rows, "n_kept": int(ki.size),
             "tare_mts_N": float(geom.get("tare_mts_N", np.nan)),
             "P_max_N": float(np.nanmax(force[ki]))}

    i_max = int(np.nanargmax(sigma[ki]))

    # Modulus per strain channel, each toe-compensated on its own.
    for name, eps in (("curvature", eps_curv), ("deflection", eps_defl),
                      ("crosshead", eps_xhead)):
        slope, toe, n_fit = fit_modulus(eps[ki], sigma[ki], MODULUS_STRAIN_RANGE)
        eps_corr = eps.copy()
        if np.isfinite(toe):
            eps_corr[ki] = eps[ki] - toe                        # §12.1
        props[f"Ef_{name}_GPa"] = slope / 1000.0 if np.isfinite(slope) else np.nan
        props[f"Ef_{name}_chord_GPa"] = chord_modulus(eps_corr[ki], sigma[ki],
                                                      CHORD_STRAIN_RANGE) / 1000.0
        props[f"toe_{name}"] = toe
        props[f"n_fit_{name}"] = n_fit
        props[f"eps_at_max_{name}"] = float(eps_corr[ki][i_max])
        frames[f"eps_{name}"] = eps_corr

    # Curvature modulus: E = (M/I) / kappa. No span, no machine compliance.
    kap = frames["kappa_1pmm"].to_numpy()
    fw = np.flatnonzero((eps_curv >= MODULUS_STRAIN_RANGE[0])
                        & (eps_curv <= MODULUS_STRAIN_RANGE[1]) & kept)
    props["Ef_MI_kappa_GPa"] = (float(np.polyfit(kap[fw], m_over_i[fw], 1)[0]) / 1000.0
                                if fw.size >= 3 else np.nan)
    props["n_fit_MI_kappa"] = int(fw.size)

    # Strength (§3.2.7: flexural strength = maximum flexural stress).
    props["sigma_fM_MPa"] = float(np.nanmax(sigma[ki]))
    props["sigma_fB_MPa"] = float(sigma[ki][-1])       # last correlated frame = break
    props["broke_before_5pct"] = bool(props["eps_at_max_curvature"] < 0.05)

    # Bending-quality diagnostics. The through-depth R^2 is meaningless at low
    # load (the strain range is below the DIC noise floor), so it is summarised
    # only above 25 % of peak.
    y_c = 0.5 * (float(geom["y_bot_mm"]) + float(geom["y_top_mm"]))
    loaded = ki[force[ki] >= 0.25 * props["P_max_N"]]
    low = ki[force[ki] <= 0.25 * props["P_max_N"]]
    r2 = frames["profile_r2"].to_numpy()
    na = frames["na_Y_mm"].to_numpy() - y_c
    props["profile_r2_median"] = float(np.nanmedian(r2[loaded]))
    props["profile_r2_at_max"] = float(r2[ki][i_max])
    props["na_offset_at_max_mm"] = float(na[ki][i_max])
    props["na_offset_at_low_load_mm"] = float(np.nanmedian(na[low]))
    props["na_drift_mm"] = (props["na_offset_at_max_mm"]
                            - props["na_offset_at_low_load_mm"])
    props["eps_membrane_at_max"] = float(frames["eps_membrane"].to_numpy()[ki][i_max])
    with np.errstate(invalid="ignore"):
        props["defl_over_span_at_max"] = float(frames["defl_mm"].to_numpy()[ki][i_max] / L)
    for k in ("frame_rate_hz", "mts_offset_s", "break_frame", "disp_check_rmse_mm"):
        props[k] = float(geom.get(k, np.nan))

    frames["stress_MPa"] = sigma
    frames["M_over_I_MPa_per_mm"] = m_over_i
    return frames, props


def print_coupon(p, ki0, ki1):
    print(f"[{p['coupon']}]  window frames {ki0}..{ki1} "
          f"({p['n_kept']} of {p['n_frames']} kept)")
    print(f"    sigma_fM = {p['sigma_fM_MPa']:.1f} MPa at "
          f"{p['eps_at_max_curvature'] * 100:.2f} % strain "
          f"(peak {p['P_max_N']:.0f} N net, tare {p['tare_mts_N']:.0f} N removed)")
    print(f"    E_f      = {p['Ef_curvature_GPa']:.2f} GPa curvature (reference) / "
          f"{p['Ef_deflection_GPa']:.2f} deflection / "
          f"{p['Ef_crosshead_GPa']:.2f} crosshead / "
          f"{p['Ef_MI_kappa_GPa']:.2f} M-I-vs-kappa")
    ec = p["Ef_curvature_GPa"]
    if np.isfinite(ec) and ec != 0:
        print(f"               deflection {100 * (p['Ef_deflection_GPa'] - ec) / ec:+.1f} %, "
              f"crosshead {100 * (p['Ef_crosshead_GPa'] - ec) / ec:+.1f} % "
              f"against the curvature value")
    print(f"    bending  R^2 {p['profile_r2_median']:.4f} median above 25 % load, "
          f"neutral axis {p['na_offset_at_low_load_mm']:+.3f} mm at low load "
          f"drifting {p['na_drift_mm']:+.3f} mm to peak")
    if not p["broke_before_5pct"]:
        print(f"    [!] did NOT break before 5 % strain — D790 §12.2 asks for the "
              f"stress at 5 % strain, not the maximum")
    if p["defl_over_span_at_max"] > 0.10:
        print(f"    [!] D/L = {p['defl_over_span_at_max']:.3f} at peak — above 0.10 "
              f"D790 §12.3 calls for the large-deflection correction")


# =============================================================================
# MAIN
# =============================================================================
def main():
    t0 = time.time()
    print("=" * 74)
    print("FlexuralDIC_Level2 — truncate and compute D790 properties")
    print("=" * 74)
    print(f"DIC dir : {DIC_DIR}")
    print(f"Span    : {FLEX_SPAN_MM:.1f} mm ({FLEX_SPAN_MM / IN2MM:.2f} in)")
    print()

    scalars = load_coupon_scalars()
    if scalars.empty:
        print("[!] no coupon_scalars.csv — run FlexuralDIC_Level1 first")
        return

    rows = []
    for cid in selected_coupons():
        frames_fp = DIC_DIR / f"{cid}.csv"
        if not frames_fp.exists():
            print(f"[{cid}] no per-coupon CSV — run Level 1 first")
            continue
        if cid not in scalars.index:
            print(f"[{cid}] not in coupon_scalars.csv — re-run Level 1 for it")
            continue

        frames = pd.read_csv(frames_fp)
        result = compute_properties(cid, frames, scalars.loc[cid])
        if result is None:
            continue
        frames, props = result

        ki = np.flatnonzero(frames["kept"].to_numpy(dtype=bool))
        print_coupon(props, int(ki[0]), int(ki[-1]))

        frames.to_csv(frames_fp, index=False, float_format="%.6g")
        rows.append(props)

    if rows:
        write_specimen_sheet(rows)
        stats_fp = write_group_stats(rows)      # D790 §12.9 summary statistics
        print(f"\n{len(rows)} coupon(s) → DIC/*.csv, {SPECIMEN_CSV.name}, "
              f"DIC/{stats_fp.name}")

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
