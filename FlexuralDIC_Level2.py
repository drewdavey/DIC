#!/usr/bin/env python3
"""
FlexuralDIC_Level2.py  —  FSR Flexural Coupons (ASTM D790, 3-point bend)
=========================================================================
Reads Level-1's per-coupon frame CSV, applies failure truncation, and computes
ASTM D790 flexural properties. No plotting here — Level 3 doesn't exist yet
(see the README); plot from the per-coupon CSV and flexural_properties.csv
until it does.

Standards compliance — what each calculation cites
  Flexural stress      : D790 §12.2 Eq.3   sigma = 3PL / (2 b d^2)
  Flexural strain      : D790 §12.3 Eq.4   eps   = 6 D d / L^2
  Tangent modulus      : D790 §12.4 Eq.5   E_B   = L^3 m / (4 b d^3)
  Chord modulus        : D790 §12.5
  Toe compensation     : D790 §12.1 / D638 Annex A1
  Flexural strength    : D790 §3.2.7 (max flexural stress, if it breaks < 5 %)
  Curvature modulus    : Euler-Bernoulli, E = (M/I) / kappa. Not in D790 — the
                         DIC-native equivalent of E_B, and the reference here.

THREE STRAIN CHANNELS, AND THEY ARE SUPPOSED TO DISAGREE
---------------------------------------------------------
All three are the strain of the same outer fibre, and they share one force and
one specimen, so any separation between them is measurement method, not
material:

  eps_curvature   kappa * d/2                 DIC only. No span, no fixture, no
                                              machine compliance. THE REFERENCE.
  eps_deflection  6 * D * d / L^2  (Eq.4)     DIC deflection + D790's
                                              small-deflection kinematics +
                                              the nominal span L.
  eps_crosshead   6 * D_mts * d / L^2         all of the above, plus machine
                                              compliance and support indentation.

eps_deflection - eps_curvature is the error in D790's small-deflection
kinematics, plus the fact that the moment actually crosses zero a little
inboard of the rollers (see Level 1); eps_crosshead - eps_deflection is
machine compliance and support indentation. Both gaps are printed.
eps_curvature is the reference because it is the only one with no span in it
at all -- it stays right whatever the fixture is doing.

WHAT THE ROI-CENTRING ASSUMPTION TOUCHES
-----------------------------------------
The correlation subset radius insets the ROI from both specimen faces, so
Level 1 has to *assume* the ROI is centred on the depth to place them. That
assumption sets where the faces are, so it affects eps_bot/eps_top and it fully
determines the reported neutral-axis offset and membrane strain: an ROI
mis-centred by delta reads as a neutral axis delta off mid-depth and a membrane
strain of exactly kappa*delta, indistinguishable from a genuinely asymmetric
specimen. It does NOT touch kappa, and therefore does not touch eps_curvature
or any modulus. That is precisely why the reference channel is built from
kappa rather than from eps_bot.

The two causes separate by their load dependence — a mis-centred ROI is a
constant present from the first frame, while asymmetric material nonlinearity
only develops under load — so the low-load baseline and the drift above it are
reported separately. Only the drift is a statement about the material.

PROCESSING NOTE
  Level 1 writes the full, untruncated per-frame record. Level 2 reads that
  SAME file and appends its own columns to it in place — no separate
  per-coupon output file, same convention as the tensile pipeline. Failure
  truncation marks rows out of the analysis window rather than dropping them:
  every L2-derived column is NaN outside the window and the boolean 'kept'
  column marks which rows are inside it. The window additionally requires DIC
  validity — a bend specimen loses correlation at fracture, several frames
  before the load channel finishes falling.

  No smoothing pass. D790 does not call for filtering the stress-strain
  record, the MTS force is the primary channel and is already clean (~0.7 % of
  span), and the DIC-derived channels are the cleanest signals in the test.

INPUT per coupon
  <DIC_DIR>/<coupon_id>.csv        Level-1's per-frame record
  <DIC_DIR>/flexural_geometry.csv  per-coupon b, d, fixture, tare, alignment

OUTPUT per coupon
  <DIC_DIR>/<coupon_id>.csv          Level-1's columns, unchanged, plus: kept,
                                     stress_MPa, M_over_I_MPa_per_mm,
                                     eps_curvature, eps_deflection,
                                     eps_crosshead (all toe-corrected). All NaN
                                     where kept is False.
  <DIC_DIR>/flexural_properties.csv  one row per coupon: the D790 scalars.
  <DIC_DIR>/flexural_group_stats.csv mean/std/count per exposure × direction.

  FSR-SpecimenTesting.xlsx is NOT written. Whether the D790 scalars belong in
  that workbook alongside the tensile ones, and under what column names, is
  still undecided — see the README. Until it is, the properties CSV above is
  the single source of truth for flexural scalars.
"""

from __future__ import annotations

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

# =============================================================================
# SWITCHES
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "IS": True}
DIRECTIONS = {"00": True, "90": True}
REPLICATES = ["01", "02", "03"]

# =============================================================================
# GEOMETRY
# Must match FlexuralDIC_Level1.FLEX_SPAN_MM. 8.00 in, confirmed against the
# fixture. D790 defines stress and strain on this nominal span, which is what
# is used here even though the moment crosses zero a little inboard of the
# rollers (see Level 1). eps_curvature has no span in it either way.
# =============================================================================
IN2MM        = 25.4
FLEX_SPAN_MM = 8.00 * IN2MM

# =============================================================================
# FAILURE TRUNCATION  — mirrors TensileDIC_Level2's, plus a DIC-validity term
# =============================================================================
LOAD_START_FRAC = 0.02   # pre-load: outside the window until force exceeds this × peak
LOAD_END_FRAC   = 0.50   # post-peak: cut at the first frame below this × peak

# =============================================================================
# PROPERTY SETTINGS
# =============================================================================
# Modulus fit window, in flexural strain. D790 §12.4: "initial straight-line
# portion". These coupons break near 2 % strain, so 0.05-0.25 % is well inside
# the linear region and clear of the toe.
MODULUS_STRAIN_RANGE = (0.0005, 0.0025)
# D790 §12.5 chord modulus endpoints.
CHORD_STRAIN_RANGE   = (0.0005, 0.0025)


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
    """Return (exposure, direction_str) e.g. ('CL', '00')."""
    part = cid.split("-")[1]
    return part[1:-2], part[-2:]

def find_frames_csv(cid):
    p = DIC_DIR / f"{cid}.csv"
    return p if p.exists() else None

def load_geometry() -> pd.DataFrame:
    """flexural_geometry.csv, indexed by coupon — b, d, fixture, tare and
    alignment, written once per coupon by Level 1."""
    fp = DIC_DIR / "flexural_geometry.csv"
    if not fp.exists():
        return pd.DataFrame()
    return pd.read_csv(fp).set_index("coupon")


def truncation_mask(force: np.ndarray, has_dic: np.ndarray) -> np.ndarray:
    """Valid analysis window: after the load picks up, up to the peak's decay,
    and only where DIC actually correlated.

    The DIC-validity term is the difference from the tensile version — a bend
    specimen loses correlation at fracture, several frames before the load
    channel finishes falling.
    """
    n = len(force)
    f = np.where(has_dic, force, np.nan)
    peak = float(np.nanmax(np.abs(f))) if np.isfinite(f).any() else 0.0
    if peak <= 0:
        return np.zeros(n, dtype=bool)
    i_pk = int(np.nanargmax(np.abs(f)))
    # Start on the *rising edge*, not the first sample over threshold: with the
    # tare removed the pre-touchdown baseline sits on zero with noise either
    # side of it, so "first sample above 2 % of peak" can trigger on noise tens
    # of seconds early. Take the last sample below threshold before the peak.
    below = np.flatnonzero(np.abs(np.nan_to_num(f[:i_pk + 1])) < LOAD_START_FRAC * peak)
    i0 = int(below[-1]) + 1 if below.size else 0
    post = np.flatnonzero(np.abs(np.nan_to_num(f[i_pk:], nan=0.0)) < LOAD_END_FRAC * peak)
    i1 = int(i_pk + post[0]) - 1 if post.size else n - 1
    mask = np.zeros(n, dtype=bool)
    mask[i0:max(i0, i1) + 1] = True
    return mask & has_dic & np.isfinite(force)


def fit_modulus(eps: np.ndarray, sig: np.ndarray,
                window: tuple[float, float]) -> tuple[float, float, int]:
    """Least-squares slope of sigma against eps inside `window`.

    Returns (slope MPa, toe offset in strain, n points). The toe offset is the
    fit line's x-intercept — D790 §12.1 / D638 Annex A1: extend the straight
    line back to zero stress and take that as the true strain origin.

    Fitted twice. The window is specified in *corrected* strain, but the
    correction is what the first fit produces, so the second pass re-selects the
    window after shifting the origin. Without it the tangent and the chord
    modulus below are quietly measured over two different strain ranges — offset
    by the toe — and are not comparable.
    """
    toe = 0.0
    slope = np.nan
    n = 0
    for _ in range(2):
        e = eps - toe
        m = (e >= window[0]) & (e <= window[1]) & np.isfinite(e) & np.isfinite(sig)
        n = int(m.sum())
        if n < 3:
            return np.nan, np.nan, n
        slope, icept = np.polyfit(e[m], sig[m], 1)
        if slope == 0:
            return np.nan, np.nan, n
        toe = toe + (-icept / slope)      # x-intercept, accumulated
    return float(slope), float(toe), n


def chord_modulus(eps: np.ndarray, sig: np.ndarray,
                  window: tuple[float, float]) -> float:
    """D790 §12.5 chord modulus: secant between two strain points."""
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
def compute_properties(cid: str, frames: pd.DataFrame,
                       geom: pd.Series) -> tuple[pd.DataFrame, dict] | None:
    """D790 property extraction for one coupon's truncated record.

    Appends the L2 columns to `frames` in place and returns (frames, props).
    """
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

    nan = lambda: np.full(len(frames), np.nan)

    # --- D790 §12.2 Eq.3: sigma = 3 P L / (2 b d^2) ---
    sigma = nan()
    sigma[ki] = 3.0 * force[ki] * L / (2.0 * b * d ** 2)
    # --- bending moment per unit second moment, for the curvature modulus ---
    m_over_i = nan()
    m_over_i[ki] = (force[ki] * L / 4.0) / I_mm4

    # --- the three strain channels (see module docstring) ---
    eps_curv = nan()
    eps_curv[ki] = frames["kappa_1pmm"].to_numpy()[ki] * (d / 2.0)
    eps_defl = nan()
    eps_defl[ki] = 6.0 * frames["defl_mm"].to_numpy()[ki] * d / L ** 2
    eps_xhead = nan()
    dx = frames["disp_mts_mm"].to_numpy()[ki]
    dx = dx - np.nanmin(dx) if np.isfinite(dx).any() else dx
    eps_xhead[ki] = 6.0 * dx * d / L ** 2

    props: dict = {"coupon": cid, "b_mm": b, "d_mm": d, "span_mm": L,
                   "I_mm4": I_mm4, "span_to_depth": L / d,
                   "n_frames": len(frames), "n_kept": int(ki.size),
                   "tare_mts_N": float(geom.get("tare_mts_N", np.nan)),
                   "P_max_N": float(np.nanmax(force[ki]))}

    # --- modulus per strain channel, each toe-compensated on its own ---
    for name, eps in (("curvature", eps_curv), ("deflection", eps_defl),
                      ("crosshead", eps_xhead)):
        slope, toe, n_fit = fit_modulus(eps[ki], sigma[ki], MODULUS_STRAIN_RANGE)
        eps_corr = eps.copy()
        if np.isfinite(toe):
            eps_corr[ki] = eps[ki] - toe          # D790 §12.1 toe compensation
        props[f"Ef_{name}_GPa"] = slope / 1000.0 if np.isfinite(slope) else np.nan
        props[f"Ef_{name}_chord_GPa"] = chord_modulus(eps_corr[ki], sigma[ki],
                                                      CHORD_STRAIN_RANGE) / 1000.0
        props[f"toe_{name}"] = toe
        props[f"n_fit_{name}"] = n_fit
        frames[f"eps_{name}"] = eps_corr
        i_max = int(np.nanargmax(sigma[ki]))
        props[f"eps_at_max_{name}"] = float(eps_corr[ki][i_max])

    # --- curvature modulus: E = (M/I) / kappa, no span, no compliance -------
    kap = frames["kappa_1pmm"].to_numpy()
    fw = np.flatnonzero((eps_curv >= MODULUS_STRAIN_RANGE[0])
                        & (eps_curv <= MODULUS_STRAIN_RANGE[1]) & kept)
    props["Ef_MI_kappa_GPa"] = (float(np.polyfit(kap[fw], m_over_i[fw], 1)[0]) / 1000.0
                                if fw.size >= 3 else np.nan)
    props["n_fit_MI_kappa"] = int(fw.size)

    # --- strength (D790 §3.2.7: flexural strength = max flexural stress) ----
    i_max = int(np.nanargmax(sigma[ki]))
    props["sigma_fM_MPa"] = float(np.nanmax(sigma[ki]))
    props["sigma_fB_MPa"] = float(sigma[ki][-1])       # last correlated frame = break
    props["broke_before_5pct"] = bool(props["eps_at_max_curvature"] < 0.05)

    # --- bending-quality diagnostics ---
    # R^2 of the through-depth fit is meaningless at low load (the strain range
    # is below the DIC noise floor there), so it is summarised only over the
    # upper part of the ramp where the measurement has something to measure.
    y_c = 0.5 * (float(geom["y_bot_mm"]) + float(geom["y_top_mm"]))
    loaded = ki[force[ki] >= 0.25 * props["P_max_N"]]
    props["profile_r2_median"] = float(
        np.nanmedian(frames["profile_r2"].to_numpy()[loaded]))
    props["profile_r2_at_max"] = float(frames["profile_r2"].to_numpy()[ki][i_max])
    na = frames["na_Y_mm"].to_numpy() - y_c
    props["na_offset_at_max_mm"] = float(na[ki][i_max])
    props["na_offset_at_low_load_mm"] = float(np.nanmedian(
        na[ki[force[ki] <= 0.25 * props["P_max_N"]]]))
    props["na_drift_mm"] = (props["na_offset_at_max_mm"]
                            - props["na_offset_at_low_load_mm"])
    props["eps_membrane_at_max"] = float(frames["eps_membrane"].to_numpy()[ki][i_max])
    with np.errstate(invalid="ignore"):
        props["defl_over_span_at_max"] = float(
            frames["defl_mm"].to_numpy()[ki][i_max] / L)
    for k in ("frame_rate_hz", "mts_offset_s", "break_frame", "disp_check_rmse_mm"):
        props[k] = float(geom.get(k, np.nan))

    frames["stress_MPa"] = sigma
    frames["M_over_I_MPa_per_mm"] = m_over_i
    return frames, props


def print_coupon(props: dict, ki0: int, ki1: int) -> None:
    p = props
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
def main() -> None:
    t0 = time.time()
    print("=" * 74)
    print("FlexuralDIC_Level2 — truncate and compute D790 properties")
    print("=" * 74)
    print(f"DIC dir : {DIC_DIR}")
    print(f"Span    : {FLEX_SPAN_MM:.1f} mm ({FLEX_SPAN_MM / IN2MM:.2f} in) — "
          f"confirmed fixture setting, D790 nominal")
    print()

    geometry = load_geometry()
    if geometry.empty:
        print("[!] no flexural_geometry.csv — run FlexuralDIC_Level1 first")
        return

    rows = []
    for cid in selected_coupons():
        frames_fp = find_frames_csv(cid)
        if frames_fp is None:
            print(f"[{cid}] no per-coupon CSV — run Level 1 first")
            continue
        if cid not in geometry.index:
            print(f"[{cid}] not in flexural_geometry.csv — re-run Level 1 for it")
            continue

        frames = pd.read_csv(frames_fp)
        result = compute_properties(cid, frames, geometry.loc[cid])
        if result is None:
            continue
        frames, props = result

        ki = np.flatnonzero(frames["kept"].to_numpy(dtype=bool))
        print_coupon(props, int(ki[0]), int(ki[-1]))

        frames.to_csv(frames_fp, index=False, float_format="%.6g")
        rows.append(props)

    if rows:
        prop_fp = DIC_DIR / "flexural_properties.csv"
        pd.DataFrame(rows).to_csv(prop_fp, index=False, float_format="%.6g")

        # ---- mean & std per (exposure, direction), as D638 §11.7 does for
        # tensile. D790 §12.9 asks for the same summary statistics.
        df_sum = pd.DataFrame(rows)
        df_sum["exposure"] = df_sum["coupon"].map(lambda c: parse_id(c)[0])
        df_sum["direction"] = df_sum["coupon"].map(lambda c: parse_id(c)[1])
        agg_cols = ["sigma_fM_MPa", "Ef_curvature_GPa", "Ef_deflection_GPa",
                    "Ef_crosshead_GPa", "eps_at_max_curvature", "P_max_N"]
        stats_fp = DIC_DIR / "flexural_group_stats.csv"
        (df_sum.groupby(["exposure", "direction"])[agg_cols]
               .agg(["mean", "std", "count"]).to_csv(stats_fp))

        print(f"\n{len(rows)} coupon(s) → DIC/*.csv, {prop_fp.name}, {stats_fp.name}")

    print(f"\nDone. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
