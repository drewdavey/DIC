#!/usr/bin/env python
"""
mts_plots.py — MTS-only curves for tensile, flexural, and pin-bearing coupons.

One figure per test type, two panels each:
  (left)  force vs. crosshead displacement  — the machine channels
  (right) stress vs. strain                 — the ASTM reduction of the same
                                               record, using the coupon
                                               geometry from the specimen sheet

Colour    = exposure group (CL / UV / SW / IS)
Linestyle = print orientation (0 / 45 / 90 deg)

Standards:
  Tensile      ASTM D638-22   (Tensile Properties of Plastics)
  Flexural     ASTM D790-17   (Flexural Properties, Procedure A, 3-point)
  Pin-bearing  ASTM D953-19   (Bearing Strength, Procedure A)

Every formula this script applies is written out in full, with its ASTM
source, in the ASTM EQUATIONS block below and again at the point of use.
The code is deliberately written to be read next to the standards rather
than to be fast.

Outputs (per print folder):
  figs/mts_tensile.png
  figs/mts_flexural.png
  figs/mts_bearing.png

This is a first-look script — no filtering, no property extraction. For the
D638 tensile properties use TensileDIC_Level1-3, for the D790 flexural
properties use FlexuralDIC_Level1-2, and for the D953 bearing strengths use
bearing_group_plots.py.

Usage:
  python mts_plots.py                      # P01, all three tests
  python mts_plots.py --print P02
  python mts_plots.py --tests flexural
  python mts_plots.py --no-toe             # disable D638 Annex A1 toe compensation
  python mts_plots.py --keep-tare          # leave the flexural load-cell offset in
  python mts_plots.py --large-deflection   # apply D790 §12.3 to the flexural stress
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IN2MM = 25.4

# The ASTM section marker (§) appears in the console output and in --help.
# A real Windows console handles it, but a piped or redirected stdout defaults
# to cp1252, which cannot encode it. Make stdout UTF-8 so the citations survive
# `python mts_plots.py > log.txt`.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# =============================================================================
# SCALAR PARAMETERS
#
# Every number this script can be tuned by, gathered in one place so nothing
# further down has to be edited to change the fixture geometry, the reference
# lengths, the ASTM reference ordinates or the toe-compensation search.
# =============================================================================

# --- Fixture geometry --------------------------------------------------------
# FLEX_SPAN_IN — support span L of the 3-point bend fixture. CONFIRMED at
#   8.00 in against the fixture, which is D790 §7.2's 16:1 span-to-depth for the
#   nominal 0.50 in coupon depth. Must match FLEX_SPAN_MM in FlexuralDIC_Level1.py
#   and FlexuralDIC_Level2.py; if you ever change it, change it in all three.
FLEX_SPAN_IN = 8.00
FLEX_SPAN_MM = FLEX_SPAN_IN * IN2MM                  # L = 203.2 mm

# TENSILE_GRIP_SEPARATION_IN — the ORIGINAL GRIP SEPARATION of the tensile
#   coupons, measured at 9.00 in. This is the length D638 §3.2.5 refers its
#   nominal strain to, so the right-hand tensile panel is a D638 nominal strain:
#
#       eps_nom = (change in grip separation) / (original grip separation)
#
#   It is still crosshead travel, so it also carries grip slip and load-frame
#   compliance and will read higher than an extensometer strain; for a D638
#   §11.4 strain use the DIC pipeline. Set this to None if the grip separation
#   is ever in doubt, and the axis falls back to TENSILE_DIC_GAUGE_IN and is
#   labelled as not being a D638 strain.
TENSILE_GRIP_SEPARATION_IN: float | None = 9.00

# TENSILE_DIC_GAUGE_IN — fallback reference length, the DIC axial gauge length
#   (AXIAL_GAUGE_IN in TensileDIC_Level1.py). Used only when the grip separation
#   above is None, so that the axis is at least on the same scale as the DIC
#   record and the two can be laid side by side.
TENSILE_DIC_GAUGE_IN = 4.36

TENSILE_GRIP_SEPARATION_MM = (None if TENSILE_GRIP_SEPARATION_IN is None
                              else TENSILE_GRIP_SEPARATION_IN * IN2MM)
TENSILE_DIC_GAUGE_MM  = TENSILE_DIC_GAUGE_IN * IN2MM           # 110.7 mm
TENSILE_REF_LENGTH_MM = (TENSILE_GRIP_SEPARATION_MM
                         if TENSILE_GRIP_SEPARATION_MM is not None
                         else TENSILE_DIC_GAUGE_MM)            # 228.6 mm
TENSILE_REF_IS_D638   = TENSILE_GRIP_SEPARATION_MM is not None

# --- ASTM reference ordinates drawn on the stress-strain panels ---------------
FLEX_STRAIN_LIMIT_PCT   = 5.0    # D790 §12.2 / §3.2.7 — 5 % flexural strain
BEARING_DEFORM_PCT      = 4.0    # D953 §13.3 — 4 % hole deformation
FLEX_DEFL_SPAN_LIMIT    = 0.10   # D790 §12.3 — D/L above which the correction applies

# --- MTS file format ---------------------------------------------------------
MTS_HEADERS = 8                  # header rows in the MTS .txt before row 1 of data

# --- Toe compensation search (ASTM D638 Annex A1; see the block below) --------
TOE_SEARCH_LOAD_BAND = (0.05, 0.60)   # fraction of peak load to search within
TOE_WINDOW_FRAC      = 0.15           # fit-window width, as a fraction of that band
TOE_WINDOW_MIN_PTS   = 25             # never fit a window shorter than this


# =============================================================================
# ASTM EQUATIONS
#
# Written out here once, in the standards' own symbols, so the reductions
# further down can be checked against them line by line. Each is repeated as a
# comment where it is actually evaluated.
#
# A note on the section/equation numbers: they follow the same convention as
# FlexuralDIC_Level2.py so the two scripts agree with each other. Different
# revisions of D790 renumber §12 (some insert the large-support-span stress as
# its own numbered equation, which shifts everything after it), so treat the
# section *names* as authoritative and re-check the numbers against your copy.
#
# -----------------------------------------------------------------------------
# ASTM D638-22 — TENSILE
# -----------------------------------------------------------------------------
#   §11.2  Tensile stress
#
#             sigma = P / A_0
#
#          P    = load carried by the specimen                          [N]
#          A_0  = *original* average cross-sectional area of the gauge
#                 section, = width x thickness                          [mm^2]
#          The original area is used at every point on the curve; D638 does not
#          ask for a true-stress correction.
#
#   §3.2.5 Nominal strain
#
#             eps_nom = (change in grip separation) / (original grip separation)
#
#          D638's crosshead-based strain is referred to the ORIGINAL GRIP
#          SEPARATION, measured at 9.00 in for this batch — see
#          TENSILE_GRIP_SEPARATION_IN in SCALAR PARAMETERS above. The
#          right-hand tensile axis is therefore a D638 nominal strain, but it
#          is still crosshead travel and so carries grip slip and load-frame
#          compliance; a D638 §11.4 strain needs the extensometer/DIC record.
#
#   Annex A1 (A1.2/A1.3)  Toe compensation
#
#          The toe region is take-up of slack plus seating of the specimen, not
#          a property of the material. Extend the linear (Hookean) portion of
#          the curve back to the zero-load axis; that intercept is the corrected
#          zero point, and all extensions are measured from it.
#
# -----------------------------------------------------------------------------
# ASTM D790-17 — FLEXURAL, PROCEDURE A, 3-POINT LOADING
# -----------------------------------------------------------------------------
#   §12.2 Eq.3  Flexural stress (stress in the outer fibre at midspan)
#
#             sigma_f = 3 P L / (2 b d^2)
#
#          P = load at a point on the load-deflection curve              [N]
#          L = support span                                              [mm]
#          b = width of beam tested                                      [mm]
#          d = depth of beam tested                                      [mm]
#
#   §12.3       Flexural stress for beams tested at large support spans
#
#             sigma_f = (3 P L / (2 b d^2)) * [1 + 6 (D/L)^2 - 4 (d/L)(D/L)]
#
#          D = deflection of the centreline of the specimen at midspan   [mm]
#          Required when the span-to-depth ratio is greater than 16:1 such that
#          deflections exceed 10 % of the support span. This fixture is 16:1 and
#          the flexural records reach D/L = 0.067 at peak, so the correction is
#          NOT required here; --large-deflection applies it anyway if you want
#          to see the size of it. The script prints max D/L either way.
#
#   §12.4 Eq.5  Flexural strain (nominal fractional change in length of the
#               outer surface at midspan)
#
#             eps_f = 6 D d / L^2
#
#          This is small-deflection beam kinematics. It is exactly what the
#          flexural pipeline's eps_deflection channel computes, so the right-hand
#          panel here is directly comparable to it (the difference is that D is
#          crosshead travel here, DIC midspan deflection there).
#
#   §3.2.7      Flexural strength = the maximum flexural stress the specimen
#               sustains. Valid only when the specimen breaks (or the load
#               drops) before 5 % strain; if it does not, D790 §12.2 asks for
#               the stress AT 5 % strain instead. The 5 % line is drawn on the
#               flexural panel and the per-specimen table flags any coupon that
#               reaches it.
#
#   §12.1       Toe compensation — refers to Annex A1 of D638, as above.
#
#   §7.2        Support span-to-depth ratio 16:1 for this material class.
#
# -----------------------------------------------------------------------------
# ASTM D953-19 — PIN-BEARING, PROCEDURE A
# -----------------------------------------------------------------------------
#   §13.2  Bearing stress — load divided by the PROJECTED bearing area, i.e.
#          the hole diameter times the specimen thickness:
#
#             sigma_b = P / (D_hole * t)
#
#          P      = load carried by the specimen                         [N]
#          D_hole = diameter of the loaded hole                          [mm]
#          t      = thickness of the specimen at the hole                [mm]
#
#   §13.3  Bearing stress at 4 % hole deformation — the bearing stress reached
#          when the hole has deformed by 4 % of its original diameter. The
#          abscissa here is therefore normalised as
#
#             hole deformation (%) = 100 * (deformation) / D_hole
#
#          and the 4 % ordinate is drawn as a reference line on the panel.
#
#   §3.2.5 Bearing strength = maximum bearing stress carried by the specimen.
# =============================================================================


# ---------------------------------------------------------------------------
# Paths
#
# The specimen geometry lives in FSR-SpecimenTesting.csv one level above the
# coupon folders. It used to be an .xlsx with a CSV export beside it; the
# workbook is retired (its formula columns kept coming back blank once a Level
# 2 had saved it through openpyxl), so the CSV is now the sheet itself — read
# here, and written by the Level 2s.
# ---------------------------------------------------------------------------
SCRIPT_DIR     = Path(__file__).resolve().parent
COUPONS_ROOT   = SCRIPT_DIR.parent
SPECIMEN_STEM  = COUPONS_ROOT.parent / "FSR-SpecimenTesting"
SPECIMEN_CSV   = SPECIMEN_STEM.with_suffix(".csv")

# The CSV started life as a Windows Excel export, so its degree/superscript
# characters ("Computed Area (in^2)") may still be cp1252 rather than UTF-8.
# Try in this order; the Level 2s re-write the file as utf-8-sig.
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")


# ---------------------------------------------------------------------------
# Toe compensation (ASTM D638 Annex A1)
#
# A1.3 says to construct a continuation of the LINEAR (Hookean) region of the
# curve through the zero-load axis and take that intercept as the corrected
# zero. The linear region is found here as the steepest straight segment of the
# load-deflection record, searched over the load band below — restricted so the
# search cannot land on the noise floor at the very start or on the roll-over
# near peak, both of which are outside the Hookean region by definition.
# The band and window are TOE_SEARCH_LOAD_BAND / TOE_WINDOW_* in SCALAR
# PARAMETERS above.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Display
#
# Colours are the dataviz categorical slots 1/2/3/7. The tab10 set used by the
# older figures (#1f77b4 / #ff7f0e / #17becf / #2ca02c) fails colour-vision
# separation: UV orange vs. IS green sit at OKLab dE 0.7 under protanopia, i.e.
# the same colour. This set clears the all-pairs CVD and normal-vision floors.
# Aqua (SW) is below 3:1 contrast on white, so the per-specimen stdout table is
# the required non-colour fallback.
# ---------------------------------------------------------------------------
EXPOSURE_ORDER  = ["CL", "UV", "SW", "IS"]
EXPOSURE_LABELS = {"CL": "Control", "UV": "UV", "SW": "Seawater", "IS": "SW+UV"}
EXPOSURE_COLORS = {"CL": "#2a78d6", "UV": "#eb6834", "SW": "#1baf7a", "IS": "#4a3aa7"}
# Long gaps on purpose: three replicates per group sit almost on top of each
# other, and a tight dash pattern (1.5 on / 1.8 off) has its gaps filled in by
# the neighbouring reps' dashes, so the group reads as a solid band.
DIR_STYLES      = {0: "-", 45: (0, (6.5, 3.5)), 90: (0, (1.5, 4.0))}

TEST_LETTERS = {"tensile": "T", "flexural": "F", "bearing": "B"}
TEST_TITLES  = {
    "tensile":  "Tensile",
    "flexural": "Flexural",
    "bearing":  "Bearing",
}

# Fixture malfunction on P01-BCL00-01 (coupon expanded inside the hole, no
# bushings); P01-BGM00 is the replacement Control 0 deg run.
BAD_SPECIMENS  = {"P01-BCL00-01"}
SUBSTITUTIONS  = {"P01-BGM00": ("B", "CL", 0)}   # stem -> (test letter, exposure, direction)

plt.rcParams.update({
    "font.size":        10,
    "axes.titlesize":   11,
    "axes.labelsize":   11,
    "axes.edgecolor":   "#8a8a85",
    "axes.linewidth":   0.8,
    "grid.color":       "#c9c9c4",
    "grid.linewidth":   0.6,
    "legend.frameon":   True,
    "legend.framealpha": 0.9,
    "legend.edgecolor": "#c9c9c4",
    "figure.facecolor": "white",
})


# ===========================================================================
# SPECIMEN SHEET
# ===========================================================================
def read_specimen_table() -> tuple[pd.DataFrame, Path]:
    """Read FSR-SpecimenTesting.csv.

    Returns (dataframe, path actually read). Raises SystemExit with every
    failure reported if it cannot be read, since without geometry there is no
    stress-strain panel to draw.
    """
    problems = []

    if SPECIMEN_CSV.exists():
        for enc in CSV_ENCODINGS:
            try:
                return pd.read_csv(SPECIMEN_CSV, encoding=enc), SPECIMEN_CSV
            except UnicodeDecodeError:
                continue
            except Exception as exc:                      # malformed CSV, locked file
                problems.append(f"{SPECIMEN_CSV.name} [{enc}]: {exc}")
                break
        else:
            problems.append(f"{SPECIMEN_CSV.name}: not decodable as "
                            f"{'/'.join(CSV_ENCODINGS)}")
    else:
        problems.append(f"{SPECIMEN_CSV.name}: not found")

    raise SystemExit("Could not read the specimen sheet:\n  "
                     + "\n  ".join(problems))


def load_specimens(print_id: str) -> tuple[pd.DataFrame, Path]:
    """Geometry columns for one print, indexed by Specimen ID, in mm.

    The sheet's "Width / Dia. (in)" column means a different dimension in each
    test, which is why it is carried through unrenamed and interpreted per test
    in the stress-strain functions:

        tensile   width of the gauge section              (b, 1.500 in)
        flexural  width of the beam                       (b, 1.000 in)
        bearing   diameter of the loaded hole             (D, 0.5625 in)

    "Measured Gauge Thickness (in)" is likewise the tensile/bearing thickness
    and the flexural beam DEPTH (d), and "Computed Area (in^2)" is thickness x
    width, which is the D638 original area A_0 for the tensile coupons only.
    """
    df, source = read_specimen_table()

    if "Print ID" not in df.columns or "Specimen ID" not in df.columns:
        raise SystemExit(f"{source.name} has no 'Print ID'/'Specimen ID' column "
                         f"(found: {list(df.columns)[:6]}...)")
    df = df[df["Print ID"] == print_id].copy()

    # The area header carries a superscript 2 whose byte value depends on the
    # export encoding, so match it loosely rather than by exact string.
    area_col = next((c for c in df.columns if "computed area" in c.lower()), None)
    if area_col is None:
        raise SystemExit(f"{source.name} has no 'Computed Area' column")

    df = df.rename(columns={
        "Measured Gauge Thickness (in)": "t_in",
        "Width / Dia. (in)":             "w_in",
        area_col:                        "area_in2",
    })
    for col in ("t_in", "w_in", "area_in2"):
        if col not in df.columns:
            raise SystemExit(f"{source.name} has no column mapping to '{col}'")
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["t_mm"]     = df["t_in"] * IN2MM
    df["w_mm"]     = df["w_in"] * IN2MM
    df["area_mm2"] = df["area_in2"] * IN2MM ** 2
    return df.set_index("Specimen ID"), source


# ===========================================================================
# FILENAME -> (exposure, direction)
# ===========================================================================
def parse_stem(stem: str, test_letter: str) -> tuple[str, int] | None:
    """P01-TCL45-02 -> ('CL', 45). Returns None if the stem is not this test."""
    if stem in SUBSTITUTIONS:
        sub_letter, exp_code, d_int = SUBSTITUTIONS[stem]
        return (exp_code, d_int) if sub_letter == test_letter.upper() else None

    m = re.fullmatch(rf"[A-Z]\d+-{test_letter}([A-Z]{{2}})(\d{{2}})(?:-(\d+))?",
                     stem, flags=re.IGNORECASE)
    if not m:
        return None
    exp_code = m.group(1).upper()
    d_int    = int(m.group(2))
    if exp_code not in EXPOSURE_COLORS or d_int not in DIR_STYLES:
        return None
    return exp_code, d_int


# ===========================================================================
# LOAD-CELL TARE OFFSET
# ===========================================================================
def find_force_baseline(d: np.ndarray, f: np.ndarray) -> float:
    """Level of the leading flat run in `f`, or 0.0 if there isn't one.

    Not an ASTM correction — an instrumentation one, and it has to happen
    before the D638 Annex A1 toe compensation because a constant offset moves
    the Hookean line's zero-load intercept.

    The flexural records open on a constant 862 +/- 3 N (sd 12 N) that holds for
    ~1.6 mm of crosshead travel at 10-18 N/mm, against loading slopes of
    ~164 N/mm (0 deg) and ~85 N/mm (90 deg). Force cannot stay constant while
    the crosshead advances, and the same level shows up on all 12 specimens in
    both orientations, so it is the weight of the loading nose / compression
    platen hanging on an un-tared load cell, not specimen load. Left in, it
    inflates flexural strength by ~1.6x. Tensile and bearing records do not
    show it (they open at 43-184 N, and the -1745 N seen on some tensile runs is
    real grip-seating compression, not an offset), so this returns 0.0 for them
    and the caller only applies it where a plateau is actually detected.
    """
    n = len(f)
    if n < 200:
        return 0.0

    # Reference loading slope, taken at 40-70% of the rise. That band is inside
    # the ramp whether or not the record carries an offset, and a slope is
    # immune to the offset itself.
    f_max = f[-1]
    rise  = f_max - f[0]
    band  = (f >= f[0] + 0.40 * rise) & (f <= f[0] + 0.70 * rise)
    if band.sum() < 10:
        return 0.0
    load_slope = float(np.polyfit(d[band], f[band], 1)[0])
    if load_slope <= 0:
        return 0.0

    # Touchdown = first place the local slope reaches 40% of the loading slope.
    # Testing slope (not scatter about the first few hundred samples) is what
    # separates the two cases: a tensile record ramps from its very first
    # sample, so touchdown lands at index 0 and no offset is claimed, whereas
    # the flexural plateau creeps at 10-18 N/mm against 85-164 N/mm of loading.
    # Block MEANS, not every w'th sample: the force channel carries ~12 N of
    # noise over 0.07 mm sample spacing, which is ~125 N/mm of slope noise on
    # raw samples and would trip the test on its own.
    w  = max(25, n // 100)
    nb = n // w
    if nb < 8:
        return 0.0
    fs = f[:nb * w].reshape(nb, w).mean(axis=1)
    ds = d[:nb * w].reshape(nb, w).mean(axis=1)

    above = np.flatnonzero(np.gradient(fs, ds) > 0.40 * load_slope)
    if above.size == 0:
        return 0.0
    kb = int(above[0])
    k  = kb * w

    # A plateau must be a real run of samples spanning real travel, and be flat.
    if k < 50 or (d[k] - d[0]) < 0.2 or kb < 2:
        return 0.0
    baseline = float(np.median(f[:k]))
    if np.ptp(fs[:kb]) > 0.10 * (f_max - baseline):
        return 0.0

    # Finally, the offset has to be big enough to matter. Fixture take-up also
    # produces a slow leading rise, but there the force really is ~zero at the
    # origin, so its "plateau" level is a rounding error on the rise (0.8% on
    # the bearing runs) whereas a true tare is a large share of it (56% on the
    # flexural runs). Below this cut the correction would be noise, and take-up
    # is the toe correction's job, not this one's.
    if abs(baseline) < 0.05 * (f_max - baseline):
        return 0.0
    return baseline


# ===========================================================================
# TOE COMPENSATION  —  ASTM D638 Annex A1  (D790 §12.1 refers to it)
# ===========================================================================
def toe_compensation_shift(d: np.ndarray, f: np.ndarray) -> tuple[float, float]:
    """Corrected zero point on the displacement axis, per D638 Annex A1.3.

    A1.3, for a material with a Hookean region: "a continuation of the linear
    (CD) region of the curve is constructed through the zero-stress axis. This
    intersection (B) is the corrected zero-strain point from which all
    extensions or strains must be measured."

    So: find the linear region, fit a straight line to it, and return that
    line's zero-load intercept

        d_zero = -intercept / slope        for  f = slope * d + intercept

    The linear region is taken to be the steepest straight segment of the
    load-deflection curve, searched only within TOE_SEARCH_LOAD_BAND of peak
    load. The band matters: below it the record is on the load cell's noise
    floor and a spuriously steep window can appear there, and above it the
    curve has already rolled over, which is the non-Hookean part A1.3 excludes.

    Returns (d_zero, slope). d_zero is 0.0 and slope NaN if no usable linear
    region is found, which leaves the record untouched.
    """
    f_max  = float(f[-1])           # the record is truncated at peak, so f[-1] is peak
    f_span = f_max - float(f[0])
    if f_span <= 0:
        return 0.0, float("nan")

    # Search band, expressed on the load axis. Floored above f[0] so a record
    # that opens already loaded cannot collapse the band onto a near-flat run.
    lo_frac, hi_frac = TOE_SEARCH_LOAD_BAND
    lo = max(lo_frac * f_max, f[0] + lo_frac * f_span)
    hi = max(hi_frac * f_max, f[0] + hi_frac * f_span)
    band = np.flatnonzero((f >= lo) & (f <= hi))
    if band.size < TOE_WINDOW_MIN_PTS:
        return 0.0, float("nan")

    i0, i1 = int(band[0]), int(band[-1])
    n_band = i1 - i0 + 1
    width  = max(TOE_WINDOW_MIN_PTS, int(TOE_WINDOW_FRAC * n_band))
    if width > n_band:
        width = n_band

    # Slide the fit window across the band and keep the steepest fit. Step is a
    # quarter window, which is fine resolution for a segment this long and keeps
    # the loop obvious.
    step = max(1, width // 4)
    best_slope, best_intercept = float("nan"), float("nan")
    for start in range(i0, i1 - width + 2, step):
        seg = slice(start, start + width)
        slope, intercept = np.polyfit(d[seg], f[seg], 1)
        if not np.isfinite(best_slope) or slope > best_slope:
            best_slope, best_intercept = float(slope), float(intercept)

    if not np.isfinite(best_slope) or best_slope <= 0:
        return 0.0, float("nan")
    return -best_intercept / best_slope, best_slope


# ===========================================================================
# MTS RECORD -> (displacement, force), signs made positive
# ===========================================================================
def read_mts(path: Path, test: str, toe_correct: bool, remove_baseline: bool):
    """Read one MTS .txt and return the reduced machine channels.

    Returns (d_mm, f_N, d_raw_mm, f_raw_N, d_zero_mm, tare_N) or None.
    The _raw arrays are the same record with only the sign convention applied,
    kept so the tare and toe shifts can be reported against something.
    """
    raw = pd.read_csv(path, sep="\t", skiprows=MTS_HEADERS, header=None,
                      names=["disp_mm", "force_N", "output_V", "time_s"],
                      encoding="utf-8-sig", on_bad_lines="skip")
    raw = (raw.apply(pd.to_numeric, errors="coerce")
              .dropna(subset=["disp_mm", "force_N"]))
    if len(raw) < 10:
        return None

    d = raw["disp_mm"].to_numpy()
    f = raw["force_N"].to_numpy()

    # Flexural runs in compression: crosshead travels down and force is
    # negative. Flip both so every test type plots in the first quadrant.
    if test == "flexural":
        d, f = -d, -f
    d = d - d[0]

    # Truncate at peak force — everything past it is fracture + crosshead return.
    # D790 §3.2.7 and D953 §3.2.5 both define their strength at the maximum, and
    # D638's tensile strength is likewise the maximum load, so the peak is the
    # last point any of the three standards needs from these curves.
    i_peak = int(np.argmax(f))
    if i_peak < 5:
        return None
    d, f = d[:i_peak + 1], f[:i_peak + 1]
    d_raw, f_raw = d.copy(), f.copy()

    # 1. Instrument tare (see find_force_baseline) — must precede the toe fit.
    tare = find_force_baseline(d, f) if remove_baseline else 0.0
    f = f - tare

    # 2. D638 Annex A1 toe compensation — the corrected zero on the
    #    displacement axis. With the tare removed this shift is the approach
    #    travel before the loading nose touches down (~1.6 mm on the flexural
    #    runs) plus fixture take-up.
    d_zero = 0.0
    if toe_correct:
        d_zero, _slope = toe_compensation_shift(d, f)
    d = d - d_zero

    # Drop the non-contact approach travel that the toe compensation just pushed
    # to the left of the corrected origin: A1.3 measures all extensions FROM the
    # corrected zero, so negative extensions are not part of the record.
    keep = d >= 0
    if keep.sum() < 5:
        return None
    return d[keep], f[keep], d_raw, f_raw, d_zero, tare


# ===========================================================================
# STRESS / STRAIN  —  one function per standard, equations written out
# ===========================================================================
def tensile_stress_strain(d: np.ndarray, f: np.ndarray, row: pd.Series):
    """ASTM D638-22.  Returns (x, x_label, stress_MPa, y_label, extras)."""
    # --- D638 §11.2:  sigma = P / A_0,  A_0 = ORIGINAL area (width x thickness)
    A_0    = row["area_mm2"]
    stress = f / A_0

    # --- D638 §3.2.5 nominal strain is (change in grip separation) / (original
    #     grip separation). The grip separation was measured at 9.00 in, so
    #     this axis is a D638 nominal strain. If TENSILE_GRIP_SEPARATION_IN is
    #     set back to None, the axis divides crosshead travel by the DIC axial
    #     gauge length instead and becomes a machine-frame extension ratio,
    #     NOT a D638 strain — the label says which one it is.
    L_ref  = TENSILE_REF_LENGTH_MM
    strain = 100.0 * d / L_ref

    x_label = (f"Nominal strain (%)  [D638 §3.2.5, grip sep. {L_ref:.0f} mm]"
               if TENSILE_REF_IS_D638 else
               f"Crosshead extension / L$_{{ref}}$ (%)  "
               f"[L$_{{ref}}$ = {L_ref:.0f} mm — not a D638 strain]")
    return strain, x_label, stress, "Tensile stress (MPa)", {}


def flexural_stress_strain(d: np.ndarray, f: np.ndarray, row: pd.Series,
                           large_deflection: bool):
    """ASTM D790-17 Procedure A, 3-point.  Returns (x, x_label, stress, y_label, extras)."""
    P = f                       # load                                      [N]
    D = d                       # midspan deflection ~ crosshead travel     [mm]
    L = FLEX_SPAN_MM            # support span, confirmed 8.00 in           [mm]
    b = row["w_mm"]             # width of beam tested                      [mm]
    depth = row["t_mm"]         # depth of beam tested (d in D790)          [mm]

    # --- D790 §12.2 Eq.3:  sigma_f = 3 P L / (2 b d^2) -----------------------
    stress = 3.0 * P * L / (2.0 * b * depth ** 2)

    # --- D790 §12.3 large-support-span correction ----------------------------
    #     sigma_f = (3 P L / 2 b d^2) * [1 + 6 (D/L)^2 - 4 (d/L)(D/L)]
    #     Required only for spans greater than 16:1 with D/L > 0.10. This
    #     fixture is exactly 16:1 and D/L peaks at ~0.067, so it is off by
    #     default; --large-deflection turns it on to show the size of it.
    if large_deflection:
        stress = stress * (1.0 + 6.0 * (D / L) ** 2
                           - 4.0 * (depth / L) * (D / L))

    # --- D790 §12.4 Eq.5:  eps_f = 6 D d / L^2 -------------------------------
    strain = 100.0 * 6.0 * D * depth / L ** 2

    extras = {
        "span_to_depth":  L / depth,
        "defl_over_span": float(D[-1] / L),          # D/L at peak, for §12.4
        "reached_5pct":   bool(strain[-1] >= FLEX_STRAIN_LIMIT_PCT),
    }
    return strain, "Flexural strain (%)", stress, "Flexural stress (MPa)", extras


def bearing_stress_deformation(d: np.ndarray, f: np.ndarray, row: pd.Series):
    """ASTM D953-19 Procedure A.  Returns (x, x_label, stress, y_label, extras)."""
    D_hole = row["w_mm"]        # diameter of the loaded hole, 0.5625 in    [mm]
    t      = row["t_mm"]        # specimen thickness at the hole            [mm]

    # --- D953 §13.2:  sigma_b = P / (D_hole * t), the PROJECTED bearing area --
    stress = f / (D_hole * t)

    # --- D953 §13.3 works in hole deformation as a fraction of hole diameter;
    #     the 4 % ordinate is what the standard's bearing stress is quoted at.
    deform = 100.0 * d / D_hole

    extras = {"reached_4pct": bool(deform[-1] >= BEARING_DEFORM_PCT)}
    return (deform, f"Hole deformation (% of D = {D_hole:.2f} mm)",
            stress, "Bearing stress (MPa)", extras)


def stress_strain(test: str, d: np.ndarray, f: np.ndarray, row: pd.Series,
                  large_deflection: bool):
    """Dispatch to the right standard."""
    if test == "tensile":
        return tensile_stress_strain(d, f, row)
    if test == "flexural":
        return flexural_stress_strain(d, f, row, large_deflection)
    return bearing_stress_deformation(d, f, row)


# ===========================================================================
# ONE FIGURE FOR ONE TEST TYPE
# ===========================================================================
def collect_coupons(test: str, mts_dir: Path, spec: pd.DataFrame, print_id: str,
                    toe_correct: bool, remove_baseline: bool,
                    large_deflection: bool) -> list[dict]:
    letter  = TEST_LETTERS[test]
    coupons = []

    for fp in sorted(mts_dir.glob(f"{print_id}-*.txt")):
        stem   = fp.stem
        parsed = parse_stem(stem, letter)
        if parsed is None:
            continue
        exp_code, d_int = parsed

        if stem in BAD_SPECIMENS:
            print(f"  [skip] {stem} - flagged bad in spreadsheet")
            continue
        if stem not in spec.index:
            print(f"  [skip] {stem} - no row in specimen sheet")
            continue
        row = spec.loc[stem]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        if pd.isna(row["t_mm"]) or pd.isna(row["w_mm"]):
            print(f"  [skip] {stem} - missing geometry")
            continue
        if test == "tensile" and pd.isna(row["area_mm2"]):
            print(f"  [skip] {stem} - missing computed area (D638 A_0)")
            continue

        rec = read_mts(fp, test, toe_correct, remove_baseline)
        if rec is None:
            print(f"  [skip] {stem} - unusable record")
            continue
        d, f, _d_raw, _f_raw, d_zero, tare = rec

        strain, x_lab, stress, y_lab, extras = stress_strain(
            test, d, f, row, large_deflection)

        # Strength is the MAXIMUM stress the specimen carried — D638 §11.2
        # (tensile strength at maximum load), D790 §3.2.7 (flexural strength)
        # and D953 §3.2.5 (bearing strength) all define it that way. That is
        # the last point of the record for every reduction except the D790
        # §12.3 corrected stress, which is not monotonic in load, so take the
        # argmax rather than assuming.
        i_max = int(np.nanargmax(stress))

        coupons.append({
            "id": stem, "exp": exp_code, "dir": d_int,
            "d": d, "f": f, "strain": strain, "stress": stress,
            "x_lab": x_lab, "y_lab": y_lab, "i_max": i_max,
            "d_zero": d_zero, "tare": tare, **extras,
        })

    return coupons


def print_table(test: str, coupons: list[dict]) -> None:
    """Per-specimen numbers — the required non-colour fallback for the figure,
    and where the ASTM validity flags are reported."""
    strength_label = {"tensile":  "sigma_max (MPa)",     # D638 §11.2 at P_max
                      "flexural": "sigma_fM (MPa)",      # D790 §3.2.7
                      "bearing":  "sigma_b,max (MPa)"}[test]
    x_label = {"tensile":  "x at peak (%)",
               "flexural": "eps_f at peak (%)",
               "bearing":  "def. at peak (%)"}[test]

    print(f"\n  {'specimen':<16}{'group':<10}{'dir':>4}"
          f"{'P_max (kN)':>12}{strength_label:>18}{x_label:>18}"
          f"{'toe (mm)':>10}{'tare (N)':>10}  flags")
    print("  " + "-" * 104)
    for c in sorted(coupons, key=lambda c: (EXPOSURE_ORDER.index(c["exp"]),
                                            c["dir"], c["id"])):
        flags = []
        if test == "flexural":
            # D790 §12.2: flexural strength is the maximum stress only if the
            # specimen breaks before 5 % strain; otherwise report stress at 5 %.
            if c.get("reached_5pct"):
                flags.append("reached 5% strain — D790 §12.2: quote sigma at 5%")
            # D790 §12.3: correction needed only above D/L = 0.10.
            if c.get("defl_over_span", 0.0) > FLEX_DEFL_SPAN_LIMIT:
                flags.append(f"D/L={c['defl_over_span']:.3f} > "
                             f"{FLEX_DEFL_SPAN_LIMIT:.2f} — D790 §12.3")
        if test == "bearing" and not c.get("reached_4pct", True):
            # D953 §13.3's 4 % bearing stress does not exist for this coupon.
            flags.append("failed before 4% deformation — D953 §13.3 n/a")

        i = c["i_max"]
        print(f"  {c['id']:<16}{EXPOSURE_LABELS[c['exp']]:<10}{c['dir']:>4}"
              f"{c['f'][-1] / 1000:>12.2f}{c['stress'][i]:>18.1f}"
              f"{c['strain'][i]:>18.2f}{c['d_zero']:>10.2f}"
              f"{c['tare']:>10.0f}  {'; '.join(flags)}")

    if test == "flexural":
        s_over_d = coupons[0].get("span_to_depth")
        dl_max   = max(c.get("defl_over_span", 0.0) for c in coupons)
        print(f"\n  D790 §7.2 span-to-depth L/d = {s_over_d:.1f} "
              f"(nominal 16:1, span {FLEX_SPAN_MM:.1f} mm confirmed)")
        print(f"  D790 §12.3 max D/L over this set = {dl_max:.3f} "
              f"({'ABOVE' if dl_max > FLEX_DEFL_SPAN_LIMIT else 'below'} the "
              f"{FLEX_DEFL_SPAN_LIMIT:.2f} threshold)")


def plot_test(test: str, coupons: list[dict], figs_dir: Path, print_id: str,
              toe_correct: bool, remove_baseline: bool,
              large_deflection: bool) -> Path:
    x_lab, y_lab = coupons[0]["x_lab"], coupons[0]["y_lab"]

    fig, (ax_fd, ax_ss) = plt.subplots(1, 2, figsize=(12, 5))

    for c in coupons:
        col = EXPOSURE_COLORS[c["exp"]]
        ls  = DIR_STYLES[c["dir"]]
        for ax, x, y in ((ax_fd, c["d"], c["f"] / 1000),
                         (ax_ss, c["strain"], c["stress"])):
            ax.plot(x, y, color=col, ls=ls, lw=1.4, alpha=0.9,
                    solid_capstyle="round")
            ax.plot(x[-1], y[-1], marker="o", ms=5, color=col,
                    mec="white", mew=1.0, zorder=5)

    # ---- left panel: the machine channels ----------------------------------
    ax_fd.set_xlabel("Crosshead displacement (mm)")
    ax_fd.set_ylabel("Force (kN)")
    processing = [s for s in (("tare removed" if remove_baseline else None),
                              ("toe-compensated" if toe_correct else None)) if s]
    ax_fd.set_title("Force – displacement"
                    + (f"\n({', '.join(processing)})" if processing else "\n(as recorded)"))

    # ---- right panel: the ASTM reduction -----------------------------------
    ax_ss.set_xlabel(x_lab)
    ax_ss.set_ylabel(y_lab)
    ss_title = {"tensile":  "Stress – strain",
                "flexural": "Flexural stress – strain",
                "bearing":  "Bearing stress – hole deformation"}[test]
    if test == "flexural" and large_deflection:
        ss_title += "\n(large-deflection correction applied)"
    ax_ss.set_title(ss_title)

    # ---- ASTM reference ordinate --------------------------------------------
    #   flexural  D790 §12.2 / §3.2.7 — 5 % strain, past which the maximum
    #             stress is no longer the flexural strength
    #   bearing   D953 §13.3        — 4 % hole deformation, the ordinate the
    #             standard's bearing stress is quoted at
    #
    # Drawn as a line only when the data gets near it. If every coupon fails
    # well short, a line at 5 % would squash all the curves into the left third
    # of the panel to show an empty region, so the fact is stated in a corner
    # note instead — the information is the same either way.
    ref_x, ref_txt = {
        "tensile":  (None, None),
        "flexural": (FLEX_STRAIN_LIMIT_PCT, "D790 5 % strain"),
        "bearing":  (BEARING_DEFORM_PCT,    "D953 4 % hole deformation"),
    }[test]
    if ref_x is not None:
        x_max = max(float(c["strain"][-1]) for c in coupons)
        if x_max >= 0.5 * ref_x:
            ax_ss.axvline(ref_x, color="#8a8a85", lw=1.0, ls=(0, (2, 2)), zorder=1)
            ax_ss.annotate(ref_txt, xy=(ref_x, 0.98),
                           xycoords=("data", "axes fraction"), xytext=(4, -2),
                           textcoords="offset points", fontsize=8.5,
                           color="#52514e", ha="left", va="top", rotation=90)
        else:
            ax_ss.annotate(f"{ref_txt} not reached (max {x_max:.2f} %)",
                           xy=(0.02, 0.98), xycoords="axes fraction",
                           fontsize=8.5, color="#52514e", ha="left", va="top")

    for ax in (ax_fd, ax_ss):
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.5, ls="--")
        ax.spines[["top", "right"]].set_visible(False)

    # ---- legend: colour = group, linestyle = orientation --------------------
    exps_here = [e for e in EXPOSURE_ORDER if any(c["exp"] == e for c in coupons)]
    dirs_here = sorted({c["dir"] for c in coupons})
    exp_handles = [mlines.Line2D([], [], color=EXPOSURE_COLORS[e], lw=2.4,
                                 label=EXPOSURE_LABELS[e]) for e in exps_here]
    dir_handles = [mlines.Line2D([], [], color="#52514e", lw=1.6,
                                 ls=DIR_STYLES[d], label=f"{d}°")
                   for d in dirs_here]
    leg_exp = ax_ss.legend(handles=exp_handles, title="Exposure", fontsize=9,
                           title_fontsize=9, loc="lower right",
                           bbox_to_anchor=(1.0, 0.0))
    leg_exp._legend_box.align = "left"
    ax_ss.add_artist(leg_exp)
    leg_dir = ax_ss.legend(handles=dir_handles, title="Orientation", fontsize=9,
                           title_fontsize=9, loc="lower right",
                           bbox_to_anchor=(0.76, 0.0))
    leg_dir._legend_box.align = "left"

    fig.suptitle(f"{print_id} — {TEST_TITLES[test]}", fontsize=12.5)
    fig.tight_layout()

    out = figs_dir / f"mts_{test}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


# ===========================================================================
# MAIN
# ===========================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--print", dest="print_id", default="P01",
                    help="print ID to plot (default P01)")
    ap.add_argument("--no-toe", dest="toe", action="store_false",
                    help="skip the D638 Annex A1 toe compensation")
    ap.add_argument("--keep-tare", dest="baseline", action="store_false",
                    help="keep the load-cell tare offset (see find_force_baseline); "
                         "inflates flexural strength ~1.6x")
    ap.add_argument("--large-deflection", action="store_true",
                    help="apply the D790 §12.3 large-support-span stress "
                         "correction (not required at this fixture's 16:1 span "
                         "and D/L < 0.10 — see the ASTM EQUATIONS block)")
    ap.add_argument("--tests", nargs="+", default=list(TEST_LETTERS),
                    choices=list(TEST_LETTERS), help="subset of tests to plot")
    args = ap.parse_args()

    print_dirs = sorted(COUPONS_ROOT.glob(f"{args.print_id}-*"))
    print_dirs = [p for p in print_dirs if (p / "MTS").is_dir()]
    if not print_dirs:
        raise SystemExit(f"No '{args.print_id}-*/MTS' folder under {COUPONS_ROOT}")
    print_dir = print_dirs[0]

    mts_dir  = print_dir / "MTS"
    figs_dir = print_dir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)

    spec, source = load_specimens(args.print_id)
    print(f"{args.print_id}: {len(spec)} specimen rows from {source.name}")
    print(f"MTS dir: {mts_dir}")
    print(f"flexural span L = {FLEX_SPAN_MM:.1f} mm ({FLEX_SPAN_MM / IN2MM:.2f} in) "
          f"— confirmed fixture setting")
    print(f"tensile reference length = {TENSILE_REF_LENGTH_MM:.1f} mm "
          f"({TENSILE_REF_LENGTH_MM / IN2MM:.2f} in) — "
          + ("measured grip separation, D638 §3.2.5 nominal strain"
             if TENSILE_REF_IS_D638 else
             "DIC gauge length; NOT a D638 strain"))
    print(f"toe compensation (D638 Annex A1): {'on' if args.toe else 'off'}   "
          f"tare removal: {'on' if args.baseline else 'off'}   "
          f"D790 §12.3 correction: {'on' if args.large_deflection else 'off'}")

    saved = []
    for test in args.tests:
        print(f"\n=== {test} — {TEST_TITLES[test]} ===")
        coupons = collect_coupons(test, mts_dir, spec, args.print_id,
                                  args.toe, args.baseline, args.large_deflection)
        if not coupons:
            print(f"  no {test} records found")
            continue
        print_table(test, coupons)
        saved.append(plot_test(test, coupons, figs_dir, args.print_id,
                               args.toe, args.baseline, args.large_deflection))

    print()
    for p in saved:
        print(f"Saved -> {p}")


if __name__ == "__main__":
    main()
