#!/usr/bin/env python3
"""
FlexuralDIC_Level1.py  —  FSR Flexural Coupons (ASTM D790, 3-point bend)
=========================================================================
Step A — exports each coupon's raw VIC-3D .out files to per-frame CSVs
(next to the .out files, on the raw data drive), exactly as TensileDIC_Level1
Step A does and with the same columns. ~1 GB per coupon.
Step B — reduces every .out frame to the bending kinematics (midspan
deflection, curvature, neutral axis, extreme-fibre strains), pairs them with
the MTS force/displacement record, and writes one per-coupon CSV (full,
untruncated record) to DIC_DIR.
Step C — reports which variables each coupon's .out files actually contain.
Cheap and vicpyx-only; run it after a VIC-3D reprocess to confirm new
inspector items are reaching the .out export before spending an hour on
Step B.

No truncation and no property extraction here — that's Level 2's job, so it
can be tuned without re-running the slow per-.out pass.

WHY THIS ISN'T TensileDIC_Level1 WITH DIFFERENT CONSTANTS
---------------------------------------------------------
The tensile ROI is the gauge section and each frame reduces to two virtual
point extensometers. The flexural ROI is the *side profile* — the specimen
depth seen edge-on over the whole span — so each frame carries the entire
bending field, and the useful reduction is different:

  defl_mm       midspan deflection, measured against the chord through the
                two supports. Referencing to the supports rather than the
                world frame subtracts rigid settling and fixture tilt, which
                the crosshead cannot separate from real deflection.
  kappa_1pmm    curvature, from a straight-line fit of exx against Y across
                the depth at midspan. DIC-native: no span, no fixture, no
                machine compliance in it.
  na_Y_mm       that same fit's zero crossing — the neutral axis. Its drift
                away from mid-depth is the check on whether the specimen is
                really in pure bending.
  eps_bot/top   the fit extrapolated to the two faces. The ROI is inset from
                both faces by the correlation subset radius (~61 % of the
                depth is covered here), so the surface strain D790 asks for
                is never directly measured — it has to come off the fit.

THE SYNC CSV IS NOT A DATA SOURCE HERE — ONLY A CLOCK
------------------------------------------------------
Every P01 coupon (tensile, flexural and bearing alike) has bad analog data in
its per-coupon VIC sync CSV: there was a connection issue during the batch.
On the tensile side that was worked around by anchoring the two peak forces
and back-tracking to the start of test (see TensileDIC_Level2_tmp.py). The
flexural files are worse — there is no load channel at all. Dev1/ai2, the
input the DAQ config labels LOAD, is a scaled copy of Dev1/ai1 (the
displacement input) at r > 0.99: nothing was wired to it, so what appears on
it is multiplexer crosstalk from the input beside it. It shows no touchdown
knee where the MTS sits flat on its tare, and it ramps at 81-103 % of its
loaded rate during the approach travel. A load cell cannot do that.

So this script reads exactly one thing out of the sync CSV — Time_0_0, the
frame trigger clock — and takes force and crosshead travel entirely from the
MTS record. None of the analog channels are carried into the per-coupon CSV.

ALIGNMENT: THE LAST CORRELATED FRAME
-------------------------------------
The MTS record and the DIC frames have independent clocks and independent
start triggers, so an offset has to be found, and with no usable sync load
channel only a *feature* can supply one. Scanning the offset to minimise a
regression residual cannot work even in principle: force is very nearly
linear in time in both records, and a time shift between two straight lines
is absorbed exactly by the regression's intercept.

Fracture is the one sharp event in the test, and on this specimen type it can
be read straight off the images. A bend specimen decorrelates the instant it
cracks — every point in the frame goes invalid at once, not gradually — so the
last correlated frame IS the fracture instant, measured with no analog channel
in it. That frame is aligned to the MTS force peak. See find_break_frame() and
align_mts_to_frames().

The crosshead-travel cross-check is reported: the sync CSV's displacement
channel (ai1, the one input that is really connected) and the MTS crosshead
record measure the same motion through separate paths, so after alignment they
should agree. That does not independently confirm the time offset — a
constant-rate ramp cannot — but a large residual is a loud signal that the two
files are not the same test.

LOAD-CELL TARE
--------------
The flexural records open on a constant ~860 N held through the approach
travel: the loading nose hanging on an un-tared cell. It shows on all 12
flexural specimens in both orientations, and left in it inflates flexural
strength ~1.6x. find_force_baseline() detects and removes it (ported from
mts_plots.py). Tensile and bearing records don't have it.

SPAN
----
FLEX_SPAN_MM = 203.2 (8.00 in, D790's 16:1 for the nominal 0.50 in depth).
CONFIRMED against the fixture — it was an open assumption until then, since
nothing in the repo or the workbook recorded the setting.

Two related numbers that are not the same thing, and both are right:

  203.2 mm   the fixture setting: where the support rollers are. D790 defines
             flexural stress and strain on this nominal span, so it is what
             every formula here uses.
  ~200 mm    where the bending moment actually crosses zero at low load,
             measured off the DIC curvature diagram. The contact points sit a
             little inboard of the roller centres and ride further inboard as
             the beam rotates on them (~-6 mm/kN), so the effective span
             shrinks under load.

The gap between them is contact geometry, not a fixture error. It is one of the
things feeding the eps_deflection/eps_curvature disagreement Level 2 reports.
Level 2's curvature strain channel has no span in it at all and stays the
reference regardless.

USAGE
-----
1. Toggle the SWITCHES section below to pick which coupons to run.
2. Make sure vicpyx is installed in the active Python environment:
       pip install vicpyx
3. Run:
       python FlexuralDIC_Level1.py

INPUT per coupon
  <coupon_dir>/*.out            VIC-3D full-field export, one per DIC frame
  <coupon_dir>/<folder>.csv     VIC sync CSV — read for its frame clock only
  <MTS_DIR>/<coupon_id>.txt     MTS raw: disp_mm, force_N, output_V, time_s
  FSR-SpecimenTesting.csv       depth d and width b (read only, never written)

OUTPUTS
  <coupon_dir>/<out_filename>.csv   Step A: one CSV per .out, written next to
                                    it, same EXPORT_VARS columns TensileDIC_
                                    Level1 writes.
  <DIC_DIR>/<coupon_id>.csv         full per-frame record: step, time_s,
                                    force_N, disp_mts_mm, n_pts, defl_mm,
                                    kappa_1pmm, na_Y_mm, profile_r2, eps_bot,
                                    eps_top, eps_membrane. Level 2 reads this
                                    same file and appends its own columns.
  <DIC_DIR>/coupon_scalars.csv      one row per coupon: b, d, the located
                                    fixture, the tare, the alignment offset
                                    and its cross-check. The only place these
                                    per-coupon scalars are stored.

ONE FILE LAYOUT FOR BOTH TEST TYPES
-----------------------------------
This script used to write its per-coupon scalars to a flexural_geometry.csv of
its own. It no longer does: they go into <DIC_DIR>/coupon_scalars.csv, the same
file TensileDIC_Level1 writes, keyed by coupon. The two test types contribute
different columns — tensile mts_peak_N and area_mm2, flexural the fixture and
alignment block — so the table is sparse, and that is the point: one place to
look up any coupon's scalars, tensile or flexural, instead of one file per test
type. The upsert is by coupon ID, so neither script disturbs the other's rows.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from vicpyx import VicDataSet as VICDataSet

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# PATHS
# =============================================================================
MTS_DIR = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\MTS"
)
DIC_DIR  = MTS_DIR.parent / "DIC"   # per-coupon CSVs land here, next to MTS/
RAW_ROOT = DIC_DIR / "raw" / "2026_FSR_Flexural_FCL_FIS"
# Specimen geometry. THE CSV *IS* THE SHEET — there is no .xlsx any more.
#
# 'Width / Dia. (in)' used to be a FORMULA cell in FSR-SpecimenTesting.xlsx.
# openpyxl does not evaluate formulas, so every time a Level 2 wrote its
# scalars back into the workbook it saved the formula and dropped the cached
# value; pandas and openpyxl both read the CACHE, so that column read blank
# until somebody opened the workbook in Excel and saved it, and
# specimen_geometry() then returned b_mm = NaN with every stress downstream
# NaN too. So the workbook is retired. FSR-SpecimenTesting.csv holds evaluated
# values, needs no Excel engine, and is not locked while something else has it
# open — mts_plots.py and TensileDIC_Level1.py read it the same way, and the
# Level 2s write their scalars back into it.
SPECIMEN_STEM  = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\FSR-SpecimenTesting"
)
SPECIMEN_CSV   = SPECIMEN_STEM.with_suffix(".csv")

# The CSV started life as a Windows Excel export, so it may still be cp1252
# rather than UTF-8; the Level 2s re-write it as utf-8-sig.
CSV_ENCODINGS = ("utf-8-sig", "cp1252", "latin-1")

# =============================================================================
# SWITCHES — toggle which coupons to process
# =============================================================================
PRINTS     = ["P01"]
EXPOSURES  = {"CL": True, "IS": False}      # only CL and IS were bend-tested
DIRECTIONS = {"00": True, "90": True}
REPLICATES = ["01", "02", "03"]

DO_LIST_VARS     = True     # Step C: report the variables in each coupon's .out files
DO_EXPORT_FRAMES = True     # Step A: dump each .out to a CSV next to it (~1 GB/coupon)
OVERWRITE_FRAMES = False    # if False, skip .out files whose .csv already exists
DO_BUILD_L1      = True     # Step B: reduce frames + pair with MTS
OVERWRITE_L1     = False    # if False, skip coupons whose per-coupon CSV already exists

# =============================================================================
# .OUT VARIABLES
#
# REQUIRED_VARS must be present or the frame is unusable. OPTIONAL_VARS are
# read when they exist and silently skipped when they don't — that is where
# any new VIC-3D inspector items land once the projects are reprocessed, so
# adding one here costs nothing on coupons that predate it. Step C prints what
# each coupon actually has; check there first after a reprocess.
#
# EXPORT_VARS is a separate, wider list used only by Step A. Step B needs five
# variables; the per-frame CSVs are a general-purpose export that other tools
# read, so they carry everything the .out has — the same list TensileDIC_Level1
# exports, so a flexural per-frame CSV and a tensile one have the same columns.
# =============================================================================
REQUIRED_VARS = ["sigma", "X", "Y", "V", "exx"]
OPTIONAL_VARS = ["Z", "U", "W", "eyy", "exy"]

EXPORT_VARS = [
    "sigma",                              # correlation confidence (filter on this)
    "X", "Y", "Z",                        # world coords (mm)
    "U", "V", "W",                        # displacements (mm)
    "exx", "eyy", "exy",                  # in-plane strains
    "e1", "e2", "gamma",                  # principal & max-shear strains
    "x", "y", "u", "v",                   # pixel coords / pixel disps
    "q", "r", "q_ref", "r_ref",           # subset coords
]

# =============================================================================
# FIXTURE / SPECIMEN GEOMETRY
# =============================================================================
IN2MM        = 25.4
FLEX_SPAN_MM = 8.00 * IN2MM     # L, support span — CONFIRMED against the fixture

# Set to a number to pin midspan in DIC X coordinates instead of locating it
# from the deflected shape. None = auto (recommended; it is a fixture check).
MIDSPAN_X_MM = None

# Fraction of peak load at which to grab the frame the fixture is located on.
# Deflected enough to measure, early enough to be undamaged.
FIXTURE_PROBE_LOAD_FRAC = 0.60

# Sampling windows on the specimen, in mm.
MIDSPAN_HALF_WIDTH_MM = 10.0   # +/- about midspan, for the exx-vs-Y curvature fit
SUPPORT_HALF_WIDTH_MM = 4.0    # +/- about each support, for the deflection chord
PROFILE_BIN_MM        = 4.0    # X bin width when reducing V(X) to a profile

# Drop points this close to the top/bottom edge of the ROI before fitting exx
# vs Y. The outermost row of subsets is only partly on the specimen and
# consistently pulls low.
PROFILE_EDGE_TRIM_MM = 0.30

# The ideal 3-point deflected shape corrects for the support sampling windows
# sitting a few mm inboard of the supports themselves, where the true
# deflection is not quite zero. ~2-5 % at the default window width.
APPLY_SUPPORT_OFFSET_CORRECTION = True

MIN_VALID_POINTS = 200   # frames with fewer correlated points count as uncorrelated

# ROI ORIENTATION CHECK
#
# Every reduction in this script assumes the VIC-3D world frame is oriented the
# same way on every coupon: X runs ALONG the specimen (the span) and Y runs
# THROUGH THE DEPTH. That is not automatic — it comes from the calibration and
# the alignment plane chosen when the project was built, and three coupons in
# this batch (FCL9003, FIS0001, FIS9001) were reconstructed on a different
# frame: on two of them the span and the depth are swapped outright.
#
# Nothing downstream can detect that on its own. The curvature fit happily
# regresses exx against a coordinate that isn't depth, returns R^2 ~ 0 and a
# neutral axis tens of mm off the specimen, and Level 2 turns that into a
# plausible-looking row of numbers. So it is checked here and the coupon is
# skipped rather than written — a missing coupon is recoverable, a silently
# wrong one is not.
#
# The test: the ROI's Y extent must be a sensible fraction of the specimen
# depth (the subset radius insets it from both faces, ~61 % here), and its X
# extent must be much larger, since the ROI spans most of a 203 mm span against
# a 12.5 mm depth.
ROI_DEPTH_FRAC_RANGE = (0.30, 1.10)   # ROI Y extent / specimen depth d
ROI_MIN_ASPECT       = 3.0            # ROI X extent / ROI Y extent

# =============================================================================
# SYNC / MTS FILE FORMAT
# =============================================================================
HEADERS = 8                 # MTS .txt header lines

# The ONLY column read out of the sync CSV. Everything else in that file is
# either unconnected or crosstalk — see the module docstring.
SYNC_TIME_COL = "Time_0_0"

# All flexural DIC recordings should be 5 Hz. Checked per coupon rather than
# assumed; a coupon off by more than the tolerance gets a warning and is still
# processed on its own measured clock.
EXPECTED_FRAME_RATE_HZ = 5.0
FRAME_RATE_TOL_HZ      = 0.05


# =============================================================================
# HELPERS — coupon selection and paths
# =============================================================================
def coupon_id(p, e, d, r): return f"{p}-F{e}{d}-{r}"

def selected_coupons():
    return [coupon_id(p, e, d, r)
            for p in PRINTS
            for e, on in EXPOSURES.items() if on
            for d, on2 in DIRECTIONS.items() if on2
            for r in REPLICATES]

def raw_folder(cid: str) -> Path:
    """P01-FCL00-01 -> <RAW_ROOT>/FCL0001.

    The flexural VIC-3D projects drop the print prefix and the dashes, so
    unlike tensile the coupon ID cannot be used to find the folder directly.
    """
    _, group, rep = cid.split("-")        # "P01", "FCL00", "01"
    return RAW_ROOT / f"{group}{rep}"

def find_out_files(cdir: Path) -> list[Path]:
    return sorted(cdir.rglob("*.out"))

def find_sync_csv(cdir: Path) -> Path | None:
    p = cdir / f"{cdir.name}.csv"
    return p if p.exists() else None

def find_mts_txt(cid: str) -> Path | None:
    p = MTS_DIR / f"{cid}.txt"
    if p.exists():
        return p
    hits = sorted(MTS_DIR.glob(f"{cid}*.txt"))
    return hits[0] if hits else None

def load_mts_txt(fp: Path) -> pd.DataFrame:
    """MTS .txt: 8-line header, tab separated, cols disp_mm force_N output_V time_s."""
    return (pd.read_csv(fp, sep="\t", skiprows=HEADERS, header=None,
                        names=["disp_mm", "force_N", "output_V", "time_s"],
                        encoding="utf-8-sig", on_bad_lines="skip")
              .apply(pd.to_numeric, errors="coerce")
              .dropna(subset=["force_N", "disp_mm"]))

def _read_specimen_table() -> tuple[pd.DataFrame, Path]:
    """Specimen sheet as (dataframe, path actually read).

    See the SPECIMEN_CSV block near the top for why this reads a CSV and not a
    workbook. Same reader TensileDIC_Level1 and mts_plots use; duplicated
    rather than imported so each script stays runnable on its own, as
    FLEX_SPAN_MM already is.
    """
    problems = []
    if SPECIMEN_CSV.exists():
        for enc in CSV_ENCODINGS:
            try:
                return pd.read_csv(SPECIMEN_CSV, encoding=enc), SPECIMEN_CSV
            except UnicodeDecodeError:
                continue
            except Exception as exc:              # malformed CSV, locked file
                problems.append(f"{SPECIMEN_CSV.name} [{enc}]: {exc}")
                break
        else:
            problems.append(f"{SPECIMEN_CSV.name}: not decodable as "
                            f"{'/'.join(CSV_ENCODINGS)}")
    else:
        problems.append(f"{SPECIMEN_CSV.name}: not found")

    raise RuntimeError("Could not read the specimen sheet:\n  "
                       + "\n  ".join(problems))


def specimen_geometry(cid: str) -> tuple[float, float]:
    """Return (b, d) in mm: width and depth of the beam, from the specimen sheet.

    For a bend test the DIC ROI views the b x L face edge-on, so 'depth' d is
    the sheet's Measured Gauge Thickness (0.50 in nominal) and 'width' b is
    Width / Dia. (1.00 in nominal).
    """
    df, src = _read_specimen_table()
    t_col = next(c for c in df.columns if "thickness" in c.lower())
    w_col = next(c for c in df.columns if "width" in c.lower() and "dia" in c.lower())
    row = df.loc[df["Specimen ID"] == cid]
    if row.empty:
        raise RuntimeError(f"{cid} not in {src.name}")
    row = row.iloc[0]
    b_in = pd.to_numeric(row[w_col], errors="coerce")
    d_in = pd.to_numeric(row[t_col], errors="coerce")
    # Refuse rather than returning NaN: a NaN b or d makes every stress NaN
    # further down with nothing saying why.
    if not (np.isfinite(b_in) and np.isfinite(d_in)):
        raise RuntimeError(
            f"{cid}: missing geometry in {src.name} "
            f"(b={row[w_col]!r}, d={row[t_col]!r}).\n"
            f"    Fix: fill '{w_col}' and '{t_col}' in {SPECIMEN_CSV.name}. "
            f"That file is a plain CSV — edit it directly.")
    return float(b_in) * IN2MM, float(d_in) * IN2MM


# =============================================================================
# HELPERS — reading one .out frame
# =============================================================================
def frame_variables(fp: Path) -> list[str]:
    """The variable names one .out actually carries. Step C's whole content."""
    ds = VICDataSet()
    try:
        ds.load(str(fp))
        return list(ds.variables())
    except Exception as ex:
        print(f"    [!] could not read variables from {fp.name}: {ex}")
        return []


def read_frame(fp: Path) -> pd.DataFrame | None:
    """Load one .out and return its valid points, or None if uncorrelated.

    Post-fracture frames come back with every point invalid (sigma < 0) —
    correlation is lost the instant the specimen breaks — so they are dropped
    here rather than producing garbage downstream. That same property is what
    makes the last correlated frame a usable fracture timestamp.
    """
    ds = VICDataSet()
    try:
        ds.load(str(fp))
    except Exception as ex:
        print(f"    [!] could not load {fp.name}: {ex}")
        return None
    try:
        available = set(ds.variables())
    except Exception:
        available = set()
    missing = [v for v in REQUIRED_VARS if available and v not in available]
    if missing:
        print(f"    [!] {fp.name} is missing required variable(s) {missing}")
        return None
    wanted = REQUIRED_VARS + [v for v in OPTIONAL_VARS
                              if (not available) or (v in available)]
    try:
        vals = ds.get_values(wanted)
    except Exception as ex:
        print(f"    [!] get_values failed on {fp.name}: {ex}")
        return None
    mask = np.asarray(vals["sigma"]) >= 0
    if mask.sum() < MIN_VALID_POINTS:
        return None
    return pd.DataFrame({v: np.asarray(vals[v])[mask]
                         for v in wanted if v != "sigma"})


def check_roi_orientation(ref: pd.DataFrame, d_mm: float) -> str | None:
    """Return a complaint if the ROI isn't a side profile in the X-Y frame.

    See the ROI ORIENTATION CHECK note in the settings block above for why this
    exists. Returns None when the ROI looks right.
    """
    x_ext = float(ref["X"].max() - ref["X"].min())
    y_ext = float(ref["Y"].max() - ref["Y"].min())
    frac = y_ext / d_mm if d_mm else np.nan
    lo, hi = ROI_DEPTH_FRAC_RANGE
    if not (lo <= frac <= hi):
        return (f"ROI Y extent is {y_ext:.1f} mm = {frac:.1f}x the specimen depth "
                f"d = {d_mm:.2f} mm (expected {lo:g}-{hi:g}x)")
    if y_ext > 0 and x_ext / y_ext < ROI_MIN_ASPECT:
        return (f"ROI is {x_ext:.1f} mm along X by {y_ext:.1f} mm along Y — "
                f"aspect {x_ext / y_ext:.1f}, expected > {ROI_MIN_ASPECT:g}")
    return None


def find_break_frame(out_files: list[Path]) -> int:
    """Index of the last frame that still correlated — i.e. fracture.

    Scans backwards rather than bisecting: correlation is not guaranteed to be
    monotonic in the frames just before failure.
    """
    for i in range(len(out_files) - 1, -1, -1):
        if read_frame(out_files[i]) is not None:
            return i
    raise RuntimeError("no frame in this coupon correlated at all")


# =============================================================================
# HELPERS — load-cell tare
# =============================================================================
def find_force_baseline(d: np.ndarray, f: np.ndarray) -> float:
    """Level of the leading flat run in `f`, or 0.0 if there isn't one.

    Same detector as mts_plots.find_force_baseline. The flexural records
    open on a constant ~862 +/- 3 N holding over ~1.6 mm of crosshead travel,
    against loading slopes of ~164 N/mm (0 deg) and ~85 N/mm (90 deg): force
    cannot stay constant while the crosshead advances, and the same level shows
    on all 12 specimens in both orientations, so it is the loading nose hanging
    on an un-tared cell, not specimen load. Expects `f` already truncated at
    peak and sign-flipped positive.
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

    # Touchdown = first block where the local slope reaches 40% of the loading
    # slope. Block means, not raw samples: the force channel's own noise is
    # worth more slope than the plateau creep we are trying to detect.
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
# HELPERS — the frame clock and the MTS alignment
# =============================================================================
def frame_clock(sync_fp: Path, n_frames: int) -> tuple[np.ndarray, float]:
    """Frame trigger times (s, zeroed to the first frame) and the median rate.

    The sync CSV's analog channels are unusable on this batch, but Time_0_0 is
    a real epoch clock stamped on each frame trigger and is the only thing that
    puts the DIC frames on a time axis at all.
    """
    sync = pd.read_csv(sync_fp)
    if SYNC_TIME_COL not in sync.columns:
        raise RuntimeError(f"sync CSV missing {SYNC_TIME_COL!r}; "
                           f"has {list(sync.columns)}")
    t = pd.to_numeric(sync[SYNC_TIME_COL], errors="coerce").to_numpy()
    t = t[:n_frames] - t[0]
    dt = float(np.median(np.diff(t)))
    return t, (1.0 / dt if dt > 0 else np.nan)


def align_mts_to_frames(t_frames: np.ndarray, i_break: int,
                        t_mts: np.ndarray, f_mts: np.ndarray) -> float:
    """Seconds the MTS clock leads the frame clock by.

    Fracture is the anchor: the last correlated DIC frame against the MTS force
    peak. See the module docstring for why nothing else works here.
    """
    i_pk_mts = int(np.nanargmax(f_mts))
    return float(t_mts[i_pk_mts] - t_frames[i_break])


def _median_smooth(x: np.ndarray, win: int) -> np.ndarray:
    """Rolling median, odd window, edge-replicating. Used only on the sync
    displacement channel, whose noise is ~2 % of its full range."""
    win = win - 1 if win % 2 == 0 else win
    if win < 3 or len(x) < win:
        return np.asarray(x, dtype=float).copy()
    pad = win // 2
    xp = np.pad(np.asarray(x, dtype=float), pad, mode="edge")
    return np.median(np.lib.stride_tricks.sliding_window_view(xp, win), axis=-1)


def check_displacement_agreement(sync_fp: Path, t_frames, t_mts, d_mts,
                                 delta: float) -> dict:
    """Cross-check the alignment on the crosshead-travel channels.

    The sync CSV's displacement input (ai1) is the one channel on this batch
    that is genuinely connected, and it records the same crosshead motion the
    MTS does by a separate path. Fitting its arbitrary units onto the MTS
    millimetres and reporting the residual says whether the two records really
    are the same test. It does NOT independently confirm the time offset — a
    constant-rate ramp cannot — but a large residual is a loud signal that
    something is wrong with the pairing.
    """
    sync = pd.read_csv(sync_fp)
    col = next((c for c in sync.columns if "disp" in c.lower() and "scaled" in c.lower()), None)
    if col is None:
        return {"rmse_mm": np.nan, "n": 0, "span_mm": np.nan}
    ds = -pd.to_numeric(sync[col], errors="coerce").to_numpy()[:len(t_frames)]
    ds = _median_smooth(ds, SYNC_SMOOTH_FRAMES)
    d_at_frames = np.interp(t_frames + delta, t_mts, d_mts,
                            left=np.nan, right=np.nan)
    ok = np.isfinite(d_at_frames) & np.isfinite(ds)
    if ok.sum() < 30:
        return {"rmse_mm": np.nan, "n": int(ok.sum()), "span_mm": np.nan}
    scale, offset = np.polyfit(ds[ok], d_at_frames[ok], 1)
    resid = (ds[ok] * scale + offset) - d_at_frames[ok]
    return {"rmse_mm": float(np.sqrt(np.mean(resid ** 2))),
            "n": int(ok.sum()),
            "span_mm": float(np.nanmax(d_at_frames) - np.nanmin(d_at_frames))}


# =============================================================================
# HELPERS — fixture location from the deflected shape
# =============================================================================
def deflection_profile(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Reduce a frame's point cloud to a binned V(X) profile."""
    xb = np.round(df["X"].to_numpy() / PROFILE_BIN_MM) * PROFILE_BIN_MM
    prof = pd.DataFrame({"xb": xb, "V": df["V"].to_numpy()}).groupby("xb")["V"]
    counts = prof.size().to_numpy()
    means = prof.mean()
    x = means.index.to_numpy(dtype=float)
    v = means.to_numpy(dtype=float)
    keep = counts >= 5           # ignore sparse end bins
    return x[keep], v[keep]


def chord_at(x_query, x_l, v_l, x_r, v_r):
    """Value at x_query of the straight line through (x_l, v_l) and (x_r, v_r)."""
    if x_r == x_l:
        return np.nan
    t = (x_query - x_l) / (x_r - x_l)
    return v_l + t * (v_r - v_l)


def locate_fixture(df: pd.DataFrame, span_mm: float) -> dict:
    """Find midspan and the two supports in DIC X coordinates from one loaded frame.

    Midspan is the one place along a 3-point beam where dV/dx = 0, by symmetry,
    so it is the interior extremum of the deflected shape — and finding it that
    way needs no prior knowledge of the sign of V, of where the specimen sits in
    the world frame, or of how much of it the ROI covers.

    Note what does NOT work: taking the largest departure from a straight line
    fitted through the whole profile. For a symmetric bow the fit line sits at
    the mean, and the two ends are then twice as far off it as the middle is,
    so that lands on the ROI edge every time.

    After the interior extremum initialises it, two passes refine midspan as the
    largest departure from the chord through the two supports, which is the
    right criterion once you are already near midspan.
    """
    x, v = deflection_profile(df)
    if len(x) < 10:
        raise RuntimeError("deflected-shape profile too sparse to locate the fixture")

    # Pass 0 — interior extremum. Excluding a sixth of a span at each end keeps
    # the overhang tips, which travel furthest of all, out of the running.
    margin = span_mm / 6.0
    interior = (x > x.min() + margin) & (x < x.max() - margin)
    if interior.sum() < 5:
        interior = np.ones_like(x, dtype=bool)
    xi, vi = x[interior], v[interior]
    v_end = 0.5 * (v[0] + v[-1])
    i_lo, i_hi = int(np.argmin(vi)), int(np.argmax(vi))
    x_mid = float(xi[i_lo] if abs(vi[i_lo] - v_end) >= abs(vi[i_hi] - v_end)
                  else xi[i_hi])

    # Passes 1-2 — departure from the support chord.
    for _ in range(2):
        x_l, x_r = x_mid - span_mm / 2, x_mid + span_mm / 2
        v_l = np.interp(x_l, x, v)
        v_r = np.interp(x_r, x, v)
        dev = v - chord_at(x, x_l, v_l, x_r, v_r)
        near = np.abs(x - x_mid) < span_mm / 4       # midspan cannot wander far
        if not near.any():
            break
        idx = np.flatnonzero(near)
        x_mid = float(x[idx[np.argmax(np.abs(dev[near]))]])

    x_l, x_r = x_mid - span_mm / 2, x_mid + span_mm / 2
    v_l, v_r = np.interp(x_l, x, v), np.interp(x_r, x, v)
    # Sign convention: positive deflection = the way the beam actually moved.
    sign = 1.0 if (chord_at(x_mid, x_l, v_l, x_r, v_r)
                   - np.interp(x_mid, x, v)) > 0 else -1.0

    return {"x_mid": x_mid, "x_left": x_l, "x_right": x_r, "defl_sign": sign,
            "roi_x_min": float(df["X"].min()), "roi_x_max": float(df["X"].max())}


def window_mean_V(df: pd.DataFrame, x_c: float, half: float) -> tuple[float, float, int]:
    """Mean V and mean X of the points within +/- half of x_c. (nan, nan, 0) if empty."""
    m = np.abs(df["X"].to_numpy() - x_c) <= half
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    return (float(df["V"].to_numpy()[m].mean()),
            float(df["X"].to_numpy()[m].mean()),
            int(m.sum()))


def support_offset_factor(x_l_eff: float, x_r_eff: float, x_mid: float,
                          span_mm: float) -> float:
    """Correction for support windows sampled inboard of the supports themselves.

    A deflection referenced to the chord through two points that are Delta
    inboard of the supports under-reads, because the true deflection at those
    points is not zero. Using the ideal 3-point shape,
        w(u)/w_mid = u (3L^2 - 4u^2) / L^3      (u = distance from a support)
    the measured chord-referenced deflection is w_mid * (1 - mean of that at the
    two windows), so this returns its reciprocal. Returns 1.0 if the correction
    is switched off or looks unphysical.
    """
    if not APPLY_SUPPORT_OFFSET_CORRECTION:
        return 1.0
    L = span_mm
    u_l = max(0.0, x_l_eff - (x_mid - L / 2))
    u_r = max(0.0, (x_mid + L / 2) - x_r_eff)
    shape = lambda u: u * (3 * L**2 - 4 * u**2) / L**3
    lost = 0.5 * (shape(u_l) + shape(u_r))
    if not np.isfinite(lost) or not (0.0 <= lost < 0.25):
        return 1.0
    return 1.0 / (1.0 - lost)


# =============================================================================
# HELPERS — the through-depth strain fit
# =============================================================================
def midspan_strain_fit(df: pd.DataFrame, x_mid: float,
                       y_lo: float, y_hi: float) -> dict:
    """Straight-line fit of exx against Y over the depth at midspan.

    Euler-Bernoulli says exx varies linearly through the depth; the slope is
    the curvature (up to sign), the zero crossing is the neutral axis, and the
    value at mid-depth is the membrane (pure axial) part, which should be ~0 in
    clean bending. The fit's R^2 is the test of whether any of that holds.
    """
    x = df["X"].to_numpy()
    y = df["Y"].to_numpy()
    e = df["exx"].to_numpy()
    m = (np.abs(x - x_mid) <= MIDSPAN_HALF_WIDTH_MM) & np.isfinite(e) & np.isfinite(y)
    m &= (y >= y_lo + PROFILE_EDGE_TRIM_MM) & (y <= y_hi - PROFILE_EDGE_TRIM_MM)
    if m.sum() < 20:
        return {"slope": np.nan, "icept": np.nan, "r2": np.nan, "n": int(m.sum())}

    yy, ee = y[m], e[m]
    slope, icept = np.polyfit(yy, ee, 1)
    resid = ee - (slope * yy + icept)
    ss_tot = float(np.sum((ee - ee.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else np.nan
    return {"slope": float(slope), "icept": float(icept), "r2": r2, "n": int(m.sum())}


# =============================================================================
# STEP A — .out → per-frame CSV
# =============================================================================
def export_out_to_csv(out_path: Path, csv_path: Path,
                      var_names: list[str]) -> bool:
    """Convert a single .out to CSV. Byte-for-byte the same routine as
    TensileDIC_Level1.export_out_to_csv, so the two test types produce
    interchangeable per-frame CSVs.

    Writes one row per valid AOI point, columns named after var_names, with
    sigma used to filter invalid points and then dropped from the output.
    Returns True on success.
    """
    ds = VICDataSet()
    try:
        ds.load(str(out_path))
    except Exception as ex:
        print(f"    [!] could not load {out_path.name}: {ex}")
        return False

    try:
        available = list(ds.variables())
    except Exception:
        available = []
    wanted = [v for v in var_names if (not available) or (v in available)]
    if not wanted:
        print(f"    [!] none of the requested variables found in {out_path.name}")
        return False
    try:
        values = ds.get_values(wanted)
    except Exception as ex:
        print(f"    [!] get_values failed on {out_path.name}: {ex}")
        return False

    mask = np.ones(len(values), dtype=bool)
    if "sigma" in values.dtype.names:
        mask &= values["sigma"] >= 0

    export_cols = [v for v in wanted if v != "sigma"]
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write(",".join(export_cols) + "\n")
        arrs = [np.asarray(values[v])[mask] for v in export_cols]
        for row in zip(*arrs):
            fh.write(",".join(f"{x}" for x in row) + "\n")
    return True


def export_frames(cid: str, cdir: Path) -> None:
    """Step A for one coupon: export every .out to a sibling CSV.

    Costs ~1 GB per coupon and Step B does not read these — it loads each .out
    itself. They exist so the flexural raw data is laid out exactly like the
    tensile raw data, and so anything that wants a frame (the heatmap animation,
    an ad-hoc look at one frame) can read a CSV instead of needing vicpyx.
    """
    out_files = find_out_files(cdir)
    if not out_files:
        print(f"  Step A: [skip] no .out files in {cdir}")
        return
    n_written = n_skipped = 0
    for fp in out_files:
        csv_fp = fp.with_suffix(".csv")
        if csv_fp.exists() and not OVERWRITE_FRAMES:
            n_skipped += 1
            continue
        if export_out_to_csv(fp, csv_fp, EXPORT_VARS):
            n_written += 1
    print(f"  Step A: {len(out_files)} .out files, wrote {n_written}, "
          f"skipped {n_skipped} (already existed)")


# =============================================================================
# STEP B — reduce frames + pair with MTS, write the per-coupon CSV
# =============================================================================
def build_l1(cid: str, cdir: Path) -> dict | None:
    """Step B for one coupon. Returns its coupon_scalars row, or None if skipped."""
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    out_fp = DIC_DIR / f"{cid}.csv"
    if out_fp.exists() and not OVERWRITE_L1:
        print(f"  Step B: [skip] {out_fp.name} exists")
        return None

    out_files = find_out_files(cdir)
    sync_fp = find_sync_csv(cdir)
    mts_fp = find_mts_txt(cid)
    if not out_files:  print(f"  Step B: [skip] no .out files in {cdir}");   return None
    if sync_fp is None: print(f"  Step B: [skip] no sync CSV in {cdir}");    return None
    if mts_fp is None:  print(f"  Step B: [skip] no MTS .txt for {cid}");    return None

    b_mm, d_mm = specimen_geometry(cid)
    print(f"  Step B: {len(out_files)} frames | sync {sync_fp.name} | MTS {mts_fp.name}")
    print(f"    geometry: b = {b_mm:.2f} mm, d = {d_mm:.3f} mm "
          f"(L/d = {FLEX_SPAN_MM / d_mm:.1f})")

    # ---- is this project even in the expected coordinate frame? -----------
    # Checked before anything expensive runs, and treated as fatal for the
    # coupon: everything below assumes X is the span and Y is the depth, and
    # nothing downstream can notice when that is wrong.
    ref = read_frame(out_files[0])
    if ref is None:
        raise RuntimeError("reference frame has no correlated points")
    complaint = check_roi_orientation(ref, d_mm)
    if complaint is not None:
        print(f"    [!] SKIPPED — {complaint}.")
        print(f"        This script needs X along the span and Y through the depth. "
              f"Re-align this\n        VIC-3D project's coordinate system and re-export, "
              f"then re-run. Writing nothing:\n        a missing coupon is recoverable, "
              f"a silently wrong one is not.")
        return None

    # ---- the frame clock, and the 5 Hz check ------------------------------
    t_frames, rate_hz = frame_clock(sync_fp, len(out_files))
    if not np.isfinite(rate_hz) or abs(rate_hz - EXPECTED_FRAME_RATE_HZ) > FRAME_RATE_TOL_HZ:
        print(f"    [!] frame rate {rate_hz:.3f} Hz, expected "
              f"{EXPECTED_FRAME_RATE_HZ:g} Hz — check the DIC acquisition settings")
    else:
        print(f"    frame rate: {rate_hz:.3f} Hz over {t_frames[-1]:.1f} s")

    # ---- MTS: flip to first quadrant, zero displacement, strip the tare ----
    mts = load_mts_txt(mts_fp)
    # Bending runs in compression: crosshead travels down, force reads negative.
    d_mts = -mts["disp_mm"].to_numpy()
    f_mts = -mts["force_N"].to_numpy()
    t_mts = mts["time_s"].to_numpy()
    d_mts = d_mts - d_mts[0]
    i_pk = int(np.argmax(f_mts))
    tare = find_force_baseline(d_mts[:i_pk + 1], f_mts[:i_pk + 1])
    f_mts = f_mts - tare
    peak_net_N = float(f_mts[i_pk])
    print(f"    MTS: tare {tare:.0f} N removed | net peak {peak_net_N:.0f} N "
          f"at t = {t_mts[i_pk]:.1f} s on the machine clock")

    # ---- align on fracture, read off the images ---------------------------
    i_break = find_break_frame(out_files)
    delta = align_mts_to_frames(t_frames, i_break, t_mts, f_mts)
    dchk = check_displacement_agreement(sync_fp, t_frames, t_mts, d_mts, delta)
    print(f"    break: correlation lost after frame {i_break} of {len(out_files)} "
          f"(t = {t_frames[i_break]:.1f} s)")
    print(f"    align: MTS clock leads the frames by {delta:.2f} s | "
          f"crosshead cross-check {dchk['rmse_mm']:.2f} mm RMS over "
          f"{dchk['span_mm']:.1f} mm of travel")

    t_query = t_frames + delta
    force_N = np.interp(t_query, t_mts, f_mts, left=np.nan, right=np.nan)
    disp_mts_mm = np.interp(t_query, t_mts, d_mts, left=np.nan, right=np.nan)
    n_off = int(np.sum(~np.isfinite(force_N)))
    if n_off:
        print(f"    [i] {n_off} frame(s) fall outside the MTS record and carry "
              f"no force — Level 2 truncates them out")

    # ---- locate the fixture on a well-loaded, undamaged frame -------------
    target = FIXTURE_PROBE_LOAD_FRAC * np.nanmax(force_N)
    cand = np.flatnonzero(np.isfinite(force_N) & (force_N >= target))
    probe_i = min(int(cand[0]) if cand.size else len(out_files) // 2,
                  len(out_files) - 1)
    probe = None
    while probe is None and probe_i > 0:
        probe = read_frame(out_files[probe_i])
        if probe is None:
            probe_i -= 5
    if probe is None:
        raise RuntimeError("no loaded frame with correlation, cannot locate the fixture")

    fix = (locate_fixture(probe, FLEX_SPAN_MM) if MIDSPAN_X_MM is None else
           {"x_mid": MIDSPAN_X_MM,
            "x_left": MIDSPAN_X_MM - FLEX_SPAN_MM / 2,
            "x_right": MIDSPAN_X_MM + FLEX_SPAN_MM / 2,
            "defl_sign": 1.0,
            "roi_x_min": float(probe["X"].min()),
            "roi_x_max": float(probe["X"].max())})
    x_mid = fix["x_mid"]
    print(f"    fixture (frame {probe_i}): midspan X = {x_mid:.1f} mm, supports X = "
          f"{fix['x_left']:.1f} / {fix['x_right']:.1f} mm, ROI X = "
          f"{fix['roi_x_min']:.1f} .. {fix['roi_x_max']:.1f} mm")

    _, x_l_eff, n_l = window_mean_V(probe, fix["x_left"], SUPPORT_HALF_WIDTH_MM)
    _, x_r_eff, n_r = window_mean_V(probe, fix["x_right"], SUPPORT_HALF_WIDTH_MM)
    corr = support_offset_factor(x_l_eff, x_r_eff, x_mid, FLEX_SPAN_MM)
    print(f"    deflection windows: left {n_l} pts at {x_l_eff:.1f} mm, right {n_r} pts "
          f"at {x_r_eff:.1f} mm -> D scaled by {corr:.4f}")
    for side, npts in (("left", n_l), ("right", n_r)):
        if npts < 20:
            print(f"    [!] only {npts} points in the {side} support window — "
                  f"midspan deflection is weakly referenced on that side")

    # ROI depth coverage. The subset radius insets the ROI from both faces, so
    # the surfaces D790 cares about are never measured directly.
    y_lo, y_hi = float(ref["Y"].min()), float(ref["Y"].max())
    y_c = 0.5 * (y_lo + y_hi)
    y_bot, y_top = y_c - d_mm / 2, y_c + d_mm / 2   # assumed specimen faces
    print(f"    depth: ROI Y = {y_lo:.2f} .. {y_hi:.2f} mm ({y_hi - y_lo:.2f} mm of "
          f"d = {d_mm:.2f}, {100 * (y_hi - y_lo) / d_mm:.0f}% covered); faces assumed "
          f"at Y = {y_bot:.2f} / {y_top:.2f} mm")

    # ---- per-frame pass ---------------------------------------------------
    n = min(len(out_files), len(t_frames))
    rows = []
    t0 = time.time()
    for i in range(n):
        df = read_frame(out_files[i])
        rec = {"step": i, "time_s": t_frames[i], "force_N": force_N[i],
               "disp_mts_mm": disp_mts_mm[i]}
        if df is None:
            rec.update({"n_pts": 0, "defl_mm": np.nan, "kappa_1pmm": np.nan,
                        "na_Y_mm": np.nan, "profile_r2": np.nan,
                        "eps_bot": np.nan, "eps_top": np.nan,
                        "eps_membrane": np.nan})
            rows.append(rec)
            continue

        # --- midspan deflection, referenced to the support chord ---
        v_mid, _, n_mid = window_mean_V(df, x_mid, SUPPORT_HALF_WIDTH_MM)
        v_l, x_l_i, n_l_i = window_mean_V(df, fix["x_left"], SUPPORT_HALF_WIDTH_MM)
        v_r, x_r_i, n_r_i = window_mean_V(df, fix["x_right"], SUPPORT_HALF_WIDTH_MM)
        if min(n_mid, n_l_i, n_r_i) >= 3:
            chord = chord_at(x_mid, x_l_i, v_l, x_r_i, v_r)
            defl = (fix["defl_sign"] * (chord - v_mid)
                    * support_offset_factor(x_l_i, x_r_i, x_mid, FLEX_SPAN_MM))
        else:
            defl = np.nan

        # --- curvature and neutral axis from exx(Y) at midspan ---
        fit = midspan_strain_fit(df, x_mid, y_lo, y_hi)
        slope = fit["slope"]
        if np.isfinite(slope) and slope != 0:
            kappa = abs(slope)                       # 1/mm
            na_y = -fit["icept"] / slope
            e_bot = slope * y_bot + fit["icept"]
            e_top = slope * y_top + fit["icept"]
            e_mem = slope * y_c + fit["icept"]
        else:
            kappa = na_y = e_bot = e_top = e_mem = np.nan

        rec.update({"n_pts": len(df), "defl_mm": defl, "kappa_1pmm": kappa,
                    "na_Y_mm": na_y, "profile_r2": fit["r2"], "eps_bot": e_bot,
                    "eps_top": e_top, "eps_membrane": e_mem})
        rows.append(rec)

        if (i + 1) % 100 == 0 or i == n - 1:
            el = time.time() - t0
            print(f"      [{i+1}/{n}] {el:.0f} s elapsed, "
                  f"ETA {el / (i + 1) * (n - i - 1):.0f} s")

    out = pd.DataFrame(rows)
    out.to_csv(out_fp, index=False, float_format="%.6g")
    n_bad = int((out["n_pts"] == 0).sum())
    print(f"  → DIC/{out_fp.name} ({len(out)} rows, {n_bad} uncorrelated/post-fracture)")

    return {"coupon": cid, "b_mm": b_mm, "d_mm": d_mm,
            "x_mid_mm": x_mid, "x_left_mm": fix["x_left"],
            "x_right_mm": fix["x_right"], "defl_sign": fix["defl_sign"],
            "y_roi_lo_mm": y_lo, "y_roi_hi_mm": y_hi,
            "y_bot_mm": y_bot, "y_top_mm": y_top,
            "tare_mts_N": tare, "peak_net_N": peak_net_N,
            "frame_rate_hz": rate_hz, "mts_offset_s": delta,
            "break_frame": i_break, "fixture_probe_frame": probe_i,
            "disp_check_rmse_mm": dchk["rmse_mm"],
            "n_frames": len(out), "n_uncorrelated": n_bad}


# =============================================================================
# STEP C — what is actually in the .out files
# =============================================================================
def list_vars(cid: str, cdir: Path) -> None:
    """Print the variable set of this coupon's first .out.

    Run after reprocessing a project in VIC-3D: any new inspector item that is
    genuinely being written to the export shows up here. Anything listed here
    can be added to OPTIONAL_VARS and it will be carried through Step B.
    """
    out_files = find_out_files(cdir)
    if not out_files:
        print(f"  Step C: [skip] no .out files in {cdir}")
        return
    have = frame_variables(out_files[0])
    missing = [v for v in REQUIRED_VARS if v not in have]
    unused = [v for v in have if v not in REQUIRED_VARS + OPTIONAL_VARS]
    print(f"  Step C: {len(have)} variable(s) in {out_files[0].name}")
    print(f"    {', '.join(have)}")
    if missing:
        print(f"    [!] REQUIRED but absent: {', '.join(missing)} — Step B cannot run")
    if unused:
        print(f"    [i] present but not read: {', '.join(unused)} "
              f"(add to OPTIONAL_VARS to carry them through)")


# =============================================================================
# MAIN
# =============================================================================
def process_coupon(cid: str) -> dict | None:
    print(f"[{cid}]")
    cdir = raw_folder(cid)
    if not cdir.is_dir():
        print(f"  [skip] directory not found: {cdir}")
        return None
    if DO_LIST_VARS:
        list_vars(cid, cdir)
    if DO_EXPORT_FRAMES:
        export_frames(cid, cdir)
    if DO_BUILD_L1:
        return build_l1(cid, cdir)
    return None


def write_coupon_scalars(rows: list[dict]) -> Path:
    """Upsert this run's per-coupon scalars into coupon_scalars.csv — the single
    place these are stored, so they are never repeated down every row of a
    per-frame CSV. Coupons skipped this run keep whatever row is already on disk
    from a previous run.

    The same file TensileDIC_Level1 writes. The merge is keyed on coupon ID and
    the two test types never share an ID, so tensile rows pass through this
    function untouched and keep their own columns; pandas fills the columns each
    test type does not have. Identical routine to
    TensileDIC_Level1.write_coupon_scalars — either script can be run first.
    """
    out_fp = DIC_DIR / "coupon_scalars.csv"
    DIC_DIR.mkdir(parents=True, exist_ok=True)
    merged = {}
    if out_fp.exists():
        for row in pd.read_csv(out_fp).to_dict("records"):
            merged[row["coupon"]] = row
    for row in rows:
        merged[row["coupon"]] = row
    pd.DataFrame(sorted(merged.values(), key=lambda r: r["coupon"])
                ).to_csv(out_fp, index=False, float_format="%.6g")
    print(f"Coupon scalars: {out_fp}")
    return out_fp


def main() -> None:
    t0 = time.time()
    print("=" * 74)
    print("FlexuralDIC_Level1 — reduce .out frames to bending kinematics + pair with MTS")
    print("=" * 74)
    print(f"Raw root  : {RAW_ROOT}")
    print(f"DIC dir   : {DIC_DIR}")
    print(f"Span      : {FLEX_SPAN_MM:.1f} mm ({FLEX_SPAN_MM / IN2MM:.2f} in) — confirmed fixture setting")
    print()

    coupons = selected_coupons()
    print(f"Processing {len(coupons)} coupon(s)\n")

    scalar_rows = []
    for i, cid in enumerate(coupons, start=1):
        t_c = time.time()
        try:
            row = process_coupon(cid)
            if row is not None:
                scalar_rows.append(row)
        except Exception as ex:
            print(f"  [error] {cid}: {ex}")
        elapsed = time.time() - t0
        print(f"  [{i}/{len(coupons)}] {time.time() - t_c:.1f} s this coupon | "
              f"elapsed {elapsed:.1f} s | ETA {elapsed / i * (len(coupons) - i):.1f} s\n")

    if scalar_rows:
        write_coupon_scalars(scalar_rows)

    print(f"Done. {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
