#!/usr/bin/env python3
"""
dic_heatmap_animation.py  —  FSR Tensile + Flexural Coupons
============================================================
Renders a full-field DIC animation for one or more coupons, sampling every
FRAME_STRIDE-th frame. One script covers both test types: they share every hard
part (striding frames, holding one colour scale across a test, writing an MP4)
and differ only in where the frames live, what the axes mean, and what is drawn
on top. Those differences live in TEST_TYPES, frame_set() and the two overlay
classes.

Which field is coloured is COLOR_VAR. Left on "auto" it picks eyy for tensile
(axial strain, the loading direction) and exx for flexural (the ROI is the side
profile, so exx is the through-depth bending strain D790 is written about).

USAGE
    python dic_heatmap_animation.py <coupon_id> [<coupon_id> ...] [--var exx]

    e.g.  python dic_heatmap_animation.py P01-TSW00-01
          python dic_heatmap_animation.py P01-FCL00-01 P01-FCL90-01
          python dic_heatmap_animation.py P01-FCL00-01 --var eyy

    --var overrides COLOR_VAR for this run only.

INPUT per coupon
  <coupon_dir>/<stem>-????????_0.csv   Level-1 Step A per-frame CSVs
  <coupon_dir>/*.out                   fallback when Step A has not been run
                                       (needs vicpyx; slower)
  <DIC_DIR>/coupon_scalars.csv         flexural fixture, when Level 1 wrote a row

OUTPUT
  <FIGS_ROOT>/<coupon_id>/<var>_field_animation.mp4   (or .gif)
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
import matplotlib.animation as animation

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# SWITCHES
# =============================================================================
# "auto" resolves per test type through AUTO_COLOR_VAR. Anything else is used
# for both types as given: exx, eyy, exy, U, V, W (tensile CSVs also carry
# e1, e2, gamma). A field the frames lack is reported and the coupon skipped.
COLOR_VAR = "auto"
AUTO_COLOR_VAR = {"tensile": "eyy", "flexural": "exx"}

FRAME_STRIDE  = 5    # use every Nth frame
FPS           = 5    # frames per second in the output video
N_CBAR_SAMPLE = 10   # frames sampled when range_from is "sampled"

# Flexural view. The ROI is ~180 mm along the span by ~8 mm through the depth;
# at equal aspect the through-depth gradient — the point of the plot — is
# invisible, so Y is stretched.
#   "fill"    stretch to fill the figure and report the factor on the axis
#             label. Not a fixed number because the deformed view's Y extent
#             grows with the deflection (~8 mm of depth becomes ~28 mm).
#   <number>  a fixed Y:X stretch, for comparing two coupons at one scale.
#   "equal"   true shape. Only readable with FLEX_X_WINDOW_MM cropped to midspan.
FLEX_ASPECT      = "fill"
FLEX_X_WINDOW_MM = None    # crop to +/- this about midspan; None = whole ROI

# =============================================================================
# PATHS
# =============================================================================
FIGS_ROOT = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\figs"
)
DIC_DIR = FIGS_ROOT.parent / "DIC"
RAW_ROOTS = {
    "tensile": {
        "CL": DIC_DIR / "raw" / "2026_FSR_TensileTest_TCL",
        "SW": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
        "UV": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
        "IS": DIC_DIR / "raw" / "2026_FSR_TensileTest_TSW_TIS_TUV",
    },
    "flexural": DIC_DIR / "raw" / "2026_FSR_Flexural_FCL_FIS",
}
COUPON_SCALARS = DIC_DIR / "coupon_scalars.csv"

# Tensile point extensometers — must match TensileDIC_Level1.py.
IN2MM          = 25.4
AXIAL_GAUGE_IN = 4.36    # E0 axial
TRANS_GAUGE_IN = 1.0     # E1 transverse

# =============================================================================
# PER-TEST-TYPE RENDER SETTINGS
#
# range_from — which frames set the colour range: "sampled" pools
# N_CBAR_SAMPLE frames spread across the test, "last" uses the last correlated
# frame alone. The tensile gauge section strains roughly uniformly, so pooling
# describes it well. The flexural field is neither uniform nor stationary —
# moment is triangular along the span and the extreme fibre grows ~8x over the
# test — so pooling hands the scale to the early frames where nothing is
# happening (+/-0.0023 against a real midspan strain of +/-0.019, i.e. every
# later frame saturates).
# =============================================================================
TEST_TYPES = {
    "tensile": {
        "cmap":       "plasma",
        "symmetric":  False,        # one-signed field
        "pctile":     (10, 90),
        "range_from": "sampled",
        "deformed":   False,
        "figsize":    (5.0, None),  # None height = derived from the ROI aspect
        "marker_s":   6,
    },
    "flexural": {
        # Diverging, because the field is signed about the neutral axis — but it
        # needs a DARK centre, not coolwarm/RdBu/bwr. Their near-white zero makes
        # the low-moment ends of the beam the brightest thing on a #1a1a1a figure,
        # competing with the white fixture overlay. berlin's centre is near-black,
        # so zero falls into the background and the fibres glow out of it.
        # Needs matplotlib >= 3.10; older versions fall back to "coolwarm".
        "cmap":       "berlin",
        "symmetric":  True,         # so the neutral axis is the midpoint colour
        "pctile":     (2, 98),
        "range_from": "last",
        # Plot each point at X+U, Y+V so the beam bends as it colours. VIC-3D's
        # X and Y are reference coordinates, so the default view is a rectangle
        # that only changes colour — which throws away the largest signal in the
        # test, since midspan travels ~12 mm (about a full specimen depth).
        "deformed":   True,
        "figsize":    (11.0, 3.4),
        "marker_s":   5,
    },
}

FIELD_LABELS = {
    ("tensile",  "eyy"): "$\\epsilon_{yy}$  axial strain [1]",
    ("tensile",  "exx"): "$\\epsilon_{xx}$  transverse strain [1]",
    ("flexural", "exx"): "$\\epsilon_{xx}$  bending strain [1]",
    ("flexural", "eyy"): "$\\epsilon_{yy}$  through-depth strain [1]",
}
GENERIC_LABELS = {
    "exy":   "$\\epsilon_{xy}$  shear strain [1]",
    "e1":    "$\\epsilon_1$  major principal strain [1]",
    "e2":    "$\\epsilon_2$  minor principal strain [1]",
    "gamma": "$\\gamma$  max shear strain [1]",
    "U":     "U  displacement X (mm)",
    "V":     "V  displacement Y (mm)",
    "W":     "W  displacement Z (mm)",
}

BG = "#1a1a1a"
KIND_OF = {"T": "tensile", "F": "flexural"}


# =============================================================================
# COUPON ID AND LABELS
# =============================================================================
def parse_id(cid):
    """P01-FCL00-01 -> kind flexural, exposure CL, direction 00, rep 01."""
    parts = cid.split("-")
    if len(parts) != 3:
        raise ValueError(f"'{cid}' is not a coupon ID like P01-TSW00-01")
    group = parts[1]
    kind = KIND_OF.get(group[0])
    if kind is None:
        raise ValueError(f"'{cid}': test type '{group[0]}' is not tensile (T) or "
                         f"flexural (F) — bearing coupons have no DIC field here")
    return {"print": parts[0], "kind": kind, "exposure": group[1:-2],
            "direction": group[-2:], "rep": parts[2], "group": group}


def field_label(kind, var):
    return FIELD_LABELS.get((kind, var), GENERIC_LABELS.get(var, f"{var} [1]"))


# =============================================================================
# FRAME SOURCES
# Both paths return a list of frame files and a reader that turns one into a
# DataFrame of valid points (or None). Everything downstream is blind to which.
# =============================================================================
def flex_level1():
    """Import FlexuralDIC_Level1 for its .out reader and bending helpers.
    Lazy, and only for flexural coupons: it pulls in vicpyx at module level and
    a tensile run has no business requiring that."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import FlexuralDIC_Level1 as L1
    except ImportError as ex:
        raise RuntimeError(
            f"cannot import FlexuralDIC_Level1 ({ex}). Flexural frames are read "
            f"from .out files through vicpyx — pip install vicpyx") from ex
    return L1


def read_csv_frame(fp):
    """One per-frame CSV, valid points only."""
    try:
        df = pd.read_csv(fp)
    except Exception:
        return None
    df = df.dropna(subset=["X", "Y"])
    return df if len(df) else None


def flex_scalars(cid):
    """This coupon's row of coupon_scalars.csv, or None. A row with no fixture
    in it (a tensile row, or a flexural coupon written before the fixture
    columns existed) counts as absent."""
    if not COUPON_SCALARS.exists():
        return None
    try:
        g = pd.read_csv(COUPON_SCALARS)
    except Exception:
        return None
    row = g.loc[g["coupon"] == cid]
    if row.empty:
        return None
    rec = row.iloc[0].to_dict()
    return rec if np.isfinite(rec.get("x_mid_mm", np.nan)) else None


def coupon_raw_dir(cid, meta):
    """The VIC-3D project folder holding this coupon's frames. A tensile folder
    is the coupon ID; the flexural projects drop the print prefix and the dashes
    (P01-FCL00-01 -> FCL0001), same rule as FlexuralDIC_Level1.raw_folder."""
    if meta["kind"] == "tensile":
        root = RAW_ROOTS["tensile"].get(meta["exposure"])
        if root is None:
            raise RuntimeError(f"unknown tensile exposure '{meta['exposure']}' "
                               f"(known: {sorted(RAW_ROOTS['tensile'])})")
        return root / cid
    return RAW_ROOTS["flexural"] / f"{meta['group']}{meta['rep']}"


def find_break_frame(paths, reader):
    """Index of the last frame that still correlated — i.e. fracture. Scans
    backwards rather than bisecting: correlation is not monotonic in the frames
    just before failure."""
    for i in range(len(paths) - 1, -1, -1):
        df = reader(paths[i])
        if df is not None and len(df) >= 200:
            return i
    raise RuntimeError("no frame in this coupon correlated at all")


def frame_set(cid, meta):
    """Locate a coupon's frames and hand back a reader for them.

    Step A exports a per-frame CSV next to every .out for both test types, so
    the CSVs are used when present and the .out files only when they are not.
    Flexural paths are truncated at fracture: a bend specimen decorrelates the
    instant it cracks, and the blank frames after it would otherwise be what
    the colour range is taken from.
    """
    cdir = coupon_raw_dir(cid, meta)
    if not cdir.is_dir():
        raise RuntimeError(f"{cdir} not found")

    # Frame files are "<stem>-00000000_0"; the sync CSV alongside them is just
    # "<stem>.csv" and must not be picked up as a frame.
    stem = cid if meta["kind"] == "tensile" else cdir.name
    paths = sorted(cdir.glob(f"{stem}-????????_0.csv"))
    reader, note = read_csv_frame, "per-frame CSVs"
    if not paths:
        paths = sorted(cdir.rglob("*.out"))
        if not paths:
            raise RuntimeError(f"no per-frame CSVs or .out files in {cdir}")
        reader, note = flex_level1().read_frame, ".out files (Step A not run)"

    n_total = len(paths)
    if meta["kind"] == "flexural":
        scal = flex_scalars(cid)
        if scal is not None and np.isfinite(scal.get("break_frame", np.nan)):
            paths = paths[:int(scal["break_frame"]) + 1]
        else:
            try:
                paths = paths[:find_break_frame(paths, reader) + 1]
            except RuntimeError:
                pass      # nothing correlated; let the caller report it

    return {"paths": paths, "read": reader, "n_total": n_total,
            "source": f"{n_total} {note} in {cdir}"}


def stride_frames(paths, stride):
    """Every stride-th frame, always including the last."""
    idx = list(range(0, len(paths), stride))
    if idx[-1] != len(paths) - 1:
        idx.append(len(paths) - 1)
    return idx, [paths[i] for i in idx]


# =============================================================================
# COLOUR RANGE
# =============================================================================
def frame_values(fs, i, var):
    """The finite values of `var` in frame i, or None if there are none."""
    df = fs["read"](fs["paths"][i])
    if df is None or var not in df.columns:
        return None
    v = pd.to_numeric(df[var], errors="coerce").to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return v if v.size else None


def compute_color_range(fs, var, cfg):
    """The colour range, held fixed for the whole animation — a scale that moved
    frame to frame would make every frame look alike and hide the growth of the
    field, which is the thing being animated. Which frames set it is per test
    type; see TEST_TYPES."""
    lo_p, hi_p = cfg["pctile"]
    n = len(fs["paths"])
    vals = []

    if cfg["range_from"] == "last":
        # Walk back from the end: the last path is the last correlated frame for
        # flexural, but a tensile coupon read this way may end on a bad frame.
        for i in range(n - 1, max(n - 21, -1), -1):
            v = frame_values(fs, i, var)
            if v is not None:
                vals.append(v)
                break
    else:
        for i in np.round(np.linspace(0, n - 1, min(N_CBAR_SAMPLE, n))).astype(int):
            v = frame_values(fs, i, var)
            if v is not None:
                vals.append(v)

    if not vals:
        return 0.0, 0.02

    combined = np.concatenate(vals)
    vmin = float(np.percentile(combined, lo_p))
    vmax = float(np.percentile(combined, hi_p))
    if cfg["symmetric"]:
        vmax = max(abs(vmin), abs(vmax))
        vmin = -vmax
    if vmin >= vmax:
        vmax = vmin + 1e-6
    return vmin, vmax


# =============================================================================
# TENSILE OVERLAY — the two virtual point extensometers
# (mirrors TensileDIC_Level1.ext_endpoints / point_extensometer)
# =============================================================================
def ext_endpoints(ref_df, axial_mm, trans_mm):
    ref = ref_df.dropna(subset=["X", "Y"])
    Yc, Xc = float(ref["Y"].median()), float(ref["X"].median())
    Y_top = min(Yc + axial_mm / 2, float(ref["Y"].max()))
    Y_bot = max(Yc - axial_mm / 2, float(ref["Y"].min()))
    X_rgt = min(Xc + trans_mm / 2, float(ref["X"].max()))
    X_lft = max(Xc - trans_mm / 2, float(ref["X"].min()))
    return Xc, Y_bot, Y_top, X_lft, X_rgt, Yc


def point_ext_strain(df, x0, y0, x1, y1):
    if "U" not in df.columns or "V" not in df.columns:
        return np.nan
    L0 = float(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
    if L0 == 0:
        return np.nan
    r0 = df.loc[((df["X"] - x0) ** 2 + (df["Y"] - y0) ** 2).idxmin()]
    r1 = df.loc[((df["X"] - x1) ** 2 + (df["Y"] - y1) ** 2).idxmin()]
    dx = (x1 + float(r1["U"])) - (x0 + float(r0["U"]))
    dy = (y1 + float(r1["V"])) - (y0 + float(r0["V"]))
    return (np.sqrt(dx ** 2 + dy ** 2) - L0) / L0


class TensileOverlay:
    """E0 axial and E1 transverse gauges, redrawn with their live strain."""

    def __init__(self, ax, ref):
        self.ax = ax
        (self.Xc, self.Y_bot, self.Y_top,
         self.X_lft, self.X_rgt, self.Yc) = ext_endpoints(
            ref, AXIAL_GAUGE_IN * IN2MM, TRANS_GAUGE_IN * IN2MM)
        self.artists = []

    def describe(self):
        return [f"E0 axial      : Xc = {self.Xc:.1f} mm, Y {self.Y_bot:.1f} -> "
                f"{self.Y_top:.1f} mm  (L = {self.Y_top - self.Y_bot:.1f} mm)",
                f"E1 transverse : Yc = {self.Yc:.1f} mm, X {self.X_lft:.1f} -> "
                f"{self.X_rgt:.1f} mm  (L = {self.X_rgt - self.X_lft:.1f} mm)"]

    def gauge(self, x0, y0, x1, y1, label, eps):
        ax = self.ax
        line, = ax.plot([x0, x1], [y0, y1], color="white", lw=1.0, zorder=6)
        mk0, = ax.plot(x0, y0, "s", color="white", ms=5, zorder=7)
        mk1, = ax.plot(x1, y1, "s", color="white", ms=5, zorder=7)
        val = f"{eps:.5e}" if np.isfinite(eps) else "—"
        ann = ax.annotate(
            f"{label}: {val} [1]",
            xy=((x0 + x1) / 2, (y0 + y1) / 2), xytext=(4, 4),
            textcoords="offset points", fontsize=5.5, color="white", zorder=8,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#555555", alpha=0.80, lw=0),
        )
        return [line, mk0, mk1, ann]

    def update(self, df):
        for a in self.artists:
            a.remove()
        self.artists = []
        e0 = point_ext_strain(df, self.Xc, self.Y_bot, self.Xc, self.Y_top)
        e1 = point_ext_strain(df, self.X_lft, self.Yc, self.X_rgt, self.Yc)
        self.artists += self.gauge(self.Xc, self.Y_bot, self.Xc, self.Y_top, "E0", e0)
        self.artists += self.gauge(self.X_lft, self.Yc, self.X_rgt, self.Yc, "E1", e1)


# =============================================================================
# FLEXURAL FIXTURE
# =============================================================================
def probe_frame(fs):
    """A well-loaded, undamaged frame: ~60 % of the way to the last frame."""
    i = int(0.6 * (len(fs["paths"]) - 1))
    while i > 0:
        df = fs["read"](fs["paths"][i])
        if df is not None:
            return df
        i -= 5
    raise RuntimeError("no correlated frame to locate the fixture on")


def resolve_fixture(geom, fs, L1):
    """Midspan, the supports and the specimen faces in DIC coordinates.

    Level 1's row when there is one, so the animation and the reduced CSV share
    one fixture rather than two independent guesses at it. Without a row it is
    located here from the deflected shape, the same way Level 1 would.
    Resolved once and shared: the overlay and the deformed-view transform must
    not disagree about where midspan is.
    """
    if geom is not None:
        fix = {"x_mid": float(geom["x_mid_mm"]),
               "x_left": float(geom["x_left_mm"]),
               "x_right": float(geom["x_right_mm"]),
               "sign": float(geom["defl_sign"]),
               "y_lo": float(geom["y_roi_lo_mm"]),
               "y_hi": float(geom["y_roi_hi_mm"]),
               "y_bot": float(geom["y_bot_mm"]),
               "y_top": float(geom["y_top_mm"])}
        fix["note"] = (f"fixture from DIC/{COUPON_SCALARS.name}: midspan X = "
                       f"{fix['x_mid']:.1f} mm, supports {fix['x_left']:.1f} / "
                       f"{fix['x_right']:.1f} mm")
        return fix

    probe = probe_frame(fs)
    f = L1.locate_fixture(probe, L1.FLEX_SPAN_MM)
    y_lo, y_hi = float(probe["Y"].min()), float(probe["Y"].max())
    # With no geometry row the faces are taken as the ROI edges rather than
    # guessed at from a nominal depth.
    return {"x_mid": f["x_mid"], "x_left": f["x_left"], "x_right": f["x_right"],
            "sign": f["defl_sign"], "y_lo": y_lo, "y_hi": y_hi,
            "y_bot": y_lo, "y_top": y_hi,
            "note": (f"no geometry row — fixture located here: midspan X = "
                     f"{f['x_mid']:.1f} mm, faces taken as the ROI edges")}


def support_chord(df, fix, L1):
    """The rigid line through the two support windows, as (x_l, v_l, x_r, v_r),
    or None if either window is too sparse to place."""
    v_l, x_l, n_l = L1.window_mean_V(df, fix["x_left"], L1.SUPPORT_HALF_WIDTH_MM)
    v_r, x_r, n_r = L1.window_mean_V(df, fix["x_right"], L1.SUPPORT_HALF_WIDTH_MM)
    if min(n_l, n_r) < 3:
        return None
    return x_l, v_l, x_r, v_r


# =============================================================================
# FLEXURAL OVERLAY — fixture, faces, live neutral axis and kinematics
# =============================================================================
class FlexuralOverlay:
    """Midspan and the supports (static), plus the faces, the neutral axis and a
    live curvature/deflection readout (per frame)."""

    def __init__(self, ax, fix, L1, deformed=False):
        self.ax, self.L1, self.fix = ax, L1, fix
        self.span = L1.FLEX_SPAN_MM
        # When the cloud is drawn deformed, everything measured in reference Y —
        # the faces and the neutral axis — has to travel with it or it floats off
        # the specimen it describes. Midspan and the supports are referenced in X,
        # which barely moves, so they stay put.
        self.deformed = deformed
        self.notes = [fix["note"]]
        self.x_mid, self.x_left = fix["x_mid"], fix["x_left"]
        self.x_right, self.sign = fix["x_right"], fix["sign"]
        self.y_lo, self.y_hi = fix["y_lo"], fix["y_hi"]
        self.y_bot, self.y_top = fix["y_bot"], fix["y_top"]
        self.y_c = 0.5 * (self.y_bot + self.y_top)
        self.draw_static()
        self.artists = []

    def describe(self):
        return self.notes

    def draw_static(self):
        """The X-referenced fixture, fixed in the world frame for the whole test."""
        ax = self.ax
        ax.axvline(self.x_mid, color="white", lw=0.8, ls="--", alpha=0.7, zorder=5)
        ax.annotate("midspan", xy=(self.x_mid, 0.98), xycoords=("data", "axes fraction"),
                    xytext=(3, -8), textcoords="offset points",
                    fontsize=6, color="white", zorder=8)
        for x_s in (self.x_left, self.x_right):
            ax.axvline(x_s, color="#9ecbff", lw=0.8, ls=":", alpha=0.8, zorder=5)
            ax.annotate("support", xy=(x_s, 0.98), xycoords=("data", "axes fraction"),
                        xytext=(3, -8), textcoords="offset points",
                        fontsize=6, color="#9ecbff", zorder=8)

    def update(self, df):
        for a in self.artists:
            a.remove()
        self.artists = []
        L1, ax = self.L1, self.ax

        # Midspan deflection referenced to the support chord — exactly Level 1's
        # reduction, so the number on screen is the number in the CSV.
        v_mid, _, n_mid = L1.window_mean_V(df, self.x_mid, L1.SUPPORT_HALF_WIDTH_MM)
        ch = support_chord(df, self.fix, L1)
        if ch is not None and n_mid >= 3:
            x_l, v_l, x_r, v_r = ch
            chord = L1.chord_at(self.x_mid, x_l, v_l, x_r, v_r)
            defl = (self.sign * (chord - v_mid)
                    * L1.support_offset_factor(x_l, x_r, self.x_mid, self.span))
        else:
            chord = defl = np.nan

        fit = L1.midspan_strain_fit(df, self.x_mid, self.y_lo, self.y_hi)
        slope = fit["slope"]
        if np.isfinite(slope) and slope != 0:
            kappa = abs(slope)
            na_y = -fit["icept"] / slope
            e_bot = slope * self.y_bot + fit["icept"]
            e_top = slope * self.y_top + fit["icept"]
        else:
            kappa = na_y = e_bot = e_top = np.nan

        # The deformed view is drawn chord-referenced (see xy() in
        # make_animation), so midspan's on-screen offset is its V measured
        # against that chord — the deflection itself. In the reference view the
        # shift is zero and this draws them as measured.
        dy = (v_mid - chord if (self.deformed and np.isfinite(v_mid)
                                and np.isfinite(chord)) else 0.0)
        half = L1.MIDSPAN_HALF_WIDTH_MM

        # The ROI is inset from both faces by the subset radius, so the surfaces
        # D790 asks about are drawn but were never measured — worth seeing.
        for y_f in (self.y_bot, self.y_top):
            fl, = ax.plot([self.x_mid - 3 * half, self.x_mid + 3 * half],
                          [y_f + dy, y_f + dy], color="#888888", lw=0.7,
                          ls="-.", alpha=0.8, zorder=5)
            self.artists.append(fl)

        if np.isfinite(na_y) and self.y_lo - 2 <= na_y <= self.y_hi + 2:
            ln, = ax.plot([self.x_mid - half, self.x_mid + half],
                          [na_y + dy, na_y + dy], color="lime", lw=1.4, zorder=7)
            self.artists.append(ln)

        # Returned rather than drawn: the ROI fills the axes edge to edge, so an
        # in-plot readout would sit on the specimen. The caller puts it in the title.
        return (f"$\\kappa$={kappa:.2e} 1/mm   NA={na_y:+.2f} mm   "
                f"$\\delta$={defl:.2f} mm   "
                f"$\\epsilon_{{bot}}$={e_bot:+.4f}   $\\epsilon_{{top}}$={e_top:+.4f}   "
                f"$R^2$={fit['r2']:.3f}")


# =============================================================================
# ANIMATION
# =============================================================================
def make_animation(cid, meta, fs, sel_idx, sel_paths, var, vmin, vmax,
                   cfg, geom, fig_dir, overlay_ok):
    kind = meta["kind"]
    L1 = flex_level1() if kind == "flexural" else None
    fix = resolve_fixture(geom, fs, L1) if (kind == "flexural" and overlay_ok) else None

    ref = fs["read"](sel_paths[0])
    if ref is None:
        raise RuntimeError("first sampled frame has no correlated points")

    # The deformed view is only drawn when it can be referenced to the support
    # chord. Without one, what is left in V is mostly rigid-body motion, and
    # drawing that reads as deflection when it isn't.
    deformed = cfg.get("deformed", False)
    if deformed and kind == "flexural":
        if fix is None:
            deformed = False
            print("  [!] no fixture — drawing reference coordinates, not the "
                  "deformed shape")
        elif support_chord(ref, fix, L1) is None:
            deformed = False
            print("  [!] support windows carry too few points to form a chord — "
                  "drawing reference\n      coordinates, not the deformed shape")

    def xy(df):
        """Where to draw each point: deformed position, or reference position.

        The deformed view is referenced to the support chord, not the DIC world
        frame. The raw V on these coupons carries ~15 mm of rigid-body motion
        against ~12 mm of actual midspan deflection, so drawn raw the beam sails
        up the screen while it bends. Subtracting the chord leaves bending alone
        and matches what the readout reports and what Level 1 writes.
        """
        x = df["X"].to_numpy(dtype=float)
        y = df["Y"].to_numpy(dtype=float)
        if not (deformed and "U" in df.columns and "V" in df.columns):
            return x, y
        xd = x + df["U"].to_numpy(dtype=float)
        yd = y + df["V"].to_numpy(dtype=float)
        if fix is not None:
            ch = support_chord(df, fix, L1)
            if ch is not None:
                x_l, v_l, x_r, v_r = ch
                yd = yd - L1.chord_at(x, x_l, v_l, x_r, v_r)
        return xd, yd

    x0, y0 = xy(ref)
    x_lo, x_hi = float(np.nanmin(x0)), float(np.nanmax(x0))
    y_lo, y_hi = float(np.nanmin(y0)), float(np.nanmax(y0))

    # The deformed view grows through the test — midspan travels about a full
    # specimen depth — so limits from frame 0 would clip the beam out of frame
    # exactly when it gets interesting. Take the union with the last frame.
    if deformed:
        last = fs["read"](sel_paths[-1])
        if last is not None:
            xn, yn = xy(last)
            x_lo, x_hi = min(x_lo, float(np.nanmin(xn))), max(x_hi, float(np.nanmax(xn)))
            y_lo, y_hi = min(y_lo, float(np.nanmin(yn))), max(y_hi, float(np.nanmax(yn)))

    if kind == "flexural" and FLEX_X_WINDOW_MM is not None and geom is not None:
        x_c = float(geom["x_mid_mm"])
        x_lo = max(x_lo, x_c - FLEX_X_WINDOW_MM)
        x_hi = min(x_hi, x_c + FLEX_X_WINDOW_MM)

    x_pad = (x_hi - x_lo) * 0.04
    y_pad = (y_hi - y_lo) * 0.08

    fig_w, fig_h = cfg["figsize"]
    if fig_h is None:
        fig_h = fig_w * ((y_hi - y_lo) / max(x_hi - x_lo, 1e-3)) + 1.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    valid0 = ref.dropna(subset=[var])
    vx0, vy0 = xy(valid0)
    scat = ax.scatter(vx0, vy0, c=valid0[var], cmap=cfg["cmap"],
                      norm=norm, s=cfg["marker_s"], linewidths=0, zorder=2)

    sm = plt.cm.ScalarMappable(cmap=cfg["cmap"], norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label(field_label(kind, var), fontsize=8, color="white")
    cb.ax.tick_params(labelsize=7, colors="white")
    plt.setp(cb.ax.yaxis.get_ticklines(), color="white")

    ax.set_xlim(x_lo - x_pad, x_hi + x_pad)
    ax.set_ylim(y_lo - y_pad, y_hi + y_pad)

    y_label = "Y (mm)"
    if kind == "flexural":
        if FLEX_ASPECT == "equal":
            ax.set_aspect("equal", adjustable="box")
            y_label = "Y — through depth (mm)"
        elif FLEX_ASPECT == "fill":
            ax.set_aspect("auto")     # y_label filled in after layout, below
        else:
            ax.set_aspect(float(FLEX_ASPECT), adjustable="box")
            y_label = f"Y — through depth (mm), drawn {FLEX_ASPECT:g}$\\times$ tall"
        ax.set_xlabel("X — along span (mm)"
                      + ("   [deformed shape]" if deformed else ""),
                      fontsize=8, color="white")
    else:
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("X (mm)", fontsize=8, color="white")
    ax.set_ylabel(y_label, fontsize=8, color="white")
    ax.tick_params(labelsize=7, colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#555555")

    overlay = None
    if overlay_ok:
        overlay = (TensileOverlay(ax, ref) if kind == "tensile"
                   else FlexuralOverlay(ax, fix, L1, deformed))
        for line in overlay.describe():
            print(f"  {line}")

    warn = "" if overlay_ok else "    [ROI orientation not verified]"
    title = ax.set_title("", fontsize=9, color="white", pad=6)

    def set_title(gi, readout):
        """Line 1 identifies the frame; line 2 is what the overlay measured on
        it. Keeping the live numbers in the title is what stops them landing on
        top of the specimen."""
        head = (f"{cid}   —   {field_label(kind, var)}   |   "
                f"frame {gi + 1} / {fs['n_total']}{warn}")
        title.set_text(head if not readout else f"{head}\n{readout}")

    set_title(sel_idx[0], overlay.update(ref) if overlay is not None else None)
    fig.tight_layout()

    # "fill" lets the axes take the shape of the figure, then reports the
    # exaggeration that produced — measured off the laid-out axes, so the number
    # on the label is the one the reader is looking at.
    if kind == "flexural" and FLEX_ASPECT == "fill":
        fig.canvas.draw()
        bb = ax.get_window_extent()
        (xa, xb), (ya, yb) = ax.get_xlim(), ax.get_ylim()
        exagg = (bb.height / (yb - ya)) / (bb.width / (xb - xa))
        ax.set_ylabel(f"Y — through depth (mm), drawn {exagg:.1f}$\\times$ tall",
                      fontsize=8, color="white")
        fig.tight_layout()

    def update(k):
        df = fs["read"](sel_paths[k])
        if df is None or var not in df.columns:
            return
        df = df.copy()
        df[var] = pd.to_numeric(df[var], errors="coerce")
        valid = df.dropna(subset=[var, "X", "Y"])
        if not len(valid):
            return
        vx, vy = xy(valid)
        scat.set_offsets(np.column_stack([vx, vy]))
        scat.set_array(valid[var].to_numpy())
        set_title(sel_idx[k], overlay.update(df) if overlay is not None else None)

    n_frames = len(sel_paths)
    mp4_path = fig_dir / f"{var}_field_animation.mp4"
    gif_path = fig_dir / f"{var}_field_animation.gif"

    # 1. imageio — imageio-ffmpeg bundles its own binary, no system ffmpeg needed
    try:
        import imageio
        with imageio.get_writer(str(mp4_path), fps=FPS, codec="libx264",
                                quality=8) as w:
            t0 = time.time()
            for k in range(n_frames):
                update(k)
                fig.canvas.draw()
                w.append_data(np.asarray(fig.canvas.buffer_rgba())[..., :3])
                if (k + 1) % 25 == 0 or k == n_frames - 1:
                    el = time.time() - t0
                    print(f"    [{k + 1}/{n_frames}] {el:.0f} s elapsed, "
                          f"ETA {el / (k + 1) * (n_frames - k - 1):.0f} s")
        plt.close(fig)
        return mp4_path
    except ImportError:
        print("  imageio not installed — run: pip install imageio imageio-ffmpeg")
    except Exception as ex:
        print(f"  imageio MP4 failed ({ex})")

    # 2. system ffmpeg
    ani = animation.FuncAnimation(fig, update, frames=n_frames,
                                  interval=1000 // FPS, blit=False)
    if "ffmpeg" in animation.writers.list():
        ani.save(str(mp4_path), dpi=120,
                 writer=animation.FFMpegWriter(fps=FPS, bitrate=1800,
                                               extra_args=["-vcodec", "libx264"]))
        plt.close(fig)
        return mp4_path

    # 3. GIF fallback
    print("  ffmpeg not found — saving GIF (slower/larger)")
    ani.save(str(gif_path), writer=animation.PillowWriter(fps=FPS), dpi=100)
    plt.close(fig)
    return gif_path


# =============================================================================
# ONE COUPON
# =============================================================================
def process(cid, var_override):
    t0 = time.time()
    print("=" * 70)
    print(f"DIC heatmap animation  —  {cid}")
    print("=" * 70)

    meta = parse_id(cid)
    kind = meta["kind"]
    cfg = TEST_TYPES[kind]

    requested = var_override or COLOR_VAR
    var = AUTO_COLOR_VAR[kind] if requested == "auto" else requested
    print(f"  Test type      : {kind}")
    print(f"  Colour field   : {var}{'  (auto)' if requested == 'auto' else ''}")

    fs = frame_set(cid, meta)
    print(f"  Frames         : {fs['source']}")

    ref0 = fs["read"](fs["paths"][0])
    if ref0 is None:
        raise RuntimeError("first frame has no correlated points")

    geom = flex_scalars(cid) if kind == "flexural" else None
    overlay_ok = True
    if kind == "flexural":
        print(f"  Truncated at fracture: frame {len(fs['paths'])} of {fs['n_total']}"
              + ("" if geom is not None else "  (located here, no Level-1 row)"))
        if geom is None:
            # Before locating the fixture here, confirm the frame is one where
            # that means anything: X along the span, Y through the depth.
            # Midspan and a neutral axis are meaningless otherwise.
            L1 = flex_level1()
            try:
                _, d_mm = L1.specimen_geometry(cid)
                complaint = L1.check_roi_orientation(ref0, d_mm)
            except Exception as ex:
                complaint = f"could not check ROI orientation ({ex})"
            if complaint:
                print(f"  [!] no geometry row and {complaint}")
                print(f"      Rendering the raw field with no bending overlays. "
                      f"Re-align this VIC-3D\n      project and re-run Level 1 to "
                      f"get them back.")
                overlay_ok = False
            else:
                print(f"  [!] no fixture row in {COUPON_SCALARS.name} — locating it here")

    if var not in ref0.columns:
        raise RuntimeError(f"'{var}' is not in these frames "
                           f"(available: {', '.join(map(str, ref0.columns))})")

    sel_idx, sel_paths = stride_frames(fs["paths"], FRAME_STRIDE)
    print(f"  Sampling       : {len(fs['paths'])} frames -> {len(sel_paths)} at "
          f"stride {FRAME_STRIDE}, {FPS} fps "
          f"({len(sel_paths) / FPS:.0f} s of video)")

    src = ("the last correlated frame" if cfg["range_from"] == "last"
           else f"{min(N_CBAR_SAMPLE, len(fs['paths']))} sampled frames")
    print(f"  Colour range   : from {src} …")
    vmin, vmax = compute_color_range(fs, var, cfg)
    lo_p, hi_p = cfg["pctile"]
    print(f"                   {vmin:+.5f} -> {vmax:+.5f}  ({lo_p}th-{hi_p}th pctile"
          f"{', forced symmetric' if cfg['symmetric'] else ''})")

    fig_dir = FIGS_ROOT / cid
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = make_animation(cid, meta, fs, sel_idx, sel_paths, var, vmin, vmax,
                         cfg, geom, fig_dir, overlay_ok)
    print(f"  Done in {time.time() - t0:.1f} s  ->  {out}")


# =============================================================================
# MAIN
# =============================================================================
def main():
    args = list(sys.argv[1:])
    var_override = None
    if "--var" in args:
        i = args.index("--var")
        if i + 1 >= len(args):
            print("[error] --var needs a field name")
            sys.exit(1)
        var_override = args[i + 1]
        del args[i:i + 2]

    if not args:
        print(__doc__)
        sys.exit(0)

    failures = []
    for cid in args:
        try:
            process(cid.strip(), var_override)
        except Exception as ex:
            print(f"  [error] {cid}: {ex}")
            failures.append(cid)

    if failures:
        print(f"\n{len(failures)} coupon(s) failed: {', '.join(failures)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
