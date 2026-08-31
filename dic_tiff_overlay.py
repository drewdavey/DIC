#!/usr/bin/env python3
"""
dic_tiff_overlay.py  —  FSR Flexural Coupons
============================================
A before/after figure that puts the bending model back on the photograph it came
from: the raw camera-0 TIFF with the ROI, the assumed specimen faces, the
assumed centroid and the measured neutral axis drawn on top.

WHY. Level 1 reports extreme-fibre strains at the two specimen faces, but those
faces are never measured — the ROI is inset from both surfaces by the subset
radius (on FCL0001 it covers 7.4 mm of a 12.7 mm depth, so ~2.6 mm at each face
is extrapolation), and their position rests on one assumption in
FlexuralDIC_Level1.build_l1: the centroid is ROI mid-depth, faces at +/- d/2.
Nothing in the point cloud can check that — the cloud has no edges in it. The
photograph does.

HOW THE MODEL GETS ONTO THE IMAGE. The per-frame CSVs carry both coordinate
systems for every subset: world millimetres (X, Y, U, V) and camera-0 pixels
(x, y, u, v). That is a dense exact correspondence, so the map from reference
millimetres to deformed pixels is fitted from each frame's own data:

    px, py = P(X, Y),   P = degree POLY_NX in X, degree POLY_NY in Y

Degree 6 in X carries lens distortion along the 250 mm span and the deflected
shape. Degree 1 in Y is deliberate: the face lines are extrapolated ~2.6 mm past
the last measured point and a quadratic Y term is free to do anything out there
(it also buys nothing — a quadratic-Y fit moves the faces 0.3 px). Fitting to
(x + u, y + v) is what puts reference coordinates on the *deformed* photograph.
Per-panel residuals are printed.

USAGE
    python dic_tiff_overlay.py <coupon_id> [<coupon_id> ...]
                               [--frames a,b] [--no-strain] [--zoom-mm 25]

    e.g.  python dic_tiff_overlay.py P01-FCL00-01
          python dic_tiff_overlay.py P01-FCL00-01 --frames 0,600
          python dic_tiff_overlay.py P01-FCL00-01 P01-FIS00-02 --no-strain

    --frames    two frame indices, before and after. Default: the first frame
                and the last correlated one (fracture).
    --no-strain omit the exx point cloud, leaving the photograph and the lines.
    --zoom-mm   half-width of the midspan zoom panel. Default 25 mm.

INPUT per coupon
  <coupon_dir>/<stem>-????????_0.csv   per-frame CSV (needs the pixel columns)
  <coupon_dir>/<stem>-????????_0.tif   the camera-0 TIFF beside it
  <DIC_DIR>/coupon_scalars.csv         fixture and assumed faces, from Level 1
  <DIC_DIR>/<coupon_id>.csv            Level-1 per-frame record, for the load

OUTPUT
  <FIGS_ROOT>/<coupon_id>/tiff_overlay_before_after.png
"""

from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dic_heatmap_animation as A

sys.stdout.reconfigure(encoding="utf-8")

# =============================================================================
# SWITCHES
# =============================================================================
ZOOM_HALF_MM  = 25.0     # midspan zoom panel half-width along the span
DEPTH_MARGIN  = 0.55     # crop margin past the faces, as a fraction of d
STRAIN_ALPHA  = 0.55     # opacity of the exx cloud over the photograph
STRAIN_PCTILE = (2, 98)  # same range rule as the flexural animation
NA_MIN_R2     = 0.80     # below this the exx(Y) fit is noise — see na_is_meaningful
POLY_NX       = 6        # mm -> px map: degree along the span
POLY_NY       = 1        # mm -> px map: degree through the depth (see docstring)
DPI           = 200

C_ROI  = "#4ec9ff"       # measured ROI edges
C_FACE = "#ffb454"       # assumed specimen faces (extrapolated)
C_CENT = "#ff77c8"       # assumed centroid == ROI mid-depth
C_NA   = "#66ff66"       # measured neutral axis
C_EDGE = "#ff5f5f"       # specimen surfaces, measured off the image
C_MID  = "white"
C_SUPP = "#9ecbff"


# =============================================================================
# mm -> px MAP
# =============================================================================
def poly_basis(X, Y, nx=POLY_NX, ny=POLY_NY):
    """Tensor-product monomials. X is scaled by 100 mm to keep the fit conditioned."""
    Xs = np.asarray(X, dtype=float) / 100.0
    Yv = np.asarray(Y, dtype=float)
    return np.column_stack([(Xs ** i) * (Yv ** j)
                            for j in range(ny + 1) for i in range(nx + 1)])


class MMtoPX:
    """Reference millimetres -> this frame's camera-0 pixels, fitted on its own
    points. Fitted to the deformed pixel position (x + u, y + v), so a
    reference-frame quantity is drawn where that material actually sits in this
    photograph."""

    def __init__(self, df):
        need = ["X", "Y", "x", "y", "u", "v"]
        miss = [c for c in need if c not in df.columns]
        if miss:
            raise RuntimeError(
                f"frame is missing pixel column(s) {miss}. Image coordinates come "
                f"from the VIC-3D per-frame CSV export; a coupon read through the "
                f"vicpyx .out fallback has world coordinates only, so run "
                f"FlexuralDIC_Level1 Step A on it first")
        d = df[need].apply(pd.to_numeric, errors="coerce").dropna()
        if len(d) < 200:
            raise RuntimeError(f"only {len(d)} points carry both coordinate systems")
        B = poly_basis(d["X"].to_numpy(), d["Y"].to_numpy())
        px = (d["x"] + d["u"]).to_numpy()
        py = (d["y"] + d["v"]).to_numpy()
        self.cx = np.linalg.lstsq(B, px, rcond=None)[0]
        self.cy = np.linalg.lstsq(B, py, rcond=None)[0]
        self.rms = (float(np.sqrt(np.mean((px - B @ self.cx) ** 2))),
                    float(np.sqrt(np.mean((py - B @ self.cy) ** 2))))
        self.n = len(d)

    def __call__(self, X, Y):
        X = np.atleast_1d(np.asarray(X, dtype=float))
        Y = np.atleast_1d(np.asarray(Y, dtype=float))
        X, Y = np.broadcast_arrays(X, Y)
        B = poly_basis(X.ravel(), Y.ravel())
        return (B @ self.cx).reshape(X.shape), (B @ self.cy).reshape(X.shape)

    def hline(self, y_mm, x0, x1, n=200):
        """A constant-Y line across the span, sampled so it follows the map's curve."""
        xs = np.linspace(x0, x1, n)
        return self(xs, np.full(n, float(y_mm)))

    def vline(self, x_mm, y0, y1, n=60):
        ys = np.linspace(y0, y1, n)
        return self(np.full(n, float(x_mm)), ys)


# =============================================================================
# GEOMETRY AND BENDING STATE
# =============================================================================
def face_geometry(cid, fix, L1):
    """Assumed centroid and faces, plus a note saying where they came from.

    resolve_fixture() hands back y_bot/y_top from coupon_scalars.csv when Level 1
    has a row. Without one it falls back to the ROI edges, which makes the
    extrapolation this figure is about invisible — so the faces are rebuilt here
    exactly as build_l1 would: ROI mid-depth, plus and minus half the
    specimen-sheet depth.
    """
    y_lo, y_hi = fix["y_lo"], fix["y_hi"]
    y_c = 0.5 * (y_lo + y_hi)
    y_bot, y_top = fix["y_bot"], fix["y_top"]
    d_mm = float(y_top - y_bot)
    src = f"assumed faces from DIC/{A.COUPON_SCALARS.name}"

    if abs(y_bot - y_lo) < 1e-6 and abs(y_top - y_hi) < 1e-6:
        try:
            _, d_mm = L1.specimen_geometry(cid)
            y_bot, y_top = y_c - d_mm / 2, y_c + d_mm / 2
            src = "no Level-1 row — faces rebuilt here from the specimen sheet"
        except Exception as ex:
            src = (f"no Level-1 row and no specimen depth ({ex}) — faces shown at "
                   f"the ROI edges, i.e. no extrapolation drawn")

    return {"y_lo": y_lo, "y_hi": y_hi, "y_c": y_c, "y_bot": y_bot,
            "y_top": y_top, "d_mm": d_mm, "note": src}


def bending_state(df, fix, geo, L1):
    """Level 1's midspan reduction on one frame: curvature, neutral axis, deflection."""
    fit = L1.midspan_strain_fit(df, fix["x_mid"], geo["y_lo"], geo["y_hi"])
    slope, icept = fit["slope"], fit["icept"]
    if np.isfinite(slope) and slope != 0:
        na_y = -icept / slope
        e_bot = slope * geo["y_bot"] + icept
        e_top = slope * geo["y_top"] + icept
        kappa = abs(slope)
    else:
        na_y = e_bot = e_top = kappa = np.nan

    v_mid, _, n_mid = L1.window_mean_V(df, fix["x_mid"], L1.SUPPORT_HALF_WIDTH_MM)
    ch = A.support_chord(df, fix, L1)
    if ch is not None and n_mid >= 3:
        x_l, v_l, x_r, v_r = ch
        chord = L1.chord_at(fix["x_mid"], x_l, v_l, x_r, v_r)
        defl = (fix["sign"] * (chord - v_mid)
                * L1.support_offset_factor(x_l, x_r, fix["x_mid"], L1.FLEX_SPAN_MM))
    else:
        defl = np.nan

    return {"kappa": kappa, "na_y": na_y, "e_bot": e_bot, "e_top": e_top,
            "defl": defl, "r2": fit["r2"], "n": fit["n"]}


def level1_row(cid, step):
    """That frame's row of the Level-1 per-frame CSV, for the load readout."""
    fp = A.DIC_DIR / f"{cid}.csv"
    if not fp.exists():
        return None
    try:
        g = pd.read_csv(fp)
    except Exception:
        return None
    r = g.loc[g["step"] == step]
    return None if r.empty else r.iloc[0].to_dict()


def na_is_meaningful(st, geo):
    """Is this frame's neutral axis worth drawing?

    On the first frames there is no moment: exx is noise about zero, the fit
    has no slope worth speaking of, and its zero crossing lands tens of
    millimetres off the specimen. Drawing that is worse than drawing nothing —
    it reads as a measurement. R^2 is the test Level 1 already uses for whether
    Euler-Bernoulli holds on a frame.
    """
    if not np.isfinite(st["na_y"]) or not np.isfinite(st["r2"]):
        return False
    if st["r2"] < NA_MIN_R2:
        return False
    pad = 0.75 * geo["d_mm"]
    return geo["y_bot"] - pad <= st["na_y"] <= geo["y_top"] + pad


# =============================================================================
# IMAGE
# =============================================================================
def tiff_for(csv_path):
    """Camera-0 TIFF beside the per-frame CSV: <stem>_0.csv -> <stem>_0.tif."""
    tif = csv_path.with_suffix(".tif")
    if not tif.exists():
        raise RuntimeError(f"no TIFF beside {csv_path.name} (looked for {tif.name})")
    return tif


def load_tiff(fp):
    try:
        from PIL import Image
    except ImportError as ex:
        raise RuntimeError(f"reading the TIFF needs Pillow ({ex}) — "
                           f"pip install Pillow") from ex
    with Image.open(fp) as im:
        return np.asarray(im.convert("L"))


def measure_edges(img, m, fix, geo, half_mm=20.0, n_cols=201, n_rows=1201,
                  smooth_mm=1.0):
    """The specimen's two surfaces as the camera sees them, in reference mm —
    the measurement the point cloud cannot make.

    The specimen is lit and speckled and everything around it is unlit, so
    averaging along the span at midspan gives a depth profile with a dark floor
    and a bright plateau; the half-height crossings are the two surfaces. The
    plateau must be the CONTIGUOUS run through the ROI, not every row above the
    threshold: on FIS00-02 the loading nose sits a millimetre above the top
    surface and is brighter than the specimen, and taken globally it merges with
    it and puts the measured top face 4 mm out.

    Returns (y_bot, y_top) in reference mm, or None if there is no plateau to
    find — a dark or badly cropped frame.
    """
    ys = np.linspace(geo["y_c"] - 1.4 * geo["d_mm"],
                     geo["y_c"] + 1.4 * geo["d_mm"], n_rows)
    xs = np.linspace(fix["x_mid"] - half_mm, fix["x_mid"] + half_mm, n_cols)
    XX, YY = np.meshgrid(xs, ys)
    px, py = m(XX, YY)
    rr = np.clip(np.round(py).astype(int), 0, img.shape[0] - 1)
    cc = np.clip(np.round(px).astype(int), 0, img.shape[1] - 1)
    prof = img[rr, cc].mean(axis=1)

    # Averaging along the span alone does not flatten the plateau: the speckle
    # dots are ~0.5 mm and their rows line up, so the raw profile still swings
    # from 30 to 110 counts inside the specimen and a half-height threshold cuts
    # through the middle of it. Smoothing through the depth over about one dot
    # leaves the two surfaces, the only edges that survive it.
    w = max(3, int(round(smooth_mm / (ys[1] - ys[0]))) | 1)
    prof = np.convolve(prof, np.ones(w) / w, mode="same")
    prof[:w] = prof[w]
    prof[-w:] = prof[-w - 1]

    floor, plateau = np.percentile(prof, 5), np.percentile(prof, 95)
    if plateau - floor < 8:                     # nothing lit in this crop
        return None
    thr = 0.5 * (floor + plateau)

    i_c = int(np.argmin(np.abs(ys - geo["y_c"])))
    if prof[i_c] < thr:
        return None
    lo = i_c
    while lo > 0 and prof[lo - 1] >= thr:
        lo -= 1
    hi = i_c
    while hi < len(prof) - 1 and prof[hi + 1] >= thr:
        hi += 1
    if lo == 0 or hi == len(prof) - 1:          # plateau runs off the window
        return None

    def cross(i_in, i_out):
        """Sub-row half-height crossing between an inside and an outside sample."""
        a, b = prof[i_out], prof[i_in]
        f = 0.0 if b == a else (thr - a) / (b - a)
        return ys[i_out] + f * (ys[i_in] - ys[i_out])

    return float(cross(lo, lo - 1)), float(cross(hi, hi + 1))


# =============================================================================
# ONE PANEL
# =============================================================================
def draw_panel(ax, img, m, fix, geo, st, half, edges, *, x0_mm, x1_mm, df=None,
               norm=None, cmap=None, show_labels=True):
    """One photograph with the model on it, cropped to x0_mm..x1_mm of the span."""
    # Set the display range from the crop, not the whole 4112 x 3008 frame: most
    # of that frame is unlit background, and scaling to it compresses the speckle
    # — and the specimen edge this figure exists to show — into the top of the range.
    pad = DEPTH_MARGIN * geo["d_mm"]
    cx, cy = m(np.array([x0_mm, x1_mm, x0_mm, x1_mm]),
               np.array([geo["y_bot"] - pad] * 2 + [geo["y_top"] + pad] * 2))
    px0, px1 = float(cx.min()), float(cx.max())
    py0, py1 = float(cy.min()), float(cy.max())
    sub = img[max(int(py0), 0):int(py1) + 1, max(int(px0), 0):int(px1) + 1]
    lo, hi = (np.percentile(sub, (1, 99)) if sub.size else (img.min(), img.max()))
    ax.imshow(img, cmap="gray", interpolation="nearest", zorder=0,
              vmin=float(lo), vmax=float(hi))

    if df is not None and norm is not None:
        d = df[["X", "Y", "exx"]].apply(pd.to_numeric, errors="coerce").dropna()
        if len(d):
            sx, sy = m(d["X"].to_numpy(), d["Y"].to_numpy())
            ax.scatter(sx, sy, c=d["exx"].to_numpy(), cmap=cmap, norm=norm,
                       s=1.6, linewidths=0, alpha=STRAIN_ALPHA, zorder=2)

    x_mid = fix["x_mid"]

    # Measured ROI edges — the last real data through the depth.
    for y_e in (geo["y_lo"], geo["y_hi"]):
        ax.plot(*m.hline(y_e, x0_mm, x1_mm), color=C_ROI, lw=0.9, ls=":",
                alpha=0.95, zorder=4)

    # Assumed faces — extrapolated, and the point of the figure.
    for y_f in (geo["y_bot"], geo["y_top"]):
        ax.plot(*m.hline(y_f, x0_mm, x1_mm), color=C_FACE, lw=1.1, ls="-.",
                alpha=0.95, zorder=5)

    # Assumed centroid: ROI mid-depth, which is what puts the faces where they are.
    ax.plot(*m.hline(geo["y_c"], x0_mm, x1_mm), color=C_CENT, lw=0.9,
            ls=(0, (1, 2)), alpha=0.95, zorder=5)

    # The surfaces the camera actually sees — the check on that assumption.
    if edges is not None:
        for y_e in edges:
            ax.plot(*m.hline(y_e, x0_mm, x1_mm), color=C_EDGE, lw=1.0,
                    alpha=0.95, zorder=6)

    # Measured neutral axis, over the window the curvature fit actually used.
    na_ok = na_is_meaningful(st, geo)
    if na_ok:
        nx0, nx1 = max(x0_mm, x_mid - half), min(x1_mm, x_mid + half)
        if nx1 > nx0:
            ax.plot(*m.hline(st["na_y"], nx0, nx1), color=C_NA, lw=1.8,
                    alpha=0.95, zorder=7)

    # Fixture, in the world frame.
    y0v, y1v = geo["y_bot"] - 0.6 * geo["d_mm"], geo["y_top"] + 0.6 * geo["d_mm"]
    if x0_mm <= x_mid <= x1_mm:
        ax.plot(*m.vline(x_mid, y0v, y1v), color=C_MID, lw=0.8, ls="--",
                alpha=0.7, zorder=4)
    for x_s in (fix["x_left"], fix["x_right"]):
        if x0_mm <= x_s <= x1_mm:
            ax.plot(*m.vline(x_s, y0v, y1v), color=C_SUPP, lw=0.8, ls=":",
                    alpha=0.85, zorder=4)

    # Crop to the span window asked for, with enough depth to show the assumed
    # faces outside the ROI and the real specimen edge between them. Image rows
    # increase downward, so the Y limits are reversed.
    ax.set_xlim(px0, px1)
    ax.set_ylim(py1, py0)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor("#555555")

    if not show_labels:
        return

    def tag(x_mm, y_mm, txt, col, dx=3, dy=2, ha="left"):
        lx, ly = m(x_mm, y_mm)
        ax.annotate(txt, xy=(float(np.ravel(lx)[0]), float(np.ravel(ly)[0])),
                    xytext=(dx, dy), textcoords="offset points", ha=ha,
                    fontsize=6, color=col, zorder=9,
                    bbox=dict(boxstyle="square,pad=0.12", fc=A.BG, ec="none",
                              alpha=0.6))

    x_l = x0_mm + 0.015 * (x1_mm - x0_mm)
    x_r = x1_mm - 0.015 * (x1_mm - x0_mm)
    tag(x_l, geo["y_top"], "assumed face", C_FACE)
    tag(x_l, geo["y_bot"], "assumed face", C_FACE, dy=-8)
    tag(x_l, geo["y_c"], "assumed centroid", C_CENT)
    tag(x_r, geo["y_hi"], "ROI edge", C_ROI, dx=-3, ha="right")
    if edges is not None:
        tag(x_r, edges[1], "measured surface", C_EDGE, dx=-3, dy=-9, ha="right")
    if na_ok:
        tag(min(x1_mm, x_mid + half), st["na_y"], "neutral axis", C_NA, dy=-8)


def legend_handles(extrap):
    return [
        Line2D([], [], color=C_ROI, ls=":", lw=1.2,
               label="measured ROI edge — last real data"),
        Line2D([], [], color=C_FACE, ls="-.", lw=1.4,
               label=f"assumed face (ROI mid-depth $\\pm$ d/2, "
                     f"{extrap:.2f} mm extrapolated)"),
        Line2D([], [], color=C_CENT, ls=(0, (1, 2)), lw=1.2,
               label="assumed centroid = ROI mid-depth"),
        Line2D([], [], color=C_EDGE, lw=1.4,
               label="specimen surface, measured off the image"),
        Line2D([], [], color=C_NA, lw=2.0, label="measured neutral axis"),
        Line2D([], [], color=C_MID, ls="--", lw=1.0, label="midspan"),
        Line2D([], [], color=C_SUPP, ls=":", lw=1.0, label="support"),
    ]


# =============================================================================
# ONE COUPON
# =============================================================================
def process(cid, frames, show_strain, zoom_mm):
    print("=" * 70)
    print(f"TIFF overlay  —  {cid}")
    print("=" * 70)

    meta = A.parse_id(cid)
    if meta["kind"] != "flexural":
        raise RuntimeError("this figure is the flexural neutral-axis / face check; "
                           "a tensile coupon has no through-depth model to draw")

    L1 = A.flex_level1()
    half = L1.MIDSPAN_HALF_WIDTH_MM

    fs = A.frame_set(cid, meta)
    print(f"  Frames         : {fs['source']}")
    if fs["read"] is not A.read_csv_frame:
        raise RuntimeError("this coupon is being read from .out files, which carry "
                           "no image coordinates. Run FlexuralDIC_Level1 Step A on "
                           "it to export the per-frame CSVs.")

    fix = A.resolve_fixture(A.flex_scalars(cid), fs, L1)
    print(f"  {fix['note']}")

    if frames is None:
        idx = (0, len(fs["paths"]) - 1)
    else:
        idx = tuple(min(max(i, 0), fs["n_total"] - 1) for i in frames)
    print(f"  Showing frames : {idx[0]} (before) and {idx[1]} (after) "
          f"of {fs['n_total']}")

    dfs = []
    for i in idx:
        df = fs["read"](fs["paths"][i])
        if df is None or not len(df):
            raise RuntimeError(f"frame {i} has no correlated points")
        dfs.append(df)

    geo = face_geometry(cid, fix, L1)
    extrap = 0.5 * (geo["d_mm"] - (geo["y_hi"] - geo["y_lo"]))
    print(f"  {geo['note']}")
    print(f"  ROI Y = {geo['y_lo']:.2f} .. {geo['y_hi']:.2f} mm "
          f"({geo['y_hi'] - geo['y_lo']:.2f} mm measured of d = {geo['d_mm']:.2f}); "
          f"centroid assumed at Y = {geo['y_c']:+.2f} mm, faces at "
          f"{geo['y_bot']:+.2f} / {geo['y_top']:+.2f} mm")
    print(f"  Extreme-fibre strains are extrapolated {extrap:.2f} mm past the last "
          f"measured point at each face")

    # Colour range from the "after" frame, the largest field shown — the same
    # rule the flexural animation uses.
    norm = cmap = None
    if show_strain:
        e = pd.to_numeric(dfs[1]["exx"], errors="coerce").dropna().to_numpy()
        lo, hi = np.percentile(e, STRAIN_PCTILE)
        v = max(abs(lo), abs(hi))
        norm = plt.Normalize(-v, v)
        cmap = A.TEST_TYPES["flexural"]["cmap"]
        try:
            plt.get_cmap(cmap)
        except Exception:
            cmap = "coolwarm"
        print(f"  Strain overlay : exx on {cmap}, +/-{v:.4f} "
              f"({STRAIN_PCTILE[0]}th-{STRAIN_PCTILE[1]}th pctile of frame {idx[1]})")

    x_mid = fix["x_mid"]
    x_full = (float(dfs[0]["X"].min()), float(dfs[0]["X"].max()))
    x_zoom = (max(x_full[0], x_mid - zoom_mm), min(x_full[1], x_mid + zoom_mm))

    # Layout: the full-ROI strip is 250 mm by ~26 mm of crop, so it only reads at
    # full figure width — one row each for before and after. The midspan zooms are
    # where the face lines are actually checked against the specimen edge, so they
    # get the tall row, side by side and without the strain cloud over that edge.
    fig = plt.figure(figsize=(13.5, 8.2))
    fig.patch.set_facecolor(A.BG)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 2.35],
                          hspace=0.42, wspace=0.06,
                          left=0.02, right=0.98, top=0.90, bottom=0.07)
    ax_full = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, :])]
    ax_zoom = [fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])]

    states, edge_list = [], []
    for r, (i, df) in enumerate(zip(idx, dfs)):
        img = load_tiff(tiff_for(fs["paths"][i]))
        m = MMtoPX(df)
        st = bending_state(df, fix, geo, L1)
        edges = measure_edges(img, m, fix, geo)
        states.append(st)
        edge_list.append(edges)

        drawn = "" if na_is_meaningful(st, geo) else "  [not drawn: fit is noise]"
        print(f"  frame {i:>4}: mm->px fit on {m.n} pts, rms "
              f"{m.rms[0]:.2f} / {m.rms[1]:.2f} px   NA = {st['na_y']:+.2f} mm "
              f"({st['na_y'] - geo['y_c']:+.2f} mm off the assumed centroid), "
              f"R2 = {st['r2']:.3f}{drawn}")
        if edges is None:
            print(f"            surfaces not measurable off this image — no "
                  f"half-height plateau through the ROI")
        else:
            e_c = 0.5 * (edges[0] + edges[1])
            print(f"            surfaces on image Y = {edges[0]:+.2f} .. "
                  f"{edges[1]:+.2f} mm ({edges[1] - edges[0]:.2f} mm vs sheet d = "
                  f"{geo['d_mm']:.2f}); assumed centroid is {geo['y_c'] - e_c:+.2f} "
                  f"mm off their midpoint")

        row = level1_row(cid, i)
        load = ""
        if row is not None and np.isfinite(row.get("force_N", np.nan)):
            load = f"   P = {row['force_N']:.0f} N"
        when = "before" if r == 0 else "after"
        stamp = f"{when}   frame {i} / {fs['n_total'] - 1}{load}"

        for ax, xw, lbl, strain in ((ax_full[r], x_full, "full ROI", show_strain),
                                    (ax_zoom[r], x_zoom,
                                     f"midspan $\\pm${zoom_mm:g} mm", False)):
            ax.set_facecolor(A.BG)
            draw_panel(ax, img, m, fix, geo, st, half, edges,
                       x0_mm=xw[0], x1_mm=xw[1],
                       df=df if strain else None, norm=norm, cmap=cmap,
                       show_labels=(ax is ax_zoom[r]))
            head = f"{stamp}   —   {lbl}"
            if ax is ax_zoom[r]:
                head += (f"\n$\\kappa$={st['kappa']:.2e} 1/mm   "
                         f"NA={st['na_y']:+.2f} mm   $\\delta$={st['defl']:.2f} mm"
                         f"\n$\\epsilon_{{bot}}$={st['e_bot']:+.4f}   "
                         f"$\\epsilon_{{top}}$={st['e_top']:+.4f}   "
                         f"$R^2$={st['r2']:.3f}")
            ax.set_title(head, fontsize=7, color="white", pad=4)

    leg = fig.legend(handles=legend_handles(extrap), loc="lower center", ncol=3,
                     fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, -0.004))
    for t in leg.get_texts():
        t.set_color("white")

    st_after = states[1]
    off = (f"NA {st_after['na_y'] - geo['y_c']:+.2f} mm off it under load"
           if na_is_meaningful(st_after, geo)
           else "no usable NA on the loaded frame")
    e0 = edge_list[0]
    cen = ("surfaces not measurable off the reference image" if e0 is None else
           f"centroid assumption {geo['y_c'] - 0.5 * (e0[0] + e0[1]):+.2f} mm off "
           f"the measured surfaces")
    fig.suptitle(
        f"{cid}   —   assumed faces and measured neutral axis on the raw camera-0 "
        f"image\nROI covers {geo['y_hi'] - geo['y_lo']:.2f} mm of d = "
        f"{geo['d_mm']:.2f} mm, so each face is {extrap:.2f} mm extrapolated   |   "
        f"{cen}   |   {off}",
        fontsize=9, color="white", y=0.985)

    fig_dir = A.FIGS_ROOT / cid
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = fig_dir / "tiff_overlay_before_after.png"
    fig.savefig(out, dpi=DPI, facecolor=A.BG)
    plt.close(fig)
    print(f"  ->  {out}")
    return out


# =============================================================================
# MAIN
# =============================================================================
def main():
    args = list(sys.argv[1:])
    frames = None
    show_strain = True
    zoom_mm = ZOOM_HALF_MM

    if "--frames" in args:
        i = args.index("--frames")
        frames = tuple(int(v) for v in args[i + 1].split(","))
        if len(frames) != 2:
            print("[error] --frames needs two indices, e.g. --frames 0,600")
            sys.exit(1)
        del args[i:i + 2]
    if "--no-strain" in args:
        show_strain = False
        args.remove("--no-strain")
    if "--zoom-mm" in args:
        i = args.index("--zoom-mm")
        zoom_mm = float(args[i + 1])
        del args[i:i + 2]

    if not args:
        print(__doc__)
        sys.exit(0)

    failures = []
    for cid in args:
        try:
            process(cid.strip(), frames, show_strain, zoom_mm)
        except Exception as ex:
            print(f"  [error] {cid}: {ex}")
            failures.append(cid)

    if failures:
        print(f"\n{len(failures)} coupon(s) failed: {', '.join(failures)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
