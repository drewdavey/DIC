#!/usr/bin/env python3
"""
DIC Analysis Driver — FSR Tensile Coupons (P01)
================================================
Walks each selected coupon's directory, finds VIC-3D .out files,
exports each one's AOI data to a sibling .csv, and (later) makes
plots.

USAGE
-----
1. Toggle the SWITCHES section below to pick which coupons to run.
2. Make sure vicpyx is installed in the active Python environment:
       pip install vicpyx
3. Run:
       python DIC_Level1.py

OUTPUTS
-------
- One .csv per .out file, written *next to* the .out file inside
  each coupon's directory.
"""

from __future__ import annotations

import os
import glob
import time
from pathlib import Path
from typing import Iterable

from vicpyx import VicDataSet as VICDataSet

# =============================================================================
# PATHS
# =============================================================================
DATA_ROOTS = {
    "CL": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TCL"),
    "SW": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
    "UV": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
    "IS": Path(r"G:\DrewDavey\2026_FSR_TensileTest_TSW_TIS_TUV"),
}
FIGS_ROOT = Path(
    r"Z:\2023_07_SIO_Functional_Surfing_Reef\04_Drew"
    r"\01_MaterialTesting\02_Mechanical Testing\04_TestCoupons"
    r"\P01-LT150-LH4.5\figs"
)

# =============================================================================
# SWITCHES — toggle which coupons to process
# =============================================================================

# Print ID(s). Only P01 is processed for now.
PRINTS = ["P01"]

# Exposure types: "CL" (Control), "SW" (Seawater), "UV", "IS" (In-Situ)
EXPOSURES = {
    "CL": True,
    "SW": True,
    "UV": True,
    "IS": True,
}

# Print directions (degrees relative to load)
DIRECTIONS = {
    "00": True,
    "45": True,
    "90": True,
}

# Replicate numbers (1-3, zero-padded)
REPLICATES = ["01", "02", "03"]
# REPLICATES = ["01"]

# Pipeline steps — toggle each independently
DO_EXPORT_CSV = True    # export each .out to a CSV next to it
OVERWRITE_CSV = False   # if False, skip .out files whose .csv already exists

# Variables to export from each .out
# Standard full-field DIC variables; sigma is needed to filter invalid points.
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
# HELPERS
# =============================================================================

def coupon_id(print_id: str, exposure: str, direction: str, replicate: str) -> str:
    """Build coupon ID like P01-TCL00-01."""
    return f"{print_id}-T{exposure}{direction}-{replicate}"


def selected_coupons() -> list[str]:
    """Return list of coupon IDs based on SWITCHES."""
    out = []
    for p in PRINTS:
        for exp, exp_on in EXPOSURES.items():
            if not exp_on:
                continue
            for direction, dir_on in DIRECTIONS.items():
                if not dir_on:
                    continue
                for rep in REPLICATES:
                    out.append(coupon_id(p, exp, direction, rep))
    return out


def find_out_files(coupon_dir: Path) -> list[Path]:
    """Find all .out files inside a coupon's project directory.

    VIC-3D writes them in the project working directory; depending on
    your project structure they may be in the coupon dir directly or
    in a subdir. We search recursively.
    """
    return sorted(coupon_dir.rglob("*.out"))


def export_out_to_csv(out_path: Path, csv_path: Path,
                      var_names: Iterable[str]) -> bool:
    """Convert a single .out file to CSV using vicpyx

    Writes one row per valid AOI data point with columns named after
    var_names (sigma column dropped from the output, only used to
    filter invalid points).

    Returns True on success, False otherwise.
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
        values = ds.get_values(wanted)  # numpy structured array
    except Exception as ex:
        print(f"    [!] get_values failed on {out_path.name}: {ex}")
        return False

    import numpy as np
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


def process_coupon(cid: str) -> dict:
    """Run all enabled pipeline steps on one coupon.

    Returns a small dict of stats.
    """
    stats = {"coupon": cid, "n_out": 0, "n_csv_written": 0, "n_csv_skipped": 0}

    exposure = cid.split("-")[1][1:-2]  # "P01-TCL00-01" → "CL"
    coupon_dir = DATA_ROOTS[exposure] / cid
    if not coupon_dir.is_dir():
        print(f"  [skip] directory not found: {coupon_dir}")
        return stats

    out_files = find_out_files(coupon_dir)
    stats["n_out"] = len(out_files)
    if not out_files:
        print(f"  [skip] no .out files in {coupon_dir}")
        return stats

    print(f"  found {len(out_files)} .out files")

    # ---- CSV export ----
    if DO_EXPORT_CSV:
        for out_path in out_files:
            csv_path = out_path.with_suffix(".csv")
            if csv_path.exists() and not OVERWRITE_CSV:
                stats["n_csv_skipped"] += 1
                continue
            ok = export_out_to_csv(out_path, csv_path, EXPORT_VARS)
            if ok:
                stats["n_csv_written"] += 1

        print(f"  CSV: wrote {stats['n_csv_written']}, "
              f"skipped {stats['n_csv_skipped']} (already existed)")

    return stats

# =============================================================================
# MAIN
# =============================================================================

def main():
    t0 = time.time()
    print("=" * 70)
    print("FSR DIC Analysis Driver")
    print("=" * 70)
    for exp, root in DATA_ROOTS.items():
        print(f"Data root ({exp}): {root}")
    print(f"Figs root : {FIGS_ROOT}")
    print()

    coupons = selected_coupons()
    print(f"Processing {len(coupons)} coupon(s):")
    for c in coupons:
        print(f"  - {c}")
    print()

    if DO_PLOTS:
        FIGS_ROOT.mkdir(parents=True, exist_ok=True)

    summary = []
    for cid in coupons:
        print(f"[{cid}]")
        stats = process_coupon(cid)
        summary.append(stats)
        print()

    # ---- Final summary ----
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Coupon':<20} {'#out':>6} {'#csv':>6} {'#skip':>6}")
    for s in summary:
        print(f"{s['coupon']:<20} {s['n_out']:>6} "
              f"{s['n_csv_written']:>6} {s['n_csv_skipped']:>6}")

    dt = time.time() - t0
    print(f"\nTotal elapsed: {dt:.1f} s")


if __name__ == "__main__":
    main()
