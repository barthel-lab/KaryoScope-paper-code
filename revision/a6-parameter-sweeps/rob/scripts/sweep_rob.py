#!/usr/bin/env python3
"""A6 Sweep 2 (ROB): window / coverage / dominance, and what they actually control.

Reads the EXISTING k=31 scaffolded chromosome BEDs from the A2 Step 3 sweep, so
there is no annotation step -- this is pure interpretation layer.

WHAT THE CODE ACTUALLY DOES, which the A6 brief's premise got wrong:

  The published Fig. 4 breakpoint is `cumulative_breakpoint`, and in the shipped
  pipeline (`cumulative_refine_report`) its window is the WHOLE SCAFFOLD, chosen
  deliberately -- a narrow window around the content-dominance interval truncates
  the upstream chromosome's mass and misplaces the junction. So the brief's three
  parameters (window / min_window_fraction / specificity) are inputs to
  `refine_breakpoint`, which produces the content-dominance INTERVAL, and they do
  not enter the published breakpoint at all.

  Three sweeps therefore, not one:
    A  refine_breakpoint       <- the brief's three parameters, on the quantity
                                  they control (the ISCN interval / centre)
    B  cumulative_breakpoint   <- resolution and plateau_threshold, the parameters
                                  that DO produce the published number
    C  coupling test           <- run cumulative over the LEGACY window
                                  (refine interval +/- 1 Mb, DEFAULT_WINDOW_PAD)
                                  while varying A's parameters, to measure how
                                  much the whole-scaffold choice insulates the
                                  published result

Held fixed throughout: k=31, the k=31 scaffold layout, the detected event, the
database version, pivot=chromosome, fine_bin, fallback_pad, max_overlap_windows.
All coordinates are SCAFFOLDED frame (never mixed with iscn_regions.tsv, which is
original-assembly frame).

The pinned per-sample values come from the ROB configuration ($ROBS_CONFIG):
section a6_rob names the section that holds them (a2_step3).
"""
import os
import sys
from pathlib import Path

import yaml

# KaryoScope-ISCN analysis modules: the copy vendored in this repository at commit 5507369
# (vendored/KaryoScope-ISCN/src). These scripts are copied to a working directory before they
# run, so point $KS_ISCN_SRC at that src/ directory.
sys.path.insert(0, os.environ.get("KS_ISCN_SRC", "$KS_ISCN_SRC"))

from karyoscope_iscn.core.bedio import read_bed
from karyoscope_iscn.core.detect import Breakpoint
from karyoscope_iscn.core.refine import refine_breakpoint
from karyoscope_iscn.core.cumulative import (
    cumulative_breakpoint,
    DEFAULT_RESOLUTION,
    DEFAULT_PLATEAU_THRESHOLD,
    DEFAULT_WINDOW_PAD,
)
from karyoscope_iscn.reference import load_reference

# Identical to a2_step3/sweep_refine.py -- same pinned detection, same SST1
# midpoints re-derived from our own acrocentric annotations ($ROBS_CONFIG).
def _fixed():
    c = yaml.safe_load(open(os.environ["ROBS_CONFIG"]))
    out = {}
    for s, d in c[c["a6_rob"]["samples_from"]]["samples"].items():
        g = c["breakpoint_groups"][d["breakpoint_group"]][0]
        out[s] = dict(scaffold=g["scaffold"], pos=d["detected_pos"], up=g["chrom1"],
                      down=g["chrom2"], up_band=d["up_band"], down_band=d["down_band"],
                      length=d["scaffold_length"], sst1=d["sst1_mid"])
    return out


FIXED = _fixed()

# Defaults, from the shipped module constants (asserted below, so a future
# upstream change to a default breaks this loudly instead of silently).
D_WINDOW, D_FRAC, D_SPEC = 100_000, 0.25, 0.50
D_RES, D_PLAT = DEFAULT_RESOLUTION, DEFAULT_PLATEAU_THRESHOLD
assert (D_RES, D_PLAT) == (100, 0.99), (D_RES, D_PLAT)

WINDOWS = [25_000, 50_000, 100_000, 200_000, 500_000]
FRACS = [0.10, 0.25, 0.40, 0.50]
SPECS = [0.50, 0.60, 0.75, 0.90]
RESOLUTIONS = [100, 500, 1_000, 5_000, 10_000]
PLATEAUS = [0.90, 0.95, 0.99, 0.995, 1.0]

ROOT = Path(os.environ.get("KS_WORK", "$KS_WORK")) / "a2_step3" / "sweep"
OUT = Path(os.environ.get("KS_WORK", "$KS_WORK")) / "a6_sweeps" / "rob"

COLS = ["sweep", "sample", "parameter", "value", "is_default", "scaffold",
        "refine_center", "refine_start", "refine_end", "refine_width", "refined",
        "cum_breakpoint", "plateau_start", "plateau_end", "plateau_width",
        "sst1_mid", "delta_sst1", "delta_baseline", "window_lo", "window_hi", "status"]


def row(**kw):
    return "\t".join(str(kw.get(c, "NA")) for c in COLS)


def main():
    ref = load_reference("chm13v2.0")
    out_rows = []

    for sample, f in FIXED.items():
        bed = ROOT / sample / "k31" / "chromosome.scaffolded.bed.gz"
        if not bed.exists():
            sys.exit(f"FATAL: missing {bed}")
        ivs = [i for i in read_bed(bed) if i.chrom == f["scaffold"]]
        if not ivs:
            sys.exit(f"FATAL: no intervals on {f['scaffold']} in {bed}")
        bp = Breakpoint(scaffold=f["scaffold"], pos=f["pos"], up_chrom=f["up"],
                        down_chrom=f["down"], up_band=f["up_band"], down_band=f["down_band"])
        print(f"[{sample}] {len(ivs)} intervals on {f['scaffold']}", flush=True)

        # ---- baseline: the shipped configuration, whole-scaffold window ----
        base_r = refine_breakpoint(bp, ivs, ref, window=D_WINDOW,
                                   min_window_fraction=D_FRAC, specificity=D_SPEC,
                                   pivot="chromosome")
        base_c = cumulative_breakpoint(bp, ivs, ref, (0, f["length"]),
                                       resolution=D_RES, plateau_threshold=D_PLAT,
                                       pivot="chromosome")
        assert base_c is not None, sample
        baseline = base_c.inferred_breakpoint
        print(f"[{sample}] baseline cumulative={baseline} refine_center={base_r.center}", flush=True)

        def emit(sweep, param, value, r, c, win):
            return row(sweep=sweep, sample=sample, parameter=param, value=value,
                       is_default=int(value == {"window": D_WINDOW, "min_window_fraction": D_FRAC,
                                                "specificity": D_SPEC, "resolution": D_RES,
                                                "plateau_threshold": D_PLAT}.get(param)),
                       scaffold=f["scaffold"],
                       refine_center=r.center if r else "NA",
                       refine_start=r.start if r else "NA",
                       refine_end=r.end if r else "NA",
                       refine_width=(r.end - r.start) if r else "NA",
                       refined=int(r.refined) if r else "NA",
                       cum_breakpoint=c.inferred_breakpoint if c else "NA",
                       plateau_start=c.plateau_start if c else "NA",
                       plateau_end=c.plateau_end if c else "NA",
                       plateau_width=(c.plateau_end - c.plateau_start) if c else "NA",
                       sst1_mid=f["sst1"],
                       delta_sst1=(c.inferred_breakpoint - f["sst1"]) if c else "NA",
                       delta_baseline=(c.inferred_breakpoint - baseline) if c else "NA",
                       window_lo=win[0], window_hi=win[1],
                       status="ok" if c else "no_chromosome_content")

        # ---- Sweep A: refine parameters, one at a time ----
        for param, values in (("window", WINDOWS), ("min_window_fraction", FRACS),
                              ("specificity", SPECS)):
            for v in values:
                kw = dict(window=D_WINDOW, min_window_fraction=D_FRAC, specificity=D_SPEC)
                kw[param] = v
                r = refine_breakpoint(bp, ivs, ref, pivot="chromosome", **kw)
                # cumulative held at the SHIPPED whole-scaffold window, so this row
                # shows the published number is untouched by the refine parameters.
                out_rows.append(emit("A_refine", param, v, r, base_c, (0, f["length"])))

                # ---- Sweep C: same refine parameters, but the LEGACY window ----
                lo = max(0, r.start - DEFAULT_WINDOW_PAD)
                hi = min(f["length"], r.end + DEFAULT_WINDOW_PAD)
                c = cumulative_breakpoint(bp, ivs, ref, (lo, hi), resolution=D_RES,
                                          plateau_threshold=D_PLAT, pivot="chromosome")
                out_rows.append(emit("C_legacy_window", param, v, r, c, (lo, hi)))

        # ---- Sweep B: cumulative parameters, one at a time ----
        for param, values in (("resolution", RESOLUTIONS), ("plateau_threshold", PLATEAUS)):
            for v in values:
                kw = dict(resolution=D_RES, plateau_threshold=D_PLAT)
                kw[param] = v
                c = cumulative_breakpoint(bp, ivs, ref, (0, f["length"]),
                                          pivot="chromosome", **kw)
                out_rows.append(emit("B_cumulative", param, v, base_r, c, (0, f["length"])))

    dest = OUT / "rob_param_sweep.tsv"
    dest.write_text("\t".join(COLS) + "\n" + "\n".join(out_rows) + "\n")
    print(f"\nwrote {dest} ({len(out_rows)} rows)")


main()
