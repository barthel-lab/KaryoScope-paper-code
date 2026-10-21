#!/usr/bin/env python3
"""One row of the A2 Step 3 k-sweep: re-run the SHIPPED refinement at one k.

Usage: sweep_refine.py SAMPLE K SCAFFOLDED_BED OUT_TSV

Design -- what is held fixed and what varies:
  FIXED  the scaffold layout (k=31 scaffold_map, so every k reports in ONE
         coordinate frame -- a per-k layout could flip a contig and produce a
         spurious ~contig-length jump), and the DETECTED event (the k=31
         Breakpoint: scaffold, pos, up/down chromosome). Detection needs the
         cytoband feature set, which the k-sweep databases do not ship.
  VARIES only the chromosome-identity evidence fed to refinement -- exactly the
         quantity R3-major-3 asks about.

Runs the shipped `refine_breakpoint` + `cumulative_breakpoint` (pivot=chromosome),
NOT a reimplementation, so the sweep measures the tool. `cumulative_refine_report`
spans the whole scaffold, so the window here is (0, scaffold_length) to match.

All coordinates are SCAFFOLDED frame -- deliberately not back-mapped, so the
rows are directly comparable to each other and to conf/robs.yaml.

The pinned per-sample values come from the ROB configuration ($ROBS_CONFIG,
section a2_step3; see robs/conf/config.template.yaml).
"""
import os
import sys
# KaryoScope-ISCN analysis modules: the copy vendored in this repository at commit 5507369
# (vendored/KaryoScope-ISCN/src). These scripts are copied to a working directory before they
# run, so point $KS_ISCN_SRC at that src/ directory.
sys.path.insert(0, os.environ.get("KS_ISCN_SRC", "$KS_ISCN_SRC"))
from pathlib import Path

from karyoscope_iscn.core.bedio import read_bed
from karyoscope_iscn.core.detect import Breakpoint
from karyoscope_iscn.core.refine import refine_breakpoint
from karyoscope_iscn.core.cumulative import cumulative_breakpoint
from karyoscope_iscn.reference import load_reference
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

# k=31 detection, pinned (see module docstring). sst1 = midpoint re-derived from
# our own acrocentric annotations; block sizes reproduce the published ones
# exactly and the midpoints match robs.yaml to <= 27 bp.
FIXED = {s: dict(scaffold=d["scaffold"], pos=d["detected_pos"], up=d["left"], down=d["right"],
                 up_band=d["up_band"], down_band=d["down_band"], length=d["scaffold_length"],
                 sst1=d["sst1_mid"])
         for s, d in samples().items()}


def main():
    sample, k, bed, out = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    f = FIXED[sample]
    ref = load_reference("chm13v2.0")
    ivs = [i for i in read_bed(bed) if i.chrom == f["scaffold"]]
    bp = Breakpoint(scaffold=f["scaffold"], pos=f["pos"], up_chrom=f["up"],
                    down_chrom=f["down"], up_band=f["up_band"], down_band=f["down_band"])

    r = refine_breakpoint(bp, ivs, ref, pivot="chromosome")
    c = cumulative_breakpoint(bp, ivs, ref, (0, f["length"]), pivot="chromosome")

    if c is None:
        row = [sample, k, f["scaffold"], "NA", "NA", "NA", "NA", f["sst1"], "NA",
               r.center, int(r.refined), len(ivs), "no_chromosome_content"]
    else:
        row = [sample, k, f["scaffold"], c.inferred_breakpoint, c.plateau_start,
               c.plateau_end, c.plateau_end - c.plateau_start, f["sst1"],
               c.inferred_breakpoint - f["sst1"], r.center, int(r.refined),
               len(ivs), "ok"]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\t".join(str(x) for x in row) + "\n")
    print("\t".join(str(x) for x in row))


main()
