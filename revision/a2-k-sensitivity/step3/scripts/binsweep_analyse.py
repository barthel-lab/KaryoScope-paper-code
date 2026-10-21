#!/usr/bin/env python3
"""A2 Step 3 k-sweep in the manuscript's own (binned) input representation.

Two methods per (sample, k), because they are not guaranteed to agree once the
input is winner-take-all binned:

  published : the legacy breakpoint_density.py computation -- focus window =
              _robertsonian ISCN CI +/- 1 Mb, 100 bp, plateau_center, CI at
              0.99*max. This is what main.tex:198 reports, so it is the primary.
  shipped   : the shipped refine_breakpoint + cumulative_breakpoint at
              pivot=chromosome over the whole scaffold, i.e. what the tool emits.

Emits sst1_breakpoints_by_k_binned.tsv in $KS_WORK/a2_step3. Per-sample values come
from the ROB configuration ($ROBS_CONFIG, section a2_step3).
"""
import gzip
import os
import sys
from pathlib import Path

import numpy as np

# KaryoScope-ISCN analysis modules: the copy vendored in this repository at commit 5507369
# (vendored/KaryoScope-ISCN/src). These scripts are copied to a working directory before they
# run, so point $KS_ISCN_SRC at that src/ directory.
sys.path.insert(0, os.environ.get("KS_ISCN_SRC", "$KS_ISCN_SRC"))
from karyoscope_iscn.core.bedio import iter_bed_filtered
from karyoscope_iscn.core.cumulative import cumulative_breakpoint
from karyoscope_iscn.core.detect import Breakpoint
from karyoscope_iscn.core.refine import refine_breakpoint
from karyoscope_iscn.reference import load_reference
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

RES = 100
KS = [11, 21, 31, 41, 51, 61]
W = os.path.join(os.environ.get("KS_WORK", "$KS_WORK"), "a2_step3")
ROBS = {s: dict(scaf=d["scaffold"], left=d["left"], right=d["right"], pos=d["detected_pos"],
                up_band=d["up_band"], dn_band=d["down_band"], length=d["scaffold_length"],
                sst1=d["sst1_mid"], win=(d["win_start"], d["win_end"]))
        for s, d in samples().items()}
BED = os.path.join(W, "binsweep/{s}/k{k}/binned100.plain.bed.gz")


def published(path, scaffold, left, right, fs, fe):
    n = int((fe - fs) // RES)
    L, R = np.zeros(n), np.zeros(n)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != scaffold:
                continue
            lab = f[3]
            if lab not in (left, right):
                continue
            s, e = int(f[1]), int(f[2])
            if e <= fs or s >= fe:
                continue
            i0, i1 = max(0, (s - fs) // RES), min(n, -(-(e - fs) // RES))
            tgt = L if lab == left else R
            for i in range(i0, i1):
                bs = fs + i * RES
                tgt[i] += max(0, min(e, bs + RES) - max(s, bs))
    if L.sum() == 0 or R.sum() == 0:
        return None
    pos = fs + np.arange(n) * RES
    lc, rc = np.cumsum(L), np.cumsum(R[::-1])[::-1]
    comb = lc / (lc[-1] or 1) + rc / (rc[0] or 1)
    mx = comb.max()
    at = np.where(np.isclose(comb, mx, rtol=1e-10))[0]
    ci = np.where(comb >= mx * 0.99)[0]
    return (int(pos[at[0]]) + int(pos[at[-1]])) // 2, int(pos[ci[0]]), int(pos[ci[-1]])


def main():
    ref = load_reference("chm13v2.0")
    out = Path(W) / "sst1_breakpoints_by_k_binned.tsv"
    rows = [("sample", "k", "scaffold", "pub_center", "pub_delta_bp", "pub_ci_start",
             "pub_ci_end", "pub_ci_width", "shipped_breakpoint", "shipped_delta_bp",
             "shipped_plateau_width", "sst1_mid", "status")]
    for s, d in ROBS.items():
        for k in KS:
            p = BED.format(s=s, k=k)
            try:
                pr = published(p, d["scaf"], d["left"], d["right"],
                               d["win"][0] - 1_000_000, d["win"][1] + 1_000_000)
            except (FileNotFoundError, EOFError, OSError):
                rows.append((s, k, d["scaf"], *["NA"] * 9, d["sst1"], "missing"))
                continue
            # iter_bed_filtered, NOT read_bed: read_bed materialises every record
            # in a ~20 M-line BED and OOM-killed this script at exit 137. This is
            # the same memory-bounded path analyze.py uses for exactly this reason.
            ivs = list(iter_bed_filtered(Path(p), chroms={d["scaf"]}))
            bp = Breakpoint(scaffold=d["scaf"], pos=d["pos"], up_chrom=d["left"],
                            down_chrom=d["right"], up_band=d["up_band"], down_band=d["dn_band"])
            refine_breakpoint(bp, ivs, ref, pivot="chromosome")
            cum = cumulative_breakpoint(bp, ivs, ref, (0, d["length"]), pivot="chromosome")
            if pr is None:
                rows.append((s, k, d["scaf"], *["NA"] * 4,
                             cum.inferred_breakpoint if cum else "NA",
                             (cum.inferred_breakpoint - d["sst1"]) if cum else "NA",
                             (cum.plateau_end - cum.plateau_start) if cum else "NA",
                             d["sst1"], "no_partner_content"))
                continue
            c, a, b = pr
            rows.append((s, k, d["scaf"], c, c - d["sst1"], a, b, b - a,
                         cum.inferred_breakpoint if cum else "NA",
                         (cum.inferred_breakpoint - d["sst1"]) if cum else "NA",
                         (cum.plateau_end - cum.plateau_start) if cum else "NA",
                         d["sst1"], "ok"))
    out.write_text("\n".join("\t".join(str(x) for x in r) for r in rows) + "\n")
    print(f"[ok] {out}\n")
    hdr = f"{'sample':9}{'k':>4}{'pub centre':>13}{'Δ SST1':>11}{'99% CI':>9}" \
          f"{'shipped bp':>13}{'Δ SST1':>11}{'plateau Mb':>12}"
    print(hdr)
    for r in rows[1:]:
        if r[-1] != "ok":
            print(f"{r[0]:9}{r[1]:>4}   {r[-1]}")
            continue
        print(f"{r[0]:9}{r[1]:>4}{r[3]:>13,}{r[4]/1000:>10.1f}k{r[7]/1e6:>9.2f}"
              f"{r[8]:>13,}{r[9]/1000:>10.1f}k{r[10]/1e6:>12.2f}")


main()
