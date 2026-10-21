#!/usr/bin/env python3
"""Published method on WINNER-TAKE-ALL 100 bp binned input, vs unbinned.

Tests the paper-repo session's addendum-4 hypothesis: the published scripts read
`...presmoothed.scaffolded.binned100.features.bed.gz`, and every reproduction so
far used the raw scaffolded BED with no `bin` step.

Targets, recovered from the published figure's vector geometry (addendum 3): the
distances and 99% CI widths in a2_step3 (published_delta_kb, published_ci_width_mb)
of the ROB configuration ($ROBS_CONFIG), which also gives the per-sample values.
"""
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

RES = 100
W = os.path.join(os.environ.get("KS_WORK", "$KS_WORK"), "a2_step3")
ROBS = [(s, d["scaffold"], d["left"], d["right"], d["sst1_mid"], d["win_start"], d["win_end"],
         d["published_delta_kb"], d["published_ci_width_mb"]) for s, d in samples().items()]
SRC = {
    "unbinned-presmoothed": os.path.join(W, "presm_{s}/{s}.haplotype1.HKS_human_CHM13_v2.chromosome.presmoothed.scaffolded.bed.gz"),
    "binned100 (prio)":     os.path.join(W, "bin100_{s}/binned100.prio.bed.gz"),
    "binned100 (plain)":    os.path.join(W, "bin100_{s}/binned100.plain.bed.gz"),
}


def analyse(path, scaffold, left, right, fs, fe):
    n = int((fe - fs) // RES)
    L, R = np.zeros(n), np.zeros(n)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != scaffold:
                continue
            lab = f[3].replace("_specific", "")
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


print(f"{'sample':9}{'input':24}{'centre':>13}{'Δ SST1':>11}{'target Δ':>11}"
      f"{'99% CI w':>10}{'target':>8}{'match?':>9}")
for s, scaf, l, r, mid, i0, i1, tgtd, tgtci in ROBS:
    fs, fe = i0 - 1_000_000, i1 + 1_000_000
    for tag, tmpl in SRC.items():
        try:
            res = analyse(tmpl.format(s=s), scaf, l, r, fs, fe)
        except (FileNotFoundError, EOFError, OSError):
            print(f"{s:9}{tag:24}  (not ready)")
            continue
        if res is None:
            print(f"{s:9}{tag:24}  (a partner has no content in the window)")
            continue
        c, a, b = res
        d = (c - mid) / 1000
        w = (b - a) / 1e6
        ok = "YES" if abs(d - tgtd) < 25 and abs(w - tgtci) < 0.25 else "no"
        print(f"{s:9}{tag:24}{c:>13,}{d:>10.1f}k{tgtd:>10.1f}k{w:>10.2f}{tgtci:>8.2f}{ok:>9}")
    print()
