#!/usr/bin/env python3
"""Faithful reimplementation of KaryoScope-ROBs/scripts/breakpoint_density.py.

Same method (symmetric cumulative k-mer density), same 100 bp resolution, same
reported statistic (plateau_center), same CI definition (combined >= 0.99*max),
and the same focus windows the published run used (the `_robertsonian` ISCN CI
from conf/robs.yaml, padded +/- 1 Mb -- these reproduce the axis limits printed
on the published supplementary figure).

Two windows are evaluated per sample:
  published  -- the window the manuscript used
  corrected  -- the whole inter-segment span, i.e. what a fixed pipeline would
                hand the refiner instead of a window centred on the wrong side
                of the unassignable acrocentric gap.

Difference from the original: it reads the *smoothed* scaffolded chromosome BED
(what the current pipeline uses) rather than the legacy presmoothed binned-100
BED, whose data/ directory is gitignored and gone. The CI-width comparison
against the published widths is the check on whether that matters.

Per-sample values (scaffold, partners, SST1 midpoint, ISCN window, published
width, corrected window) come from the ROB configuration ($ROBS_CONFIG,
section a2_step3).
"""
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

RES = 100
W = os.path.join(os.environ.get("KS_WORK", "$KS_WORK"), "a2_step3")
ROBS = [
    # sample, scaffold, left label, right label, sst1 mid, iscn_start, iscn_end
    (s, d["scaffold"], d["left"], d["right"], d["sst1_mid"], d["win_start"], d["win_end"])
    for s, d in samples().items()]
PUBLISHED_CI = {s: d["published_ci_width_mb"] for s, d in samples().items()}
BED = os.path.join(W, "gate_{s}/work/"
                   "{s}.haplotype1.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz")


def coverage(path, scaffold, left, right, fs, fe):
    n = int((fe - fs) // RES)
    L = np.zeros(n)
    R = np.zeros(n)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, lab = line.rstrip("\n").split("\t")[:4]
            if c != scaffold or (lab != left and lab != right):
                continue
            s = int(s); e = int(e)
            if e <= fs or s >= fe:
                continue
            i0 = max(0, (s - fs) // RES)
            i1 = min(n, -(-(e - fs) // RES))
            tgt = L if lab == left else R
            for i in range(i0, i1):
                bs = fs + i * RES
                tgt[i] += max(0, min(e, bs + RES) - max(s, bs))
    return L, R


def analyse(L, R, fs):
    n = len(L)
    pos = fs + np.arange(n) * RES
    lc = np.cumsum(L)
    rc = np.cumsum(R[::-1])[::-1]
    tl = lc[-1] or 1
    tr = rc[0] or 1
    comb = lc / tl + rc / tr
    mx = comb.max()
    idx = np.where(np.isclose(comb, mx, rtol=1e-10))[0]
    ps, pe = int(pos[idx[0]]), int(pos[idx[-1]])
    ci = np.where(comb >= mx * 0.99)[0]
    return (ps + pe) // 2, ps, pe, int(pos[ci[0]]), int(pos[ci[-1]]), mx


def main():
    # corrected windows: last left-dominant .. first sustained right-dominant,
    # from gap.sh, padded 1 Mb each side to match the published padding convention
    CORRECTED = {s: tuple(d["corrected_window"]) for s, d in samples().items()}
    print(f"{'sample':8} {'window':10} {'refined center':>15} {'d(SST1)':>10} "
          f"{'99% CI':>22} {'width':>8} {'pub width':>10}")
    for s, scaf, left, right, mid, i0, i1 in ROBS:
        path = BED.format(s=s)
        for tag, (fs, fe) in (("published", (i0 - 1000000, i1 + 1000000)),
                              ("corrected", CORRECTED[s])):
            L, R = coverage(path, scaf, left, right, fs, fe)
            c, ps, pe, c0, c1, mx = analyse(L, R, fs)
            print(f"{s:8} {tag:10} {c:15,} {(c-mid)/1000:+9.1f}k "
                  f"{c0/1e6:9.3f}-{c1/1e6:<8.3f}{(c1-c0)/1e6:8.2f} "
                  f"{PUBLISHED_CI[s] if tag=='published' else 0:10.2f}")


main()
