#!/usr/bin/env python3
"""v1 vs v2 chromosome annotation on the ROB fusion scaffolds, and the published
method's delta-to-SST1 computed on v1.

Answers: is the `chromosome` feature set actually the same between the two
databases (the trees are structurally identical), and does the published
"200-525 kb" range come out of v1?

Per-sample values come from the ROB configuration ($ROBS_CONFIG, section a2_step3).
"""
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

RES = 100
W = os.path.join(os.environ.get("KS_WORK", "$KS_WORK"), "a2_step3")
ROBS = [(s, d["scaffold"], d["left"], d["right"], d["sst1_mid"], d["win_start"], d["win_end"])
        for s, d in samples().items()]
V2 = os.path.join(W, "gate_{s}/work/"
                  "{s}.haplotype1.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz")
V1 = os.path.join(W, "v1_{s}/chromosome.scaffolded.bed.gz")


def norm(lab):
    """Map v1's vocabulary onto v2's so the two are comparable."""
    return lab.replace("_specific", "").replace("_multigroup1", "")


def load(path, scaffold):
    bp, iv = {}, 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != scaffold:
                continue
            bp[norm(f[3])] = bp.get(norm(f[3]), 0) + int(f[2]) - int(f[1])
            iv += 1
    return bp, iv


def published_method(path, scaffold, left, right, fs, fe):
    n = int((fe - fs) // RES)
    L, R = np.zeros(n), np.zeros(n)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != scaffold:
                continue
            lab = norm(f[3])
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
    pos = fs + np.arange(n) * RES
    lc, rc = np.cumsum(L), np.cumsum(R[::-1])[::-1]
    comb = lc / (lc[-1] or 1) + rc / (rc[0] or 1)
    mx = comb.max()
    at = np.where(np.isclose(comb, mx, rtol=1e-10))[0]
    ci = np.where(comb >= mx * 0.99)[0]
    return (int(pos[at[0]]) + int(pos[at[-1]])) // 2, int(pos[ci[0]]), int(pos[ci[-1]])


print("=" * 96)
print("PART 1 -- is the chromosome feature set the same? (bp per label, fusion scaffold)")
print("=" * 96)
for s, scaf, l, r, mid, i0, i1 in ROBS:
    b1, n1 = load(V1.format(s=s), scaf)
    b2, n2 = load(V2.format(s=s), scaf)
    print(f"\n{s}  ({scaf})   v1 {n1:,} intervals | v2 {n2:,} intervals")
    keys = sorted(set(b1) | set(b2), key=lambda k: -(b1.get(k, 0) + b2.get(k, 0)))
    print(f"  {'label':<16}{'v1 bp':>14}{'v2 bp':>14}{'diff':>14}{'v1-v2 %pt':>11}")
    t1, t2 = sum(b1.values()), sum(b2.values())
    for k in keys[:9]:
        x, y = b1.get(k, 0), b2.get(k, 0)
        print(f"  {k:<16}{x:>14,}{y:>14,}{x-y:>14,}{100*x/t1 - 100*y/t2:>10.2f}")
    ident = b1 == b2
    print(f"  IDENTICAL: {ident}")

print()
print("=" * 96)
print("PART 2 -- published method (100 bp, plateau_center, published focus window) on v1")
print("=" * 96)
print(f"{'sample':9}{'v1 center':>13}{'Δ SST1':>11}{'v2 center':>13}{'Δ SST1':>11}{'v1 99% CI':>22}{'width':>8}{'pub w':>7}")
PUB = {s: d["published_ci_width_mb"] for s, d in samples().items()}
for s, scaf, l, r, mid, i0, i1 in ROBS:
    fs, fe = i0 - 1_000_000, i1 + 1_000_000
    c1, a1, b1_ = published_method(V1.format(s=s), scaf, l, r, fs, fe)
    c2, _, _ = published_method(V2.format(s=s), scaf, l, r, fs, fe)
    print(f"{s:9}{c1:>13,}{(c1-mid)/1000:>10.1f}k{c2:>13,}{(c2-mid)/1000:>10.1f}k"
          f"{a1/1e6:>13.3f}-{b1_/1e6:<8.3f}{(b1_-a1)/1e6:7.2f}{PUB[s]:7.2f}")
