#!/usr/bin/env python3
"""Does the 99% CI contain the SST1 block midpoint -- and is the published
CI-width gap explained by presmoothed-vs-smoothed input?

Runs the published method (100 bp, plateau_center, CI = combined >= 0.99*max,
focus window = `_robertsonian` ISCN CI +/- 1 Mb) over BOTH the smoothed
scaffolded chromosome BED (what the current pipeline reads) and the presmoothed
scaffolded one (what the legacy breakpoint_density.py read), and reports where
the SST1 midpoint falls relative to the interval.

All coordinates SCAFFOLDED frame. Per-sample values come from the ROB
configuration ($ROBS_CONFIG, section a2_step3).
"""
import gzip, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rob_samples import samples  # noqa: E402  per-sample values from $ROBS_CONFIG (a2_step3)

RES = 100
W = os.path.join(os.environ.get("KS_WORK", "$KS_WORK"), "a2_step3")
ROBS = [(s, d["scaffold"], d["left"], d["right"], d["sst1_mid"], d["win_start"], d["win_end"],
         d["published_ci_width_mb"]) for s, d in samples().items()]
SMOOTH = os.path.join(W, "gate_{s}/work/"
                      "{s}.haplotype1.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz")
PRESM  = os.path.join(W, "presm_{s}/"
                      "{s}.haplotype1.HKS_human_CHM13_v2.chromosome.presmoothed.scaffolded.bed.gz")


def analyse(path, scaffold, left, right, fs, fe):
    n = int((fe - fs) // RES)
    L = np.zeros(n); R = np.zeros(n)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != scaffold:
                continue
            lab = f[3]
            if lab != left and lab != right:
                continue
            s = int(f[1]); e = int(f[2])
            if e <= fs or s >= fe:
                continue
            i0 = max(0, (s - fs) // RES); i1 = min(n, -(-(e - fs) // RES))
            tgt = L if lab == left else R
            for i in range(i0, i1):
                bs = fs + i * RES
                tgt[i] += max(0, min(e, bs + RES) - max(s, bs))
    pos = fs + np.arange(n) * RES
    lc = np.cumsum(L); rc = np.cumsum(R[::-1])[::-1]
    comb = lc / (lc[-1] or 1) + rc / (rc[0] or 1)
    mx = comb.max()
    at = np.where(np.isclose(comb, mx, rtol=1e-10))[0]
    ci = np.where(comb >= mx * 0.99)[0]
    return (int(pos[at[0]]) + int(pos[at[-1]])) // 2, int(pos[ci[0]]), int(pos[ci[-1]])


def main():
    print(f"{'sample':8} {'input':11} {'center':>12} {'d(SST1)':>10} "
          f"{'99% CI (Mb)':>20} {'width':>7} {'pub w':>6} {'SST1 in CI?':>12} {'d to CI':>9}")
    for s, scaf, l, r, mid, i0, i1, pubw in ROBS:
        fs, fe = i0 - 1_000_000, i1 + 1_000_000
        for tag, tmpl in (("smoothed", SMOOTH), ("presmoothed", PRESM)):
            path = tmpl.format(s=s)
            try:
                c, a, b = analyse(path, scaf, l, r, fs, fe)
            except (FileNotFoundError, EOFError, OSError):
                print(f"{s:8} {tag:11}  (not ready)")
                continue
            inside = a <= mid <= b
            d = 0 if inside else (a - mid if mid < a else mid - b)
            print(f"{s:8} {tag:11} {c:12,} {(c-mid)/1000:+9.1f}k "
                  f"{a/1e6:9.3f}-{b/1e6:<10.3f}{(b-a)/1e6:7.2f} {pubw:6.2f} "
                  f"{('YES' if inside else 'no'):>12} {d/1000:8.1f}k")


main()
