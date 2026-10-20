#!/usr/bin/env python3
"""Where do Ey15-2's presmoothed-novel bases go after smoothing? Per feature set: bp of presmoothed-novel
that is smoothed to each label (leaf / internal / still novel). Same for presmoothed-multilabel bases."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import sys, os, collections, gzip
sys.path.insert(0, os.path.dirname(__file__))
from common_mask import *
A3 = (KS_WORK + "/a3_crossgenome"); RES = (KS_WORK + "/a7_review_checks/results")
H = parse_hierarchy(f"{A3}/db/HKS_arabidopsis_ColCEN/hierarchy.tsv")
def read(f):
    out = collections.defaultdict(list)
    with opener(f) as fh:
        for line in fh:
            c, s, e, lab = line.rstrip("\n").split("\t")[:4]; out[c].append((int(s), int(e), lab))
    return out
rows = []
for fs in ("chromosome", "region", "repeat", "gene"):
    nodes, leaves, _ = H[fs]; lvl = make_level(nodes, leaves)
    P = read(f"{A3}/annot_ey15_2/ey15-2.HKS_arabidopsis_ColCEN.{fs}.presmoothed.bed.gz")
    S = read(f"{A3}/annot_ey15_2/ey15-2.HKS_arabidopsis_ColCEN.{fs}.smoothed.bed.gz")
    fate = collections.defaultdict(collections.Counter)  # pre level -> post label
    for c in P:
        A, B = P[c], S[c]; i = j = 0; pos = 0; L = A[-1][1]
        while pos < L:
            end = min(A[i][1], B[j][1]); seg = end - pos
            fate[lvl(A[i][2])][B[j][2]] += seg
            pos = end
            if A[i][1] == pos: i += 1
            if B[j][1] == pos: j += 1
    for pre in ("novel", "internal", "leaf"):
        tot = sum(fate[pre].values())
        agg = collections.Counter()
        for lab, bp in fate[pre].items(): agg[lvl(lab)] += bp
        for lab, bp in fate[pre].most_common():
            rows.append((fs, pre, tot, lab, lvl(lab), bp, 100*bp/tot if tot else 0))
        print(f"{fs:11} presmoothed {pre:8} {tot:>11,} bp -> smoothed leaf {100*agg['leaf']/tot:6.2f}%  internal {100*agg['internal']/tot:6.2f}%  novel {100*agg['novel']/tot:6.2f}%")
with open(f"{RES}/ey15_2_novel_fate.tsv", "w") as o:
    o.write("feature_set\tpresmoothed_level\tpresmoothed_level_bp\tsmoothed_label\tsmoothed_level\tbp\tpct_of_presmoothed_level\n")
    for r in rows: o.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
