#!/usr/bin/env python3
"""Per-label bp per stage for one feature set, A4 N-mask convention (N runs +30 bp left; chrM dropped)."""
import os as _ks_os
KS_SOFTWARE = _ks_os.environ.get("KS_SOFTWARE", "$KS_SOFTWARE")  # see the repository README, Paths
import sys, os, argparse, collections
sys.path.insert(0, os.path.dirname(__file__))
from common_mask import *
ap = argparse.ArgumentParser()
ap.add_argument("--pre"); ap.add_argument("--post"); ap.add_argument("--nbed"); ap.add_argument("--db"); ap.add_argument("--out")
ap.add_argument("--assembly", default="hg002hap1"); ap.add_argument("--feature-set", default="gene")
a = ap.parse_args()
HIER = f"{KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/{a.db}/hierarchy.tsv"
if not os.path.exists(HIER): HIER = a.db  # allow a path
nodes, leaves, parent_of = parse_hierarchy(HIER)[a.feature_set]
lvl = make_level(nodes, leaves)
mask = n_runs_from_bed(a.nbed) if a.nbed else {}
def tally(path):
    assert ("presmoothed" in os.path.basename(path)) == (path == a.pre), path
    by = collections.Counter(); den = 0; masked = 0
    with opener(path) as fh:
        for line in fh:
            c, s, e, lab = line.rstrip("\n").split("\t")[:4]
            if c == "chrM": continue
            s, e = int(s), int(e)
            m = overlap(mask.get(c, ()), s, e); masked += m
            eff = (e - s) - m
            if eff <= 0: continue
            by[lab] += eff; den += eff
    return by, den, masked
pre, den1, m1 = tally(a.pre); post, den2, m2 = tally(a.post)
assert den1 == den2, (den1, den2)
with open(a.out, "w") as o:
    o.write("assembly\tfeature_set\tstage\tlabel\tlevel\tbp\tpct\tdenominator_bp\n")
    for stage, d in (("presmoothed", pre), ("smoothed", post)):
        for lab in sorted(d, key=lambda l: -d[l]):
            o.write(f"{a.assembly}\t{a.feature_set}\t{stage}\t{lab}\t{lvl(lab)}\t{d[lab]}\t{100*d[lab]/den1:.4f}\t{den1}\n")
print(f"denominator {den1:,} (N+chrM masked pre {m1:,} / post {m2:,})")
for stage, d in (("presmoothed", pre), ("smoothed", post)):
    print(stage, "  ".join(f"{l}={d[l]:,} ({100*d[l]/den1:.4f}%)" for l in sorted(d, key=lambda l: -d[l])))
