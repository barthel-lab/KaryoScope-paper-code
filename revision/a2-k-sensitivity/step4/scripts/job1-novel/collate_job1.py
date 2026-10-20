#!/usr/bin/env python3
"""Collate Job 1 per-(assembly,k) TSVs into novel_by_k.tsv (brief columns first)."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import glob, os
R = (KS_WORK + "/a2_step4/job1_novel/results")
files = sorted(glob.glob(f"{R}/novel_*_k*.tsv"))
rows = []; hdr = None
for f in files:
    with open(f) as fh:
        h = next(fh).rstrip("\n").split("\t")
        hdr = hdr or h
        for line in fh:
            rows.append(line.rstrip("\n").split("\t"))
korder = {"11": 0, "21": 1, "31": 2, "41": 3, "51": 4, "61": 5}
sets = ["chromosome", "region", "repeat", "subtelomeric", "gene", "acrocentric"]
rows.sort(key=lambda r: (r[0], korder[r[1]], sets.index(r[2]), r[3]))
out = (KS_WORK + "/a2_step4/novel_by_k.tsv")
with open(out, "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows: o.write("\t".join(r) + "\n")
print(f"wrote {out}: {len(rows)} rows from {len(files)} files\n")
# headline table: repeat set (novel is identical across sets pre-smoothing; post differs)
print(f"{'assembly':10} {'k':>3} {'stage':12} {'novel_bp':>12} {'evaluated_bp':>14} {'novel%':>8} {'<=k%':>7} {'<=31%':>7} {'>1kb%':>7} {'>10kb%':>7} {'n_runs':>9} {'med':>5} {'max':>7}")
for r in rows:
    if r[2] != "repeat": continue
    d = dict(zip(hdr, r))
    print(f"{d['assembly']:10} {d['k']:>3} {d['stage']:12} {int(d['novel_bp']):>12,} {int(d['evaluated_bp']):>14,} "
          f"{100*float(d['novel_frac']):>8.4f} {100*float(d['frac_runs_le_k']):>7.2f} {100*float(d['frac_runs_le_31']):>7.2f} "
          f"{100*float(d['frac_runs_gt_1kb']):>7.3f} {100*float(d['frac_runs_gt_10kb']):>7.3f} {int(d['n_runs']):>9,} {float(d['median_run']):>5.0f} {int(d['max_run']):>7,}")
# is pre-smoothing novel identical across the six sets at every k?
import collections
pre = collections.defaultdict(set)
for r in rows:
    if r[3] == "presmoothed": pre[(r[0], r[1])].add(r[4])
bad = {k: v for k, v in pre.items() if len(v) > 1}
print("\npre-smoothing novel_bp identical across the six sets:", "YES at every (assembly,k)" if not bad else f"NO: {bad}")
