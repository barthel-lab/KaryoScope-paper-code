#!/usr/bin/env python3
"""Collate the A6 centromere-clustering sweep into one TSV.

Parses the per-chromosome Major/Outlier counts KS_allchr_barplot.py prints, so
Yi-An's scripts run unmodified. Reported quantities, per (sil_threshold,
centroid_sd, centroid_scan):

  n_chroms_split   chromosomes with at least one Outlier -- the "retains a
                   split" count, which is the published major/minor structure
  n_outlier        total Outlier haplotypes
  n_major          total Major haplotypes
  per-chromosome Outlier counts, so a chromosome that gains or loses its split
                   can be named rather than just counted

The --no-centroid-scan rows isolate the silhouette stage: both stages create and
destroy outliers, so a combined count cannot attribute a change to one of them.
"""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import re
from pathlib import Path

W = Path((KS_WORK + "/a6_sweeps/cent/sweep"))
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
ROW = re.compile(r"^\s*(chr[0-9XY]+):\s*(\d+)\s+Major,\s*(\d+)\s+Outlier")
TAG = re.compile(r"^sil([0-9.]+)_sd(\d+)_(on|off)\.stdout$")

rows = []
for f in sorted(W.glob("*.stdout")):
    m = TAG.match(f.name)
    if not m:
        continue
    sil, sd, scan = float(m.group(1)), int(m.group(2)), m.group(3)
    counts = {}
    for line in f.read_text().splitlines():
        r = ROW.match(line)
        if r:
            counts[r.group(1)] = (int(r.group(2)), int(r.group(3)))
    if not counts:
        print(f"WARNING: no counts parsed from {f.name} -- run failed or output changed")
        continue
    n_major = sum(v[0] for v in counts.values())
    n_out = sum(v[1] for v in counts.values())
    split = [c for c in CHROMS if counts.get(c, (0, 0))[1] > 0]
    rows.append(dict(sil=sil, sd=sd, scan=scan, n_chroms_split=len(split),
                     n_major=n_major, n_outlier=n_out,
                     chroms_split=",".join(split),
                     per_chrom={c: counts.get(c, (0, 0))[1] for c in CHROMS}))

rows.sort(key=lambda r: (r["scan"], r["sil"], r["sd"]))
cols = ["sil_threshold", "centroid_sd", "centroid_scan", "n_chroms_split",
        "n_major", "n_outlier", "chroms_split"] + [f"out_{c}" for c in CHROMS]
out = ["\t".join(cols)]
for r in rows:
    out.append("\t".join(str(x) for x in
        [r["sil"], r["sd"], r["scan"], r["n_chroms_split"], r["n_major"],
         r["n_outlier"], r["chroms_split"]] + [r["per_chrom"][c] for c in CHROMS]))

dest = W.parent / "cent_param_sweep.tsv"
dest.write_text("\n".join(out) + "\n")
print(f"wrote {dest} ({len(rows)} rows)\n")

print("=== centroid scan ON (the published configuration) ===")
print(f"{'sil':>5} {'sd':>3} | {'chroms split':>12} {'major':>7} {'outlier':>8}")
for r in rows:
    if r["scan"] == "on":
        star = "  <- published" if (r["sil"] == 0.5 and r["sd"] == 5) else ""
        print(f"{r['sil']:>5} {r['sd']:>3} | {r['n_chroms_split']:>12} "
              f"{r['n_major']:>7} {r['n_outlier']:>8}{star}")
print("\n=== silhouette stage only (centroid scan OFF) ===")
for r in rows:
    if r["scan"] == "off":
        print(f"{r['sil']:>5} {'--':>3} | {r['n_chroms_split']:>12} "
              f"{r['n_major']:>7} {r['n_outlier']:>8}")
