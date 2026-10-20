#!/usr/bin/env python3
"""A2 Step 4 / Job 1 -- novel fraction and novel-run lengths at one (assembly, k).

Reads the EXISTING A2 Step 2 annotation BEDs (presmoothed + smoothed, six feature
sets) for one assembly at one k. No re-annotation. Method mirrors
a4_smoothing/run_a4_novel_v3.sh so the k=31 row must reproduce its anchor
(HG002 hap1: 75,396,753 novel bp of 2,940,260,101 evaluated = 2.5643%):

  * N runs are masked out entirely (absence, not novelty), extended k-1 bp to the
    LEFT exactly as fill_gaps.py does for the k-mer that spans into the N run.
  * chrM is excluded (the database excludes it).
  * evaluated_bp = every remaining annotated base.
  * a novel RUN = maximal set of coordinate-adjacent 'novel' intervals on one contig.

Per (feature_set, stage) it reports novel_bp, evaluated_bp, novel_frac, and the
fraction of novel bp in runs <= k, <= 31 (fixed, for cross-k comparability),
> 1 kb and > 10 kb, plus run counts and the longest run.
"""
import argparse, collections, gzip, os, statistics, sys

ap = argparse.ArgumentParser()
ap.add_argument("--annot-dir", required=True)
ap.add_argument("--assembly", required=True)
ap.add_argument("--k", type=int, required=True)
ap.add_argument("--nbed", default="", help="N-region BED for the assembly ('' if gapless)")
ap.add_argument("--exclude", default="chrM", help="comma-separated contigs to drop")
ap.add_argument("--out", required=True)
a = ap.parse_args()

SETS = ("chromosome", "region", "repeat", "subtelomeric", "gene", "acrocentric")
EXCL = set(x for x in a.exclude.split(",") if x)
ext = max(0, a.k - 1)
N = collections.defaultdict(list)
if a.nbed:
    for line in open(a.nbed):
        c, s, e = line.split("\t")[:3]
        N[c].append((max(0, int(s) - ext), int(e)))
    for c in N:
        N[c].sort()

def n_overlap(c, s, e):
    tot = 0
    for x, y in N.get(c, ()):
        lo, hi = max(s, x), min(e, y)
        if hi > lo:
            tot += hi - lo
    return tot

def scan(path):
    op = gzip.open if path.endswith(".gz") else open
    tot = 0; novel = 0
    runs = []
    cur_c = None; cur_end = None; cur_len = 0
    def close():
        nonlocal cur_len
        if cur_len > 0:
            runs.append(cur_len)
        cur_len = 0
    with op(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            c = p[0]
            if c in EXCL:
                close(); cur_c = None; continue
            s = int(p[1]); e = int(p[2]); L = e - s
            eff = L - (n_overlap(c, s, e) if c in N else 0)
            if eff <= 0:
                close(); cur_c = None; continue
            tot += eff
            if p[3] == "novel":
                novel += eff
                if cur_c == c and cur_end == s:
                    cur_len += eff
                else:
                    close(); cur_len = eff
                cur_c = c; cur_end = e
            else:
                close(); cur_c = None
    close()
    return tot, novel, runs

def bucket_stats(runs, k, novel):
    d = {}
    def frac(pred):
        return sum(l for l in runs if pred(l)) / novel if novel else 0.0
    d["frac_runs_le_k"] = frac(lambda l: l <= k)
    d["frac_runs_le_31"] = frac(lambda l: l <= 31)
    d["frac_runs_32_100"] = frac(lambda l: 31 < l <= 100)
    d["frac_runs_101_1k"] = frac(lambda l: 100 < l <= 1000)
    d["frac_runs_gt_1kb"] = frac(lambda l: l > 1000)
    d["frac_runs_gt_10kb"] = frac(lambda l: l > 10000)
    d["n_runs"] = len(runs)
    d["median_run"] = statistics.median(runs) if runs else 0
    d["max_run"] = max(runs) if runs else 0
    return d

cols = ["assembly", "k", "feature_set", "stage", "novel_bp", "evaluated_bp", "novel_frac",
        "frac_runs_le_k", "frac_runs_gt_1kb", "frac_runs_le_31", "frac_runs_32_100",
        "frac_runs_101_1k", "frac_runs_gt_10kb", "n_runs", "median_run", "max_run", "source_file"]
rows = []
for fs in SETS:
    for stage in ("presmoothed", "smoothed"):
        cands = [f for f in os.listdir(a.annot_dir) if f.endswith(f".{fs}.{stage}.bed") or f.endswith(f".{fs}.{stage}.bed.gz")]
        if len(cands) != 1:
            print(f"ERROR: expected exactly one {fs}.{stage} BED in {a.annot_dir}, found {cands}", file=sys.stderr)
            sys.exit(1)
        path = os.path.join(a.annot_dir, cands[0])
        tot, novel, runs = scan(path)
        b = bucket_stats(runs, a.k, novel)
        frac = novel / tot if tot else 0.0
        rows.append([a.assembly, a.k, fs, stage, novel, tot, f"{frac:.8f}",
                     f"{b['frac_runs_le_k']:.6f}", f"{b['frac_runs_gt_1kb']:.6f}", f"{b['frac_runs_le_31']:.6f}",
                     f"{b['frac_runs_32_100']:.6f}", f"{b['frac_runs_101_1k']:.6f}", f"{b['frac_runs_gt_10kb']:.6f}",
                     b["n_runs"], b["median_run"], b["max_run"], cands[0]])
        print(f"{a.assembly} k={a.k} {fs:13s} {stage:12s} novel={novel:>12,} eval={tot:>14,} "
              f"frac={100*frac:.4f}%  <=k {100*b['frac_runs_le_k']:.2f}%  >1kb {100*b['frac_runs_gt_1kb']:.3f}%", flush=True)

with open(a.out, "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print("wrote", a.out)
