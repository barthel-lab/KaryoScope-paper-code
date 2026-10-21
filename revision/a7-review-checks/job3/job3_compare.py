#!/usr/bin/env python3
"""Compare forward vs reverse-complement annotation per base.

Coordinates: a KaryoScope BED labels k-mer START positions, so a sequence of length L tiles
[0, L-k+1). k-mer index p on the RC strand is k-mer index (L' - 1 - p) on the forward strand,
L' = L-k+1, so an RC interval [s, e) maps to forward [L'-e, L'-s). The presmoothed stage is the
control: with both strands in the index the lookup is orientation-free, so any difference there
would be a mapping error, not a finding.

The third sequence, the Robertsonian fusion contig, is read from the a7_job3 section of $ROBS_CONFIG
(file stem and description); with ROBS_CONFIG unset only the two HG002 sequences are compared."""
import argparse, collections, os, statistics, sys
ap = argparse.ArgumentParser()
ap.add_argument("--annot-dir", required=True); ap.add_argument("--db", required=True)
ap.add_argument("--backend", required=True); ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()
SEQS = [("chr21_hap1", "HG002 hap1 chr21 (acrocentric)"), ("chr20_hap1", "HG002 hap1 chr20 (ordinary autosome)")]
if os.environ.get("ROBS_CONFIG"):
    import yaml
    _rob = yaml.safe_load(open(os.environ["ROBS_CONFIG"]))["a7_job3"]
    SEQS.append((_rob["stem"], _rob["description"]))
SETS = ["chromosome", "region", "repeat", "subtelomeric", "gene", "acrocentric"]
BUCKETS = [(10, "1-10"), (31, "11-31"), (100, "32-100"), (1000, "101-1k"), (10000, "1k-10k"), (float("inf"), ">10k")]
def bucket(n):
    for lim, name in BUCKETS:
        if n <= lim: return name
def read(path):
    iv = []
    for line in open(path):
        c, s, e, lab = line.rstrip("\n").split("\t")[:4]
        iv.append((int(s), int(e), lab))
    return iv
summary, pairs, hist = [], [], []
for stem, desc in SEQS:
    for fs in SETS:
        for stage in ("presmoothed", "smoothed"):
            fwd_f = f"{a.annot_dir}/{stem}.fwd.{a.db}.{fs}.{stage}.bed"; rc_f = f"{a.annot_dir}/{stem}.rc.{a.db}.{fs}.{stage}.bed"
            if not (os.path.exists(fwd_f) and os.path.exists(rc_f)):
                print(f"MISSING {fwd_f} or {rc_f}"); continue
            F = read(fwd_f); R = read(rc_f)
            Lf = F[-1][1]; Lr = R[-1][1]
            assert F[0][0] == 0 and R[0][0] == 0 and Lf == Lr, (stem, fs, stage, Lf, Lr)
            # tiling check
            for X in (F, R):
                for i in range(1, len(X)): assert X[i][0] == X[i-1][1], (stem, fs, stage, X[i-1], X[i])
            M = [(Lf - e, Lf - s, lab) for s, e, lab in reversed(R)]  # RC mapped to forward, ascending
            same = diff = 0; pair = collections.Counter(); runs = []; cur = 0
            i = j = 0; pos = 0
            while pos < Lf:
                fs_, fe, fl = F[i]; ms, me, ml = M[j]
                end = min(fe, me); seg = end - pos
                if fl == ml:
                    same += seg
                    if cur: runs.append(cur); cur = 0
                else:
                    diff += seg; pair[(fl, ml)] += seg; cur += seg
                pos = end
                if fe == pos: i += 1
                if me == pos: j += 1
            if cur: runs.append(cur)
            assert same + diff == Lf
            summary.append((a.backend, stem, desc, Lf, fs, stage, Lf, same, diff, 100*diff/Lf, len(runs),
                            min(runs) if runs else 0, statistics.median(runs) if runs else 0, max(runs) if runs else 0, len(pair)))
            for (fl, ml), bp in sorted(pair.items(), key=lambda x: -x[1]):
                pairs.append((a.backend, stem, fs, stage, fl, ml, bp, 100*bp/Lf))
            hb = collections.Counter(); hn = collections.Counter()
            for r in runs: hb[bucket(r)] += r; hn[bucket(r)] += 1
            for _, name in BUCKETS:
                hist.append((a.backend, stem, fs, stage, name, hn[name], hb[name]))
            print(f"{stem:20} {fs:13} {stage:12} L'={Lf:>11,} identical={same:>11,} differing={diff:>9,} ({100*diff/Lf:.5f}%) runs={len(runs):,}")
def w(path, hdr, rows):
    with open(path, "w") as o:
        o.write("\t".join(hdr) + "\n")
        for r in rows: o.write("\t".join(f"{x:.5f}" if isinstance(x, float) else str(x) for x in r) + "\n")
w(a.out_prefix + ".tsv", ["backend", "sequence", "description", "sequence_kmer_positions", "feature_set", "stage", "bp_compared", "bp_identical", "bp_differing", "pct_differing", "n_differing_runs", "run_min_bp", "run_median_bp", "run_max_bp", "n_label_pairs"], summary)
w(a.out_prefix + "_label_pairs.tsv", ["backend", "sequence", "feature_set", "stage", "forward_label", "rc_label_mapped_back", "bp", "pct_of_sequence"], pairs)
w(a.out_prefix + "_run_hist.tsv", ["backend", "sequence", "feature_set", "stage", "run_length_bucket", "n_runs", "bp"], hist)
