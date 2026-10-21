#!/usr/bin/env python3
"""Pre-release (gap-masked) vs released HPRC Release 2 QC calls, centromere region only.

For every (sample, hap, chrom) centromere region in --regions, count the calls of each tool that overlap the
region (bedtools-intersect semantics: >=1 bp overlap, one record per call, bp clipped to the region), sum the
nine flags across tools (as aggregate_qc_v2 does), apply Err=0, COLLAPSE=0, COLLAPSE_VAR=0, and compare the
pass sets.

Pre-release: the files of ../prerelease_qc_index.tsv under --qc-dir (<qc-dir>/<tool>/<file>).
Released:    <released-dir>/{sample}.flagger_{hifi,ont}.bed (released_flagger_index.tsv) + the same NucFlag
             files (the released NucFlag BEDs are byte-identical to the pre-release ones, by sha256).
Validation:  the pre-release recomputation must equal every centromere row of --curated (all nine counts).
Pass counts use the curated keys present in --clustering-bed, as QCfilter_apply.R counts them.

Writes <out-dir>/comparison.tsv and <out-dir>/passfail_diff.tsv; prints a summary to stderr.
"""
import argparse
import collections, csv, os, sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--qc-index", required=True, help="../prerelease_qc_index.tsv")
ap.add_argument("--qc-dir", required=True, help="pre-release QC BEDs as <qc-dir>/<tool>/<file>")
ap.add_argument("--released-dir", required=True, help="released Flagger BEDs as {sample}.flagger_{hifi,ont}.bed")
ap.add_argument("--regions", required=True, help="../inputs/regions.tsv")
ap.add_argument("--curated", required=True, help="the curated QC table")
ap.add_argument("--clustering-bed", required=True,
                help="pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed (6,138 contigs)")
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()
CUR, CLUST, OUT = args.curated, args.clustering_bed, args.out_dir
QC = collections.defaultdict(list)       # (sample, tool) -> [path, ...]
with open(args.qc_index) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        QC[(r["sample"], r["tool"])].append(os.path.join(args.qc_dir, r["tool"], r["file"]))
REG = collections.defaultdict(list)      # (sample, hap) -> [(chrom, contig, start, stop), ...]
with open(args.regions) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if r["region"] == "centromere":
            REG[(r["sample"], r["haplotype"])].append((r["chromosome"], r["contig"], r["start"], r["stop"]))
FLAGS = "Col Dup Err NNN COLLAPSE COLLAPSE_OTHER COLLAPSE_VAR MISJOIN HET".split()
TOOLS = ["flagger_hifi", "flagger_ont", "nucflag"]


def read_calls(paths, tool, into):
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith("track") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                into[c[0]].append((int(c[1]), int(c[2]), c[3], tool))


# curated table (centromere rows)
cur = {}
with open(CUR) as fh:
    r = csv.reader(fh, delimiter="\t")
    hdr = next(r)
    for c in r:
        if c[3] == "centromere":
            cur[(c[0], c[1], c[2], c[4])] = dict(zip(FLAGS, map(int, c[5:14])))
samples = sorted({k[0] for k in cur})

clust_keys = set()
with open(CLUST) as fh:
    for line in fh:
        c = line.split("\t")
        clust_keys.add((c[0], c[4].strip()))

rows = []          # per region results
calls_out = collections.Counter()   # (set, tool, flag) -> n ; bp
bp_out = collections.Counter()
missing = []
for s in samples:
    pre = collections.defaultdict(list)
    rel = collections.defaultdict(list)
    read_calls(QC[(s, "flagger_hifi")], "flagger_hifi", pre)
    read_calls(QC[(s, "flagger_ont")], "flagger_ont", pre)
    nuc = QC[(s, "nucflag")]
    read_calls(nuc, "nucflag", pre)
    for t in ("hifi", "ont"):
        p = os.path.join(args.released_dir, f"{s}.flagger_{t}.bed")
        if os.path.exists(p):
            read_calls([p], "flagger_" + t, rel)
        else:
            missing.append((s, "flagger_" + t))
    read_calls(nuc, "nucflag", rel)
    for hap in ("hap1", "hap2"):
        # chromosomes in string order, one interval each (the order the region BEDs were globbed in)
        for chrom, contig, rs, re_ in sorted(REG[(s, hap)]):
            regs = [(contig, rs, re_)]
            for contig, rs, re_ in regs:
                rs, re_ = int(rs), int(re_)
                res = {}
                for name, calls in (("pre", pre), ("rel", rel)):
                    cnt = collections.Counter()
                    for (a, b, f, t) in calls.get(contig, ()):
                        if a < re_ and b > rs:
                            cnt[f] += 1
                            cnt[(t, f)] += 1
                            calls_out[(name, t, f)] += 1
                            bp_out[(name, t, f)] += min(b, re_) - max(a, rs)
                    res[name] = cnt
                rows.append((s, hap, chrom, contig, rs, re_, res["pre"], res["rel"]))

# curated rows whose region BED is empty (no interval): the master table wrote '.' for every tool, so the
# same (empty) region gives zero calls in both call sets.
have = {(r[0], r[1], r[2], r[3]) for r in rows}
for k in sorted(set(cur) - have):
    print("  empty region BED, zero counts in both sets:", k, file=sys.stderr)
    rows.append((k[0], k[1], k[2], k[3], -1, -1, collections.Counter(), collections.Counter()))


def passes(c):
    return c["Err"] == 0 and c["COLLAPSE"] == 0 and c["COLLAPSE_VAR"] == 0

# validation against curated table
match = mism = notcur = 0
mism_ex = []
seen = set()
for s, hap, chrom, contig, rs, re_, p, r in rows:
    k = (s, hap, chrom, contig)
    if k not in cur:
        notcur += 1
        continue
    seen.add(k)
    if all(p[f] == cur[k][f] for f in FLAGS):
        match += 1
    else:
        mism += 1
        mism_ex.append((k, {f: (p[f], cur[k][f]) for f in FLAGS if p[f] != cur[k][f]}))
print(f"regions computed: {len(rows)}; match curated: {match}; mismatch: {mism}; not in curated: {notcur}; "
      f"curated rows without region: {len(set(cur) - seen)}", file=sys.stderr)
for e in mism_ex[:20]:
    print("  MISMATCH", e, file=sys.stderr)
print("missing released files:", missing, file=sys.stderr)

# pass counts (count = curated keys present in clustering input, as QCfilter_apply.R does)
cur_pass = {k for k, v in cur.items() if passes(v) and (k[3], k[2]) in clust_keys}
pre_pass = {(s, h, c, ct) for s, h, c, ct, _, _, p, r in rows if passes(p) and (ct, c) in clust_keys
            and (s, h, c, ct) in cur}
rel_pass = {(s, h, c, ct) for s, h, c, ct, _, _, p, r in rows if passes(r) and (ct, c) in clust_keys
            and (s, h, c, ct) in cur}
print(f"pass: curated={len(cur_pass)} pre_recomputed={len(pre_pass)} released={len(rel_pass)}", file=sys.stderr)

with open(os.path.join(OUT, "comparison.tsv"), "w") as out:
    out.write("tool\tflag\tn_calls_prerelease\tn_calls_released\tdelta_n\tbp_prerelease\tbp_released\tdelta_bp\t"
              "n_regions_count_changed\n")
    for t in TOOLS:
        for f in FLAGS + ["Hap"]:
            a, b = calls_out[("pre", t, f)], calls_out[("rel", t, f)]
            if a == 0 and b == 0:
                continue
            nchg = sum(1 for row in rows if row[6][(t, f)] != row[7][(t, f)])
            out.write(f"{t}\t{f}\t{a}\t{b}\t{b-a}\t{bp_out[('pre',t,f)]}\t{bp_out[('rel',t,f)]}\t"
                      f"{bp_out[('rel',t,f)]-bp_out[('pre',t,f)]}\t{nchg}\n")
    # merged (tools summed) per-flag region-level change
    out.write("ALL_TOOLS_MERGED\tregions_with_any_count_change\t\t\t"
              f"{sum(1 for row in rows if any(row[6][f] != row[7][f] for f in FLAGS))}\t\t\t\t\n")
    for f in FLAGS:
        out.write(f"ALL_TOOLS_MERGED\t{f}\t{sum(r[6][f] for r in rows)}\t{sum(r[7][f] for r in rows)}\t"
                  f"{sum(r[7][f]-r[6][f] for r in rows)}\t\t\t\t{sum(1 for r in rows if r[6][f]!=r[7][f])}\n")

with open(os.path.join(OUT, "passfail_diff.tsv"), "w") as out:
    out.write("sample\thaplotype\tchromosome\tcontig\tregion_start\tregion_end\tpre_Err\tpre_COLLAPSE\t"
              "pre_COLLAPSE_VAR\tpre_NNN\trel_Err\trel_COLLAPSE\trel_COLLAPSE_VAR\tpre_pass\trel_pass\tin_curated_pass\n")
    for s, h, c, ct, rs, re_, p, r in rows:
        k = (s, h, c, ct)
        if k not in cur or (ct, c) not in clust_keys:
            continue
        if passes(p) != passes(r):
            out.write(f"{s}\t{h}\t{c}\t{ct}\t{rs}\t{re_}\t{p['Err']}\t{p['COLLAPSE']}\t{p['COLLAPSE_VAR']}\t"
                      f"{p['NNN']}\t{r['Err']}\t{r['COLLAPSE']}\t{r['COLLAPSE_VAR']}\t{passes(p)}\t{passes(r)}\t"
                      f"{k in cur_pass}\n")

# per-tool Err region change detail for the report
chg = collections.Counter()
for s, h, c, ct, rs, re_, p, r in rows:
    for f in FLAGS:
        if p[f] != r[f]:
            chg[f] += 1
print("regions with merged-count change by flag:", dict(chg), file=sys.stderr)
