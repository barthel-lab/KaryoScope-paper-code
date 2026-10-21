#!/usr/bin/env python3
"""Gap-containing centromeres in the analysed 4,551, gap vs Outlier, and pre-release NNN vs released Err.

(1) NNN > 0 in the centromere region (pre-release Flagger, either tool; regions of ../inputs/regions.tsv) for
    the chromosomes in --clustered.
(2) 2x2 NNN>0 x Major/Outlier with Fisher exact p:
    (a) gated per-haplotype labels, --allele-pairs (allchr_allele_pairs.tsv, step 7b; chr3/8/11/12);
    (b) ungated stage-3 labels, --assignments (allchr_structure.sequence_assignments.tsv, step 3);
    (c) gated minor clusters (--gated, allchr_outlier_annotations.tsv, step 5) mapped to stage-3 Outlier
        clusters (the cluster whose own representative is the gated representative), used only where the
        stage-3 cluster size equals the gated n.
(3) bp of pre-release NNN (centromere regions) under released Err, per tool; and assembly N runs (--gaps, e.g.
    `seqtk gap -l 1` on the assemblies) in centromere regions: bp labelled NNN (pre) / under Err (released).

--clustered: TSV with columns sample, haplotype, chrom, contig, cluster, cluster_type: the step-3 assignments of
the 4,551 analysed chromosomes (the QC-passed contigs clustered in step 3, minus chr13, chr14, chr15, chr21,
chr22 and chrY). Steps refer to ../../scripts/HPRC2_centromere_analysis.sh.
Writes addendum.txt, addendum_regions.tsv, addendum_nnn_per_chrom.tsv, addendum_gap_vs_outlier.tsv to --out-dir.
"""
import argparse, os
import collections, sys
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--qc-index", required=True, help="../prerelease_qc_index.tsv")
ap.add_argument("--qc-dir", required=True, help="pre-release QC BEDs as <qc-dir>/<tool>/<file>")
ap.add_argument("--released-dir", required=True, help="released Flagger BEDs as {sample}.flagger_{hifi,ont}.bed")
ap.add_argument("--regions", required=True, help="../inputs/regions.tsv")
ap.add_argument("--clustered", required=True)
ap.add_argument("--gaps", required=True, help="BED of assembly N runs")
ap.add_argument("--assignments", required=True)
ap.add_argument("--allele-pairs", required=True)
ap.add_argument("--gated", required=True)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()
OUT = args.out_dir
QC = collections.defaultdict(list)       # (sample, hifi|ont) -> [path, ...]
with open(args.qc_index) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if r["tool"].startswith("flagger_"):
            QC[(r["sample"], r["tool"][len("flagger_"):])].append(os.path.join(args.qc_dir, r["tool"], r["file"]))
REG = {}                                 # (sample, hap, chrom) -> [(contig, start, stop)]
with open(args.regions) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        if r["region"] == "centromere":
            REG.setdefault((r["sample"], r["haplotype"], r["chromosome"]), []).append((r["contig"], r["start"], r["stop"]))
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
out = open(os.path.join(OUT, "addendum.txt"), "w")


def say(*a):
    print(*a); print(*a, file=out)


def load(paths, keep=None):
    d = collections.defaultdict(list)
    for p in paths:
        for line in open(p):
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or line.startswith("track"):
                continue
            if keep is None or c[3] in keep:
                d[(c[0], c[3])].append((int(c[1]), int(c[2])))
    return d


def union(iv):
    m = []
    for a, b in sorted(iv):
        if m and a <= m[-1][1]:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    return m


def clip(iv, s, e):
    return [(max(a, s), min(b, e)) for a, b in iv if a < e and b > s]


def bp(iv):
    return sum(b - a for a, b in iv)


def inter_bp(x, y):  # both unions
    i = j = t = 0
    while i < len(x) and j < len(y):
        a, b = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
        if a < b:
            t += b - a
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return t


cl = pd.read_csv(args.clustered, sep="\t")
say(f"clustered.tsv rows: {len(cl)}")
gaps = collections.defaultdict(list)
for line in open(args.gaps):
    c = line.split("\t")
    gaps[c[0]].append((int(c[1]), int(c[2])))

# region-level measures for every region in the 4,650 contig list (rows of gaps_telomeres.tsv) -- compute for the 4,551
rec = []
tot = collections.Counter()
for s, g in cl.groupby("sample"):
    pre = {t: load(QC[(s, t)], {"NNN", "Err"}) for t in ("hifi", "ont")}
    rel = {t: load([os.path.join(args.released_dir, f"{s}.flagger_{t}.bed")], {"Err"}) for t in ("hifi", "ont")}
    for r in g.itertuples():
        reg = REG.get((s, r.haplotype, r.chrom), [])
        assert len(reg) == 1 and reg[0][0] == r.contig, (s, r.haplotype, r.chrom)
        rs, re_ = int(reg[0][1]), int(reg[0][2])
        gp = union(clip(gaps.get(r.contig, []), rs, re_))
        row = dict(sample=s, haplotype=r.haplotype, chrom=r.chrom, contig=r.contig, cluster=r.cluster,
                   ungated=r.cluster_type, fasta_gap_n=len(gp), fasta_gap_bp=bp(gp))
        anynnn = False
        for t in ("hifi", "ont"):
            nraw = clip(pre[t].get((r.contig, "NNN"), []), rs, re_)
            n = union(nraw)
            e = union(clip(rel[t].get((r.contig, "Err"), []), rs, re_))
            row[f"{t}_NNN_n"] = len(nraw)
            row[f"{t}_NNN_bp"] = bp(n)
            anynnn |= len(nraw) > 0
            tot[(t, "NNN_bp")] += bp(n)
            tot[(t, "NNN_bp_under_relErr")] += inter_bp(n, e)
            tot[(t, "gap_bp")] += bp(gp)
            tot[(t, "gap_bp_labelled_NNN")] += inter_bp(gp, n)
            tot[(t, "gap_bp_under_relErr")] += inter_bp(gp, e)
            # gaps with no NNN overlap at all
            tot[(t, "gaps_without_NNN")] += sum(1 for a, b in gp if inter_bp([[a, b]], n) == 0)
            tot[(t, "gaps_without_relErr")] += sum(1 for a, b in gp if inter_bp([[a, b]], e) == 0)
        row["NNN_pos"] = anynnn
        rec.append(row)
d = pd.DataFrame(rec)
d.to_csv(os.path.join(OUT, "addendum_regions.tsv"), sep="\t", index=False)
say(f"joined regions: {len(d)}")

# (1)
say("\n(1) NNN>0 in centromere region, analysed 4,551")
say(f"overall: {int(d.NNN_pos.sum())} / {len(d)} = {100*d.NNN_pos.mean():.1f}%")
say(f"cross-check, any assembly N run in the region (seqtk): {int((d.fasta_gap_n>0).sum())}; "
    f"NNN>0 but no N run: {int((d.NNN_pos & (d.fasta_gap_n==0)).sum())}; N run but NNN=0: {int((~d.NNN_pos & (d.fasta_gap_n>0)).sum())}")
pc = d.groupby("chrom").agg(n=("NNN_pos", "size"), nnn=("NNN_pos", "sum"), gap_any=("fasta_gap_n", lambda x: int((x > 0).sum())))
pc = pc.reindex([c for c in CHROMS if c in pc.index])
pc["pct"] = (100 * pc.nnn / pc.n).round(1)
say(pc.to_string())
pc.to_csv(os.path.join(OUT, "addendum_nnn_per_chrom.tsv"), sep="\t")


def table(lbl, outlier_mask, sub, name):
    a = int((sub.NNN_pos & outlier_mask).sum()); b = int((sub.NNN_pos & ~outlier_mask).sum())
    c = int((~sub.NNN_pos & outlier_mask).sum()); e = int((~sub.NNN_pos & ~outlier_mask).sum())
    orr, p = fisher_exact([[a, b], [c, e]])
    return dict(set=name, chrom=lbl, NNN_Outlier=a, NNN_Major=b, noNNN_Outlier=c, noNNN_Major=e,
                pct_Outlier_if_NNN=round(100 * a / (a + b), 1) if a + b else float("nan"),
                pct_Outlier_if_noNNN=round(100 * c / (c + e), 1) if c + e else float("nan"),
                odds_ratio=round(orr, 3), fisher_p=p)


rows = []
# (2b) ungated
m = d.ungated == "Outlier"
rows.append(table("ALL", m, d, "ungated_stage3"))
for ch in pc.index:
    s = d[d.chrom == ch]
    rows.append(table(ch, s.ungated == "Outlier", s, "ungated_stage3"))
# (2a) gated per-haplotype for chr3/8/11/12
pairs = pd.read_csv(args.allele_pairs, sep="\t")
lab = {}
for r in pairs.itertuples():
    lab[(r.chromosome, r.sample, "hap1")] = r.h1_cluster
    lab[(r.chromosome, r.sample, "hap2")] = r.h2_cluster
d["gated"] = [lab.get((c, s, h)) for c, s, h in zip(d.chrom, d["sample"], d.haplotype)]
g4 = d[d.chrom.isin(["chr3", "chr8", "chr11", "chr12"])]
say(f"\n(2a) gated labels for chr3/8/11/12: {g4.gated.notna().sum()} of {len(g4)} haplotypes labelled; "
    f"Major/Outlier agree with ungated for {int(((g4.gated=='Major')==(g4.ungated=='Major'))[g4.gated.notna()].sum())}")
g4 = g4[g4.gated.notna()]
gm = g4.gated != "Major"
rows.append(table("chr3+8+11+12", gm, g4, "gated_allele_pairs"))
for ch in ["chr3", "chr8", "chr11", "chr12"]:
    s = g4[g4.chrom == ch]
    rows.append(table(ch, s.gated != "Major", s, "gated_allele_pairs"))
# (2c) gated minors mapped to stage-3 clusters (build_fig5.py rule), only where sizes agree
asg = pd.read_csv(args.assignments, sep="\t")
asg["sample_id"] = asg.sequence.str.split("#").str[0]
size = asg.cluster.value_counts()
o = asg[asg.cluster_type == "Outlier"]
orep = o.loc[o.groupby("cluster").raw_divergence.idxmax()]
gated = pd.read_csv(args.gated, sep="\t")
matched = []
for gg in gated.itertuples():
    reps = orep[(orep.chromosome == gg.chrom) & (orep.sample_id == gg.sample)]
    reps = reps[reps.cluster.map(size) == gg.n]
    if len(reps) == 1:
        matched.append(reps.cluster.iloc[0])
say(f"\n(2c) gated minor clusters whose stage-3 cluster (own representative) has the same n: {len(matched)} of {len(gated)}; "
    f"members {int(size[matched].sum())} of gated minor total {int(gated.n.sum())}")
mm = d.cluster.isin(matched)
dd = d[mm | (d.ungated == "Major")]  # matched gated minors vs ungated Major (gated Major is a superset; see report)
rows.append(table("ALL", dd.cluster.isin(matched), dd, "gated_minor_matched_vs_ungated_Major"))
# also: matched gated-minor members vs everything else in the 4,551
rows.append(table("ALL", mm, d, "gated_minor_matched_vs_rest"))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, "addendum_gap_vs_outlier.tsv"), sep="\t", index=False)
say("\n(2) 2x2 tables NNN>0 x Outlier (Fisher exact, two-sided)")
say(T.to_string(index=False))

say("\n(3) NNN bp under released Err, and assembly N runs (centromere regions, the 4,551)")
for t in ("hifi", "ont"):
    say(f"{t}: pre-release NNN {tot[(t,'NNN_bp')]:,} bp, under released Err {tot[(t,'NNN_bp_under_relErr')]:,} "
        f"({100*tot[(t,'NNN_bp_under_relErr')]/tot[(t,'NNN_bp')]:.2f}%); assembly N-run bp {tot[(t,'gap_bp')]:,}: "
        f"labelled NNN {100*tot[(t,'gap_bp_labelled_NNN')]/tot[(t,'gap_bp')]:.2f}%, under released Err "
        f"{100*tot[(t,'gap_bp_under_relErr')]/tot[(t,'gap_bp')]:.2f}%; N runs with no NNN {tot[(t,'gaps_without_NNN')]}, "
        f"with no released Err {tot[(t,'gaps_without_relErr')]} (of {int(d.fasta_gap_n.sum())})")
out.close()
