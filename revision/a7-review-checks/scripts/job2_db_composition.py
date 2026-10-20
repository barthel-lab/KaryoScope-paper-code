#!/usr/bin/env python3
"""CHM13 v2 gene set: composition of the build INPUT by reference position (recipe's chm13_gene.bed,
sha256 03ae30a1..., 612,909 rows) and of the built index by k=31 node count (A4 `hks node-stats`)."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import collections
BED = (KS_WORK + "/recipe_repro/chm13_gene.bed")
NS = (KS_WORK + "/a4_smoothing/raw_nodestats/nodestats_shipped_gene.txt")
OUT = (KS_WORK + "/a7_review_checks/results/chm13_gene_db_composition.tsv")
bp_all = collections.Counter(); bp_nuc = collections.Counter()
for line in open(BED):
    c, s, e, lab = line.rstrip("\n").split("\t")[:4]
    L = int(e) - int(s); bp_all[lab] += L
    if c != "chrM": bp_nuc[lab] += L
kmer = {}
for line in open(NS):
    p = line.rstrip("\n").split("\t")
    if p[0] == "31": kmer[p[1]] = int(p[2])
rows = []
for measure, d, note in (("input_bed_bp_all_sequences", bp_all, "T2T-CHM13v2.0 bases per label in the build input BED (before k-mer resolution), chrM included"),
                         ("input_bed_bp_excluding_chrM", bp_nuc, "same, chrM dropped (the build excludes chrM)"),
                         ("index_nodes_k31", kmer, "hks node-stats count at k=31 per label after priority resolution (exon>intron>intergenic); index stores both strands (hks add-feature-set default, no --forward-only), so counts are ~2x distinct canonical k-mers")):
    tot = sum(d.values())
    for lab in ("exon", "intron", "intergenic", "categorized"):
        rows.append((measure, lab, d.get(lab, 0), 100*d.get(lab, 0)/tot, tot, note))
with open(OUT, "w") as o:
    o.write("measure\tlabel\tcount\tpct\ttotal\tnote\n")
    for r in rows: o.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
for r in rows: print(*r[:5], sep="\t")
