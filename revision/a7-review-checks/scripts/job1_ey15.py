#!/usr/bin/env python3
"""Job 1 -- Ey15-2 vs the Col-CEN database: composition by stage (1a) and per-sequence
chromosome matrix (1b). Reads the A3 BEDs as left by job 35378568 (both stages were kept)."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import sys, os, collections
sys.path.insert(0, os.path.dirname(__file__))
from common_mask import *
A3 = (KS_WORK + "/a3_crossgenome")
RES = (KS_WORK + "/a7_review_checks/results")
HIER = f"{A3}/db/HKS_arabidopsis_ColCEN/hierarchy.tsv"
SETS = ["chromosome", "region", "repeat", "gene"]
RUNS = {  # assembly: (annot dir, file prefix, fasta)
    "ey15_2": (f"{A3}/annot_ey15_2", "ey15-2.HKS_arabidopsis_ColCEN", f"{A3}/fa/ey15-2.fa"),
    "colcen_nuclear": (f"{A3}/annot_colcen_nuclear", "colcen_nuclear.HKS_arabidopsis_ColCEN", f"{A3}/fa/colcen_nuclear.fa"),
}
# published smoothed novel bp (a3_crossgenome/results/label_composition.tsv, 2026-08-24)
PUBLISHED_SMOOTHED_NOVEL = {"chromosome": 105, "region": 760084, "repeat": 2497066, "gene": 4824275}
H = parse_hierarchy(HIER)

def read_headers(fa):
    """sequence -> chromosome named in the assembly record's FASTA header (ENA/NCBI)."""
    m = {}
    for line in open(fa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            chrom = None
            if "chromosome:" in line: chrom = "Chr" + line.split("chromosome:")[1].strip().split()[0]
            elif name.startswith("Chr"): chrom = name
            m[name] = chrom
    return m

rows_comp, rows_lab, rows_seq, rows_san, prov = [], [], [], [], []
for asm, (d, pre, fa) in RUNS.items():
    mask = n_runs_from_fasta(fa)
    ena = read_headers(fa)
    prov.append(f"{asm}: N-mask (N runs extended 30 bp left) = " + (", ".join(f"{c}:{a}-{b}" for c, iv in mask.items() for a, b in iv) or "none"))
    prov.append(f"{asm}: sequence -> chromosome from FASTA headers: " + ", ".join(f"{k}={v}" for k, v in ena.items()))
    for fs in SETS:
        nodes, leaves, parent_of = H[fs]
        lvl = make_level(nodes, leaves)
        for stage in ("presmoothed", "smoothed"):
            f = f"{d}/{pre}.{fs}.{stage}.bed.gz"
            assert stage in os.path.basename(f) and (stage != "smoothed" or "presmoothed" not in os.path.basename(f))
            by_lvl = collections.Counter(); by_lab = collections.Counter(); by_seq = collections.defaultdict(collections.Counter)
            unmasked_novel = 0; unmasked_total = 0
            with opener(f) as fh:
                for line in fh:
                    c, s, e, lab = line.rstrip("\n").split("\t")[:4]
                    s, e = int(s), int(e); L = e - s
                    unmasked_total += L
                    if lab.lower() == "novel": unmasked_novel += L
                    eff = L - overlap(mask.get(c, ()), s, e)
                    if eff <= 0: continue
                    lv = lvl(lab); assert lv != "other", (f, lab)
                    by_lvl[lv] += eff; by_lab[lab] += eff; by_seq[c][lab] += eff
            den = sum(by_lvl.values())
            Ls, M, V = by_lvl["leaf"], by_lvl["internal"], by_lvl["novel"]
            rows_comp.append((asm, fs, stage, Ls, M, V, 100*Ls/den, 100*M/den, 100*V/den, den, 100*(Ls+M+V)/den))
            for lab, bp in sorted(by_lab.items(), key=lambda x: -x[1]):
                rows_lab.append((asm, fs, stage, lab, lvl(lab), bp, 100*bp/den))
            if stage == "smoothed" and asm == "ey15_2":
                exp = PUBLISHED_SMOOTHED_NOVEL[fs]
                rows_san.append((asm, fs, unmasked_total, unmasked_novel, exp, "MATCH" if unmasked_novel == exp else "MISMATCH", V))
            if fs == "chromosome":
                labels = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "nuclear", "novel"]
                for c in sorted(by_seq):  # fai order == accession order here
                    tot = sum(by_seq[c].values()); own = ena.get(c)
                    own_bp = by_seq[c].get(own, 0) if own else 0
                    other_chr = sum(v for k, v in by_seq[c].items() if k.startswith("Chr") and k != own)
                    rows_seq.append((asm, stage, c, own or "NA", *[by_seq[c].get(l, 0) for l in labels], tot,
                                     100*own_bp/tot, 100*other_chr/tot, 100*by_seq[c].get("nuclear", 0)/tot))
    prov.append(f"{asm}: hierarchy leaves per set: " + "; ".join(f"{fs}={sorted(H[fs][1])}" for fs in SETS))
    prov.append(f"{asm}: hierarchy internal nodes per set: " + "; ".join(f"{fs}={sorted(H[fs][0]-H[fs][1])}" for fs in SETS))

with open(f"{RES}/ey15_2_novel_by_stage.tsv", "w") as o:
    o.write("assembly\tfeature_set\tstage\tlabel_specific_bp\tmultilabel_bp\tnovel_bp\tlabel_specific_pct\tmultilabel_pct\tnovel_pct\tdenominator_bp\tsum_pct\n")
    for r in rows_comp: o.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
with open(f"{RES}/ey15_2_labels_by_stage.tsv", "w") as o:
    o.write("assembly\tfeature_set\tstage\tlabel\tlevel\tbp\tpct\n")
    for r in rows_lab: o.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
with open(f"{RES}/ey15_2_chromosome_by_sequence.tsv", "w") as o:
    o.write("assembly\tstage\tsequence\tchromosome_per_assembly_record\tChr1\tChr2\tChr3\tChr4\tChr5\tnuclear\tnovel\ttotal_bp\tpct_own_chromosome\tpct_other_chromosome\tpct_nuclear\n")
    for r in rows_seq: o.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")
with open(f"{RES}/ey15_2_novel_sanity.tsv", "w") as o:
    o.write("assembly\tfeature_set\tsmoothed_tiled_bp_unmasked\tsmoothed_novel_bp_unmasked\tpublished_label_composition_novel_bp\tcheck\tsmoothed_novel_bp_N_masked\n")
    for r in rows_san: o.write("\t".join(str(x) for x in r) + "\n")
with open(f"{RES}/job1_provenance.txt", "w") as o:
    o.write("Job 1 inputs: A3 annotate outputs (karyoscope 2.3.0, hks 0.4.0, HKS_arabidopsis_ColCEN 1.2.0, job 35378568, 2026-08-11); both stages were kept, no re-run.\n")
    o.write("\n".join(prov) + "\n")
for r in rows_comp: print(*r, sep="\t")
print("--- sanity"); [print(*r, sep="\t") for r in rows_san]
print("--- per sequence (chromosome set)"); [print(*r[:4], *[f"{x:.2f}" if isinstance(x, float) else x for x in r[4:]], sep="\t") for r in rows_seq]
