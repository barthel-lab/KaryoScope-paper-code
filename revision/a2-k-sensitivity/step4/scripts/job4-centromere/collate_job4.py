#!/usr/bin/env python3
"""Collate Job 4: centromere_by_k.tsv + centromere_by_k_summary.md."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import pandas as pd, numpy as np, os
W = (KS_WORK + "/a2_step4/job4_centromere"); OUT = (KS_WORK + "/a2_step4")
KS = ["v1_21", "v1_31", "v1_51", "published"]
V2 = ["v2_21", "v2_31", "v2_51"]
import os
files = [f"{W}/labels/labels_{k}.tsv" for k in KS + V2]
lab = pd.concat([pd.read_csv(f, sep="\t", dtype={"k": str}) for f in files if os.path.exists(f)], ignore_index=True)
KS = [k for k in KS if k in set(lab["k"])]
ALLK = KS + [k for k in V2 if k in set(lab["k"])]
lab.to_csv(f"{OUT}/centromere_by_k.tsv", sep="\t", index=False)
CH = [c for c in [f"chr{i}" for i in range(1, 23)] + ["chrX"] if c in set(lab["chromosome"])]
VALID = [("chr3", "NA21144", "HSat1A deletion (manuscript: GM21144)"), ("chr5", "HG00558", "HSat3 loss"),
         ("chr9", "HG02630", "HSat3 inversion"), ("chr12", "HG03816", "D12Z3 variant"),
         ("chr19", "NA18960", "D19Z3 variant (manuscript: GM18960); Major at every A6 setting")]
L = []
def p(s=""): L.append(s)
p("# Centromere subtype clustering at k = 21, 31, 51 -- summary\n")
p("Pipeline: KaryoScope-analysis@pangenome_structure_V2 (90e2874), unchanged, gate (0.5, 5) with its own "
  "--exclude-features default (novel) and its own EXCLUDED_CHROMS-before-gating order. "
  "Inputs: the published centromere windows re-annotated with the v1-like region rebuilds "
  "HKS_human_CHM13_v1like_region_k{21,31,51} (rows v1_*; the v2 six-k family is kept as rows v2_* for the "
  "database-generation contrast only -- its vocabulary differs from v1 so it is not comparable to published). "
  "(region set, "
  "presmoothed, 100 kb bins, trimmed to centromere-class extents as the published script did). "
  "'published' = the v1.1.0 inputs (A6 run, job 37136128) re-gated identically.\n")
p("## Per k: haplotypes, minor calls, chromosomes retaining a split\n")
p("| k | haplotypes | chromosomes | minor | major | chromosomes with a split | collapsed |"); p("|---|---|---|---|---|---|---|")
for k in ALLK:
    d = lab[lab["k"] == k]
    split = d.groupby("chromosome")["major_minor"].apply(lambda s: (s == "minor").any())
    coll = ",".join(c for c in CH if c in split.index and not split[c])
    p(f"| {k} | {len(d)} | {d['chromosome'].nunique()} | {(d['major_minor']=='minor').sum()} | {(d['major_minor']=='major').sum()} | {int(split.sum())} | {coll or '-'} |")
p("\n## Major / minor counts per chromosome\n")
p("| chromosome | " + " | ".join(f"k={k} major/minor" for k in KS) + " | k2 silhouette 21/31/51/pub |"); p("|---|" + "---|" * (len(KS) + 1))
for c in CH:
    cells = []
    sils = []
    for k in KS:
        d = lab[(lab["k"] == k) & (lab["chromosome"] == c)]
        cells.append(f"{(d['major_minor']=='major').sum()}/{(d['major_minor']=='minor').sum()}")
        sils.append(f"{d['chrom_silhouette_k2'].iloc[0]:.3f}" if len(d) and pd.notna(d['chrom_silhouette_k2'].iloc[0]) else "na")
    p(f"| {c} | " + " | ".join(cells) + " | " + "/".join(sils) + " |")
def concord(ka, kb):
    a = lab[lab["k"] == ka].set_index("haplotype"); b = lab[lab["k"] == kb].set_index("haplotype")
    common = a.index.intersection(b.index)
    same = (a.loc[common, "major_minor"] == b.loc[common, "major_minor"])
    per = same.groupby(a.loc[common, "chromosome"]).agg(["mean", "sum", "count"])
    return common, same, per, a, b
p("\n## Per-haplotype major/minor concordance vs k = 31 (new)\n")
for ka in ["v1_21", "v1_51", "published", "v2_31"]:
    common, same, per, a, b = concord(ka, "v1_31")
    only_a = len(set(a.index) - set(b.index)); only_b = len(set(b.index) - set(a.index))
    p(f"### k = {ka} vs k = 31(new): {same.sum()}/{len(common)} same call = {100*same.mean():.2f}%  "
      f"(haplotypes only in k={ka}: {only_a}; only in k=31: {only_b})\n")
    p("| chromosome | same | of | % |"); p("|---|---|---|---|")
    for c in CH:
        if c in per.index:
            p(f"| {c} | {int(per.loc[c,'sum'])} | {int(per.loc[c,'count'])} | {100*per.loc[c,'mean']:.1f} |")
    changed = common[~same.values]
    if len(changed):
        rows = [f"{h} ({a.loc[h,'chromosome']}: {a.loc[h,'major_minor']} at k={ka} -> {b.loc[h,'major_minor']} at k=31)" for h in changed]
        pd.DataFrame({"haplotype": changed, "chromosome": a.loc[changed, "chromosome"].values,
                      f"call_k{ka}": a.loc[changed, "major_minor"].values, "call_k31": b.loc[changed, "major_minor"].values}
                     ).to_csv(f"{OUT}/centromere_changes_{ka}_vs_v1_31.tsv", sep="\t", index=False)
        p(f"\nChanged calls ({len(changed)}; full list in centromere_changes_{ka}_vs_v1_31.tsv):\n")
        for r in rows[:40]: p(f"- {r}")
        if len(rows) > 40: p(f"- ... {len(rows)-40} more")
    p()
p("## k = 21 vs k = 51 directly (the brief's instability flag)\n")
common, same, per, a, b = concord("v1_21", "v1_51")
p(f"{same.sum()}/{len(common)} same call = {100*same.mean():.2f}%; {int((~same).sum())} haplotypes change major/minor between k=21 and k=51.\n")
p("## FISH-validated haplotypes\n")
p("| chromosome | sample | note | " + " | ".join(f"k={k}" for k in KS) + " |"); p("|---|---|---|" + "---|" * len(KS))
for c, s, note in VALID:
    cells = []
    for k in KS:
        d = lab[(lab["k"] == k) & (lab["chromosome"] == c) & lab["haplotype"].str.startswith(s + "#")]
        cells.append("absent" if len(d) == 0 else ",".join(f"{h.split('#')[1]}:{m}" for h, m in zip(d["haplotype"], d["major_minor"])))
    p(f"| {c} | {s} | {note} | " + " | ".join(cells) + " |")
p("\n(hap index:call; both haplotypes of the sample on that chromosome are listed when both passed QC)")
open(f"{OUT}/centromere_by_k_summary.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
