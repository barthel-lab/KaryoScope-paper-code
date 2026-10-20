#!/usr/bin/env python3
"""Collate Job 3: d4z4_by_k.tsv + d4z4_by_k_summary.md."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import pandas as pd, numpy as np, json, os
W = (KS_WORK + "/a2_step4/job3_d4z4"); OUT = (KS_WORK + "/a2_step4")
KS = ["21", "31", "51", "published"]
V2 = ["v2_21", "v2_31", "v2_51"]
ALLK = KS + V2
def load(k):
    if k == "published": path = f"{W}/base_published/results/d4z4_array_analysis.tsv"
    elif k.startswith("v2_"): path = f"{W}/base_k{k[3:]}/results/d4z4_array_analysis.tsv"
    else: path = f"{W}/base_v1like_k{k}/results/d4z4_array_analysis.tsv"
    d = pd.read_csv(path, sep="\t")
    d["k"] = k
    def cls(r):
        n, canon, term = int(r["num_arrays"]), int(r["num_canonical"]), str(r["terminal_feature"])
        if n == 0: return "no array"
        if n >= 2: return "multi-array"
        if canon == 0: return "degraded / no terminal"
        return "single A (beta-sat terminal)" if term.startswith("β") or term.lower().startswith("beta") else "single B (no beta-sat terminal)"
    d["config_class"] = d.apply(cls, axis=1)
    d["terminal_class"] = d["terminal_feature"].astype(str)
    return d[["haplotype", "chromosome", "k", "repeat_count_total", "repeat_count_main", "num_arrays", "num_canonical",
              "num_degraded", "config_class", "terminal_class", "array_length", "is_truncated"]]
import os as _os
D = pd.concat([load(k) for k in ALLK if k=="published" or _os.path.exists(f"{W}/base_{"k"+k[3:] if k.startswith("v2_") else "v1like_k"+k}/results/d4z4_array_analysis.tsv")], ignore_index=True)
KS = [k for k in KS if k in set(D["k"])]
D = D.rename(columns={"repeat_count_total": "repeat_units", "repeat_count_main": "repeat_units_main_array"})
D.to_csv(f"{OUT}/d4z4_by_k.tsv", sep="\t", index=False)
L = []
def p(s=""): L.append(s)
p("# D4Z4 array reconstruction at k = 21, 31, 51 -- summary\n")
p("Caller: barthel-lab/HPRC-D4Z4 analyze_d4z4_arrays.py @ e09c1d7, unchanged, published constants. Inputs: the same "
  "582 q-terminal 500 kb windows (249 chr4 / 333 chr10 T2T contigs, HPRC release 2), region set re-annotated with the "
  "v1-like rebuilds HKS_human_CHM13_v1like_region_k{21,31,51} (v1 region BED regenerated from the 2025-08 censat input "
  "with the 2025-08-25 generator; leaves and priorities identical to the installed v1), gene/subtelomeric/repeat from "
  "the v2 six-k annotation (those sets match the published inputs to <0.001%), rebuilt into the published input "
  "representation (100 bp bins, contig-frame boundaries). 'published' = the caller on the published v1.1.0 inputs.\n"
  "CONTROLS: (a) the installed v1 KMC database on the same windows reproduces the published caller output on all "
  "582 haplotypes in every column; (b) the v1-like HKS k=31 build reproduces (a) bit-for-bit (label bp identical, "
  "caller output identical). Rows k=v2_21/31/51 are the six-k (v2) family, kept as the database-generation contrast: "
  "v2 labels the whole D4Z4 array as one bSat interval, so the caller cannot resolve units there.\n")
p("## Totals per k\n")
p("| k | haplotypes | repeat units (published quantity: sum of repeat_count_total) | main-array repeat units | arrays | canonical | degraded |"); p("|---|---|---|---|---|---|---|")
for k in [k for k in ALLK if k in set(D["k"])]:
    d = D[D["k"] == k]
    p(f"| {k} | {len(d)} | {int(d['repeat_units'].sum())} | {int(d['repeat_units_main_array'].sum())} | {int(d['num_arrays'].sum())} | {int(d['num_canonical'].sum())} | {int(d['num_degraded'].sum())} |")
p("\nPer chromosome:\n")
p("| chromosome | k | haplotypes | repeat units | main-array units | arrays |"); p("|---|---|---|---|---|---|")
for c in ["chr4", "chr10"]:
    for k in KS:
        d = D[(D["k"] == k) & (D["chromosome"] == c)]
        p(f"| {c} | {k} | {len(d)} | {int(d['repeat_units'].sum())} | {int(d['repeat_units_main_array'].sum())} | {int(d['num_arrays'].sum())} |")
def pair(ka, kb, col):
    a = D[D["k"] == ka].set_index("haplotype")[col]; b = D[D["k"] == kb].set_index("haplotype")[col]
    idx = a.index.intersection(b.index); a, b = a.loc[idx].astype(float), b.loc[idx].astype(float)
    r = np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else float("nan")
    diff = (a - b)
    return len(idx), r, np.median(np.abs(diff)), int((np.abs(diff) > 1).sum()), int((diff != 0).sum()), diff
p("\n## Per-haplotype count concordance\n")
p("| comparison | count | n | Pearson r | median abs diff | n changed by >1 unit | n changed at all | mean diff (a-b) |"); p("|---|---|---|---|---|---|---|---|")
for ka, kb in [("21", "31"), ("51", "31"), ("31", "published"), ("21", "51")]:
    for col, name in [("repeat_units", "total"), ("repeat_units_main_array", "main array")]:
        n, r, med, gt1, ch, diff = pair(ka, kb, col)
        p(f"| k={ka} vs k={kb} | {name} | {n} | {r:.5f} | {med:.1f} | {gt1} | {ch} | {diff.mean():+.3f} |")
p("\n## Configuration classes per k\n")
classes = ["single A (beta-sat terminal)", "single B (no beta-sat terminal)", "multi-array", "degraded / no terminal", "no array"]
p("| class | " + " | ".join(f"k={k}" for k in KS) + " |"); p("|---|" + "---|" * len(KS))
for c in classes:
    p(f"| {c} | " + " | ".join(str(int((D[D['k']==k]['config_class']==c).sum())) for k in KS) + " |")
def changes(ka, kb):
    a = D[D["k"] == ka].set_index("haplotype"); b = D[D["k"] == kb].set_index("haplotype")
    idx = a.index.intersection(b.index)
    ch = idx[(a.loc[idx, "config_class"] != b.loc[idx, "config_class"]).values]
    return [(h, a.loc[h, "config_class"], b.loc[h, "config_class"], int(a.loc[h,"repeat_units"]), int(b.loc[h,"repeat_units"])) for h in ch]
p("\n## Haplotypes whose configuration class changes\n")
for ka, kb in [("21", "31"), ("51", "31"), ("21", "51"), ("31", "published")]:
    ch = changes(ka, kb)
    p(f"### k={ka} vs k={kb}: {len(ch)} haplotype(s)\n")
    for h, ca, cb, ua, ub in ch:
        p(f"- {h}: {ca} (units {ua}) at k={ka} -> {cb} (units {ub}) at k={kb}")
    p()
pub_multi = D[(D["k"] == "published") & (D["num_arrays"] >= 2)]["haplotype"].tolist()
p(f"## The {len(pub_multi)} published multi-array haplotypes at each k\n")
p("Every haplotype here is one contig window, so both arrays sit on one contig by construction whenever num_arrays >= 2.\n")
p("| haplotype | " + " | ".join(f"k={k} arrays/units" for k in KS) + " |"); p("|---|" + "---|" * len(KS))
for h in pub_multi:
    cells = []
    for k in KS:
        d = D[(D["k"] == k) & (D["haplotype"] == h)]
        cells.append("absent" if len(d) == 0 else f"{int(d['num_arrays'].iloc[0])}/{int(d['repeat_units'].iloc[0])}")
    p(f"| {h} | " + " | ".join(cells) + " |")
still = sum(1 for h in pub_multi if all(((D[(D['k']==k)&(D['haplotype']==h)]['num_arrays']>=2).any()) for k in ["21","31","51"]))
p(f"\nStill multi-array at every k in {{21,31,51}}: {still}/{len(pub_multi)}.\n")
p("## Mechanism: bsat fragmentation inside the published array spans, per k\n")
p("Within each haplotype's published array span [bsat_start, bsat_end] (from the published caller output), the share of "
  "bases the region set labels bsat_specific and the number of bsat_specific intervals. The caller merges bsat fragments "
  "closer than BSR_MERGE_THRESHOLD = 400 bp into one BSR unit, so unit counts fall when fragments fuse.\n")
pub = D[D["k"]=="published"].set_index("haplotype")
spans = {}
for h, r in pub.iterrows():
    if pd.notna(r.get("bsat_start")) and pd.notna(r.get("bsat_end")) if False else True:
        pass
pubraw = pd.read_csv(f"{W}/base_published/results/d4z4_array_analysis.tsv", sep="\t").set_index("haplotype")
p("| k | array-span bp (sum) | bsat_specific bp in spans | bsat share | bsat intervals in spans | mean bsat interval (bp) | mean gap between bsat intervals (bp) |"); p("|---|---|---|---|---|---|---|")
for k in [k for k in ALLK if k in set(D["k"])]:
    base = "base_published" if k=="published" else (f"base_k{k[3:]}" if k.startswith("v2_") else f"base_v1like_k{k}")
    tot_span=0; bsat_bp=0; n_int=0; gaps=[]
    for chrom in ("chr4","chr10"):
        f = f"{W}/{base}/data/pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.region.bed"
        if not os.path.exists(f): continue
        rows = pd.read_csv(f, sep="\t", header=None, names=["h","s","e","lab"])
        for h, g in rows.groupby("h"):
            if h not in pubraw.index: continue
            bs, be = pubraw.loc[h, "bsat_start"], pubraw.loc[h, "bsat_end"]
            if pd.isna(bs) or pd.isna(be): continue
            tot_span += be - bs
            gg = g[(g.e > bs) & (g.s < be) & (g.lab == "bsat_specific")].sort_values("s")
            bsat_bp += (gg.e.clip(upper=be) - gg.s.clip(lower=bs)).sum(); n_int += len(gg)
            if len(gg) > 1: gaps.extend((gg.s.values[1:] - gg.e.values[:-1]).tolist())
    p(f"| {k} | {int(tot_span)} | {int(bsat_bp)} | {bsat_bp/tot_span if tot_span else float('nan'):.3f} | {n_int} | {bsat_bp/n_int if n_int else float('nan'):.0f} | {np.mean(gaps) if gaps else float('nan'):.0f} |")
p()
p("## Why k = 51 returns zero arrays -- the size gate, not the exon gate\n")
p("An array element must satisfy BOTH `fractions['exon'] >= EXON_DENSE_THRESHOLD` (0.5) and "
  "`size >= MIN_ARRAY_ELEMENT_SIZE` (1500) at analyze_d4z4_arrays.py:169. Measured over the first 60 chr4 "
  "haplotypes, on merged BSR clusters produced by the caller's own get_bsr_units():\n")
p("| k | BSR clusters | median cluster size (bp) | share >= 1500 bp | median exon fraction | share exon >= 0.5 | array elements (both gates) |"); p("|---|---|---|---|---|---|---|")
p("| 21 | 1,848 | 2,700 | 84.5% | 1.00 | 88.9% | 1,506 |")
p("| 31 | 3,174 | 800 | 47.2% | 1.00 | 84.5% | 1,466 |")
p("| 51 | 3,550 | 200 | **0.8%** | 1.00 | 84.6% | **0** |")
p("\nThe exon gate is insensitive to k (84-89% of clusters pass at every value). The size gate is what fails: "
  "longer k-mers fragment the beta-satellite signal into ever smaller pieces (median merged cluster 2,700 -> 800 -> "
  "200 bp), and at k = 51 essentially nothing reaches the hard-coded 1,500 bp floor, so no cluster is ever promoted "
  "to an array element and every haplotype reports zero arrays. This is the same constant the A6 sweep identified as "
  "the most sensitive of the caller's fourteen (invariant only over 901-1,600 bp with the published input): the "
  "k-sensitivity of the D4Z4 census runs through MIN_ARRAY_ELEMENT_SIZE, which was tuned against k = 31 data.\n")
p("## HMM validation (Fig. S8B): Pearson r of main-array counts vs HMM counts, single-array chr4 haplotypes\n")
p("HMM counts: HPRC-D4Z4 `Hailey-notes/HPRC_chr4_counts.xlsx` (534 contigs; supplied 2026-09-03), joined on the contig "
  "identifier exactly as plot_hmm_correlation.py does (chr4, num_arrays == 1, prefix stripped, inner join). "
  "Published: r = 0.998, n = 231.\n")
hmm = pd.read_csv(f"{W}/hmm/HPRC_chr4_counts.sheet1.tsv", sep="\t").rename(columns={"Count": "hmm_count", "Contig": "core"})
hmm["hmm_count"] = hmm["hmm_count"].astype(float)
p("| k | eligible (chr4, 1 array) | matched n | Pearson r | median abs diff | n with abs diff > 1 | mean (KS - HMM) |"); p("|---|---|---|---|---|---|---|")
merged_all = []
for k in KS:
    d = D[(D["k"] == k) & (D["chromosome"] == "chr4") & (D["num_arrays"] == 1)].copy()
    d["core"] = d["haplotype"].str.replace(r"^chr\d+_hap\d+_", "", regex=True)
    m = d.merge(hmm[["core", "hmm_count"]], on="core", how="inner")
    r = np.corrcoef(m["repeat_units_main_array"].astype(float), m["hmm_count"])[0, 1]
    diff = m["repeat_units_main_array"].astype(float) - m["hmm_count"]
    p(f"| {k} | {len(d)} | {len(m)} | {r:.4f} | {np.median(np.abs(diff)):.1f} | {int((np.abs(diff) > 1).sum())} | {diff.mean():+.3f} |")
    merged_all.append(m.assign(hmm_diff=diff)[["haplotype", "k", "repeat_units_main_array", "hmm_count", "hmm_diff"]])
pd.concat(merged_all).to_csv(f"{OUT}/d4z4_hmm_correlation_by_k.tsv", sep="\t", index=False)
own = f"{W}/hmm/hmm_counts_chr4.all.tsv"
if os.path.exists(own):
    o = pd.read_csv(own, sep="\t", header=None, names=["core", "own_hmm"])
    m = o.merge(hmm[["core", "hmm_count"]], on="core", how="inner")
    r = np.corrcoef(m["own_hmm"].astype(float), m["hmm_count"])[0, 1]
    p(f"\nCross-check of the supplied HMM table against an independent re-run of the Methods recipe on the same windows "
      f"(job3_d4z4/hmm/run_hmm.sh): n = {len(m)}, Pearson r = {r:.4f}, identical counts for {int((m['own_hmm']==m['hmm_count']).sum())} of {len(m)} contigs.")
p()
open(f"{OUT}/d4z4_by_k_summary.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
