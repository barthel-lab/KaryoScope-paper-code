#!/usr/bin/env python3
"""Corrected D4Z4 sweep: salvage test at k = 21/51, plus the A6 re-test.

BUG FOUND 2026-09-03. Four functions bind module constants as DEFAULT ARGUMENTS,
evaluated at def time:
    get_bsr_units(bsat_regions, merge_threshold=BSR_MERGE_THRESHOLD)   # :81, called bare at :1763
    get_bsr_clusters(..., merge_threshold=BSR_MERGE_THRESHOLD)          # :990
    count_bsr_clusters(..., merge_threshold=BSR_MERGE_THRESHOLD)        # :1664
    detect_multiple_tracts(..., gap_threshold=TRACT_GAP_THRESHOLD)      # :999
Setting the module attribute therefore never reaches them, so BSR_MERGE_THRESHOLD
and TRACT_GAP_THRESHOLD were never actually varied -- in this sweep's first pass OR
in the A6 sweep, which reported both "completely inert". They are not. Fix: patch
f.__defaults__ as well. Floris's logic is still unmodified.

Part A -- salvage: can retuning recover the census at k = 21 / 51? Scored against
the HMM counts (external, k-independent), not against the published totals.
Part B -- A6 correction: the true sensitivity of the two affected constants on the
PUBLISHED inputs.
"""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import os, sys, shutil, importlib, logging
import numpy as np, pandas as pd
from pathlib import Path

W = Path((KS_WORK + "/a2_step4/job3_d4z4")); S = W / "salvage"
sys.path.insert(0, (KS_WORK + "/a6_sweeps/d4z4_bisect/code"))
HMM = pd.read_csv(W / "hmm/HPRC_chr4_counts.sheet1.tsv", sep="\t").rename(columns={"Count":"hmm_count","Contig":"core"})
HMM["hmm_count"] = HMM["hmm_count"].astype(float)
PUB = pd.read_csv(W / "base_published/results/d4z4_array_analysis.tsv", sep="\t").set_index("haplotype")

def run(datadir, tag, merge=None, minel=None, tract=None):
    base = S / "points2" / tag
    if base.exists(): shutil.rmtree(base)
    (base/"results").mkdir(parents=True); (base/"data").symlink_to(datadir)
    for m in [m for m in list(sys.modules) if m == "analyze_d4z4_arrays"]: del sys.modules[m]
    for h in logging.root.handlers[:]: logging.root.removeHandler(h)
    mod = importlib.import_module("analyze_d4z4_arrays")
    if merge is not None:
        mod.BSR_MERGE_THRESHOLD = merge
        for fn in ("get_bsr_units", "get_bsr_clusters", "count_bsr_clusters"):
            getattr(mod, fn).__defaults__ = (merge,)
    if tract is not None:
        mod.TRACT_GAP_THRESHOLD = tract
        mod.detect_multiple_tracts.__defaults__ = (tract,)
    if minel is not None:
        mod.MIN_ARRAY_ELEMENT_SIZE = minel
    os.environ["D4Z4_BASE_DIR"] = str(base); saved = sys.stdout
    try:
        with open(base/"stdout.log","w") as fh: sys.stdout = fh; mod.main()
    finally: sys.stdout = saved
    d = pd.read_csv(base/"results/d4z4_array_analysis.tsv", sep="\t")
    shutil.rmtree(base, ignore_errors=True)
    return d

def score(d):
    c4 = d[(d.chromosome=="chr4") & (d.num_arrays==1)].copy()
    c4["core"] = c4["haplotype"].str.replace(r"^chr\d+_hap\d+_","",regex=True)
    m = c4.merge(HMM[["core","hmm_count"]], on="core", how="inner")
    if len(m) >= 3 and m["repeat_count_main"].std() > 0:
        x = m["repeat_count_main"].astype(float)
        r = float(np.corrcoef(x, m["hmm_count"])[0,1]); md = float(np.median(np.abs(x-m["hmm_count"]))); bias = float((x-m["hmm_count"]).mean())
    else: r = md = bias = float("nan")
    dd = d.set_index("haplotype"); idx = dd.index.intersection(PUB.index)
    return dict(repeat_total=int(d.repeat_count_total.sum()), repeat_main=int(d.repeat_count_main.sum()),
                arrays=int(d.num_arrays.sum()), canonical=int(d.num_canonical.sum()), hmm_n=len(m), hmm_r=r,
                hmm_median_abs_diff=md, hmm_bias=bias,
                haps_eq_published=int((dd.loc[idx,"repeat_count_total"]==PUB.loc[idx,"repeat_count_total"]).sum()))

rows = []
print("=== PART A: salvage grid (v1-like annotations) ===", flush=True)
for kd in (21, 51, 31):
    merges = [200,300,400,600,800,1200,1600] if kd != 31 else [400]
    minels = [100,200,400,800,1200,1500,2000] if kd != 31 else [1500]
    for mg in merges:
        for mn in minels:
            d = run(W/f"base_v1like_k{kd}/data", f"A_k{kd}_m{mg}_s{mn}", merge=mg, minel=mn)
            s = score(d); s.update(part="salvage", k=kd, bsr_merge=mg, min_array_element=mn, tract_gap="")
            rows.append(s)
            print(f"k={kd} merge={mg:>5} minel={mn:>5} | total={s['repeat_total']:>6} main={s['repeat_main']:>6} "
                  f"arrays={s['arrays']:>4} | HMM r={s['hmm_r']:.4f} med|d|={s['hmm_median_abs_diff']:.1f} bias={s['hmm_bias']:+.2f}", flush=True)

print("\n=== PART B: A6 correction on the PUBLISHED inputs ===", flush=True)
PD = W/"base_published/data"
for mg in [200,300,400,600,800,1200]:
    d = run(PD, f"B_pub_m{mg}", merge=mg); s = score(d)
    s.update(part="a6_bsr_merge", k="published", bsr_merge=mg, min_array_element="", tract_gap=""); rows.append(s)
    print(f"published BSR_MERGE={mg:>5} | total={s['repeat_total']:>6} main={s['repeat_main']:>6} arrays={s['arrays']:>4} HMM r={s['hmm_r']:.4f}", flush=True)
for tg in [5000,7500,10000,15000,20000]:
    d = run(PD, f"B_pub_t{tg}", tract=tg); s = score(d)
    s.update(part="a6_tract_gap", k="published", bsr_merge="", min_array_element="", tract_gap=tg); rows.append(s)
    print(f"published TRACT_GAP={tg:>6} | total={s['repeat_total']:>6} main={s['repeat_main']:>6} arrays={s['arrays']:>4} HMM r={s['hmm_r']:.4f}", flush=True)

cols = ["part","k","bsr_merge","min_array_element","tract_gap","repeat_total","repeat_main","arrays","canonical",
        "hmm_n","hmm_r","hmm_median_abs_diff","hmm_bias","haps_eq_published"]
pd.DataFrame(rows)[cols].to_csv(S/"salvage_sweep_corrected.tsv", sep="\t", index=False)
print("\nwrote", S/"salvage_sweep_corrected.tsv")
