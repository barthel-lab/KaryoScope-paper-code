#!/usr/bin/env python3
"""Collate Job 2 -> divergence_by_k.tsv (brief columns first) + control check vs A5."""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import os, glob, math, pandas as pd
W = (KS_WORK + "/a2_step4/job2_divergence"); A5 = (KS_WORK + "/a5_divergence")
OUT = (KS_WORK + "/a2_step4/divergence_by_k.tsv")
conds = ["sub_sub0.001_ind0_seed42", "sub_sub0.01_ind0_seed42", "sub_sub0.05_ind0_seed42", "indel_sub0_ind0.005_seed42"]
rows = []
for tag in conds:
    f = f"{W}/{tag}/results/divergence_rows.tsv"
    if not os.path.exists(f): print("MISSING", f); continue
    d = pd.read_csv(f, sep="\t")
    a5 = pd.read_csv(f"{A5}/{tag}/results/divergence_rows.tsv", sep="\t")
    for (k, fs), g in d.groupby(["k", "feature_set"]):
        def frac(stage, cls, src=g):
            m = src[(src["stage"] == stage) & (src["class"] == cls)]
            return float(m["fraction_pct"].iloc[0]) / 100 if len(m) else float("nan")
        sub, ind = float(g["sub_rate"].iloc[0]), float(g["indel_rate"].iloc[0])
        mtype = "substitution" if sub > 0 else "indel"; rate = sub if sub > 0 else ind
        pred = 1 - (1 - rate) ** k
        a5g = a5[a5["feature_set"] == fs]
        rows.append(dict(k=k, mutation_type=mtype, rate=rate, seed=42, feature_set=fs,
            novel_frac_measured=frac("presmoothed", "novel"), novel_frac_predicted=pred,
            exact_pre=frac("presmoothed", "equivalent"), exact_post=frac("smoothed", "equivalent"),
            compatible_pre=frac("presmoothed", "compatible"), compatible_post=frac("smoothed", "compatible"),
            incompatible_pre=frac("presmoothed", "incompatible"), incompatible_post=frac("smoothed", "incompatible"),
            novel_post=frac("smoothed", "novel"),
            a5_k31_shipped_db_novel_pre=frac("presmoothed", "novel", a5g), a5_k31_shipped_db_exact_pre=frac("presmoothed", "equivalent", a5g),
            a5_k31_shipped_db_exact_post=frac("smoothed", "equivalent", a5g)))
df = pd.DataFrame(rows).sort_values(["feature_set", "mutation_type", "rate", "k"], ascending=[False, False, True, True])
df.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
print(f"wrote {OUT} ({len(df)} rows)\n")
pd.set_option("display.width", 250)
print(df[df.feature_set == "repeat"][["k", "mutation_type", "rate", "novel_frac_measured", "novel_frac_predicted", "exact_pre", "exact_post", "incompatible_post"]].to_string(index=False))
c = df[(df.k == 31)]
print("\nCONTROL k=31 six-k DB vs A5 shipped DB (must be identical):")
for _, r in c.iterrows():
    same = abs(r.novel_frac_measured - r.a5_k31_shipped_db_novel_pre) < 5e-7 and abs(r.exact_pre - r.a5_k31_shipped_db_exact_pre) < 5e-7 and abs(r.exact_post - r.a5_k31_shipped_db_exact_post) < 5e-7
    print(f"  {r.feature_set:7} {r.mutation_type:12} r={r.rate}: novel {r.novel_frac_measured:.6f} vs {r.a5_k31_shipped_db_novel_pre:.6f}; exact_pre {r.exact_pre:.6f} vs {r.a5_k31_shipped_db_exact_pre:.6f}; exact_post {r.exact_post:.6f} vs {r.a5_k31_shipped_db_exact_post:.6f}  -> {'IDENTICAL' if same else 'DIFFERS'}")
print("\nRegion set:")
print(df[df.feature_set == "region"][["k", "mutation_type", "rate", "novel_frac_measured", "exact_pre", "exact_post", "incompatible_post"]].to_string(index=False))
