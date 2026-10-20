#!/usr/bin/env python3
"""Per-haplotype Major/Outlier labels at (sil, sd) from a step-3 assignments table,
using Yi-An's own apply_sil_and_centroid() (KS_allchr_barplot.py @ 90e2874), plus
the per-chromosome k=2 Ward silhouette that the gate tests (recomputed with the
same calls the gate makes). No reimplementation of the gating itself.
"""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import argparse, sys, types
sys.path.insert(0, (KS_WORK + "/a6_sweeps/cent/scripts"))
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score
from KS_allchr_barplot import apply_sil_and_centroid
from KS_allchr_dendrogram import load_bed_features, build_global_matrix, CHROM_ORDER, EXCLUDED_CHROMS

ap = argparse.ArgumentParser()
ap.add_argument("--assignments", required=True)
ap.add_argument("--bed", required=True, help="the QC pass.bed")
ap.add_argument("--k", required=True, help="label for the k column (21/31/51/published)")
ap.add_argument("--sil", type=float, default=0.5)
ap.add_argument("--sd", type=float, default=5.0)
ap.add_argument("--exclude-features", default="novel", help="must match KS_allchr_barplot.py's default")
ap.add_argument("--out", required=True)
a = ap.parse_args()

base = pd.read_csv(a.assignments, sep="\t")
# KS_allchr_barplot.py filters EXCLUDED_CHROMS BEFORE gating (line 154 then 158); match that.
base = base[~base["chromosome"].isin(EXCLUDED_CHROMS)].copy()
args = types.SimpleNamespace(sil_threshold=a.sil, centroid_sd=a.sd, centroid_scan=True,
                             bed=[a.bed], exclude_features=a.exclude_features, edge_mode="directional",
                             matrix_type="count_log1p_zscore")
gated = apply_sil_and_centroid(base.copy(), args)

# per-chromosome k=2 silhouette, exactly as the gate computes it
feat, _, _ = load_bed_features([a.bed], set(base["sequence"]), a.exclude_features)
sil_by_chrom = {}
for chrom in CHROM_ORDER:
    cdf = base[base["chromosome"] == chrom]
    if cdf["cluster"].nunique() <= 1:
        sil_by_chrom[chrom] = float("nan"); continue
    seqs = [s for s in sorted(cdf["sequence"]) if s in feat]
    if len(seqs) < 5:
        sil_by_chrom[chrom] = float("nan"); continue
    feats = sorted(set(f for s in seqs for f, _ in feat.get(s, [])))
    M, _ = build_global_matrix(seqs, feat, feats, "directional", "count_log1p_zscore")
    lab = fcluster(linkage(pdist(M, metric="euclidean"), method="ward"), 2, criterion="maxclust")
    sil_by_chrom[chrom] = silhouette_score(M, lab) if len(set(lab)) >= 2 else 0.0

gated = gated.copy()
gated["k"] = a.k
gated["major_minor"] = gated["cluster_type"].map(lambda t: "major" if t == "Major" else "minor")
gated["chrom_silhouette_k2"] = gated["chromosome"].map(sil_by_chrom)
gated["step3_cluster"] = base.set_index("sequence").loc[gated["sequence"], "cluster"].values
out = gated[["sequence", "chromosome", "k", "cluster", "major_minor", "chrom_silhouette_k2", "norm_divergence", "step3_cluster"]]
out = out.rename(columns={"sequence": "haplotype", "cluster": "cluster_label"})
out.to_csv(a.out, sep="\t", index=False)
n_split = out.groupby("chromosome")["major_minor"].apply(lambda s: (s == "minor").any()).sum()
print(f"k={a.k}: {len(out)} haplotypes, {out['chromosome'].nunique()} chromosomes, "
      f"{(out['major_minor']=='minor').sum()} minor, {n_split} chromosomes retain a split")
