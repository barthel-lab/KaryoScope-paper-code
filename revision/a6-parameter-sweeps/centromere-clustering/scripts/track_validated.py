#!/usr/bin/env python3
"""Does each FISH-validated haplotype stay an Outlier across the (sil, sd) grid?

The counts sweep answers "how many outliers", which is not the published claim.
The published claim is that five specific minor subtypes are real, each backed by
FISH. This asks whether those five haplotypes are still called Outlier at every
setting in the neighbourhood of the published (0.5, 5).

Calls Yi-An's own apply_sil_and_centroid() from KS_allchr_barplot -- no
reimplementation of the gating.
"""
import os as _ks_os
KS_WORK = _ks_os.environ.get("KS_WORK", "$KS_WORK")  # see the repository README, Paths
import sys, types
sys.path.insert(0, (KS_WORK + "/a6_sweeps/cent/scripts"))
import pandas as pd
from KS_allchr_barplot import apply_sil_and_centroid

W = (KS_WORK + "/a6_sweeps/cent")
ASGN = f"{W}/out/allchr_structure.sequence_assignments.tsv"
PASS = f"{W}/out/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed"

# The five FISH-validated minor subtypes, from main.tex's validation paragraph,
# plus the chr9 trio parents (the inheritance test).
# TRAP: the manuscript names these by Coriell catalogue number (GM#####) but the
# assemblies use the 1000G/HPRC name (NA#####). GM21144 IS NA21144 and GM18960 IS
# NA18960 -- looking them up as written returns "absent", not an error. The table
# holds 3,476 HG and 1,170 NA samples and no GM at all.
VALIDATED = [
    ("chr3",  "NA21144", "HSat1A deletion (GM21144)"),
    ("chr5",  "HG00558", "HSat3 loss"),
    ("chr9",  "HG02630", "HSat3 megabase inversion"),
    ("chr12", "HG03816", "D12Z3 variant"),
    ("chr19", "NA18960", "D19Z3 variant (GM18960)"),
]

SILS = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7]
SDS = [3, 4, 5, 6, 7]

base = pd.read_csv(ASGN, sep="\t")
print(f"assignments: {len(base)} rows, columns: {list(base.columns)}")
samples = set(base["sequence"].str.split("#").str[0])
for chrom, samp, what in VALIDATED:
    present = samp in samples
    n = ((base["chromosome"] == chrom) & base["sequence"].str.startswith(samp + "#")).sum()
    print(f"  {chrom:>6} {samp:<9} {what:<28} in table: {present}, {n} haplotype(s) on {chrom}")

rows = []
for sil in SILS:
    for sd in SDS:
        args = types.SimpleNamespace(
            sil_threshold=sil, centroid_sd=float(sd), centroid_scan=True,
            bed=[PASS], exclude_features=None, edge_mode="directional",
            matrix_type="count_log1p_zscore")
        a = apply_sil_and_centroid(base.copy(), args)
        rec = dict(sil=sil, sd=sd)
        for chrom, samp, _ in VALIDATED:
            m = (a["chromosome"] == chrom) & a["sequence"].str.startswith(samp + "#")
            sub = a[m]
            if len(sub) == 0:
                rec[f"{chrom}:{samp}"] = "absent"
            else:
                types_ = set(sub["cluster_type"])
                rec[f"{chrom}:{samp}"] = "Outlier" if types_ - {"Major"} else "Major"
        rows.append(rec)
        print(f"  sil={sil} sd={sd} done", flush=True)

df = pd.DataFrame(rows)
dest = f"{W}/validated_haplotypes_by_setting.tsv"
df.to_csv(dest, sep="\t", index=False)
print(f"\nwrote {dest}\n")
cols = [c for c in df.columns if ":" in c]
print(f"{'sil':>5} {'sd':>3} | " + " | ".join(f"{c:>20}" for c in cols))
for _, r in df.iterrows():
    star = "  <- published" if (r["sil"] == 0.5 and r["sd"] == 5) else ""
    print(f"{r['sil']:>5} {r['sd']:>3} | " + " | ".join(f"{r[c]:>20}" for c in cols) + star)
print("\n=== survival summary ===")
for c in cols:
    vals = df[c].tolist()
    n_out = sum(v == "Outlier" for v in vals)
    print(f"  {c:<22} Outlier at {n_out}/{len(vals)} settings"
          + ("  <- INVARIANT" if n_out == len(vals) else ""))
