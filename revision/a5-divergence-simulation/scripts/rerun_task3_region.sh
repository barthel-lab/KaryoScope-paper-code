#!/bin/bash
#SBATCH --job-name=a5_t3_region
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=a5_t3_region_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
D=${KS_WORK}/a5_divergence/sub_sub0.01_ind0_seed42
RES=$D/results
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
bash ${KS_WORK}/ks_concord.sh \
  --truth "$D/mut.region.projected.bed" \
  --pre  "$(ls $D/annot/*.region.presmoothed.bed | head -1)" \
  --post "$(ls $D/annot/*.region.smoothed.bed | head -1)" \
  --hierarchy $ROOT/HKS_human_CHM13_v2/hierarchy.tsv \
  --feature-set region --fai "$D/mut.fa.fai" --gaps "$D/mut.n_regions.bed" \
  --kmer-size 31 --out "$RES/region"
echo "rc=$?"
# rebuild the tidy rows now that both sets exist
python3 - "$RES" sub 0.01 0 42 <<'PY'
import sys, os
RES, series, sub, ind, seed = sys.argv[1:6]
with open(os.path.join(RES,"divergence_rows.tsv"),"w") as o:
    o.write("series\tsub_rate\tindel_rate\tseed\tfeature_set\tstage\tclass\tbases\tfraction_pct\n")
    for fs in ("repeat","region"):
        p=os.path.join(RES,f"{fs}.summary.tsv")
        if not os.path.exists(p): continue
        with open(p) as fh:
            next(fh)
            for line in fh:
                stage,c,b,frac=line.rstrip("\n").split("\t")
                o.write(f"{series}\t{sub}\t{ind}\t{seed}\t{fs}\t{stage}\t{c}\t{b}\t{frac}\n")
print(open(os.path.join(RES,"divergence_rows.tsv")).read())
PY
if [ -s "$RES/repeat.summary.tsv" ] && [ -s "$RES/region.summary.tsv" ]; then
  echo "both sets present -> freeing intermediates"
  rm -f "$D"/annot/*.presmoothed.bed "$D"/annot/*.smoothed.bed "$D/mut.fa"
fi
