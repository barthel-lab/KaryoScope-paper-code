#!/bin/bash
#SBATCH --job-name=a5_reconcord
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=a5_reconcord_%A_%a.log
#SBATCH --array=1-2
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
MERGE="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"
case $SLURM_ARRAY_TASK_ID in
  1) TAG=baseline_sub0_ind0_seed42;   SET=repeat;;
  2) TAG=sub_sub0.02_ind0_seed42;     SET=region;;
esac
D=${KS_WORK}/a5_divergence/$TAG; RES=$D/results
echo "re-running $SET concordance for $TAG on $(hostname)"
ARGS=(--truth "$D/mut.${SET}.projected.bed"
      --pre  "$(ls $D/annot/*.${SET}.presmoothed.bed | head -1)"
      --post "$(ls $D/annot/*.${SET}.smoothed.bed | head -1)"
      --hierarchy $ROOT/HKS_human_CHM13_v2/hierarchy.tsv
      --feature-set "$SET" --fai "$D/mut.fa.fai" --gaps "$D/mut.n_regions.bed"
      --kmer-size 31 --out "$RES/${SET}")
[ "$SET" = repeat ] && ARGS+=(--merge "$MERGE")
bash ${KS_WORK}/ks_concord.sh "${ARGS[@]}"
echo "rc=$?"
[ -s "$RES/repeat.summary.tsv" ] && [ -s "$RES/region.summary.tsv" ] \
  && { echo "both sets now present -> freeing intermediates"; rm -f $D/annot/*.bed $D/mut.fa; } \
  || echo "still incomplete"
