#!/bin/bash
# A2 Step 2b: per-base concordance of each k against k=31, BOTH DIRECTIONS.
# Direction matters: contrast_HG002.py's 'compatible' means the prediction is an
# ANCESTOR of the truth (correct but less specific). With k=31 always as truth, a
# k>31 call that is MORE specific would score incompatible though it is not an
# error -- hence both directions.
#
# --no-fill: truth here is a KaryoScope annotation, which TILES the input
# (invariant) and is already in genome order, so there is nothing to fill and
# fill_gaps' 'nonrepeat' filler (a repeat-set concept) must not be injected.
# N regions need no special handling: CHM13 is gapless, and for HG002 both sides
# are novel in the same places, so they cannot bias an annotation-vs-annotation
# comparison.
#SBATCH --job-name=a2_concord
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=1-40%6
#SBATCH --output=a2conc_%A_%a.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
W=${KS_WORK}/a2_ksweep; RES=$W/results/concord; mkdir -p "$RES"
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
HIER=$ROOT/HKS_human_CHM13_v2/hierarchy.tsv
MERGE="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"
ASMS=(chm13 hg002hap1); KLIST=(11 21 41 51 61); SETS=(repeat region); DIRS=(fwd rev)
i=$((SLURM_ARRAY_TASK_ID-1))
D=${DIRS[$((i%2))]}; i=$((i/2))
S=${SETS[$((i%2))]}; i=$((i/2))
K=${KLIST[$((i%5))]}; i=$((i/5))
A=${ASMS[$((i%2))]}
case $A in
  chm13) FAI=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta.fai;;
  *)     FAI=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz.fai;;
esac
pick () { ls $W/annot_${1}_k${2}/*.${3}.${4}.bed 2>/dev/null | head -1; }
if [ "$D" = fwd ]; then
  TRUTH=$(pick $A 31 $S smoothed); PRE=$(pick $A $K $S presmoothed); POST=$(pick $A $K $S smoothed)
else
  TRUTH=$(pick $A $K $S smoothed); PRE=$(pick $A 31 $S presmoothed); POST=$(pick $A 31 $S smoothed)
fi
OUT=$RES/${A}_k${K}_${S}_${D}
echo "task=$SLURM_ARRAY_TASK_ID asm=$A k=$K set=$S dir=$D node=$(hostname)"
[ -n "$TRUTH" ] && [ -n "$PRE" ] && [ -n "$POST" ] || { echo "MISSING inputs"; exit 1; }
ARGS=(--truth "$TRUTH" --pre "$PRE" --post "$POST" --hierarchy "$HIER"
      --feature-set "$S" --fai "$FAI" --no-fill --out "$OUT")
[ "$S" = repeat ] && ARGS+=(--merge "$MERGE")
bash ${KS_WORK}/ks_concord.sh "${ARGS[@]}"
rc=$?
echo "rc=$rc"
[ -s "${OUT}.summary.tsv" ] && awk -v a=$A -v k=$K -v s=$S -v d=$D 'NR>1{print a"\t"k"\t"s"\t"d"\t"$0}' \
   "${OUT}.summary.tsv" >> $RES/concord_rows.tsv
exit $rc
