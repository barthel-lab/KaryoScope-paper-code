#!/bin/bash
# A2 Step 4 / Job 1 -- novel fraction vs k, from the EXISTING A2 Step 2 BEDs.
# One array task per (assembly x k); no annotation is run. k=11 is included
# (its BEDs exist and finished in 24 min), given more time because its
# presmoothed BEDs are 8-37 GB each.
#   tasks 1-6  = chm13     x k(11,21,31,41,51,61)
#   tasks 7-12 = hg002hap1 x k(11,21,31,41,51,61)
#SBATCH --job-name=a2s4_novel
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --array=1-12%12
#SBATCH --output=novel_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job1_novel
SRC=${KS_WORK}/a2_ksweep
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
T=${SLURM_ARRAY_TASK_ID}
KS_LIST=(11 21 31 41 51 61)
if [ "$T" -le 6 ]; then ASM=chm13; K=${KS_LIST[$((T-1))]}; NBED=""
else ASM=hg002hap1; K=${KS_LIST[$((T-7))]}; NBED=${KS_WORK}/a4_smoothing/HG002.hap1.n_regions.renamed.bed; fi
mkdir -p $W/results
PROV=$W/results/prov_${ASM}_k${K}.txt
{
  echo "=== A2 step4 job1 ${ASM} k=${K} ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[$T]  node: $(hostname)"
  echo "input annotation dir: $SRC/annot_${ASM}_k${K} (A2 Step 2, job 35382788 / k11 35382789, db HKS_human_CHM13_v2_k${K})"
  cat $SRC/results/prov_${ASM}_k${K}.txt 2>/dev/null | grep -E "karyoscope|db:" 
  echo "n-mask: ${NBED:-none (gapless)}  exclude: chrM  left-extension: $((K-1)) bp"
} | tee "$PROV"
/usr/bin/time -v $PY $W/novel_stats.py --annot-dir $SRC/annot_${ASM}_k${K} --assembly $ASM --k $K \
    --nbed "$NBED" --exclude chrM --out $W/results/novel_${ASM}_k${K}.tsv 2> $W/results/time_${ASM}_k${K}.txt
RC=$?
grep -E "Elapsed|Maximum resident" $W/results/time_${ASM}_k${K}.txt | tee -a "$PROV"
echo "rc=$RC date_end: $(date -Is)" | tee -a "$PROV"
exit $RC
