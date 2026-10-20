#!/bin/bash
# A4 anchor check: does ks_concord.sh reproduce the PUBLISHED HG002 hap1 repeat
# numbers (Supp Note 1: exact 84.6 -> 95.25%, compatible 9.2 -> ~1%,
# novel 2.5 -> <0.05%)? If not, the wrapper is wrong and A5 inherits the error.
#SBATCH --job-name=a4_anchor
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=a4_anchor_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
A=${KS_SOFTWARE}/github/barthel/KaryoScope-archive
W=${KS_WORK}/a4_smoothing
RES=$W/results; mkdir -p "$RES"
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
echo "node: $(hostname)  start: $(date -Is)"
bash ${KS_WORK}/ks_concord.sh \
  --truth ${KS_WORK}/a4_smoothing/HG002.hap1.RepeatMasker.renamed.bed \
  --pre  $W/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.repeat.presmoothed.bed \
  --post $W/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.repeat.smoothed.bed \
  --hierarchy ${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv \
  --feature-set repeat \
  --fai ${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz.fai \
  --gaps ${KS_WORK}/a4_smoothing/HG002.hap1.n_regions.renamed.bed \
  --kmer-size 31 \
  --merge "Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat" \
  --out $RES/anchor_hg002hap1_repeat
rc=$?
echo "concord_rc=$rc"
echo "=============== ANCHOR COMPARISON ==============="
echo "published (Supp Note 1): presmoothed exact 84.6 / compat 9.2 / novel 2.5"
echo "                         smoothed    exact 95.25 / compat ~1 / novel <0.05"
echo "measured:"
# STALE GUARD: a killed wrapper used to leave the previous run's summary in
# place, and this script cat'd it as if it were fresh (that is how job
# 35379578 reported run 1's numbers). Refuse to print anything unless the
# file was written by THIS invocation.
[ -s $RES/anchor_hg002hap1_repeat.summary.tsv ] || { echo "FATAL: no summary produced this run"; exit 3; }
cat $RES/anchor_hg002hap1_repeat.summary.tsv 2>/dev/null
echo "================================================="
exit $rc
