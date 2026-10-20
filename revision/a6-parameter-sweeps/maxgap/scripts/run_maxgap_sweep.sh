#!/bin/bash
# A6 §3 item 2/3 — smoothing max_gap concordance sweep.
#
# main.tex:366 says 1000 bp "was chosen by evaluating the concordance between
# the smoothed chromosome feature annotations of the HG002 assembly and the
# corresponding reference-derived ground truth". Those numbers appear nowhere.
# This produces them.
#
# Cheap by construction: annotation is NOT re-run. `hks smooth` takes the
# existing presmoothed BED and --max-gap, so each point is a re-smooth plus a
# concordance pass. Feature set is `repeat`, because that is the one whose
# ground truth (RepeatMasker) is in hand AND whose anchor is proven exact
# (job 35385021 reproduced the published Supp Note 1 numbers to the base).
#
# max_gap=0 is the control: smoothing off. It must reproduce the PRESMOOTHED
# numbers. If it does not, the harness is wrong before any threshold conclusion.
#SBATCH --job-name=a6_maxgap
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --array=1-8
#SBATCH --output=a6mg_%A_%a.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
W=${KS_WORK}/a4_smoothing
A6=${KS_WORK}/a6_sweeps
RES=$A6/results/maxgap; mkdir -p "$RES"
REG=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2
GAPS=(0 100 250 500 1000 2000 5000 10000)
G=${GAPS[$((SLURM_ARRAY_TASK_ID-1))]}
FS=repeat
PRE=$W/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.${FS}.presmoothed.bed
HIER_SET=$REG/index/features.${FS}.hierarchy.txt
OUTBED=$RES/smoothed_maxgap${G}.${FS}.bed

echo "=== A6 max_gap sweep: G=$G set=$FS node=$(hostname) start=$(date -Is) ==="
[ -s "$PRE" ] || { echo "FATAL missing presmoothed $PRE"; exit 1; }
[ -s "$HIER_SET" ] || { echo "FATAL missing hierarchy $HIER_SET"; exit 1; }

echo "[$(date +%T)] hks smooth --max-gap $G"
/usr/bin/time -v hks smooth --feature-hierarchy "$HIER_SET" -i "$PRE" -o "$OUTBED" \
    --max-gap "$G" -t 8 --miss-label novel --no-header \
    > $RES/smooth_g${G}.log 2>&1
rc=$?; echo "smooth_rc=$rc  wall=$(grep 'Elapsed (wall' $RES/smooth_g${G}.log | awk '{print $NF}')"
[ $rc -eq 0 ] || { tail -20 $RES/smooth_g${G}.log; exit $rc; }
echo "  smoothed bytes: $(stat -c%s "$OUTBED")"

echo "[$(date +%T)] concordance vs RepeatMasker truth"
bash ${KS_WORK}/ks_concord.sh \
  --truth $W/HG002.hap1.RepeatMasker.renamed.bed \
  --pre  "$PRE" \
  --post "$OUTBED" \
  --hierarchy $REG/hierarchy.tsv \
  --feature-set $FS \
  --fai ${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz.fai \
  --gaps $W/HG002.hap1.n_regions.renamed.bed \
  --kmer-size 31 \
  --merge "Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat" \
  --out $RES/maxgap${G}_${FS}
crc=$?; echo "concord_rc=$crc"
# STALE GUARD (same failure mode as the A4 anchor): refuse to report a summary
# this invocation did not write.
[ -s $RES/maxgap${G}_${FS}.summary.tsv ] || { echo "FATAL: no summary produced this run"; exit 3; }
echo "=== max_gap=$G SUMMARY ==="
cat $RES/maxgap${G}_${FS}.summary.tsv
rm -f "$OUTBED"          # regenerable in minutes; these are multi-GB
echo "done $(date -Is)"
exit $crc
