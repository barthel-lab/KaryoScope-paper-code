#!/bin/bash
# A6 §6 — TAR1 / ITS merge-distance sweep.
#
# main.tex:338 (TAR1) and :340 (ITS) both use a 1 kb merge distance with the
# same stated rationale: consolidate RepeatMasker's fragmented per-element hits
# into one feature per array. So the evidence wanted is a KNEE -- interval count
# falling steeply then flattening, with 1 kb past the knee -- and total bases
# labelled staying nearly FLAT. If total bases move much, the threshold is
# absorbing real sequence, which is the problematic outcome.
#
# BED layout (verified): chrom start end name score strand class family div id
#   TAR1 = field 4 contains "TAR1"
#   ITS  = field 7 == Simple_repeat AND field 4 contains a TTAGGG/CCCTAA rotation,
#          THEN dropped if it overlaps a seqtk telo region, THEN merged.
# The filter-before-merge ordering is load-bearing and is asserted below: merging
# first would let an interior hit bridge into a terminal array.
#SBATCH --job-name=a6_tar1_sweep
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=sweep_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
W=${KS_WORK}/a6_sweeps/tar1_its
RM=$W/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
TELO=$W/chm13v2.0.telo.bed
RES=$W/results; mkdir -p "$RES"
echo "=== TAR1/ITS sweep  $(date -Is)  node $(hostname) ==="

# --- candidate sets -------------------------------------------------------
awk -F'\t' '$4 ~ /TAR1/ {print $1"\t"$2"\t"$3}' "$RM" | sort -k1,1 -k2,2n > $RES/tar1.raw.bed
ROT='TTAGGG|TAGGGT|AGGGTT|GGGTTA|GGTTAG|GTTAGG|CCCTAA|CCTAAC|CTAACC|TAACCC|AACCCT|ACCCTA'
awk -F'\t' -v rot="$ROT" '$7=="Simple_repeat" && $4 ~ rot {print $1"\t"$2"\t"$3}' "$RM" \
  | sort -k1,1 -k2,2n > $RES/its.candidates.bed
sort -k1,1 -k2,2n "$TELO" | cut -f1-3 > $RES/telo.sorted.bed
$BT intersect -v -a $RES/its.candidates.bed -b $RES/telo.sorted.bed > $RES/its.raw.bed

echo "  TAR1 hits            : $(wc -l < $RES/tar1.raw.bed)"
echo "  ITS candidates       : $(wc -l < $RES/its.candidates.bed)"
echo "  ITS after telo filter: $(wc -l < $RES/its.raw.bed)"
echo "  (removed by filter   : $(( $(wc -l < $RES/its.candidates.bed) - $(wc -l < $RES/its.raw.bed) )))"

# --- ORDERING CHECK: filter-then-merge vs merge-then-filter ---------------
$BT merge -d 1000 -i $RES/its.candidates.bed | $BT intersect -v -a - -b $RES/telo.sorted.bed > $RES/its.wrong_order.bed
echo "=== ordering check at d=1000 ==="
echo "  filter-then-merge: $($BT merge -d 1000 -i $RES/its.raw.bed | wc -l) intervals, $($BT merge -d 1000 -i $RES/its.raw.bed | awk '{s+=$3-$2} END{print s+0}') bp"
echo "  merge-then-filter: $(wc -l < $RES/its.wrong_order.bed) intervals, $(awk '{s+=$3-$2} END{print s+0}' $RES/its.wrong_order.bed) bp"

# --- sweep ----------------------------------------------------------------
printf "feature\tmerge_d\tintervals\ttotal_bp\tmin\tq1\tmedian\tq3\tmax\n" > $RES/tar1_its_sweep.tsv
for SET in tar1 its; do
  for D in 0 100 500 1000 2000 5000 10000; do
    $BT merge -d $D -i $RES/${SET}.raw.bed > $RES/${SET}.merged_d${D}.bed
    awk -F'\t' -v set=$SET -v d=$D '
      {n++; L[n]=$3-$2; s+=$3-$2}
      END{ if(n==0){printf "%s\t%s\t0\t0\t0\t0\t0\t0\t0\n",set,d; exit}
           asort(L)
           printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", set, d, n, s,
                  L[1], L[int(n*0.25)+((n*0.25)==int(n*0.25)?0:1)],
                  L[int(n*0.5)+((n*0.5)==int(n*0.5)?0:1)],
                  L[int(n*0.75)+((n*0.75)==int(n*0.75)?0:1)], L[n] }' \
      $RES/${SET}.merged_d${D}.bed >> $RES/tar1_its_sweep.tsv
  done
done
echo "=== SWEEP ==="
column -t -s$'\t' $RES/tar1_its_sweep.tsv
echo "done $(date -Is)"
