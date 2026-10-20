#!/bin/bash
# A2 Step 4 / Job 4a -- annotate the published centromere windows at k = 21, 31, 51.
# One array task per chromosome (24). Windows = the per-contig extent of the
# published pangenome.<chr>.centromere.*.region.bed (forward frame; the published
# scripts dropped reverse-complemented contigs), restricted to contigs that pass
# the NucFlag QC filter and whose span is <= 50 Mb -- both are hard filters the
# published pipeline applies downstream, so dropping them here changes nothing
# and saves 27 Gb of annotation per k. Region set only, presmoothed only.
#SBATCH --job-name=a2s4_cenv1
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-24%12
#SBATCH --output=cenv1_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job4_centromere; C=${KS_WORK}/a2_step4/common
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
DBROOT=${KS_WORK}/ks_build_sources/v1like_ksweep/db
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
CH=${CHROMS[$((SLURM_ARRAY_TASK_ID-1))]}
mkdir -p $W/fasta $W/prov $W/extents
PROV=$W/prov/prov_v1like_$CH.txt
{
  echo "=== A2 step4 job4a $CH ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[$SLURM_ARRAY_TASK_ID]  node: $(hostname)  $(lscpu | grep -m1 'Model name' | sed 's/ \+/ /g')"
  echo "karyoscope: $(karyoscope --version)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')  threads: 16"
} | tee "$PROV"
FA=$W/fasta/windows_$CH.fa
if [ ! -s "$FA.done" ]; then
  echo "[$(date +%T)] extracting $CH windows" | tee -a "$PROV"; : > "$FA"
  awk -F'\t' -v c=$CH '$2==c && $9==1' $C/cen_windows.tsv | while IFS=$'\t' read -r pansn ch mn mx span fasta L qc use; do
    samtools faidx "$fasta" "${pansn}:$((mn+1))-${mx}" | sed "1s/.*/>${pansn}/" >> "$FA"
  done
  samtools faidx "$FA" && touch "$FA.done"
fi
echo "  $CH: $(grep -c '>' $FA) windows, $(awk '{s+=$2}END{print s}' $FA.fai) bp" | tee -a "$PROV"
for K in 21 31 51; do
  DB=HKS_human_CHM13_v1like_region_k$K; A=$W/annot_v1like_k$K/$CH; mkdir -p "$A" $W/beds_v1like_k$K
  OUTBED=$W/beds_v1like_k$K/pangenome.$CH.centromere.KS_human_CHM13.presmoothed.region.bed
  if [ -s "$OUTBED.done" ]; then echo "  k=$K already done"; continue; fi
  echo "[$(date +%T)] annotate $CH k=$K" | tee -a "$PROV"
  /usr/bin/time -v karyoscope annotate -i "$FA" -o "$A" --db $DB --db-root $DBROOT --feature-set region \
      -t 16 --no-bgzip --no-smooth --force > $W/logs/annotate_v1like_${CH}_k$K.log 2> $W/logs/time_v1like_${CH}_k$K.txt
  RC=$?
  echo "  annotate rc=$RC $(grep Elapsed $W/logs/time_v1like_${CH}_k$K.txt | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $W/logs/time_v1like_${CH}_k$K.txt | awk '{print $NF}')" | tee -a "$PROV"
  [ $RC -eq 0 ] || { echo "FATAL annotate $CH k=$K"; exit 1; }
  BED=$A/windows_$CH.$DB.region.presmoothed.bed
  echo "[$(date +%T)] bin 100 kb + trim $CH k=$K" | tee -a "$PROV"
  python3 $W/cen_postprocess.py --bed "$BED" --windows $C/cen_windows.tsv --chrom $CH \
      --bin-script $C/bin_features_feb2026.py --out "$OUTBED" --extents-out $W/extents/extents_v1like_${CH}_k$K.tsv | tee -a "$PROV" \
    && touch "$OUTBED.done" && rm -f "$BED"
done
echo "date_end: $(date -Is)" | tee -a "$PROV"
