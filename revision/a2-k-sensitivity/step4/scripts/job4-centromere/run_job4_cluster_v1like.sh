#!/bin/bash
# A2 Step 4 / Job 4b -- Yi-An's centromere pipeline, unchanged, on each k's inputs.
# Steps 1-3 (concatenate, NucFlag QC filter, per-chromosome clustering) and the
# published gate (--sil-threshold 0.5 --centroid-sd 5) via KS_allchr_barplot.py,
# then per-haplotype labels. Scripts: KaryoScope-analysis@pangenome_structure_V2
# (90e2874), exported to a6_sweeps/cent/scripts. The PUBLISHED-input control is
# the A6 run (job 37136128) whose assignments are reused and re-gated identically.
#SBATCH --job-name=a2s4_cencluv1
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=cencluv1_%j.log
: "${KS_HOME:?set KS_HOME (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job4_centromere
S=${KS_WORK}/a6_sweeps/cent/scripts
QC="${KS_WORK}/aggregate_qc_v4.manual_curation (2).tsv"
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
export PATH=${KS_HOME}/software/bin:$PATH
NAME=pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed
PROV=$W/prov/prov_cluster_v1like.txt; mkdir -p $W/prov $W/labels
{
  echo "=== A2 step4 job4b clustering ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}  node: $(hostname)"
  echo "scripts: KaryoScope-analysis@pangenome_structure_V2 90e2874 (a6_sweeps/cent/scripts)"
  echo "gate: --sil-threshold 0.5 --centroid-sd 5 (centroid scan on)"
} | tee "$PROV"
for K in 21 31 51; do KTAG=v1like_$K
  D=$W/cluster_v1like_k$K; rm -rf "$D"; mkdir -p $D/centromere_region_beds $D/run $D/out
  n=$(ls $W/beds_v1like_k$K/pangenome.chr*.centromere.KS_human_CHM13.presmoothed.region.bed 2>/dev/null | wc -l)
  [ "$n" -eq 24 ] || { echo "FATAL: k=$K has $n/24 per-chromosome BEDs" | tee -a "$PROV"; exit 1; }
  ln -sf $W/beds_v1like_k$K/pangenome.chr*.centromere.KS_human_CHM13.presmoothed.region.bed $D/centromere_region_beds/
  echo "[$(date +%T)] k=$K step 1 concatenate" | tee -a "$PROV"
  ( cd $D/run && Rscript $S/concatenate_centromere_regions.R && mv $NAME $D/out/$NAME ) || { echo FATAL step1; exit 1; }
  echo "  rows: $(wc -l < $D/out/$NAME)  contigs: $(cut -f1 $D/out/$NAME | sort -u | wc -l)" | tee -a "$PROV"
  echo "[$(date +%T)] k=$K step 2 QC filter" | tee -a "$PROV"
  ( cd $D/out && Rscript $S/QCfilter_apply.R $D/out/$NAME "$QC" centromere Err=0 COLLAPSE=0 COLLAPSE_VAR=0 ) || { echo FATAL step2; exit 1; }
  PASS=$D/out/${NAME%.bed}.pass.bed
  echo "  pass rows: $(wc -l < $PASS)  contigs: $(cut -f1 $PASS | sort -u | wc -l)" | tee -a "$PROV"
  [ "$(wc -l < $PASS)" -gt 1000 ] || { echo "FATAL: empty pass.bed (contig join failed)"; exit 1; }
  echo "[$(date +%T)] k=$K step 3 clustering" | tee -a "$PROV"
  /usr/bin/time -v $PY $S/KaryoScope_cluster_analysis.py --bed "$PASS" --output-prefix $D/out/allchr_structure \
      --analysis-mode structure --edges directional --no-abundance --max-sequence-length 50000000 \
      --exclude-features "novel" --matrix-type count_log1p_zscore_blockweight --k-selection silhouette \
      --background white > $D/out/step3.log 2> $D/out/step3.time
  echo "  step3 rc=$? $(grep Elapsed $D/out/step3.time | sed 's/.*: //')  assignments: $(wc -l < $D/out/allchr_structure.sequence_assignments.tsv)" | tee -a "$PROV"
  $PY $S/KS_allchr_barplot.py --assignments $D/out/allchr_structure.sequence_assignments.tsv --bed "$PASS" \
      --output $D/out/barplot_sil0.5_sd5.svg --matrix-type count_log1p_zscore --sil-threshold 0.5 --centroid-sd 5 \
      > $D/out/barplot_sil0.5_sd5.stdout 2> $D/out/barplot.stderr
  grep -E "Grand total|chr[0-9XY]+:" $D/out/barplot_sil0.5_sd5.stdout | tee -a "$PROV"
  $PY $W/cen_labels.py --assignments $D/out/allchr_structure.sequence_assignments.tsv --bed "$PASS" --k v1_$K --out $W/labels/labels_v1_$K.tsv | tee -a "$PROV"
done
# published-input control: A6's step-3 output on the v1.1.0 annotations, re-gated identically
A6=${KS_WORK}/a6_sweeps/cent/out
$PY $W/cen_labels.py --assignments $A6/allchr_structure.sequence_assignments.tsv \
    --bed $A6/${NAME%.bed}.pass.bed --k published --out $W/labels/labels_published.tsv | tee -a "$PROV"
echo "[$(date +%T)] collate" | tee -a "$PROV"
echo 'labels written; collate separately'
echo "date_end: $(date -Is)" | tee -a "$PROV"
