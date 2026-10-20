#!/bin/bash
# Job 2a -- HG002 hap1 gene set: per-label bp, both stages, A4 N-mask, denominator 2,940,260,101.
# Plus the database's own gene composition (input BED by reference position; k=31 node counts).
#SBATCH --job-name=a7_job2a
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=job2a_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a7_review_checks
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
echo "job=${SLURM_JOB_ID:-NA} node=$(hostname) start=$(date -Is)"
$PY $W/scripts/job2_gene_labels.py \
  --pre ${KS_WORK}/a4_smoothing/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.gene.presmoothed.bed \
  --post ${KS_WORK}/a4_smoothing/annot_hg002hap1/hg002v1.1_hap1.HKS_human_CHM13_v2.gene.smoothed.bed \
  --nbed ${KS_WORK}/a4_smoothing/HG002.hap1.n_regions.renamed.bed \
  --db HKS_human_CHM13_v2 --out $W/results/hg002_gene_labels_by_stage.tsv
$PY $W/scripts/job2_db_composition.py
echo "end=$(date -Is)"
