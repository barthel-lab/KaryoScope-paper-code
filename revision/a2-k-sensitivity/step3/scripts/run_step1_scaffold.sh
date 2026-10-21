#!/bin/bash
# A2 Step 3, step 1 (reproduction path): the k=31 scaffold layout that the k-sweep holds
# fixed, made with public KaryoScope.
#
# The original runs got this layout from `karyoscope-iscn analyze` (run_step1_baseline.sh),
# which annotates each haplotype with `karyoscope annotate` and then calls
# `karyoscope scaffold --mode bed --feature-set region`. This script runs those two
# KaryoScope commands directly for haplotype 1, with the same input label (hap1), so the
# scaffold map is byte-identical to the one analyze wrote. run_step3_sweep.sh reads
# $KS_WORK/a2_step3/<sample>/work/<sample>.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv.
#SBATCH --job-name=a2s3_scaffold
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-3
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_GUARRACINO:?set KS_GUARRACINO (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
SRC=${KS_GUARRACINO}/ks-robs/results
SAMPLES=(GM03417 GM03786 GM04890)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
OUT=$W/$S/work; mkdir -p "$OUT"
FA=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
[ -r "$FA" ] || { echo "FATAL: unreadable $FA"; exit 1; }

echo "=== A2 Step 3 scaffold layout: $S ==="
echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
echo "karyoscope: $(karyoscope --version 2>&1 | tail -1)"

for FS in region chromosome; do
  karyoscope annotate --input "$FA" --db HKS_human_CHM13_v2 --feature-set $FS \
      --outdir "$OUT" --threads 16 --db-root "$ROOT" --no-keep-presmoothed
  echo "annotate $FS rc=$?"
done
karyoscope scaffold --input "hap1=$FA" --db HKS_human_CHM13_v2 --feature-set region --mode bed \
    --outdir "$OUT" --threads 16 --db-root "$ROOT"
echo "scaffold rc=$?"
ls -l "$OUT/$S.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv"
