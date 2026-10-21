#!/bin/bash
# Test the claim that the 99% CI-width gap (ours vs the published widths,
# a2_step3 published_ci_width_mb in the ROB config) is explained by presmoothed-vs-smoothed input.
# The legacy breakpoint_density.py read the PRESMOOTHED scaffolded chromosome
# BED; the current pipeline reads the SMOOTHED one. We still have the
# presmoothed BEDs (unscaffolded) and the scaffold maps, so remap and re-run
# the identical analysis on both inputs.
#SBATCH --job-name=a2s3_presm
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_presm_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$KSBIN:$PATH
SAMPLES=(GM03417 GM03786 GM04890)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
SRC=$W/$S/work
OUT=$W/presm_$S; mkdir -p "$OUT"

BED=$SRC/$S.haplotype1.HKS_human_CHM13_v2.chromosome.presmoothed.bed.gz
MAP=$SRC/$S.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv
DST=$OUT/$S.haplotype1.HKS_human_CHM13_v2.chromosome.presmoothed.scaffolded.bed.gz

echo "=== $S  $(date -Is)  node $(hostname) ==="
ls -la "$BED" "$MAP"
if [ -s "$DST" ]; then
  echo "remap already present, reusing"
else
  echo "[$(date +%T)] remap-bed presmoothed -> scaffolded"
  /usr/bin/time -v $KSBIN/karyoscope remap-bed --bed "$BED" --map "$MAP" --output "$DST" \
      > "$OUT/remap.log" 2> "$OUT/remap.time"
  echo "remap_rc=$? wall=$(grep 'Elapsed (wall' "$OUT/remap.time" | awk '{print $NF}')"
fi
ls -la "$DST"
echo "[$(date +%T)] done"
