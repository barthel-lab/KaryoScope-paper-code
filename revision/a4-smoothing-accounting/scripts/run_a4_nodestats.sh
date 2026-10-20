#!/bin/bash
#SBATCH --job-name=a4_nodestats
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=a4_nodestats_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing; RAW=$W/raw_nodestats; mkdir -p "$RAW"
H=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
IDX=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/index
echo "node: $(hostname)"; ls $IDX | head
for SET in chromosome region repeat subtelomeric gene acrocentric; do
  fs=$IDX/features.$SET.hksf
  [ -r "$fs" ] || { echo "  $SET: MISSING"; continue; }
  echo "[$(date +%T)] node-stats $SET"
  /usr/bin/time -v $H node-stats --index $IDX/features.hksb --feature-set-file "$fs" -t 16 \
      > "$RAW/nodestats_shipped_${SET}.txt" 2> "$RAW/time_${SET}.txt"
  echo "  rc=$? bytes=$(wc -c < "$RAW/nodestats_shipped_${SET}.txt")"
done
echo "=== head of repeat output ==="; grep -v '^\[' "$RAW/nodestats_shipped_repeat.txt" 2>/dev/null | head -8
echo "=== DONE ==="
