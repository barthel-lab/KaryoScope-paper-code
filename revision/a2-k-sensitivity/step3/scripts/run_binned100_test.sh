#!/bin/bash
# =============================================================================
# Does the published 200-525 kb come from WINNER-TAKE-ALL 100 bp BINNED input?
# =============================================================================
# The paper-repo session found that both published scripts read
#   {assembly}.{db}.chromosome.presmoothed.scaffolded.binned100.features.bed.gz
# (breakpoint_density.py:133, breakpoint_density_stacked.py:351; robs.yaml:58
# declares the same template). CONFIRMED on this side: every reproduction so far
# used the raw scaffolded BED (smoothed, and presmoothed in the CI test) with NO
# `karyoscope bin` step. So this axis has never been tested.
#
# Mechanism: `bin` is winner-take-all, so a 100 bp bin holding 60 bp chr14 +
# 40 bp chr13 contributes 100/0 instead of 60/40. The discarded minority signal
# is concentrated in the fusion transition zone -- exactly what localises the
# breakpoint -- so it should move BOTH the argmax and the plateau flatness,
# i.e. both the delta and the 99% CI width. One cause, both discrepancies.
#
# Prediction if true (figure-derived, addendum 3):
#   deltas and 99% CI widths as read from the submitted figure (a2_step3
#   published_delta_kb, published_ci_width_mb in the ROB config), against the
#   unbinned per-base values
#SBATCH --job-name=a2s3_bin100
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_bin100_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
SAMPLES=(GM03417 GM03786 GM04890)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
OUT=$W/bin100_$S; mkdir -p "$OUT"

# The presmoothed SCAFFOLDED BEDs already exist (run_presmoothed_ci_test.sh).
SRCBED=$W/presm_$S/$S.haplotype1.HKS_human_CHM13_v2.chromosome.presmoothed.scaffolded.bed.gz
[ -s "$SRCBED" ] || { echo "FATAL: missing $SRCBED"; exit 1; }
echo "=== binned100 test: $S  $(date -Is)  node $(hostname) ==="
echo "source: $SRCBED ($(stat -c%s "$SRCBED") bytes)"

# Variant A: with leaf prioritisation (--db/--feature-set), as the docs describe.
if [ -s "$OUT/binned100.prio.bed.gz" ]; then echo "variant A present, reusing"; else
echo "[$(date +%T)] bin 100 bp, WITH leaf prioritisation"
/usr/bin/time -v $KSBIN/karyoscope bin -i "$SRCBED" -o "$OUT/binned100.prio.bed.gz" \
    -b 100 --db HKS_human_CHM13_v2 --db-root "$ROOT" --feature-set chromosome -t 8 \
    > "$OUT/binA.log" 2> "$OUT/binA.time"
echo "binA_rc=$? wall=$(grep 'Elapsed (wall' "$OUT/binA.time" | awk '{print $NF}')"
fi

# Variant B: no prioritisation flags -- pure winner-take-all by overlap bp.
if [ -s "$OUT/binned100.plain.bed.gz" ]; then echo "variant B present, reusing"; else
echo "[$(date +%T)] bin 100 bp, WITHOUT prioritisation flags"
/usr/bin/time -v $KSBIN/karyoscope bin -i "$SRCBED" -o "$OUT/binned100.plain.bed.gz" \
    -b 100 -t 8 > "$OUT/binB.log" 2> "$OUT/binB.time"
echo "binB_rc=$? wall=$(grep 'Elapsed (wall' "$OUT/binB.time" | awk '{print $NF}')"
fi

for v in prio plain; do
  f="$OUT/binned100.$v.bed.gz"
  [ -s "$f" ] && echo "  $v: $(stat -c%s "$f") bytes" || echo "  $v: MISSING"
done
echo "[$(date +%T)] done"
