#!/bin/bash
# =============================================================================
# A2 Step 3 k-sweep, re-run in the MANUSCRIPT'S OWN input representation
# =============================================================================
# The published 200-525 kb range comes from winner-take-all 100 bp binned input
# (an earlier job reproduced the published distances and 99% interval widths
# from it). The decision at the time was to KEEP the published numbers
# in the manuscript rather than switch to the tighter unbinned ones; the
# manuscript was later switched to the per-base numbers.
#
# Consequence: the published deliverable `sst1_breakpoints_by_k.tsv` currently
# carries UNBINNED deltas (46-252 kb), which a referee would recompute and find
# inconsistent with the paper. So the sweep must move to the manuscript's
# configuration, not the other way round -- and that is also the right test,
# since k should be the only variable within the configuration being reported.
#
# Per (sample, k): presmoothed -> clip to the k=31 map -> remap onto the FIXED
# k=31 layout -> `karyoscope bin -b 100` PLAIN (no --db/--feature-set; the
# prioritised variant does NOT reproduce the published numbers).
#SBATCH --job-name=a2s3_binsweep
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=1-18%6
#SBATCH --output=a2s3_binsweep_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$KSBIN:$PATH
SAMPLES=(GM03417 GM03786 GM04890)
KS=(11 21 31 41 51 61)
I=$((SLURM_ARRAY_TASK_ID-1))
S=${SAMPLES[$((I / 6))]}
K=${KS[$((I % 6))]}
DB=HKS_human_CHM13_v2_k$K

SRC=$W/sweep/$S/k$K/ann/$S.haplotype1.$DB.chromosome.presmoothed.bed.gz
MAP=$W/$S/work/$S.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv
OUT=$W/binsweep/$S/k$K; mkdir -p "$OUT"

echo "=== binned sweep: $S k=$K  $(date -Is)  node $(hostname) ==="
[ -s "$SRC" ] || { echo "FATAL: missing $SRC"; exit 1; }
echo "source: $SRC ($(stat -c%s "$SRC") bytes)"

echo "[$(date +%T)] clip to the k=31 map lengths"
zcat "$SRC" | awk -F'\t' -v OFS='\t' -v map="$MAP" '
  BEGIN{ while((getline line < map) > 0){ split(line,f,"\t"); if(f[2]!="original_name") L[f[2]]=f[7]+0 } }
  { if(!($1 in L)){print; next} if($2>=L[$1]) next; if($3>L[$1]) $3=L[$1]; print }' \
  | gzip > "$OUT/presmoothed.clipped.bed.gz"

echo "[$(date +%T)] remap onto the fixed k=31 layout"
$KSBIN/karyoscope remap-bed --bed "$OUT/presmoothed.clipped.bed.gz" --map "$MAP" \
    --output "$OUT/presmoothed.scaffolded.bed.gz" > "$OUT/remap.log" 2>&1
RC=$?; echo "remap_rc=$RC"; [ $RC -eq 0 ] || { tail -20 "$OUT/remap.log"; exit $RC; }

echo "[$(date +%T)] bin 100 bp, PLAIN winner-take-all (the published configuration)"
/usr/bin/time -v $KSBIN/karyoscope bin -i "$OUT/presmoothed.scaffolded.bed.gz" \
    -o "$OUT/binned100.plain.bed.gz" -b 100 -t 8 > "$OUT/bin.log" 2> "$OUT/bin.time"
RC=$?; echo "bin_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/bin.time" | awk '{print $NF}')"
[ $RC -eq 0 ] || { tail -20 "$OUT/bin.log"; exit $RC; }

rm -f "$OUT/presmoothed.clipped.bed.gz"          # reconstructible, and large
echo "  binned: $(stat -c%s "$OUT/binned100.plain.bed.gz") bytes"
echo "[$(date +%T)] done"
