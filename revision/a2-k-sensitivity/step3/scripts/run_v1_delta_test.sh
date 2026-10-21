#!/bin/bash
# =============================================================================
# Does the published "200-525 kb" SST1 range come from the v1 database?
# =============================================================================
# conf/robs.yaml declares `database: KS_human_CHM13` -- the **v1** database --
# while all of A2 Step 3 runs on HKS_human_CHM13_v2. v1 is a KMC index at k=31
# whose chromosome labels are `chr13_specific`/`chr14_specific`/`chr21_specific`,
# exactly the vocabulary the legacy breakpoint_density.py hardcodes. That is a
# strong hint the published sentence was computed on v1 annotations.
#
# This re-runs the published method (100 bp, plateau_center, focus window =
# _robertsonian ISCN CI +/- 1 Mb) on v1 annotations, so the resulting three
# distances are directly comparable to "200-525 kb".
#SBATCH --job-name=a2s3_v1
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_v1_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_GUARRACINO:?set KS_GUARRACINO (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
SRC=${KS_GUARRACINO}/ks-robs/results
SAMPLES=(GM03417 GM03786 GM04890)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
OUT=$W/v1_$S; mkdir -p "$OUT"
FA=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
MAP=$W/$S/work/$S.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv

echo "=== v1 delta test: $S  $(date -Is)  node $(hostname) ==="
ANN="$OUT/ann/$S.haplotype1.KS_human_CHM13.chromosome.smoothed.bed.gz"
if [ -s "$ANN" ]; then
  echo "reusing annotate output ($(stat -c%s "$ANN") bytes)"
else
  echo "[$(date +%T)] annotate against KS_human_CHM13 (v1, KMC, k=31, chromosome only)"
  /usr/bin/time -v $KSBIN/karyoscope annotate \
      --input "$FA" --outdir "$OUT/ann" --db KS_human_CHM13 --db-root "$ROOT" \
      --feature-set chromosome --smooth --threads 16 \
      > "$OUT/annotate.log" 2> "$OUT/annotate.time"
  RC=$?
  echo "annotate_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/annotate.time" | awk '{print $NF}')"
  [ $RC -eq 0 ] || { echo "FATAL: annotate failed"; tail -30 "$OUT/annotate.log"; exit $RC; }
fi
[ -s "$ANN" ] || { echo "FATAL: no smoothed BED at $ANN"; ls -la "$OUT/ann"; exit 1; }

echo "[$(date +%T)] clip to the k=31 map lengths, then remap"
zcat "$ANN" | awk -F'\t' -v OFS='\t' -v map="$MAP" '
  BEGIN{ while((getline line < map) > 0){ split(line,f,"\t"); if(f[2]!="original_name") L[f[2]]=f[7]+0 } }
  { if(!($1 in L)){print; next} if($2>=L[$1]) next; if($3>L[$1]) $3=L[$1]; print }' \
  | gzip > "$OUT/chromosome.smoothed.clipped.bed.gz"
$KSBIN/karyoscope remap-bed --bed "$OUT/chromosome.smoothed.clipped.bed.gz" \
    --map "$MAP" --output "$OUT/chromosome.scaffolded.bed.gz" > "$OUT/remap.log" 2>&1
RC=$?; echo "remap_rc=$RC"
[ $RC -eq 0 ] || { echo "FATAL: remap failed"; tail -20 "$OUT/remap.log"; exit $RC; }

echo "[$(date +%T)] label inventory on the fusion scaffold"
zcat "$OUT/chromosome.scaffolded.bed.gz" | awk -F'\t' '{c[$4]++} END{for(l in c) print "  "l"\t"c[l]}' | sort -k2,2nr | head -12
echo "[$(date +%T)] done"
