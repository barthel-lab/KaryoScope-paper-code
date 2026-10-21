#!/bin/bash
# =============================================================================
# A2 Step 3 -- the k-sweep (R3-major-3: are the BIOLOGICAL results robust to k?)
# =============================================================================
# 18 tasks = 3 ROB samples x k in {11,21,31,41,51,61}.
#
# k=31 is run through the SAME driver as every other k, as a control: it must
# reproduce the shipped pipeline's gate answers (run_step3_gate.sh) exactly, to
# the base, in the scaffolded frame.
# If it does not, the driver is wrong and no other row can be trusted. It also
# checks the k=31 sweep rebuild against the shipped v2 database.
#
# Per task: annotate hap1 with HKS_human_CHM13_v2_k$K (chromosome set only),
# remap onto the FIXED k=31 scaffold map, then run the shipped
# refine_breakpoint + cumulative_breakpoint at pivot=chromosome.
# See sweep_refine.py for what is held fixed and why.
#SBATCH --job-name=a2s3_sweep
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --array=1-18%6
#SBATCH --output=a2s3_sweep_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_ISCN_SRC:?set KS_ISCN_SRC to vendored/KaryoScope-ISCN/src in this repository}"
export KS_ISCN_SRC
: "${KS_GUARRACINO:?set KS_GUARRACINO (see the repository README, Paths)}"
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
set -uo pipefail
W=${KS_WORK}/a2_step3
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
SWEEPDB=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db
SRC=${KS_GUARRACINO}/ks-robs/results
ISCN=${KS_SOFTWARE}/github/barthel/KaryoScope-ISCN

SAMPLES=(GM03417 GM03786 GM04890)
KS=(11 21 31 41 51 61)
I=$((SLURM_ARRAY_TASK_ID-1))
S=${SAMPLES[$((I / 6))]}
K=${KS[$((I % 6))]}
DB=HKS_human_CHM13_v2_k$K

OUT=$W/sweep/$S/k$K; mkdir -p "$OUT"
# Keep ann/ -- annotate is the expensive step and all 18 completed it on the
# first attempt (the failures were downstream). Everything else is regenerated.
rm -f "$OUT"/{row.tsv,chromosome.scaffolded.bed.gz,chromosome.smoothed.clipped.bed.gz,clip.log,remap.log}
FA=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
MAP=$W/$S/work/$S.haplotype1.HKS_human_CHM13_v2.scaffold_map.tsv   # FIXED k=31 layout

{
  echo "=== A2 Step 3 sweep: $S k=$K ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  echo "db   : $DB"
  echo "map  : $MAP  (k=31 layout, held fixed across k)"
  echo "iscn : $(git -C $ISCN rev-parse --short HEAD) on $(git -C $ISCN rev-parse --abbrev-ref HEAD)"
} | tee "$OUT/provenance.txt"

[ -d "$SWEEPDB/$DB" ] || { echo "FATAL: no database $SWEEPDB/$DB"; exit 1; }
[ -s "$MAP" ]         || { echo "FATAL: no scaffold map $MAP"; exit 1; }

ANN="$OUT/ann/$S.haplotype1.$DB.chromosome.smoothed.bed.gz"
if [ -s "$ANN" ]; then
  echo "[$(date +%T)] annotate output already present, reusing: $(stat -c%s "$ANN") bytes"
else
echo "[$(date +%T)] annotate (chromosome feature set only)"
/usr/bin/time -v $KSBIN/karyoscope annotate \
    --input "$FA" --outdir "$OUT/ann" --db "$DB" --db-root "$SWEEPDB" \
    --feature-set chromosome --smooth --threads 8 \
    > "$OUT/annotate.log" 2> "$OUT/annotate.time"
RC=$?
echo "annotate_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/annotate.time" | awk '{print $NF}') peakRSS=$(grep 'Maximum resident' "$OUT/annotate.time" | awk '{print $NF}')KB"
[ $RC -eq 0 ] || { echo "FATAL: annotate failed"; exit $RC; }
fi

# NAME THE FILE EXACTLY. A glob like *chromosome*smoothed*.bed.gz sorts
# "presmoothed" BEFORE "smoothed", so `| head -1` silently selects the wrong
# input and every downstream number is quietly wrong (the k=31 control caught
# exactly this on the first attempt). Same class as the ls|head -1 trap.
SM="$OUT/ann/$S.haplotype1.$DB.chromosome.smoothed.bed.gz"
[ -s "$SM" ] || { echo "FATAL: no smoothed chromosome BED at $SM"; ls -la "$OUT/ann"; exit 1; }
echo "smoothed: $SM  ($(stat -c%s "$SM") bytes)"

# The BED's extent is k-dependent: k-mer starts run 0..L-k, so a k=11 BED reaches
# (31-11)=20 bp further than the k=31 scaffold map's recorded contig length and
# remap-bed rejects it as "annotated against a different assembly". Clip to the
# map's lengths so every k shares ONE layout. The effect is <= 30 bp at contig
# ends, nowhere near a breakpoint.
echo "[$(date +%T)] clip to the k=31 map's contig lengths"
CLIP="$OUT/chromosome.smoothed.clipped.bed.gz"
zcat "$SM" | awk -F'\t' -v OFS='\t' -v map="$MAP" '
  BEGIN{ while((getline line < map) > 0){ split(line,f,"\t"); if(f[2]!="original_name") L[f[2]]=f[7]+0 } }
  { if(!($1 in L)) { print; next }
    if($2 >= L[$1]) { dropped++; next }
    if($3 >  L[$1]) { $3 = L[$1]; trimmed++ }
    print }
  END{ printf "clip: dropped %d, trimmed %d\n", dropped, trimmed > "/dev/stderr" }' \
  2> "$OUT/clip.log" | gzip > "$CLIP"
cat "$OUT/clip.log"

echo "[$(date +%T)] remap onto the fixed k=31 scaffold layout"
$KSBIN/karyoscope remap-bed --bed "$CLIP" --map "$MAP" \
    --output "$OUT/chromosome.scaffolded.bed.gz" > "$OUT/remap.log" 2>&1
RC=$?
echo "remap_rc=$RC"
[ $RC -eq 0 ] || { echo "FATAL: remap failed"; tail -20 "$OUT/remap.log"; exit $RC; }

echo "[$(date +%T)] shipped refine + cumulative at pivot=chromosome"
$KSBIN/python $W/sweep_refine.py "$S" "$K" "$OUT/chromosome.scaffolded.bed.gz" "$OUT/row.tsv"
RC=$?
echo "refine_rc=$RC"
echo "[$(date +%T)] done"
exit $RC
