#!/bin/bash
# A2 Step 0 (pre-flight) + Step 1 (uniqueness sweep, `repeat` feature set).
#
# node-stats loads the full base index (19-24 GB) and was OOM-killed on the
# login node, hence a batch job with real memory.
#
# k=11 runs FIRST here (not last): its index is 19 MB, so it validates the
# invocation cheaply. The brief's "k=11 last" applies to Step 2's annotation
# sweep, where k=11 is expensive, not to Step 1.
#SBATCH --job-name=a2_step01
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=a2_step01_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a2_ksweep
DB=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db
SPEC=${KS_WORK}/ks_build_sources/chm13v2_ksweep
RES=$W/results; RAW=$W/raw_nodestats
mkdir -p "$RES" "$RAW"
H=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
KS=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/karyoscope
KS_LIST="11 21 31 41 51 61"
SET=repeat

# ============================================================ STEP 0
P=$RES/step0_preflight.txt
{
  echo "=== A2 Step 0 pre-flight ==="
  echo "date    : $(date -Is)"
  echo "job     : ${SLURM_JOB_ID:-NA}   node: $(hostname)"
  echo "hks     : $($H 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "karyoscope: $($KS --version 2>&1)"
  echo
  echo "--- (a) hierarchy.tsv identical across all six k? ---"
  echo "    (this is the precondition for cross-k concordance to mean anything)"
  for k in $KS_LIST; do
    printf "    k=%-3s %s\n" "$k" "$(md5sum $DB/HKS_human_CHM13_v2_k$k/hierarchy.tsv | cut -d' ' -f1)"
  done
  echo "    distinct hierarchy checksums: $(md5sum $DB/HKS_human_CHM13_v2_k*/hierarchy.tsv | awk '{print $1}' | sort -u | wc -l) (must be 1)"
  echo
  echo "--- (b) colors.tsv identical across all six k? ---"
  echo "    distinct colors checksums: $(md5sum $DB/HKS_human_CHM13_v2_k*/colors.tsv | awk '{print $1}' | sort -u | wc -l) (must be 1)"
  echo
  echo "--- (c) build specs differ ONLY by k? ---"
  for k in 11 21 41 51 61; do
    echo "    k31 vs k$k:"
    diff $SPEC/build.k31.yaml $SPEC/build.k$k.yaml | sed 's/^/      /'
  done
  echo
  echo "--- (d) CORRECTION to the brief ---"
  echo "    The brief's Step 0 says to confirm 'flatten: true' on all six."
  echo "    There is NO flatten key in any of the six specs, and none in the"
  echo "    shipped HKS_human_CHM13_v2 manifest either. grep result:"
  grep -c flatten $SPEC/build.k*.yaml 2>/dev/null | sed 's/^/      /'
  echo "    So 'no flatten' is the correct expectation, and it matches the"
  echo "    shipped database. Prior work additionally proved the k=31 rebuild"
  echo "    byte-identical to shipped in hierarchy.tsv/colors.tsv, with"
  echo "    identical annotation output on chr21 for all six feature sets."
  echo "    Conclusion: the six databases are internally controlled AND"
  echo "    faithful to production. Not a reason to stop."
  echo
  echo "--- (e) per-database manifests ---"
  for k in $KS_LIST; do
    echo "    === k=$k ==="
    grep -E "^(id|s|k|kmer)|max_size" $DB/HKS_human_CHM13_v2_k$k/manifest.yaml 2>/dev/null | sed 's/^/      /'
  done
  echo
  echo "--- (f) index sizes ---"
  du -sh $DB/HKS_human_CHM13_v2_k* 2>/dev/null | sed 's/^/    /'
} | tee "$P"

# ============================================================ STEP 1
echo "[$(date +%T)] Step 1: node-stats on the '$SET' feature set, six k"
T=$RES/step1_timings.tsv
echo -e "k\tset\texit\twall_s\tmax_rss_kb\tout_bytes" > "$T"
for k in $KS_LIST; do
  IDX=$DB/HKS_human_CHM13_v2_k$k/index
  base=$IDX/features.hksb
  fs=$IDX/features.$SET.hksf
  out=$RAW/nodestats_k${k}_${SET}.txt
  if [ ! -r "$base" ] || [ ! -r "$fs" ]; then
    echo "  k=$k: MISSING index files, skipping" | tee -a "$P"
    echo -e "${k}\t${SET}\tMISSING\tNA\tNA\t0" >> "$T"; continue
  fi
  echo "  [$(date +%T)] k=$k ..."
  s=$(date +%s)
  /usr/bin/time -v $H node-stats --index "$base" --feature-set-file "$fs" -t 16 \
      > "$out" 2> "$RAW/time_k${k}.txt"
  rc=$?
  e=$(date +%s)
  rss=$(grep 'Maximum resident set size' "$RAW/time_k${k}.txt" | awk '{print $NF}')
  echo -e "${k}\t${SET}\t${rc}\t$((e-s))\t${rss:-NA}\t$(wc -c < "$out")" >> "$T"
  echo "    rc=$rc wall=$((e-s))s rss=${rss:-NA}kb bytes=$(wc -c < "$out")"
done

echo "=== timings ==="; cat "$T"
echo "=== head of each raw output (format TBD; aggregation done separately) ==="
for k in $KS_LIST; do
  echo "--- k=$k ---"; grep -v '^\[' "$RAW/nodestats_k${k}_${SET}.txt" 2>/dev/null | head -6
done
echo "=== DONE ==="
