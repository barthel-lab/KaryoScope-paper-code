#!/bin/bash
# A2 Step 3, step 1 RE-RUN to verify the two ISCN false-positive fixes.
#
# The baseline called del(Y)(p11.1q12) on GM03786, an XX sample. Root cause: a
# 2.48 Mb hap2 contig of X pseudo-autosomal + X-transposed sequence is laid out as
# chrY, and two isolated 100 kb fine bins labelled Yp11.2 (shared k-mers with the
# XTR) made chrY look "present, then truncated" for the remaining 51.9 Mb.
#   fix 1  detect_terminal_deletions: sustained-retention gate (>= 500 kb of the
#          chromosome's OWN bands), the terminal analogue of DEFAULT_MIN_FLANK_BP.
#   fix 2  build_sample_report: a deletion covering >= 50% of a chromosome the
#          copy-number table puts at 0 copies is the numerical call restated as an
#          SV -- suppressed, with a warning. Local events are untouched.
#
# REUSES THE CACHED INTERMEDIATES. Both fixes live at the very end of the pipeline
# (SV detection + report assembly); annotate/smooth/scaffold/bin are unchanged, and
# every pipeline step is skip-if-present. The baseline workdirs are symlinked into a
# fresh directory so the reused artefacts stay byte-identical AND the baseline
# outputs cannot be touched -- anything that did get recomputed would land here.
# ~30 min/sample of k-mer work becomes a couple of minutes of report assembly.
#
# Expected: GM03786 loses ONLY the del(Y); rob(13;14) and dic(10;X) survive with
# identical breakpoints; GM03417/GM04890 unchanged.
#SBATCH --job-name=a2s3_verify
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_verify_%A_%a.log
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
ISCN=${KS_SOFTWARE}/github/barthel/KaryoScope-ISCN

SAMPLES=(GM03417 GM03786 GM04890)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
OUT=$W/verify_$S; rm -rf "$OUT"; mkdir -p "$OUT/work"
H1=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
H2=$SRC/$S/assembly/haplotype2/$S.haplotype2.fasta.gz

# Symlink every cached intermediate under the name the pipeline predicts.
BASE=$W/$S/work
[ -d "$BASE" ] || { echo "FATAL: no cached workdir $BASE"; exit 1; }
n=0
for f in "$BASE"/*; do ln -s "$f" "$OUT/work/$(basename "$f")" && n=$((n+1)); done
echo "linked $n cached intermediates from $BASE"

{
  echo "=== A2 Step 3 verify (retention + copy-number gates): $S ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  echo "karyoscope      : $($KSBIN/karyoscope --version 2>&1)"
  echo "karyoscope-iscn : $($KSBIN/karyoscope-iscn --version 2>&1 | tail -1)"
  echo "iscn git        : $(git -C $ISCN rev-parse --short HEAD) on $(git -C $ISCN rev-parse --abbrev-ref HEAD)"
  echo "iscn dirty      : $(git -C $ISCN status --porcelain -- src | tr '\n' ' ')"
  echo "reused workdir  : $BASE ($n files, symlinked)"
} | tee "$OUT/provenance.txt"

echo "[$(date +%T)] running karyoscope-iscn analyze (cached intermediates)"
/usr/bin/time -v $KSBIN/karyoscope-iscn -v analyze \
    --sample "$S" \
    --input "hap1=$H1" \
    --input "hap2=$H2" \
    --output "$OUT/${S}.iscn.txt" \
    --workdir "$OUT/work" \
    --genome chm13 \
    --db-root "$ROOT" \
    --sex XX \
    --threads 8 \
    > "$OUT/analyze.log" 2> "$OUT/time.txt"
RC=$?
echo "analyze_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/time.txt" | awk '{print $NF}')" | tee -a "$OUT/provenance.txt"

# Nothing should have been recomputed: a real file in work/ means a cache miss.
REAL=$(find "$OUT/work" -maxdepth 1 -type f | wc -l)
echo "cache check: $REAL recomputed file(s) in the verify workdir (0 = full reuse)"

echo "=== ISCN karyotype: baseline -> verify ==="
echo "  baseline: $(grep -m1 -E '^[0-9]+,X' "$W/$S/${S}.iscn.txt" 2>/dev/null)"
echo "  verify  : $(grep -m1 -E '^[0-9]+,X' "$OUT/${S}.iscn.txt" 2>/dev/null)"
echo "=== warnings ==="
sed -n '/WARNING/,/^$/p' "$OUT/${S}.iscn.txt" 2>/dev/null | head -12
echo "=== intrachrom events (baseline -> verify) ==="
echo "--- baseline ---"; cat "$W/$S/${S}.iscn_intrachrom.tsv" 2>/dev/null
echo "--- verify ---";   cat "$OUT/${S}.iscn_intrachrom.tsv" 2>/dev/null
echo "=== regions diff (empty = breakpoints identical) ==="
diff "$W/$S/${S}.iscn_regions.tsv" "$OUT/${S}.iscn_regions.tsv" && echo "  IDENTICAL"
exit $RC
