#!/bin/bash
# A2 Step 3, step 1: BASELINE run of the current ISCN pipeline on the three ROB
# assemblies at k=31, with the shipped databases and NO code changes.
#
# Purpose is not the k-sweep. It is to establish, before touching anything:
#   * which haplotype/scaffold each Robertsonian fusion actually lives on
#     (an open unknown - the legacy ROBs data dir is gitignored and long gone),
#   * the detected breakpoint coordinates from the maintained implementation,
#   * that the pipeline runs on these assemblies at all.
#
# Both haplotypes are passed per sample: analyze reports a single haplotype
# against a HAPLOID baseline, which would misread these XX samples.
#SBATCH --job-name=a2s3_baseline
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_baseline_%A_%a.log
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
OUT=$W/$S; rm -rf "$OUT"; mkdir -p "$OUT"
H1=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
H2=$SRC/$S/assembly/haplotype2/$S.haplotype2.fasta.gz
for f in "$H1" "$H2"; do [ -r "$f" ] || { echo "FATAL: unreadable $f"; exit 1; }; done

{
  echo "=== A2 Step 3 baseline: $S ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "karyoscope      : $($KSBIN/karyoscope --version 2>&1)"
  echo "karyoscope-iscn : $($KSBIN/karyoscope-iscn --version 2>&1 | tail -1)"
  echo "hap1: $H1"
  echo "hap2: $H2"
  echo "cytoband_db=HKS_human_CHM13_cytoband  region_db=HKS_human_CHM13_v2  (k=31, shipped)"
} | tee "$OUT/provenance.txt"

echo "[$(date +%T)] running karyoscope-iscn analyze"
/usr/bin/time -v $KSBIN/karyoscope-iscn -v analyze \
    --sample "$S" \
    --input "hap1=$H1" \
    --input "hap2=$H2" \
    --output "$OUT/${S}.iscn.txt" \
    --workdir "$OUT/work" \
    --genome chm13 \
    --db-root "$ROOT" \
    --sex XX \
    --threads 16 \
    > "$OUT/analyze.log" 2> "$OUT/time.txt"
RC=$?
echo "analyze_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/time.txt" | awk '{print $NF}')" | tee -a "$OUT/provenance.txt"

echo "=== outputs ==="; ls -lh "$OUT"/*.txt "$OUT"/*.tsv 2>/dev/null | awk '{print "  "$9,$5}'
echo "=== ISCN report ==="; sed -n '1,40p' "$OUT/${S}.iscn.txt" 2>/dev/null
echo "=== detected regions (breakpoint coordinates + which scaffold) ==="
column -t -s$'\t' "$OUT/${S}.iscn_regions.tsv" 2>/dev/null | head -25 || \
  ls "$OUT"/*regions* 2>/dev/null
exit $RC
