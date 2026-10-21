#!/bin/bash
# =============================================================================
# A2 Step 3 -- THE GATE
# =============================================================================
# Question: does the chromosome pivot, at k=31, reproduce the PUBLISHED
# Robertsonian breakpoints before we spend k=21/51 compute?
#
# Published target (manuscript ROB section, machine-readable in
# KaryoScope-ROBs/conf/robs.yaml as the *_robertsonian SST1 overrides):
#   "In all three cases the refined density center fell within approximately
#    200--525 kb of the SST1 block midpoint"
#   99% CI widths: a2_step3 published_ci_width_mb in the ROB config
#
# SST1 midpoints, INDEPENDENTLY re-derived from our own acrocentric
# annotations (block sizes reproduce the published ones exactly; midpoints
# agree with robs.yaml to <= 27 bp): a2_step3 sst1_mid in the ROB config, on
# the fusion scaffold of each sample's *_robertsonian breakpoint group.
#
# Baseline for contrast -- the CURRENT DEFAULT (cytoband pivot) at k=31: two
# samples within the published range, the third (the flipped fusion contig)
# several megabases away, so the gate is a real test, not a formality.
#
# REUSES THE CACHED INTERMEDIATES, same mechanism as run_step1_verify.sh: the
# baseline workdirs are symlinked in, every pipeline step is skip-if-present,
# and anything actually recomputed lands in the gate workdir where the cache
# check will count it. Expected recomputation: ONE remap per haplotype, for
# `<stem>.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz`, which the
# cytoband-pivot baseline never needed.
#
# Per-sample values (fusion scaffold, SST1 midpoint) come from the ROB
# configuration in $ROBS_CONFIG (robs/conf/config.template.yaml, a2_step3).
# =============================================================================
#SBATCH --job-name=a2s3_gate
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=1-3
#SBATCH --output=a2s3_gate_%A_%a.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_GUARRACINO:?set KS_GUARRACINO (see the repository README, Paths)}"
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
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
OUT=$W/gate_$S; rm -rf "$OUT"; mkdir -p "$OUT/work"
H1=$SRC/$S/assembly/haplotype1/$S.haplotype1.fasta.gz
H2=$SRC/$S/assembly/haplotype2/$S.haplotype2.fasta.gz

BASE=$W/$S/work
[ -d "$BASE" ] || { echo "FATAL: no cached workdir $BASE"; exit 1; }
n=0
for f in "$BASE"/*; do ln -s "$f" "$OUT/work/$(basename "$f")" && n=$((n+1)); done
echo "linked $n cached intermediates from $BASE"

{
  echo "=== A2 Step 3 GATE (chromosome pivot, k=31): $S ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  echo "karyoscope      : $($KSBIN/karyoscope --version 2>&1)"
  echo "karyoscope-iscn : $($KSBIN/karyoscope-iscn --version 2>&1 | tail -1)"
  echo "iscn git        : $(git -C $ISCN rev-parse --short HEAD) on $(git -C $ISCN rev-parse --abbrev-ref HEAD)"
  echo "iscn dirty      : $(git -C $ISCN status --porcelain -- src | tr '\n' ' ')"
  echo "reused workdir  : $BASE ($n files, symlinked)"
  echo "refine pivot    : chromosome"
} | tee "$OUT/provenance.txt"

echo "[$(date +%T)] running karyoscope-iscn analyze --refine-pivot chromosome"
/usr/bin/time -v $KSBIN/karyoscope-iscn -v analyze \
    --sample "$S" \
    --input "hap1=$H1" \
    --input "hap2=$H2" \
    --output "$OUT/${S}.iscn.txt" \
    --workdir "$OUT/work" \
    --genome chm13 \
    --db-root "$ROOT" \
    --sex XX \
    --refine-pivot chromosome \
    --threads 8 \
    > "$OUT/analyze.log" 2> "$OUT/time.txt"
RC=$?
echo "analyze_rc=$RC wall=$(grep 'Elapsed (wall' "$OUT/time.txt" | awk '{print $NF}')" | tee -a "$OUT/provenance.txt"

# A real file in work/ is a cache miss. Expect exactly the chromosome remaps.
echo "=== cache check (expect only *.chromosome.smoothed.scaffolded.bed.gz) ==="
find "$OUT/work" -maxdepth 1 -type f -printf '  %f\n' | sort

echo "=== karyotype: cytoband pivot -> chromosome pivot ==="
echo "  cytoband  : $(grep -m1 -E '^[0-9]+,X' "$W/verify_$S/${S}.iscn.txt" 2>/dev/null)"
echo "  chromosome: $(grep -m1 -E '^[0-9]+,X' "$OUT/${S}.iscn.txt" 2>/dev/null)"

echo "=== regions table (chromosome pivot) ==="
cat "$OUT/${S}.iscn_regions.tsv"

echo "=== GATE: inferred_breakpoint vs SST1 midpoint ==="
# fusion scaffold (from the sample's breakpoint group) and SST1 midpoint, from $ROBS_CONFIG
read -r SCAF MID < <($KSBIN/python -c '
import os, sys, yaml
c = yaml.safe_load(open(os.environ["ROBS_CONFIG"]))
d = c["a2_step3"]["samples"][sys.argv[1]]
print(c["breakpoint_groups"][d["breakpoint_group"]][0]["scaffold"], d["sst1_mid"])' "$S")
awk -F'\t' -v scaf="$SCAF" -v mid="$MID" -v s="$S" '
  NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  $h["scaffold"]==scaf && $h["region_type"]=="breakpoint" {
    bp=$h["inferred_breakpoint"]+0
    if(bp==0){printf "  %s  NO inferred_breakpoint emitted -- GATE FAIL\n", s; next}
    d=bp-mid; ad=(d<0?-d:d)
    ci=($h["plateau_end"]+0)-($h["plateau_start"]+0)
    printf "  %s  bp=%d  sst1_mid=%d  delta=%+.1f kb  |delta|=%.1f kb  plateau=%.2f Mb  -> %s\n", \
      s, bp, mid, d/1000, ad/1000, ci/1e6, (ad<=525000 ? "PASS (within published 200-525 kb)" : "FAIL")
  }' "$OUT/${S}.iscn_regions.tsv"
exit $RC
