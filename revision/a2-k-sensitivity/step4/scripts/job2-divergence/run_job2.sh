#!/bin/bash
# A2 Step 4 / Job 2 -- the A5 divergence series at k = 21, 31, 51.
# One array task per mutation condition. Each task:
#   1. regenerates the mutated CHM13 with mutate_fasta.py (deterministic, seed 42)
#      and VERIFIES it against the .fai and mutation summary A5 left on disk
#      (the FASTAs themselves were deleted after A5's concordance succeeded);
#   2. annotates it at k = 21, 31, 51 with the six-k database family
#      (HKS_human_CHM13_v2_k<K>; k = 31 is the CONTROL and must reproduce the
#      A5 rows, which used the shipped HKS_human_CHM13_v2);
#   3. runs the published concordance chain (ks_concord.sh, --kmer-size K) on the
#      repeat set (the brief's set) and the region set (cheap, for context)
#      against the PROJECTED truth, pre- and post-smoothing.
# Conditions (brief): sub 0.001, 0.01, 0.05; indel 0.005; all seed 42.
#SBATCH --job-name=a2s4_div
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-4%4
#SBATCH --output=div_%A_%a.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job2_divergence
A5=${KS_WORK}/a5_divergence
ARCH=${KS_SOFTWARE}/github/barthel/KaryoScope-archive
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
KS=$KSBIN/karyoscope
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
DBROOT=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db
HIER=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv
REF=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta
TRUTH_REPEAT=$ARCH/resources/KS_human_CHM13/CHM13_rep.bed
TRUTH_REGION=$ARCH/resources/KS_human_CHM13/CHM13_reg.bed
MERGE_REPEAT="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) SERIES=sub;   SUB=0.001; IND=0;     SEED=42;;
  2) SERIES=sub;   SUB=0.01;  IND=0;     SEED=42;;
  3) SERIES=sub;   SUB=0.05;  IND=0;     SEED=42;;
  4) SERIES=indel; SUB=0;     IND=0.005; SEED=42;;
  *) echo "bad task"; exit 1;;
esac
TAG="${SERIES}_sub${SUB}_ind${IND}_seed${SEED}"
D=$W/$TAG; RES=$D/results; mkdir -p "$RES"
PROV=$RES/provenance.txt
{
  echo "=== A2 step4 job2 $TAG ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  lscpu | grep -E 'Model name|Core\(s\) per socket|Socket'
  echo "karyoscope: $($KS --version)   hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "db_root: $DBROOT  dbs: HKS_human_CHM13_v2_k{21,31,51}  threads: 16"
  echo "series=$SERIES sub=$SUB indel=$IND seed=$SEED"
} | tee "$PROV"

# ---- 1. regenerate the mutated genome and verify it is the A5 one
if [ ! -s "$D/mut.fa" ] || [ ! -s "$D/mut.fa.fai" ]; then
  echo "[$(date +%T)] mutate (seed $SEED)" | tee -a "$PROV"
  /usr/bin/time -v python3 ${KS_WORK}/mutate_fasta.py \
      -i "$REF" -o "$D/mut" --sub-rate "$SUB" --indel-rate "$IND" --seed "$SEED" \
      --bed repeat="$TRUTH_REPEAT" --bed region="$TRUTH_REGION" > "$RES/mutate.log" 2> "$RES/time_mutate.txt"
  echo "mutate_rc: $?" | tee -a "$PROV"
  samtools faidx "$D/mut.fa"
fi
if cmp -s "$D/mut.fa.fai" "$A5/$TAG/mut.fa.fai" && cmp -s "$D/mut.summary.tsv" "$A5/$TAG/mut.summary.tsv" \
   && cmp -s "$D/mut.repeat.projected.bed" "$A5/$TAG/mut.repeat.projected.bed"; then
  echo "VERIFIED: regenerated .fai, mutation summary and projected truth are byte-identical to A5's" | tee -a "$PROV"
else
  echo "FATAL: regenerated mutant differs from A5's on-disk .fai/summary/projected truth -- stopping" | tee -a "$PROV"
  cmp "$D/mut.fa.fai" "$A5/$TAG/mut.fa.fai" | head -2; exit 3
fi
python3 "$ARCH/workflow/scripts/find_n_regions.py" "$D/mut.fa" > "$D/mut.n_regions.bed" 2>>"$RES/prep.log"
echo "  n_gap_intervals: $(wc -l < "$D/mut.n_regions.bed")" | tee -a "$PROV"

# ---- 2+3. annotate at each k, then concordance
OK_ALL=1
for K in 21 31 51; do
  DB=HKS_human_CHM13_v2_k${K}; OUT=$D/annot_k${K}; rm -rf "$OUT"; mkdir -p "$OUT"
  echo "[$(date +%T)] annotate k=$K ($DB)" | tee -a "$PROV"
  /usr/bin/time -v $KS -v annotate -i "$D/mut.fa" -o "$OUT" --db $DB --db-root "$DBROOT" \
      --feature-set repeat --feature-set region -t 16 --no-bgzip \
      > "$RES/annotate_k${K}.log" 2> "$RES/time_annotate_k${K}.txt"
  RC=$?
  echo "  annotate_rc=$RC $(grep -E 'Elapsed' "$RES/time_annotate_k${K}.txt" | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' "$RES/time_annotate_k${K}.txt" | awk '{print $NF}')" | tee -a "$PROV"
  [ $RC -eq 0 ] || { OK_ALL=0; continue; }
  for SET in repeat region; do
    PRE=$OUT/mut.${DB}.${SET}.presmoothed.bed; POST=$OUT/mut.${DB}.${SET}.smoothed.bed
    PROJ="$D/mut.${SET}.projected.bed"
    [ -s "$PRE" ] && [ -s "$POST" ] || { echo "  MISSING $SET BEDs for k=$K" | tee -a "$PROV"; OK_ALL=0; continue; }
    ARGS=(--truth "$PROJ" --pre "$PRE" --post "$POST" --hierarchy "$HIER" --feature-set "$SET" \
          --fai "$D/mut.fa.fai" --gaps "$D/mut.n_regions.bed" --kmer-size "$K" --out "$RES/${SET}_k${K}")
    [ "$SET" = repeat ] && ARGS+=(--merge "$MERGE_REPEAT")
    echo "[$(date +%T)] concordance $SET k=$K" | tee -a "$PROV"
    bash ${KS_WORK}/ks_concord.sh "${ARGS[@]}" >> "$RES/concord_${SET}_k${K}.log" 2>&1
    echo "  $SET k=$K rc=$?" | tee -a "$PROV"
    [ -s "$RES/${SET}_k${K}.summary.tsv" ] && { echo "--- $SET k=$K ---"; cat "$RES/${SET}_k${K}.summary.tsv"; } | tee -a "$PROV"
    [ -s "$RES/${SET}_k${K}.summary.tsv" ] || OK_ALL=0
  done
  # presmoothed novel fraction straight from the BED as well (independent of the concordance chain)
  awk -v k=$K -v set=repeat '$4=="novel"{n+=$3-$2} {t+=$3-$2} END{printf "  raw presmoothed repeat BED k=%d: novel=%d of %d = %.6f\n", k, n, t, n/t}' \
      "$OUT/mut.${DB}.repeat.presmoothed.bed" | tee -a "$PROV"
done

# ---- tidy rows
python3 - "$RES" "$SERIES" "$SUB" "$IND" "$SEED" <<'PY' | tee -a "$PROV"
import sys, os
RES, series, sub, ind, seed = sys.argv[1:6]
out = os.path.join(RES, "divergence_rows.tsv")
with open(out, "w") as o:
    o.write("series\tsub_rate\tindel_rate\tseed\tk\tfeature_set\tstage\tclass\tbases\tfraction_pct\n")
    for k in (21, 31, 51):
        for fs in ("repeat", "region"):
            p = os.path.join(RES, f"{fs}_k{k}.summary.tsv")
            if not os.path.exists(p): continue
            with open(p) as fh:
                next(fh)
                for line in fh:
                    stage, c, b, frac = line.rstrip("\n").split("\t")
                    o.write(f"{series}\t{sub}\t{ind}\t{seed}\t{k}\t{fs}\t{stage}\t{c}\t{b}\t{frac}\n")
print(open(out).read())
PY
if [ $OK_ALL -eq 1 ]; then
  echo "all concordance OK -> freeing BEDs and mut.fa" | tee -a "$PROV"
  rm -f "$D"/annot_k*/*.bed "$D/mut.fa"
else
  echo "INCOMPLETE -> keeping intermediates" | tee -a "$PROV"
fi
echo "date_end: $(date -Is)" >> "$PROV"
