#!/bin/bash
# A5: SNP/indel divergence series on CHM13. One array task per mutated genome:
# mutate -> annotate (HKS, keeping pre-smoothing) -> concordance vs PROJECTED truth.
#
# Ordering matters and is easy to get backwards: truth must be PROJECTED into
# mutated coordinates FIRST, then put into k-mer space against the MUTATED .fai
# with N-regions re-derived from the MUTATED FASTA. N coordinates shift with
# indels, so projecting a pre-filled truth would put the k-1 trim and the
# exclude regions in the wrong place.
#
# Run matrix (12 tasks), from the brief:
#   1-5   substitution-only  0.001 0.005 0.01 0.02 0.05     seed 42
#   6-9   indel-only         0.0001 0.0005 0.001 0.005      seed 42
#   10-12 variance check     sub 0.01                       seeds 42 43 44
#SBATCH --job-name=a5_div
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-12%4
#SBATCH --output=a5_%A_%a.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a5_divergence
ARCH=${KS_SOFTWARE}/github/barthel/KaryoScope-archive
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
KS=$KSBIN/karyoscope
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
DB=HKS_human_CHM13_v2
HIER=$ROOT/$DB/hierarchy.tsv
REF=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta
# truth = the BEDs the database was built from
TRUTH_REPEAT=$ARCH/resources/KS_human_CHM13/CHM13_rep.bed
TRUTH_REGION=$ARCH/resources/KS_human_CHM13/CHM13_reg.bed
MERGE_REPEAT="Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat"

# ---- task -> parameters
case "${SLURM_ARRAY_TASK_ID:-1}" in
   1) SERIES=sub;   SUB=0.001;  IND=0;      SEED=42;;
   2) SERIES=sub;   SUB=0.005;  IND=0;      SEED=42;;
   3) SERIES=sub;   SUB=0.01;   IND=0;      SEED=42;;
   4) SERIES=sub;   SUB=0.02;   IND=0;      SEED=42;;
   5) SERIES=sub;   SUB=0.05;   IND=0;      SEED=42;;
   6) SERIES=indel; SUB=0;      IND=0.0001; SEED=42;;
   7) SERIES=indel; SUB=0;      IND=0.0005; SEED=42;;
   8) SERIES=indel; SUB=0;      IND=0.001;  SEED=42;;
   9) SERIES=indel; SUB=0;      IND=0.005;  SEED=42;;
  10) SERIES=var;   SUB=0.01;   IND=0;      SEED=42;;
  11) SERIES=var;   SUB=0.01;   IND=0;      SEED=43;;
  12) SERIES=var;   SUB=0.01;   IND=0;      SEED=44;;
  # 13: ZERO-DIVERGENCE BASELINE. Not in the brief, but required: the
  # region set carries a ~6.5% 'arm_as_centromeric' incompatible floor that
  # is a property of the priority-built database (shared k-mers between the
  # arms and the centromeric transition), NOT of divergence -- it is flat
  # across smoothing and cannot be caused by 3% k-mer loss. Without this row
  # every curve is read against the wrong origin.
  13) SERIES=baseline; SUB=0;   IND=0;      SEED=42;;
  *) echo "bad task id"; exit 1;;
esac
TAG="${SERIES}_sub${SUB}_ind${IND}_seed${SEED}"
D=$W/$TAG; RES=$D/results
rm -rf "$D"; mkdir -p "$RES"
PROV=$RES/provenance.txt
{
  echo "=== A5 $TAG ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]  node: $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)"
  echo "karyoscope: $($KS --version)   hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "series=$SERIES sub=$SUB indel=$IND seed=$SEED db=$DB"
  echo "expected surviving 31-mer fraction (1-sub)^31 = $(python3 -c "print(f'{(1-$SUB)**31:.4f}')")"
} | tee "$PROV"

# ---- 1. mutate (+ project truth into mutated coordinates)
echo "[$(date +%T)] mutate"
/usr/bin/time -v python3 ${KS_WORK}/mutate_fasta.py \
    -i "$REF" -o "$D/mut" --sub-rate "$SUB" --indel-rate "$IND" --seed "$SEED" \
    --bed repeat="$TRUTH_REPEAT" --bed region="$TRUTH_REGION" \
    > "$RES/mutate.log" 2> "$RES/time_mutate.txt"
echo "mutate_rc: $?" | tee -a "$PROV"
cp "$D/mut.summary.tsv" "$RES/" 2>/dev/null
awk 'NR>1{s+=$3;i+=$4;d+=$5;ib+=$6;db+=$7}END{printf "  totals: subs=%d ins=%d del=%d ins_bp=%d del_bp=%d\n",s,i,d,ib,db}' \
    "$D/mut.summary.tsv" | tee -a "$PROV"

# ---- 2. mutated-genome .fai and N-regions (MUST be re-derived, not projected)
echo "[$(date +%T)] faidx + N regions on the MUTATED genome"
samtools faidx "$D/mut.fa"
python3 "$ARCH/workflow/scripts/find_n_regions.py" "$D/mut.fa" > "$D/mut.n_regions.bed" 2>>"$RES/prep.log"
echo "  n_gap_intervals: $(wc -l < "$D/mut.n_regions.bed")" | tee -a "$PROV"

# ---- 3. annotate (keeps *.presmoothed.bed)
echo "[$(date +%T)] annotate"
/usr/bin/time -v $KS -v annotate -i "$D/mut.fa" -o "$D/annot" --db $DB --db-root "$ROOT" \
    -t 16 --no-bgzip > "$RES/annotate.log" 2> "$RES/time_annotate.txt"
echo "annotate_rc: $?" | tee -a "$PROV"
grep -E "Elapsed \(wall|Maximum resident" "$RES/time_annotate.txt" | tee -a "$PROV"

# ---- 4. concordance vs projected truth, pre and post smoothing
for SET in repeat region; do
  PRE=$(ls "$D/annot"/*.${SET}.presmoothed.bed* 2>/dev/null | head -1)
  POST=$(ls "$D/annot"/*.${SET}.smoothed.bed* 2>/dev/null | head -1)
  PROJ="$D/mut.${SET}.projected.bed"
  if [ -z "$PRE" ] || [ -z "$POST" ] || [ ! -s "$PROJ" ]; then
    echo "  $SET: MISSING inputs (pre='$PRE' post='$POST' proj='$PROJ')" | tee -a "$PROV"; continue
  fi
  ARGS=(--truth "$PROJ" --pre "$PRE" --post "$POST" --hierarchy "$HIER" \
        --feature-set "$SET" --fai "$D/mut.fa.fai" --gaps "$D/mut.n_regions.bed" \
        --kmer-size 31 --out "$RES/${SET}")
  [ "$SET" = repeat ] && ARGS+=(--merge "$MERGE_REPEAT")
  echo "[$(date +%T)] concordance: $SET"
  bash ${KS_WORK}/ks_concord.sh "${ARGS[@]}" >> "$RES/concord_${SET}.log" 2>&1
  echo "  $SET rc=$?" | tee -a "$PROV"
  [ -s "$RES/${SET}.summary.tsv" ] && { echo "--- $SET ---"; cat "$RES/${SET}.summary.tsv"; } | tee -a "$PROV"
done

# ---- 5. tidy row for the combined table
python3 - "$RES" "$SERIES" "$SUB" "$IND" "$SEED" <<'PY' | tee -a "$PROV"
import sys, os, csv
RES, series, sub, ind, seed = sys.argv[1:6]
out=os.path.join(RES,"divergence_rows.tsv")
with open(out,"w") as o:
    o.write("series\tsub_rate\tindel_rate\tseed\tfeature_set\tstage\tclass\tbases\tfraction_pct\n")
    for fs in ("repeat","region"):
        p=os.path.join(RES,f"{fs}.summary.tsv")
        if not os.path.exists(p): continue
        with open(p) as fh:
            next(fh)
            for line in fh:
                stage,c,b,frac=line.rstrip("\n").split("\t")
                o.write(f"{series}\t{sub}\t{ind}\t{seed}\t{fs}\t{stage}\t{c}\t{b}\t{frac}\n")
print(open(out).read() if os.path.exists(out) else "no rows")
PY

# Free bulk intermediates ONLY if concordance actually produced results —
# otherwise keep them, so a wrapper fix can be re-run without re-annotating
# (annotate is ~9 min; regenerating all 12 would be hours).
if [ -s "$RES/repeat.summary.tsv" ] && [ -s "$RES/region.summary.tsv" ]; then
  echo "concordance OK -> freeing intermediates" | tee -a "$PROV"
  rm -f "$D"/annot/*.presmoothed.bed "$D"/annot/*.smoothed.bed "$D/mut.fa" 2>/dev/null
else
  echo "concordance INCOMPLETE -> keeping intermediates for re-run" | tee -a "$PROV"
fi
echo "date_end: $(date -Is)" >> "$PROV"
echo "=== DONE $TAG ==="
