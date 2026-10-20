#!/bin/bash
# A4 step 1: annotate CHM13 and HG002 hap1 on the HKS backend, KEEPING the
# pre-smoothing output (annotate writes *.presmoothed.bed.gz by default).
#
# Fresh runs are required, confirmed not assumed: the v2.3.0 benchmark job
# 35247555 deleted its intermediates (`rm -rf $BENCH/out` per iteration AND at
# the end of run_bench.sh), so nothing from it can be harvested.
#SBATCH --job-name=a4_annotate
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=a4_annotate_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing
RES=$W/results; mkdir -p "$RES"
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
KS=$KSBIN/karyoscope
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export KARYOSCOPE_GET_FEATUREIDS=${KS_SOFTWARE}/github/barthel/KaryoScope/native/get_featureIDs/build/get_featureIDs
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry

{ echo "=== A4 annotate ==="; echo "date: $(date -Is)"; echo "job: ${SLURM_JOB_ID:-NA} node: $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)"
  echo "karyoscope: $($KS --version)"; echo "threads: 16"
} | tee "$RES/provenance_annotate.txt"

run () { # $1 label $2 fasta $3 db
  local out=$W/annot_$1; rm -rf "$out"; mkdir -p "$out"
  echo "[$(date +%T)] annotate $1 ($3)"
  /usr/bin/time -v $KS -v annotate -i "$2" -o "$out" --db "$3" --db-root "$ROOT" \
      -t 16 --no-bgzip > "$RES/annotate_$1.log" 2> "$RES/time_$1.txt"
  local rc=$?
  echo "  $1 rc=$rc wall=$(grep 'Elapsed (wall' "$RES/time_$1.txt" | awk '{print $NF}')" \
      | tee -a "$RES/provenance_annotate.txt"
  [ $rc -eq 0 ] || return 1
}

run chm13      ${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta          HKS_human_CHM13_v2
run hg002hap1  ${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz    HKS_human_CHM13_v2
# one KMC run on a single feature set, to verify the fractions are
# backend-independent as A4.2 requires before reporting HKS throughout
echo "[$(date +%T)] KMC cross-check (repeat set only)"
out=$W/annot_chm13_kmc_repeat; rm -rf "$out"; mkdir -p "$out"
/usr/bin/time -v $KS -v annotate -i ${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta \
    -o "$out" --db KS_human_CHM13_v2 --db-root "$ROOT" --feature-set repeat -t 16 --no-bgzip \
    > "$RES/annotate_chm13_kmc.log" 2> "$RES/time_chm13_kmc.txt"
echo "  kmc rc=$? wall=$(grep 'Elapsed (wall' "$RES/time_chm13_kmc.txt" | awk '{print $NF}')" \
    | tee -a "$RES/provenance_annotate.txt"

echo "=== outputs ==="; ls -lh $W/annot_*/ | head -40
echo "=== DONE ==="
