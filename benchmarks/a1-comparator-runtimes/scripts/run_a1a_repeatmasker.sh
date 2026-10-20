#!/bin/bash
# A1a: RepeatMasker on HG002 v1.1 hap1, for the head-to-head runtime comparison
# against KaryoScope (brief: $KS_WORK/a1-head-to-head-brief.md).
#
# Deviations from the brief, both forced by this cluster:
#   - --constraint instead of --nodelist: same CPU model is what makes the
#     comparison valid; pinning a 48h job to one node just queues it behind
#     whatever is resident (and g-c-1-4-20 is running our thread sweep).
#   - /scratch instead of /local/scratch: no local disk here (TmpDisk=0).
#     annbench3 also ran on /scratch, so storage is matched between the two
#     tools, which is what comparability actually requires.
#SBATCH --job-name=a1a_rm_hg002hap1
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=a1a_rm_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a1_bench/a1a_repeatmasker
RES=$W/results
rm -rf "$W"; mkdir -p "$W/out" "$RES"
RM=${KS_SOFTWARE}/RepeatMasker
A1=${KS_WORK}/envs/a1tools
# Private rmblast + a python that has h5py (famdb.py needs it), WITHOUT
# editing the shared RepeatMasker config.
export RMBLAST_DIR=$A1/bin
export PATH=$RM:$RM/util:$A1/bin:${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz

PROV=$RES/provenance.txt
{
  echo "=== A1a RepeatMasker run ==="
  echo "date_start   : $(date -Is)"
  echo "slurm_job    : ${SLURM_JOB_ID:-NA}"
  echo "node         : $(hostname)"
  lscpu | grep -E 'Model name|^CPU\(s\)|Socket|Core\(s\) per socket|Thread\(s\) per core'
  echo "cpus_alloc   : ${SLURM_CPUS_PER_TASK:-NA}  (-pa 8 = 8 jobs x 4 rmblast threads = 32 cores)"
  echo "mem_alloc    : ${SLURM_MEM_PER_NODE:-NA} MB"
  echo "input_path   : $IN"
  echo "input_md5    : $(md5sum $IN | cut -d' ' -f1)"
  echo "repeatmasker : $($RM/RepeatMasker -v 2>&1 | head -1)"
  echo "rmblastn     : $($A1/bin/rmblastn -version 2>&1 | head -1)"
  echo "trf          : ${KS_SOFTWARE}/trf409.linux64"
  echo "--- famdb library ---"
  $RM/famdb.py -i $RM/Libraries/famdb info 2>&1 | grep -viE "_warn|UserWarning"
  echo "--- NOTE ---"
  echo "The paper's KaryoScope repeat feature set derives from RepeatMasker"
  echo "4.1.2p1 / Dfam 2022Apr14 (UCSC CHM13 track). This run is 4.1.6 /"
  echo "Dfam 3.8 - fine for a RUNTIME comparison, but record it before any"
  echo "concordance analysis against this output."
} | tee "$PROV"

echo "[$(date +%T)] decompressing input (RepeatMasker cannot read .gz)"
zcat "$IN" > "$W/hg002.hap1.fa"
echo "decompressed_bases: $(grep -v '^>' "$W/hg002.hap1.fa" | tr -d '\n' | wc -c)" | tee -a "$PROV"

CMD=($RM/RepeatMasker -engine rmblast -species human -pa 8 -dir "$W/out" "$W/hg002.hap1.fa")
echo "command      : ${CMD[*]}" | tee -a "$PROV"

echo "[$(date +%T)] starting RepeatMasker"
/usr/bin/time -v "${CMD[@]}" 2> "$RES/time.repeatmasker.txt"
RC=$?
echo "[$(date +%T)] RepeatMasker rc=$RC"
echo "repeatmasker_rc: $RC" >> "$PROV"

{
  echo "=== headline numbers ==="
  grep -E "Elapsed \(wall clock\)|Maximum resident set size|Percent of CPU" "$RES/time.repeatmasker.txt"
  echo "=== sanity check: bases masked (expect 45-55%) ==="
  sed -n '1,10p' "$W/out"/*.tbl 2>/dev/null
} | tee -a "$PROV"

# Non-redundant per-base class map; raw .out double-counts nested repeats.
echo "[$(date +%T)] RM2Bed"
python3 $RM/util/RM2Bed.py --ovlp_resolution higher_score \
    -o "$RES" -d "$W/out"/*.out > "$RES/rm2bed.log" 2>&1 \
    && echo "rm2bed_rc: 0" >> "$PROV" \
    || echo "rm2bed_rc: FAILED (see rm2bed.log; can be redone offline)" >> "$PROV"

cp "$W/out"/*.tbl "$RES/" 2>/dev/null
gzip -c "$W/out"/*.out > "$RES/repeatmasker.out.gz" 2>/dev/null
echo "date_end     : $(date -Is)" >> "$PROV"
echo "=== DONE ==="; cat "$PROV"
