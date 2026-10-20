#!/bin/bash
# A1c: Tiberius (gene prediction) on HG002 v1.1 hap1, GPU.
#
# Deviation from the brief's example command: it shows
# --model_cfg mammalia_softmasking_v2, but that config sets softmasking:True
# and our input is a RAW UNMASKED assembly. Running a softmasking model on
# unmasked sequence silently degrades predictions. The brief's own text says to
# prefer the FASTA-only path, so this uses mammalia_nosofttmasking_v2
# (softmasking:False), and records the choice.
#SBATCH --job-name=a1c_tiberius
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=a1c_tiberius_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

B=${KS_WORK}/a1_bench
W=$B/a1c_tiberius; RES=$W/results
rm -rf "$W"; mkdir -p "$RES"
TIB=$B/Tiberius
PY=${KS_WORK}/envs/tiberius/bin/python
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
MODEL=mammalia_nosofttmasking_v2

# Keep the image build off /tmp: a 2.5 GB SIF build OOM-killed on a login node
# earlier in this session with its temp on /tmp.
export APPTAINER_CACHEDIR=${KS_WORK}/.apptainer_cache
export APPTAINER_TMPDIR=${KS_WORK}/.apptainer_tmp
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
# apptainer mounts $HOME, so the host ~/.local numpy 2.4.3 shadowed the
# container numpy and broke ml_dtypes (verified: job 35377201 failed this way,
# job 35378573 passed with these set).
export APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_PYTHONNOUSERSITE=1

PROV=$RES/provenance.txt
{
  echo "=== A1c Tiberius run ==="
  echo "date_start : $(date -Is)"
  echo "slurm_job  : ${SLURM_JOB_ID:-NA}"
  echo "node       : $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "--- GPU ---"
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>&1 | head -5
  echo "cuda_visible: ${CUDA_VISIBLE_DEVICES:-NA}"
  echo "--- software ---"
  echo "tiberius_pkg : $($PY -c 'import importlib.metadata as m; print(m.version("tiberius"))' 2>/dev/null)"
  echo "repo_commit  : $(git -C $TIB rev-parse --short HEAD 2>/dev/null)"
  echo "model_cfg    : $MODEL (softmasking:$(grep -m1 softmasking $TIB/model_cfg/${MODEL}.yaml | awk '{print $2}'))"
  echo "input        : $IN"
  echo "input_md5    : $(md5sum $IN | cut -d' ' -f1)"
} | tee "$PROV"

echo "[$(date +%T)] decompressing input"
zcat "$IN" > "$W/hg002.hap1.fa"

cd "$TIB"
CMD=($PY tiberius.py --singularity --genome "$W/hg002.hap1.fa" \
     --model_cfg "$MODEL" --out "$RES/tiberius.gtf" --threads 8)
echo "command    : ${CMD[*]}" | tee -a "$PROV"

echo "[$(date +%T)] running Tiberius (first run also pulls the container image)"
/usr/bin/time -v "${CMD[@]}" > "$RES/tiberius.log" 2> "$RES/time.tiberius.txt"
RC=$?
echo "tiberius_rc: $RC" | tee -a "$PROV"

{
  echo "=== headline numbers ==="
  grep -E "Elapsed \(wall clock\)|Maximum resident set size|Percent of CPU" "$RES/time.tiberius.txt"
  echo "=== output sanity ==="
  echo "gtf_lines : $(wc -l < "$RES/tiberius.gtf" 2>/dev/null || echo 0)"
  echo "genes     : $(awk '$3=="gene"' "$RES/tiberius.gtf" 2>/dev/null | wc -l)"
  echo "transcripts: $(awk '$3=="transcript"' "$RES/tiberius.gtf" 2>/dev/null | wc -l)"
  echo "(expect roughly 20k genes for a human haplotype; wildly different means"
  echo " the model or the masking assumption is wrong)"
} | tee -a "$PROV"

echo "date_end   : $(date -Is)" >> "$PROV"
# make the job state reflect the tool, not the wrapper
echo "=== DONE ==="; cat "$PROV"
[ "${RC:-1}" -eq 0 ] && [ -s "$RES/tiberius.gtf" ] || { echo "FAILED: rc=$RC or empty GTF"; exit 1; }
