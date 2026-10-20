#!/bin/bash
#SBATCH --job-name=tib_smoke
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --output=tib_smoke_%j.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
B=${KS_WORK}/a1_bench
W=$B/tib_smoke; rm -rf "$W"; mkdir -p "$W"
PY=${KS_WORK}/envs/tiberius/bin/python
export APPTAINER_CACHEDIR=${KS_WORK}/.apptainer_cache
export APPTAINER_TMPDIR=${KS_WORK}/.apptainer_tmp
# THE FIX: apptainer mounts $HOME, so the host's ~/.local numpy 2.4.3 shadowed
# the container's numpy and broke ml_dtypes. Block user site-packages inside.
export APPTAINERENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_PYTHONNOUSERSITE=1
nvidia-smi --query-gpu=name --format=csv,noheader
cp $B/smoke2/test.fa "$W/tiny.fa"
cd $B/Tiberius
$PY tiberius.py --singularity --genome "$W/tiny.fa" \
    --model_cfg mammalia_nosofttmasking_v2 --out "$W/tiny.gtf" --threads 8 \
    > "$W/tib.log" 2>&1
echo "tiberius_rc=$?"
echo "--- import errors? ---"; grep -ciE "ImportError|numpy.core" "$W/tib.log"
echo "--- gtf lines ---"; wc -l < "$W/tiny.gtf" 2>/dev/null || echo "NO GTF"
tail -6 "$W/tib.log"
