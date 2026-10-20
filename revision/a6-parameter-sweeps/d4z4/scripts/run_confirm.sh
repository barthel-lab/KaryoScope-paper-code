#!/bin/bash
# A6 Sweep 4 refinement -- resolve the MIN_ARRAY_ELEMENT_SIZE plateau to the base pair.
# Step-50 scan over [1000,1750] plus a data-breakpoint walk to the exact edges.
#SBATCH --job-name=a6_d4z4_confirm
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=confirm_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
echo "=== A6 D4Z4 MIN_ARRAY_ELEMENT_SIZE plateau refinement ==="
echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}  node: $(hostname)"
echo "code: HPRC-D4Z4:analyze_d4z4_arrays.py (base_dir line env-overridable only)"
$PY ${KS_WORK}/a6_sweeps/d4z4_bisect/confirm_array_edge.py
RC=$?
echo "rc=$RC  done $(date -Is)"
exit $RC
