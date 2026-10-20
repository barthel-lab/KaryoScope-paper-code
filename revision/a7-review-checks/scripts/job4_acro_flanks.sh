#!/bin/bash
#SBATCH --job-name=a7_job4
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=job4_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
echo "job=${SLURM_JOB_ID:-NA} node=$(hostname) start=$(date -Is)"
${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python ${KS_WORK}/a7_review_checks/scripts/job4_acro_flanks.py
echo "end=$(date -Is)"
