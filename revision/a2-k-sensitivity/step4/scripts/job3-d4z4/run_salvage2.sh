#!/bin/bash
#SBATCH --job-name=a2s4_salv2
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=salvage_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
echo "date: $(date -Is) job: ${SLURM_JOB_ID} node: $(hostname)"
${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python \
  ${KS_WORK}/a2_step4/job3_d4z4/salvage/sweep_salvage2.py
echo "rc=$? done $(date -Is)"
