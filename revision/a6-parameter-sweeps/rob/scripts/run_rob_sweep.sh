#!/bin/bash
# A6 Sweep 2 (ROB) -- window / coverage / dominance. Interpretation layer only:
# reads the existing k=31 scaffolded chromosome BEDs, so no annotate step.
# sweep_rob.py is run from $KS_WORK/a6_sweeps/rob (copy it there).
#SBATCH --job-name=a6_rob
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=a6rob_%j.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_ISCN_SRC:?set KS_ISCN_SRC to vendored/KaryoScope-ISCN/src in this repository}"
export KS_ISCN_SRC
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
set -uo pipefail
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
ISCN=${KS_SOFTWARE}/github/barthel/KaryoScope-ISCN
echo "=== A6 ROB parameter sweep ==="
echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}  node: $(hostname)"
echo "iscn: $(git -C $ISCN rev-parse --short HEAD) on $(git -C $ISCN rev-parse --abbrev-ref HEAD)"
$KSBIN/python ${KS_WORK}/a6_sweeps/rob/sweep_rob.py
RC=$?
echo "rc=$RC  done $(date -Is)"
exit $RC
