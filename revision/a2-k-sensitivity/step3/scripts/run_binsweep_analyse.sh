#!/bin/bash
# Runs binsweep_analyse.py (copied into $KS_WORK/a2_step3) over the binned k-sweep.
#SBATCH --job-name=a2s3_binan
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=a2s3_binan_%j.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_ISCN_SRC:?set KS_ISCN_SRC to vendored/KaryoScope-ISCN/src in this repository}"
export KS_ISCN_SRC
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
set -uo pipefail
cd ${KS_WORK}/a2_step3
/usr/bin/time -v ${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python \
    binsweep_analyse.py 2> >(tail -25 >&2)
echo "rc=$?"
