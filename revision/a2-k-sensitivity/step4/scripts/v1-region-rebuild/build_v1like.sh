#!/bin/bash
#SBATCH --job-name=v1like_build
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=build_k%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
K="${K:?}"; W=${KS_WORK}/ks_build_sources/v1like_ksweep
KS=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/karyoscope
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
echo "=== v1-like region build k=$K job ${SLURM_JOB_ID} node $(hostname) mem=${SLURM_MEM_PER_NODE}M ==="; $KS --version
/usr/bin/time -v $KS -v build --spec $W/build.k$K.yaml --db-root $W/db --force 2> $W/metrics/k$K.time.txt
RC=$?; grep -E "Elapsed|Maximum resident" $W/metrics/k$K.time.txt; echo "rc=$RC"
ls $W/db/HKS_human_CHM13_v1like_region_k$K/ 2>/dev/null; cat $W/db/HKS_human_CHM13_v1like_region_k$K/hierarchy.tsv 2>/dev/null | head -30
exit $RC
