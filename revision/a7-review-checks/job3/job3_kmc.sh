#!/bin/bash
# Job 3 (second backend) -- same inputs through the KMC database KS_human_CHM13_v2, whose smoothing
# runs in karyoscope/core/smooth.py (the "Advance i" scan the brief reasons about).
# Run after job3_hks.sh, which extracts the FASTAs. The Robertsonian contig is included when
# ROBS_CONFIG is set (a7_job3 section), as in job3_hks.sh.
#SBATCH --job-name=a7_job3_kmc
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=job3_kmc_%j.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a7_review_checks; J=$W/job3
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
REG=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
STEMS="chr21_hap1 chr20_hap1"
if [ -n "${ROBS_CONFIG:-}" ]; then
  STEMS="$STEMS $($PY -c 'import os, yaml; print(yaml.safe_load(open(os.environ["ROBS_CONFIG"]))["a7_job3"]["stem"])')"
fi
PROV=$J/provenance_kmc.txt
mkdir -p $J/annot_kmc
{ echo "=== Job 3 (KMC backend, python smoother) ==="; echo "job=${SLURM_JOB_ID:-NA} node=$(hostname) start=$(date -Is)"
  echo "karyoscope: $(karyoscope --version 2>&1)"; } | tee $PROV
/usr/bin/time -v karyoscope annotate $(for s in $STEMS; do echo -i $J/fa/$s.fwd.fa -i $J/fa/$s.rc.fa; done) \
   -o $J/annot_kmc --db KS_human_CHM13_v2 --db-root $REG -t 16 --no-bgzip > $J/annotate_kmc.log 2> $J/annotate_kmc.time
echo "  rc=$? $(grep Elapsed $J/annotate_kmc.time | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $J/annotate_kmc.time | awk '{print $NF}')" | tee -a $PROV
ls $J/annot_kmc | wc -l | tee -a $PROV
$PY $W/scripts/job3_compare.py --annot-dir $J/annot_kmc --db KS_human_CHM13_v2 --backend kmc_v2.3.0 --out-prefix $W/results/rc_symmetry_kmc | tee -a $PROV
echo "end=$(date -Is)" | tee -a $PROV
