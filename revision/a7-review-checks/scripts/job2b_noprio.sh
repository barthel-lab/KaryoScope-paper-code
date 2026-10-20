#!/bin/bash
# Job 2b -- build the gene set WITHOUT the priority file, annotate HG002 hap1, tally labels.
#SBATCH --job-name=a7_job2b
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=job2b_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a7_review_checks; J=$W/job2_noprio
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
PROV=$J/provenance.txt
{ echo "=== Job 2b: gene set without priorities ==="; echo "job=${SLURM_JOB_ID:-NA} node=$(hostname) start=$(date -Is)"
  echo "karyoscope: $(karyoscope --version 2>&1)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "input bed sha256: $(sha256sum ${KS_WORK}/recipe_repro/chm13_gene.bed | cut -c1-16)..."; } | tee $PROV
mkdir -p $J/db $J/extmem $J/annot
echo "[$(date +%T)] build" | tee -a $PROV
/usr/bin/time -v karyoscope build --spec $J/build.yaml --db-root $J/db --force > $J/build.log 2> $J/build.time
echo "  build rc=$? $(grep Elapsed $J/build.time | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $J/build.time | awk '{print $NF}')" | tee -a $PROV
echo "[$(date +%T)] annotate HG002 hap1, gene set" | tee -a $PROV
/usr/bin/time -v karyoscope annotate -i ${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz -o $J/annot \
   --db HKS_human_CHM13_v2_gene_noprio --db-root $J/db --feature-set gene -t 16 --no-bgzip > $J/annotate.log 2> $J/annotate.time
echo "  annotate rc=$? $(grep Elapsed $J/annotate.time | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $J/annotate.time | awk '{print $NF}')" | tee -a $PROV
ls $J/annot | tee -a $PROV
$PY $W/scripts/job2_gene_labels.py --pre $J/annot/hg002v1.1_hap1.HKS_human_CHM13_v2_gene_noprio.gene.presmoothed.bed \
   --post $J/annot/hg002v1.1_hap1.HKS_human_CHM13_v2_gene_noprio.gene.smoothed.bed \
   --nbed ${KS_WORK}/a4_smoothing/HG002.hap1.n_regions.renamed.bed \
   --db $J/db/HKS_human_CHM13_v2_gene_noprio/hierarchy.tsv --assembly hg002hap1 --out $W/results/hg002_gene_labels_by_stage_NOPRIORITY.tsv | tee -a $PROV
echo "end=$(date -Is)" | tee -a $PROV
