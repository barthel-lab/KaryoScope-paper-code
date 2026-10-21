#!/bin/bash
# Job 3 -- reverse-complement symmetry of smoothing. Extract 3 query contigs, RC them with seqtk,
# annotate both orientations (HKS_human_CHM13_v2, k=31, hks smooth max-gap 1000 = annotate's default),
# then compare per base after mapping the RC annotation back to forward k-mer coordinates.
# The two HG002 contigs are public. The third, the Robertsonian fusion contig, is read from the
# a7_job3 section of $ROBS_CONFIG (sample, haplotype, contig name, file stem); with ROBS_CONFIG
# unset, only the two HG002 contigs are run.
#SBATCH --job-name=a7_job3_hks
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=job3_hks_%j.log
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a7_review_checks; J=$W/job3
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
REG=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
HG002=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
STEMS="chr21_hap1 chr20_hap1"
if [ -n "${ROBS_CONFIG:-}" ]; then
  : "${KS_GUARRACINO:?set KS_GUARRACINO (see the repository README, Paths)}"
  read -r ROB_SAMPLE ROB_HAP ROB_CONTIG ROB_STEM < <($PY -c '
import os, yaml
d = yaml.safe_load(open(os.environ["ROBS_CONFIG"]))["a7_job3"]
print(d["sample"], d["haplotype"], d["contig"], d["stem"])')
  ROB=${KS_GUARRACINO}/ks-robs/results/$ROB_SAMPLE/assembly/$ROB_HAP/$ROB_SAMPLE.$ROB_HAP.fasta.gz
  STEMS="$STEMS $ROB_STEM"
  LAST=$ROB_STEM
else
  LAST=chr20_hap1
fi
PROV=$J/provenance_hks.txt
mkdir -p $J/fa $J/annot_hks
{ echo "=== Job 3 (HKS backend) ==="; echo "job=${SLURM_JOB_ID:-NA} node=$(hostname) start=$(date -Is)"
  echo "karyoscope: $(karyoscope --version 2>&1)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')  seqtk: $(seqtk 2>&1 | grep -m1 Version)"
  echo "smoother on this path: hks smooth (src/smooth.rs), invoked by karyoscope.core.io.hks.run_hks_smooth with max_gap default 1000"; } | tee $PROV
if [ ! -s $J/fa/$LAST.fwd.fa ]; then
  echo "[$(date +%T)] extracting" | tee -a $PROV
  samtools faidx $HG002 chr21_hap1 > $J/fa/chr21_hap1.fwd.fa
  samtools faidx $HG002 chr20_hap1 > $J/fa/chr20_hap1.fwd.fa
  if [ -n "${ROBS_CONFIG:-}" ]; then
    zcat $ROB | awk -v name=">$ROB_CONTIG" '/^>/{p=($1==name)} p' > $J/fa/$ROB_STEM.fwd.fa
  fi
  for s in $STEMS; do
    seqtk seq -r $J/fa/$s.fwd.fa > $J/fa/$s.rc.fa
    samtools faidx $J/fa/$s.fwd.fa; samtools faidx $J/fa/$s.rc.fa
    echo "  $s fwd: $(cut -f1,2 $J/fa/$s.fwd.fa.fai | tr '\t' ' ')   rc: $(cut -f1,2 $J/fa/$s.rc.fa.fai | tr '\t' ' ')" | tee -a $PROV
  done
  # RC self-check: seqtk seq -r applied twice must return the original sequence
  seqtk seq -r $J/fa/chr21_hap1.rc.fa | seqtk seq -l0 - | tail -n +2 | md5sum > $J/fa/rr.md5
  seqtk seq -l0 $J/fa/chr21_hap1.fwd.fa | tail -n +2 | md5sum > $J/fa/fwd.md5
  cmp -s $J/fa/rr.md5 $J/fa/fwd.md5 && echo "  rc(rc(chr21)) == chr21: OK" | tee -a $PROV || echo "  rc(rc(chr21)) != chr21: FAIL" | tee -a $PROV
fi
echo "[$(date +%T)] annotate HKS ($(( $(echo $STEMS | wc -w) * 2 )) inputs, all feature sets)" | tee -a $PROV
/usr/bin/time -v karyoscope annotate $(for s in $STEMS; do echo -i $J/fa/$s.fwd.fa -i $J/fa/$s.rc.fa; done) \
   -o $J/annot_hks --db HKS_human_CHM13_v2 --db-root $REG -t 16 --no-bgzip > $J/annotate_hks.log 2> $J/annotate_hks.time
echo "  rc=$? $(grep Elapsed $J/annotate_hks.time | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $J/annotate_hks.time | awk '{print $NF}')" | tee -a $PROV
ls $J/annot_hks | wc -l | tee -a $PROV
echo "[$(date +%T)] compare" | tee -a $PROV
$PY $W/scripts/job3_compare.py --annot-dir $J/annot_hks --db HKS_human_CHM13_v2 --backend hks_v2.3.0 --out-prefix $W/results/rc_symmetry_hks | tee -a $PROV
echo "end=$(date -Is)" | tee -a $PROV
