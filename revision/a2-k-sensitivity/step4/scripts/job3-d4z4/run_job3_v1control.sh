#!/bin/bash
# D4Z4 control with the ACTUAL v1 database (KS_human_CHM13, KMC, k=31 only):
# annotate the same 582 windows, rebuild the input representation, run the caller.
# Must reproduce the published counts (582 / 629 / 14433 / 28580) closely.
#SBATCH --job-name=a2s4_d4z4v1
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=d4z4_v1ctl_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job3_d4z4; C=${KS_WORK}/a2_step4/common
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
DBROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry; DB=KS_human_CHM13
CODE=${KS_WORK}/a6_sweeps/d4z4_bisect/code
PROV=$W/provenance_v1control.txt
{ echo "=== D4Z4 v1 control ($DB, KMC, k=31) ==="; echo "date: $(date -Is) job: ${SLURM_JOB_ID:-NA} node: $(hostname)"; echo "karyoscope: $(karyoscope --version)"; } | tee $PROV
A=$W/annot_v1_k31; rm -rf $A; mkdir -p $A
/usr/bin/time -v karyoscope annotate -i $W/fasta/windows_chr4.fa -i $W/fasta/windows_chr10.fa -o $A --db $DB --db-root $DBROOT \
    --feature-set region --feature-set subtelomeric --feature-set gene --feature-set repeat -t 16 --no-bgzip --no-smooth \
    > $W/logs/annotate_v1.log 2> $W/logs/time_annotate_v1.txt
echo "annotate rc=$? $(grep Elapsed $W/logs/time_annotate_v1.txt | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $W/logs/time_annotate_v1.txt | awk '{print $NF}')" | tee -a $PROV
ls $A | head -12 | tee -a $PROV
B=$W/base_v1_k31; rm -rf $B; mkdir -p $B/data $B/results
python3 $W/d4z4_postprocess.py --annot-dir $A --db $DB --windows $C/d4z4_windows.tsv --bin-script $C/bin_features_feb2026.py --outdir $B/data | tee -a $PROV
( cd $CODE && D4Z4_BASE_DIR=$B python analyze_d4z4_arrays.py > $B/results/stdout.log 2>&1 ); echo "caller rc=$?" | tee -a $PROV
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{a+=$h["num_arrays"]; m+=$h["repeat_count_main"]; t+=$h["repeat_count_total"]; n++} END{printf "v1 control: haplotypes=%d arrays=%d repeat_main=%d repeat_total=%d (published 582/629/14433/28580)\n", n, a, m, t}' $B/results/d4z4_array_analysis.tsv | tee -a $PROV
echo "date_end: $(date -Is)" | tee -a $PROV
