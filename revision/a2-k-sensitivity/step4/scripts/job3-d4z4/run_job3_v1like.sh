#!/bin/bash
# D4Z4 at k = 21, 31, 51 with the v1-like REGION rebuilds (HKS_human_CHM13_v1like_region_k<K>).
# Region set from the v1-like build; gene/subtelomeric/repeat from the v2 six-k annotation
# already on disk (those sets match the published inputs to <0.001%). Caller unchanged.
#SBATCH --job-name=a2s4_d4z4v1k
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=d4z4_v1like_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job3_d4z4; C=${KS_WORK}/a2_step4/common
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
DBROOT=${KS_WORK}/ks_build_sources/v1like_ksweep/db
CODE=${KS_WORK}/a6_sweeps/d4z4_bisect/code
PROV=$W/provenance_v1like.txt
{ echo "=== D4Z4 by k, v1-like region builds ==="; echo "date: $(date -Is) job: ${SLURM_JOB_ID:-NA} node: $(hostname)"; echo "karyoscope: $(karyoscope --version)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"; } | tee $PROV
for K in ${KLIST:-21 31 51}; do
  DB=HKS_human_CHM13_v1like_region_k$K; A=$W/annot_v1like_k$K; rm -rf $A; mkdir -p $A
  [ -d $DBROOT/$DB ] || { echo "MISSING DB $DB" | tee -a $PROV; continue; }
  echo "[$(date +%T)] annotate k=$K ($DB, region only)" | tee -a $PROV
  /usr/bin/time -v karyoscope annotate -i $W/fasta/windows_chr4.fa -i $W/fasta/windows_chr10.fa -o $A --db $DB --db-root $DBROOT \
      --feature-set region -t 16 --no-bgzip --no-smooth > $W/logs/annotate_v1like_k$K.log 2> $W/logs/time_annotate_v1like_k$K.txt
  echo "  rc=$? $(grep Elapsed $W/logs/time_annotate_v1like_k$K.txt | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $W/logs/time_annotate_v1like_k$K.txt | awk '{print $NF}')" | tee -a $PROV
  B=$W/base_v1like_k$K; rm -rf $B; mkdir -p $B/data $B/results
  python3 $W/d4z4_postprocess.py --annot-dir $W/annot_k$K --db HKS_human_CHM13_v2_k$K --region-annot-dir $A --region-db $DB \
      --windows $C/d4z4_windows.tsv --bin-script $C/bin_features_feb2026.py --outdir $B/data | tee -a $PROV
  ( cd $CODE && D4Z4_BASE_DIR=$B python analyze_d4z4_arrays.py > $B/results/stdout.log 2>&1 ); echo "  caller rc=$?" | tee -a $PROV
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{a+=$h["num_arrays"]; m+=$h["repeat_count_main"]; t+=$h["repeat_count_total"]; n++} END{printf "  v1-like k=%s: haplotypes=%d arrays=%d repeat_main=%d repeat_total=%d (published 582/629/14433/28580)\n", K, n, a, m, t}' K=$K $B/results/d4z4_array_analysis.tsv | tee -a $PROV
done
echo "date_end: $(date -Is)" | tee -a $PROV
