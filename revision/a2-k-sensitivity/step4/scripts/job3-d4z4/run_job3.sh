#!/bin/bash
# A2 Step 4 / Job 3 -- D4Z4 array reconstruction at k = 21, 31, 51.
# 1. Extract the SAME 500 kb q-terminal windows the published pipeline used
#    (last 500 kb of each of the 582 T2T chr4/chr10 contigs, forward frame --
#    the published scripts dropped every reverse-complemented contig).
# 2. Annotate them with the six-k family (region, subtelomeric, gene, repeat).
# 3. Rebuild the published input representation (100 bp bins, contig-frame
#    boundaries, normalised coordinates) and run Floris's caller UNCHANGED at
#    published constants. Also run it on the PUBLISHED inputs as the baseline.
#SBATCH --job-name=a2s4_d4z4
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=d4z4_%j.log
: "${KS_PROJECTS:?set KS_PROJECTS (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job3_d4z4; C=${KS_WORK}/a2_step4/common
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
DBROOT=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db
PUB=${KS_PROJECTS}/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed
CODE=${KS_WORK}/a6_sweeps/d4z4_bisect/code   # upstream HEAD e09c1d7 + base_dir env patch (diff-verified)
PROV=$W/provenance.txt
{
  echo "=== A2 step4 job3 D4Z4 by k ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}  node: $(hostname)"
  lscpu | grep -E 'Model name'
  echo "karyoscope: $(karyoscope --version)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')  threads: 16"
  echo "dbs: HKS_human_CHM13_v2_k{21,31,51} at $DBROOT"
  echo "caller: barthel-lab/HPRC-D4Z4 analyze_d4z4_arrays.py @ e09c1d7 (+1-line base_dir patch)"
  echo "binning: KaryoScope-archive bin_features.py as of 18dccc6 (Feb 2026), plain -b 100"
} | tee "$PROV"

# ---- 1. extract windows (once)
mkdir -p $W/fasta
for chrom in chr4 chr10; do
  FA=$W/fasta/windows_$chrom.fa
  if [ ! -s "$FA.done" ]; then
    echo "[$(date +%T)] extracting $chrom windows" | tee -a "$PROV"; : > "$FA"
    awk -F'\t' -v c=$chrom '$3==c' $C/d4z4_windows.tsv | while IFS=$'\t' read -r pub pansn ch hap fasta L start end; do
      samtools faidx "$fasta" "${pansn}:$((start+1))-${end}" | sed "1s/.*/>${pansn}/" >> "$FA"
    done
    samtools faidx "$FA" && touch "$FA.done"
  fi
  echo "  $chrom: $(grep -c '>' $FA) windows, $(awk '{s+=$2}END{print s}' $FA.fai) bp" | tee -a "$PROV"
done

# ---- 2+3. annotate + rebuild inputs + run caller, per k
for K in 21 31 51; do
  DB=HKS_human_CHM13_v2_k$K; A=$W/annot_k$K; rm -rf "$A"; mkdir -p "$A"
  echo "[$(date +%T)] annotate k=$K" | tee -a "$PROV"
  /usr/bin/time -v karyoscope annotate -i $W/fasta/windows_chr4.fa -i $W/fasta/windows_chr10.fa -o "$A" \
      --db $DB --db-root $DBROOT --feature-set region --feature-set subtelomeric --feature-set gene --feature-set repeat \
      -t 16 --no-bgzip --no-smooth > $W/logs/annotate_k$K.log 2> $W/logs/time_annotate_k$K.txt
  echo "  rc=$? $(grep Elapsed $W/logs/time_annotate_k$K.txt | sed 's/.*: //') peakRSS_kb=$(grep 'Maximum resident' $W/logs/time_annotate_k$K.txt | awk '{print $NF}')" | tee -a "$PROV"
  B=$W/base_k$K; rm -rf "$B"; mkdir -p "$B/data" "$B/results"
  python3 $W/d4z4_postprocess.py --annot-dir "$A" --db $DB --windows $C/d4z4_windows.tsv \
      --bin-script $C/bin_features_feb2026.py --outdir "$B/data" | tee -a "$PROV"
  echo "[$(date +%T)] caller k=$K" | tee -a "$PROV"
  ( cd $CODE && D4Z4_BASE_DIR=$B python analyze_d4z4_arrays.py > $B/results/stdout.log 2>&1 ); echo "  caller rc=$?" | tee -a "$PROV"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{a+=$h["num_arrays"]; m+=$h["repeat_count_main"]; t+=$h["repeat_count_total"]; n++} END{printf "  k=%s: haplotypes=%d arrays=%d repeat_main=%d repeat_total=%d\n", K, n, a, m, t}' K=$K $B/results/d4z4_array_analysis.tsv | tee -a "$PROV"
done

# ---- baseline on the PUBLISHED inputs (v1.1.0 annotations)
B=$W/base_published; rm -rf "$B"; mkdir -p "$B/results"; ln -sfn $PUB $B/data
( cd $CODE && D4Z4_BASE_DIR=$B python analyze_d4z4_arrays.py > $B/results/stdout.log 2>&1 ); echo "published-input caller rc=$?" | tee -a "$PROV"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{a+=$h["num_arrays"]; m+=$h["repeat_count_main"]; t+=$h["repeat_count_total"]; n++} END{printf "  published inputs: haplotypes=%d arrays=%d repeat_main=%d repeat_total=%d (expect 582/629/14433/28580)\n", n, a, m, t}' $B/results/d4z4_array_analysis.tsv | tee -a "$PROV"
echo "date_end: $(date -Is)" | tee -a "$PROV"
