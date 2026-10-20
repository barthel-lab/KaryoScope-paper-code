#!/bin/bash
# Re-derive HMM-based D4Z4 repeat-unit counts for the 249 chr4 windows, following
# the manuscript Methods ("HMM-based D4Z4 repeat-unit counting") step by step:
# CHM13 chr4:193,538,193-193,541,504 -> hmmbuild -> nhmmer --cpu 12 --notextw --noali
# --tblout -> hmmertblout2bed.awk th=0.7 -> keep hits > 2 kb -> count per contig.
# This is NOT Hailey Loucks's count table (HPRC_chr4_counts.xlsx, not on the
# cluster); it is the same recipe applied to the same windows, so it is the
# k-independent reference the k-vs-HMM correlation needs.
#SBATCH --job-name=a2s4_hmm
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=hmm_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
H=${KS_WORK}/a2_step4/job3_d4z4/hmm; C=${KS_WORK}/a2_step4/common
A1=${KS_WORK}/envs/a1tools/bin
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
cd $H
echo "date: $(date -Is) job: ${SLURM_JOB_ID:-NA} node: $(hostname)"; $A1/nhmmer -h | sed -n '2p'
FA=$H/windows_chr4.fa
if [ ! -s "$FA.fai" ]; then
  : > "$FA"
  awk -F'\t' '$3=="chr4"' $C/d4z4_windows.tsv | while IFS=$'\t' read -r pub pansn ch hap fasta L start end; do
    samtools faidx "$fasta" "${pansn}:$((start+1))-${end}" | sed "1s/.*/>${pansn}/" >> "$FA"
  done
  samtools faidx "$FA"
fi
echo "windows: $(wc -l < $FA.fai)"
$A1/hmmbuild D4Z4rep.hmm D4Z4rep.fa > hmmbuild.log 2>&1; echo "hmmbuild rc=$?"; grep -E "^D4Z4|LENG|NSEQ" D4Z4rep.hmm | head -3
/usr/bin/time -v $A1/nhmmer --cpu 12 --notextw --noali --tblout windows_chr4.out D4Z4rep.hmm $FA > nhmmer.stdout 2> time_nhmmer.txt
echo "nhmmer rc=$? $(grep Elapsed time_nhmmer.txt)"
gawk -v th=0.7 -f hmmertblout2bed.awk windows_chr4.out > windows_chr4.bed
echo "bed rows: $(wc -l < windows_chr4.bed)"
# hits > 2 kb per contig (Methods: "hits larger than 2 kb were counted")
awk -F'\t' 'BEGIN{OFS="\t"} ($3-$2)>2000 {n[$1]++} END{for(c in n) print c, n[c]}' windows_chr4.bed | sort > hmm_counts_chr4.tsv
# every window present, zero if no hits
awk -F'\t' 'NR==FNR{n[$1]=$2; next} {print $1"\t"(($1 in n)?n[$1]:0)}' hmm_counts_chr4.tsv <(cut -f1 $FA.fai) | sort > hmm_counts_chr4.all.tsv
echo "contigs with counts: $(wc -l < hmm_counts_chr4.tsv) / $(wc -l < hmm_counts_chr4.all.tsv); total hits>2kb: $(awk '{s+=$2}END{print s}' hmm_counts_chr4.all.tsv)"
awk -F'\t' '($3-$2)>2000{print $3-$2}' windows_chr4.bed | sort -n | awk '{a[NR]=$1} END{print "hit length >2kb: min",a[1],"median",a[int(NR/2)],"max",a[NR]}'
echo "done $(date -Is)"
