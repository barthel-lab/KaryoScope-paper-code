#!/bin/bash
#SBATCH --job-name=rm_smoke
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=smoke_rm_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a1_bench/smoke
rm -rf "$W"; mkdir -p "$W/out"
RM=${KS_SOFTWARE}/RepeatMasker
A1=${KS_WORK}/envs/a1tools
# Point RepeatMasker at our private rmblast WITHOUT editing the shared config.
export RMBLAST_DIR=$A1/bin
export PATH=$RM:$RM/util:$A1/bin:${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz

echo "node: $(hostname)"; echo "rmblastn: $($A1/bin/rmblastn -version | head -1)"
echo "RepeatMasker: $($RM/RepeatMasker -v 2>&1 | head -1)"

# ~2 Mb of a single contig: enough to exercise engine + library load.
CONTIG=$(samtools faidx "$IN" 2>/dev/null && samtools idxstats "$IN" 2>/dev/null | head -1 | cut -f1)
CONTIG=${CONTIG:-$(zcat "$IN" | head -1 | tr -d '>' | awk '{print $1}')}
echo "contig: $CONTIG"
samtools faidx "$IN" "${CONTIG}:1-2000000" > "$W/test.fa"
grep -c "" "$W/test.fa"

/usr/bin/time -v $RM/RepeatMasker -engine rmblast -species human -pa 2 \
    -dir "$W/out" "$W/test.fa" 2> "$W/out/time.txt"
echo "--- RepeatMasker rc=$? ---"

echo "=== .tbl summary (the sanity check) ==="
cat "$W/out"/*.tbl 2>/dev/null | head -40
echo "=== did rmblast actually get used? ==="
grep -iE "search engine|rmblast|abblast|crossmatch" "$W/out"/*.log "$W/out"/*.tbl 2>/dev/null | head
tail -20 "$W/out/time.txt"
