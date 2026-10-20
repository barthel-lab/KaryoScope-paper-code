#!/bin/bash
#SBATCH --job-name=rm_smoke2
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=smoke_rm2_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a1_bench/smoke2
rm -rf "$W"; mkdir -p "$W/out"
RM=${KS_SOFTWARE}/RepeatMasker
A1=${KS_WORK}/envs/a1tools
# Point RepeatMasker at our private rmblast WITHOUT editing the shared config.
export RMBLAST_DIR=$A1/bin
export PATH=$RM:$RM/util:$A1/bin:${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz

echo "node: $(hostname)"
echo "rmblastn: $($A1/bin/rmblastn -version | head -1)"
echo "RepeatMasker: $($RM/RepeatMasker -v 2>&1 | head -1)"
echo "input md5: $(md5sum $IN | cut -d' ' -f1)  path: $IN"

# 5 Mb of chr1_hap1, away from the telomere and centromere so the masked
# fraction is broadly representative of bulk euchromatin.
REGION="chr1_hap1:20000000-25000000"
samtools faidx "$IN" "$REGION" > "$W/test.fa"
echo "test.fa bases: $(grep -v '^>' "$W/test.fa" | tr -d '\n' | wc -c)"

/usr/bin/time -v $RM/RepeatMasker -engine rmblast -species human -pa 2 \
    -dir "$W/out" "$W/test.fa" 2> "$W/out/time.txt"
echo "--- RepeatMasker rc=$? ---"

echo "=== engine actually used ==="
grep -iE "search engine" "$W/out"/*.tbl 2>/dev/null
echo "=== masked fraction + families (the library-load check) ==="
sed -n '1,32p' "$W/out"/*.tbl 2>/dev/null
echo "=== wall-clock / RSS ==="
grep -E "Elapsed \(wall|Maximum resident" "$W/out/time.txt"
