#!/bin/bash
# A6 §6 inputs: the RepeatMasker release the SUBTELOMERIC set was built from,
# plus seqtk telo regions. Note this is the S3 BED, NOT the UCSC bigBed the
# `repeat` set uses -- the recipe is explicit that they are different
# representations and not interchangeable.
#SBATCH --job-name=a6_tar1_fetch
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fetch_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export PATH=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH
W=${KS_WORK}/a6_sweeps/tar1_its
RM=$W/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
FA=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta
echo "=== fetch  $(date -Is)  node $(hostname) ==="
if [ ! -s "$RM" ]; then
  curl -fL --retry 3 -o "$RM" \
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
  echo "curl_rc=$?"
fi
echo "bytes=$(stat -c%s "$RM")  (expected 342862439)"
[ "$(stat -c%s "$RM")" -eq 342862439 ] || { echo "FATAL size mismatch"; exit 1; }
echo "sha256: $(sha256sum "$RM" | cut -d' ' -f1)"
echo "=== column layout (first 3 rows, tab-delimited) ==="
head -3 "$RM" | cat -A | sed 's/\$$//' | head -3
echo "=== field count ==="
head -1 "$RM" | awk -F'\t' '{print "  NF =", NF}'
echo "=== candidate name/class fields ==="
head -200000 "$RM" | awk -F'\t' '{for(i=4;i<=NF && i<=12;i++) c[i"|"$i]++} END{for(k in c) if(c[k]>500) print "  ", k, c[k]}' | sort -t'|' -k1,1n | head -25
echo "=== TAR1 hits (any field) ==="
grep -c "TAR1" "$RM"
echo "=== seqtk telo ==="
seqtk telo "$FA" > $W/chm13v2.0.telo.bed 2> $W/telo.err
echo "telo_rc=$?  regions=$(wc -l < $W/chm13v2.0.telo.bed)  bp=$(awk '{s+=$3-$2} END{print s+0}' $W/chm13v2.0.telo.bed)"
head -3 $W/chm13v2.0.telo.bed
echo "done $(date -Is)"
