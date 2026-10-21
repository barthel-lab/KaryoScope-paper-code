#!/bin/bash
# Extract SST1 blocks on each ROB fusion scaffold from our own acrocentric
# annotations, merging intervals separated by < GAP bp into one block.
# The fusion scaffolds come from $ROBS_CONFIG (a2_step3), through rob_samples.py
# next to this script.
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
HERE=$(cd "$(dirname "$0")" && pwd)
GAP=${GAP:-20000}
W=${KS_WORK}/a2_step3
while read -r S SCAF; do
  echo "=== $S  $SCAF (gap<${GAP}) ==="
  zcat $W/$S/work/$S.haplotype1.HKS_human_CHM13_v2.acrocentric.smoothed.scaffolded.bed.gz \
  | awk -F'\t' -v scaf="$SCAF" '$1==scaf && $4=="SST1"{print $2"\t"$3}' \
  | sort -k1,1n \
  | awk -v gap="$GAP" 'NR==1{s=$1;e=$2;n=$2-$1;next}
       {if($1-e<gap){e=$2;n+=$2-$1}else{printf "  block %12d - %12d  span %8.1f kb  SST1 %8.1f kb  mid %12d\n",s,e,(e-s)/1000,n/1000,int((s+e)/2); s=$1;e=$2;n=$2-$1}}
       END{printf "  block %12d - %12d  span %8.1f kb  SST1 %8.1f kb  mid %12d\n",s,e,(e-s)/1000,n/1000,int((s+e)/2)}'
done < <(${PYTHON:-python3} "$HERE/rob_samples.py" sample scaffold)
