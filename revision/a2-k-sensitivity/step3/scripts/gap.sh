#!/bin/bash
# How wide is the stretch between the two chromosome-specific segments -- i.e.
# sequence the chromosome featureset can assign to NEITHER parent?
# Per-sample values (fusion scaffold, partner chromosomes) come from $ROBS_CONFIG
# (a2_step3), through rob_samples.py next to this script.
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
HERE=$(cd "$(dirname "$0")" && pwd)
run(){ S=$1; SCAF=$2; C1=$3; C2=$4
  zcat ${KS_WORK}/a2_step3/gate_$S/work/$S.haplotype1.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz \
  | awk -F'\t' -v RES=500000 -v scaf="$SCAF" -v c1="$C1" -v c2="$C2" -v s="$S" '
     $1==scaf{ b=int($2/RES); n[b]+=$3-$2; v[b"|"$4]+=$3-$2; if(b>mx)mx=b }
     END{ for(k in v){split(k,a,"|"); if(v[k]>best[a[1]]){best[a[1]]=v[k]; lab[a[1]]=a[2]}}
          last1=-1; first2=-1
          for(i=0;i<=mx;i++){ if(lab[i]==c1) last1=i; if(lab[i]==c2 && first2<0 && i>0) first2=i }
          # first sustained c2: 3 consecutive c2-dominant windows after last1
          sus=-1; for(i=last1+1;i<=mx-2;i++){ if(lab[i]==c2&&lab[i+1]==c2&&lab[i+2]==c2){sus=i;break} }
          printf "%-8s last %s-dominant %7.1f Mb | first sustained %s %7.1f Mb | GAP %6.2f Mb\n", \
            s, c1, last1*RES/1e6, c2, sus*RES/1e6, (sus-last1)*RES/1e6 }'
}
while read -r S SCAF C1 C2; do
  run "$S" "$SCAF" "$C1" "$C2"
done < <(${PYTHON:-python3} "$HERE/rob_samples.py" sample scaffold left right)
