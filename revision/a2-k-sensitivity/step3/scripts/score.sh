#!/bin/bash
# Independent implementation of the published cumulative k-mer density method,
# run over the WHOLE fusion scaffold (no detection window), from our own
# chromosome-featureset annotations.
# Per-sample values (fusion scaffold, partner chromosomes, SST1 midpoint and the
# breakpoint the cytoband-pivot pipeline reported) come from $ROBS_CONFIG
# (a2_step3), through rob_samples.py next to this script.
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
: "${ROBS_CONFIG:?set ROBS_CONFIG to your filled ROB configuration (robs/README.md)}"
HERE=$(cd "$(dirname "$0")" && pwd)
run() { S=$1; SCAF=$2; C1=$3; C2=$4; MID=$5; CALLED=$6
  zcat ${KS_WORK}/a2_step3/gate_$S/work/$S.haplotype1.HKS_human_CHM13_v2.chromosome.smoothed.scaffolded.bed.gz \
  | awk -F'\t' -v RES=10000 -v scaf="$SCAF" -v c1="$C1" -v c2="$C2" -v mid="$MID" -v called="$CALLED" -v s="$S" '
     $1==scaf && ($4==c1||$4==c2){ b=int($2/RES)
       if($4==c1){a[b]+=$3-$2; ta+=$3-$2} else {d[b]+=$3-$2; td+=$3-$2}; if(b>mx)mx=b }
     END{ cum=0; for(i=0;i<=mx;i++) cumA[i]=cum+=a[i]
          cum=0; for(i=mx;i>=0;i--) cumD[i]=cum+=d[i]
          best=-1; for(i=0;i<=mx;i++){sc=cumA[i]/ta+cumD[i]/td; if(sc>best){best=sc;bi=i}}
          gm=bi*RES; dm=gm-mid; adm=(dm<0?-dm:dm); dc=called-mid; adc=(dc<0?-dc:dc)
          printf "%-8s argmax %9.3f Mb  (%+7.1f kb from SST1)  %-4s   |  pipeline %9.3f Mb  (%+8.1f kb)  %s\n", \
            s, gm/1e6, dm/1000, (adm<=525000?"PASS":"FAIL"), called/1e6, dc/1000, (adc<=525000?"PASS":"FAIL") }'
}
printf "%-8s %-45s   |  %s\n" "sample" "independent method, full scaffold" "pipeline (gate run)"
while read -r S SCAF C1 C2 MID CALLED; do
  run "$S" "$SCAF" "$C1" "$C2" "$MID" "$CALLED"
done < <(${PYTHON:-python3} "$HERE/rob_samples.py" sample scaffold left right sst1_mid cytoband_pivot_breakpoint)
