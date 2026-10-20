#!/bin/bash
# Per-phase timings for Fig. S6 from runs/index.tsv (same regexes as the old
# sweep's extract_phases.sh; the code is the same v2.3.0).
#   assignment = KMC `get_featureIDs` / HKS `hks lookup`   -> Mb/s
#   smoothing  = KMC Python pool ("smoothing pass:" -> "concat pass complete")
#                / HKS `hks smooth`
# Usage: extract_phases.sh [cold|warm]   (default cold) -> stdout TSV in the
# column layout plot_figS6.py reads.

MODE=${1:-cold}
W=$(cd "$(dirname "$0")" && pwd)   # this directory: runs/index.tsv + runs/*.log
GB=2950   # HG002 hap1, Mb (as in the published figure)
printf "backend\tthreads\trep\tassign_s\tassign_Mb_s\tsmooth_s\tpeak_rss_GiB\n"
awk -F'\t' -v m="$MODE" -v w="$W" 'NR>1 && $4==m {print $1, $2, $3, w "/" $10}' $W/runs/index.tsv |
while read -r b t r f; do
  rss=$(awk '/Maximum resident/{print $NF}' "$f")
  if [ "$b" = hks ]; then
    a=$(sed -nE "s/.*hks lookup for 'repeat' wrote .* in ([0-9.]+)s.*/\1/p" "$f" | head -1)
    s=$(sed -nE "s/.*hks smooth for 'repeat' wrote .* in ([0-9.]+)s.*/\1/p" "$f" | head -1)
  else
    a=$(sed -nE 's/.*ran get_featureIDs in ([0-9.]+)s.*/\1/p' "$f" | head -1)
    s=$(awk '/INFO: smoothing pass:/{split($1,x,":"); st=x[1]*3600+x[2]*60+x[3]}
             /concat pass complete/{split($1,x,":"); en=x[1]*3600+x[2]*60+x[3]; print en-st; exit}' "$f")
  fi
  [ -n "$a" ] || { echo "no assignment time in $f" >&2; continue; }
  printf "%s\t%s\t%s\t%s\t%.2f\t%s\t%.2f\n" "$b" "$t" "$r" "$a" \
      "$(echo "$GB/$a" | bc -l)" "${s:-NA}" "$(echo "${rss:-0}/1048576" | bc -l)"
done | sort -k1,1 -k2,2n -k3,3n
