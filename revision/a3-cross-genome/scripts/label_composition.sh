#!/bin/bash
# Per-feature-set label composition for the A3 cross-genome panels.
# Reads the SAME scaffolded+smoothed BEDs the karyotype renders were drawn from,
# so the fractions describe what is actually on the panel.
# Annotations tile the input, so per-label bp sums to the genome length.
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
D=${KS_WORK}/a3_crossgenome
run() {
  local tag=$1 dir=$2 pre=$3
  for fs in chromosome region repeat gene; do
    f=$dir/$pre.$fs.smoothed.scaffolded.bed.gz
    [ -s "$f" ] || { echo "MISSING $f"; continue; }
    zcat "$f" | awk -F'\t' -v tag="$tag" -v fs="$fs" '
      {bp[$4]+=$3-$2; tot+=$3-$2}
      END{for(l in bp) printf "%s\t%s\t%s\t%d\t%.4f\n", tag, fs, l, bp[l], 100*bp[l]/tot}'
  done
}
{ printf "run\tfeature_set\tlabel\tbp\tpct\n"
  run ey15_2 $D/karyotype_ey15_2 ey15-2.HKS_arabidopsis_ColCEN
  run colcen $D/karyotype_colcen Col-CEN_v1.2.HKS_arabidopsis_ColCEN
} > $D/results/label_composition.tsv
wc -l $D/results/label_composition.tsv
