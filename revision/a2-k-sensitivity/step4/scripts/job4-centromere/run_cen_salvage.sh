#!/bin/bash
# Can the centromere major/minor calls be made k-stable by changing the FEATURE
# REPRESENTATION rather than the gate? The gate parameters are already known to be
# invariant (A6, 30-point grid) while k moves everything, which says the
# k-dependence is in the matrix, not the threshold. Five representations x k in
# {21,31,51}, on the v1-like annotations, gated exactly as published (0.5, 5).
#
# NOTE: KS_allchr_barplot's gate always rebuilds an EDGE-ONLY matrix regardless of
# what step 3 clustered on, so a representation change cannot remove k-sensitivity
# from the gate itself -- only from the cluster assignment. Measured, not assumed.
#SBATCH --job-name=a2s4_censalv
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=censalv_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_step4/job4_centromere
S=${KS_WORK}/a6_sweeps/cent/scripts
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
N=pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region
mkdir -p $W/salvage $W/labels
echo "date: $(date -Is) job: ${SLURM_JOB_ID} node: $(hostname)"
run_cfg () {  # $1=name $2=step3 extra flags $3=gate matrix type
  local name=$1 flags=$2 gmt=$3
  for K in 21 31 51; do
    local D=$W/cluster_v1like_k$K/out; local PASS=$D/$N.pass.bed
    local O=$W/salvage/${name}_k$K; mkdir -p $O
    $PY $S/KaryoScope_cluster_analysis.py --bed "$PASS" --output-prefix $O/st \
        --analysis-mode structure --max-sequence-length 50000000 --exclude-features "novel" \
        --k-selection silhouette --background white --no-umap $flags > $O/step3.log 2>&1
    local rc=$?
    [ -s $O/st.sequence_assignments.tsv ] || { echo "  $name k=$K step3 FAILED rc=$rc"; tail -3 $O/step3.log; continue; }
    $PY $W/cen_labels.py --assignments $O/st.sequence_assignments.tsv --bed "$PASS" \
        --k ${name}_$K --out $W/labels/salv_${name}_$K.tsv 2>/dev/null | sed "s/^/  /"
  done
}
run_cfg published  "--edges directional --no-abundance --matrix-type count_log1p_zscore_blockweight" count_log1p_zscore
run_cfg lengthw    "--edges directional --no-abundance --matrix-type length_weighted"                length_weighted
run_cfg abundonly  "--no-include-edges --abundance --matrix-type count_log1p_zscore_blockweight"     count_log1p_zscore
run_cfg abundlw    "--no-include-edges --abundance --matrix-type length_weighted"                    length_weighted
run_cfg edgesabund "--edges directional --abundance --matrix-type count_log1p_zscore_blockweight"    count_log1p_zscore
echo "done $(date -Is)"
