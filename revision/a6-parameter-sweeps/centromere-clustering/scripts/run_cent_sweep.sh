#!/bin/bash
# =============================================================================
# A6 Sweep 3 -- centromere clustering (--sil-threshold x --centroid-sd)
# =============================================================================
# Stage 1 runs steps 1-3 of Yi-An's HPRC2_centromere_analysis.sh ONCE. Step 3
# (KaryoScope_cluster_analysis.py) takes NEITHER swept parameter -- they are
# consumed by steps 4-7, which read allchr_structure.sequence_assignments.tsv.
# So the clustering is done once and the sweep runs over the interpretation
# layer only.
#
# Stage 2 runs KS_allchr_barplot.py at each (sil, sd) pair. That script carries
# apply_sil_and_centroid(), "same logic as the dendrogram", and prints the
# per-chromosome Major/Outlier counts -- which is exactly the published quantity
# (the major/minor allele structure). Parsed from stdout; Yi-An's scripts are
# used unmodified.
#
#   --sil-threshold {0.0, 0.3, 0.4, 0.5, 0.6, 0.7}   0.0 = ungated control
#   --centroid-sd   {3, 4, 5, 6, 7}                  published (0.5, 5)
#   + 6 runs with --no-centroid-scan, to separate the silhouette stage from
#     the centroid stage (both create/destroy outliers, so the combined count
#     alone cannot attribute a change to one of them)
#
# Stage 3 runs the FULL steps 4-7 chain at the published (0.5, 5) as a control:
# if it does not reproduce the published figure, no sweep row is trustworthy.
#
# Inputs: see A6-DATA-LOCATIONS-for-YiAn.md. COLORS_DIR is the registry path,
# not the driver's default resources/databases path (old layout).
#SBATCH --job-name=a6_cent
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=a6cent_%j.log
: "${KS_HOME:?set KS_HOME (see the repository README, Paths)}"
: "${KS_PROJECTS:?set KS_PROJECTS (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a6_sweeps/cent
BEDDIR=${KS_PROJECTS}/KaryoScope_guarracino/results/pangenome/KaryoScope/KS_human_CHM13/presmoothed
QC="${KS_WORK}/aggregate_qc_v4.manual_curation (2).tsv"
COLORS=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/KS_human_CHM13
PY=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/python
export PATH=${KS_HOME}/software/bin:$PATH   # Rscript
S=$W/scripts
OUT=$W/out
NAME=pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed
PASS=$OUT/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed
ASGN=$OUT/allchr_structure.sequence_assignments.tsv

echo "=== A6 centromere-clustering sweep ==="
echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}  node: $(hostname)"
echo "beds:   $BEDDIR"
echo "qc:     $QC"
echo "colors: $COLORS"
echo "scripts: origin/pangenome_structure_V2 @ $(git -C ${KS_SOFTWARE}/github/barthel/KaryoScope-analysis rev-parse --short origin/pangenome_structure_V2)"
mkdir -p "$OUT" "$W/logs" "$W/sweep"

# ---- Stage 1: steps 1-3, once ----------------------------------------------
if [ -s "$PASS" ]; then
  echo "[$(date +%T)] steps 1-2 already done ($(wc -l < "$PASS") rows), reusing"
else
  echo "[$(date +%T)] step 1: concatenate 24 per-chromosome BEDs"
  T=$(mktemp -d -t HPRC2_concat.XXXXXX)
  ln -sfn "$BEDDIR" "$T/centromere_region_beds"; mkdir -p "$T/run"
  ( cd "$T/run" && Rscript "$S/concatenate_centromere_regions.R" && mv "$NAME" "$OUT/$NAME" ) || { echo "FATAL step 1"; exit 1; }
  rm -rf "$T"
  echo "  -> $(wc -l < "$OUT/$NAME") rows"

  echo "[$(date +%T)] step 2: NucFlag QC filter (Err=0 COLLAPSE=0 COLLAPSE_VAR=0)"
  ( cd "$OUT" && Rscript "$S/QCfilter_apply.R" "$OUT/$NAME" "$QC" centromere Err=0 COLLAPSE=0 COLLAPSE_VAR=0 ) || { echo "FATAL step 2"; exit 1; }
  # An empty pass.bed is the documented silent failure here (plain string join on
  # contig name), so fail loudly instead of clustering nothing.
  N=$(wc -l < "$PASS"); echo "  -> $N rows, $(cut -f1 "$PASS" | sort -u | wc -l) contigs"
  [ "$N" -gt 1000 ] || { echo "FATAL: pass.bed has $N rows -- contig join almost certainly failed"; exit 1; }
fi

if [ -s "$ASGN" ]; then
  echo "[$(date +%T)] step 3 already done, reusing $ASGN"
else
  echo "[$(date +%T)] step 3: per-chromosome clustering (Ward, structure, k in [2,10])"
  /usr/bin/time -v $PY "$S/KaryoScope_cluster_analysis.py" \
      --bed "$PASS" --output-prefix "$OUT/allchr_structure" \
      --analysis-mode structure --edges directional --no-abundance \
      --max-sequence-length 50000000 --exclude-features "novel" \
      --matrix-type count_log1p_zscore_blockweight \
      --k-selection silhouette --background white \
      > "$W/logs/step3.log" 2> "$W/logs/step3.time"
  RC=$?
  echo "  step3_rc=$RC wall=$(grep 'Elapsed (wall' "$W/logs/step3.time" | awk '{print $NF}') peakRSS=$(grep 'Maximum resident' "$W/logs/step3.time" | awk '{print $NF}')KB"
  [ $RC -eq 0 ] && [ -s "$ASGN" ] || { echo "FATAL step 3"; tail -30 "$W/logs/step3.log"; exit 1; }
fi
echo "  assignments: $(wc -l < "$ASGN") rows"

# ---- Stage 2: the sweep, steps 7a only -------------------------------------
run_point () {   # $1=sil $2=sd $3=scan(on|off)
  local sil=$1 sd=$2 scan=$3
  local tag="sil${sil}_sd${sd}_${scan}"
  local flag=""; [ "$scan" = "off" ] && flag="--no-centroid-scan"
  $PY "$S/KS_allchr_barplot.py" \
      --assignments "$ASGN" --bed "$PASS" \
      --output "$W/sweep/barplot_${tag}.svg" \
      --matrix-type count_log1p_zscore \
      --sil-threshold "$sil" --centroid-sd "$sd" $flag \
      > "$W/sweep/${tag}.stdout" 2> "$W/sweep/${tag}.stderr"
  local rc=$?
  printf "  %-28s rc=%s %s\n" "$tag" "$rc" \
    "$(grep -m1 'Grand total' "$W/sweep/${tag}.stdout" 2>/dev/null | sed 's/^ *//')"
  return $rc
}

echo
echo "[$(date +%T)] Stage 2: 30 combinations + 6 silhouette-only"
FAIL=0
for sil in 0.0 0.3 0.4 0.5 0.6 0.7; do
  for sd in 3 4 5 6 7; do run_point $sil $sd on || FAIL=$((FAIL+1)); done
  run_point $sil 5 off || FAIL=$((FAIL+1))
done
echo "stage2 failures: $FAIL"

# ---- Stage 3: the published-value control, full steps 4-7 ------------------
echo
echo "[$(date +%T)] Stage 3: full steps 4-7 at the published (0.5, 5)"
CTL=$W/control; mkdir -p "$CTL"
$PY "$S/KS_allchr_dendrogram.py" --assignments "$ASGN" --bed "$PASS" \
    --colors "$COLORS" --output "$CTL/allchr_dendrogram.svg" \
    --matrix-type count_log1p_zscore --sil-threshold 0.5 --centroid-sd 5 \
    --row-height 12 --bar-height 10 > "$W/logs/ctl_step4.log" 2>&1
echo "  step4 rc=$?"
$PY "$S/KS_allchr_annotate.py" --svg "$CTL/allchr_dendrogram.svg" \
    --assignments "$ASGN" --bed "$PASS" \
    --output "$CTL/allchr_outlier_annotations.tsv" > "$W/logs/ctl_step5.log" 2>&1
echo "  step5 rc=$?"
$PY "$S/KS_allchr_dendrogram.py" --assignments "$ASGN" --bed "$PASS" \
    --colors "$COLORS" --output "$CTL/allchr_dendrogram_annotated.svg" \
    --matrix-type count_log1p_zscore --sil-threshold 0.5 --centroid-sd 5 \
    --row-height 12 --bar-height 10 \
    --annotations "$CTL/allchr_outlier_annotations.tsv" > "$W/logs/ctl_step6.log" 2>&1
echo "  step6 rc=$?"
$PY "$S/KS_allchr_allele_heatmap.py" --assignments "$ASGN" --bed "$PASS" \
    --output-dir "$CTL" --matrix-type count_log1p_zscore \
    --sil-threshold 0.5 --centroid-sd 5 > "$W/logs/ctl_step7b.log" 2>&1
echo "  step7b rc=$?"

echo
echo "[$(date +%T)] done. Collate with: $PY $W/collate.py"
