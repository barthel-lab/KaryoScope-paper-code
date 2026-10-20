#!/bin/bash
# Table S1 re-run under the Fig. S6 cold/warm protocol (../threadsweep_cold).
#
# Why: the published Table S1 (job 35247555, g-c-1-7-24, 2026-08-07) ran its 18
# runs back-to-back in one job with no cache control, in rep-major order; rep 1
# was the slowest run for every input x backend (8-20 %), i.e. cache state was
# mixed across replicates. Its per-run logs were purged from scratch.
#
# Design (identical protocol to the Fig. S6 sweep):
#  * ONE job, ONE node, 42 CPUs + 490 GB, so no other job can meaningfully share
#    the node; co-tenant count logged per run (others=);
#  * 18 points = {chm13, hg002_hap1, hg002_diploid} x {hks, kmc} x rep{1,2,3},
#    fixed random order (manifest.tsv, seed 20261013);
#  * before each COLD run: DB + input (+.fai/.gzi) + binaries + source evicted
#    from the page cache (posix_fadvise DONTNEED), fincore residency must be
#    < 64 MiB or the job aborts; then an identical WARM rerun;
#  * Table S1 design: all six feature sets (no --feature-set), -t 16;
#  * pinned: KaryoScope v2.3.0 worktree, hks 0.4.0 and get_featureIDs copies
#    (the same files the Fig. S6 sweep used);
#  * outputs to node-local /tmp, deleted between runs;
#  * resumable: points already in runs/index.tsv (warm row) are skipped.
#SBATCH --job-name=ks_tableS1_cold
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --mem=490G
#SBATCH --time=48:00:00
#SBATCH --output=tableS1_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -u
W=${KS_WORK}/ks_slurm/tableS1_cold
PIN=${KS_WORK}/ks_slurm/threadsweep_cold
ENVBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$ENVBIN:$PATH
export PYTHONPATH=$PIN/ks-v2.3.0/src
export KARYOSCOPE_HKS=$PIN/hks-0.4.0
export KARYOSCOPE_GET_FEATUREIDS=$PIN/get_featureIDs
KS=$ENVBIN/karyoscope
THREADS=16
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
declare -A DBS=( [hks]=HKS_human_CHM13_v2 [kmc]=KS_human_CHM13_v2 )
declare -A INPUTS=(
  [chm13]=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta
  [hg002_hap1]=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
  [hg002_diploid]=${KS_REFERENCES}/HG002/hg002v1.1.fasta.gz
)
LOCAL=/tmp/ks_tableS1_${SLURM_JOB_ID}
trap 'rm -rf "$LOCAL"' EXIT
INDEX=$W/runs/index.tsv
[ -s $INDEX ] || printf "input\tbackend\trep\tmode\tnode\tjob\tresident_bytes_before\tstart\tloadavg1\tlog\n" > $INDEX

echo "node: $(hostname)  job: $SLURM_JOB_ID  $(date -Is)"
lscpu | grep -E 'Model name|^CPU\(s\)|Socket|MHz'
echo "affinity: $(taskset -cp $$)"
echo "karyoscope: $($KS --version) from $(python -c 'import karyoscope;print(karyoscope.__file__)')"
echo "hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version.*')"
df -h /tmp | tail -1

while read -r idx name backend rep; do
  if grep -q "^${name}	${backend}	${rep}	warm	" $INDEX; then
    echo "skip $idx ($name $backend rep=$rep): done"; continue
  fi
  IN=${INPUTS[$name]}
  CACHED=( "$ROOT/${DBS[$backend]}" "$IN" "$KARYOSCOPE_HKS" "$KARYOSCOPE_GET_FEATUREIDS" "$PIN/ks-v2.3.0/src" )
  for ext in fai gzi; do [ -e "$IN.$ext" ] && CACHED+=( "$IN.$ext" ); done
  for mode in cold warm; do
    [ $mode = cold ] && python3 $PIN/evict.py "${CACHED[@]}"
    res=$($PIN/resident.sh "${CACHED[@]}")
    if [ $mode = cold ] && [ "$res" -gt $((64*1024*1024)) ]; then
      echo "ABORT: eviction failed for $idx ($res bytes resident)"; exit 1
    fi
    rm -rf "$LOCAL"; mkdir -p "$LOCAL"
    log=$W/runs/${name}_${backend}_r${rep}_${mode}.log
    start=$(date +%T); la=$(cut -d' ' -f1 /proc/loadavg)
    echo "=== [$idx] $name $backend rep=$rep $mode start $start resident=$res loadavg=$la others=$(squeue -h -w "$(hostname)" | grep -vc "$SLURM_JOB_ID")"
    /usr/bin/time -v "$KS" -v annotate -i "$IN" -o "$LOCAL" \
        --db "${DBS[$backend]}" --db-root "$ROOT" -t $THREADS > "$log" 2>&1
    rc=$?
    echo "=== [$idx] $mode done rc=$rc wall=$(grep 'Elapsed (wall' "$log" | awk '{print $NF}') out=$(du -sh "$LOCAL" | cut -f1)"
    if [ $rc -ne 0 ]; then echo "ABORT: run failed"; exit $rc; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$name" "$backend" "$rep" "$mode" \
        "$(hostname)" "$SLURM_JOB_ID" "$res" "$start" "$la" "$log" >> $INDEX
  done
done < $W/manifest.tsv
echo "ALL DONE $(date -Is)"
