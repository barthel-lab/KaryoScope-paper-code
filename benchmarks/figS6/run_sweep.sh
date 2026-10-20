#!/bin/bash
# Fig. S6 thread sweep, redesigned so cache state cannot confound thread count.
#
# Why the old sweep (jobs 35376001/35376043) is not reused:
#  * thread counts ran in DESCENDING order, so the one cold-index run landed on
#    t=32 (the shortest point) -> the 104 s replicate;
#  * rep 3 came from a second job and was 5-30% faster at all 12 points;
#  * both jobs shared their node with other tenants (job 2's CPU affinity was a
#    scattered 32-of-48 set; the node showed 25 GB free of 515 GB), so memory
#    pressure and other tenants could evict/compete at random.
#
# This design:
#  * ONE job, ONE node, 42 CPUs (the most a job may request here; MaxCPUsPerNode=44 is not reachable) and 490 GB, so
#    no other job can meaningfully share the node; every point sees the same
#    hardware and there is no job-to-job offset;
#  * the 36 (backend, threads, rep) points run in a FIXED RANDOM order
#    (manifest.tsv, seed 20261012), so any drift is not aligned with threads;
#  * before each timed run every DB / input / binary / source file is evicted
#    from the page cache (posix_fadvise DONTNEED) and residency is checked with
#    fincore (must be ~0): every point is a COLD start, as a user's run would be
#    if it were the only thing they were running;
#  * each cold run is followed by an identical WARM run (cache left as the cold
#    run left it), so the cold-start cost is measured, not guessed;
#  * code pinned: KaryoScope v2.3.0 worktree + copied hks 0.4.0 binary;
#  * outputs go to node-local disk and are deleted between runs;
#  * resumable: points already in runs/index.tsv are skipped.
#SBATCH --job-name=ks_coldsweep
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --mem=490G
#SBATCH --time=48:00:00
#SBATCH --output=coldsweep_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -u
W=${KS_WORK}/ks_slurm/threadsweep_cold
ENVBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$ENVBIN:$PATH
export PYTHONPATH=$W/ks-v2.3.0/src
export KARYOSCOPE_HKS=$W/hks-0.4.0
export KARYOSCOPE_GET_FEATUREIDS=$W/get_featureIDs
KS=$ENVBIN/karyoscope
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
INPUT=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
declare -A DBS=( [hks]=HKS_human_CHM13_v2 [kmc]=KS_human_CHM13_v2 )
LOCAL=/tmp/ks_coldsweep_${SLURM_JOB_ID}
trap 'rm -rf "$LOCAL"' EXIT
mkdir -p $W/runs
INDEX=$W/runs/index.tsv
[ -s $INDEX ] || printf "backend\tthreads\trep\tmode\tnode\tjob\tresident_bytes_before\tstart\tloadavg1\tlog\n" > $INDEX

echo "node: $(hostname)  job: $SLURM_JOB_ID  $(date -Is)"
lscpu | grep -E 'Model name|^CPU\(s\)|Socket|MHz'
echo "affinity: $(taskset -cp $$)"
echo "karyoscope: $($KS --version) from $(python -c 'import karyoscope;print(karyoscope.__file__)')"
echo "hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version.*')"

while read -r idx backend t rep; do
  if grep -q "^${backend}	${t}	${rep}	warm	" $INDEX; then
    echo "skip $idx ($backend t=$t rep=$rep): done"; continue
  fi
  DB=$ROOT/${DBS[$backend]}
  CACHED=( "$DB" "$INPUT" "$KARYOSCOPE_HKS" "$KARYOSCOPE_GET_FEATUREIDS" "$W/ks-v2.3.0/src" )
  for mode in cold warm; do
    [ $mode = cold ] && python3 $W/evict.py "${CACHED[@]}"
    res=$($W/resident.sh "${CACHED[@]}")
    if [ $mode = cold ] && [ "$res" -gt $((64*1024*1024)) ]; then
      echo "ABORT: eviction failed for $idx ($res bytes resident)"; exit 1
    fi
    rm -rf "$LOCAL"; mkdir -p "$LOCAL"
    log=$W/runs/${backend}_t${t}_r${rep}_${mode}.log
    start=$(date +%T); la=$(cut -d' ' -f1 /proc/loadavg)
    echo "=== [$idx] $backend t=$t rep=$rep $mode start $start resident=$res loadavg=$la others=$(squeue -h -w "$(hostname)" | grep -vc "$SLURM_JOB_ID")"
    /usr/bin/time -v "$KS" -v annotate -i "$INPUT" -o "$LOCAL" \
        --db "${DBS[$backend]}" --db-root "$ROOT" \
        --feature-set repeat -t "$t" > "$log" 2>&1
    rc=$?
    echo "=== [$idx] $mode done rc=$rc wall=$(grep 'Elapsed (wall' "$log" | awk '{print $NF}')"
    if [ $rc -ne 0 ]; then echo "ABORT: run failed"; exit $rc; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$backend" "$t" "$rep" "$mode" \
        "$(hostname)" "$SLURM_JOB_ID" "$res" "$start" "$la" "$log" >> $INDEX
  done
done < $W/manifest.tsv
echo "ALL DONE $(date -Is)"
