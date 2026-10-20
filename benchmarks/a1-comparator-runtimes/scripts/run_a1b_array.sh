#!/bin/bash
# A1b, restructured: ONE ARRAY TASK PER CONTIG/REGION.
#
# The brief said splitting the assembly one sequence per file makes nhmmer
# parallelise across files. It does not: hmmer-run.sh line 34 is a plain
# sequential `for filename in ...` loop, so splitting only chunks the work.
# Measured sequentially: chr10_hap1 (135.7 Mb) took 2h16 (~1.0 Mb/min), which
# projects to ~49 h for the full assembly (48 h limit) and ~46 h for the
# centromeric set (12 h limit). Both would have timed out.
#
# Each task runs the SAME nhmmer configuration (16 threads) as the sequential
# run, so per-item wall-clocks are unchanged and their SUM is the honest
# single-machine cost to report against KaryoScope's single sequential run.
# The array only buys us the answer sooner; it does not change the number.
#SBATCH --job-name=a1b_arr
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=a1b_arr_%A_%a.log
: "${KS_HOME:?set KS_HOME (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
B=${KS_WORK}/a1_bench
REPO=$B/HumAS-HMMER_for_AnVIL
HMM=$REPO/AS-HORs-hmmer3.3.2-120124.hmm
A1=${KS_WORK}/envs/a1tools
export PATH=$A1/bin:${KS_HOME}/software/bin:${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin:$PATH

MODE=${MODE:?set MODE=full or MODE=cen}
case $MODE in
  full) SRC=$B/a1b_humas/fa_full;        OUT=$B/a1b_arr_full;;
  cen)  SRC=$B/a1b_centromeric/fa_cen;   OUT=$B/a1b_arr_cen;;
esac
mkdir -p "$OUT/results"
FA=$(ls "$SRC"/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$FA" ] || { echo "no input for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }
NAME=$(basename "$FA" .fa)
WD="$OUT/work_$NAME"; rm -rf "$WD"; mkdir -p "$WD"; cd "$WD"
cp "$REPO"/hmmertblout2bed.awk "$REPO"/overlap_filter.py .
mkdir -p one; cp "$FA" one/

BP=$(grep -v '^>' "$FA" | tr -d '\n' | wc -c)
echo "task=${SLURM_ARRAY_TASK_ID} name=$NAME bases=$BP node=$(hostname) threads=16"
S=$(date +%s)
/usr/bin/time -v bash "$REPO/hmmer-run.sh" one "$HMM" 16 > "$OUT/results/${NAME}.log" 2> "$OUT/results/time_${NAME}.txt"
RC=$?
E=$(date +%s)
HOR=$(cat AS-HOR-vs-*.bed 2>/dev/null | wc -l)
printf "%s\t%s\t%d\t%d\t%d\t%s\n" "$MODE" "$NAME" "$BP" "$((E-S))" "$HOR" "$RC" >> "$OUT/results/timings.tsv"
echo "  rc=$RC wall=$((E-S))s HOR_lines=$HOR"
cp AS-HOR*-vs-*.bed "$OUT/results/" 2>/dev/null
rm -f "$WD"/*.out
# Propagate hmmer-run.sh's status. Without this the task exits 0 and SLURM
# records COMPLETED even when nhmmer failed, which is how chr13/21/22 were
# lost silently on the first pass.
exit $RC
