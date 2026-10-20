#!/bin/bash
#SBATCH --job-name=a2s2_k11
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --array=1-2
#SBATCH --output=a2s2_k11_%A_%a.log
# HARD 3h CAP by design. Step 1 showed the k=11 index contains 99.997% of all
# possible 11-mers, so 99.80% of s-mers resolve only to internal nodes. A
# timeout here is itself the reportable result ("k=11 is unusable in practice"),
# not a failure to chase.
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
export SLURM_ARRAY_TASK_ID_ORIG=$SLURM_ARRAY_TASK_ID
if [ "$SLURM_ARRAY_TASK_ID" = "1" ]; then export ASM_OVERRIDE=chm13; else export ASM_OVERRIDE=hg002hap1; fi
W=${KS_WORK}/a2_ksweep
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
DBROOT=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db
if [ "$ASM_OVERRIDE" = chm13 ]; then FA=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta; else FA=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz; fi
OUT=$W/annot_${ASM_OVERRIDE}_k11; RES=$W/results; mkdir -p "$RES"; rm -rf "$OUT"; mkdir -p "$OUT"
echo "node=$(hostname) asm=$ASM_OVERRIDE k=11 (capped at 3h by design)"
/usr/bin/time -v $KSBIN/karyoscope -v annotate -i "$FA" -o "$OUT" --db HKS_human_CHM13_v2_k11 \
    --db-root "$DBROOT" -t 16 --no-bgzip > "$RES/annotate_${ASM_OVERRIDE}_k11.log" 2> "$RES/time_${ASM_OVERRIDE}_k11.txt"
rc=$?
echo "annotate_rc=$rc wall=$(grep 'Elapsed (wall' "$RES/time_${ASM_OVERRIDE}_k11.txt" | awk '{print $NF}')"
printf "%s\t11\t%s\t%s\t%s\n" "$ASM_OVERRIDE" "$rc" \
  "$(grep 'Elapsed (wall' "$RES/time_${ASM_OVERRIDE}_k11.txt" | awk '{print $NF}')" \
  "$(grep 'Maximum resident' "$RES/time_${ASM_OVERRIDE}_k11.txt" | awk '{print $NF}')" >> "$RES/step2_runtimes.tsv"
du -sh "$OUT" 2>/dev/null | sed 's/^/output size: /'
