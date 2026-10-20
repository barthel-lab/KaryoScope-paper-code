#!/bin/bash
# A1b, part 1 (centromere-restricted) — REDO.
#
# The in-job attempt produced 0 regions for two reasons, both found by
# inspecting output that already existed:
#   1. `karyoscope centromeres` writes the BED GZIPPED
#      (*.centromeres.bed.gz); the extraction globbed *.bed and matched nothing.
#   2. Contig names in that BED are SCAFFOLDED names (chr1_hap1_chr1_hap1),
#      which do not exist in the source FASTA. Translation must go through
#      scaffold_map.tsv (new_name -> original_name), not string surgery.
# Verified before running: all 23 contigs are flipped=no, so coordinates need
# no reversal. The map/fai length delta is a constant -30 = k-1 for k=31 (the
# map counts k-mer positions), negligible over multi-Mb regions.
#
# The full-assembly scan is running in job 35377064 and is untouched by this.
#SBATCH --job-name=a1b_cen
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=a1b_cen_%j.log
: "${KS_HOME:?set KS_HOME (see the repository README, Paths)}"
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

B=${KS_WORK}/a1_bench
SRC=$B/a1b_humas/results/centromeres
W=$B/a1b_centromeric; RES=$W/results
rm -rf "$W"; mkdir -p "$RES" "$W/fa_cen" "$W/work"
REPO=$B/HumAS-HMMER_for_AnVIL
HMM=$REPO/AS-HORs-hmmer3.3.2-120124.hmm
A1=${KS_WORK}/envs/a1tools
KS=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$A1/bin:${KS_HOME}/software/bin:$KS:$PATH
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
THREADS=16

PROV=$RES/provenance.txt
{
  echo "=== A1b part 1: HumAS-HMMER on centromeric regions ==="
  echo "date_start : $(date -Is)"
  echo "slurm_job  : ${SLURM_JOB_ID:-NA}"
  echo "node       : $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket|Thread\(s\) per core'
  echo "governor   : $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo NA)"
  echo "no_turbo   : $(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo NA)"
  echo "threads    : $THREADS   (matches Table S1's thread count)"
  echo "input      : $IN"
  echo "input_md5  : $(md5sum $IN | cut -d' ' -f1)"
  echo "nhmmer     : $($A1/bin/nhmmer -h 2>&1 | grep -m1 HMMER)"
  echo "hmm        : $(basename $HMM)  md5=$(md5sum $HMM | cut -d' ' -f1)"
  echo "hmm_models : $($A1/bin/hmmstat $HMM 2>/dev/null | grep -vc '^#')"
  echo "cen_source : $SRC (from karyoscope centromeres, job 35377064)"
} | tee "$PROV"

MAP=$(ls $SRC/*scaffold_map.tsv)
CENGZ=$(ls $SRC/*.centromeres.bed.gz)
echo "flipped_contigs: $(awk -F'\t' 'NR>1 && $6!="no"' $MAP | wc -l) (must be 0)" | tee -a "$PROV"

# scaffolded name -> original FASTA name, then extract each region
zcat "$CENGZ" | awk -F'\t' 'NF>=3' | while IFS=$'\t' read -r sc s e rest; do
  orig=$(awk -F'\t' -v n="$sc" 'NR>1 && $1==n{print $2; exit}' "$MAP")
  [ -n "$orig" ] || { echo "WARN: no mapping for $sc" >> "$PROV"; continue; }
  reg="${orig}:$((s+1))-${e}"
  samtools faidx "$IN" "$reg" > "$W/fa_cen/${orig}_${s}_${e}.fa" 2>/dev/null
done

N=$(ls "$W/fa_cen"/*.fa 2>/dev/null | wc -l)
BASES=$(cat "$W/fa_cen"/*.fa 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c)
{
  echo "regions_extracted : $N (expect 23)"
  echo "centromeric_bases : $BASES"
  echo "fraction_of_asm   : $(awk -v b=$BASES 'BEGIN{printf "%.1f%%", 100*b/2947895164}')"
} | tee -a "$PROV"
[ "$N" -gt 0 ] || { echo "FATAL: no regions extracted"; exit 1; }

cd "$W/work"
cp "$REPO"/hmmertblout2bed.awk "$REPO"/overlap_filter.py .
echo "[$(date +%T)] nhmmer over centromeric regions"
/usr/bin/time -v bash "$REPO/hmmer-run.sh" "$W/fa_cen" "$HMM" "$THREADS" \
    > "$RES/hmmer_centromeric.log" 2> "$RES/time.hmmer_centromeric.txt"
echo "hmmer_rc: $?" | tee -a "$PROV"

{
  echo "=== headline numbers (centromeric regions only) ==="
  grep -E "Elapsed \(wall clock\)|Maximum resident set size|Percent of CPU" "$RES/time.hmmer_centromeric.txt"
  echo "HOR bed lines : $(cat AS-HOR-vs-*.bed 2>/dev/null | wc -l)"
  echo "distinct HORs : $(cat AS-HOR-vs-*.bed 2>/dev/null | awk '{print $4}' | sort -u | wc -l)"
  echo "(a haplotype should yield thousands of HOR calls across ~20 centromeres;"
  echo " zero means the models or the regions are wrong)"
} | tee -a "$PROV"
cp AS-HOR*-vs-*.bed "$RES/" 2>/dev/null
echo "date_end   : $(date -Is)" >> "$PROV"
echo "=== DONE ==="; cat "$PROV"
