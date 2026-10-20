#!/bin/bash
# A1b: HumAS-HMMER (alpha-satellite HOR annotation) on HG002 v1.1 hap1.
# Brief asks for TWO numbers: (1) centromeric regions only, (2) full assembly.
# Centromere-restricted runs FIRST so the cheap number lands even if the
# overnight full-assembly scan is cut short.
#SBATCH --job-name=a1b_humas
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=a1b_humas_%j.log
: "${KS_HOME:?set KS_HOME (see the repository README, Paths)}"
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

B=${KS_WORK}/a1_bench
W=$B/a1b_humas; RES=$W/results
rm -rf "$W"; mkdir -p "$RES"
REPO=$B/HumAS-HMMER_for_AnVIL
HMM=$REPO/AS-HORs-hmmer3.3.2-120124.hmm     # per brief; repo also ships a 3.4 build
A1=${KS_WORK}/envs/a1tools
KS=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export PATH=$A1/bin:${KS_HOME}/software/bin:$KS:$PATH
IN=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz
THREADS=16

PROV=$RES/provenance.txt
{
  echo "=== A1b HumAS-HMMER run ==="
  echo "date_start : $(date -Is)"
  echo "slurm_job  : ${SLURM_JOB_ID:-NA}"
  echo "node       : $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket|Thread\(s\) per core'
  echo "threads    : $THREADS"
  echo "input      : $IN"
  echo "input_md5  : $(md5sum $IN | cut -d' ' -f1)"
  echo "nhmmer     : $($A1/bin/nhmmer -h 2>&1 | grep -m1 HMMER)"
  echo "bedmap     : $(bedmap --version 2>&1 | head -1)"
  echo "hmm        : $HMM"
  echo "hmm_md5    : $(md5sum $HMM | cut -d' ' -f1)"
  echo "hmm_models : $($A1/bin/hmmstat $HMM 2>/dev/null | grep -vc '^#')"
  echo "repo_commit: $(git -C $REPO rev-parse --short HEAD 2>/dev/null)"
} | tee "$PROV"

# ---------------------------------------------------------------- (1) centromeres
echo "[$(date +%T)] deriving centromere coordinates with karyoscope"
CENBED=$RES/centromeres_hap1.bed
/usr/bin/time -v $KS/karyoscope centromeres -i "$IN" \
    --db HKS_human_CHM13_v2 \
    --db-root ${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry \
    -o "$RES/centromeres" > "$RES/centromeres.log" 2>"$RES/time.centromeres.txt"
echo "karyoscope_centromeres_rc: $?" >> "$PROV"
find "$RES/centromeres" -name '*.bed' 2>/dev/null | head -5 >> "$PROV"
# take whatever BED it produced as the region list
cat "$RES/centromeres"/*.bed 2>/dev/null | awk 'NF>=3' > "$CENBED" || true
echo "centromere_regions: $(wc -l < "$CENBED" 2>/dev/null || echo 0)" | tee -a "$PROV"

run_scan () {   # $1 = label, $2 = dir of per-contig .fa
  local label=$1 dir=$2
  local n=$(ls "$dir"/*.fa 2>/dev/null | wc -l)
  local bases=$(cat "$dir"/*.fa 2>/dev/null | grep -v '^>' | tr -d '\n' | wc -c)
  echo "[$(date +%T)] scan '$label': $n file(s), $bases bases" | tee -a "$PROV"
  [ "$n" -gt 0 ] || { echo "  SKIP (no input)" | tee -a "$PROV"; return; }
  mkdir -p "$W/work_$label"; cd "$W/work_$label"
  cp "$REPO"/hmmertblout2bed.awk "$REPO"/overlap_filter.py . 2>/dev/null
  /usr/bin/time -v bash "$REPO/hmmer-run.sh" "$dir" "$HMM" "$THREADS" \
      > "$RES/hmmer_$label.log" 2>"$RES/time.hmmer_$label.txt"
  echo "  rc=$?" | tee -a "$PROV"
  {
    echo "--- $label ---"
    echo "  files=$n bases=$bases"
    grep -E "Elapsed \(wall clock\)|Maximum resident set size|Percent of CPU" "$RES/time.hmmer_$label.txt"
    echo "  HOR bed lines: $(cat AS-HOR-vs-*.bed 2>/dev/null | wc -l)"
  } | tee -a "$PROV"
  cp AS-HOR*-vs-*.bed "$RES/" 2>/dev/null
  cd "$W"
}

# per-contig FASTA for the centromeric regions
if [ -s "$CENBED" ]; then
  mkdir -p "$W/fa_cen"
  awk '{printf "%s:%d-%d\n",$1,($2+1),$3}' "$CENBED" > "$W/cen.regions"
  while read -r reg; do
    safe=$(echo "$reg" | tr ':' '_' | tr '-' '_')
    samtools faidx "$IN" "$reg" > "$W/fa_cen/${safe}.fa" 2>/dev/null
  done < "$W/cen.regions"
  run_scan centromeric "$W/fa_cen"
else
  echo "centromeric scan SKIPPED: no centromere BED produced" | tee -a "$PROV"
fi

# ---------------------------------------------------------------- (2) full assembly
echo "[$(date +%T)] splitting full assembly one sequence per file"
mkdir -p "$W/fa_full"
samtools faidx "$IN" 2>/dev/null
cut -f1 "${IN}.fai" | while read -r c; do
  samtools faidx "$IN" "$c" > "$W/fa_full/${c}.fa"
done
run_scan fullassembly "$W/fa_full"

echo "date_end   : $(date -Is)" >> "$PROV"
echo "=== DONE ==="; cat "$PROV"
