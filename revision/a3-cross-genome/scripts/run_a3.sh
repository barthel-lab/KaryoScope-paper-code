#!/bin/bash
# A3: cross-genome demonstration — Ey15-2 annotated against the Col-CEN database.
#
# THREE annotate runs, not two. The brief says to strip organelles from Ey15-2,
# but the data is the other way round:
#   - Ey15-2 (GCA_946499665.1) has 5 sequences, all nuclear -> nothing to strip.
#   - Col-CEN v1.2 HAS ChrM (366,924 bp) and ChrC (154,478 bp), and the database
#     excludes both, so annotating the full file puts ~521 kb of unmatched
#     sequence into the BASELINE's novel fraction.
# The database excluding organelles is why the query needs filtering, not a
# substitute for it. So we run the baseline BOTH ways and measure the effect
# rather than asserting it. (Not necessarily 100% novel: Arabidopsis NUPT/NUMT
# insertions mean some organellar k-mers match nuclear labels and come back
# MISLABELLED instead of novel. Both outcomes pollute the baseline.)
#SBATCH --job-name=a3_crossgenome
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=a3_%j.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail

W=${KS_WORK}/a3_crossgenome
RES=$W/results; DBROOT=$W/db
mkdir -p "$RES" "$DBROOT" "$W/fa"
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
KS=$KSBIN/karyoscope
DATASETS=${KS_WORK}/envs/ncbi/bin/datasets
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:${KS_WORK}/envs/ncbi/bin:$PATH
COLCEN=${KS_WORK}/ks_build_sources/arabidopsis_colcen/Col-CEN_v1.2.fasta.gz
DB=HKS_arabidopsis_ColCEN

PROV=$RES/provenance.txt
{
  echo "=== A3 cross-genome ==="
  echo "date_start : $(date -Is)"
  echo "job        : ${SLURM_JOB_ID:-NA}   node: $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "governor   : $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo NA)"
  echo "karyoscope : $($KS --version 2>&1)"
  echo "hks        : $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "threads    : 16"
} | tee "$PROV"

# ---------------------------------------------------------------- database
echo "[$(date +%T)] fetching the database through the RELEASED path"
$KS download $DB --db-root "$DBROOT" > "$RES/db_download.log" 2>&1
echo "db_download_rc: $?" | tee -a "$PROV"
$KS info --db-root "$DBROOT" 2>&1 | head -20 | tee -a "$PROV"

# ---------------------------------------------------------------- Ey15-2
echo "[$(date +%T)] downloading GCA_946499665.1 (Ey15-2)"
cd "$W/fa"
[ -s ey15-2.zip ] || $DATASETS download genome accession GCA_946499665.1 --include genome \
    --filename ey15-2.zip > "$RES/datasets.log" 2>&1
echo "datasets_rc: $?" | tee -a "$PROV"
# `unzip` is not installed on these nodes; use python's zipfile (always present).
python3 -c "
import zipfile,sys
zipfile.ZipFile('ey15-2.zip').extractall('ey15-2_raw')
print('extracted', len(zipfile.ZipFile('ey15-2.zip').namelist()), 'entries')
" 2>>"$RES/datasets.log"
EY=$(ls ey15-2_raw/ncbi_dataset/data/GCA_946499665.1/*.fna 2>/dev/null | head -1)
if [ -z "$EY" ]; then echo "FATAL: Ey15-2 FASTA not found"; exit 1; fi
cp "$EY" "$W/fa/ey15-2.fa"; samtools faidx "$W/fa/ey15-2.fa"
{
  echo "--- Ey15-2 sequences (brief: keep nuclear only) ---"
  awk '{printf "    %-22s %12d\n",$1,$2}' "$W/fa/ey15-2.fa.fai"
  echo "    n_seqs: $(wc -l < "$W/fa/ey15-2.fa.fai")  total_bp: $(awk '{s+=$2}END{print s}' "$W/fa/ey15-2.fa.fai")"
  echo "    organelles present: $(grep -icE 'mitochond|chloroplast|ChrM|ChrC|plastid' "$W/fa/ey15-2.fa.fai" || true) (expected 0 -> filtering is a no-op)"
} | tee -a "$PROV"

# ---------------------------------------------------------------- Col-CEN nuclear-only
echo "[$(date +%T)] building nuclear-only Col-CEN (Chr1-Chr5)"
samtools faidx "$COLCEN" Chr1 Chr2 Chr3 Chr4 Chr5 > "$W/fa/colcen_nuclear.fa" 2>>"$RES/prep.log"
samtools faidx "$W/fa/colcen_nuclear.fa"
{
  echo "--- Col-CEN full vs nuclear-only ---"
  echo "    full    : $(awk '{s+=$2}END{print s}' ${COLCEN}.fai) bp over $(wc -l < ${COLCEN}.fai) seqs"
  echo "    nuclear : $(awk '{s+=$2}END{print s}' "$W/fa/colcen_nuclear.fa.fai") bp over $(wc -l < "$W/fa/colcen_nuclear.fa.fai") seqs"
  echo "    dropped : ChrM+ChrC = $(awk '$1=="ChrM"||$1=="ChrC"{s+=$2}END{print s}' ${COLCEN}.fai) bp"
} | tee -a "$PROV"

# ---------------------------------------------------------------- annotate x3
annotate_run () {   # $1=label  $2=fasta
  local label=$1 fa=$2 out=$W/annot_$1
  rm -rf "$out"; mkdir -p "$out"
  echo "[$(date +%T)] annotate: $label"
  /usr/bin/time -v $KS -v annotate -i "$fa" -o "$out" --db $DB --db-root "$DBROOT" \
      -t 16 > "$RES/annotate_$label.log" 2> "$RES/time_$label.txt"
  local rc=$?
  local wall=$(grep 'Elapsed (wall clock)' "$RES/time_$label.txt" | awk '{print $NF}')
  local rss=$(grep 'Maximum resident set size' "$RES/time_$label.txt" | awk '{print $NF}')
  echo "  $label rc=$rc wall=$wall rss_kb=$rss" | tee -a "$PROV"
  echo -e "${label}\t${rc}\t${wall}\t${rss}" >> "$RES/runtimes.tsv"
}
echo -e "run\texit\twall\tmax_rss_kb" > "$RES/runtimes.tsv"
annotate_run ey15_2          "$W/fa/ey15-2.fa"
annotate_run colcen_nuclear  "$W/fa/colcen_nuclear.fa"
annotate_run colcen_full     "$COLCEN"

# ---------------------------------------------------------------- novel fractions
echo "[$(date +%T)] novel fractions"
python3 - "$W" "$RES" <<'PY' | tee -a "$PROV"
import glob, gzip, os, sys, collections
W, RES = sys.argv[1], sys.argv[2]
rows=[]
for run in ("ey15_2","colcen_nuclear","colcen_full"):
    for f in sorted(glob.glob(os.path.join(W, f"annot_{run}", "*.smoothed.bed*"))):
        base=os.path.basename(f)
        if ".binned" in base or ".scaffolded" in base: continue
        fs=base.split(".")[-3] if base.endswith(".gz") else base.split(".")[-2]
        tot=collections.Counter(); nov=collections.Counter()
        op=gzip.open if f.endswith(".gz") else open
        with op(f,"rt") as fh:
            for line in fh:
                p=line.rstrip("\n").split("\t")
                if len(p)<4: continue
                try: ln=int(p[2])-int(p[1])
                except ValueError: continue
                seq=p[0]; lab=p[3]
                grp = "organelle" if seq in ("ChrM","ChrC") else "nuclear"
                tot[grp]+=ln
                if lab.lower()=="novel": nov[grp]+=ln
        t=sum(tot.values()); n=sum(nov.values())
        rows.append((run,fs,t,n,100.0*n/t if t else 0.0,
                     tot["organelle"], nov["organelle"]))
print("\n=== NOVEL FRACTION per run x feature set ===")
print(f"{'run':16}{'feature_set':14}{'total_bp':>14}{'novel_bp':>13}{'novel_%':>9}{'organelle_bp':>14}{'org_novel_bp':>13}")
for r in rows:
    print(f"{r[0]:16}{r[1]:14}{r[2]:14d}{r[3]:13d}{r[4]:9.3f}{r[5]:14d}{r[6]:13d}")
with open(os.path.join(RES,"novel_fraction.tsv"),"w") as o:
    o.write("run\tfeature_set\ttotal_bp\tnovel_bp\tnovel_pct\torganelle_bp\torganelle_novel_bp\n")
    for r in rows: o.write("\t".join(str(x) for x in r)+"\n")
PY

# ---------------------------------------------------------------- renders (BOTH fresh)
render () {  # $1=label $2=fasta $3=samplelabel
  local out=$W/karyotype_$1; rm -rf "$out"; mkdir -p "$out"
  echo "[$(date +%T)] karyotype: $1"
  $KS karyotype -i "$3=$2" --db $DB --db-root "$DBROOT" \
      --sample-label "$3" --telo-motif CCCTAAA --sex unknown \
      --format svg --format png -t 16 -o "$out" > "$RES/karyotype_$1.log" 2>&1
  echo "  $1 rc=$? files=$(ls "$out"/*.png 2>/dev/null | wc -l)" | tee -a "$PROV"
  cp "$out"/*.png "$RES/" 2>/dev/null
}
render ey15_2 "$W/fa/ey15-2.fa" Ey15-2
render colcen "$COLCEN" Col-CEN_v1.2

echo "date_end : $(date -Is)" >> "$PROV"
echo "=== DONE ==="; cat "$PROV"
