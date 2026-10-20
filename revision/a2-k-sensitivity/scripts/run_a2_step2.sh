#!/bin/bash
# A2 Step 2: annotation sweep. One array task per (assembly x k).
# Annotates on the HKS backend at 16 threads keeping the pre-smoothing output,
# then computes the cheap per-base metrics (leaf / internal-node / novel by
# feature set). Concordance vs k=31 is a SEPARATE phase, because it needs every
# k=31 run to exist first.
#
# k=11 is deliberately NOT in this array. Step 1 showed the k=11 index holds
# 4,194,198 of the 4,194,304 possible 11-mers (99.997% of the entire 11-mer
# space), so 99.80% of s-mers resolve only to internal nodes. It is submitted
# separately with a hard cap: informative for panel A, expected to be
# degenerate here, and not allowed to block the rest.
#
# Tasks: 1-5 = CHM13 x k(21,31,41,51,61); 6-10 = HG002 hap1 x same.
#SBATCH --job-name=a2_step2
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=1-10%5
#SBATCH --output=a2s2_%A_%a.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a2_ksweep
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
KS=$KSBIN/karyoscope
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
DBROOT=${KS_WORK}/ks_build_sources/chm13v2_ksweep/db

T=${SLURM_ARRAY_TASK_ID:-1}
KS_LIST=(21 31 41 51 61)
if [ "$T" -le 5 ]; then
  ASM=chm13; FA=${KS_REFERENCES}/CHM13v2/chm13v2.0.fasta; K=${KS_LIST[$((T-1))]}
else
  ASM=hg002hap1; FA=${KS_REFERENCES}/HG002/hg002v1.1_hap1.fasta.gz; K=${KS_LIST[$((T-6))]}
fi
DB=HKS_human_CHM13_v2_k${K}
OUT=$W/annot_${ASM}_k${K}; RES=$W/results
mkdir -p "$RES"; rm -rf "$OUT"; mkdir -p "$OUT"
PROV=$RES/prov_${ASM}_k${K}.txt
{
  echo "=== A2 step2 ${ASM} k=${K} ==="
  echo "date: $(date -Is)  job: ${SLURM_JOB_ID:-NA}[$T]  node: $(hostname)"
  lscpu | grep -E 'Model name|Socket|Core\(s\) per socket'
  echo "governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null)"
  echo "karyoscope: $($KS --version)  hks: $($KARYOSCOPE_HKS 2>&1 | grep -m1 -o 'version [0-9.]*')"
  echo "db: $DB   db_root: $DBROOT   threads: 16"
} | tee "$PROV"

echo "[$(date +%T)] annotate"
/usr/bin/time -v $KS -v annotate -i "$FA" -o "$OUT" --db "$DB" --db-root "$DBROOT" \
    -t 16 --no-bgzip > "$RES/annotate_${ASM}_k${K}.log" 2> "$RES/time_${ASM}_k${K}.txt"
RC=$?
WALL=$(grep 'Elapsed (wall clock)' "$RES/time_${ASM}_k${K}.txt" | awk '{print $NF}')
RSS=$(grep 'Maximum resident set size' "$RES/time_${ASM}_k${K}.txt" | awk '{print $NF}')
echo "annotate_rc=$RC wall=$WALL rss_kb=$RSS" | tee -a "$PROV"
printf "%s\t%s\t%s\t%s\t%s\n" "$ASM" "$K" "$RC" "$WALL" "$RSS" >> "$RES/step2_runtimes.tsv"
[ $RC -eq 0 ] || exit 1

# ---- per-base leaf / internal-node / novel, per feature set
echo "[$(date +%T)] label-level accounting"
python3 - "$OUT" "$W" "$ASM" "$K" "$RES" <<'PY' | tee -a "$PROV"
import glob, gzip, os, sys, collections
OUT, W, ASM, K, RES = sys.argv[1:6]
HIER=os.environ["KS_SOFTWARE"]+"/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv"
# leaves are per feature set: nodes that are never a parent within that set
par=collections.defaultdict(set); nodes=collections.defaultdict(set)
for line in open(HIER):
    p=line.rstrip("\n").split("\t")
    if len(p)<3 or p[0]=="feature_set": continue
    fs,c,pa=p[0],p[1],p[2]
    nodes[fs].add(c); nodes[fs].add(pa); par[fs].add(pa)
leaves={fs:nodes[fs]-par[fs] for fs in nodes}

rows=[]
for stage in ("presmoothed","smoothed"):
    for f in sorted(glob.glob(os.path.join(OUT, f"*.{stage}.bed*"))):
        b=os.path.basename(f)
        if ".binned" in b or ".scaffolded" in b: continue
        parts=b.split(".")
        fs=parts[parts.index(stage)-1]
        cnt=collections.Counter()
        op=gzip.open if f.endswith(".gz") else open
        with op(f,"rt") as fh:
            for line in fh:
                p=line.rstrip("\n").split("\t")
                if len(p)<4: continue
                try: L=int(p[2])-int(p[1])
                except ValueError: continue
                lab=p[3]
                if lab.lower()=="novel": lvl="novel"
                elif lab in leaves.get(fs,()): lvl="leaf"
                elif lab in nodes.get(fs,()): lvl="internal"
                else: lvl="other"
                cnt[lvl]+=L
        tot=sum(cnt.values())
        for lvl in ("leaf","internal","novel","other"):
            rows.append((ASM,K,fs,stage,lvl,cnt[lvl],100.0*cnt[lvl]/tot if tot else 0.0))
out=os.path.join(RES,f"labellevel_{ASM}_k{K}.tsv")
with open(out,"w") as o:
    o.write("assembly\tk\tfeature_set\tstage\tlevel\tbases\tfraction_pct\n")
    for r in rows: o.write("\t".join(str(x) for x in r)+"\n")
print(f"{'set':14}{'stage':13}{'leaf%':>8}{'internal%':>11}{'novel%':>9}")
seen=set()
for a,k,fs,stage,lvl,b,frac in rows:
    key=(fs,stage)
    if key in seen: continue
    seen.add(key)
    d={l:f for (_,_,f2,s2,l,_,f) in rows if f2==fs and s2==stage for l,f in [(l,f)]}
    print(f"{fs:14}{stage:13}{d.get('leaf',0):>8.2f}{d.get('internal',0):>11.2f}{d.get('novel',0):>9.3f}")
print(f"wrote {out}")
PY
echo "date_end: $(date -Is)" >> "$PROV"
echo "=== DONE ${ASM} k=${K} ==="
