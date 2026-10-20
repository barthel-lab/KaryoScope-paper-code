#!/bin/bash
# A4.2 + A4.3: smoothing characterization from the annotations already on disk.
#   A4.2  per-base {leaf, internal-node, novel} fractions pre- and post-smoothing,
#         per feature set, PLUS the movement matrix (pre level -> post level).
#   A4.3  persistent vs transient novel: pre-smoothing novel intervals classified
#         by post-smoothing fate, with counts, bp and length medians.
# Task 13 is the single-feature-set KMC cross-check A4.2 requires before
# reporting HKS throughout.
#
# No sorting anywhere: annotations tile the input and are already in genome
# order (invariant). unionbedg needs both inputs in the SAME order, which they
# are, coming from one annotate run.
#SBATCH --job-name=a4_analysis
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-13%6
#SBATCH --output=a4an_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing; RES=$W/results/smoothing; mkdir -p "$RES"
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
HIER=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv
SETS=(chromosome region repeat subtelomeric gene acrocentric)
T=${SLURM_ARRAY_TASK_ID}
if [ "$T" -eq 13 ]; then
  ASM=chm13_kmc; SET=repeat; DIR=$W/annot_chm13_kmc_repeat
elif [ "$T" -le 6 ]; then
  ASM=chm13; SET=${SETS[$((T-1))]}; DIR=$W/annot_chm13
else
  ASM=hg002hap1; SET=${SETS[$((T-7))]}; DIR=$W/annot_hg002hap1
fi
PRE=$(ls $DIR/*.${SET}.presmoothed.bed 2>/dev/null | head -1)
POST=$(ls $DIR/*.${SET}.smoothed.bed 2>/dev/null | head -1)
echo "task=$T asm=$ASM set=$SET node=$(hostname)"
[ -n "$PRE" ] && [ -n "$POST" ] || { echo "MISSING inputs in $DIR"; exit 1; }

TMP=$(mktemp -d $RES/tmp.${ASM}_${SET}.XXXX); trap 'rm -rf $TMP' EXIT
echo "[$(date +%T)] unionbedg pre vs post"
$BT unionbedg -i "$PRE" "$POST" > $TMP/join.bed
echo "  rows: $(wc -l < $TMP/join.bed)"

python3 - "$TMP/join.bed" "$HIER" "$SET" "$ASM" "$RES" "$PRE" <<'PY'
import sys, collections, statistics
join, hier, FS, ASM, RES, PREBED = sys.argv[1:7]
par=collections.defaultdict(set); nodes=collections.defaultdict(set)
for line in open(hier):
    p=line.rstrip("\n").split("\t")
    if len(p)<3 or p[0]=="feature_set": continue
    nodes[p[0]].add(p[1]); nodes[p[0]].add(p[2]); par[p[0]].add(p[2])
leaves=nodes[FS]-par[FS]
def lvl(lab):
    if lab.lower()=="novel": return "novel"
    if lab in leaves: return "leaf"
    if lab in nodes[FS]: return "internal"
    return "other"

move=collections.Counter(); pre_tot=collections.Counter(); post_tot=collections.Counter()
novel_fate=collections.Counter()          # bases
run_len=[]; run_novel_post=[]             # per-interval persistence
cur=None; cur_len=0; cur_still=0
for line in open(join):
    c,s,e,a,b=line.rstrip("\n").split("\t")
    if a=="0" or b=="0": continue         # unionbedg filler
    L=int(e)-int(s)
    la,lb=lvl(a),lvl(b)
    move[(la,lb)]+=L; pre_tot[la]+=L; post_tot[lb]+=L
    if la=="novel":
        novel_fate["persistent" if lb=="novel" else "transient"]+=L

# interval-level novel fate, from the ORIGINAL presmoothed intervals
post_by_pos={}
pt=sum(pre_tot.values())
with open(RES+f"/smoothing_accounting_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tstage\tlevel\tbases\tfraction_pct\n")
    for stage,d in (("presmoothed",pre_tot),("smoothed",post_tot)):
        t=sum(d.values())
        for l in ("leaf","internal","novel","other"):
            o.write(f"{ASM}\t{FS}\t{stage}\t{l}\t{d[l]}\t{100*d[l]/t if t else 0:.4f}\n")
with open(RES+f"/movement_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tpre_level\tpost_level\tbases\tfraction_of_pre_pct\n")
    for (a,b),v in sorted(move.items(), key=lambda x:-x[1]):
        o.write(f"{ASM}\t{FS}\t{a}\t{b}\t{v}\t{100*v/pre_tot[a] if pre_tot[a] else 0:.4f}\n")
nf=sum(novel_fate.values())
with open(RES+f"/novel_fate_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tclass\tbases\tfraction_of_novel_pct\n")
    for k in ("transient","persistent"):
        o.write(f"{ASM}\t{FS}\t{k}\t{novel_fate[k]}\t{100*novel_fate[k]/nf if nf else 0:.4f}\n")

print(f"\n=== {ASM} / {FS} ===")
print(f"{'stage':13}{'leaf%':>9}{'internal%':>11}{'novel%':>9}")
for stage,d in (("presmoothed",pre_tot),("smoothed",post_tot)):
    t=sum(d.values())
    print(f"{stage:13}{100*d['leaf']/t:>9.3f}{100*d['internal']/t:>11.3f}{100*d['novel']/t:>9.4f}")
print(f"novel fate: transient {novel_fate['transient']:,} bp ({100*novel_fate['transient']/nf if nf else 0:.2f}%), "
      f"persistent {novel_fate['persistent']:,} bp ({100*novel_fate['persistent']/nf if nf else 0:.2f}%)")
print("top movements (pre -> post, bases):")
for (a,b),v in sorted(move.items(), key=lambda x:-x[1])[:6]:
    print(f"   {a:9} -> {b:9} {v:>14,}")
PY
echo "[$(date +%T)] done"
