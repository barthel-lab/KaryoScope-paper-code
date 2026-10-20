#!/bin/bash
# A4 follow-up: smoothing composition per feature set on ONE denominator.
#
# Answers a4-followup-novel-imputation-ask.md. A4.2 (run_a4_analysis.sh) counted
# every evaluated base; A4.3 (run_a4_novel_v3.sh) masked N runs and chrM. Mixing
# their percentages is why the gene rows did not sum to 100. This redoes the
# A4.2 level accounting with A4.3's exact mask so all three categories share the
# 2,940,260,101 denominator and sum to 100% per set per stage.
#
# Mask is copied verbatim from run_a4_novel_v3.sh: N intervals extended 30 bp
# left (fill_gaps' k-1 for k=31), plus chrM.
#SBATCH --job-name=a4_comp
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-6%6
#SBATCH --output=a4comp_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing; RES=$W/results/composition; mkdir -p "$RES"
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
HIER=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv
SETS=(chromosome region repeat subtelomeric gene acrocentric)
T=${SLURM_ARRAY_TASK_ID}
ASM=hg002hap1; SET=${SETS[$((T-1))]}; DIR=$W/annot_hg002hap1
NBED=$W/HG002.hap1.n_regions.renamed.bed; EXCL=chrM
# .${SET}. anchor keeps this off presmoothed -- do NOT relax it to *smoothed*
PRE=$(ls $DIR/*.${SET}.presmoothed.bed 2>/dev/null | head -1)
POST=$(ls $DIR/*.${SET}.smoothed.bed 2>/dev/null | head -1)
echo "task=$T asm=$ASM set=$SET node=$(hostname)"
echo "  PRE =$PRE"
echo "  POST=$POST"
[ -n "$PRE" ] && [ -n "$POST" ] || { echo "MISSING inputs in $DIR"; exit 1; }
case "$PRE" in *presmoothed*) ;; *) echo "PRE glob wrong"; exit 1;; esac
case "$POST" in *presmoothed*) echo "POST matched presmoothed"; exit 1;; esac

TMP=$(mktemp -d $RES/tmp.${ASM}_${SET}.XXXX); trap 'rm -rf $TMP' EXIT
echo "[$(date +%T)] unionbedg pre vs post"
$BT unionbedg -i "$PRE" "$POST" > $TMP/join.bed
echo "  rows: $(wc -l < $TMP/join.bed)"

python3 - "$TMP/join.bed" "$HIER" "$SET" "$ASM" "$RES" "$NBED" "$EXCL" <<'PY'
import sys, collections
join, hier, FS, ASM, RES, nbed, excl = sys.argv[1:8]
EXCL=set(x for x in excl.split(",") if x)

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

# N intervals, extended 30 bp left exactly as run_a4_novel_v3.sh does
N=collections.defaultdict(list)
for line in open(nbed):
    c,s,e=line.split("\t")[:3]
    N[c].append((max(0,int(s)-30), int(e)))
def n_overlap(c,s,e):
    tot=0
    for a,b in N.get(c,()):
        lo,hi=max(s,a),min(e,b)
        if hi>lo: tot+=hi-lo
    return tot

pre_tot=collections.Counter(); post_tot=collections.Counter(); masked=0
for line in open(join):
    c,s,e,a,b=line.rstrip("\n").split("\t")
    if a=="0" or b=="0": continue          # unionbedg filler
    if c in EXCL: continue
    s,e=int(s),int(e)
    nl=n_overlap(c,s,e); masked+=nl
    eff=(e-s)-nl
    if eff<=0: continue
    pre_tot[lvl(a)]+=eff; post_tot[lvl(b)]+=eff

den=sum(pre_tot.values()); assert den==sum(post_tot.values())
out=f"{RES}/composition_{ASM}_{FS}.tsv"
with open(out,"w") as o:
    o.write("assembly\tfeature_set\tstage\tlabel_specific_bp\tmultilabel_bp\tnovel_bp\t"
            "label_specific_pct\tmultilabel_pct\tnovel_pct\tdenominator_bp\tsum_pct\n")
    for stage,d in (("presmoothed",pre_tot),("smoothed",post_tot)):
        assert d["other"]==0, f"unexpected 'other' bases: {d['other']}"
        L,M,V=d["leaf"],d["internal"],d["novel"]
        pl,pm,pv=(100*L/den,100*M/den,100*V/den)
        o.write(f"{ASM}\t{FS}\t{stage}\t{L}\t{M}\t{V}\t{pl:.4f}\t{pm:.4f}\t{pv:.4f}\t{den}\t{pl+pm+pv:.4f}\n")
print(f"=== {ASM} / {FS}   denominator {den:,}   N+chrM masked {masked:,}")
for stage,d in (("presmoothed",pre_tot),("smoothed",post_tot)):
    L,M,V=d["leaf"],d["internal"],d["novel"]
    print(f"  {stage:12} label_specific {100*L/den:8.4f}  multilabel {100*M/den:7.4f}  novel {100*V/den:7.4f}  sum {100*(L+M+V)/den:8.4f}")
PY
echo "[$(date +%T)] done"
