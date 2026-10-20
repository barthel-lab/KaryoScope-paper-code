#!/bin/bash
# A6 §3 item 1 — MULTILABEL (ambiguous) k-mer run-length distribution.
#
# The gap Rhyker flagged: main.tex:366 says novel k-mers "are treated uniformly
# with multilabel k-mers during smoothing", so max_gap bridges gaps arising from
# BOTH populations -- but only the novel one has ever been characterised (A4).
# This is the matching measurement for the multilabel population.
#
# "Multilabel" = the presmoothed label is an INTERNAL node of the feature set's
# hierarchy. A k-mer matching several leaves resolves to their LCA, so an
# internal-node call IS the ambiguous call. Leaves are children that never
# appear as a parent within that feature set. `novel` is excluded -- it is the
# other population and A4 already measured it.
#
# Deliberately mirrors run_a4_novel_v3.sh: same unionbedg join, same N masking
# (30 bp left extension, matching fill_gaps at k=31), same chrM exclusion, and
# the same run-length buckets, so the two tables are directly comparable.
#SBATCH --job-name=a6_multilabel
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-12%6
#SBATCH --output=a6ml_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing
RES=${KS_WORK}/a6_sweeps/results/multilabel; mkdir -p "$RES"
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
HIER=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry/HKS_human_CHM13_v2/hierarchy.tsv
SETS=(chromosome region repeat subtelomeric gene acrocentric)
T=${SLURM_ARRAY_TASK_ID}
if [ "$T" -le 6 ]; then ASM=chm13; SET=${SETS[$((T-1))]}; DIR=$W/annot_chm13; NBED=""; EXCL=chrM
else ASM=hg002hap1; SET=${SETS[$((T-7))]}; DIR=$W/annot_hg002hap1; NBED=$W/HG002.hap1.n_regions.renamed.bed; EXCL=chrM; fi
export EXCL
PRE=$(ls $DIR/*.${SET}.presmoothed.bed 2>/dev/null | head -1)
POST=$(ls $DIR/*.${SET}.smoothed.bed 2>/dev/null | head -1)
echo "task=$T asm=$ASM set=$SET node=$(hostname) start=$(date -Is)"
[ -n "$PRE" ] && [ -n "$POST" ] || { echo "MISSING inputs"; exit 1; }
TMP=$(mktemp -d $RES/tmp.${ASM}_${SET}.XXXX); trap 'rm -rf $TMP' EXIT
$BT unionbedg -i "$PRE" "$POST" > $TMP/join.bed

python3 - "$TMP/join.bed" "${NBED:-}" "$ASM" "$SET" "$RES" "${EXCL:-}" "$HIER" <<'PY'
import sys, statistics, collections
join, nbed, ASM, FS, RES = sys.argv[1:6]
EXCL=set(x for x in (sys.argv[6] or "").split(",") if x)
HIER=sys.argv[7]

# Leaves of THIS feature set: children that are never a parent.
children, parents = set(), set()
with open(HIER) as f:
    hdr = f.readline()
    for line in f:
        fs, child, parent = line.rstrip("\n").split("\t")[:3]
        if fs != FS: continue
        children.add(child); parents.add(parent)
LEAVES = children - parents
def is_multilabel(lab):
    l = lab.lower()
    if l in ("novel", "unknown", ""): return False
    return lab not in LEAVES
print(f"  hierarchy: {len(children)} nodes, {len(LEAVES)} leaves for set '{FS}'")

N=collections.defaultdict(list)
if nbed:
    for line in open(nbed):
        c,s,e=line.split("\t")[:3]
        N[c].append((max(0,int(s)-30), int(e)))
def n_overlap(c,s,e):
    tot=0
    for a,b in N.get(c,()):
        lo,hi=max(s,a),min(e,b)
        if hi>lo: tot+=hi-lo
    return tot

tot_bases=0; ml_pre=0; ml_post=0
runs=[]
cur_c=None; cur_end=None; cur_len=0; cur_still=0
def close():
    global cur_len,cur_still
    if cur_len>0: runs.append((cur_len,cur_still))
    cur_len=0; cur_still=0
for line in open(join):
    c,s,e,a,b=line.rstrip("\n").split("\t")
    if a=="0" or b=="0": continue
    if c in EXCL: close(); cur_c=None; continue
    s,e=int(s),int(e); L=e-s
    eff=L-n_overlap(c,s,e)
    if eff<=0:
        close(); cur_c=None; continue
    tot_bases+=eff
    if is_multilabel(a):
        ml_pre+=eff
        # "still ambiguous" = smoothing did NOT resolve it to a leaf
        still = eff if is_multilabel(b) else 0
        ml_post+=still
        if cur_c==c and cur_end==s: cur_len+=eff; cur_still+=still
        else: close(); cur_len=eff; cur_still=still
        cur_c=c; cur_end=e
    else:
        close(); cur_c=None
close()

resolved=[l for l,st in runs if st==0]
unresolved=[l for l,st in runs if st==l]
partial=[l for l,st in runs if 0<st<l]
def med(x): return statistics.median(x) if x else 0
with open(f"{RES}/multilabel_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tmetric\tvalue\n")
    for k,v in (("evaluated_bases_N_excluded",tot_bases),
                ("multilabel_pre_bases",ml_pre),
                ("multilabel_pre_pct",round(100*ml_pre/tot_bases,4) if tot_bases else 0),
                ("multilabel_post_bases",ml_post),
                ("multilabel_post_pct",round(100*ml_post/tot_bases,4) if tot_bases else 0),
                ("runs_total",len(runs)),
                ("runs_resolved_to_leaf",len(resolved)),("bp_resolved",sum(resolved)),("median_len_resolved",med(resolved)),
                ("runs_unresolved",len(unresolved)),("bp_unresolved",sum(unresolved)),("median_len_unresolved",med(unresolved)),
                ("runs_partial",len(partial)),("bp_partial",sum(partial)),("median_len_partial",med(partial)),
                ("max_run_len",max((l for l,_ in runs), default=0))):
        o.write(f"{ASM}\t{FS}\t{k}\t{v}\n")

buck=collections.Counter(); buckbp=collections.Counter()
for l,_ in runs:
    b = "1-10" if l<=10 else "11-31" if l<=31 else "32-100" if l<=100 else "101-1k" if l<=1000 else "1k-10k" if l<=10000 else ">10k"
    buck[b]+=1; buckbp[b]+=l
cum=0
with open(f"{RES}/multilabel_runlen_hist_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tbucket\tn_runs\tbases\tpct_of_multilabel_bp\tcum_pct\n")
    for b in ("1-10","11-31","32-100","101-1k","1k-10k",">10k"):
        p = 100*buckbp[b]/ml_pre if ml_pre else 0; cum += p
        o.write(f"{ASM}\t{FS}\t{b}\t{buck[b]}\t{buckbp[b]}\t{p:.4f}\t{cum:.4f}\n")

print(f"\n=== MULTILABEL {ASM} / {FS}  (N-excluded) ===")
print(f"  evaluated bases          : {tot_bases:,}")
print(f"  multilabel PRE           : {ml_pre:,} bp = {100*ml_pre/tot_bases if tot_bases else 0:.4f}%")
print(f"  multilabel POST (unresolved): {ml_post:,} bp = {100*ml_post/tot_bases if tot_bases else 0:.4f}%")
print(f"  runs: {len(runs):,}  resolved {len(resolved):,} (median {med(resolved)} bp)  "
      f"unresolved {len(unresolved):,} (median {med(unresolved)} bp)  partial {len(partial):,}")
print(f"  longest multilabel run   : {max((l for l,_ in runs), default=0):,} bp")
print("  run-length histogram (the distribution max_gap=1000 has to bridge):")
cum=0
for b in ("1-10","11-31","32-100","101-1k","1k-10k",">10k"):
    p = 100*buckbp[b]/ml_pre if ml_pre else 0; cum += p
    print(f"     {b:8} {buck[b]:>10,} runs  {buckbp[b]:>14,} bp  {p:6.2f}%  (cum {cum:6.2f}%)")
PY
echo "done $(date -Is)"
