#!/bin/bash
# A4.3 v2 — novel sequence characterisation, corrected framing.
#
# `novel` is the one label that is not a claim about what the sequence IS: the
# k-mers are absent from the database. Smoothing does NOT remove novelty, it
# IMPUTES a label from flanking context. So:
#   * pre-smoothing novel  = the measure of novel sequence content
#   * post-smoothing novel = novel sequence smoothing could not even label
#   * the difference       = novel sequence with an IMPUTED label
#
# Two corrections over v1:
#   1. N runs are MASKED OUT. A k-mer spanning N matches nothing and reports as
#      novel, but that is ABSENCE, not novelty - categorically different. Ns are
#      9.2% of HG002 hap1's raw novel (7,634,373 of 83,031,126 bp). The
#      manuscript excludes them too (hence its 2,940,260,101 evaluated
#      positions). Masking uses the same 30 bp left extension fill_gaps applies.
#   2. RUN-LENGTH DISTRIBUTIONS, which v1 declared and never computed. Run
#      length is the discriminator: a substitution breaks ~k=31 bp of k-mers, so
#      SHORT novel runs are divergence artifacts, while LONG runs are genuinely
#      absent sequence (insertions, divergent haplotypes).
#SBATCH --job-name=a4_novel2
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-12%6
#SBATCH --output=a4nv_%A_%a.log
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
W=${KS_WORK}/a4_smoothing; RES=$W/results/novel_v2; mkdir -p "$RES"
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
SETS=(chromosome region repeat subtelomeric gene acrocentric)
T=${SLURM_ARRAY_TASK_ID}
if [ "$T" -le 6 ]; then ASM=chm13; SET=${SETS[$((T-1))]}; DIR=$W/annot_chm13; NBED=""
else ASM=hg002hap1; SET=${SETS[$((T-7))]}; DIR=$W/annot_hg002hap1; NBED=$W/HG002.hap1.n_regions.renamed.bed; fi
PRE=$(ls $DIR/*.${SET}.presmoothed.bed 2>/dev/null | head -1)
POST=$(ls $DIR/*.${SET}.smoothed.bed 2>/dev/null | head -1)
echo "task=$T asm=$ASM set=$SET nmask='${NBED:-none (gapless)}' node=$(hostname)"
[ -n "$PRE" ] && [ -n "$POST" ] || { echo "MISSING inputs"; exit 1; }
TMP=$(mktemp -d $RES/tmp.${ASM}_${SET}.XXXX); trap 'rm -rf $TMP' EXIT
$BT unionbedg -i "$PRE" "$POST" > $TMP/join.bed

python3 - "$TMP/join.bed" "${NBED:-}" "$ASM" "$SET" "$RES" <<'PY'
import sys, statistics, collections
join, nbed, ASM, FS, RES = sys.argv[1:6]
# N intervals, extended 30 bp left exactly as fill_gaps does for k=31
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

tot_bases=0; novel_pre=0; novel_post=0
runs=[]                       # (length, bases_still_novel)
cur_c=None; cur_end=None; cur_len=0; cur_still=0
def close():
    global cur_len,cur_still
    if cur_len>0: runs.append((cur_len,cur_still))
    cur_len=0; cur_still=0
for line in open(join):
    c,s,e,a,b=line.rstrip("\n").split("\t")
    if a=="0" or b=="0": continue
    s,e=int(s),int(e); L=e-s
    nl=n_overlap(c,s,e)
    eff=L-nl                                   # N bases masked out entirely
    if eff<=0:
        close(); cur_c=None; continue
    tot_bases+=eff
    if a.lower()=="novel":
        novel_pre+=eff
        still = eff if b.lower()=="novel" else 0
        novel_post+=still
        if cur_c==c and cur_end==s: cur_len+=eff; cur_still+=still
        else: close(); cur_len=eff; cur_still=still
        cur_c=c; cur_end=e
    else:
        close(); cur_c=None
close()

imputed=[l for l,st in runs if st==0]
unlabeled=[l for l,st in runs if st==l]
partial=[l for l,st in runs if 0<st<l]
def med(x): return statistics.median(x) if x else 0
def tot(x): return sum(x)
out=f"{RES}/novel_v2_{ASM}_{FS}.tsv"
with open(out,"w") as o:
    o.write("assembly\tfeature_set\tmetric\tvalue\n")
    for k,v in (("evaluated_bases_N_excluded",tot_bases),
                ("novel_pre_bases",novel_pre),
                ("novel_pre_pct",round(100*novel_pre/tot_bases,4) if tot_bases else 0),
                ("novel_post_bases",novel_post),
                ("novel_post_pct",round(100*novel_post/tot_bases,4) if tot_bases else 0),
                ("runs_total",len(runs)),
                ("runs_fully_imputed",len(imputed)),("bp_fully_imputed",tot(imputed)),("median_len_imputed",med(imputed)),
                ("runs_fully_unlabeled",len(unlabeled)),("bp_fully_unlabeled",tot(unlabeled)),("median_len_unlabeled",med(unlabeled)),
                ("runs_partial",len(partial)),("bp_partial",tot(partial)),("median_len_partial",med(partial)),
                ("max_run_len",max((l for l,_ in runs), default=0))):
        o.write(f"{ASM}\t{FS}\t{k}\t{v}\n")
# length histogram, log-ish buckets: run length vs k=31
buck=collections.Counter()
for l,_ in runs:
    b = "1-10" if l<=10 else "11-31" if l<=31 else "32-100" if l<=100 else "101-1k" if l<=1000 else "1k-10k" if l<=10000 else ">10k"
    buck[b]+=1
with open(f"{RES}/novel_runlen_hist_{ASM}_{FS}.tsv","w") as o:
    o.write("assembly\tfeature_set\tbucket\tn_runs\n")
    for b in ("1-10","11-31","32-100","101-1k","1k-10k",">10k"):
        o.write(f"{ASM}\t{FS}\t{b}\t{buck[b]}\n")

print(f"\n=== {ASM} / {FS}  (N-excluded) ===")
print(f"  evaluated bases      : {tot_bases:,}")
print(f"  novel PRE  (novel sequence content) : {novel_pre:,} bp = {100*novel_pre/tot_bases if tot_bases else 0:.4f}%")
print(f"  novel POST (left unlabeled)         : {novel_post:,} bp = {100*novel_post/tot_bases if tot_bases else 0:.4f}%")
print(f"  runs: {len(runs):,}  imputed {len(imputed):,} (median {med(imputed)} bp)  "
      f"unlabeled {len(unlabeled):,} (median {med(unlabeled)} bp)  partial {len(partial):,}")
print(f"  longest novel run: {max((l for l,_ in runs), default=0):,} bp")
print("  run-length histogram (k=31 is the substitution-artifact scale):")
for b in ("1-10","11-31","32-100","101-1k","1k-10k",">10k"):
    print(f"     {b:8} {buck[b]:>10,}")
PY
