#!/bin/bash
# HG002 hap2 reproduction check. The manuscript states, for hap2:
#   evaluated positions 3,028,094,375 ; smoothed exact 95.21% ; compatible 1.08%
# and the two-haplotype ranges: pre exact 84.6-84.8, pre compatible 9.2-9.5,
# pre novel 2.5-2.6, smoothed novel <0.05.
# hap1 gave the UPPER end of pre-exact (84.7787), so hap2 should give ~84.6.
#SBATCH --job-name=hap2_check
#SBATCH --partition=compute
#SBATCH --constraint=2.80GHz&Ice_Lake
#SBATCH --exclude=g-c-1-4-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=hap2_check_%j.log
: "${KS_REFERENCES:?set KS_REFERENCES (see the repository README, Paths)}"
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
: "${KS_WORK:?set KS_WORK (see the repository README, Paths)}"
set -uo pipefail
A=${KS_SOFTWARE}/github/barthel/KaryoScope-archive
W=${KS_WORK}/a4_smoothing; RES=$W/results
KSBIN=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin
export KARYOSCOPE_HKS=${KS_SOFTWARE}/github/barthel/HKS/target/release/hks
export PATH=$KSBIN:$PATH
ROOT=${KS_SOFTWARE}/github/barthel/KaryoScope/resources/registry
FA=${KS_REFERENCES}/HG002/hg002v1.1_hap2.fasta.gz
OUT=$W/annot_hg002hap2; rm -rf "$OUT"; mkdir -p "$OUT"

echo "node=$(hostname) start=$(date -Is)"
# rename BOTH truth and gaps: _MATERNAL -> _hap2 (the gaps rename is the step
# whose omission cost 7.6 Mb on hap1)
sed 's/_MATERNAL\t/_hap2\t/' $A/resources/KS_human_HG002/HG002.hap2.RepeatMasker.processed.sorted.merged.bed > $W/HG002.hap2.RepeatMasker.renamed.bed
sed 's/_MATERNAL\t/_hap2\t/' $A/resources/KS_human_HG002/HG002.hap2.n_regions.bed > $W/HG002.hap2.n_regions.renamed.bed
echo "truth contigs: $(cut -f1 $W/HG002.hap2.RepeatMasker.renamed.bed | head -1)  gaps contigs: $(cut -f1 $W/HG002.hap2.n_regions.renamed.bed | head -1)"
samtools faidx "$FA" 2>/dev/null

echo "[$(date +%T)] annotate hap2 (repeat set only)"
/usr/bin/time -v $KSBIN/karyoscope -v annotate -i "$FA" -o "$OUT" --db HKS_human_CHM13_v2 \
    --db-root "$ROOT" --feature-set repeat -t 16 --no-bgzip \
    > "$RES/annotate_hap2.log" 2> "$RES/time_hap2.txt"
echo "annotate_rc=$? wall=$(grep 'Elapsed (wall' "$RES/time_hap2.txt" | awk '{print $NF}')"

echo "[$(date +%T)] concordance"
bash ${KS_WORK}/ks_concord.sh \
  --truth $W/HG002.hap2.RepeatMasker.renamed.bed \
  --pre  "$(ls $OUT/*.repeat.presmoothed.bed | head -1)" \
  --post "$(ls $OUT/*.repeat.smoothed.bed | head -1)" \
  --hierarchy $ROOT/HKS_human_CHM13_v2/hierarchy.tsv \
  --feature-set repeat \
  --fai ${FA}.fai \
  --gaps $W/HG002.hap2.n_regions.renamed.bed \
  --kmer-size 31 \
  --merge "Simple_repeat,Satellite,Low_complexity,Noninterspersed:tandem_repeat" \
  --out $RES/hap2_repeat
rc=$?; echo "concord_rc=$rc"
[ -s $RES/hap2_repeat.summary.tsv ] || { echo "FATAL: no summary"; exit 3; }

python3 - $RES/hap2_repeat.summary.tsv <<'PY'
import sys
rows={}
for i,l in enumerate(open(sys.argv[1])):
    if i==0: continue
    st,c,b,f=l.split("\t"); rows[(st,c)]=(int(b),float(f))
tot=sum(rows[("presmoothed",c)][0] for c in ("equivalent","compatible","novel","incompatible"))
print("\n=========== HG002 hap2 vs MANUSCRIPT ===========")
print(f"  positions evaluated : {tot:,}   manuscript: 3,028,094,375   {'EXACT MATCH' if tot==3028094375 else '*** MISMATCH: delta %+d ***'%(tot-3028094375)}")
def chk(st,c,lo,hi,label):
    v=rows[(st,c)][1]
    ok = lo<=v<=hi
    print(f"  {label:22}: {v:8.4f}%   manuscript: {lo}-{hi}%   {'IN RANGE' if ok else '*** OUT OF RANGE ***'}")
chk("presmoothed","equivalent",84.6,84.8,"pre exact")
chk("presmoothed","compatible",9.2,9.5,"pre compatible")
chk("presmoothed","novel",2.5,2.6,"pre novel")
v=rows[("smoothed","equivalent")][1]; print(f"  {'smoothed exact':22}: {v:8.4f}%   manuscript: 95.21%   rounds to {round(v,2)}  {'MATCH' if round(v,2)==95.21 else '*** MISMATCH ***'}")
v=rows[("smoothed","compatible")][1]; print(f"  {'smoothed compatible':22}: {v:8.4f}%   manuscript: 1.08%    rounds to {round(v,2)}  {'MATCH' if round(v,2)==1.08 else '*** MISMATCH ***'}")
v=rows[("smoothed","novel")][1];      print(f"  {'smoothed novel':22}: {v:8.4f}%   manuscript: <0.05%   {'SATISFIES' if v<0.05 else '*** EXCEEDS ***'}")
print("================================================")
PY
