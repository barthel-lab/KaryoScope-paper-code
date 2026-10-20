#!/bin/bash
# ks_concord.sh — per-base concordance of a KaryoScope annotation against a truth
# BED, pre- AND post-smoothing, using the SAME logic that produced the paper's
# Figure 2 / Supp Note 1 numbers.
#
# It is a thin wrapper on the published chain, deliberately adding no logic of
# its own:
#     fill_gaps.py  ->  bedtools unionbedg  ->  contrast_HG002.py
#
# Why each step matters (all four are required for the numbers to be comparable
# to the published ones):
#   1. fill_gaps.py --kmer-size K trims each contig's effective end by K-1,
#      because KaryoScope annotates n-K+1 positions, not n. It also labels
#      N-runs 'exclude'. Without it the last K-1 bp of every contig and every
#      N-run score as disagreement.
#   2. unionbedg puts truth/presmoothed/smoothed on one coordinate grid; its
#      default "0" filler marks bases missing from some input, and
#      contrast_HG002.py skips those rows.
#   3. contrast_HG002.py is hierarchy-aware: 'compatible' means the prediction
#      is an ANCESTOR of the truth (correct but less specific). Raw per-base
#      equality would score those as errors and cannot reproduce the anchors.
#   4. --merge is feature-set specific and load-bearing: the published repeat
#      numbers use the tandem_repeat merge; region uses none.
#
# Output: <out>.contrast.tsv  (stage, truth, prediction, length, compatibility)
#     and <out>.summary.tsv   (stage, compatibility class, bases, fraction)
: "${KS_SOFTWARE:?set KS_SOFTWARE (see the repository README, Paths)}"
set -euo pipefail

ARCH=${KS_SOFTWARE}/github/barthel/KaryoScope-archive/workflow/scripts
BT=${KS_SOFTWARE}/miniforge3_new/envs/karyoscope/bin/bedtools
KMER=31; MERGE=""; GAPS=""; KEEP=0; NOFILL=0
usage () { sed -n '2,30p' "$0"; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --truth)       TRUTH=$2; shift 2;;
    --pre)         PRE=$2; shift 2;;
    --post)        POST=$2; shift 2;;
    --hierarchy)   HIER=$2; shift 2;;
    --feature-set) FS=$2; shift 2;;
    --fai)         FAI=$2; shift 2;;
    --gaps)        GAPS=$2; shift 2;;
    --kmer-size)   KMER=$2; shift 2;;
    --merge)       MERGE=$2; shift 2;;
    --out)         OUT=$2; shift 2;;
    --keep-temp)   KEEP=1; shift;;
    --no-fill)     NOFILL=1; shift;;
    -h|--help)     usage;;
    *) echo "unknown arg: $1" >&2; usage;;
  esac
done
: "${TRUTH:?--truth required}" "${PRE:?--pre required}" "${POST:?--post required}"
: "${HIER:?--hierarchy required}" "${FS:?--feature-set required}"
: "${FAI:?--fai required}" "${OUT:?--out required}"

TMP=$(mktemp -d "${OUT}.tmp.XXXXXX")
cleanup () { [ "$KEEP" -eq 1 ] || rm -rf "$TMP"; }
trap cleanup EXIT

decompress () { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

echo "[concord] feature_set=$FS kmer=$KMER merge='${MERGE:-none}'" >&2

# 0. CONTIG-NAME ASSERTION.
# Learned the hard way (job 35378922): the HG002 truth BED uses chr1_PATERNAL
# while the annotated assembly uses chr1_hap1. fill_gaps then filled EVERY base
# with 'nonrepeat', the truth intervals were dropped as unionbedg filler, and
# the run exited 0 with plausible-looking percentages (88% incompatible,
# compatible exactly 0.0000). Only comparison against a published number caught
# it. Never let that pass silently again.
contigs () { decompress "$1" | awk '{if(!($1 in s)){s[$1]; print $1}}' | sort; }
echo "[concord] checking contig names agree" >&2
contigs "$TRUTH" > "$TMP/c_truth"
contigs "$POST"  > "$TMP/c_pred"
cut -f1 "$FAI" | sort > "$TMP/c_fai"
N_TRUTH=$(wc -l < "$TMP/c_truth"); N_PRED=$(wc -l < "$TMP/c_pred")
N_SHARED=$(comm -12 "$TMP/c_truth" "$TMP/c_pred" | wc -l)
N_FAI_SHARED=$(comm -12 "$TMP/c_truth" "$TMP/c_fai" | wc -l)
echo "[concord]   truth=$N_TRUTH contigs, prediction=$N_PRED, shared=$N_SHARED, truth∩fai=$N_FAI_SHARED" >&2
if [ "$N_SHARED" -eq 0 ]; then
  echo "[concord] FATAL: truth and prediction share NO contig names." >&2
  echo "[concord]   truth examples:      $(head -3 "$TMP/c_truth" | tr '\n' ' ')" >&2
  echo "[concord]   prediction examples: $(head -3 "$TMP/c_pred"  | tr '\n' ' ')" >&2
  echo "[concord]   Rename one side to match (e.g. _PATERNAL -> _hap1) and re-run." >&2
  exit 2
fi
if [ "$N_FAI_SHARED" -eq 0 ]; then
  echo "[concord] FATAL: truth shares no contig names with --fai ($FAI)." >&2
  echo "[concord]   fill_gaps would label the whole genome 'nonrepeat'." >&2
  exit 2
fi
if [ "$N_SHARED" -lt "$N_TRUTH" ] || [ "$N_SHARED" -lt "$N_PRED" ]; then
  echo "[concord] WARNING: contig sets only partially overlap ($N_SHARED shared);" >&2
  echo "[concord]   non-shared contigs contribute nothing and are silently dropped." >&2
  comm -23 "$TMP/c_truth" "$TMP/c_pred" | head -5 | sed 's/^/[concord]   truth-only: /' >&2
  comm -13 "$TMP/c_truth" "$TMP/c_pred" | head -5 | sed 's/^/[concord]   pred-only:  /' >&2
fi

# 1. truth -> k-mer space, N-runs excluded
# INVARIANT: KaryoScope annotations, fill_gaps/trim_bed_ends output, and the
# pipeline's truth BEDs are ALL already in genome (fai) order. Do NOT sort:
# `sort -k1,1` is lexicographic (chr1, chr10, chr11...) and would destroy genome
# order, leaving unionbedg with inputs in mismatched chromosome order. The
# published smk (compare_HG002/compare_CHM13) sorts nothing for this reason.
echo "[concord] using inputs in genome order (no re-sorting)" >&2
decompress "$TRUTH" > "$TMP/truth.sorted.bed"
if [ "$NOFILL" -eq 1 ]; then
  # --no-fill: the "truth" is itself a KaryoScope annotation that already tiles
  # k-mer space with no gaps. fill_gaps would inject the literal label
  # 'nonrepeat' into those non-existent gaps, which is meaningful only for a
  # RepeatMasker-style truth and would corrupt e.g. the region set. Differences
  # in k-mer-space extent between two k values are handled by unionbedg's "0"
  # filler, which contrast_HG002.py already skips.
  echo "[concord] --no-fill: using truth as-is (annotation-vs-annotation)" >&2
  cp "$TMP/truth.sorted.bed" "$TMP/truth.filled.bed"
else
  FILL_ARGS=(--fai "$FAI" --input "$TMP/truth.sorted.bed" --output "$TMP/truth.filled.bed" --kmer-size "$KMER")
  [ -n "$GAPS" ] && FILL_ARGS+=(--gaps "$GAPS")
  echo "[concord] fill_gaps (trim k-1 = $((KMER-1)) bp per contig; mark N-runs exclude)" >&2
  python3 "$ARCH/fill_gaps.py" "${FILL_ARGS[@]}" >&2
fi

# 2. one coordinate grid
decompress "$PRE"  > "$TMP/pre.bed"
decompress "$POST" > "$TMP/post.bed"
echo "[concord] unionbedg" >&2
$BT unionbedg -i "$TMP/truth.filled.bed" "$TMP/pre.bed" "$TMP/post.bed" > "$TMP/comparison.bed"

# 3. hierarchy-aware contrast, per base, both stages
echo "[concord] contrast_HG002.py" >&2
CARGS=(--input "$TMP/comparison.bed" --hierarchy "$HIER" --feature-set "$FS" --output "${OUT}.contrast.tsv")
[ -n "$MERGE" ] && CARGS+=(--merge "$MERGE")
python3 "$ARCH/contrast_HG002.py" "${CARGS[@]}" >&2

# 4. collapse to per-stage class fractions
python3 - "${OUT}.contrast.tsv" "${OUT}.summary.tsv" <<'PY'
import sys, collections
inp, out = sys.argv[1], sys.argv[2]
tot=collections.defaultdict(int); cls=collections.defaultdict(int)
with open(inp) as fh:
    next(fh)
    for line in fh:
        stage, truth, pred, length, compat = line.rstrip("\n").split("\t")
        length=int(length)
        # collapse the wrong_*/X_as_Y families into one 'incompatible' class,
        # keeping equivalent/compatible/novel exactly as the paper reports them
        c = compat if compat in ("equivalent","compatible","novel") else "incompatible"
        tot[stage]+=length; cls[(stage,c)]+=length
with open(out,"w") as o:
    o.write("stage\tclass\tbases\tfraction_pct\n")
    for stage in ("presmoothed","smoothed"):
        for c in ("equivalent","compatible","novel","incompatible"):
            b=cls[(stage,c)]; t=tot[stage]
            o.write(f"{stage}\t{c}\t{b}\t{100.0*b/t if t else 0:.4f}\n")
print(open(out).read())
PY
# 5. COVERAGE ASSERTION.
# A chromosome-order mismatch between unionbedg inputs makes bedtools emit "0"
# filler for whole chromosomes, which contrast_HG002.py skips. That silently
# dropped 40.6% of the genome while leaving all the PERCENTAGES looking right --
# it even reproduced the published anchors. Percentages alone cannot detect it;
# total bases scored can.
EXPECTED=$(awk '{s+=$2}END{print s}' "$FAI")
SCORED=$(awk -F'\t' 'NR>1 && $1=="presmoothed"{s+=$4}END{print s+0}' "${OUT}.contrast.tsv")
PCT=$(awk -v a="$SCORED" -v b="$EXPECTED" 'BEGIN{printf "%.1f", (b?100*a/b:0)}')
echo "[concord] coverage: ${SCORED} / ${EXPECTED} bp (${PCT}% of the .fai genome)" >&2
LOW=$(awk -v p="$PCT" 'BEGIN{print (p<90)?1:0}')
if [ "$LOW" -eq 1 ]; then
  echo "[concord] FATAL: only ${PCT}% of the genome was scored (expected >90%)." >&2
  echo "[concord]   Usual cause: truth and predictions in different chromosome ORDER," >&2
  echo "[concord]   so unionbedg emitted \"0\" filler that contrast then skipped." >&2
  echo "[concord]   Everything must be in genome (fai) order - do NOT sort -k1,1." >&2
  rm -f "${OUT}.summary.tsv"
  exit 4
fi
echo "[concord] wrote ${OUT}.contrast.tsv and ${OUT}.summary.tsv" >&2
