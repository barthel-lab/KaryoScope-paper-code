#!/usr/bin/env bash
# Download the HPRC Release 2 pre-release QC calls listed in prerelease_qc_index.tsv and check their sha256.
#
# Usage: bash download_prerelease_qc.sh <out_dir> [prerelease_qc_index.tsv]
#
# Writes <out_dir>/<tool>/<file> (1,386 BED files, 462 per tool), the layout that
# build_master_qc_table.py --qc-dir reads. The bucket is public: the AWS CLI is used without
# credentials when it is installed, otherwise curl over HTTPS. Files already present with the
# right sha256 are skipped. Exits non-zero if any checksum differs.
set -euo pipefail
out=${1:?usage: download_prerelease_qc.sh <out_dir> [index]}
index=${2:-$(dirname "$0")/prerelease_qc_index.tsv}
tail -n +2 "$index" | while IFS=$'\t' read -r sample hap tool file sha idx idx_sample idx_hap idx_asm url; do
    mkdir -p "$out/$tool"
    dest="$out/$tool/$file"
    if [[ -f "$dest" ]] && echo "$sha  $dest" | sha256sum -c --quiet - 2>/dev/null; then
        continue
    fi
    if command -v aws >/dev/null 2>&1; then
        aws s3 --no-sign-request cp --only-show-errors "$url" "$dest"
    else
        curl -fsSL -o "$dest" "https://s3-us-west-2.amazonaws.com/${url#s3://}"
    fi
    echo "$sha  $dest" | sha256sum -c --quiet - || { echo "sha256 mismatch: $dest" >&2; exit 1; }
done
echo "all files in $index present in $out with the listed sha256"
