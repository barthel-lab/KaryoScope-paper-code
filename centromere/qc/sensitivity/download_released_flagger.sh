#!/usr/bin/env bash
# Download the released HMM-Flagger v1.2.0 conservative calls listed in released_flagger_index.tsv
# (one BED per sample and platform, both haplotypes) to <out_dir>/<file> and check their sha256.
# Usage: bash download_released_flagger.sh <out_dir> [released_flagger_index.tsv]
set -euo pipefail
out=${1:?usage: download_released_flagger.sh <out_dir> [index]}
index=${2:-$(dirname "$0")/released_flagger_index.tsv}
mkdir -p "$out"
tail -n +2 "$index" | while IFS=$'\t' read -r sample tool file sha url; do
    dest="$out/$file"
    if [[ -f "$dest" ]] && echo "$sha  $dest" | sha256sum -c --quiet - 2>/dev/null; then
        continue
    fi
    curl -fsSL -o "$dest" "$url"
    echo "$sha  $dest" | sha256sum -c --quiet - || { echo "sha256 mismatch: $dest" >&2; exit 1; }
done
echo "all files in $index present in $out with the listed sha256"
