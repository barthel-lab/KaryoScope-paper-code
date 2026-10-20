#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# check_data.sh
#
# Verify that data/ contains the expected per-sample KaryoScope BEDs for v1
# and v2. Checks per-version:
#   - centromere BED(s)
#   - scaffolded feature BEDs across featuresets × bin sizes
#
# Does NOT check database.info.txt files -- those are bundled with the
# KaryoScope software repo at
#   $KS_RESOURCES/databases/KS_human_CHM13{,_v2}/
# and are not expected inside this project's data/ directory.
#
# Exits non-zero if anything required is missing. Intended to be run after
# the data were copied into data/.
# -----------------------------------------------------------------------------
set -u

LOCAL_BASE="$(cd "$(dirname "$0")/.." && pwd)/data"
# v1 ships 5 scaffolded featuresets (no acrocentric beds); v2 adds acrocentric.
FEATURESETS_V1=(chromosome region repeat gene subtelomeric)
FEATURESETS_V2=(chromosome region repeat gene acrocentric subtelomeric)
BINSIZES=(binned100 binned100000 binned1000000)

missing=0
report_missing() { echo "  MISSING: $1"; missing=$((missing + 1)); }

check_one() {
    local pattern="$1"
    # shellcheck disable=SC2086
    if ! compgen -G "$pattern" > /dev/null; then
        report_missing "$pattern"
    fi
}

check_version() {
    local version="$1"
    local prefix="$2"         # e.g. HG002.assembly.KS_human_CHM13 or ..._v2
    # Remaining args are the featureset list (avoids bash 3.2 nameref issues on macOS)
    shift 2
    local featuresets=("$@")
    local dir="$LOCAL_BASE/$version"

    echo "==> $version ($dir)"
    if [[ ! -d "$dir" ]]; then
        echo "  MISSING: directory $dir"
        missing=$((missing + 1))
        return
    fi

    # centromere coords
    check_one "$dir/${prefix}.region.presmoothed.scaffolded.centromeres.bed"
    if [[ "$version" == "KS_human_CHM13_v2" ]]; then
        check_one "$dir/${prefix}.region.presmoothed.scaffolded.centromeres_v1coords.bed"
    fi

    # scaffolded feature BEDs across featuresets × bin sizes
    for fs in "${featuresets[@]}"; do
        for bin in "${BINSIZES[@]}"; do
            check_one "$dir/${prefix}.${fs}.presmoothed.scaffolded.${bin}.features.bed*"
        done
    done
}

check_version "KS_human_CHM13_v1" "HG002.assembly.KS_human_CHM13"    "${FEATURESETS_V1[@]}"
check_version "KS_human_CHM13_v2" "HG002.assembly.KS_human_CHM13_v2" "${FEATURESETS_V2[@]}"

echo
if (( missing == 0 )); then
    echo "OK: all expected files present."
    exit 0
else
    echo "FAIL: $missing expected files/patterns missing. Copy the missing files into data/."
    exit 1
fi
