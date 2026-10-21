#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# check_data.sh
#
# Verify that data/robs/ contains the expected per-sample KaryoScope BEDs for
# all three Robertsonian samples (GM03417, GM03786, GM04890). Checks per sample:
#   - 6 featuresets x 3 bin sizes scaffolded feature BEDs
#   - centromere BED
#
# Does NOT check database.info.txt / color files -- those live in
# vendored/KaryoScope-archive-resources/ (resources/databases/KS_human_CHM13/)
# of KaryoScope-paper-code.
#
# Exits non-zero if anything required is missing. Intended to be run after
# the data were copied into data/robs/.
# -----------------------------------------------------------------------------
set -u

LOCAL_BASE="$(cd "$(dirname "$0")/.." && pwd)/data/robs"
SAMPLES=(GM03417 GM03786 GM04890)
FEATURESETS=(chromosome region repeat gene acrocentric subtelomeric)
BINSIZES=(binned100 binned100000 binned1000000)
DATABASE="KS_human_CHM13"

missing=0
report_missing() { echo "  MISSING: $1"; missing=$((missing + 1)); }

check_one() {
    local pattern="$1"
    # shellcheck disable=SC2086
    if ! compgen -G "$pattern" > /dev/null; then
        report_missing "$pattern"
    fi
}

echo "==> data/robs/ ($LOCAL_BASE)"
if [[ ! -d "$LOCAL_BASE" ]]; then
    echo "  MISSING: directory $LOCAL_BASE"
    exit 1
fi

for sample in "${SAMPLES[@]}"; do
    prefix="${sample}.assembly.${DATABASE}"
    echo "  -- ${sample}"

    # centromere coords
    check_one "${LOCAL_BASE}/${prefix}.region.presmoothed.scaffolded.centromeres.bed*"

    # scaffolded feature BEDs across featuresets x bin sizes
    for fs in "${FEATURESETS[@]}"; do
        for bin in "${BINSIZES[@]}"; do
            check_one "${LOCAL_BASE}/${prefix}.${fs}.presmoothed.scaffolded.${bin}.features.bed*"
        done
    done
done

echo
if (( missing == 0 )); then
    echo "OK: all expected files present."
    exit 0
else
    echo "FAIL: $missing expected files/patterns missing. Copy the missing files into data/robs/."
    exit 1
fi
