#!/usr/bin/env python3
"""Annotate and scaffold one ROB assembly with KaryoScope v1.1.0 (Fig. 3, Extended Data Fig. 7 inputs).

Runs KaryoScope v1.1.0's scaffold pipeline in BED mode on the two haplotype FASTAs of one assembly
against the `KS_human_CHM13` (v1, KMC) database. Missing annotations are made first (`annotate`
with smoothing, keeping the presmoothed BEDs); the scaffolded BEDs are written from the
PRESMOOTHED annotations, which is what both figures read. The `karyoscope scaffold` command line
of v1.1.0 writes only smoothed scaffolded BEDs, hence this call into the same function.

Outputs in OUT_DIR, per haplotype H (1, 2):
  <stem H>.KS_human_CHM13.<set>.presmoothed.bed.gz             per-base annotation
  <stem H>.KS_human_CHM13.<set>.presmoothed.scaffolded.bed.gz  the same, on scaffold names
  <stem H>.KS_human_CHM13.scaffold_map.tsv                     contig -> scaffold name, orientation
where <stem H> is the FASTA file name without .fasta.gz / .fa.gz.

Needs KaryoScope v1.1.0 importable (tag 4c09ce5) and its get_featureIDs binary
(KARYOSCOPE_GET_FEATUREIDS), seqtk on PATH, and about 40 GB of memory (KMC database).

usage: annotate_scaffold_v110.py --hap1 H1.fasta.gz --hap2 H2.fasta.gz --db-root ROOT -o OUT_DIR [-t 16]
"""
import argparse
import logging
from pathlib import Path

from karyoscope.core.scaffold_run import InputSpec, scaffold_run

FEATURE_SETS = ["chromosome", "region", "repeat", "gene", "acrocentric", "subtelomeric"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hap1", type=Path, required=True)
    ap.add_argument("--hap2", type=Path, required=True)
    ap.add_argument("--db-root", type=Path, required=True)
    ap.add_argument("-o", "--out-dir", type=Path, required=True)
    ap.add_argument("-t", "--threads", type=int, default=16)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    scaffold_run(
        [InputSpec("hap1", a.hap1), InputSpec("hap2", a.hap2)],
        db_root=a.db_root, db_id="KS_human_CHM13", feature_sets=FEATURE_SETS,
        mode="bed", threads=a.threads, bgzip=True, output_dir=a.out_dir,
        annotation_variant="presmoothed")
