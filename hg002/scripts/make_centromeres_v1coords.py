#!/usr/bin/env python3
"""Write the Fig. 1b centromere window file: v2 scaffold names with v1 centromere coordinates.

conf/hg002.yaml reads
  data/KS_human_CHM13_v2/HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres_v1coords.bed
for the centromere zoom windows. It is built from the two scaffolded centromere BEDs of the HG002 v1.1
assembly, which KaryoScope writes one row per scaffolded chromosome:
  v1: HG002.assembly.KS_human_CHM13.region.presmoothed.scaffolded.centromeres.bed      (KS_human_CHM13)
  v2: HG002.assembly.KS_human_CHM13_v2.region.presmoothed.scaffolded.centromeres.bed.gz (KS_human_CHM13_v2)
The two files list the same 46 chromosomes in the same order. Names differ only in the last field
(v1 `chr1_hap1_chr1_hap1`, v2 `chr1_hap1_chr1_PATERNAL`). Row i of the output is the name of v2 row i with
the coordinates of v1 row i; the script stops if the rows do not pair up by chromosome and haplotype.

Usage: make_centromeres_v1coords.py <v1.bed> <v2.bed[.gz]> > <out.bed>
"""
import gzip
import sys


def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    v1, v2 = rows(sys.argv[1]), rows(sys.argv[2])
    if len(v1) != len(v2):
        sys.exit(f"row counts differ: v1 {len(v1)}, v2 {len(v2)}")
    for a, b in zip(v1, v2):
        # chrN_hapH prefix (chromosome and haplotype of the scaffold) must agree
        if a[0].split("_")[:2] != b[0].split("_")[:2]:
            sys.exit(f"rows do not pair up: {a[0]} vs {b[0]}")
        sys.stdout.write("\t".join([b[0]] + a[1:]) + "\n")


if __name__ == "__main__":
    main()
