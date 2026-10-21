#!/usr/bin/env python3
"""Step 2 (+ step 3) of the centromere QC chain: count QC flags per region, then drop the one bad contig.

A Python port of aggregate_qc_v2.awk (sha256 266dc648...4914), which needs gawk for arrays of arrays.

Step 2 (the awk logic, unchanged):
  * skip the header and every N/A row (not a chromosome-level contig);
  * key = sample, haplotype, chromosome, region, contig (the tool column is ignored, so the three
    tools are merged);
  * count the nine flags below, ignoring "." (a region/tool with no call);
  * print a key when the key changes (a control break), so the input order matters and a key that
    comes back after a different key is printed again.
  The header keeps the awk's column names: `printf "%s_n_%s", OFS, flag` gives "_n_<flag>".
  On master_qc_table.full.tsv this gives a header + 18,416 rows (aggregate_qc_v4.tsv).

Step 3 (--curate; Yi-An's manual edit, made explicit):
  drop the rows of HG03486 hap1 chr8 on contig HG03486#1#CM091971.1. That contig is 90,381 bp, a
  chromosome-level (CM) accession that KaryoScope assigned to chr8, but it is not chromosome-length;
  the real chr8 is HG03486#1#CM091976.1. It is the only chromosome with two CM contigs, and it appears
  in the master table only as the contig name of the "." rows (no call). This gives a header +
  18,414 rows (tables/aggregate_qc_v4.curated.tsv). Yi-An Chen's hand-edited file
  (aggregate_qc_v4.manual_curation.tsv) holds the same rows; it differs only in line endings (CRLF,
  and no terminator after its last row).

Usage: aggregate_qc.py master_qc_table.full.tsv [--curate] > out.tsv
"""
import argparse
import collections
import sys

FLAGS = ["Col", "Dup", "Err", "NNN", "COLLAPSE", "COLLAPSE_OTHER", "COLLAPSE_VAR", "MISJOIN", "HET"]
KEY_COLS = ["sample", "haplotype", "chromosome", "region", "contig"]

# (sample, haplotype, chromosome, contig) rows removed by --curate
EXCLUDE = {
    # 90,381 bp CM contig labelled chr8; not chromosome-length (real chr8 = HG03486#1#CM091976.1)
    ("HG03486", "hap1", "chr8", "HG03486#1#CM091971.1"),
}


def aggregate(lines):
    """Yield one output row (list of str) per control-break group."""
    prev, counts = None, collections.Counter()
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if f[8] == "N/A":
            continue
        key = f[0:4] + [f[5]]
        if prev is not None and key != prev:
            yield prev + [str(counts[x]) for x in FLAGS]
            counts = collections.Counter()
        prev = key
        if f[8] != ".":
            counts[f[8]] += 1
    if prev is not None:
        yield prev + [str(counts[x]) for x in FLAGS]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table")
    ap.add_argument("--curate", action="store_true", help="apply step 3 (drop the EXCLUDE rows)")
    a = ap.parse_args()
    out = sys.stdout
    out.write("\t".join(KEY_COLS + ["_n_" + x for x in FLAGS]) + "\n")
    with open(a.master_table) as fh:
        next(fh)
        for row in aggregate(fh):
            if a.curate and (row[0], row[1], row[2], row[4]) in EXCLUDE:
                continue
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
