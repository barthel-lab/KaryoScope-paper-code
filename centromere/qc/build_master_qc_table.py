#!/usr/bin/env python3
"""Build master_qc_table.full.tsv (step 1 of the centromere QC chain) without Snakemake.

This is a standalone re-implementation of three rules of KaryoScope-archive workflow/rules/pangenome.smk
(commit 0bf3a7d, run 2025-11-12/13 and 2025-11-19):
  get_regions               -> the region intervals are an INPUT here (--regions)
  find_errors               -> `bedtools intersect -a <QC calls> -b <region>`, re-implemented below
  aggregate_master_qc_table -> workflow/scripts/aggregate_master_qc.py, re-implemented below
It reproduces the original file byte for byte (sha256 2b7faf0a...f218, 225,079 lines = header + 225,078 rows).

Inputs
  --qc-index   TSV with columns sample, haplotype, tool, file (prerelease_qc_index.tsv).
               The table's SAMPLE ORDER is the order of first appearance here (the original took it from
               the Snakemake config: HG005 first, then the rest).
  --qc-dir     directory holding the QC files as <qc-dir>/<tool>/<file> (download_prerelease_qc.sh
               writes this layout).
  --contigs    TSV sample, haplotype, chromosome, contig: the chromosome-level (GenBank CM) contigs that
               KaryoScope assigned to each chromosome, in scaffold_stats order. A (sample, haplotype,
               chromosome) listed here is "chromosome-level"; one that is not gets N/A rows.
  --regions    TSV sample, haplotype, chromosome, region, contig, start, stop (0-based, half-open):
               the subtelo_p / centromere / subtelo_q intervals (one each; see README.md for how they
               were defined). A chromosome-level combination with no interval gets "." rows.

Conventions reproduced from the original
  * rows are written for every sample x hap(hap1, hap2) x chrom(chr1..22, X, Y) x region(subtelo_p,
    centromere, subtelo_q) x tool(flagger_hifi, flagger_ont, nucflag), in that nesting order;
  * not chromosome-level: one row `contig start stop flag` = `N/A N/A N/A N/A`;
  * chromosome-level, no call in the region (or no region): one row `<contig> . . .`, where <contig> is
    the LAST contig listed for that chromosome in --contigs (the original dict was overwritten per
    stats line). Only HG03486 hap1 chr8 lists two CM contigs; this is why its "." rows name
    CM091971.1 while its calls (on the region contig, the FIRST listed, CM091976.1) do;
  * otherwise one row per QC call that overlaps the region by >= 1 bp, in QC-file order, with start and
    stop CLIPPED to the region (bedtools intersect without -wa reports only the overlap).
"""
import argparse, collections, os, sys

HAPS = ["hap1", "hap2"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
REGIONS = ["subtelo_p", "centromere", "subtelo_q"]
TOOLS = ["flagger_hifi", "flagger_ont", "nucflag"]


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def load_calls(qc_path):
    """contig -> [(start, stop, flag), ...] in file order; header lines skipped, as bedtools does."""
    calls = collections.defaultdict(list)
    with open(qc_path) as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            calls[f[0]].append((int(f[1]), int(f[2]), f[3]))
    return calls


def intersect(calls, contig, start, stop):
    """bedtools intersect -a <calls> -b <one interval>, default output, first four columns:
    every call with >= 1 bp overlap, in input order, clipped to the interval."""
    return [(contig, max(a, start), min(b, stop), flag)
            for a, b, flag in calls.get(contig, ()) if a < stop and b > start]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qc-index", required=True)
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()

    qc, samples = {}, []
    for r in read_tsv(a.qc_index):
        if r["sample"] not in samples:
            samples.append(r["sample"])
        qc[(r["sample"], r["haplotype"], r["tool"])] = os.path.join(a.qc_dir, r["tool"], r["file"])

    contigs = collections.defaultdict(list)
    for r in read_tsv(a.contigs):
        contigs[(r["sample"], r["haplotype"], r["chromosome"])].append(r["contig"])

    region = {}
    for r in read_tsv(a.regions):
        k = (r["sample"], r["haplotype"], r["chromosome"], r["region"])
        assert k not in region, f"two intervals for {k}"
        region[k] = (r["contig"], int(r["start"]), int(r["stop"]))

    n = 0
    with open(a.output, "w") as out:
        out.write("sample\thaplotype\tchromosome\tregion\ttool\tcontig\tstart\tstop\tflag\n")
        for s in samples:
            for h in HAPS:
                cache = {}
                for c in CHROMS:
                    listed = contigs.get((s, h, c))
                    for z in REGIONS:
                        for t in TOOLS:
                            key = f"{s}\t{h}\t{c}\t{z}\t{t}"
                            if not listed:
                                out.write(f"{key}\tN/A\tN/A\tN/A\tN/A\n"); n += 1
                                continue
                            hits = []
                            if (s, h, c, z) in region:
                                if t not in cache:
                                    cache[t] = load_calls(qc[(s, h, t)])
                                hits = intersect(cache[t], *region[(s, h, c, z)])
                            if not hits:
                                out.write(f"{key}\t{listed[-1]}\t.\t.\t.\n"); n += 1
                            for ctg, x, y, flag in hits:
                                out.write(f"{key}\t{ctg}\t{x}\t{y}\t{flag}\n"); n += 1
    print(f"wrote {n} rows (+ header) to {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
