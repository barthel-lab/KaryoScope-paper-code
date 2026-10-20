#!/usr/bin/env python3
"""
trim_bed_ends.py

Trim the last k-1 bases from each chromosome in a complete-tiling BED
file.  This is used when benchmarking KaryoScope against truth annotations
(e.g., from get_repeat_featureset.py or get_region_featureset.py) to
exclude positions that KaryoScope cannot annotate because there are not
enough bases to form a complete k-mer at contig ends.

The input BED is assumed to be a complete tiling (every base covered
exactly once, as produced by the gap-filling scripts).

Usage:
    python trim_bed_ends.py --input truth.bed --fai genome.fai --kmer-size 31 --output trimmed.bed
"""

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Trim the last k-1 bases from each chromosome in a tiled BED."
    )
    parser.add_argument("--input", required=True,
                        help="Input BED file (complete tiling).")
    parser.add_argument("--fai", required=True,
                        help="Chromosome sizes file (.fai).")
    parser.add_argument("--kmer-size", type=int, required=True,
                        help="K-mer size used by KaryoScope.")
    parser.add_argument("--output", required=True,
                        help="Output trimmed BED file.")
    args = parser.parse_args()

    trim = args.kmer_size - 1
    if trim <= 0:
        print("kmer-size must be > 1; nothing to trim.", file=sys.stderr)
        sys.exit(1)

    # Read chromosome lengths
    chrom_lengths = {}
    with open(args.fai) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])

    # Compute effective end per chromosome
    effective_end = {
        chrom: length - trim for chrom, length in chrom_lengths.items()
    }

    print(f"Trimming last {trim} bp (kmer-size={args.kmer_size}) "
          f"from {len(chrom_lengths)} chromosomes.", file=sys.stderr)

    with open(args.input) as fin, open(args.output, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, start, stop = parts[0], int(parts[1]), int(parts[2])
            label = parts[3]

            eff = effective_end.get(chrom, stop)

            # Skip features entirely beyond the effective end
            if start >= eff:
                continue

            # Clamp features that extend past the effective end
            stop = min(stop, eff)

            fout.write(f"{chrom}\t{start}\t{stop}\t{label}\n")

    print("Trimming complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
