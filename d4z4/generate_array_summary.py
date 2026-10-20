#!/usr/bin/env python3
"""
Generate D4Z4 Array Summary Table

Creates a summary table showing for each haplotype:
- Number of arrays
- Number of canonical vs degraded arrays
- Array types in order (C+=canonical with β-satellite, C-=canonical without, D=degraded)
- Number of elements in each array

Usage:
    python generate_array_summary.py
    python generate_array_summary.py --input results/d4z4_array_analysis.tsv --output results/d4z4_array_summary.tsv
"""

import argparse
import json
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Generate D4Z4 array summary table")
    parser.add_argument("--input", default="results/d4z4_array_analysis.tsv",
                        help="Input analysis TSV file")
    parser.add_argument("--output", default="results/d4z4_array_summary.tsv",
                        help="Output summary TSV file")
    return parser.parse_args()


def generate_summary(input_file, output_file):
    """Generate summary table from analysis results."""

    # Load the analysis results
    df = pd.read_csv(input_file, sep='\t')

    rows = []

    for _, row in df.iterrows():
        chromosome = row['chromosome']
        sample = row['sample']
        hap_num = row['hap_num']
        haplotype = row['haplotype']

        # Extract contig from haplotype (e.g., "chr10_hap2_HG00146#2#CM090038.1" -> "CM090038.1")
        contig = haplotype.split('#')[-1] if '#' in haplotype else haplotype

        # Parse arrays_info JSON
        arrays_info_str = row.get('arrays_info')
        if pd.isna(arrays_info_str) or not arrays_info_str:
            # Single array fallback
            repeat_count = int(row.get('repeat_count_main', 0)) if pd.notna(row.get('repeat_count_main')) else 0
            has_plam = row.get('plam_detected', False)
            if pd.isna(has_plam):
                has_plam = False

            # Canonical = has terminal, Degraded = everything else
            # C+ = has β-satellite terminal, C- = no β-satellite (but has ct_specific)
            arr_type = 'canonical' if (repeat_count >= 3 and has_plam) else 'degraded'
            arr_code = 'C+' if arr_type == 'canonical' else 'D'

            rows.append({
                'chromosome': chromosome,
                'sample': sample,
                'hap': hap_num,
                'contig': contig,
                'num_arrays': 1,
                'num_canonical': 1 if arr_type == 'canonical' else 0,
                'num_degraded': 1 if arr_type == 'degraded' else 0,
                'array_types': arr_code,
                'array_sizes': str(repeat_count)
            })
            continue

        try:
            arrays_info = json.loads(arrays_info_str)
        except (json.JSONDecodeError, TypeError):
            continue

        num_canonical = 0
        num_degraded = 0
        array_types = []
        array_sizes = []

        for arr in arrays_info:
            arr_type_str = arr.get('type', 'unknown')
            has_plam = arr.get('has_plam', False)
            element_count = arr.get('element_count', 0) or 0

            # C+ = canonical with β-satellite, C- = canonical without, D = degraded
            if arr_type_str == 'canonical':
                arr_code = 'C+' if has_plam else 'C-'
                num_canonical += 1
            else:
                arr_code = 'D'
                num_degraded += 1

            array_types.append(arr_code)
            array_sizes.append(str(int(element_count)))

        rows.append({
            'chromosome': chromosome,
            'sample': sample,
            'hap': hap_num,
            'contig': contig,
            'num_arrays': len(arrays_info),
            'num_canonical': num_canonical,
            'num_degraded': num_degraded,
            'array_types': ','.join(array_types),
            'array_sizes': ','.join(array_sizes)
        })

    # Create DataFrame and save
    summary_df = pd.DataFrame(rows)
    summary_df = summary_df.sort_values(['chromosome', 'sample', 'hap'])

    # Save to TSV
    summary_df.to_csv(output_file, sep='\t', index=False)

    return summary_df


def main():
    args = parse_args()

    print(f"Generating D4Z4 array summary...")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")

    summary_df = generate_summary(args.input, args.output)

    print(f"\n=== Legend ===")
    print(f"C+ = Canonical array with β-satellite terminal (has pLAM)")
    print(f"C- = Canonical array without β-satellite terminal (has ct_specific only)")
    print(f"D  = Degraded array (no terminal or <3 elements)")

    print(f"\n=== Summary ===")
    print(f"Total haplotypes: {len(summary_df)}")
    print(f"Canonical arrays: {summary_df['num_canonical'].sum()}")
    print(f"Degraded arrays: {summary_df['num_degraded'].sum()}")

    print(f"\n=== Array Count Distribution ===")
    print(summary_df['num_arrays'].value_counts().sort_index().to_string())

    print(f"\nSaved to: {args.output}")


if __name__ == '__main__':
    main()
