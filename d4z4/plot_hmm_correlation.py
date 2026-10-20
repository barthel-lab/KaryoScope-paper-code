#!/usr/bin/env python3
"""
Plot correlation between KaryoScope D4Z4 repeat counts and HMM-based counts.

Compares repeat counts from KaryoScope analysis with HMM-based counts from
Hailey's analysis, for single-array chr4 samples only.

Usage:
    python plot_hmm_correlation.py
    python plot_hmm_correlation.py --output-dir results
"""

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def parse_args():
    parser = argparse.ArgumentParser(description="Plot HMM vs KaryoScope correlation")
    parser.add_argument("--analysis-file",
                        default="results/d4z4_array_analysis.tsv",
                        help="D4Z4 analysis results TSV")
    parser.add_argument("--hmm-file",
                        default="Hailey-notes/HPRC_chr4_counts.xlsx",
                        help="HMM counts Excel file")
    parser.add_argument("--output-dir", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--output-prefix", default="hmm_correlation",
                        help="Output filename prefix (default: hmm_correlation)")
    return parser.parse_args()


def load_and_match_data(analysis_file, hmm_file):
    """Load KaryoScope and HMM data, match samples."""

    # Load KaryoScope analysis
    df = pd.read_csv(analysis_file, sep='\t')

    # Filter to single-array chr4 samples
    chr4_single = df[(df['chromosome'] == 'chr4') & (df['num_arrays'] == 1)].copy()
    print(f"Single-array chr4 samples in KaryoScope: {len(chr4_single)}")

    # Extract core haplotype ID (remove chr4_hap1/2_ prefix)
    # KaryoScope format: chr4_hap2_HG00099#2#CM087359.1 -> HG00099#2#CM087359.1
    chr4_single['core_haplotype'] = chr4_single['haplotype'].str.replace(
        r'^chr\d+_hap\d+_', '', regex=True
    )

    # Load HMM counts
    hmm_df = pd.read_excel(hmm_file)
    hmm_df = hmm_df.rename(columns={'Count': 'hmm_count', 'Contig': 'core_haplotype'})
    print(f"Samples in HMM file: {len(hmm_df)}")

    # Merge on core haplotype
    merged = chr4_single.merge(
        hmm_df[['core_haplotype', 'hmm_count']],
        on='core_haplotype',
        how='inner'
    )
    print(f"Matched samples: {len(merged)}")

    return merged


def plot_correlation(data, output_dir, output_prefix):
    """Create correlation scatter plot."""

    # Calculate correlation
    r, p = stats.pearsonr(data['repeat_count_main'], data['hmm_count'])
    print(f"\nCorrelation: r = {r:.3f}, p = {p:.2e}")

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # Scatter plot
    ax.scatter(
        data['hmm_count'],
        data['repeat_count_main'],
        alpha=0.6,
        s=50,
        c='steelblue',
        edgecolors='white',
        linewidths=0.5
    )

    # Perfect agreement line (identity line)
    max_val = max(data['hmm_count'].max(), data['repeat_count_main'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y = x')

    # Labels and title
    ax.set_xlabel('HMM Repeat Count', fontsize=12)
    ax.set_ylabel('KaryoScope Repeat Count', fontsize=12)
    ax.set_title(f'D4Z4 Repeat Count Correlation (chr4, single-array)\n'
                 f'n = {len(data)}, r = {r:.3f}', fontsize=14)

    # Equal aspect ratio
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, max_val * 1.05)
    ax.set_ylim(0, max_val * 1.05)

    # Grid
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower right')

    # Footer with timestamp
    fig.text(0.5, 0.02, f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
             ha='center', fontsize=8, color='gray')

    plt.tight_layout(rect=[0, 0.03, 1, 1])

    # Save outputs
    png_path = os.path.join(output_dir, f'{output_prefix}.png')
    pdf_path = os.path.join(output_dir, f'{output_prefix}.pdf')

    plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")

    return r, p


def print_outliers(data, threshold=10):
    """Print samples with large discrepancies."""

    data = data.copy()
    data['diff'] = abs(data['repeat_count_main'] - data['hmm_count'])
    outliers = data[data['diff'] >= threshold].sort_values('diff', ascending=False)

    if len(outliers) > 0:
        print(f"\nSamples with discrepancy >= {threshold}:")
        print("-" * 80)
        for _, row in outliers.head(10).iterrows():
            truncated = " (truncated)" if row.get('is_truncated', False) else ""
            print(f"  {row['haplotype']}: KS={row['repeat_count_main']}, "
                  f"HMM={row['hmm_count']}, diff={row['diff']}{truncated}")


def main():
    args = parse_args()

    print("=" * 60)
    print("HMM vs KaryoScope Correlation Analysis")
    print("=" * 60)

    # Load and match data
    data = load_and_match_data(args.analysis_file, args.hmm_file)

    # Plot correlation
    r, p = plot_correlation(data, args.output_dir, args.output_prefix)

    # Show outliers
    print_outliers(data)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Matched samples: {len(data)}")
    print(f"  Pearson r: {r:.3f}")
    print(f"  p-value: {p:.2e}")


if __name__ == '__main__':
    main()
