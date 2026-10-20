#!/usr/bin/env python3
"""
Supplementary Figure 1: Representative D4Z4 haplotypes, validation, and pipeline.

Panel (a): Representative 4qA vs 4qB diploid (HG01358 chr4)
Panel (b): KaryoScope vs HMM repeat count correlation
Panel (c): Classification pipeline diagram (existing)
Panel (d): Classification decision tree

Outputs individual panel SVGs to results/karyoscope_companion/
"""

import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'results' / 'karyoscope_companion'

# Register Barthel fonts
FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
FONT_NAME = 'sans-serif'
plt.rcParams['svg.fonttype'] = 'none'
if FONT_DIR.exists():
    for font_file in FONT_DIR.glob("BasicSans-*.otf"):
        fm.fontManager.addfont(str(font_file))
    FONT_NAME = 'Arial'
    plt.rcParams['font.family'] = FONT_NAME

# Barthel colors
COLORS = {
    'green': '#40D392',
    'blue': '#60A5FA',
    'coral': '#F07167',
    'gray': '#545454',
    'black': '#000000',
    'lavender': '#C4A9E8',
}


def generate_panel_a():
    """Representative 4qA vs 4qB diploid (HG01358 chr4)."""
    output_path = OUTPUT_DIR / 'supp1_panel_a_representative.svg'
    cmd = [
        sys.executable, str(BASE_DIR / 'plot_d4z4_arrays.py'),
        '--samples', 'HG01358',
        '--chromosome', 'chr4',
        '--output', str(output_path),
        '--background', 'white',
        '--width', '1200',
        '--adaptive-window',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
    else:
        print(f"  Saved: {output_path}")


def generate_panel_b():
    """KaryoScope vs HMM repeat count correlation."""
    output_path = OUTPUT_DIR / 'supp1_panel_b_hmm_correlation'

    # Load analysis results
    analysis_df = pd.read_csv(BASE_DIR / 'results' / 'd4z4_array_analysis.tsv', sep='\t')

    # Load HMM counts from Hailey's data
    hmm_file = BASE_DIR / 'Hailey-notes' / 'HPRC_chr4_counts.xlsx'
    if not hmm_file.exists():
        print(f"  WARNING: HMM counts file not found: {hmm_file}")
        return

    hmm_df = pd.read_excel(hmm_file)

    # Match samples: single-array chr4 only
    chr4_single = analysis_df[
        (analysis_df['chromosome'] == 'chr4') &
        (analysis_df['num_arrays'] == 1) &
        (~analysis_df['is_truncated'])
    ].copy()

    # Build HMM lookup by contig name
    hmm_lookup = {}
    for _, hmm_row in hmm_df.iterrows():
        hmm_lookup[str(hmm_row['Contig'])] = int(hmm_row['Count'])

    # Match by contig name (strip chr4_hap1_ prefix from haplotype names)
    merged = []
    for _, row in chr4_single.iterrows():
        hap = row['haplotype']
        contig = hap.split('_', 2)[-1] if '_' in hap else hap
        if contig in hmm_lookup:
            merged.append({
                'haplotype': hap,
                'karyoscope_count': row['repeat_count_main'],
                'hmm_count': hmm_lookup[contig],
            })

    if not merged:
        print("  WARNING: No matching samples for HMM correlation")
        return

    merged_df = pd.DataFrame(merged)

    # Calculate correlation
    r, p = stats.pearsonr(merged_df['karyoscope_count'], merged_df['hmm_count'])

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    ax.scatter(merged_df['hmm_count'], merged_df['karyoscope_count'],
               alpha=0.6, s=30, color=COLORS['blue'], edgecolor='white', linewidth=0.5)

    # Identity line
    max_val = max(merged_df['hmm_count'].max(), merged_df['karyoscope_count'].max()) + 5
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=1)

    # Regression line
    slope, intercept = np.polyfit(merged_df['hmm_count'], merged_df['karyoscope_count'], 1)
    x_line = np.linspace(0, max_val, 100)
    ax.plot(x_line, slope * x_line + intercept, color=COLORS['coral'], linewidth=1.5, alpha=0.8)

    ax.set_xlabel('HMM repeat count', fontsize=12)
    ax.set_ylabel('KaryoScope repeat count', fontsize=12)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')

    ax.text(0.05, 0.95, f'r = {r:.3f}\nn = {len(merged_df)}',
            transform=ax.transAxes, fontsize=11, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(str(output_path) + '.svg', dpi=300, bbox_inches='tight')
    fig.savefig(str(output_path) + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}.svg/.png (r={r:.3f}, n={len(merged_df)})")


def generate_panel_c():
    """Classification pipeline diagram — use existing script."""
    output_path = OUTPUT_DIR / 'supp1_panel_c_pipeline.svg'
    cmd = [
        sys.executable, str(BASE_DIR / 'create_pipeline_diagram.py'),
        '--output', str(output_path),
        '--background', 'white',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
    else:
        print(f"  Saved: {output_path}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Supplementary Figure 1: Representative haplotypes, validation, pipeline")
    print("=" * 60)

    print("\nPanel (a): Representative 4qA vs 4qB (HG01358)")
    generate_panel_a()

    print("\nPanel (b): KaryoScope vs HMM correlation")
    generate_panel_b()

    print("\nPanel (c): Pipeline diagram")
    generate_panel_c()

    print(f"\nAll panels saved to {OUTPUT_DIR}/supp1_panel_*.svg")
    print("Panel (d) decision tree: create manually or in Illustrator.")


if __name__ == '__main__':
    main()
