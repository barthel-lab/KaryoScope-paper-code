#!/usr/bin/env python3
"""
D4Z4 Array Variability Analysis

Investigates whether D4Z4 array variability in HPRC samples is due to:
1. Normal Mendelian inheritance
2. Genetic remodeling (recombination, gene conversion, 4q-10q exchanges)

Analysis components:
- Trio inheritance analysis (using 1000 Genomes pedigree data)
- Population stratification
- Remodeling signature detection (4q-10q exchanges, fragmented arrays, asymmetry)

Usage:
    python analyze_d4z4_variability.py
    python analyze_d4z4_variability.py --background black --output-dir results/variability
"""

import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy import stats

# =============================================================================
# Configuration
# =============================================================================

PEDIGREE_URL = "https://raw.githubusercontent.com/human-pangenomics/HPRC_metadata/main/data/IGSR/20130606_g1k.ped"
D4Z4_SUMMARY_FILE = "results/d4z4_array_summary.tsv"

# Color palette from CLAUDE.md
COLORS = {
    'dark1': '#000000',
    'light1': '#FFFFFF',
    'dark2': '#545454',
    'light2': '#C4A9E8',
    'accent1': '#40D392',  # Green
    'accent2': '#60A5FA',  # Blue
    'accent3': '#F07167',  # Coral
    'accent4': '#FBBF24',  # Yellow
    'accent5': '#10B981',  # Emerald
    'accent6': '#3B82F6',  # Royal Blue
}

# Population colors
POP_COLORS = {
    'AFR': '#F07167',  # Coral
    'EUR': '#60A5FA',  # Blue
    'EAS': '#40D392',  # Green
    'AMR': '#FBBF24',  # Yellow
    'SAS': '#C4A9E8',  # Lavender
}

# Superpopulation mapping
SUPERPOP_MAP = {
    'ACB': 'AFR', 'ASW': 'AFR', 'ESN': 'AFR', 'GWD': 'AFR', 'LWK': 'AFR', 'MSL': 'AFR', 'YRI': 'AFR',
    'CDX': 'EAS', 'CHB': 'EAS', 'CHS': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS',
    'CEU': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 'TSI': 'EUR',
    'CLM': 'AMR', 'MXL': 'AMR', 'PEL': 'AMR', 'PUR': 'AMR',
    'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS',
}

# Configure matplotlib
plt.rcParams['svg.fonttype'] = 'none'

# Register custom fonts
FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
if FONT_DIR.exists():
    for font_file in FONT_DIR.glob("BasicSans-*.otf"):
        fm.fontManager.addfont(str(font_file))
    plt.rcParams['font.family'] = 'Basic Sans'

# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(output_dir):
    """Configure logging to file and console."""
    log_file = output_dir / f"variability_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_pedigree(url=PEDIGREE_URL):
    """
    Fetch and parse the 1000 Genomes pedigree file.

    Returns DataFrame with columns:
    - family_id, individual_id, father_id, mother_id, sex, population
    """
    logging.info(f"Fetching pedigree from {url}")

    response = requests.get(url)
    response.raise_for_status()

    lines = response.text.strip().split('\n')
    records = []

    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 7:
            records.append({
                'family_id': parts[0],
                'individual_id': parts[1],
                'father_id': parts[2] if parts[2] != '0' else None,
                'mother_id': parts[3] if parts[3] != '0' else None,
                'sex': 'M' if parts[4] == '1' else 'F' if parts[4] == '2' else 'U',
                'population': parts[6] if len(parts) > 6 else None,
            })

    df = pd.DataFrame(records)
    df['superpopulation'] = df['population'].map(SUPERPOP_MAP)

    logging.info(f"Loaded {len(df)} individuals from pedigree")
    logging.info(f"Populations: {df['population'].nunique()}, Superpopulations: {df['superpopulation'].nunique()}")

    return df


def load_d4z4_data(filepath=D4Z4_SUMMARY_FILE):
    """Load D4Z4 array summary data."""
    logging.info(f"Loading D4Z4 data from {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    logging.info(f"Loaded {len(df)} records for {df['sample'].nunique()} samples")
    return df


def identify_trios(pedigree_df, hprc_samples):
    """
    Identify complete trios where child AND both parents are in HPRC data.

    Returns list of dicts with child, father, mother sample IDs.
    """
    hprc_set = set(hprc_samples)
    trios = []

    for _, row in pedigree_df.iterrows():
        child = row['individual_id']
        father = row['father_id']
        mother = row['mother_id']

        if child in hprc_set and father in hprc_set and mother in hprc_set:
            trios.append({
                'child': child,
                'father': father,
                'mother': mother,
                'family_id': row['family_id'],
                'population': row['population'],
            })

    logging.info(f"Found {len(trios)} complete trios in HPRC data")
    return trios


# =============================================================================
# Analysis Functions
# =============================================================================

def analyze_trio_inheritance(d4z4_df, trios):
    """
    Compare child array sizes to parental ranges.

    For each trio:
    - Get array sizes for child, father, mother (both haplotypes)
    - Check if child arrays fall within parental range (±2 RU tolerance)
    - Flag potential de novo events
    """
    logging.info("Analyzing trio inheritance patterns...")
    results = []

    for trio in trios:
        child_data = d4z4_df[d4z4_df['sample'] == trio['child']]
        father_data = d4z4_df[d4z4_df['sample'] == trio['father']]
        mother_data = d4z4_df[d4z4_df['sample'] == trio['mother']]

        for chrom in ['chr4', 'chr10']:
            child_chr = child_data[child_data['chromosome'] == chrom]
            father_chr = father_data[father_data['chromosome'] == chrom]
            mother_chr = mother_data[mother_data['chromosome'] == chrom]

            if child_chr.empty or father_chr.empty or mother_chr.empty:
                continue

            # Parse array sizes
            def get_array_sizes(df):
                sizes = []
                for _, row in df.iterrows():
                    if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                        for s in str(row['array_sizes']).split(','):
                            try:
                                sizes.append(int(s.strip()))
                            except ValueError:
                                pass
                return sizes

            child_sizes = get_array_sizes(child_chr)
            father_sizes = get_array_sizes(father_chr)
            mother_sizes = get_array_sizes(mother_chr)

            if not child_sizes or not (father_sizes or mother_sizes):
                continue

            # Get parental range
            parent_sizes = father_sizes + mother_sizes
            parent_min = min(parent_sizes) - 2  # Allow ±2 RU tolerance
            parent_max = max(parent_sizes) + 2

            # Check each child array
            for i, child_size in enumerate(child_sizes):
                within_range = parent_min <= child_size <= parent_max

                # Find closest parent
                min_dist_father = min([abs(child_size - f) for f in father_sizes]) if father_sizes else float('inf')
                min_dist_mother = min([abs(child_size - m) for m in mother_sizes]) if mother_sizes else float('inf')

                results.append({
                    'family_id': trio['family_id'],
                    'child': trio['child'],
                    'father': trio['father'],
                    'mother': trio['mother'],
                    'population': trio['population'],
                    'chromosome': chrom,
                    'child_hap': i + 1,
                    'child_array_size': child_size,
                    'father_sizes': ','.join(map(str, father_sizes)),
                    'mother_sizes': ','.join(map(str, mother_sizes)),
                    'parent_min': parent_min + 2,  # Report actual min
                    'parent_max': parent_max - 2,  # Report actual max
                    'within_parental_range': within_range,
                    'min_distance_to_father': min_dist_father,
                    'min_distance_to_mother': min_dist_mother,
                    'likely_parent': 'father' if min_dist_father <= min_dist_mother else 'mother',
                    'potential_de_novo': not within_range and min(min_dist_father, min_dist_mother) > 5,
                })

    df = pd.DataFrame(results)

    if not df.empty:
        mendelian_rate = df['within_parental_range'].mean() * 100
        de_novo_count = df['potential_de_novo'].sum()
        logging.info(f"Mendelian inheritance rate: {mendelian_rate:.1f}%")
        logging.info(f"Potential de novo events: {de_novo_count}")

    return df


def analyze_population_stratification(d4z4_df, pedigree_df):
    """
    Compare array size distributions across populations.
    """
    logging.info("Analyzing population stratification...")

    # Merge population info
    pop_map = pedigree_df.set_index('individual_id')[['population', 'superpopulation']].to_dict('index')

    d4z4_df = d4z4_df.copy()
    d4z4_df['population'] = d4z4_df['sample'].map(lambda x: pop_map.get(x, {}).get('population'))
    d4z4_df['superpopulation'] = d4z4_df['sample'].map(lambda x: pop_map.get(x, {}).get('superpopulation'))

    results = []

    for chrom in ['chr4', 'chr10']:
        chrom_data = d4z4_df[d4z4_df['chromosome'] == chrom]

        for superpop in ['AFR', 'EUR', 'EAS', 'AMR', 'SAS']:
            pop_data = chrom_data[chrom_data['superpopulation'] == superpop]

            if pop_data.empty:
                continue

            # Get all array sizes
            sizes = []
            for _, row in pop_data.iterrows():
                if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                    for s in str(row['array_sizes']).split(','):
                        try:
                            sizes.append(int(s.strip()))
                        except ValueError:
                            pass

            if not sizes:
                continue

            results.append({
                'chromosome': chrom,
                'superpopulation': superpop,
                'n_samples': pop_data['sample'].nunique(),
                'n_arrays': len(sizes),
                'mean_size': np.mean(sizes),
                'median_size': np.median(sizes),
                'std_size': np.std(sizes),
                'min_size': min(sizes),
                'max_size': max(sizes),
                'range_size': max(sizes) - min(sizes),
            })

    df = pd.DataFrame(results)

    # Statistical test (Kruskal-Wallis) for each chromosome
    for chrom in ['chr4', 'chr10']:
        chrom_data = d4z4_df[d4z4_df['chromosome'] == chrom]
        groups = []

        for superpop in ['AFR', 'EUR', 'EAS', 'AMR', 'SAS']:
            pop_data = chrom_data[chrom_data['superpopulation'] == superpop]
            sizes = []
            for _, row in pop_data.iterrows():
                if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                    for s in str(row['array_sizes']).split(','):
                        try:
                            sizes.append(int(s.strip()))
                        except ValueError:
                            pass
            if sizes:
                groups.append(sizes)

        if len(groups) >= 2:
            stat, p_value = stats.kruskal(*groups)
            logging.info(f"{chrom} Kruskal-Wallis test: H={stat:.2f}, p={p_value:.4f}")

    return df, d4z4_df  # Return merged df for plotting


def detect_4q10q_exchanges(d4z4_df):
    """
    Detect potential 4q-10q exchanges based on C+/C- patterns.

    Expected patterns:
    - chr4: predominantly C- (no pLAM proxy)
    - chr10: predominantly C+ (has pLAM proxy)

    Flag samples with "wrong" patterns.
    """
    logging.info("Detecting potential 4q-10q exchanges...")
    results = []

    for sample in d4z4_df['sample'].unique():
        sample_data = d4z4_df[d4z4_df['sample'] == sample]

        for hap in ['hap1', 'hap2']:
            hap_data = sample_data[sample_data['hap'] == hap]

            chr4_data = hap_data[hap_data['chromosome'] == 'chr4']
            chr10_data = hap_data[hap_data['chromosome'] == 'chr10']

            # Get array types
            chr4_types = []
            chr10_types = []

            for _, row in chr4_data.iterrows():
                if pd.notna(row['array_types']) and row['array_types'] != '-':
                    chr4_types.extend(row['array_types'].split(','))

            for _, row in chr10_data.iterrows():
                if pd.notna(row['array_types']) and row['array_types'] != '-':
                    chr10_types.extend(row['array_types'].split(','))

            # Check for unexpected patterns
            chr4_has_cplus = 'C+' in chr4_types
            chr10_has_cminus = 'C-' in chr10_types

            if chr4_types or chr10_types:
                results.append({
                    'sample': sample,
                    'haplotype': hap,
                    'chr4_types': ','.join(chr4_types) if chr4_types else '-',
                    'chr10_types': ','.join(chr10_types) if chr10_types else '-',
                    'chr4_has_pLAM_proxy': chr4_has_cplus,
                    'chr10_missing_pLAM_proxy': chr10_has_cminus,
                    'potential_4q10q_exchange': chr4_has_cplus or chr10_has_cminus,
                })

    df = pd.DataFrame(results)

    if not df.empty:
        exchange_rate = df['potential_4q10q_exchange'].mean() * 100
        logging.info(f"Samples with potential 4q-10q exchange signature: {exchange_rate:.1f}%")

    return df


def detect_fragmented_arrays(d4z4_df):
    """
    Detect samples with fragmented arrays (>1 array per haplotype).

    Multiple arrays on a single contig suggest recombination breakpoints.
    """
    logging.info("Detecting fragmented arrays...")
    results = []

    for _, row in d4z4_df.iterrows():
        if row['num_arrays'] > 1:
            # Parse array sizes
            sizes = []
            if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                for s in str(row['array_sizes']).split(','):
                    try:
                        sizes.append(int(s.strip()))
                    except ValueError:
                        pass

            results.append({
                'sample': row['sample'],
                'haplotype': row['hap'],
                'chromosome': row['chromosome'],
                'contig': row['contig'],
                'num_arrays': row['num_arrays'],
                'array_sizes': row['array_sizes'],
                'total_repeat_units': sum(sizes) if sizes else 0,
                'largest_array': max(sizes) if sizes else 0,
                'fragmentation_ratio': len(sizes) / sum(sizes) if sizes and sum(sizes) > 0 else 0,
            })

    df = pd.DataFrame(results)

    if not df.empty:
        logging.info(f"Found {len(df)} haplotypes with fragmented arrays")
        logging.info(f"Max fragments: {df['num_arrays'].max()}")

    return df


def analyze_asymmetry(d4z4_df):
    """
    Analyze intra-individual haplotype asymmetry.

    Large differences between hap1 and hap2 could indicate de novo events.
    """
    logging.info("Analyzing intra-individual asymmetry...")
    results = []

    for sample in d4z4_df['sample'].unique():
        sample_data = d4z4_df[d4z4_df['sample'] == sample]

        for chrom in ['chr4', 'chr10']:
            chrom_data = sample_data[sample_data['chromosome'] == chrom]

            hap1_data = chrom_data[chrom_data['hap'] == 'hap1']
            hap2_data = chrom_data[chrom_data['hap'] == 'hap2']

            if hap1_data.empty or hap2_data.empty:
                continue

            # Get total repeat units per haplotype
            def get_total_ru(df):
                total = 0
                for _, row in df.iterrows():
                    if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                        for s in str(row['array_sizes']).split(','):
                            try:
                                total += int(s.strip())
                            except ValueError:
                                pass
                return total

            hap1_ru = get_total_ru(hap1_data)
            hap2_ru = get_total_ru(hap2_data)

            if hap1_ru == 0 and hap2_ru == 0:
                continue

            asymmetry = abs(hap1_ru - hap2_ru)
            fold_diff = max(hap1_ru, hap2_ru) / max(min(hap1_ru, hap2_ru), 1)

            results.append({
                'sample': sample,
                'chromosome': chrom,
                'hap1_total_ru': hap1_ru,
                'hap2_total_ru': hap2_ru,
                'asymmetry_ru': asymmetry,
                'fold_difference': fold_diff,
                'high_asymmetry': asymmetry > 20 or fold_diff > 3,
            })

    df = pd.DataFrame(results)

    if not df.empty:
        high_asym_rate = df['high_asymmetry'].mean() * 100
        logging.info(f"Samples with high haplotype asymmetry: {high_asym_rate:.1f}%")

    return df


# =============================================================================
# Visualization Functions
# =============================================================================

def setup_plot_style(background='white'):
    """Configure plot style based on background color."""
    if background == 'black':
        plt.style.use('dark_background')
        text_color = COLORS['light1']
        grid_color = COLORS['dark2']
    else:
        plt.style.use('default')
        text_color = COLORS['dark1']
        grid_color = '#CCCCCC'

    return text_color, grid_color


def plot_population_distributions(d4z4_merged, output_dir, background='white'):
    """Create box plots of array sizes by population."""
    text_color, grid_color = setup_plot_style(background)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[idx]
        chrom_data = d4z4_merged[d4z4_merged['chromosome'] == chrom]

        # Collect data by superpopulation
        plot_data = []
        labels = []
        colors = []

        for superpop in ['AFR', 'EUR', 'EAS', 'AMR', 'SAS']:
            pop_data = chrom_data[chrom_data['superpopulation'] == superpop]
            sizes = []
            for _, row in pop_data.iterrows():
                if pd.notna(row['array_sizes']) and row['array_sizes'] != '-':
                    for s in str(row['array_sizes']).split(','):
                        try:
                            sizes.append(int(s.strip()))
                        except ValueError:
                            pass
            if sizes:
                plot_data.append(sizes)
                labels.append(superpop)
                colors.append(POP_COLORS.get(superpop, COLORS['dark2']))

        if plot_data:
            bp = ax.boxplot(plot_data, labels=labels, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

        ax.set_title(f'{chrom} Array Sizes by Population', color=text_color, fontsize=12)
        ax.set_xlabel('Superpopulation', color=text_color)
        ax.set_ylabel('Repeat Units', color=text_color)
        ax.tick_params(colors=text_color)
        ax.grid(True, alpha=0.3, color=grid_color)

    plt.tight_layout()

    output_path = output_dir / 'population_distributions.png'
    fig.savefig(output_path, dpi=150, facecolor=background, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved population distributions plot to {output_path}")


def plot_asymmetry_distribution(asymmetry_df, output_dir, background='white'):
    """Plot distribution of haplotype asymmetry."""
    text_color, grid_color = setup_plot_style(background)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[idx]
        chrom_data = asymmetry_df[asymmetry_df['chromosome'] == chrom]

        if not chrom_data.empty:
            ax.hist(chrom_data['asymmetry_ru'], bins=20,
                   color=COLORS['accent2'], alpha=0.7, edgecolor='white')

            # Mark high asymmetry threshold
            ax.axvline(x=20, color=COLORS['accent3'], linestyle='--',
                      label='High asymmetry threshold (20 RU)')

        ax.set_title(f'{chrom} Haplotype Asymmetry', color=text_color, fontsize=12)
        ax.set_xlabel('|hap1 - hap2| (Repeat Units)', color=text_color)
        ax.set_ylabel('Count', color=text_color)
        ax.tick_params(colors=text_color)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, color=grid_color)

    plt.tight_layout()

    output_path = output_dir / 'asymmetry_distribution.png'
    fig.savefig(output_path, dpi=150, facecolor=background, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved asymmetry distribution plot to {output_path}")


def plot_trio_inheritance(trio_df, output_dir, background='white'):
    """Plot trio inheritance analysis results."""
    if trio_df.empty:
        logging.warning("No trio data available for plotting")
        return

    text_color, grid_color = setup_plot_style(background)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[idx]
        chrom_data = trio_df[trio_df['chromosome'] == chrom]

        if chrom_data.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                   transform=ax.transAxes, color=text_color)
            ax.set_title(f'{chrom} Trio Inheritance', color=text_color)
            continue

        # Scatter plot: child size vs distance to nearest parent
        min_dist = chrom_data[['min_distance_to_father', 'min_distance_to_mother']].min(axis=1)

        colors = [COLORS['accent1'] if x else COLORS['accent3']
                 for x in chrom_data['within_parental_range']]

        ax.scatter(chrom_data['child_array_size'], min_dist,
                  c=colors, alpha=0.7, s=50)

        ax.axhline(y=5, color=COLORS['accent4'], linestyle='--',
                  label='De novo threshold (5 RU)')

        ax.set_title(f'{chrom} Trio Inheritance', color=text_color, fontsize=12)
        ax.set_xlabel('Child Array Size (RU)', color=text_color)
        ax.set_ylabel('Min Distance to Parent (RU)', color=text_color)
        ax.tick_params(colors=text_color)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, color=grid_color)

    plt.tight_layout()

    output_path = output_dir / 'trio_inheritance.png'
    fig.savefig(output_path, dpi=150, facecolor=background, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved trio inheritance plot to {output_path}")


def plot_remodeling_summary(exchange_df, fragmented_df, asymmetry_df, output_dir, background='white'):
    """Create summary visualization of remodeling signatures."""
    text_color, grid_color = setup_plot_style(background)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. 4q-10q exchange pie chart
    ax = axes[0]
    if not exchange_df.empty:
        exchange_count = exchange_df['potential_4q10q_exchange'].sum()
        normal_count = len(exchange_df) - exchange_count
        ax.pie([normal_count, exchange_count],
               labels=['Normal', 'Potential Exchange'],
               colors=[COLORS['accent1'], COLORS['accent3']],
               autopct='%1.1f%%', startangle=90)
    ax.set_title('4q-10q Exchange Signatures', color=text_color, fontsize=12)

    # 2. Array fragmentation bar chart
    ax = axes[1]
    if not fragmented_df.empty:
        frag_counts = fragmented_df['num_arrays'].value_counts().sort_index()
        ax.bar(frag_counts.index, frag_counts.values,
               color=COLORS['accent2'], alpha=0.7)
        ax.set_xlabel('Number of Array Fragments', color=text_color)
        ax.set_ylabel('Count', color=text_color)
    ax.set_title('Array Fragmentation', color=text_color, fontsize=12)
    ax.tick_params(colors=text_color)

    # 3. Asymmetry distribution
    ax = axes[2]
    if not asymmetry_df.empty:
        high_asym = asymmetry_df['high_asymmetry'].sum()
        normal_asym = len(asymmetry_df) - high_asym
        ax.pie([normal_asym, high_asym],
               labels=['Normal', 'High Asymmetry'],
               colors=[COLORS['accent1'], COLORS['accent4']],
               autopct='%1.1f%%', startangle=90)
    ax.set_title('Haplotype Asymmetry', color=text_color, fontsize=12)

    plt.tight_layout()

    output_path = output_dir / 'remodeling_summary.png'
    fig.savefig(output_path, dpi=150, facecolor=background, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved remodeling summary plot to {output_path}")


# =============================================================================
# Report Generation
# =============================================================================

def generate_report(trio_df, pop_df, exchange_df, fragmented_df, asymmetry_df, output_dir):
    """Generate summary report as TSV and markdown."""

    # Master summary
    summary = {
        'metric': [],
        'value': [],
        'notes': []
    }

    # Trio metrics
    if not trio_df.empty:
        summary['metric'].append('Trios analyzed')
        summary['value'].append(trio_df['child'].nunique())
        summary['notes'].append('Complete trios with both parents in HPRC')

        summary['metric'].append('Mendelian inheritance rate')
        summary['value'].append(f"{trio_df['within_parental_range'].mean() * 100:.1f}%")
        summary['notes'].append('Child within parental range (±2 RU)')

        summary['metric'].append('Potential de novo events')
        summary['value'].append(trio_df['potential_de_novo'].sum())
        summary['notes'].append('Child >5 RU from both parents')

    # Population metrics
    if not pop_df.empty:
        summary['metric'].append('Populations analyzed')
        summary['value'].append(pop_df['superpopulation'].nunique())
        summary['notes'].append('Superpopulations with data')

    # Remodeling metrics
    if not exchange_df.empty:
        summary['metric'].append('4q-10q exchange signatures')
        summary['value'].append(f"{exchange_df['potential_4q10q_exchange'].mean() * 100:.1f}%")
        summary['notes'].append('Samples with unexpected C+/C- patterns')

    if not fragmented_df.empty:
        summary['metric'].append('Fragmented arrays')
        summary['value'].append(len(fragmented_df))
        summary['notes'].append('Haplotypes with >1 array')

    if not asymmetry_df.empty:
        summary['metric'].append('High haplotype asymmetry')
        summary['value'].append(f"{asymmetry_df['high_asymmetry'].mean() * 100:.1f}%")
        summary['notes'].append('>20 RU or >3-fold difference between haplotypes')

    summary_df = pd.DataFrame(summary)
    summary_path = output_dir / 'd4z4_variability_report.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)
    logging.info(f"Saved summary report to {summary_path}")

    return summary_df


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Analyze D4Z4 array variability in HPRC samples'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='results/d4z4_variability',
        help='Output directory for results (default: results/d4z4_variability)'
    )
    parser.add_argument(
        '--background', '-b',
        choices=['white', 'black'],
        default='white',
        help='Plot background color (default: white)'
    )
    parser.add_argument(
        '--d4z4-file',
        default=D4Z4_SUMMARY_FILE,
        help='Path to D4Z4 summary file'
    )

    args = parser.parse_args()

    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logger = setup_logging(output_dir)
    logger.info("Starting D4Z4 variability analysis")
    logger.info(f"Output directory: {output_dir}")

    # Load data
    d4z4_df = load_d4z4_data(args.d4z4_file)
    pedigree_df = load_pedigree()

    # Get unique samples
    hprc_samples = d4z4_df['sample'].unique().tolist()
    logger.info(f"HPRC samples: {len(hprc_samples)}")

    # Identify trios
    trios = identify_trios(pedigree_df, hprc_samples)

    # Run analyses
    trio_df = analyze_trio_inheritance(d4z4_df, trios)
    pop_df, d4z4_merged = analyze_population_stratification(d4z4_df, pedigree_df)
    exchange_df = detect_4q10q_exchanges(d4z4_df)
    fragmented_df = detect_fragmented_arrays(d4z4_df)
    asymmetry_df = analyze_asymmetry(d4z4_df)

    # Save detailed results
    if not trio_df.empty:
        trio_df.to_csv(output_dir / 'trio_inheritance_analysis.tsv', sep='\t', index=False)

    if not pop_df.empty:
        pop_df.to_csv(output_dir / 'population_array_distributions.tsv', sep='\t', index=False)

    if not exchange_df.empty:
        exchange_df.to_csv(output_dir / 'exchange_signatures.tsv', sep='\t', index=False)

    if not fragmented_df.empty:
        fragmented_df.to_csv(output_dir / 'fragmented_arrays.tsv', sep='\t', index=False)

    if not asymmetry_df.empty:
        asymmetry_df.to_csv(output_dir / 'haplotype_asymmetry.tsv', sep='\t', index=False)

    # Generate visualizations
    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)

    plot_population_distributions(d4z4_merged, figures_dir, args.background)
    plot_asymmetry_distribution(asymmetry_df, figures_dir, args.background)
    plot_trio_inheritance(trio_df, figures_dir, args.background)
    plot_remodeling_summary(exchange_df, fragmented_df, asymmetry_df, figures_dir, args.background)

    # Generate summary report
    summary_df = generate_report(trio_df, pop_df, exchange_df, fragmented_df, asymmetry_df, output_dir)

    # Print summary to console
    print("\n" + "="*60)
    print("D4Z4 VARIABILITY ANALYSIS SUMMARY")
    print("="*60)
    for _, row in summary_df.iterrows():
        print(f"{row['metric']}: {row['value']}")
    print("="*60)

    logger.info("Analysis complete!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
