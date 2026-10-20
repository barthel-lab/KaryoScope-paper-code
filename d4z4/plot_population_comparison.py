#!/usr/bin/env python3
"""
D4Z4 Population Comparison Visualizations

Generates multiple visualizations showing D4Z4 array structure differences
across HPRC populations:
1. Multi-sunburst grid (one per population)
2. Stacked bar chart (structure types by population)
3. Alluvial/Sankey diagram (population → type → repeat range)
4. Ridge plot (repeat count distributions)

Usage:
    python plot_population_comparison.py --viz all --background black
    python plot_population_comparison.py --viz sunburst --output-dir results/
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from matplotlib.patches import Wedge, Rectangle, FancyBboxPatch
from matplotlib.sankey import Sankey
from scipy import stats

# =============================================================================
# Configuration
# =============================================================================

PEDIGREE_URL = "https://raw.githubusercontent.com/human-pangenomics/HPRC_metadata/main/data/IGSR/20130606_g1k.ped"
D4Z4_ANALYSIS_FILE = "results/d4z4_array_analysis.tsv"

# Brand colors
COLORS = {
    'dark1': '#000000',
    'light1': '#FFFFFF',
    'dark2': '#545454',
    'light2': '#C4A9E8',
    'accent1': '#40D392',
    'accent2': '#60A5FA',
    'accent3': '#F07167',
    'accent4': '#FBBF24',
    'accent5': '#10B981',
    'accent6': '#3B82F6',
}

# Population colors
POP_COLORS = {
    'AFR': '#F07167',  # Coral
    'EUR': '#60A5FA',  # Blue
    'EAS': '#40D392',  # Green
    'AMR': '#FBBF24',  # Yellow
    'SAS': '#C4A9E8',  # Lavender
}

POP_ORDER = ['AFR', 'EUR', 'EAS', 'AMR', 'SAS']

# Superpopulation mapping
SUPERPOP_MAP = {
    'ACB': 'AFR', 'ASW': 'AFR', 'ESN': 'AFR', 'GWD': 'AFR', 'LWK': 'AFR', 'MSL': 'AFR', 'YRI': 'AFR',
    'CDX': 'EAS', 'CHB': 'EAS', 'CHS': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS',
    'CEU': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 'TSI': 'EUR',
    'CLM': 'AMR', 'MXL': 'AMR', 'PEL': 'AMR', 'PUR': 'AMR',
    'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS',
}

# Structure type colors (dark mode)
STRUCTURE_COLORS = {
    'qA': '#34D399',      # Emerald - C+ with β-satellite
    'qB': '#60A5FA',      # Blue - C- no β-satellite
    'D': '#F87171',       # Coral - Degraded
    'No arrays': '#6B7280',
    '2+ arrays': '#C4B5FD',
}

# Repeat range colors
RANGE_COLORS = {
    '1-10': '#F87171',    # Coral - FSHD1
    '11-20': '#FBBF24',   # Amber - FSHD2
    '21-30': '#34D399',   # Green - Normal
    '31+': '#60A5FA',     # Blue - Long
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
# Data Loading
# =============================================================================

def load_pedigree():
    """Fetch and parse the 1000 Genomes pedigree file."""
    logging.info(f"Fetching pedigree from {PEDIGREE_URL}")
    response = requests.get(PEDIGREE_URL)
    response.raise_for_status()

    lines = response.text.strip().split('\n')
    records = []

    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 7:
            records.append({
                'individual_id': parts[1],
                'population': parts[6] if len(parts) > 6 else None,
            })

    df = pd.DataFrame(records)
    df['superpopulation'] = df['population'].map(SUPERPOP_MAP)
    return df


def load_d4z4_data(filepath=D4Z4_ANALYSIS_FILE):
    """Load D4Z4 array analysis data."""
    logging.info(f"Loading D4Z4 data from {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    return df


def merge_population_data(d4z4_df, pedigree_df):
    """Merge population info into D4Z4 data."""
    pop_map = pedigree_df.set_index('individual_id')['superpopulation'].to_dict()
    d4z4_df = d4z4_df.copy()
    d4z4_df['superpopulation'] = d4z4_df['sample'].map(pop_map)
    return d4z4_df


def get_array_type(row):
    """Determine array type for single-array samples."""
    if row['num_arrays'] != 1:
        return None

    if pd.isna(row.get('arrays_info')):
        if row['terminal_feature'] == 'β-satellite terminal':
            return 'C+'
        elif row['terminal_feature'] == 'No β-satellite terminal':
            return 'C-'
        else:
            return 'D'

    try:
        arrays = json.loads(row['arrays_info'])
        if len(arrays) == 1:
            arr = arrays[0]
            if arr['type'] == 'canonical':
                return 'C+' if arr.get('has_plam', False) else 'C-'
            else:
                return 'D'
    except:
        pass
    return 'D'


def get_repeat_range(count):
    """Categorize repeat count into clinical ranges."""
    if pd.isna(count) or count == 0:
        return None
    if count <= 10:
        return '1-10'
    elif count <= 20:
        return '11-20'
    elif count <= 30:
        return '21-30'
    else:
        return '31+'


def calculate_distributions(df, chromosome, population=None, plam_df=None):
    """Calculate distributions for sunburst chart."""
    chrom_df = df[df['chromosome'] == chromosome].copy()
    if population:
        chrom_df = chrom_df[chrom_df['superpopulation'] == population]

    total = len(chrom_df)
    if total == 0:
        return None

    # Level 1: Number of arrays
    no_arrays = len(chrom_df[chrom_df['num_arrays'] == 0])
    single_array = len(chrom_df[chrom_df['num_arrays'] == 1])
    multi_array = len(chrom_df[chrom_df['num_arrays'] > 1])

    level1 = {
        'No arrays': no_arrays,
        '1 array': single_array,
        '2+ arrays': multi_array,
    }

    # Level 2: Array type for single-array samples
    single_df = chrom_df[chrom_df['num_arrays'] == 1].copy()
    single_df['array_type'] = single_df.apply(get_array_type, axis=1)

    level2 = {
        'C+': len(single_df[single_df['array_type'] == 'C+']),
        'C-': len(single_df[single_df['array_type'] == 'C-']),
        'D': len(single_df[single_df['array_type'] == 'D']),
    }

    # Level 3: pLAM status for C+ (qA) samples
    cplus_df = single_df[single_df['array_type'] == 'C+'].copy()
    cplus_df['plam_status'] = 'Absent'
    if plam_df is not None:
        chrom_plam = plam_df[plam_df['chromosome'] == chromosome]
        for idx, row in cplus_df.iterrows():
            match = chrom_plam[(chrom_plam['sample'] == row['sample']) &
                               (chrom_plam['hap'] == row['hap_num'])]
            if len(match) > 0:
                status = match.iloc[0]['plam_functional']
                if status == 'func_Plam':
                    cplus_df.at[idx, 'plam_status'] = 'Functional'
                elif status == 'nonfunc_Plam':
                    cplus_df.at[idx, 'plam_status'] = 'Non-functional'

    level3 = {
        'Functional': len(cplus_df[cplus_df['plam_status'] == 'Functional']),
        'Non-functional': len(cplus_df[cplus_df['plam_status'] == 'Non-functional']),
        'Absent': len(cplus_df[cplus_df['plam_status'] == 'Absent']),
    }

    # Level 4: Repeat count ranges for Functional pLAM only
    func_df = cplus_df[cplus_df['plam_status'] == 'Functional'].copy()
    func_df['repeat_range'] = func_df['repeat_count_main'].apply(get_repeat_range)

    level4 = {
        '1-10': len(func_df[func_df['repeat_range'] == '1-10']),
        '11-20': len(func_df[func_df['repeat_range'] == '11-20']),
        '21-30': len(func_df[func_df['repeat_range'] == '21-30']),
        '31+': len(func_df[func_df['repeat_range'] == '31+']),
    }

    return {
        'total': total,
        'level1': level1,
        'level2': level2,
        'level3': level3,
        'level4': level4,
    }


# =============================================================================
# Visualization 1: Multi-Sunburst Grid
# =============================================================================

def plot_sunburst_mini(ax, data, title, chromosome='chr4', text_color='white'):
    """Plot a compact sunburst for grid layout."""
    if data is None:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, color=text_color, fontsize=10)
        ax.set_title(title, fontsize=10, color=text_color)
        ax.axis('off')
        return

    # Ring dimensions (4 rings, compact for grid)
    inner_radius = 0.20
    middle_radius = 0.38
    plam_radius = 0.56
    outer_radius = 0.74
    ring_width = 0.16

    level1 = data['level1']
    level2 = data['level2']
    level3 = data['level3']
    level4 = data['level4']
    total = data['total']

    chrom_prefix = '4q' if chromosome == 'chr4' else '10q'
    qA_label = f'{chrom_prefix}A'
    qB_label = f'{chrom_prefix}B'

    level1_colors = {
        'No arrays': '#6B7280',
        '1 array': '#6EE7B7',
        '2+ arrays': '#C4B5FD',
    }
    level2_colors = {
        qA_label: '#34D399',
        qB_label: '#60A5FA',
        'D': '#F87171',
    }
    plam_colors = {
        'Functional': '#40D392',
        'Non-functional': '#F07167',
        'Absent': '#6B7280',
    }
    level4_colors = {
        '1-10': '#F87171',
        '11-20': '#FBBF24',
        '21-30': '#34D399',
        '31+': '#60A5FA',
    }

    level2_mapped = {
        qA_label: level2.get('C+', 0),
        qB_label: level2.get('C-', 0),
        'D': level2.get('D', 0),
    }

    # Draw level 1
    level1_order = ['No arrays', '1 array', '2+ arrays']
    level1_values = [level1.get(k, 0) for k in level1_order]
    level1_total = sum(level1_values)

    start_angle = 90
    level1_angles = []

    for label, value in zip(level1_order, level1_values):
        if value == 0:
            level1_angles.append((start_angle, start_angle, label, value))
            continue
        angle_span = 360 * value / level1_total
        end_angle = start_angle - angle_span
        level1_angles.append((start_angle, end_angle, label, value))
        start_angle = end_angle

    for start, end, label, value in level1_angles:
        if value == 0:
            continue
        color = level1_colors.get(label, '#CCCCCC')
        wedge = Wedge((0, 0), inner_radius + ring_width, end, start,
                      width=ring_width, facecolor=color, edgecolor=text_color, linewidth=0.8)
        ax.add_patch(wedge)

    # Draw level 2 (for 1 array slice)
    single_array_angles = None
    for start, end, label, value in level1_angles:
        if label == '1 array' and value > 0:
            single_array_angles = (start, end)
            break

    if single_array_angles:
        l1_start, l1_end = single_array_angles
        l1_span = l1_start - l1_end

        level2_order = [qA_label, qB_label, 'D']
        level2_values = [level2_mapped.get(k, 0) for k in level2_order]
        level2_total = sum(level2_values)

        l2_start = l1_start
        level2_angles = []

        for label, value in zip(level2_order, level2_values):
            if value == 0 or level2_total == 0:
                level2_angles.append((l2_start, l2_start, label, value))
                continue
            angle_span = l1_span * value / level2_total
            l2_end = l2_start - angle_span

            color = level2_colors.get(label, '#CCCCCC')
            wedge = Wedge((0, 0), middle_radius + ring_width, l2_end, l2_start,
                          width=ring_width, facecolor=color, edgecolor=text_color, linewidth=0.8)
            ax.add_patch(wedge)

            level2_angles.append((l2_start, l2_end, label, value))
            l2_start = l2_end

        # Draw level 3: pLAM status ring (for qA slice)
        qA_angles = None
        for start, end, label, value in level2_angles:
            if label == qA_label and value > 0:
                qA_angles = (start, end)
                break

        if qA_angles:
            l2_start_qa, l2_end_qa = qA_angles
            l2_span_qa = l2_start_qa - l2_end_qa

            level3_order = ['Functional', 'Non-functional', 'Absent']
            level3_values = [level3.get(k, 0) for k in level3_order]
            level3_total = sum(level3_values)

            l3_start = l2_start_qa
            level3_angles = []

            for label, value in zip(level3_order, level3_values):
                if value == 0 or level3_total == 0:
                    level3_angles.append((l3_start, l3_start, label, value))
                    continue
                angle_span = l2_span_qa * value / level3_total
                l3_end = l3_start - angle_span

                color = plam_colors.get(label, '#CCCCCC')
                wedge = Wedge((0, 0), plam_radius + ring_width, l3_end, l3_start,
                              width=ring_width, facecolor=color, edgecolor=text_color, linewidth=0.8)
                ax.add_patch(wedge)

                level3_angles.append((l3_start, l3_end, label, value))
                l3_start = l3_end

            # Draw level 4: repeat ranges (for Functional pLAM only)
            func_angles = None
            for start, end, label, value in level3_angles:
                if label == 'Functional' and value > 0:
                    func_angles = (start, end)
                    break

            if func_angles:
                l3_start_f, l3_end_f = func_angles
                l3_span_f = l3_start_f - l3_end_f

                level4_order = ['1-10', '11-20', '21-30', '31+']
                level4_values = [level4.get(k, 0) for k in level4_order]
                level4_total = sum(level4_values)

                l4_start = l3_start_f
                for label, value in zip(level4_order, level4_values):
                    if value == 0 or level4_total == 0:
                        continue
                    angle_span = l3_span_f * value / level4_total
                    l4_end = l4_start - angle_span

                    color = level4_colors.get(label, '#CCCCCC')
                    wedge = Wedge((0, 0), outer_radius + ring_width, l4_end, l4_start,
                                  width=ring_width, facecolor=color, edgecolor=text_color, linewidth=0.8)
                    ax.add_patch(wedge)
                    l4_start = l4_end

    # Center text
    ax.text(0, 0, f"n={total}", ha='center', va='center', fontsize=8, fontweight='bold', color=text_color)

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=10, fontweight='bold', pad=5, color=text_color)


def plot_sunburst_grid_legend(ax, text_color='white'):
    """Draw a compact legend for the sunburst grid."""
    legend_items = [
        ('Inner Ring - Array Count:', None, 'header'),
        ('No arrays', '#6B7280', 'item'),
        ('1 array', '#6EE7B7', 'item'),
        ('2+ arrays', '#C4B5FD', 'item'),
        ('', None, 'spacer'),
        ('Middle Ring - Array Type:', None, 'header'),
        ('qA (with β-satellite)', '#34D399', 'item'),
        ('qB (no β-satellite)', '#60A5FA', 'item'),
        ('D (degraded)', '#F87171', 'item'),
        ('', None, 'spacer'),
        ('pLAM Ring (qA only):', None, 'header'),
        ('Functional pLAM', '#40D392', 'item'),
        ('Non-functional pLAM', '#F07167', 'item'),
        ('Absent', '#6B7280', 'item'),
        ('', None, 'spacer'),
        ('Outer Ring - Repeat Count', None, 'header'),
        ('(Functional pLAM only):', None, 'header'),
        ('1-10 (FSHD1)', '#F87171', 'item'),
        ('11-20 (FSHD2)', '#FBBF24', 'item'),
        ('21-30 (normal)', '#34D399', 'item'),
        ('31+ (long)', '#60A5FA', 'item'),
    ]

    y = 0.95
    for label, color, item_type in legend_items:
        if item_type == 'header':
            ax.text(0.05, y, label, fontsize=8, fontweight='bold', transform=ax.transAxes, color=text_color)
            y -= 0.05
        elif item_type == 'item':
            ax.add_patch(plt.Rectangle((0.05, y - 0.015), 0.04, 0.03,
                                        facecolor=color, edgecolor='gray',
                                        linewidth=0.5, transform=ax.transAxes))
            ax.text(0.12, y, label, fontsize=7, transform=ax.transAxes, va='center', color=text_color)
            y -= 0.045
        elif item_type == 'spacer':
            y -= 0.02

    ax.axis('off')


def plot_sunburst_grid(df, output_path, background='black', plam_df=None):
    """Create 2x6 grid of sunbursts per population with legend."""
    logging.info("Generating multi-sunburst grid...")

    bg_color = '#1a1a1a' if background == 'black' else 'white'
    text_color = '#FFFFFF' if background == 'black' else '#000000'

    # Create figure with space for legend on the right
    fig = plt.figure(figsize=(24, 10), facecolor=bg_color)

    # Grid axes (2 rows x 6 cols) on the left 80%
    axes = []
    for row_idx in range(2):
        row_axes = []
        for col_idx in range(6):
            left = 0.02 + col_idx * 0.13
            bottom = 0.52 - row_idx * 0.48 + 0.05
            ax = fig.add_axes([left, bottom, 0.12, 0.42], facecolor=bg_color)
            row_axes.append(ax)
        axes.append(row_axes)

    # Legend axis on the right
    ax_legend = fig.add_axes([0.82, 0.05, 0.17, 0.88], facecolor=bg_color)

    for row_idx, chrom in enumerate(['chr4', 'chr10']):
        # First column: All populations combined
        ax = axes[row_idx][0]
        data = calculate_distributions(df, chrom, plam_df=plam_df)
        plot_sunburst_mini(ax, data, f"{chrom} - All", chromosome=chrom, text_color=text_color)

        # Remaining columns: Per population
        for col_idx, pop in enumerate(POP_ORDER):
            ax = axes[row_idx][col_idx + 1]
            data = calculate_distributions(df, chrom, population=pop, plam_df=plam_df)
            title = f"{chrom} - {pop}"
            plot_sunburst_mini(ax, data, title, chromosome=chrom, text_color=text_color)
            # Add population color indicator
            ax.add_patch(Rectangle((-0.95, 0.85), 0.15, 0.1,
                                   facecolor=POP_COLORS.get(pop, '#FFFFFF'), edgecolor='none',
                                   transform=ax.transAxes, clip_on=False))

    # Draw legend
    plot_sunburst_grid_legend(ax_legend, text_color=text_color)

    fig.suptitle('D4Z4 Array Structure by Population', fontsize=14, fontweight='bold',
                 y=0.98, color=text_color)

    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.savefig(output_path.replace('.svg', '.png'), dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    logging.info(f"Saved: {output_path}")


# =============================================================================
# Visualization 2: Stacked Bar Chart
# =============================================================================

def plot_stacked_bar(df, output_path, background='black'):
    """Create stacked bar chart of structure types by population."""
    logging.info("Generating stacked bar chart...")

    bg_color = '#000000' if background == 'black' else 'white'
    text_color = '#FFFFFF' if background == 'black' else '#000000'

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor=bg_color)

    categories = ['No arrays', 'qA (C+)', 'qB (C-)', 'Degraded', '2+ arrays']
    cat_colors = ['#6B7280', '#34D399', '#60A5FA', '#F87171', '#C4B5FD']

    for ax_idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[ax_idx]
        ax.set_facecolor(bg_color)

        chrom_df = df[df['chromosome'] == chrom].copy()
        chrom_df['array_type'] = chrom_df.apply(get_array_type, axis=1)

        data = {cat: [] for cat in categories}
        totals = []

        for pop in POP_ORDER:
            pop_df = chrom_df[chrom_df['superpopulation'] == pop]
            total = len(pop_df)
            totals.append(total)

            if total == 0:
                for cat in categories:
                    data[cat].append(0)
                continue

            no_arrays = len(pop_df[pop_df['num_arrays'] == 0]) / total * 100
            multi = len(pop_df[pop_df['num_arrays'] > 1]) / total * 100

            single_df = pop_df[pop_df['num_arrays'] == 1]
            cplus = len(single_df[single_df['array_type'] == 'C+']) / total * 100
            cminus = len(single_df[single_df['array_type'] == 'C-']) / total * 100
            degraded = len(single_df[single_df['array_type'] == 'D']) / total * 100

            data['No arrays'].append(no_arrays)
            data['qA (C+)'].append(cplus)
            data['qB (C-)'].append(cminus)
            data['Degraded'].append(degraded)
            data['2+ arrays'].append(multi)

        x = np.arange(len(POP_ORDER))
        width = 0.6

        bottom = np.zeros(len(POP_ORDER))
        for cat, color in zip(categories, cat_colors):
            values = data[cat]
            bars = ax.bar(x, values, width, bottom=bottom, label=cat, color=color, edgecolor=text_color, linewidth=0.5)
            bottom += values

        # Add sample counts above bars
        for i, (pop, total) in enumerate(zip(POP_ORDER, totals)):
            ax.text(i, 102, f'n={total}', ha='center', va='bottom', fontsize=8, color=text_color)

        ax.set_ylabel('Percentage of Haplotypes', color=text_color, fontsize=11)
        ax.set_xlabel('Superpopulation', color=text_color, fontsize=11)
        ax.set_title(f'{chrom}', color=text_color, fontsize=12, fontweight='bold', pad=12)
        ax.set_xticks(x)
        ax.set_xticklabels(POP_ORDER, fontsize=10, fontweight='bold', color=text_color)
        ax.tick_params(colors=text_color)
        ax.set_ylim(0, 115)
        ax.set_yticks([0, 20, 40, 60, 80, 100])  # Only show ticks up to 100

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_color(text_color)
        ax.spines['left'].set_color(text_color)
        ax.spines['left'].set_bounds(0, 100)  # Stop y-axis line at 100

    # Legend
    axes[1].legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False,
                   labelcolor=text_color, fontsize=9)

    fig.suptitle('D4Z4 Array Structure Types by Population', fontsize=14, fontweight='bold',
                 y=0.98, color=text_color)

    plt.tight_layout(rect=[0, 0, 0.92, 0.95])
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.savefig(output_path.replace('.svg', '.png'), dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    logging.info(f"Saved: {output_path}")


# =============================================================================
# Visualization 3: Alluvial/Sankey Diagram
# =============================================================================

def plot_alluvial(df, output_path, background='black'):
    """Create alluvial diagram: Population → Array Type → Repeat Range."""
    logging.info("Generating alluvial diagram...")

    bg_color = '#1a1a1a' if background == 'black' else 'white'
    text_color = '#FFFFFF' if background == 'black' else '#000000'

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor=bg_color)

    for ax_idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[ax_idx]
        ax.set_facecolor(bg_color)

        chrom_df = df[df['chromosome'] == chrom].copy()
        chrom_df['array_type'] = chrom_df.apply(get_array_type, axis=1)
        chrom_df['repeat_range'] = chrom_df['repeat_count_main'].apply(get_repeat_range)

        # Build flow data
        # Node positions: Left (populations), Middle (types), Right (ranges)
        pop_nodes = {pop: i for i, pop in enumerate(POP_ORDER)}
        type_nodes = {'C+': 0, 'C-': 1, 'D': 2, 'None': 3}
        range_nodes = {'1-10': 0, '11-20': 1, '21-30': 2, '31+': 3, 'None': 4}

        # Calculate flows
        # Pop -> Type flows
        pop_type_flows = {}
        for pop in POP_ORDER:
            pop_df = chrom_df[chrom_df['superpopulation'] == pop]
            for atype in ['C+', 'C-', 'D']:
                count = len(pop_df[pop_df['array_type'] == atype])
                if count > 0:
                    pop_type_flows[(pop, atype)] = count
            # Also count no arrays + multi as "None" type
            none_count = len(pop_df[pop_df['num_arrays'] != 1])
            if none_count > 0:
                pop_type_flows[(pop, 'None')] = none_count

        # Type -> Range flows (only for C+)
        type_range_flows = {}
        cplus_df = chrom_df[chrom_df['array_type'] == 'C+']
        for rng in ['1-10', '11-20', '21-30', '31+']:
            count = len(cplus_df[cplus_df['repeat_range'] == rng])
            if count > 0:
                type_range_flows[('C+', rng)] = count

        # Draw as curved flow bands
        # Left column: Populations
        left_x = 0
        mid_x = 1
        right_x = 2

        # Calculate vertical positions
        pop_heights = {}
        pop_y_starts = {}
        y_pos = 0
        total_left = sum(pop_type_flows.values())

        for pop in POP_ORDER:
            pop_total = sum(v for (p, t), v in pop_type_flows.items() if p == pop)
            pop_heights[pop] = pop_total / total_left * 10
            pop_y_starts[pop] = y_pos
            y_pos += pop_heights[pop] + 0.2

        type_heights = {'C+': 0, 'C-': 0, 'D': 0, 'None': 0}
        for (p, t), v in pop_type_flows.items():
            type_heights[t] += v

        type_total = sum(type_heights.values())
        type_y_starts = {}
        y_pos = 0
        for t in ['C+', 'C-', 'D', 'None']:
            type_y_starts[t] = y_pos
            type_heights[t] = type_heights[t] / total_left * 10
            y_pos += type_heights[t] + 0.15

        # Draw population blocks
        for pop in POP_ORDER:
            y = pop_y_starts[pop]
            h = pop_heights[pop]
            rect = FancyBboxPatch((left_x - 0.15, y), 0.3, h,
                                  boxstyle="round,pad=0.02",
                                  facecolor=POP_COLORS[pop], edgecolor=text_color, linewidth=1)
            ax.add_patch(rect)
            ax.text(left_x - 0.25, y + h/2, pop, ha='right', va='center',
                    fontsize=9, fontweight='bold', color=text_color)

        # Draw type blocks
        type_colors = {'C+': '#34D399', 'C-': '#60A5FA', 'D': '#F87171', 'None': '#6B7280'}
        type_labels = {'C+': 'qA', 'C-': 'qB', 'D': 'D', 'None': '0/2+'}
        for t in ['C+', 'C-', 'D', 'None']:
            if type_heights[t] > 0:
                y = type_y_starts[t]
                h = type_heights[t]
                rect = FancyBboxPatch((mid_x - 0.15, y), 0.3, h,
                                      boxstyle="round,pad=0.02",
                                      facecolor=type_colors[t], edgecolor=text_color, linewidth=1)
                ax.add_patch(rect)
                ax.text(mid_x, y + h/2, type_labels[t], ha='center', va='center',
                        fontsize=9, fontweight='bold', color='black')

        # Draw range blocks for C+ only
        range_heights = {'1-10': 0, '11-20': 0, '21-30': 0, '31+': 0}
        for (t, r), v in type_range_flows.items():
            range_heights[r] += v

        range_total = sum(range_heights.values())
        range_y_starts = {}
        y_pos = type_y_starts['C+']  # Align with C+ block
        for r in ['1-10', '11-20', '21-30', '31+']:
            range_y_starts[r] = y_pos
            range_heights[r] = range_heights[r] / total_left * 10 if range_total > 0 else 0
            y_pos += range_heights[r] + 0.1

        range_colors = {'1-10': '#F87171', '11-20': '#FBBF24', '21-30': '#34D399', '31+': '#60A5FA'}
        range_labels = {'1-10': 'FSHD1', '11-20': 'FSHD2', '21-30': 'Normal', '31+': 'Long'}
        for r in ['1-10', '11-20', '21-30', '31+']:
            if range_heights[r] > 0:
                y = range_y_starts[r]
                h = range_heights[r]
                rect = FancyBboxPatch((right_x - 0.15, y), 0.3, h,
                                      boxstyle="round,pad=0.02",
                                      facecolor=range_colors[r], edgecolor=text_color, linewidth=1)
                ax.add_patch(rect)
                ax.text(right_x + 0.25, y + h/2, range_labels[r], ha='left', va='center',
                        fontsize=8, fontweight='bold', color=text_color)

        # Draw flow bands (simplified bezier curves)
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path

        # Pop -> Type flows
        pop_y_offsets = {pop: 0 for pop in POP_ORDER}
        type_y_offsets = {t: 0 for t in type_heights}

        for (pop, atype), count in sorted(pop_type_flows.items(), key=lambda x: x[1], reverse=True):
            flow_height = count / total_left * 10
            y1 = pop_y_starts[pop] + pop_y_offsets[pop]
            y2 = type_y_starts[atype] + type_y_offsets[atype]

            # Bezier curve
            verts = [
                (left_x + 0.15, y1),
                (left_x + 0.5, y1),
                (mid_x - 0.5, y2),
                (mid_x - 0.15, y2),
                (mid_x - 0.15, y2 + flow_height),
                (mid_x - 0.5, y2 + flow_height),
                (left_x + 0.5, y1 + flow_height),
                (left_x + 0.15, y1 + flow_height),
                (left_x + 0.15, y1),
            ]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]

            path = Path(verts, codes)
            patch = PathPatch(path, facecolor=POP_COLORS[pop], alpha=0.4, edgecolor='none')
            ax.add_patch(patch)

            pop_y_offsets[pop] += flow_height
            type_y_offsets[atype] += flow_height

        # Type -> Range flows
        cplus_y_offset = 0
        range_y_offsets = {r: 0 for r in range_heights}

        for (atype, rng), count in sorted(type_range_flows.items(), key=lambda x: x[1], reverse=True):
            flow_height = count / total_left * 10
            y1 = type_y_starts['C+'] + cplus_y_offset
            y2 = range_y_starts[rng] + range_y_offsets[rng]

            verts = [
                (mid_x + 0.15, y1),
                (mid_x + 0.5, y1),
                (right_x - 0.5, y2),
                (right_x - 0.15, y2),
                (right_x - 0.15, y2 + flow_height),
                (right_x - 0.5, y2 + flow_height),
                (mid_x + 0.5, y1 + flow_height),
                (mid_x + 0.15, y1 + flow_height),
                (mid_x + 0.15, y1),
            ]
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                     Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]

            path = Path(verts, codes)
            patch = PathPatch(path, facecolor=range_colors[rng], alpha=0.4, edgecolor='none')
            ax.add_patch(patch)

            cplus_y_offset += flow_height
            range_y_offsets[rng] += flow_height

        # Column headers
        ax.text(left_x, -0.8, 'Population', ha='center', va='top',
                fontsize=10, fontweight='bold', color=text_color)
        ax.text(mid_x, -0.8, 'Array Type', ha='center', va='top',
                fontsize=10, fontweight='bold', color=text_color)
        ax.text(right_x, -0.8, 'Repeat Range\n(qA only)', ha='center', va='top',
                fontsize=9, fontweight='bold', color=text_color)

        ax.set_xlim(-0.8, 2.8)
        ax.set_ylim(-1.5, 13)
        ax.axis('off')
        ax.set_title(f'{chrom}', color=text_color, fontsize=12, fontweight='bold', pad=10)

    fig.suptitle('D4Z4 Array Structure Flow: Population > Type > Clinical Range',
                 fontsize=14, fontweight='bold', y=0.98, color=text_color)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.savefig(output_path.replace('.svg', '.png'), dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    logging.info(f"Saved: {output_path}")


# =============================================================================
# Visualization 4: Ridge Plot
# =============================================================================

def plot_ridge(df, output_path, background='black'):
    """Create ridge plot of repeat count distributions by population."""
    logging.info("Generating ridge plot...")

    bg_color = '#000000' if background == 'black' else 'white'
    text_color = '#FFFFFF' if background == 'black' else '#000000'

    # Use more horizontal space to prevent overlap
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), facecolor=bg_color)
    plt.subplots_adjust(wspace=0.3)  # Add space between subplots

    for ax_idx, chrom in enumerate(['chr4', 'chr10']):
        ax = axes[ax_idx]
        ax.set_facecolor(bg_color)

        chrom_df = df[(df['chromosome'] == chrom) & (df['num_arrays'] == 1)].copy()

        # Ridge plot parameters - no overlap, each ridge separate
        n_pops = len(POP_ORDER)
        x_range = np.linspace(0, 80, 200)
        ridge_height = 1.0  # Height allocated per population
        max_density_height = 0.85  # Max height of density curve within allocated space

        for i, pop in enumerate(reversed(POP_ORDER)):  # Reverse for top-to-bottom
            pop_df = chrom_df[chrom_df['superpopulation'] == pop]
            values = pop_df['repeat_count_main'].dropna().values

            if len(values) < 3:
                continue

            # KDE
            try:
                kde = stats.gaussian_kde(values, bw_method=0.3)
                density = kde(x_range)
                density = density / density.max() * max_density_height  # Normalize height
            except:
                continue

            y_offset = i * ridge_height
            actual_pop = POP_ORDER[n_pops - 1 - i]
            color = POP_COLORS[actual_pop]

            # Fill
            ax.fill_between(x_range, y_offset, y_offset + density,
                           color=color, alpha=0.7, edgecolor=text_color, linewidth=1)

            # Population label in white - positioned at x=0
            ax.text(-1, y_offset + 0.35, actual_pop,
                    ha='right', va='center', fontsize=11, fontweight='bold', color=text_color)

            # Stats
            median = np.median(values)
            min_val = np.min(values)
            max_val = np.max(values)

            # Stats annotation on right - n= on top, med/range below
            ax.text(78, y_offset + 0.55, f'n={len(values)}',
                    ha='right', va='center', fontsize=8, color=color, alpha=0.9)
            ax.text(78, y_offset + 0.25, f'med={median:.0f}, range={min_val:.0f}-{max_val:.0f}',
                    ha='right', va='center', fontsize=7, color=text_color, alpha=0.7)

        # Clinical thresholds - white lines from x-axis to top
        y_max = n_pops * ridge_height
        for threshold, label in [(10, 'FSHD1'), (20, 'FSHD2')]:
            ax.vlines(x=threshold, ymin=0, ymax=y_max + 0.25,
                     color='white', linewidth=1.5, linestyle=':', alpha=0.6)
            ax.text(threshold, y_max + 0.35, label,
                    ha='center', va='bottom', fontsize=8, color=text_color, fontweight='bold')

        ax.set_xlim(0, 80)
        ax.set_ylim(-0.1, y_max + 0.5)
        ax.set_xlabel('Repeat Units (RU)', color=text_color, fontsize=11)
        ax.set_yticks([])
        ax.set_title(f'{chrom}', color=text_color, fontsize=12, fontweight='bold', pad=10)
        ax.tick_params(colors=text_color)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color(text_color)

    fig.suptitle('D4Z4 Repeat Count Distribution by Population',
                 fontsize=14, fontweight='bold', y=0.98, color=text_color)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.savefig(output_path.replace('.svg', '.png'), dpi=150, bbox_inches='tight', facecolor=bg_color)
    plt.close()
    logging.info(f"Saved: {output_path}")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate D4Z4 population comparison visualizations'
    )
    parser.add_argument(
        '--viz', '-v',
        choices=['sunburst', 'bar', 'alluvial', 'ridge', 'all'],
        default='all',
        help='Visualization type to generate (default: all)'
    )
    parser.add_argument(
        '--background', '-b',
        choices=['white', 'black'],
        default='black',
        help='Background color (default: black)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='results',
        help='Output directory (default: results)'
    )
    parser.add_argument(
        '--analysis-file',
        default=D4Z4_ANALYSIS_FILE,
        help='D4Z4 analysis TSV file'
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logging.info("Loading data...")
    d4z4_df = load_d4z4_data(args.analysis_file)
    pedigree_df = load_pedigree()
    df = merge_population_data(d4z4_df, pedigree_df)

    # Report population coverage
    pop_counts = df[df['superpopulation'].notna()]['superpopulation'].value_counts()
    logging.info(f"Population distribution:\n{pop_counts.to_string()}")

    unknown = df['superpopulation'].isna().sum()
    if unknown > 0:
        logging.warning(f"{unknown} records have unknown population")

    # Load pLAM data
    plam_file = Path(args.analysis_file).parent / 'd4z4_array_coordinates.plam.tsv'
    plam_df = pd.read_csv(plam_file, sep='\t') if plam_file.exists() else None
    if plam_df is not None:
        logging.info(f"Loaded pLAM data: {len(plam_df)} rows from {plam_file}")
    else:
        logging.warning(f"pLAM file not found: {plam_file}")

    # Generate visualizations
    viz_type = args.viz

    if viz_type in ['sunburst', 'all']:
        plot_sunburst_grid(df, str(output_dir / 'population_sunburst_grid.svg'), args.background, plam_df=plam_df)

    if viz_type in ['bar', 'all']:
        plot_stacked_bar(df, str(output_dir / 'population_stacked_bar.svg'), args.background)

    if viz_type in ['alluvial', 'all']:
        plot_alluvial(df, str(output_dir / 'population_alluvial.svg'), args.background)

    if viz_type in ['ridge', 'all']:
        plot_ridge(df, str(output_dir / 'population_ridge.svg'), args.background)

    logging.info("Done!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
