#!/usr/bin/env python3
"""
D4Z4 Array Distribution Pie Charts

Generates hierarchical pie charts showing:
1. Number of arrays (0, 1, 2+)
2. For single-array: Array type (C+, C-, D)
3. For C+: Repeat count ranges (FSHD1, FSHD2, normal, long)

Usage:
    python plot_array_distribution.py
    python plot_array_distribution.py --output results/array_distribution.pdf
"""

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Configure matplotlib for SVG text
plt.rcParams['svg.fonttype'] = 'none'

# Color scheme matching the SVG plots
COLORS = {
    'C+': '#90EE90',      # Light green - canonical with β-satellite
    'C-': '#87CEEB',      # Light blue - canonical without β-satellite
    'D': '#FF6B6B',       # Light red - degraded
    'No arrays': '#CCCCCC',  # Gray
    'Multi-array': '#DDA0DD',  # Plum
    # Repeat count ranges
    '1-10 (FSHD1)': '#FF4444',    # Red - disease range
    '11-20 (FSHD2)': '#FFA500',   # Orange - intermediate
    '21-30': '#90EE90',           # Green - normal
    '31+': '#4444FF',             # Blue - long
}


def load_data(analysis_file):
    """Load and process analysis data."""
    df = pd.read_csv(analysis_file, sep='\t')
    return df


def get_array_type(row):
    """Determine array type for single-array samples."""
    if row['num_arrays'] != 1:
        return None

    if pd.isna(row.get('arrays_info')):
        # Fallback to terminal_feature
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
    if count <= 10:
        return '1-10 (FSHD1)'
    elif count <= 20:
        return '11-20 (FSHD2)'
    elif count <= 30:
        return '21-30'
    else:
        return '31+'


def calculate_distributions(df, chromosome, plam_df=None):
    """Calculate distributions for a chromosome."""
    chrom_df = df[df['chromosome'] == chromosome]
    total = len(chrom_df)

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

    # Match pLAM status from plam coordinates file
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
        '1-10 (FSHD1)': len(func_df[func_df['repeat_range'] == '1-10 (FSHD1)']),
        '11-20 (FSHD2)': len(func_df[func_df['repeat_range'] == '11-20 (FSHD2)']),
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


def make_autopct(values):
    """Create autopct function that shows both count and percentage."""
    def autopct(pct):
        total = sum(values)
        count = int(round(pct * total / 100.0))
        if pct > 5:  # Only show label if slice is big enough
            return f'{count}\n({pct:.0f}%)'
        elif pct > 2:
            return f'{count}'
        return ''
    return autopct


def plot_pie(ax, data, title, colors_map):
    """Plot a single pie chart."""
    # Filter out zero values
    labels = [k for k, v in data.items() if v > 0]
    values = [v for v in data.values() if v > 0]
    colors = [colors_map.get(l, '#CCCCCC') for l in labels]

    if not values:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return

    wedges, texts, autotexts = ax.pie(
        values,
        labels=labels,
        colors=colors,
        autopct=make_autopct(values),
        startangle=90,
        pctdistance=0.6,
        labeldistance=1.15,
    )

    # Style the text
    for text in texts:
        text.set_fontsize(9)
    for autotext in autotexts:
        autotext.set_fontsize(8)
        autotext.set_color('black')

    ax.set_title(title, fontsize=11, fontweight='bold')


def plot_sunburst(ax, data, title, chromosome='chr4', text_color='black', dark_mode=False):
    """
    Plot a sunburst chart with three levels:
    - Inner: Number of arrays (0, 1, 2+)
    - Middle: Array type (qA, qB, D) - only for "1 array"
    - Outer: Repeat ranges - only for qA
    """
    from matplotlib.patches import Wedge, FancyArrowPatch
    from matplotlib.collections import PatchCollection

    # Ring dimensions
    inner_radius = 0.25
    middle_radius = 0.47
    plam_radius = 0.69
    outer_radius = 0.91
    ring_width = 0.19

    # Extract data
    level1 = data['level1']  # No arrays, 1 array, 2+ arrays
    level2 = data['level2']  # C+, C-, D
    level3 = data['level3']  # pLAM status (Functional, Non-functional, Absent)
    level4 = data['level4']  # Repeat ranges (Functional pLAM only)

    total = data['total']

    # Chromosome-specific labels
    chrom_prefix = '4q' if chromosome == 'chr4' else '10q'
    qA_label = f'{chrom_prefix}A'  # C+ equivalent
    qB_label = f'{chrom_prefix}B'  # C- equivalent

    # Colors - muted pastels for dark mode, original for light mode
    if dark_mode:
        level1_colors = {
            'No arrays': '#6B7280',   # Muted gray
            '1 array': '#6EE7B7',     # Soft mint
            '2+ arrays': '#C4B5FD',   # Soft lavender
        }
        level2_colors = {
            qA_label: '#34D399',      # Emerald
            qB_label: '#60A5FA',      # Soft blue
            'D': '#F87171',           # Soft coral
        }
        plam_colors = {
            'Functional': '#40D392',
            'Non-functional': '#F07167',
            'Absent': '#6B7280',
        }
        level4_colors = {
            '1-10 (FSHD1)': '#F87171',   # Soft coral
            '11-20 (FSHD2)': '#FBBF24',  # Amber
            '21-30': '#34D399',          # Emerald
            '31+': '#60A5FA',            # Soft blue
        }
        wedge_text_color = 'black'    # Text on colored wedges
        leader_text_color = text_color  # Text on background
    else:
        level1_colors = {
            'No arrays': '#CCCCCC',
            '1 array': '#98FB98',
            '2+ arrays': '#DDA0DD',
        }
        level2_colors = {
            qA_label: '#90EE90',
            qB_label: '#87CEEB',
            'D': '#FF6B6B',
        }
        plam_colors = {
            'Functional': '#40D392',
            'Non-functional': '#F07167',
            'Absent': '#6B7280',
        }
        level4_colors = {
            '1-10 (FSHD1)': '#FF6B6B',
            '11-20 (FSHD2)': '#FFA07A',
            '21-30': '#98FB98',
            '31+': '#87CEEB',
        }
        wedge_text_color = 'black'
        leader_text_color = text_color

    # Map C+/C- to chromosome-specific labels
    level2_mapped = {
        qA_label: level2.get('C+', 0),
        qB_label: level2.get('C-', 0),
        'D': level2.get('D', 0),
    }

    # Calculate angles for level 1
    level1_order = ['No arrays', '1 array', '2+ arrays']
    level1_values = [level1.get(k, 0) for k in level1_order]
    level1_total = sum(level1_values)

    # Start from 90 degrees (top)
    start_angle = 90
    level1_angles = []

    for i, (label, value) in enumerate(zip(level1_order, level1_values)):
        if value == 0:
            level1_angles.append((start_angle, start_angle, label, value))
            continue
        angle_span = 360 * value / level1_total
        end_angle = start_angle - angle_span
        level1_angles.append((start_angle, end_angle, label, value))
        start_angle = end_angle

    # Draw level 1 (inner ring)
    for start, end, label, value in level1_angles:
        if value == 0:
            continue
        color = level1_colors.get(label, '#CCCCCC')
        wedge = Wedge(
            (0, 0), inner_radius + ring_width, end, start,
            width=ring_width, facecolor=color, edgecolor=text_color, linewidth=1.5
        )
        ax.add_patch(wedge)

        # Add label
        mid_angle = np.radians((start + end) / 2)
        label_r = inner_radius + ring_width / 2
        x = label_r * np.cos(mid_angle)
        y = label_r * np.sin(mid_angle)
        pct = value / level1_total * 100

        # "1 array" label inside, others with leader lines outside
        if label == '1 array' and pct > 3:
            display_label = f"{label}: n = {value} ({pct:.0f}%)"
            ax.text(x, y, display_label, ha='center', va='center', fontsize=7, fontweight='bold', color=wedge_text_color)
        elif value > 0:
            # Leader line to outside for "No arrays" and "2+ arrays"
            outer_r = 0.65
            outer_x = outer_r * np.cos(mid_angle)
            outer_y = outer_r * np.sin(mid_angle)
            edge_r = inner_radius + ring_width + 0.02
            edge_x = edge_r * np.cos(mid_angle)
            edge_y = edge_r * np.sin(mid_angle)
            ax.plot([edge_x, outer_x], [edge_y, outer_y], color=leader_text_color, linewidth=0.8, alpha=0.7)
            ha = 'left' if outer_x > 0 else 'right'
            display_label = f"{label}: n = {value} ({pct:.0f}%)"
            ax.text(outer_x + (0.03 if outer_x > 0 else -0.03), outer_y, display_label,
                    ha=ha, va='center', fontsize=7, fontweight='bold', color=leader_text_color)

    # Draw level 2 (middle ring) - only for "1 array" slice
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
            if value == 0:
                level2_angles.append((l2_start, l2_start, label, value))
                continue
            angle_span = l1_span * value / level2_total
            l2_end = l2_start - angle_span

            color = level2_colors.get(label, '#CCCCCC')
            wedge = Wedge(
                (0, 0), middle_radius + ring_width, l2_end, l2_start,
                width=ring_width, facecolor=color, edgecolor=text_color, linewidth=1.5
            )
            ax.add_patch(wedge)

            level2_angles.append((l2_start, l2_end, label, value))

            # Add label - inside if big enough, outside with leader line if small
            mid_angle = np.radians((l2_start + l2_end) / 2)
            label_r = middle_radius + ring_width / 2
            x = label_r * np.cos(mid_angle)
            y = label_r * np.sin(mid_angle)
            pct = value / level2_total * 100

            if pct > 8:
                # Label inside wedge
                ax.text(x, y, f"{label}: n = {value} ({pct:.0f}%)", ha='center', va='center', fontsize=7, fontweight='bold', color=wedge_text_color)
            elif value > 0:
                # Leader line to outside
                outer_r = 1.08
                outer_x = outer_r * np.cos(mid_angle)
                outer_y = outer_r * np.sin(mid_angle)
                # Draw line from wedge edge to label
                edge_r = middle_radius + ring_width + 0.02
                edge_x = edge_r * np.cos(mid_angle)
                edge_y = edge_r * np.sin(mid_angle)
                ax.plot([edge_x, outer_x], [edge_y, outer_y], color=leader_text_color, linewidth=0.8, alpha=0.7)
                # Position text
                ha = 'left' if outer_x > 0 else 'right'
                ax.text(outer_x + (0.03 if outer_x > 0 else -0.03), outer_y, f"{label}: n = {value} ({pct:.0f}%)",
                        ha=ha, va='center', fontsize=7, fontweight='bold', color=leader_text_color)

            l2_start = l2_end

        # Draw level 3 (pLAM status ring) - only for qA slice
        qA_angles = None
        for start, end, label, value in level2_angles:
            if label == qA_label and value > 0:
                qA_angles = (start, end)
                break

        if qA_angles:
            l2_start, l2_end = qA_angles
            l2_span = l2_start - l2_end

            level3_order = ['Functional', 'Non-functional', 'Absent']
            level3_values = [level3.get(k, 0) for k in level3_order]
            level3_total = sum(level3_values)

            l3_start = l2_start
            level3_angles = []

            for label, value in zip(level3_order, level3_values):
                if value == 0:
                    level3_angles.append((l3_start, l3_start, label, value))
                    continue
                angle_span = l2_span * value / level3_total
                l3_end = l3_start - angle_span

                color = plam_colors.get(label, '#CCCCCC')
                wedge = Wedge(
                    (0, 0), plam_radius + ring_width, l3_end, l3_start,
                    width=ring_width, facecolor=color, edgecolor=text_color, linewidth=1.5
                )
                ax.add_patch(wedge)

                level3_angles.append((l3_start, l3_end, label, value))

                # Add label
                mid_angle = np.radians((l3_start + l3_end) / 2)
                label_r = plam_radius + ring_width / 2
                x = label_r * np.cos(mid_angle)
                y = label_r * np.sin(mid_angle)
                pct = value / level3_total * 100

                if pct > 8:
                    ax.text(x, y, f"{label}\nn={value} ({pct:.0f}%)", ha='center', va='center', fontsize=6, fontweight='bold', color=wedge_text_color)
                elif value > 0:
                    outer_r = 1.15
                    outer_x = outer_r * np.cos(mid_angle)
                    outer_y = outer_r * np.sin(mid_angle)
                    edge_r = plam_radius + ring_width + 0.02
                    edge_x = edge_r * np.cos(mid_angle)
                    edge_y = edge_r * np.sin(mid_angle)
                    ax.plot([edge_x, outer_x], [edge_y, outer_y], color=leader_text_color, linewidth=0.8, alpha=0.7)
                    ha = 'left' if outer_x > 0 else 'right'
                    ax.text(outer_x + (0.03 if outer_x > 0 else -0.03), outer_y, f"{label}: n={value} ({pct:.0f}%)",
                            ha=ha, va='center', fontsize=6, fontweight='bold', color=leader_text_color)

                l3_start = l3_end

            # Draw level 4 (outer ring) - repeat ranges for Functional pLAM only
            func_angles = None
            for start, end, label, value in level3_angles:
                if label == 'Functional' and value > 0:
                    func_angles = (start, end)
                    break

            if func_angles:
                l3_start, l3_end = func_angles
                l3_span = l3_start - l3_end

                level4_order = ['1-10 (FSHD1)', '11-20 (FSHD2)', '21-30', '31+']
                level4_values = [level4.get(k, 0) for k in level4_order]
                level4_total = sum(level4_values)

                l4_start = l3_start

                for label, value in zip(level4_order, level4_values):
                    if value == 0:
                        continue
                    angle_span = l3_span * value / level4_total
                    l4_end = l4_start - angle_span

                    color = level4_colors.get(label, '#CCCCCC')
                    wedge = Wedge(
                        (0, 0), outer_radius + ring_width, l4_end, l4_start,
                        width=ring_width, facecolor=color, edgecolor=text_color, linewidth=1.5
                    )
                    ax.add_patch(wedge)

                    mid_angle = np.radians((l4_start + l4_end) / 2)
                    pct = value / level4_total * 100

                    outer_r = 1.22
                    outer_x = outer_r * np.cos(mid_angle)
                    outer_y = outer_r * np.sin(mid_angle)
                    edge_r = outer_radius + ring_width + 0.02
                    edge_x = edge_r * np.cos(mid_angle)
                    edge_y = edge_r * np.sin(mid_angle)
                    leader_color = color if dark_mode else text_color
                    label_color = color if dark_mode else text_color
                    ax.plot([edge_x, outer_x], [edge_y, outer_y], color=leader_color, linewidth=1.0, alpha=0.8)
                    ha_val = 'left' if outer_x > 0 else 'right'
                    short_label = label.replace(' (FSHD1)', '').replace(' (FSHD2)', '')
                    display_label = f"{short_label}: n={value} ({pct:.0f}%)"
                    ax.text(outer_x + (0.03 if outer_x > 0 else -0.03), outer_y, display_label,
                            ha=ha_val, va='center', fontsize=6, fontweight='bold', color=label_color)

                    l4_start = l4_end

    # Add center text
    ax.text(0, 0, f"n = {total}", ha='center', va='center', fontsize=10, fontweight='bold', color=text_color)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10, color=text_color)


def plot_sunburst_legend(ax, text_color='black', dark_mode=False):
    """Draw a legend for the sunburst chart."""
    from matplotlib.patches import Patch

    if dark_mode:
        legend_items = [
            # Level 1
            ('Inner Ring - Array Count:', None, 'header'),
            ('No arrays', '#6B7280', 'item'),
            ('1 array', '#6EE7B7', 'item'),
            ('2+ arrays', '#C4B5FD', 'item'),
            ('', None, 'spacer'),
            # Level 2
            ('Middle Ring - Array Type (1 array only):', None, 'header'),
            ('qA (with β-satellite)', '#34D399', 'item'),
            ('qB (no β-satellite)', '#60A5FA', 'item'),
            ('D (degraded)', '#F87171', 'item'),
            ('', None, 'spacer'),
            # Level 3 - pLAM
            ('pLAM Ring - pLAM Status (qA only):', None, 'header'),
            ('Functional pLAM', '#40D392', 'item'),
            ('Non-functional pLAM', '#F07167', 'item'),
            ('Absent', '#6B7280', 'item'),
            ('', None, 'spacer'),
            # Level 4
            ('Outer Ring - Repeat Count (Functional pLAM only):', None, 'header'),
            ('1-10 (FSHD1 range)', '#F87171', 'item'),
            ('11-20 (FSHD2 range)', '#FBBF24', 'item'),
            ('21-30 (normal)', '#34D399', 'item'),
            ('31+ (long)', '#60A5FA', 'item'),
        ]
    else:
        legend_items = [
            # Level 1
            ('Inner Ring - Array Count:', None, 'header'),
            ('No arrays', '#CCCCCC', 'item'),
            ('1 array', '#98FB98', 'item'),
            ('2+ arrays', '#DDA0DD', 'item'),
            ('', None, 'spacer'),
            # Level 2
            ('Middle Ring - Array Type (1 array only):', None, 'header'),
            ('qA (with β-satellite)', '#90EE90', 'item'),
            ('qB (no β-satellite)', '#87CEEB', 'item'),
            ('D (degraded)', '#FF6B6B', 'item'),
            ('', None, 'spacer'),
            # Level 3 - pLAM
            ('pLAM Ring - pLAM Status (qA only):', None, 'header'),
            ('Functional pLAM', '#40D392', 'item'),
            ('Non-functional pLAM', '#F07167', 'item'),
            ('Absent', '#6B7280', 'item'),
            ('', None, 'spacer'),
            # Level 4
            ('Outer Ring - Repeat Count (Functional pLAM only):', None, 'header'),
            ('1-10 (FSHD1 range)', '#FF6B6B', 'item'),
            ('11-20 (FSHD2 range)', '#FFA07A', 'item'),
            ('21-30 (normal)', '#98FB98', 'item'),
            ('31+ (long)', '#87CEEB', 'item'),
        ]

    y = 0.95
    for label, color, item_type in legend_items:
        if item_type == 'header':
            ax.text(0.05, y, label, fontsize=9, fontweight='bold', transform=ax.transAxes, color=text_color)
            y -= 0.06
        elif item_type == 'item':
            patch = Patch(facecolor=color, edgecolor='gray', linewidth=0.5)
            ax.add_patch(plt.Rectangle((0.05, y - 0.02), 0.04, 0.035,
                                        facecolor=color, edgecolor='gray',
                                        linewidth=0.5, transform=ax.transAxes))
            ax.text(0.12, y, label, fontsize=8, transform=ax.transAxes, va='center', color=text_color)
            y -= 0.05
        elif item_type == 'spacer':
            y -= 0.03

    ax.axis('off')


def plot_distribution(df, output_file, dark_mode=False):
    """Create the full distribution figure with sunburst charts."""

    # Load pLAM data
    plam_file = os.path.join(os.path.dirname(output_file) or '.', 'd4z4_array_coordinates.plam.tsv')
    plam_df = pd.read_csv(plam_file, sep='\t') if os.path.exists(plam_file) else None

    # Calculate distributions for both chromosomes
    chr4_data = calculate_distributions(df, 'chr4', plam_df=plam_df)
    chr10_data = calculate_distributions(df, 'chr10', plam_df=plam_df)

    # Set up colors based on mode
    if dark_mode:
        bg_color = '#1a1a1a'
        text_color = '#FFFFFF'
        plt.rcParams['text.color'] = text_color
        plt.rcParams['axes.labelcolor'] = text_color
        plt.rcParams['xtick.color'] = text_color
        plt.rcParams['ytick.color'] = text_color
    else:
        bg_color = 'white'
        text_color = 'black'

    # Create figure: 1 row with 2 sunbursts + legend
    fig = plt.figure(figsize=(16, 7), facecolor=bg_color)

    # Create grid: 2 sunbursts on left/center, legend on right
    ax1 = fig.add_axes([0.02, 0.1, 0.4, 0.8])  # chr4
    ax2 = fig.add_axes([0.42, 0.1, 0.4, 0.8])  # chr10
    ax_legend = fig.add_axes([0.82, 0.1, 0.17, 0.8])  # legend

    for ax in [ax1, ax2, ax_legend]:
        ax.set_facecolor(bg_color)

    # Plot sunbursts
    plot_sunburst(ax1, chr4_data, f"chr4 (n={chr4_data['total']})", chromosome='chr4', text_color=text_color, dark_mode=dark_mode)
    plot_sunburst(ax2, chr10_data, f"chr10 (n={chr10_data['total']})", chromosome='chr10', text_color=text_color, dark_mode=dark_mode)

    # Plot legend
    plot_sunburst_legend(ax_legend, text_color=text_color, dark_mode=dark_mode)

    fig.suptitle('D4Z4 Array Distribution in HPRC Samples', fontsize=14, fontweight='bold', y=0.98, color=text_color)

    # Save
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor=bg_color)
    print(f"Saved: {output_file}")

    # Also save SVG and PDF
    if output_file.endswith('.png'):
        svg_file = output_file.replace('.png', '.svg')
        plt.savefig(svg_file, bbox_inches='tight', facecolor=bg_color)
        print(f"Saved: {svg_file}")
        pdf_file = output_file.replace('.png', '.pdf')
        plt.savefig(pdf_file, bbox_inches='tight', facecolor=bg_color)
        print(f"Saved: {pdf_file}")

    plt.close()

    # Print summary
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)

    for chrom, data in [('chr4', chr4_data), ('chr10', chr10_data)]:
        print(f"\n{chrom} (n={data['total']}):")
        print(f"  Level 1 - Array count:")
        for k, v in data['level1'].items():
            pct = v / data['total'] * 100
            print(f"    {k}: {v} ({pct:.1f}%)")

        single_total = sum(data['level2'].values())
        if single_total > 0:
            print(f"  Level 2 - Single array type (n={single_total}):")
            for k, v in data['level2'].items():
                pct = v / single_total * 100 if single_total > 0 else 0
                print(f"    {k}: {v} ({pct:.1f}%)")

        plam_total = sum(data['level3'].values())
        if plam_total > 0:
            print(f"  Level 3 - pLAM status (n={plam_total}):")
            for k, v in data['level3'].items():
                pct = v / plam_total * 100 if plam_total > 0 else 0
                print(f"    {k}: {v} ({pct:.1f}%)")

        func_total = sum(data['level4'].values())
        if func_total > 0:
            print(f"  Level 4 - Functional pLAM repeat ranges (n={func_total}):")
            for k, v in data['level4'].items():
                pct = v / func_total * 100 if func_total > 0 else 0
                print(f"    {k}: {v} ({pct:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description="Plot D4Z4 array distribution")
    parser.add_argument("--analysis-file",
                        default="results/d4z4_array_analysis.tsv",
                        help="Analysis TSV file")
    parser.add_argument("--output",
                        default="results/array_distribution.png",
                        help="Output file (PNG or PDF)")
    parser.add_argument("--background", choices=["white", "black"], default="white",
                        help="Background color (default: white)")
    args = parser.parse_args()

    print("Loading data...")
    df = load_data(args.analysis_file)

    dark_mode = args.background == "black"
    print(f"Generating distribution plots (dark_mode={dark_mode})...")
    plot_distribution(df, args.output, dark_mode=dark_mode)


if __name__ == '__main__':
    main()
