#!/usr/bin/env python3
"""
Create D4Z4 Classification Pipeline Diagram

Generates an SVG diagram showing the classification pipeline steps.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

# Configure matplotlib to embed text as text (not paths) in SVG
plt.rcParams['svg.fonttype'] = 'none'

def create_pipeline_diagram(output_path='results/pipeline_diagram.svg', dark_mode=True):
    """Create the pipeline diagram."""

    fig, ax = plt.subplots(1, 1, figsize=(16, 14))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 14)
    ax.axis('off')

    # Set background color
    if dark_mode:
        fig.patch.set_facecolor('#1a1a1a')
        ax.set_facecolor('#1a1a1a')
        text_color = '#FFFFFF'
        text_secondary = '#CCCCCC'
        text_muted = '#999999'
        box_edge = '#666666'
        arrow_color = '#AAAAAA'
    else:
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        text_color = '#333333'
        text_secondary = '#444444'
        text_muted = '#666666'
        box_edge = '#666666'
        arrow_color = '#455A64'

    # Colors - matching KaryoScope/plot_d4z4_arrays.py color scheme
    colors = {
        'input': '#E3F2FD',      # Light blue
        'process': '#3d3225' if dark_mode else '#FFF3E0',    # Dark orange / Light orange
        'decision': '#F3E5F5',   # Light purple
        'output': '#E8F5E9',     # Light green
        'canonical': '#90EE90',  # Light green - canonical array (C+)
        'canonical_ct': '#87CEEB',  # Light blue - canonical without β-sat (C-)
        'degraded': '#FF6B6B',   # Light red - degraded array
        'bsat': '#fa99ff',       # Light pink/magenta - bsat_specific
        'exon': '#4444FF' if dark_mode else '#090088',       # Brighter blue for dark mode
        'intron': '#22BB66' if dark_mode else '#198450',     # Brighter green for dark mode
        'intergenic': '#E89EB8', # Light pink - intergenic_specific
        'plam': '#33CC55' if dark_mode else '#228B22',       # Brighter green - β-satellite terminal
        'ct': '#C8BFB2' if dark_mode else '#B8AFA2',         # Lighter brown - ct_specific
        'line': '#FFB366',       # Orange - LINE repeats
        'sine': '#66B3FF',       # Light blue - SINE repeats
        'ltr': '#FF6666',        # Light red - LTR repeats
        'arrow': arrow_color,
        'text': text_color,
        'text_secondary': text_secondary,
        'text_muted': text_muted,
        'box_edge': box_edge,
    }

    # Title
    ax.text(8, 13.7, 'D4Z4 Array Classification Pipeline (v4)',
            fontsize=18, fontweight='bold', ha='center', va='top', color=colors['text'])

    # =========================================================================
    # LEFT SIDE: Pipeline Steps
    # =========================================================================

    box_width = 4.0
    box_height = 0.8
    left_x = 0.3

    steps = [
        ('Step 1: Extract BSR', 'bsat_specific fragments\nfrom KaryoScope', 12.5),
        ('Step 2: Merge Clusters', 'Merge fragments\n(gap ≤ 400bp)', 11.0),
        ('Step 3: Annotate Clusters', 'Calculate exon fraction\nfrom gene track', 9.5),
        ('Step 4: Identify Arrays', 'Group exon-dense clusters\n(≥50% exon, ≥1500bp)', 8.0),
        ('Step 5: Extend Arrays', 'Add trailing fragments\n(gap ≤ 1500bp)', 6.5),
        ('Step 6: Find Terminals', 'Search for terminal β-satellite\n(pLAM proxy*) or ct_specific', 5.0),
        ('Step 7: Classify', 'Apply classification rules', 3.5),
    ]

    for i, (title, desc, y) in enumerate(steps):
        # Box
        box = FancyBboxPatch((left_x, y - box_height/2), box_width, box_height,
                             boxstyle="round,pad=0.03,rounding_size=0.1",
                             facecolor=colors['process'], edgecolor=colors['box_edge'], linewidth=1.5)
        ax.add_patch(box)

        # Text
        ax.text(left_x + box_width/2, y + 0.18, title,
                fontsize=10, fontweight='bold', ha='center', va='center', color=colors['text'])
        ax.text(left_x + box_width/2, y - 0.22, desc,
                fontsize=8, ha='center', va='center', color=colors['text_secondary'])

        # Arrow to next step
        if i < len(steps) - 1:
            ax.annotate('', xy=(left_x + box_width/2, steps[i+1][2] + box_height/2 + 0.05),
                       xytext=(left_x + box_width/2, y - box_height/2 - 0.05),
                       arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))

    # =========================================================================
    # RIGHT SIDE: Visual Examples with Multiple Tracks
    # =========================================================================

    right_x = 5.5
    track_width = 10

    # Example title
    ax.text(right_x + track_width/2, 13.3, 'Visual Examples',
            fontsize=14, fontweight='bold', ha='center', color=colors['text'])

    # ----- EXAMPLE 1: Canonical C+ Array -----
    ax.text(right_x, 12.8, 'Example 1: Canonical (C+) - β-satellite terminal',
            fontsize=11, fontweight='bold', color=colors['canonical'])

    track_height = 0.25
    track_gap = 0.4

    # --- Satellites track (bsat_specific) ---
    y = 12.3
    ax.text(right_x - 0.1, y, 'Satellites:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    # Draw bsat fragments across array region
    bsat_positions = [0.1, 0.4, 0.7, 1.1, 1.4, 1.7, 2.1, 2.4, 2.7, 3.1, 3.4, 3.7,
                      4.1, 4.4, 4.7, 5.1, 5.4, 6.0, 6.3, 6.6]
    for x in bsat_positions:
        rect = Rectangle((right_x + x, y - track_height/2), 0.2, track_height,
                         facecolor=colors['bsat'], edgecolor='none')
        ax.add_patch(rect)

    # --- Genes track (exon/intron) ---
    y = 12.3 - track_gap
    ax.text(right_x - 0.1, y, 'Genes:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    # Array region has exons (DUX4)
    gene_regions = [
        (0.1, 0.5, 'exon'), (0.5, 0.7, 'intron'), (0.7, 1.0, 'exon'),
        (1.1, 1.5, 'exon'), (1.5, 1.7, 'intron'), (1.7, 2.0, 'exon'),
        (2.1, 2.5, 'exon'), (2.5, 2.7, 'intron'), (2.7, 3.0, 'exon'),
        (3.1, 3.5, 'exon'), (3.5, 3.7, 'intron'), (3.7, 4.0, 'exon'),
        (4.1, 4.5, 'exon'), (4.5, 4.7, 'intron'), (4.7, 5.0, 'exon'),
        # pLAM region - intron-rich
        (5.8, 6.0, 'intron'), (6.0, 6.2, 'intergenic'), (6.2, 6.8, 'intron'),
    ]
    for (xs, xe, gtype) in gene_regions:
        color = colors['exon'] if gtype == 'exon' else (colors['intron'] if gtype == 'intron' else colors['intergenic'])
        rect = Rectangle((right_x + xs, y - track_height/2), xe - xs, track_height,
                         facecolor=color, edgecolor='none', alpha=0.8)
        ax.add_patch(rect)

    # --- Repeats track ---
    y = 12.3 - 2*track_gap
    ax.text(right_x - 0.1, y, 'Repeats:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    # Some LINE/SINE scattered
    repeat_regions = [
        (5.2, 5.5, 'line'), (7.0, 7.3, 'sine'), (7.5, 7.8, 'line'),
    ]
    for (xs, xe, rtype) in repeat_regions:
        color = colors['line'] if rtype == 'line' else colors['sine']
        rect = Rectangle((right_x + xs, y - track_height/2), xe - xs, track_height,
                         facecolor=color, edgecolor='none', alpha=0.7)
        ax.add_patch(rect)

    # --- Array detection result ---
    y = 12.3 - 3*track_gap
    ax.text(right_x - 0.1, y, 'Result:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])

    # Array box (C+)
    array_box = FancyBboxPatch((right_x + 0.05, y - 0.2), 5.1, 0.4,
                               boxstyle="round,pad=0.02,rounding_size=0.05",
                               facecolor=colors['canonical'], edgecolor=colors['box_edge'],
                               linewidth=2, alpha=0.5)
    ax.add_patch(array_box)

    # Array elements
    for i, xs in enumerate([0.1, 1.1, 2.1, 3.1, 4.1]):
        rect = Rectangle((right_x + xs, y - 0.12), 0.8, 0.24,
                         facecolor=colors['bsat'], edgecolor=colors['box_edge'], linewidth=0.5)
        ax.add_patch(rect)
        ax.text(right_x + xs + 0.4, y, str(i+1), fontsize=7, ha='center', va='center',
                color='#333', fontweight='bold')

    # pLAM terminal
    rect = Rectangle((right_x + 5.8, y - 0.12), 1.0, 0.24,
                     facecolor=colors['plam'], edgecolor=colors['box_edge'], linewidth=0.5)
    ax.add_patch(rect)
    ax.text(right_x + 6.3, y, 'β-sat', fontsize=7, ha='center', va='center', color='white', fontweight='bold')

    # Bracket showing pLAM criteria
    ax.annotate('intron-rich (<35% exon)\n*pLAM proxy only', xy=(right_x + 6.3, y + 0.35),
               fontsize=7, ha='center', va='bottom', color=colors['text_muted'])

    # ----- EXAMPLE 2: Canonical C- Array (no β-satellite terminal) -----
    y_base = 9.8
    ax.text(right_x, y_base + 0.5, 'Example 2: Canonical (C-) - no β-satellite terminal',
            fontsize=11, fontweight='bold', color=colors['canonical_ct'])

    # --- Satellites track ---
    y = y_base
    ax.text(right_x - 0.1, y, 'Satellites:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    # bsat in array only, then ct_specific region
    bsat_positions2 = [0.1, 0.4, 0.7, 1.1, 1.4, 1.7, 2.1, 2.4, 2.7, 3.1, 3.4, 3.7, 4.1, 4.4]
    for x in bsat_positions2:
        rect = Rectangle((right_x + x, y - track_height/2), 0.2, track_height,
                         facecolor=colors['bsat'], edgecolor='none')
        ax.add_patch(rect)
    # ct_specific region after array
    rect = Rectangle((right_x + 5.5, y - track_height/2), 2.5, track_height,
                     facecolor=colors['ct'], edgecolor='none', alpha=0.7)
    ax.add_patch(rect)
    ax.text(right_x + 6.75, y, 'ct_specific', fontsize=6, ha='center', va='center', color='#333')

    # --- Genes track ---
    y = y_base - track_gap
    ax.text(right_x - 0.1, y, 'Genes:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    # Array has exons, terminal region is mixed
    gene_regions2 = [
        (0.1, 0.5, 'exon'), (0.5, 0.7, 'intron'), (0.7, 1.0, 'exon'),
        (1.1, 1.5, 'exon'), (1.5, 1.7, 'intron'), (1.7, 2.0, 'exon'),
        (2.1, 2.5, 'exon'), (2.5, 2.7, 'intron'), (2.7, 3.0, 'exon'),
        (3.1, 3.5, 'exon'), (3.5, 3.7, 'intron'), (3.7, 4.0, 'exon'),
        (4.1, 4.5, 'exon'), (4.5, 4.7, 'intron'),
        # Terminal region - mixed (no pLAM signature)
        (5.5, 6.0, 'exon'), (6.0, 6.5, 'intergenic'), (6.5, 7.0, 'exon'), (7.0, 7.5, 'intergenic'),
    ]
    for (xs, xe, gtype) in gene_regions2:
        color = colors['exon'] if gtype == 'exon' else (colors['intron'] if gtype == 'intron' else colors['intergenic'])
        rect = Rectangle((right_x + xs, y - track_height/2), xe - xs, track_height,
                         facecolor=color, edgecolor='none', alpha=0.8)
        ax.add_patch(rect)

    # --- Repeats track ---
    y = y_base - 2*track_gap
    ax.text(right_x - 0.1, y, 'Repeats:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    repeat_regions2 = [(4.8, 5.2, 'line'), (8.0, 8.5, 'sine')]
    for (xs, xe, rtype) in repeat_regions2:
        color = colors['line'] if rtype == 'line' else colors['sine']
        rect = Rectangle((right_x + xs, y - track_height/2), xe - xs, track_height,
                         facecolor=color, edgecolor='none', alpha=0.7)
        ax.add_patch(rect)

    # --- Result track ---
    y = y_base - 3*track_gap
    ax.text(right_x - 0.1, y, 'Result:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])

    # Array box (C-)
    array_box = FancyBboxPatch((right_x + 0.05, y - 0.2), 4.6, 0.4,
                               boxstyle="round,pad=0.02,rounding_size=0.05",
                               facecolor=colors['canonical_ct'], edgecolor=colors['box_edge'],
                               linewidth=2, alpha=0.5)
    ax.add_patch(array_box)

    # Array elements
    for i, xs in enumerate([0.1, 1.1, 2.1, 3.1]):
        rect = Rectangle((right_x + xs, y - 0.12), 0.8, 0.24,
                         facecolor=colors['bsat'], edgecolor=colors['box_edge'], linewidth=0.5)
        ax.add_patch(rect)
        ax.text(right_x + xs + 0.4, y, str(i+1), fontsize=7, ha='center', va='center',
                color='#333', fontweight='bold')

    # ct_specific terminal marker
    rect = Rectangle((right_x + 5.5, y - 0.12), 2.5, 0.24,
                     facecolor=colors['ct'], edgecolor=colors['box_edge'], linewidth=0.5, alpha=0.7)
    ax.add_patch(rect)
    ax.text(right_x + 6.75, y, 'ct_specific', fontsize=7, ha='center', va='center', color='#333')

    # Annotation: no β-sat terminal
    ax.annotate('No β-sat terminal\n(no intron-rich bsat)', xy=(right_x + 4.9, y + 0.35),
               fontsize=7, ha='center', va='bottom', color=colors['text_muted'])

    # ----- EXAMPLE 3: Degraded Array -----
    y_base = 7.0
    ax.text(right_x, y_base + 0.5, 'Example 3: Degraded (D) - no terminal detected',
            fontsize=11, fontweight='bold', color=colors['degraded'])

    # --- Satellites track ---
    y = y_base
    ax.text(right_x - 0.1, y, 'Satellites:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    bsat_positions3 = [0.1, 0.4, 0.7, 1.1, 1.4, 1.7, 2.1, 2.4, 2.7]
    for x in bsat_positions3:
        rect = Rectangle((right_x + x, y - track_height/2), 0.2, track_height,
                         facecolor=colors['bsat'], edgecolor='none')
        ax.add_patch(rect)
    # Gap then scattered fragments (no terminal structure)

    # --- Genes track ---
    y = y_base - track_gap
    ax.text(right_x - 0.1, y, 'Genes:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])
    gene_regions3 = [
        (0.1, 0.5, 'exon'), (0.5, 0.7, 'intron'), (0.7, 1.0, 'exon'),
        (1.1, 1.5, 'exon'), (1.5, 1.7, 'intron'), (1.7, 2.0, 'exon'),
        (2.1, 2.5, 'exon'), (2.5, 2.7, 'intron'), (2.7, 3.0, 'exon'),
        # After array - intergenic, no structure
        (3.5, 5.0, 'intergenic'),
    ]
    for (xs, xe, gtype) in gene_regions3:
        color = colors['exon'] if gtype == 'exon' else (colors['intron'] if gtype == 'intron' else colors['intergenic'])
        rect = Rectangle((right_x + xs, y - track_height/2), xe - xs, track_height,
                         facecolor=color, edgecolor='none', alpha=0.8)
        ax.add_patch(rect)

    # --- Result track ---
    y = y_base - 2*track_gap
    ax.text(right_x - 0.1, y, 'Result:', fontsize=8, ha='right', va='center', color=colors['text_secondary'])

    # Array box (Degraded)
    array_box = FancyBboxPatch((right_x + 0.05, y - 0.2), 3.0, 0.4,
                               boxstyle="round,pad=0.02,rounding_size=0.05",
                               facecolor=colors['degraded'], edgecolor=colors['box_edge'],
                               linewidth=2, alpha=0.5)
    ax.add_patch(array_box)

    # Array elements
    for i, xs in enumerate([0.1, 1.1, 2.1]):
        rect = Rectangle((right_x + xs, y - 0.12), 0.8, 0.24,
                         facecolor=colors['bsat'], edgecolor=colors['box_edge'], linewidth=0.5)
        ax.add_patch(rect)
        ax.text(right_x + xs + 0.4, y, str(i+1), fontsize=7, ha='center', va='center',
                color='#333', fontweight='bold')

    # Annotation
    ax.annotate('No β-sat terminal, no ct_specific', xy=(right_x + 3.5, y + 0.35),
               fontsize=7, ha='center', va='bottom', color=colors['text_muted'])

    # =========================================================================
    # BOTTOM: Classification Rules (moved down)
    # =========================================================================

    y = 1.8
    ax.text(8, y + 0.6, 'Classification Rules', fontsize=13, fontweight='bold', ha='center', color=colors['text'])

    # Rule boxes
    rules = [
        ('Canonical\n(C+)', 'β-sat terminal* + ≥1 elem', colors['canonical'], 2.5),
        ('Canonical\n(C-)', 'ct_specific + ≥3 elem\n(no β-sat)', colors['canonical_ct'], 6.5),
        ('Degraded', '≥3 elements\nno terminal', colors['degraded'], 10.5),
        ('Transition\nElement', '<3 elements', '#FFD700', 14),  # Gold - matches plot_d4z4_arrays.py
    ]

    for label, desc, color, x in rules:
        box = FancyBboxPatch((x - 1.4, y - 0.5), 2.8, 0.9,
                             boxstyle="round,pad=0.02,rounding_size=0.1",
                             facecolor=color, edgecolor=colors['box_edge'], linewidth=1, alpha=0.85)
        ax.add_patch(box)
        # Use dark text for light backgrounds (green, blue, gold)
        rule_text_color = '#333' if color in [colors['canonical'], colors['canonical_ct'], '#FFD700'] else 'white'
        ax.text(x, y + 0.15, label, fontsize=11, fontweight='bold', ha='center', va='center', color=rule_text_color)
        ax.text(x, y - 0.2, desc, fontsize=7, ha='center', va='center', color=rule_text_color)

    # =========================================================================
    # LEGEND
    # =========================================================================

    y = 0.5
    ax.text(0.3, y, 'Legend:', fontsize=10, fontweight='bold', color=colors['text'])

    legend_items = [
        (colors['bsat'], 'bsat_specific'),
        (colors['exon'], 'exon'),
        (colors['intron'], 'intron'),
        (colors['intergenic'], 'intergenic'),
        (colors['ct'], 'ct_specific'),
        (colors['plam'], 'β-sat terminal*'),
        (colors['line'], 'LINE'),
        (colors['sine'], 'SINE'),
    ]

    x = 1.8
    for color, label in legend_items:
        rect = Rectangle((x, y - 0.12), 0.4, 0.24, facecolor=color, edgecolor=colors['box_edge'], linewidth=0.5)
        ax.add_patch(rect)
        ax.text(x + 0.55, y, label, fontsize=8, va='center', color=colors['text_secondary'])
        x += 1.8

    # Footnote explaining pLAM proxy limitation
    ax.text(8, 0.1, '* β-satellite terminal is a structural proxy for pLAM; actual sequence (ATTAAA vs ATCAAA) cannot currently be confirmed from k-mer data',
            fontsize=7, ha='center', color=colors['text_muted'], style='italic')

    # Save
    plt.tight_layout()
    plt.savefig(output_path, format='svg', bbox_inches='tight', dpi=150)
    plt.savefig(output_path.replace('.svg', '.png'), format='png', bbox_inches='tight', dpi=150)
    plt.close()

    print(f"Saved: {output_path}")
    print(f"Saved: {output_path.replace('.svg', '.png')}")


def main():
    parser = argparse.ArgumentParser(description="Create D4Z4 Classification Pipeline Diagram")
    parser.add_argument("--output", default="pipeline_diagram.svg",
                        help="Output SVG file path (default: pipeline_diagram.svg)")
    parser.add_argument("--background", choices=["white", "black"], default="white",
                        help="Background color (default: white)")
    args = parser.parse_args()

    dark_mode = args.background == "black"
    create_pipeline_diagram(output_path=args.output, dark_mode=dark_mode)


if __name__ == '__main__':
    main()
