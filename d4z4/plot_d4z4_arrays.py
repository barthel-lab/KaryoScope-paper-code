#!/usr/bin/env python3
"""
D4Z4 Array Visualization Script

Generates SVG plots showing D4Z4 array structure with multiple annotation layers:
- Region annotations (bsat_specific, ct_specific, etc.)
- Gene annotations (exon_specific, intergenic_specific, etc.)
- Repeat annotations (LINE_specific, SINE_specific, LTR_specific)
- Subtelomeric annotations (TAR1_specific)
- Inferred main array boundary and repeat counts

Usage:
    python plot_d4z4_arrays.py --haplotypes chr4_hap2_HG02129#2#CM085859.1 --output d4z4_plot.svg
    python plot_d4z4_arrays.py --samples HG02129 HG00099 --output d4z4_plot.svg
    python plot_d4z4_arrays.py --chromosome chr4 --max-samples 10 --output d4z4_plot.svg
"""

import argparse
import json
import math
import os
from collections import defaultdict
from datetime import datetime

import drawsvg as draw
import pandas as pd

# =============================================================================
# Constants
# =============================================================================

EXPECTED_BSR_SIZE = 2300  # Expected BSR cluster size in bp
EXPECTED_REPEAT_UNIT = 3300  # Full D4Z4 repeat unit size
BSR_MERGE_THRESHOLD = 800  # Merge bSat blocks within this distance (increased to capture within-BSR gaps)
MIN_VALID_CLUSTER_SIZE = 1500  # Clusters smaller than this may be split fragments

# pLAM detection thresholds - pLAM must be intron/intergenic (NOT exon-specific)
PLAM_SIZE_THRESHOLD = 3000  # Minimum size for pLAM (with intron/intergenic annotation)

# Color scheme for cluster types
COLOR_CANONICAL_BSAT = '#90EE90'   # Light green - canonical array with β-satellite terminal
COLOR_CANONICAL_CT = '#87CEEB'     # Light blue - canonical array lacking β-satellite terminal
COLOR_DEGRADED_ARRAY = '#FF6B6B'   # Light red - degraded array elements
COLOR_TRANSITION = '#FFD700'       # Gold - transition regions
COLOR_PLAM = '#228B22'             # Dark green - β-satellite terminal (array-associated)
COLOR_SOLITARY_PLAM = '#8B0000'    # Dark red - solitary β-satellite (not array-associated)
COLOR_PLAM_FUNCTIONAL = '#40D392'  # Green - functional pLAM sequence
COLOR_PLAM_NONFUNCTIONAL = '#F07167'  # Coral - non-functional pLAM sequence

# Default color for unknown features
DEFAULT_COLOR = '#CCCCCC'


def _style(style, key, default):
    """Look up an optional drawing-style override (see make_figures.py); default otherwise."""
    if style and style.get(key) is not None:
        return style[key]
    return default


# =============================================================================
# Color Loading Functions
# =============================================================================

def load_color_file(filepath):
    """Load color mappings from a KaryoScope color file"""
    colors = {}
    if not os.path.exists(filepath):
        return colors

    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                feature = parts[0]
                color = parts[1]
                # Skip header
                if i == 0 and feature.lower() == 'feature':
                    continue
                colors[feature] = color
    return colors


def load_all_colors(colors_dir):
    """Load all color files from KaryoScope resources"""
    color_files = {
        'region': 'KS_human_CHM13.region.colors.txt',
        'gene': 'KS_human_CHM13.gene.colors.txt',
        'repeat': 'KS_human_CHM13.repeat.colors.txt',
        'subtelomeric': 'KS_human_CHM13.subtelomeric.colors.txt',
    }

    all_colors = {}
    for annot_type, filename in color_files.items():
        filepath = os.path.join(colors_dir, filename)
        all_colors[annot_type] = load_color_file(filepath)

    return all_colors


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_bed_file(filepath):
    """Load BED file into list of tuples (chrom, start, end, annotation)"""
    regions = []
    if not os.path.exists(filepath):
        return regions
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                annotation = parts[3]
                regions.append((chrom, start, end, annotation))
    return regions


def load_analysis_results(results_file):
    """Load D4Z4 analysis results TSV"""
    if os.path.exists(results_file):
        return pd.read_csv(results_file, sep='\t')
    return None


def get_haplotype_data(hap_name, data_dir, chromosome):
    """Load all annotation data for a single haplotype"""
    chrom = chromosome if chromosome else 'chr4'

    # Load each annotation type
    region_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.region.bed')
    gene_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.gene.bed')
    repeat_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.repeat.bed')
    subtelo_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.subtelomeric.bed')

    data = {
        'region': [],
        'gene': [],
        'repeat': [],
        'subtelomeric': []
    }

    for annot_type, filepath in [('region', region_file), ('gene', gene_file),
                                   ('repeat', repeat_file), ('subtelomeric', subtelo_file)]:
        all_regions = load_bed_file(filepath)
        data[annot_type] = [(s, e, a) for (c, s, e, a) in all_regions if c == hap_name]

    return data


# =============================================================================
# Smoothing Functions (simplified from KaryoScope smooth_features.py)
# =============================================================================

def smooth_regions(regions, max_gap=500):
    """
    Apply simplified smoothing to annotation regions.

    This is a simplified version of KaryoScope's phylogeny-based smoothing that:
    1. Merges adjacent intervals with identical features
    2. Fills small gaps when flanking features match
    3. Removes tiny fragments (< 50bp) surrounded by different features

    Args:
        regions: List of (start, end, annotation) tuples, sorted by start
        max_gap: Maximum gap size to fill (default 500bp)

    Returns:
        Smoothed list of (start, end, annotation) tuples
    """
    if not regions or len(regions) < 2:
        return regions

    # Sort by start position
    sorted_regions = sorted(regions, key=lambda x: x[0])

    # Pass 1: Merge adjacent identical features
    merged = []
    curr_start, curr_end, curr_annot = sorted_regions[0]

    for start, end, annot in sorted_regions[1:]:
        if annot == curr_annot and start <= curr_end + max_gap:
            # Same annotation and close enough - extend
            curr_end = max(curr_end, end)
        else:
            merged.append((curr_start, curr_end, curr_annot))
            curr_start, curr_end, curr_annot = start, end, annot
    merged.append((curr_start, curr_end, curr_annot))

    # Pass 2: Fill small gaps when flanking features match
    if len(merged) < 3:
        return merged

    filled = [merged[0]]
    for i in range(1, len(merged) - 1):
        prev_start, prev_end, prev_annot = filled[-1]
        curr_start, curr_end, curr_annot = merged[i]
        next_start, next_end, next_annot = merged[i + 1]

        # Check if current region is a small fragment between matching features
        curr_size = curr_end - curr_start
        gap_before = curr_start - prev_end
        gap_after = next_start - curr_end

        if (curr_size < 200 and
            prev_annot == next_annot and
            gap_before < max_gap and
            gap_after < max_gap):
            # Small fragment between matching features - absorb into previous
            filled[-1] = (prev_start, curr_end, prev_annot)
        else:
            filled.append((curr_start, curr_end, curr_annot))

    # Add the last region
    filled.append(merged[-1])

    # Pass 3: Final merge of adjacent identical features
    final = []
    curr_start, curr_end, curr_annot = filled[0]

    for start, end, annot in filled[1:]:
        if annot == curr_annot and start <= curr_end + 100:  # Small gap tolerance
            curr_end = max(curr_end, end)
        else:
            final.append((curr_start, curr_end, curr_annot))
            curr_start, curr_end, curr_annot = start, end, annot
    final.append((curr_start, curr_end, curr_annot))

    return final


def smooth_all_tracks(hap_data, max_gap=500):
    """
    Apply smoothing to all annotation tracks.

    Args:
        hap_data: Dict with keys 'region', 'gene', 'repeat', 'subtelomeric'
                  Each value is a list of (start, end, annotation) tuples
        max_gap: Maximum gap to fill during smoothing

    Returns:
        Dict with same structure but smoothed regions
    """
    smoothed = {}
    for track_name, regions in hap_data.items():
        if regions:
            smoothed[track_name] = smooth_regions(regions, max_gap)
        else:
            smoothed[track_name] = regions
    return smoothed


# =============================================================================
# BSR Cluster Detection - Refactored Pipeline (v2)
# =============================================================================
#
# Core principle: Arrays are defined by exon-dense clusters. Everything else is transition.
#
# Pipeline:
#   1. get_bsr_units() - Conservative merge to reconstruct BSR units
#   2. identify_arrays() - Find contiguous exon-dense stretches
#   3. classify_arrays() - Check for terminals (pLAM or ct_specific)
#   4. Drawing uses simple logic: in array = element, after array = pLAM, else = transition
#

# Thresholds for array detection
EXON_DENSE_THRESHOLD = 0.5       # Minimum exon fraction to be considered exon-dense
MIN_ARRAY_ELEMENT_SIZE = 1500    # Minimum size for an array element
MAX_INTER_ELEMENT_GAP = 3500     # Maximum gap between elements in same array (~1 repeat unit)
MIN_ARRAY_ELEMENTS = 3           # Minimum elements for a valid array
PLAM_MIN_SIZE = 2500             # Minimum size for pLAM terminal
PLAM_MAX_EXON_FRACTION = 0.25    # pLAM must be < 25% exon
CT_TERMINAL_SEARCH_DIST = 5000   # Search distance for ct_specific terminal


def get_bsr_units(bsat_regions, merge_threshold=500):
    """
    Merge bsat fragments into BSR units using CONSERVATIVE threshold.

    Only merges fragments that are clearly part of the same BSR unit.
    Expected BSR size: ~2300bp. Within-BSR gaps from binning: typically 100-400bp.

    Args:
        bsat_regions: List of (start, end) tuples for bsat_specific regions
        merge_threshold: Maximum gap to merge (default 500bp - conservative)

    Returns:
        List of (cluster_start, cluster_end) tuples
    """
    if not bsat_regions:
        return []

    sorted_regions = sorted(bsat_regions, key=lambda x: x[0])

    clusters = []
    cluster_start, cluster_end = sorted_regions[0]

    for start, end in sorted_regions[1:]:
        gap = start - cluster_end
        if gap <= merge_threshold:
            # Extend current cluster
            cluster_end = end
        else:
            # Finalize current cluster, start new one
            clusters.append((cluster_start, cluster_end))
            cluster_start, cluster_end = start, end

    # Don't forget the last cluster
    clusters.append((cluster_start, cluster_end))

    return clusters


def get_cluster_annotation_fractions(cluster_start, cluster_end, gene_regions):
    """
    Calculate fraction of cluster covered by each gene annotation type.

    Returns:
        dict with keys 'exon', 'intron', 'intergenic' and fraction values
    """
    size = cluster_end - cluster_start
    if size <= 0 or not gene_regions:
        return {'exon': 0, 'intron': 0, 'intergenic': 0}

    fractions = {'exon': 0, 'intron': 0, 'intergenic': 0}

    for gs, ge, ga in gene_regions:
        if ge > cluster_start and gs < cluster_end:
            overlap = min(ge, cluster_end) - max(gs, cluster_start)
            if ga == 'exon_specific':
                fractions['exon'] += overlap
            elif ga == 'intron_specific':
                fractions['intron'] += overlap
            elif ga == 'intergenic_specific':
                fractions['intergenic'] += overlap

    # Normalize to fractions
    for key in fractions:
        fractions[key] /= size

    return fractions


def identify_arrays(clusters, gene_regions):
    """
    Find contiguous stretches of EXON-DENSE clusters.
    These are the D4Z4 arrays. Everything else is transition.

    Args:
        clusters: List of (start, end) tuples from get_bsr_units()
        gene_regions: List of (start, end, annotation) tuples

    Returns:
        List of array dicts with:
            - start, end: array boundaries
            - elements: list of (start, end) tuples for exon-dense clusters
            - element_count: number of elements
    """
    if not clusters:
        return []

    # Classify each cluster as exon-dense or not
    cluster_info = []
    for start, end in clusters:
        size = end - start
        fractions = get_cluster_annotation_fractions(start, end, gene_regions)
        is_exon_dense = (fractions['exon'] >= EXON_DENSE_THRESHOLD and
                         size >= MIN_ARRAY_ELEMENT_SIZE)
        cluster_info.append({
            'start': start,
            'end': end,
            'size': size,
            'fractions': fractions,
            'is_exon_dense': is_exon_dense
        })

    # Group consecutive exon-dense clusters into arrays
    arrays = []
    current_array_elements = []

    for i, info in enumerate(cluster_info):
        if info['is_exon_dense']:
            # Check gap to previous element in current array
            if current_array_elements:
                gap = info['start'] - current_array_elements[-1]['end']
                if gap > MAX_INTER_ELEMENT_GAP:
                    # Gap too large - finalize current array, start new one
                    if len(current_array_elements) >= 1:  # Keep even small arrays for now
                        arrays.append(_finalize_array(current_array_elements))
                    current_array_elements = []

            current_array_elements.append(info)
        else:
            # Not exon-dense - finalize current array if exists
            if current_array_elements:
                arrays.append(_finalize_array(current_array_elements))
                current_array_elements = []

    # Don't forget the last array
    if current_array_elements:
        arrays.append(_finalize_array(current_array_elements))

    return arrays


def _finalize_array(elements):
    """Create array dict from list of element info dicts."""
    return {
        'start': elements[0]['start'],
        'end': elements[-1]['end'],
        'elements': [(e['start'], e['end']) for e in elements],
        'element_count': len(elements)
    }


def find_plam_after_array(array_end, clusters, gene_regions, region_regions):
    """
    Look for pLAM terminal immediately after an array.

    pLAM characteristics:
    - Large (>= PLAM_MIN_SIZE)
    - Predominantly intron-rich (< PLAM_MAX_EXON_FRACTION exon)
    - Contains bsat
    - Within reasonable distance after array end

    Args:
        array_end: End position of the array
        clusters: All BSR clusters
        gene_regions: Gene annotations
        region_regions: Region annotations (for ct_specific check)

    Returns:
        dict with pLAM info or None if not found
    """
    # Find clusters immediately after array end
    candidates = [(s, e) for s, e in clusters
                  if s >= array_end - 500 and s <= array_end + 2000]

    for start, end in sorted(candidates):
        size = end - start
        if size < PLAM_MIN_SIZE:
            continue

        fractions = get_cluster_annotation_fractions(start, end, gene_regions)

        # pLAM should be intron-rich or intergenic, NOT exon-rich
        if fractions['exon'] < PLAM_MAX_EXON_FRACTION:
            return {
                'start': start,
                'end': end,
                'size': size,
                'type': 'pLAM',
                'fractions': fractions
            }

    return None


def find_ct_specific_terminal(array_end, region_regions):
    """
    Look for ct_specific terminal near array end.

    Args:
        array_end: End position of the array
        region_regions: Region annotations

    Returns:
        dict with terminal info or None if not found
    """
    for rs, re, ra in region_regions:
        if ra == 'ct_specific':
            # Check if ct_specific is near array end
            if rs <= array_end + CT_TERMINAL_SEARCH_DIST and re >= array_end - 1000:
                return {
                    'start': rs,
                    'end': re,
                    'type': 'ct_specific'
                }

    return None


def classify_arrays(arrays, clusters, gene_regions, region_regions):
    """
    Classify arrays as canonical or degraded based on terminal presence.

    Canonical = has pLAM or ct_specific terminal
    Degraded = no terminal (or too few elements)

    Args:
        arrays: List of array dicts from identify_arrays()
        clusters: All BSR clusters
        gene_regions: Gene annotations
        region_regions: Region annotations

    Returns:
        Updated arrays list with classification info
    """
    for i, arr in enumerate(arrays):
        arr['index'] = i + 1

        # Check for pLAM
        plam = find_plam_after_array(arr['end'], clusters, gene_regions, region_regions)
        arr['plam'] = plam

        # Check for ct_specific terminal
        ct_terminal = find_ct_specific_terminal(arr['end'], region_regions)
        arr['ct_terminal'] = ct_terminal

        # Determine if array has a valid terminal
        arr['has_terminal'] = plam is not None or ct_terminal is not None

        # Classify: canonical if >= MIN_ARRAY_ELEMENTS and has terminal
        if arr['element_count'] >= MIN_ARRAY_ELEMENTS and arr['has_terminal']:
            arr['type'] = 'canonical'
        elif arr['element_count'] >= MIN_ARRAY_ELEMENTS:
            arr['type'] = 'degraded'
        else:
            arr['type'] = 'fragment'  # Too few elements to be a real array

    return arrays


def get_cluster_classification(cluster_start, cluster_end, arrays, gene_regions):
    """
    Classify a single cluster based on array membership.

    Returns:
        dict with:
            - type: 'array_element', 'plam', or 'transition'
            - array_idx: index of containing array (if applicable)
            - element_num: element number within array (if applicable)
            - array_type: 'canonical', 'degraded', or 'fragment' (if in array)
    """
    # Check if cluster is part of any array
    for arr in arrays:
        # Check if this cluster is an array element
        for elem_idx, (es, ee) in enumerate(arr['elements']):
            if abs(cluster_start - es) < 100 and abs(cluster_end - ee) < 100:
                return {
                    'type': 'array_element',
                    'array_idx': arr['index'],
                    'element_num': elem_idx + 1,
                    'array_type': arr['type']
                }

        # Check if this cluster is the pLAM for this array
        if arr['plam']:
            ps, pe = arr['plam']['start'], arr['plam']['end']
            if abs(cluster_start - ps) < 100 and abs(cluster_end - pe) < 100:
                return {
                    'type': 'plam',
                    'array_idx': arr['index'],
                    'array_type': arr['type']
                }

    # Not in any array - it's a transition element
    return {
        'type': 'transition',
        'array_idx': None,
        'element_num': None,
        'array_type': None
    }


# Keep the old function name as an alias for compatibility
def get_bsr_clusters_smart(bsat_regions, merge_threshold=None, gene_regions=None):
    """
    Compatibility wrapper - now uses conservative merging.
    """
    return get_bsr_units(bsat_regions, merge_threshold=500)


# =============================================================================
# Drawing Functions
# =============================================================================

def get_color(annotation, annot_type, color_maps):
    """Get color for an annotation from loaded color maps"""
    if annot_type in color_maps and annotation in color_maps[annot_type]:
        return color_maps[annot_type][annotation]
    return DEFAULT_COLOR


def draw_annotation_track(d, regions, y_pos, track_height, x_scale, x_offset,
                          annot_type, label, text_color, left_margin, color_maps, used_colors,
                          style=None):
    """Draw a single annotation track and track which colors are used"""
    # Scale font to match track height
    label_font = _style(style, 'track_label_font', max(7, int(track_height * 0.85)))
    # Draw track label
    if _style(style, 'track_labels', True):
        d.append(draw.Text(
            label, font_size=label_font, x=left_margin - _style(style, 'track_label_gap', 10),
            y=y_pos + track_height/2 + label_font * 0.35,
            fill=text_color, font_family=_style(style, 'font_family', 'sans-serif'), text_anchor='end'
        ))

    # Draw regions
    for start, end, annotation in regions:
        x = x_offset + start * x_scale
        width = (end - start) * x_scale
        color = get_color(annotation, annot_type, color_maps)

        # Track used colors
        if annot_type not in used_colors:
            used_colors[annot_type] = {}
        if annotation not in used_colors[annot_type]:
            used_colors[annot_type][annotation] = color

        d.append(draw.Rectangle(
            x, y_pos, max(width, 1), track_height,
            fill=color, stroke='none'
        ))


def draw_bsr_clusters_from_tsv(d, clusters_info, y_pos, track_height, x_scale, x_offset,
                               text_color, left_margin, style=None):
    """Draw BSR cluster markers using pre-computed classifications from TSV.

    This is the simplified implementation that reads cluster classifications
    directly from the analysis TSV, avoiding duplicate logic.

    Args:
        clusters_info: List of dicts with pre-computed classification:
            - start, end, size
            - type: 'array_element', 'plam', or 'transition'
            - array_idx, element_num, array_type

    Colors:
    - Light green (COLOR_CANONICAL_BSAT): Canonical array with β-satellite terminal (numbered)
    - Light blue (COLOR_CANONICAL_CT): Canonical array lacking β-satellite terminal (numbered)
    - Light red (COLOR_DEGRADED_ARRAY): Degraded array elements (numbered)
    - Gold (COLOR_TRANSITION): Transition regions (marked T)
    - Dark green (COLOR_PLAM): β-satellite terminal elements (marked B)

    Returns:
        tuple: (label_y, arrays_for_arrows) for further drawing
    """
    # Scale fonts to match track height
    label_font = _style(style, 'track_label_font', max(7, int(track_height * 0.85)))
    element_font = _style(style, 'element_font', max(5, int(track_height * 0.5)))
    font_family = _style(style, 'font_family', 'sans-serif')
    # Draw track label
    if _style(style, 'track_labels', True):
        d.append(draw.Text(
            'D4Z4 arrays', font_size=label_font, x=left_margin - _style(style, 'track_label_gap', 10),
            y=y_pos + track_height/2 + label_font * 0.35,
            fill=text_color, font_family=font_family, text_anchor='end'
        ))

    if not clusters_info:
        return y_pos + track_height + int(15 * track_height / 15), []

    label_y = y_pos + track_height + _style(style, 'element_label_offset', int(10 * track_height / 15))

    # Pre-compute element counts per array for label thinning
    array_element_counts = {}
    array_extent = {}   # array_idx -> (first start, last end), for spacing-aware label thinning
    for cls in clusters_info:
        if cls.get('type') == 'array_element' and cls.get('array_idx'):
            idx = cls['array_idx']
            array_element_counts[idx] = array_element_counts.get(idx, 0) + 1
            lo, hi = array_extent.get(idx, (cls['start'], cls['end']))
            array_extent[idx] = (min(lo, cls['start']), max(hi, cls['end']))

    def _label_spacing(total):
        if total <= 10: return 1
        elif total <= 20: return 2
        elif total <= 40: return 4
        elif total <= 60: return 5
        elif total <= 100: return 10
        else: return 20

    # Track arrays for arrow drawing
    array_info = {}  # array_idx -> {start, end, type, has_terminal}

    # Track last T label position to thin out consecutive transition labels
    last_t_label_x = -999

    for cls in clusters_info:
        start = cls['start']
        end = cls['end']
        x = x_offset + start * x_scale
        width = (end - start) * x_scale

        cls_type = cls['type']
        array_type = cls.get('array_type')
        array_idx = cls.get('array_idx')
        element_num = cls.get('element_num')
        terminal_type = cls.get('terminal_type')

        # Determine color and label
        if cls_type == 'array_element':
            if array_type == 'canonical':
                # Canonical arrays: green if b-sat terminal, blue if ct_specific terminal
                if terminal_type == 'plam':
                    color = COLOR_CANONICAL_BSAT
                elif terminal_type == 'ct_specific':
                    color = COLOR_CANONICAL_CT
                else:
                    color = COLOR_CANONICAL_BSAT  # Default to green
            elif array_type == 'degraded':
                color = COLOR_DEGRADED_ARRAY
            else:
                # Fragments (< MIN_ARRAY_ELEMENTS) are colored as transition
                color = COLOR_TRANSITION
            # Thin out element labels for large arrays
            if element_num:
                total = array_element_counts.get(array_idx, 1)
                spacing = _label_spacing(total)
                # Optional: widen the label interval until neighbouring labels are at least
                # label_min_gap_em font sizes apart on the page (arrays drawn at a compressed scale)
                min_gap = _style(style, 'label_min_gap_em', 0) * element_font
                if min_gap and total > 1:
                    lo, hi = array_extent[array_idx]
                    pitch = (hi - lo) * x_scale / total
                    spacing = next((k for k in (1, 2, 4, 5, 10, 20, 50, 100)
                                    if k >= spacing and k * pitch >= min_gap), 100)
                if element_num == total or element_num % spacing == 0:
                    label_text = str(element_num)
                else:
                    label_text = ''
                # Optional per-figure removal of crowded labels, as (array_idx, element_num) pairs
                if (array_idx, element_num) in _style(style, 'hidden_element_labels', ()):
                    label_text = ''
                # Optional automatic version of the rule applied by hand in the published figures:
                # keep the array's total and drop the label next to it when they would collide
                if label_text and element_num != total and _style(style, 'auto_hide_crowded_labels', False) \
                        and total - element_num <= spacing / 2:
                    label_text = ''
            else:
                label_text = ''

            # Track array info for arrows (only valid arrays, not fragments)
            if array_idx and array_type in ('canonical', 'degraded'):
                if array_idx not in array_info:
                    array_info[array_idx] = {'start': start, 'end': end, 'type': array_type, 'has_terminal': False, 'terminal_type': terminal_type}
                else:
                    array_info[array_idx]['start'] = min(array_info[array_idx]['start'], start)
                    array_info[array_idx]['end'] = max(array_info[array_idx]['end'], end)

        elif cls_type == 'plam':
            color = COLOR_PLAM
            label_text = ''  # No label for terminus
            # Mark that this array has a terminal and extend end to include pLAM
            if array_idx and array_idx in array_info:
                array_info[array_idx]['has_terminal'] = True
                array_info[array_idx]['end'] = max(array_info[array_idx]['end'], end)

        elif cls_type == 'solitary_plam':
            color = COLOR_SOLITARY_PLAM
            label_text = ''  # No label for solitary terminus

        else:  # transition
            color = COLOR_TRANSITION
            # Thin out T labels: only show if >15px from last T label
            t_label_x = x_offset + (start + end) / 2 * x_scale
            if t_label_x - last_t_label_x > 15:
                label_text = 'T'
                last_t_label_x = t_label_x
            else:
                label_text = ''

        # Draw cluster box
        d.append(draw.Rectangle(
            x, y_pos, max(width, 2), track_height,
            fill=color, stroke='none'
        ))

        # Draw cluster label below
        if label_text:
            d.append(draw.Text(
                label_text, font_size=element_font,
                x=x + width/2, y=label_y,
                fill=text_color, font_family=font_family, text_anchor='middle'
            ))

    # Build arrays_for_arrows with element sizes
    arrays_for_arrows = []
    for idx, info in sorted(array_info.items()):
        # Get element sizes for this array
        element_sizes = [c.get('size', 0) for c in clusters_info
                        if c.get('array_idx') == idx and c.get('type') == 'array_element']
        # Calculate median element size
        median_size = 0
        if element_sizes:
            sorted_sizes = sorted(element_sizes)
            mid = len(sorted_sizes) // 2
            if len(sorted_sizes) % 2 == 0:
                median_size = (sorted_sizes[mid - 1] + sorted_sizes[mid]) // 2
            else:
                median_size = sorted_sizes[mid]
        arrays_for_arrows.append({
            'start': info['start'],
            'end': info['end'],
            'is_functional': info['type'] in ('canonical', 'degraded'),
            'has_plam': info['has_terminal'],
            'terminal_type': info.get('terminal_type'),
            'index': idx,
            'element_count': len(element_sizes),
            'median_element_size': median_size
        })

    return label_y + 5, arrays_for_arrows


# Compatibility wrapper that can use either TSV data or compute on the fly
def draw_bsr_clusters(d, clusters, y_pos, track_height, x_scale, x_offset,
                      main_array_start, text_color, left_margin, plam_detected=False, plam_size=None,
                      array_boundaries=None, gene_regions=None, repeat_regions=None,
                      region_regions=None, clusters_info=None, style=None):
    """Draw BSR clusters - uses TSV data if available, otherwise computes.

    If clusters_info is provided (from TSV), uses pre-computed classifications.
    Otherwise falls back to computing classifications on the fly.
    """
    if clusters_info:
        # Use pre-computed classifications from TSV
        return draw_bsr_clusters_from_tsv(
            d, clusters_info, y_pos, track_height, x_scale, x_offset,
            text_color, left_margin, style=style
        )
    else:
        # Compute classifications on the fly (fallback)
        arrays = identify_arrays(clusters, gene_regions)
        region_regions_list = region_regions if region_regions else []
        arrays = classify_arrays(arrays, clusters, gene_regions, region_regions_list)
        cluster_classifications = []
        for start, end in clusters:
            cls = get_cluster_classification(start, end, arrays, gene_regions)
            cluster_classifications.append({
                'start': start,
                'end': end,
                'type': cls['type'],
                'array_idx': cls['array_idx'],
                'element_num': cls['element_num'],
                'array_type': cls['array_type']
            })
        return draw_bsr_clusters_from_tsv(
            d, cluster_classifications, y_pos, track_height, x_scale, x_offset,
            text_color, left_margin, style=style
        )


def draw_array_labels(d, array_boundaries, y_pos, x_scale, x_offset, text_color, plam_data=None, track_height=15,
                      style=None):
    """Draw array labels below element indices, numbered from telomere (right)

    Args:
        array_boundaries: List of dicts with:
            - start: tract/array start position
            - end: array end position
            - is_functional: True for functional arrays, False for degenerate
            - median_element_size: median size of elements in bp
        y_pos: Y position for labels (below element indices)
        plam_data: Optional dict keyed by array_index with plam_start, plam_end, plam_functional
        style: optional overrides. style['array_labels'] == 'figure' draws the manuscript layout:
            pLAM markers directly under the array track, no label for a single array, and for
            several arrays a bar under each array's element numbers with '#N' beneath it and
            one 'Array' caption to the left.
    """
    if not array_boundaries:
        return
    if _style(style, 'array_labels', 'classic') == 'figure':
        return _draw_array_labels_figure(d, array_boundaries, y_pos, x_scale, x_offset, text_color,
                                         plam_data, track_height, style)

    # Scale fonts and arrow sizes to track height
    label_font = max(6, int(track_height * 0.6))
    line_width = 1
    arrow_half_w = max(2, track_height * 0.125)
    arrow_h = max(3, track_height * 0.25)

    # Sort arrays by position (left to right) and assign numbers from right (telomere)
    sorted_arrays = sorted(array_boundaries, key=lambda a: a.get('start', 0))
    total_arrays = len(sorted_arrays)

    for i, arr in enumerate(sorted_arrays):
        if isinstance(arr, dict):
            arr_start = arr.get('array_start') or arr['start']
            arr_end = arr.get('end', arr_start)
            is_functional = arr.get('is_functional', True)
            median_size = arr.get('median_element_size', 0)
        else:
            arr_start = arr[0]
            arr_end = arr[1] if len(arr) > 1 else arr_start
            is_functional = True
            median_size = 0

        arr_num = total_arrays - i
        arr_center = (arr_start + arr_end) / 2
        x = x_offset + arr_center * x_scale

        has_terminal = arr.get('has_plam', False) if isinstance(arr, dict) else False
        terminal_type = arr.get('terminal_type') if isinstance(arr, dict) else None
        if is_functional and has_terminal:
            if terminal_type == 'ct_specific':
                label_color = COLOR_CANONICAL_CT
            else:
                label_color = COLOR_CANONICAL_BSAT
        elif is_functional:
            label_color = COLOR_DEGRADED_ARRAY
        else:
            label_color = COLOR_TRANSITION

        # Only the leftmost array (highest number) gets "Array #N", others just "#N"
        if total_arrays > 1 and arr_num == total_arrays:
            label = f"Array #{arr_num}"
        elif total_arrays > 1:
            label = f"#{arr_num}"
        else:
            label = f"Array #{arr_num}"

        x_start = x_offset + arr_start * x_scale
        x_end = x_offset + arr_end * x_scale
        line_y = y_pos - 2
        d.append(draw.Line(
            x_start, line_y, x_end, line_y,
            stroke=text_color, stroke_width=line_width
        ))

        d.append(draw.Text(
            label, font_size=label_font,
            x=x, y=y_pos + label_font * 0.7,
            fill=text_color, font_family='sans-serif', text_anchor='middle'
        ))

        if plam_data and arr_num in plam_data:
            plam_info = plam_data[arr_num]
            plam_functional = plam_info.get('plam_functional')
            if pd.notna(plam_functional) and plam_functional in ('func_Plam', 'nonfunc_Plam'):
                plam_mid = (plam_info['plam_start'] + plam_info['plam_end']) / 2
                arrow_x = x_offset + plam_mid * x_scale
                arrow_color = COLOR_PLAM_FUNCTIONAL if plam_functional == 'func_Plam' else COLOR_PLAM_NONFUNCTIONAL
                d.append(draw.Lines(
                    arrow_x, line_y - arrow_h / 2,
                    arrow_x - arrow_half_w, line_y + arrow_h / 2,
                    arrow_x + arrow_half_w, line_y + arrow_h / 2,
                    close=True, fill=arrow_color, stroke='none'
                ))


def _draw_array_labels_figure(d, array_boundaries, y_pos, x_scale, x_offset, text_color,
                              plam_data, track_height, style):
    """Array labels in the manuscript-figure layout (see draw_array_labels)."""
    font_family = _style(style, 'font_family', 'sans-serif')
    label_font = _style(style, 'array_label_font', max(6, int(track_height * 0.6)))
    arrow_half_w = _style(style, 'plam_marker_half_width', max(2, track_height * 0.125))
    arrow_h = _style(style, 'plam_marker_height', max(3, track_height * 0.25))
    # y_pos is where the classic layout puts its bar: element-number baseline + 5
    track_bottom = y_pos - 5 - _style(style, 'element_label_offset', int(10 * track_height / 15))
    bar_y = y_pos + _style(style, 'array_bar_offset', 0)
    caption_y = bar_y + _style(style, 'array_caption_offset', label_font)

    sorted_arrays = sorted(array_boundaries, key=lambda a: a.get('start', 0))
    total_arrays = len(sorted_arrays)
    caption_left = None
    for i, arr in enumerate(sorted_arrays):
        arr_num = total_arrays - i
        arr_start = arr.get('array_start') or arr['start']
        arr_end = arr.get('end', arr_start)
        x_start = x_offset + arr_start * x_scale
        x_end = x_offset + arr_end * x_scale
        if total_arrays > 1:
            d.append(draw.Line(x_start, bar_y, x_end, bar_y, stroke=text_color,
                               stroke_width=_style(style, 'stroke_width', 1)))
            label = f'#{arr_num}'
            x_mid = (x_start + x_end) / 2
            d.append(draw.Text(label, font_size=label_font, x=x_mid, y=caption_y,
                               fill=text_color, font_family=font_family, text_anchor='middle'))
            # '#N' is about 0.56 em per character in Arial/Helvetica
            left = x_mid - len(label) * label_font * 0.28
            caption_left = left if caption_left is None else min(caption_left, left)

        if plam_data and arr_num in plam_data:
            plam_info = plam_data[arr_num]
            plam_functional = plam_info.get('plam_functional')
            if pd.notna(plam_functional) and plam_functional in ('func_Plam', 'nonfunc_Plam'):
                plam_mid = (plam_info['plam_start'] + plam_info['plam_end']) / 2
                arrow_x = x_offset + plam_mid * x_scale
                arrow_color = COLOR_PLAM_FUNCTIONAL if plam_functional == 'func_Plam' else COLOR_PLAM_NONFUNCTIONAL
                top = track_bottom + _style(style, 'plam_marker_gap', 1)
                d.append(draw.Lines(
                    arrow_x, top,
                    arrow_x - arrow_half_w, top + arrow_h,
                    arrow_x + arrow_half_w, top + arrow_h,
                    close=True, fill=arrow_color, stroke='none'
                ))

    if caption_left is not None:
        # 'Array' starts at the track's left edge, unless that would run into the first '#N'
        gap = _style(style, 'array_caption_gap', label_font * 0.5)
        caption_width = 2.39 * label_font  # 'Array' in Arial/Helvetica
        if x_offset + caption_width + gap <= caption_left:
            x, anchor = x_offset, 'start'
        else:
            x, anchor = caption_left - gap, 'end'
        d.append(draw.Text('Array', font_size=label_font, x=x, y=caption_y, fill=text_color,
                           font_family=font_family, text_anchor=anchor))


def draw_scale_bar(d, x_offset, y_pos, x_scale, text_color, track_height=15):
    """Draw a scale bar"""
    scale_bp = 10000
    scale_width = scale_bp * x_scale
    bar_weight = max(2, track_height / 8)
    tick_size = max(5, track_height / 3)
    label_font = max(8, int(track_height * 0.85))

    d.append(draw.Line(
        x_offset, y_pos, x_offset + scale_width, y_pos,
        stroke=text_color, stroke_width=bar_weight
    ))

    d.append(draw.Line(x_offset, y_pos - tick_size, x_offset, y_pos + tick_size, stroke=text_color, stroke_width=bar_weight))
    d.append(draw.Line(x_offset + scale_width, y_pos - tick_size, x_offset + scale_width, y_pos + tick_size, stroke=text_color, stroke_width=bar_weight))

    d.append(draw.Text(
        '10 kb', font_size=label_font, x=x_offset + scale_width/2, y=y_pos + tick_size + label_font,
        fill=text_color, font_family='sans-serif', text_anchor='middle'
    ))


def draw_legend(d, x_pos, y_pos, text_color, used_colors, color_maps, plot_width, track_height=15):
    """Draw vertical column legend with bold section headers and items underneath.

    Layout: columns side by side, each column has a bold header + items stacked below.
    Font sizes scale with track_height to match the annotation track labels.
    """
    # Scale all sizes to match track height
    header_font = max(8, int(track_height * 0.7))
    item_font = max(6, int(track_height * 0.55))
    swatch_size = max(6, int(track_height * 0.5))
    item_height = max(10, int(track_height * 0.8))
    col_width = max(140, int(track_height * 10))

    # Section labels for each annotation type
    type_labels = {
        'region': 'Satellites',
        'gene': 'Genes',
        'repeat': 'Repeats',
        'subtelomeric': 'Subtelomere',
    }

    # Build sections: list of (title, [(label, color, marker_type), ...])
    sections = []

    # Collect "novel" entries from all annotation types
    common_items = []
    for annot_type in used_colors:
        if 'novel' in used_colors[annot_type]:
            common_items.append(('novel', used_colors[annot_type]['novel'], 'rect'))
    if common_items:
        sections.append(('Common', common_items))

    # KaryoScope annotation types
    for annot_type in ['region', 'gene', 'repeat', 'subtelomeric']:
        if annot_type not in used_colors or not used_colors[annot_type]:
            continue
        type_colors = {k: v for k, v in used_colors[annot_type].items() if k != 'novel'}
        if not type_colors:
            continue
        items = []
        for feature, color in sorted(type_colors.items()):
            display_name = feature.replace('_specific', '').replace('_multigroup1', '')
            items.append((display_name, color, 'rect'))
        sections.append((type_labels[annot_type], items))

    # D4Z4 Arrays
    sections.append(('D4Z4 Arrays', [
        ('C+ (bsat terminal)', COLOR_CANONICAL_BSAT, 'rect'),
        ('C- (ct terminal)', COLOR_CANONICAL_CT, 'rect'),
        ('Degraded array', COLOR_DEGRADED_ARRAY, 'rect'),
    ]))

    # Other elements
    sections.append(('Other Elements', [
        ('Transition (T)', COLOR_TRANSITION, 'rect'),
        ('bsat terminal', COLOR_PLAM, 'rect'),
        ('Solitary bsat', COLOR_SOLITARY_PLAM, 'rect'),
    ]))

    # pLAM markers
    sections.append(('pLAM Sequence', [
        ('Functional pLAM', COLOR_PLAM_FUNCTIONAL, 'triangle'),
        ('Non-functional pLAM', COLOR_PLAM_NONFUNCTIONAL, 'triangle'),
    ]))

    # Draw columns: distribute sections across columns
    col_x = x_pos
    max_y = y_pos

    for section_title, items in sections:
        col_y = y_pos

        # Section header (bold)
        d.append(draw.Text(
            section_title, font_size=header_font, x=col_x, y=col_y,
            fill=text_color, font_family='sans-serif',
            text_anchor='start', font_weight='bold'
        ))
        col_y += header_font + 2

        # Items
        for label, color, marker in items:
            col_y += item_height

            if marker == 'triangle':
                cx = col_x + 3 + swatch_size / 2
                cy = col_y - swatch_size / 2
                half = swatch_size / 2
                d.append(draw.Lines(
                    cx, cy - half,
                    cx - half, cy + half,
                    cx + half, cy + half,
                    close=True, fill=color, stroke='none'
                ))
            else:
                d.append(draw.Rectangle(
                    col_x + 3, col_y - swatch_size, swatch_size, swatch_size,
                    fill=color, stroke=text_color, stroke_width=0.5
                ))

            d.append(draw.Text(
                label, font_size=item_font, x=col_x + swatch_size + 8,
                y=col_y - swatch_size / 2 + 1,
                fill=text_color, font_family='sans-serif', text_anchor='start'
            ))

        max_y = max(max_y, col_y)
        col_x += col_width

    return max_y + 15


def parse_tract_info(tract_info_str):
    """Parse tract_info string into list of (start, end) tuples in bp.

    Example input: "Tract1:140-175kb(~11rpt); Tract2:200-220kb(~6rpt)"
    Returns: [(140000, 175000), (200000, 220000)]
    """
    if not tract_info_str or pd.isna(tract_info_str):
        return []

    tracts = []
    for part in tract_info_str.split('; '):
        # Parse "Tract1:140-175kb(~11rpt)"
        if ':' in part and '-' in part:
            try:
                range_part = part.split(':')[1].split('kb')[0]  # "140-175"
                start_kb, end_kb = range_part.split('-')
                tracts.append((int(start_kb) * 1000, int(end_kb) * 1000))
            except (ValueError, IndexError):
                continue
    return tracts


def adaptive_window_size(raw_max_pos, analysis_row):
    """--adaptive-window span: bsat_start to the telomeric end of the contig, rounded up to 50 kb."""
    bsat_start = analysis_row.get('bsat_start', 0)
    if pd.notna(bsat_start) and float(bsat_start) > 0:
        span = raw_max_pos - int(bsat_start)
        return math.ceil(span / 50000) * 50000
    return math.ceil(raw_max_pos / 50000) * 50000


def plot_single_haplotype(d, hap_data, hap_name, analysis_row, y_offset,
                          left_margin, plot_width, text_color, background_color,
                          color_maps, used_colors, adaptive_window=False, plam_data=None,
                          height_ratio=1.0, style=None):
    """Plot a single haplotype with all annotation tracks.

    style: optional dict of drawing overrides (fonts, label layout); None keeps the defaults.
    """

    # Work on a copy to avoid mutating cached data when flipping coordinates
    hap_data = {k: list(v) for k, v in hap_data.items()}

    track_height = int(15 * height_ratio)
    track_spacing = max(1, int(3 * height_ratio))
    tracks = ['region', 'gene', 'repeat', 'subtelomeric']
    track_labels = {
        'region': 'Satellites',
        'gene': 'Genes',
        'repeat': 'Repeats',
        'subtelomeric': 'Subtelomere',
    }

    # Calculate coordinate range from raw data
    all_positions = []
    for annot_type in tracks:
        for start, end, _ in hap_data.get(annot_type, []):
            all_positions.extend([start, end])

    if not all_positions:
        return y_offset + 50

    raw_max_pos = max(all_positions)

    # Adaptive window: flip coordinates so 0=telomeric end, clip from bsat_start to telomeric end
    if adaptive_window and analysis_row is not None:
        window_size = adaptive_window_size(raw_max_pos, analysis_row)
        # Optional larger window, e.g. to draw several haplotypes at one scale (never smaller)
        window_size = max(window_size, _style(style, 'window_size', 0) or 0)
        clip_start = max(0, raw_max_pos - window_size)

        # Flip and clip all track data: display_pos = raw_max_pos - genomic_pos
        for annot_type in tracks:
            flipped = []
            for start, end, annotation in hap_data.get(annot_type, []):
                if end <= clip_start or start >= raw_max_pos:
                    continue  # Outside window
                cstart = max(start, clip_start)
                cend = min(end, raw_max_pos)
                new_start = raw_max_pos - cend
                new_end = raw_max_pos - cstart
                flipped.append((new_start, new_end, annotation))
            hap_data[annot_type] = flipped

        min_pos = 0
        max_pos = window_size
    else:
        min_pos = min(all_positions)
        max_pos = max(all_positions)

    coord_range = max_pos - min_pos

    x_scale = plot_width / coord_range if coord_range > 0 else 1
    x_offset = left_margin - min_pos * x_scale

    # Get BSR clusters using smart merging (pass gene_regions to split at exon→intron boundaries)
    bsat_regions = [(s, e) for s, e, a in hap_data.get('region', []) if a == 'bsat_specific']
    gene_regions = hap_data.get('gene', [])
    bsr_clusters = get_bsr_clusters_smart(bsat_regions, gene_regions=gene_regions)

    # Get main array start from analysis results
    main_array_start = None
    repeat_count_main = None
    repeat_count_total = None
    terminal_feature = None
    num_tracts = 1
    arrays_info_str = None

    plam_detected = False
    plam_size = None
    has_any_terminal = False  # True if pLAM OR ct_specific terminal

    if analysis_row is not None:
        main_array_start = analysis_row.get('main_array_start')
        if pd.notna(main_array_start):
            main_array_start = int(main_array_start)
        else:
            main_array_start = None
        repeat_count_main = analysis_row.get('repeat_count_main')
        repeat_count_total = analysis_row.get('repeat_count_total')
        terminal_feature = analysis_row.get('terminal_feature')
        # Get pLAM detection info
        plam_detected = analysis_row.get('plam_detected', False)
        if pd.isna(plam_detected):
            plam_detected = False
        plam_size = analysis_row.get('plam_size')
        if pd.notna(plam_size):
            plam_size = int(plam_size)
        else:
            plam_size = None
        # Check for ANY terminal (pLAM or ct_specific) - determines canonical vs degraded
        if plam_detected:
            has_any_terminal = True
        elif terminal_feature and pd.notna(terminal_feature):
            if 'ct_specific' in str(terminal_feature) or 'pLAM' in str(terminal_feature):
                has_any_terminal = True
        # Get multi-tract info (for older TSV format)
        num_tracts = analysis_row.get('num_tracts', 1)
        if pd.isna(num_tracts):
            num_tracts = 1
        # Get arrays_info JSON (array-level classification)
        arrays_info_str = analysis_row.get('arrays_info')
        # Get clusters_info JSON (v4 format: per-cluster classification)
        clusters_info_str = analysis_row.get('clusters_info')

    # Parse clusters_info JSON (v4 format: per-cluster classification)
    clusters_info = None
    if analysis_row is not None:
        clusters_info_str = analysis_row.get('clusters_info')
        if clusters_info_str and pd.notna(clusters_info_str):
            try:
                clusters_info = json.loads(clusters_info_str)
            except (json.JSONDecodeError, TypeError) as e:
                print(f"Warning: Could not parse clusters_info for {hap_name}: {e}")
                clusters_info = None

    # Construct array_boundaries: List of dicts with is_functional flag
    # Parse arrays_info JSON if available (new format)
    array_boundaries = []
    if arrays_info_str and pd.notna(arrays_info_str):
        try:
            arrays_info = json.loads(arrays_info_str)
            # Find the last functional array (which would have the terminal)
            last_functional_idx = None
            for i, arr in enumerate(arrays_info):
                if arr.get('is_functional', True):
                    last_functional_idx = i

            for i, arr in enumerate(arrays_info):
                # Use has_any_terminal for the last functional array (ct_specific or pLAM)
                # This overrides the JSON's has_plam if we detected a terminal from terminal_feature
                arr_has_terminal = arr.get('has_plam', False)
                if i == last_functional_idx and has_any_terminal:
                    arr_has_terminal = True

                # Use plam_end if available to include terminal in array span
                arr_end = arr.get('plam_end') or arr['end']
                # Determine terminal_type for coloring
                if arr.get('has_plam'):
                    terminal_type = 'plam'
                elif arr.get('has_ct_terminal'):
                    terminal_type = 'ct_specific'
                else:
                    terminal_type = None
                array_boundaries.append({
                    'start': arr['start'],
                    'end': arr_end,
                    'array_start': arr.get('array_start') or arr['start'],
                    'transition_end': arr.get('transition_end'),
                    'is_functional': arr.get('type') != 'degraded',
                    'has_plam': arr_has_terminal,  # Now includes ct_specific terminals
                    'terminal_type': terminal_type,
                    'repeat_count': arr.get('repeat_count', 0),
                    'index': arr.get('index', 0),
                    'is_mosaic': arr.get('is_mosaic', False),
                    'sub_arrays': arr.get('sub_arrays', [])
                })
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warning: Could not parse arrays_info for {hap_name}: {e}")
            array_boundaries = []

    # Fallback: if no arrays_info, use main_array_start
    if not array_boundaries and main_array_start is not None and bsr_clusters:
        main_end = max(e for s, e in bsr_clusters)
        array_boundaries = [{
            'start': main_array_start,
            'end': main_end,
            'array_start': main_array_start,
            'is_functional': True,
            'has_plam': has_any_terminal,  # Use any terminal (pLAM or ct_specific)
            'repeat_count': repeat_count_main if pd.notna(repeat_count_main) else 0,
            'index': 1
        }]

    # Flip derived coordinate data for adaptive window mode
    if adaptive_window and analysis_row is not None:
        # Flip main_array_start
        if main_array_start is not None:
            main_array_start = raw_max_pos - main_array_start

        # Flip clusters_info coordinates
        if clusters_info:
            for cls in clusters_info:
                old_start, old_end = cls['start'], cls['end']
                cls['start'] = raw_max_pos - old_end
                cls['end'] = raw_max_pos - old_start

        # Flip array_boundaries coordinates
        for arr in array_boundaries:
            old_start, old_end = arr['start'], arr['end']
            arr['start'] = raw_max_pos - old_end
            arr['end'] = raw_max_pos - old_start
            if 'array_start' in arr and arr['array_start'] is not None:
                arr['array_start'] = raw_max_pos - arr['array_start']
            if 'transition_end' in arr and arr['transition_end'] is not None:
                arr['transition_end'] = raw_max_pos - arr['transition_end']
            # Flip sub_arrays if present
            for sub in arr.get('sub_arrays', []):
                if 'start' in sub and 'end' in sub:
                    old_s, old_e = sub['start'], sub['end']
                    sub['start'] = raw_max_pos - old_e
                    sub['end'] = raw_max_pos - old_s

        # Flip plam_data coordinates
        if plam_data:
            for idx, pinfo in plam_data.items():
                if pd.notna(pinfo.get('plam_start')) and pd.notna(pinfo.get('plam_end')):
                    old_start, old_end = pinfo['plam_start'], pinfo['plam_end']
                    pinfo['plam_start'] = raw_max_pos - old_end
                    pinfo['plam_end'] = raw_max_pos - old_start

    # Draw haplotype title
    title_y = y_offset
    # Parse haplotype name: "chr4_hap1_HG02129#1#CM085848.1"
    parts = hap_name.split('_')
    if len(parts) >= 3:
        chrom = parts[0]  # chr4 or chr10
        chrom_arm = chrom.replace('chr', '') + 'q'  # 4q or 10q
        sample_hap_contig = parts[2]  # HG02129#1#CM085848.1
    else:
        chrom_arm = ''
        sample_hap_contig = hap_name

    title = f"{sample_hap_contig} ({chrom_arm})"

    # Scale fonts to track height
    title_font = _style(style, 'title_font', max(9, int(track_height * 0.85)))
    tick_font = _style(style, 'tick_font', max(6, int(track_height * 0.55)))
    tick_size = max(5, track_height / 3)
    font_family = _style(style, 'font_family', 'sans-serif')
    line_w = _style(style, 'stroke_width', 1)

    d.append(draw.Text(
        title, font_size=title_font, x=left_margin, y=title_y,
        fill=text_color, font_family=font_family, font_weight='bold'
    ))

    # Draw coordinate axis
    axis_y = title_y + _style(style, 'title_axis_gap', title_font + 4)
    d.append(draw.Line(
        left_margin, axis_y, left_margin + plot_width, axis_y,
        stroke=text_color, stroke_width=line_w
    ))

    # Add tick marks — use 100kb interval for windows >300kb, otherwise 50kb
    window_span = max_pos - min_pos
    tick_interval = 100000 if window_span > 300000 else 50000

    def _tick_anchor(tick_pos, all_ticks):
        """Left-align first tick, right-align last, center others."""
        if tick_pos == all_ticks[0]:
            return 'start'
        elif tick_pos == all_ticks[-1]:
            return 'end'
        return 'middle'

    if adaptive_window and analysis_row is not None:
        all_ticks = list(range(0, int(max_pos) + 1, tick_interval))
        if _style(style, 'label_window_end', False) and all_ticks[-1] < max_pos:
            # The window end is not a multiple of the tick interval: label it too, keeping
            # the last regular tick right-aligned so the two labels do not collide
            last = all_ticks[-1]
            all_ticks.append(int(max_pos))
            for tick_pos in (last, int(max_pos)):
                tick_x = x_offset + tick_pos * x_scale
                d.append(draw.Line(tick_x, axis_y, tick_x, axis_y + tick_size, stroke=text_color, stroke_width=line_w))
                d.append(draw.Text(
                    f"{tick_pos//1000}kb", font_size=tick_font, x=tick_x, y=axis_y + tick_size + tick_font,
                    fill=text_color, font_family=font_family, text_anchor='end'
                ))
            all_ticks = all_ticks[:-2] + [None]  # None: placeholder so the rest keep their anchors
        for tick_pos in all_ticks:
            if tick_pos is None:
                continue
            tick_x = x_offset + tick_pos * x_scale
            anchor = _tick_anchor(tick_pos, all_ticks)
            d.append(draw.Line(tick_x, axis_y, tick_x, axis_y + tick_size, stroke=text_color, stroke_width=line_w))
            d.append(draw.Text(
                f"{tick_pos//1000}kb", font_size=tick_font, x=tick_x, y=axis_y + tick_size + tick_font,
                fill=text_color, font_family=font_family, text_anchor=anchor
            ))
    else:
        tick_start = (min_pos // tick_interval) * tick_interval
        all_ticks = [t for t in range(int(tick_start), int(max_pos) + tick_interval, tick_interval)
                     if t >= min_pos and t <= max_pos]
        for tick_pos in all_ticks:
            tick_x = x_offset + tick_pos * x_scale
            anchor = _tick_anchor(tick_pos, all_ticks)
            d.append(draw.Line(tick_x, axis_y, tick_x, axis_y + tick_size, stroke=text_color, stroke_width=line_w))
            d.append(draw.Text(
                f"{tick_pos//1000}kb", font_size=tick_font, x=tick_x, y=axis_y + tick_size + tick_font,
                fill=text_color, font_family=font_family, text_anchor=anchor
            ))

    # Draw annotation tracks
    current_y = axis_y + tick_size + tick_font + 5

    for track in tracks:
        regions = hap_data.get(track, [])
        draw_annotation_track(d, regions, current_y, track_height, x_scale, x_offset,
                            track, track_labels[track], text_color, left_margin,
                            color_maps, used_colors, style=style)

        # Track border (solid)
        d.append(draw.Rectangle(
            left_margin, current_y, plot_width, track_height,
            fill='none', stroke=text_color, stroke_width=line_w
        ))

        current_y += track_height + track_spacing

    # Draw BSR clusters track (no border for cleaner look)
    # Use clusters_info from TSV if available, otherwise compute on the fly
    gene_regions = hap_data.get('gene', [])
    region_regions = hap_data.get('region', [])
    repeat_regions = hap_data.get('repeat', [])
    label_y, valid_arrays = draw_bsr_clusters(d, bsr_clusters, current_y, track_height, x_scale, x_offset,
                                              main_array_start, text_color, left_margin,
                                              plam_detected=plam_detected, plam_size=plam_size,
                                              array_boundaries=array_boundaries,
                                              gene_regions=gene_regions,
                                              repeat_regions=repeat_regions,
                                              region_regions=region_regions,
                                              clusters_info=clusters_info, style=style)
    current_y += track_height + track_spacing

    # Draw array labels below element indices
    if valid_arrays:
        draw_array_labels(d, valid_arrays, label_y, x_scale, x_offset, text_color, plam_data=plam_data,
                          track_height=track_height, style=style)
        current_y = max(current_y, label_y + int(15 * track_height / 15))

    return current_y + 20


# =============================================================================
# Main Script
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate D4Z4 array visualization",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Selection options
    parser.add_argument("--haplotypes", nargs='+',
                        help="Specific haplotype names to plot")
    parser.add_argument("--samples", nargs='+',
                        help="Sample names to plot (plots both haplotypes)")
    parser.add_argument("--chromosome", default="chr4",
                        help="Chromosome (chr4 or chr10, default: chr4)")
    parser.add_argument("--max-samples", type=int, default=5,
                        help="Maximum samples to plot if no specific selection (default: 5)")

    # Filtering
    parser.add_argument("--fshd-only", action="store_true",
                        help="Only plot FSHD-range samples (<=20 repeats with pLAM-like)")
    parser.add_argument("--short-arrays", action="store_true",
                        help="Only plot short arrays (<=15 repeats)")

    # Output options
    parser.add_argument("--output", default="d4z4_arrays.svg",
                        help="Output SVG file (default: d4z4_arrays.svg)")
    parser.add_argument("--background", choices=["white", "black"], default="white",
                        help="Background color (default: white)")
    parser.add_argument("--width", type=int, default=1200,
                        help="Plot width in pixels (default: 1200)")
    parser.add_argument("--height-ratio", type=float, default=1.0,
                        help="Scale factor for per-haplotype height (track heights, spacing, labels). "
                             "Values <1 compress vertically, >1 expand (default: 1.0)")
    parser.add_argument("--aspect-ratio", type=float, default=None,
                        help="Target width:height aspect ratio per haplotype (e.g. 2.5). "
                             "Overrides --height-ratio by computing the needed scaling.")

    # Window options
    parser.add_argument("--adaptive-window", action="store_true",
                        help="Use sample-specific window size based on array_length, with 0=telomeric end")

    # Smoothing options
    parser.add_argument("--smooth", action="store_true",
                        help="Apply smoothing to annotation tracks (merges adjacent identical features)")
    parser.add_argument("--smooth-gap", type=int, default=500,
                        help="Maximum gap size to fill during smoothing (default: 500bp)")

    # Data paths
    parser.add_argument("--data-dir", default="data",
                        help="Data directory with BED files")
    parser.add_argument("--results-file", default="results/d4z4_array_analysis.tsv",
                        help="D4Z4 analysis results TSV")
    parser.add_argument("--colors-dir",
                        default="KaryoScope/resources/KS_human_CHM13",
                        help="Directory with KaryoScope color files")

    return parser.parse_args()


def main():
    args = parse_args()

    # Setup
    background_color = args.background
    text_color = '#000000' if background_color == 'white' else '#FFFFFF'
    left_margin = 120
    plot_width = args.width - left_margin - 50

    # Compute height_ratio from aspect_ratio if specified
    # Default per-haplotype height at ratio=1.0 is ~128px (4 tracks × 18 + BSR + labels)
    BASE_HAP_HEIGHT = 128
    if args.aspect_ratio is not None:
        target_height = args.width / args.aspect_ratio
        args.height_ratio = target_height / BASE_HAP_HEIGHT
        print(f"Aspect ratio {args.aspect_ratio}:1 → height_ratio={args.height_ratio:.2f} "
              f"(~{int(target_height)}px per haplotype)")

    print(f"D4Z4 Array Visualization")
    print(f"=" * 50)
    print(f"Data directory: {args.data_dir}")
    print(f"Results file: {args.results_file}")
    print(f"Colors directory: {args.colors_dir}")

    # Load color maps from KaryoScope files
    color_maps = load_all_colors(args.colors_dir)
    for annot_type, colors in color_maps.items():
        print(f"  Loaded {len(colors)} colors for {annot_type}")

    # Track which colors are actually used
    used_colors = {}

    # Load pLAM coordinates
    plam_coords_file = os.path.join(os.path.dirname(args.results_file), 'd4z4_array_coordinates.plam.tsv')
    plam_df = pd.read_csv(plam_coords_file, sep='\t') if os.path.exists(plam_coords_file) else None

    # Load analysis results
    analysis_df = load_analysis_results(args.results_file)
    if analysis_df is None:
        print(f"Warning: Could not load analysis results from {args.results_file}")
        analysis_df = pd.DataFrame()
    else:
        print(f"Loaded {len(analysis_df)} haplotypes from analysis results")

    # Determine which haplotypes to plot
    haplotypes_to_plot = []

    if args.haplotypes:
        haplotypes_to_plot = args.haplotypes
    elif args.samples:
        for sample in args.samples:
            if len(analysis_df) > 0:
                sample_haps = analysis_df[
                    (analysis_df['sample'] == sample) &
                    (analysis_df['chromosome'] == args.chromosome)
                ]['haplotype'].tolist()
                haplotypes_to_plot.extend(sample_haps)
    else:
        if len(analysis_df) > 0:
            filtered_df = analysis_df[analysis_df['chromosome'] == args.chromosome].copy()

            if args.fshd_only:
                filtered_df = filtered_df[
                    (filtered_df['repeat_count_main'] <= 20) &
                    (filtered_df['terminal_feature'] == 'pLAM-like')
                ]
            elif args.short_arrays:
                filtered_df = filtered_df[filtered_df['repeat_count_main'] <= 15]

            filtered_df = filtered_df.sort_values('repeat_count_main')
            haplotypes_to_plot = filtered_df['haplotype'].head(args.max_samples).tolist()

    if not haplotypes_to_plot:
        print("No haplotypes to plot!")
        return

    print(f"Plotting {len(haplotypes_to_plot)} haplotypes:")
    for hap in haplotypes_to_plot:
        print(f"  - {hap}")

    # Load data and optionally apply smoothing
    hap_data_cache = {}
    for hap in haplotypes_to_plot:
        # Extract chromosome from haplotype name (e.g., "chr4_hap1_..." or "chr10_hap2_...")
        hap_chrom = hap.split('_')[0] if '_' in hap else args.chromosome
        hap_data = get_haplotype_data(hap, args.data_dir, hap_chrom)
        if args.smooth:
            hap_data = smooth_all_tracks(hap_data, max_gap=args.smooth_gap)
        hap_data_cache[hap] = hap_data

    if args.smooth:
        print(f"Applied smoothing with max_gap={args.smooth_gap}bp")

    # Estimate height (scaled by height_ratio)
    per_hap_height = int(150 * args.height_ratio)
    estimated_height = len(haplotypes_to_plot) * per_hap_height + 400

    # Create drawing
    d = draw.Drawing(args.width, estimated_height)
    d.append(draw.Rectangle(0, 0, args.width, estimated_height, fill=background_color))

    # Plot each haplotype (no title - starts at top)
    current_y = 20
    for hap in haplotypes_to_plot:
        hap_data = hap_data_cache[hap]

        analysis_row = None
        if len(analysis_df) > 0:
            matches = analysis_df[analysis_df['haplotype'] == hap]
            if len(matches) > 0:
                analysis_row = matches.iloc[0]

        # Build pLAM data dict for this haplotype
        hap_plam_data = None
        if plam_df is not None:
            # Parse hap name: "chr4_hap1_HG02129#1#CM085848.1"
            hap_parts = hap.split('_')
            if len(hap_parts) >= 3:
                hap_chrom = hap_parts[0]
                hap_num = hap_parts[1]  # "hap1" or "hap2"
                hap_sample = hap_parts[2].split('#')[0]
                plam_matches = plam_df[
                    (plam_df['chromosome'] == hap_chrom) &
                    (plam_df['sample'] == hap_sample) &
                    (plam_df['hap'] == hap_num)
                ]
                if len(plam_matches) > 0:
                    hap_plam_data = {}
                    for _, row in plam_matches.iterrows():
                        if pd.notna(row.get('plam_start')):
                            hap_plam_data[int(row['array_index'])] = {
                                'plam_start': row['plam_start'],
                                'plam_end': row['plam_end'],
                                'plam_functional': row['plam_functional']
                            }

        current_y = plot_single_haplotype(
            d, hap_data, hap, analysis_row, current_y,
            left_margin, plot_width, text_color, background_color,
            color_maps, used_colors, adaptive_window=args.adaptive_window,
            plam_data=hap_plam_data,
            height_ratio=args.height_ratio
        )

    # Draw legend with only used colors
    legend_y = current_y + 20
    track_height_for_legend = int(15 * args.height_ratio)
    legend_end_y = draw_legend(d, left_margin, legend_y, text_color, used_colors, color_maps, plot_width, track_height=track_height_for_legend)

    # Add footer with generation timestamp (bottom center)
    footer_y = legend_end_y + 30
    d.append(draw.Text(
        f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
        font_size=8, x=args.width / 2, y=footer_y,
        fill=text_color, font_family='sans-serif', text_anchor='middle'
    ))

    # Adjust canvas height
    final_height = max(estimated_height, footer_y + 20)
    d.height = final_height

    d.save_svg(args.output)
    print(f"\nSaved to: {args.output}")


if __name__ == '__main__':
    main()
