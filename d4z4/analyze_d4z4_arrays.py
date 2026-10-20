#!/usr/bin/env python3
"""
D4Z4 Array Length and Repeat Unit Analysis for HPRC Pangenome Data (v4)

This script analyzes KaryoScope output files to:
1. Identify D4Z4 beta-satellite array boundaries for each haplotype
2. Calculate array length and estimated repeat unit count
3. Detect terminal β-satellite elements (pLAM proxy) or ct_specific terminals
4. Flag arrays exceeding the 250kb analysis window
5. Classify each cluster as: array_element, plam, or transition

IMPORTANT - Terminal Detection Limitation:
This pipeline detects terminal β-satellite elements as a PROXY for pLAM.
It identifies the structural presence of β-satellite (large, intron-rich bsat
cluster after array), NOT the actual pLAM polyadenylation signal sequence.
- Functional pLAM sequence: ATTAAA (enables DUX4 expression)
- Non-functional variant: ATCAAA (found in some chr10 haplotypes)
- C+ classification = β-satellite terminal detected, NOT confirmed functional pLAM

v4 Architecture (Array-First approach):
- Arrays are defined by contiguous exon-dense clusters
- Everything outside arrays is transition
- Terminal β-satellite is only detected immediately after arrays
- Per-cluster classification is output for visualization
"""

import os
import sys
import json
import logging
from datetime import datetime
from collections import defaultdict
import pandas as pd

# =============================================================================
# Constants
# =============================================================================

REPEAT_UNIT_SIZE = 3300  # D4Z4 repeat unit is ~3.3kb (for reference)
EXPECTED_BSR_SIZE = 2300  # Expected BSR cluster size in bp
WINDOW_SIZE = 250000  # 250kb analysis window
TRUNCATION_THRESHOLD = 1000  # If bsat starts within 1kb of position 0, likely truncated

# v4 Array Detection Thresholds (Array-First approach)
BSR_MERGE_THRESHOLD = 400      # Conservative merge for BSR units (preserves individual repeat units)
EXON_DENSE_THRESHOLD = 0.5     # Minimum exon fraction to be considered exon-dense (for pLAM detection)
MIN_ARRAY_ELEMENT_SIZE = 1500  # Minimum size for a "full" array element (used for counting)
ARRAY_GAP_THRESHOLD = 10000    # Gap > 10kb suggests separate arrays/tracts
MIN_ARRAY_ELEMENTS = 3         # Minimum elements for a valid array
PLAM_MIN_SIZE = 2500           # Minimum size for pLAM terminal
PLAM_MAX_EXON_FRACTION = 0.35  # pLAM must be < 35% exon
SOLITARY_PLAM_MIN_SIZE = 2500  # Minimum size for solitary pLAM (same as terminal)
CT_TERMINAL_SEARCH_DIST = 5000 # Search distance for ct_specific terminal

# Multi-tract detection (same as array gap)
TRACT_GAP_THRESHOLD = 10000    # Gap > 10kb suggests separate tracts

def load_bed_file(filepath):
    """Load BED file into list of tuples (chrom, start, end, annotation)"""
    regions = []
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

def get_haplotypes(regions):
    """Extract unique haplotype names from regions"""
    return sorted(set(r[0] for r in regions))


# =============================================================================
# v4 Core Functions: Array-First Detection Pipeline
# =============================================================================

def get_bsr_units(bsat_regions, merge_threshold=BSR_MERGE_THRESHOLD):
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
            cluster_end = end
        else:
            clusters.append((cluster_start, cluster_end))
            cluster_start, cluster_end = start, end

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

    for key in fractions:
        fractions[key] /= size

    return fractions


def identify_arrays(clusters, gene_regions):
    """
    Find contiguous stretches of EXON-DENSE D4Z4 clusters (arrays).
    Small gaps between exon-dense clusters may fragment arrays, which are
    merged in a subsequent step by merge_adjacent_arrays().

    Args:
        clusters: List of (start, end) tuples from get_bsr_units()
        gene_regions: List of (start, end, annotation) tuples

    Returns:
        Tuple of (arrays, cluster_info):
        - arrays: List of array dicts
        - cluster_info: List of cluster info dicts (for merge step)
    """
    if not clusters:
        return [], []

    # Build cluster info for all clusters
    cluster_info = []
    for start, end in clusters:
        size = end - start
        fractions = get_cluster_annotation_fractions(start, end, gene_regions)
        # Exon-dense if exon > threshold AND size is sufficient
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
    # Non-exon-dense clusters break arrays (will be recovered in merge step)
    arrays = []
    current_array_elements = []

    for info in cluster_info:
        if info['is_exon_dense']:
            # Check gap to previous element
            if current_array_elements:
                gap = info['start'] - current_array_elements[-1]['end']
                if gap > ARRAY_GAP_THRESHOLD:
                    # Very large gap - finalize current array
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

    return arrays, cluster_info


def _finalize_array(elements):
    """Create array dict from list of element info dicts."""
    return {
        'start': elements[0]['start'],
        'end': elements[-1]['end'],
        'elements': [(e['start'], e['end']) for e in elements],
        'element_count': len(elements),
        '_element_info': elements  # Keep full info for merging
    }


# Maximum gap between arrays to consider merging
ARRAY_MERGE_GAP = 4000  # If gap between arrays is < 4kb, consider merging
NOVEL_MERGE_THRESHOLD = 0.7  # If gap is > 70% novel, merge regardless of size
EXON_MERGE_THRESHOLD = 0.6  # If gap is > 60% exon, merge (gap contains array material)
COMBINED_MERGE_THRESHOLD = 0.8  # If gap is > 80% (exon + novel), merge
BOUNDARY_WINDOW_SIZE = 3000  # Check 3kb windows at gap boundaries
BOUNDARY_BLOCK_THRESHOLD = 0.3  # Block merge if boundary window has > 30% intergenic/intron


def get_gap_annotation_fractions(gap_start, gap_end, gene_regions):
    """Calculate fraction of gap covered by each annotation type.

    Novel sequence is calculated as the remainder - regions with NO gene annotation.
    This captures truly novel insertions not present in the reference.
    """
    size = gap_end - gap_start
    if size <= 0:
        return {'novel': 0, 'exon': 0, 'intron': 0, 'intergenic': 0}

    # Track absolute coverage in bp
    coverage = {'exon': 0, 'intron': 0, 'intergenic': 0}

    for gs, ge, ga in gene_regions:
        if ge > gap_start and gs < gap_end:
            overlap = min(ge, gap_end) - max(gs, gap_start)
            if ga == 'exon_specific':
                coverage['exon'] += overlap
            elif ga == 'intron_specific':
                coverage['intron'] += overlap
            elif ga == 'intergenic_specific':
                coverage['intergenic'] += overlap

    # Calculate fractions
    fractions = {key: val / size for key, val in coverage.items()}

    # Novel is the unannotated remainder
    annotated = fractions['exon'] + fractions['intron'] + fractions['intergenic']
    fractions['novel'] = max(0, 1.0 - annotated)

    return fractions


def merge_adjacent_arrays(arrays, all_clusters, gene_regions=None):
    """
    Merge adjacent arrays that are separated by:
    1. Small gaps with exon-dense clusters, OR
    2. Large insertions of predominantly "novel" sequence

    Novel insertions represent sequence not in the reference, which can
    interrupt an array without representing a true array boundary.
    """
    if len(arrays) <= 1:
        return arrays

    merged = []
    current = arrays[0]

    for next_arr in arrays[1:]:
        gap = next_arr['start'] - current['end']
        should_merge = False

        # Check 1: Small gap with exon-dense clusters
        if gap <= ARRAY_MERGE_GAP:
            gap_clusters = [c for c in all_clusters
                           if c['start'] >= current['end'] and c['end'] <= next_arr['start']]
            all_exon_dense = all(c['fractions']['exon'] > 0 for c in gap_clusters) if gap_clusters else True
            if all_exon_dense:
                should_merge = True

        # Check 2: Gap is predominantly novel sequence (insertion) or exon-dense (array material)
        # Use sliding window to detect intergenic/intron-rich regions anywhere in the gap
        if not should_merge and gene_regions:
            gap_start = current['end']
            gap_end = next_arr['start']
            gap_size = gap_end - gap_start

            # Slide window across gap to detect boundary-blocking regions
            has_boundary_block = False
            window_step = BOUNDARY_WINDOW_SIZE // 2  # 50% overlap

            if gap_size > BOUNDARY_WINDOW_SIZE:
                pos = gap_start
                while pos + BOUNDARY_WINDOW_SIZE <= gap_end:
                    window_fractions = get_gap_annotation_fractions(pos, pos + BOUNDARY_WINDOW_SIZE, gene_regions)
                    non_array_content = window_fractions['intergenic'] + window_fractions['intron']
                    if non_array_content >= BOUNDARY_BLOCK_THRESHOLD:
                        has_boundary_block = True
                        break
                    pos += window_step

            if not has_boundary_block:
                gap_fractions = get_gap_annotation_fractions(gap_start, gap_end, gene_regions)
                if gap_fractions['novel'] >= NOVEL_MERGE_THRESHOLD:
                    should_merge = True
                elif gap_fractions['exon'] >= EXON_MERGE_THRESHOLD:
                    # Gap contains exon-dense material - likely array elements that got separated
                    should_merge = True
                elif (gap_fractions['exon'] + gap_fractions['novel']) >= COMBINED_MERGE_THRESHOLD:
                    # Gap is mostly array-compatible (exon + novel) even if neither alone meets threshold
                    should_merge = True

        if should_merge:
            # Merge: combine elements from both arrays
            gap_clusters = [c for c in all_clusters
                           if c['start'] >= current['end'] and c['end'] <= next_arr['start']]
            merged_elements = current['_element_info'] + gap_clusters + next_arr['_element_info']
            current = _finalize_array(merged_elements)
        else:
            # No merge - finalize current and move to next
            merged.append(current)
            current = next_arr

    # Don't forget the last array
    merged.append(current)

    return merged


# Within-array merge parameters
# Expected D4Z4 repeat unit: ~3.3kb total, ~2.3kb exon-dense BSR element, ~1kb gap
EXPECTED_ELEMENT_SIZE = 2300  # Expected BSR element size
EXPECTED_GAP_SIZE = 1000      # Expected gap between elements
SIZE_TOLERANCE = 0.5          # Allow 50% deviation from expected
FRAGMENT_SIZE = 500           # Elements smaller than this are fragments
ELEMENT_WINDOW_SIZE = 3       # Number of elements to consider in sliding window


def merge_elements_within_arrays(arrays, all_clusters):
    """
    Apply pattern-based merging WITHIN each array to create consistent elements.

    Uses sliding window to ensure elements are equally sized and equally spaced:
    1. First pass: initial merge of obvious fragments (very small clusters)
    2. Calculate typical element size and spacing from the array
    3. Second pass: validate each element against expected patterns
    4. Merge fragments that fit the expected pattern

    Returns:
        Updated arrays with merged elements and element_size stats
    """
    merged_arrays = []

    for arr in arrays:
        arr_start = arr['start']
        arr_end = arr['end']

        # Get all clusters within this array's boundaries
        arr_clusters = [c for c in all_clusters
                       if c['start'] >= arr_start and c['end'] <= arr_end]

        if not arr_clusters:
            merged_arrays.append(arr)
            continue

        # Sort by position
        arr_clusters = sorted(arr_clusters, key=lambda x: x['start'])

        # First pass: calculate initial element sizes and gaps
        initial_elements = []
        for c in arr_clusters:
            initial_elements.append({
                'start': c['start'],
                'end': c['end'],
                'size': c['end'] - c['start']
            })

        # Calculate gaps between elements
        gaps = []
        for i in range(len(initial_elements) - 1):
            gap = initial_elements[i + 1]['start'] - initial_elements[i]['end']
            gaps.append(gap)

        # Determine typical element size from larger elements (filter out fragments)
        large_elements = [e['size'] for e in initial_elements if e['size'] >= FRAGMENT_SIZE]
        if large_elements:
            typical_element_size = sorted(large_elements)[len(large_elements) // 2]
        else:
            typical_element_size = EXPECTED_ELEMENT_SIZE

        # Determine typical gap from gaps between large elements
        large_element_gaps = []
        for i in range(len(initial_elements) - 1):
            if initial_elements[i]['size'] >= FRAGMENT_SIZE and initial_elements[i + 1]['size'] >= FRAGMENT_SIZE:
                large_element_gaps.append(gaps[i])
        if large_element_gaps:
            typical_gap = sorted(large_element_gaps)[len(large_element_gaps) // 2]
        else:
            typical_gap = EXPECTED_GAP_SIZE

        # Calculate acceptable ranges
        min_element_size = typical_element_size * (1 - SIZE_TOLERANCE)
        max_element_size = typical_element_size * (1 + SIZE_TOLERANCE)
        max_merge_gap = typical_gap * (1 + SIZE_TOLERANCE)

        # Second pass: merge fragments using sliding window validation
        merged_elements = []
        current_start = initial_elements[0]['start']
        current_end = initial_elements[0]['end']

        for i, elem in enumerate(initial_elements[1:], 1):
            gap = elem['start'] - current_end
            current_size = current_end - current_start
            next_size = elem['size']

            # Decide whether to merge based on pattern consistency
            should_merge = False

            # Rule 1: Merge very small fragments (< FRAGMENT_SIZE) if gap is reasonable
            if (current_size < FRAGMENT_SIZE or next_size < FRAGMENT_SIZE) and gap <= max_merge_gap:
                # Check if merged size would be reasonable
                merged_size = elem['end'] - current_start
                if merged_size <= max_element_size * 1.5:  # Allow some slack for merged elements
                    should_merge = True

            # Rule 2: Use sliding window to check local consistency
            if not should_merge and gap <= max_merge_gap:
                # Look at nearby elements to see if this merge would create consistency
                window_start = max(0, len(merged_elements) - ELEMENT_WINDOW_SIZE + 1)
                nearby_sizes = [e['size'] for e in merged_elements[window_start:]]

                if nearby_sizes:
                    nearby_median = sorted(nearby_sizes)[len(nearby_sizes) // 2]
                    merged_size = elem['end'] - current_start

                    # Merge if it would bring size closer to nearby elements
                    current_deviation = abs(current_size - nearby_median)
                    merged_deviation = abs(merged_size - nearby_median)
                    if merged_deviation < current_deviation and merged_size <= max_element_size * 1.5:
                        should_merge = True

            if should_merge:
                # Merge: extend current element
                current_end = elem['end']
            else:
                # Finalize current element, start new one
                merged_elements.append({
                    'start': current_start,
                    'end': current_end,
                    'size': current_end - current_start
                })
                current_start = elem['start']
                current_end = elem['end']

        # Don't forget the last element
        merged_elements.append({
            'start': current_start,
            'end': current_end,
            'size': current_end - current_start
        })

        # Calculate element size statistics
        element_sizes = [e['size'] for e in merged_elements]
        median_size = sorted(element_sizes)[len(element_sizes) // 2] if element_sizes else 0
        mean_size = sum(element_sizes) / len(element_sizes) if element_sizes else 0

        # Update array with merged elements
        updated_arr = {
            'start': arr['start'],
            'end': arr['end'],
            'elements': [(e['start'], e['end']) for e in merged_elements],
            'element_count': len(merged_elements),
            'element_sizes': element_sizes,
            'median_element_size': median_size,
            'mean_element_size': mean_size,
            '_element_info': arr.get('_element_info', [])  # Preserve for compatibility
        }

        merged_arrays.append(updated_arr)

    return merged_arrays


# Fragment collection parameters
FRAGMENT_MAX_GAP = 1500  # Max gap between array and collected fragments
FRAGMENT_EXON_THRESHOLD = 0.5  # Minimum exon fraction for fragments


def extend_arrays_with_fragments(arrays, all_clusters, gene_regions, direction='trailing'):
    """
    Extend arrays to include exon-dense fragments at array boundaries.

    Small exon-dense fragments may not be included in the initial array detection.
    This function collects them if they're exon-dense and within reasonable distance.

    Args:
        arrays: List of array dicts
        all_clusters: List of cluster info dicts
        gene_regions: Gene annotation regions
        direction: 'leading' (before array) or 'trailing' (after array)

    For trailing: stops when hitting a pLAM-like cluster (intron-rich).
    For leading: stops when hitting a non-exon-dense cluster (transition).
    """
    extended_arrays = []

    for arr in arrays:
        if direction == 'trailing':
            array_boundary = arr['end']
            # Find clusters after array end
            candidates = [(c['start'], c['end'], c['size'], c['fractions'])
                          for c in all_clusters
                          if c['start'] > array_boundary]
            candidates = sorted(candidates, key=lambda x: x[0])
        else:  # leading
            array_boundary = arr['start']
            # Find clusters before array start
            candidates = [(c['start'], c['end'], c['size'], c['fractions'])
                          for c in all_clusters
                          if c['end'] < array_boundary]
            # Sort in reverse order (closest to array first)
            candidates = sorted(candidates, key=lambda x: x[1], reverse=True)

        # Collect exon-dense fragments
        collected = []
        current_boundary = array_boundary

        for start, end, size, fracs in candidates:
            # Calculate gap
            if direction == 'trailing':
                gap = start - current_boundary
            else:
                gap = current_boundary - end

            # Stop if gap too large
            if gap > FRAGMENT_MAX_GAP:
                break

            # For trailing: stop if this is a pLAM-like cluster (intron-rich)
            if direction == 'trailing':
                if fracs['exon'] < PLAM_MAX_EXON_FRACTION and size >= PLAM_MIN_SIZE:
                    break

            # Include if exon-dense (regardless of size)
            if fracs['exon'] >= FRAGMENT_EXON_THRESHOLD:
                collected.append({'start': start, 'end': end, 'size': size, 'fractions': fracs})
                current_boundary = end if direction == 'trailing' else start
            else:
                # Not exon-dense - stop collecting
                break

        # Update array if we collected fragments
        if collected:
            arr = arr.copy()

            if direction == 'trailing':
                new_end = collected[-1]['end']
                arr['end'] = new_end

                # Append collected fragments
                existing_elements = arr.get('_element_info', [])
                arr['_element_info'] = existing_elements + collected

                existing_elem_coords = list(arr.get('elements', []))
                new_elem_coords = [(c['start'], c['end']) for c in collected]
                arr['elements'] = existing_elem_coords + new_elem_coords

                existing_sizes = arr.get('element_sizes', [])
                new_sizes = [c['size'] for c in collected if c['size'] >= 500]
                arr['element_sizes'] = existing_sizes + new_sizes
            else:  # leading
                # Reverse collected to get correct order (earliest first)
                collected = list(reversed(collected))
                new_start = collected[0]['start']
                arr['start'] = new_start

                # Prepend collected fragments
                existing_elements = arr.get('_element_info', [])
                arr['_element_info'] = collected + existing_elements

                existing_elem_coords = list(arr.get('elements', []))
                new_elem_coords = [(c['start'], c['end']) for c in collected]
                arr['elements'] = new_elem_coords + existing_elem_coords

                existing_sizes = arr.get('element_sizes', [])
                new_sizes = [c['size'] for c in collected if c['size'] >= 500]
                arr['element_sizes'] = new_sizes + existing_sizes

            # Element count increases by number of substantial fragments collected
            arr['element_count'] = arr.get('element_count', 0) + len(new_sizes)

            # Recalculate stats if we have sizes
            all_sizes = arr['element_sizes']
            if all_sizes:
                arr['median_element_size'] = sorted(all_sizes)[len(all_sizes) // 2]
                arr['mean_element_size'] = sum(all_sizes) / len(all_sizes)

        extended_arrays.append(arr)

    return extended_arrays


PLAM_SEARCH_DIST = 3000  # Search up to 3kb after last BSR cluster for pLAM


def find_plam_after_array(array_end, clusters, gene_regions):
    """
    Detect terminal β-satellite element after an array (pLAM proxy).

    IMPORTANT: This detects the structural presence of terminal β-satellite,
    NOT the actual pLAM polyadenylation signal sequence. The functional pLAM
    sequence is ATTAAA; some haplotypes (especially chr10) have terminal
    β-satellite but non-functional pLAM (ATCAAA). C+ classification indicates
    β-satellite terminal detected, not confirmed functional pLAM.

    Detection criteria:
    - Large (>= PLAM_MIN_SIZE = 2500bp)
    - Intron-rich (< PLAM_MAX_EXON_FRACTION = 35% exon)
    - Within PLAM_SEARCH_DIST (3kb) of array end

    Args:
        array_end: End position of the array
        clusters: List of (start, end) BSR clusters
        gene_regions: Gene annotation data for fraction calculation

    Returns:
        dict with terminal β-satellite info or None if not found
    """
    candidates = [(s, e) for s, e in clusters
                  if s >= array_end - 500 and s <= array_end + PLAM_SEARCH_DIST]

    for start, end in sorted(candidates):
        size = end - start
        if size < PLAM_MIN_SIZE:
            continue

        fractions = get_cluster_annotation_fractions(start, end, gene_regions)

        if fractions['exon'] < PLAM_MAX_EXON_FRACTION:
            return {
                'start': start,
                'end': end,
                'size': size,
                'type': 'pLAM',
                'fractions': fractions
            }

    return None


CT_TERMINAL_MIN_SIZE = 1000  # Minimum size for a contiguous ct_specific terminal
CT_TERMINAL_MERGE_GAP = 500  # Max gap to merge ct_specific fragments within terminal region
CT_TERMINAL_MIN_DENSITY = 0.5  # At least 50% of the terminal region must be ct_specific


def find_ct_specific_terminal(array_end, region_regions):
    """
    Look for CONTIGUOUS ct_specific terminal near array end.

    The ct_specific terminal region is typically 30-40kb. We find the START
    within the search distance, then follow the full extent of the terminal.

    Scattered small ct_specific fragments don't count - we need a substantial
    contiguous block with high ct_specific density.

    Returns:
        dict with terminal info or None if not found
    """
    # First, check if ct_specific starts within search distance of array end
    ct_start_candidates = [(rs, re) for rs, re, ra in region_regions
                           if ra == 'ct_specific'
                           and rs >= array_end
                           and rs <= array_end + CT_TERMINAL_SEARCH_DIST]

    if not ct_start_candidates:
        return None

    # Find the first ct_specific region
    first_ct = min(ct_start_candidates, key=lambda x: x[0])

    # Now collect ALL ct_specific regions starting from this point
    # (the terminal can extend 30-40kb)
    all_ct = [(rs, re) for rs, re, ra in region_regions
              if ra == 'ct_specific' and rs >= first_ct[0]]

    # Merge ct_specific regions with gaps up to 500bp (terminal region)
    all_ct = sorted(all_ct)
    merged = []
    ct_bp_in_merged = []  # Track actual ct_specific bp in each merged region
    current_ct_bp = 0

    for rs, re in all_ct:
        if merged and rs - merged[-1][1] <= CT_TERMINAL_MERGE_GAP:
            merged[-1] = (merged[-1][0], re)
            current_ct_bp += re - rs
        else:
            if merged:
                ct_bp_in_merged.append(current_ct_bp)
            merged.append((rs, re))
            current_ct_bp = re - rs

    if merged:
        ct_bp_in_merged.append(current_ct_bp)

    # The terminal should be the first merged block (starting near array end)
    if not merged:
        return None

    terminal = merged[0]
    terminal_span = terminal[1] - terminal[0]
    terminal_ct_bp = ct_bp_in_merged[0]
    terminal_density = terminal_ct_bp / terminal_span if terminal_span > 0 else 0

    # Only return if the block is large enough AND has sufficient ct_specific density
    if terminal_span >= CT_TERMINAL_MIN_SIZE and terminal_density >= CT_TERMINAL_MIN_DENSITY:
        return {
            'start': terminal[0],
            'end': terminal[1],
            'size': terminal_span,
            'ct_bp': terminal_ct_bp,
            'density': round(terminal_density, 2),
            'type': 'ct_specific'
        }

    return None


def classify_arrays(arrays, clusters, gene_regions, region_regions):
    """
    Classify arrays as canonical or degraded based on terminal presence.

    Canonical = has pLAM or ct_specific terminal
    Degraded = no terminal (or too few elements)

    Returns:
        Updated arrays list with classification info
    """
    for i, arr in enumerate(arrays):
        arr['index'] = i + 1

        # Find last exon-dense cluster for THIS array specifically
        # This handles orphan fragments beyond array end but before next array/pLAM
        # For multi-array cases, we only look at clusters within this array's region
        next_array_start = arrays[i + 1]['start'] if i + 1 < len(arrays) else float('inf')

        arr_last_exon_dense_end = arr['end']
        for s, e in clusters:
            # Only consider clusters between this array's start and next array's start
            if s >= arr['start'] and s < next_array_start and s <= arr['end'] + ARRAY_GAP_THRESHOLD:
                fracs = get_cluster_annotation_fractions(s, e, gene_regions)
                if fracs['exon'] >= EXON_DENSE_THRESHOLD:
                    arr_last_exon_dense_end = max(arr_last_exon_dense_end, e)

        # Search for pLAM from last exon-dense cluster of THIS array
        search_pos = arr_last_exon_dense_end
        plam = find_plam_after_array(search_pos, clusters, gene_regions)
        arr['plam'] = plam

        ct_terminal = find_ct_specific_terminal(arr['end'], region_regions)
        arr['ct_terminal'] = ct_terminal

        arr['has_terminal'] = plam is not None or ct_terminal is not None

        # Classification logic:
        # - pLAM terminal: only need 1+ elements (β-satellite suggests functional)
        # - ct_specific terminal only: need 3+ elements
        # - No terminal: need 3+ elements for degraded, else fragment
        if plam is not None and arr['element_count'] >= 1:
            arr['type'] = 'canonical'
        elif ct_terminal is not None and arr['element_count'] >= MIN_ARRAY_ELEMENTS:
            arr['type'] = 'canonical'
        elif arr['element_count'] >= MIN_ARRAY_ELEMENTS:
            arr['type'] = 'degraded'
        else:
            arr['type'] = 'fragment'

    return arrays


def build_clusters_info_from_merged(arrays, original_cluster_info, gene_regions):
    """
    Build clusters_info from MERGED elements in arrays.

    This ensures the visualization shows the merged elements (consistent sizes)
    rather than the original pre-merge fragments.

    Returns:
        List of dicts with cluster info for visualization
    """
    TRANSITION_EXON_THRESHOLD = 0.3
    classifications = []

    # Collect all positions covered by arrays (including merged elements and pLAMs)
    array_coverage = set()

    # 1. Add merged elements from each array
    for arr in arrays:
        arr_type = arr.get('type', 'fragment')
        arr_idx = arr.get('index', 0)

        # Determine terminal type for this array (for coloring canonical arrays)
        if arr.get('plam'):
            terminal_type = 'plam'
        elif arr.get('ct_terminal'):
            terminal_type = 'ct_specific'
        else:
            terminal_type = None

        # Add merged elements
        for elem_idx, (es, ee) in enumerate(arr['elements']):
            fractions = get_cluster_annotation_fractions(es, ee, gene_regions)
            classifications.append({
                'start': es,
                'end': ee,
                'size': ee - es,
                'exon_fraction': fractions['exon'],
                'type': 'array_element',
                'array_idx': arr_idx,
                'element_num': elem_idx + 1,
                'array_type': arr_type,
                'terminal_type': terminal_type
            })
            # Mark this range as covered
            array_coverage.add((es, ee))

        # Add pLAM if present
        if arr.get('plam'):
            ps, pe = arr['plam']['start'], arr['plam']['end']
            fractions = get_cluster_annotation_fractions(ps, pe, gene_regions)
            classifications.append({
                'start': ps,
                'end': pe,
                'size': pe - ps,
                'exon_fraction': fractions['exon'],
                'type': 'plam',
                'array_idx': arr_idx,
                'element_num': None,
                'array_type': arr_type
            })
            array_coverage.add((ps, pe))

    # 2. Add clusters that are NOT part of any array (transitions or orphan elements)
    for c in original_cluster_info:
        cs, ce = c['start'], c['end']

        # Check if this cluster overlaps with any array coverage
        is_covered = False
        for (es, ee) in array_coverage:
            if cs >= es and ce <= ee:
                is_covered = True
                break
            # Also check partial overlap
            if cs < ee and ce > es:
                is_covered = True
                break

        if not is_covered:
            size = ce - cs
            exon_frac = c['fractions']['exon']
            intron_frac = c['fractions']['intron']
            intergenic_frac = c['fractions']['intergenic']

            if exon_frac >= TRANSITION_EXON_THRESHOLD:
                # High exon = orphan array element (fragment)
                cls_type = 'array_element'
                arr_type = 'fragment'
            elif (size >= SOLITARY_PLAM_MIN_SIZE and
                  exon_frac < PLAM_MAX_EXON_FRACTION and
                  intron_frac > intergenic_frac):
                # Large, low exon, intron-dominant = solitary pLAM
                cls_type = 'solitary_plam'
                arr_type = None
            else:
                # Small or intergenic-dominant = true transition
                cls_type = 'transition'
                arr_type = None

            classifications.append({
                'start': cs,
                'end': ce,
                'size': size,
                'exon_fraction': exon_frac,
                'intron_fraction': intron_frac,
                'intergenic_fraction': intergenic_frac,
                'type': cls_type,
                'array_idx': None,
                'element_num': None,
                'array_type': arr_type
            })

    # Sort by position
    classifications.sort(key=lambda x: x['start'])

    return classifications


def classify_all_clusters(clusters, arrays, gene_regions):
    """
    Classify every cluster as array_element, plam, or transition.

    Classification logic:
    - If cluster is WITHIN an array's boundaries → array_element
    - If cluster is a pLAM (matched by coordinates) → plam
    - If cluster has HIGH exon content (>= 0.3) but not in array → array_element (orphan)
    - If cluster has LOW exon content (< 0.3) → transition

    This ensures that exon-dense clusters are never classified as transition,
    since true transition regions are intergenic/intron-only (exon ~= 0).

    Returns:
        List of dicts with cluster classification info
    """
    # Threshold for considering something a true transition (intergenic/intron)
    TRANSITION_EXON_THRESHOLD = 0.3

    classifications = []

    for start, end in clusters:
        size = end - start
        fractions = get_cluster_annotation_fractions(start, end, gene_regions)

        classification = {
            'start': start,
            'end': end,
            'size': size,
            'exon_fraction': fractions['exon'],
            'type': None,  # Will be set below
            'array_idx': None,
            'element_num': None,
            'array_type': None
        }

        # First check if this is a pLAM (after any array)
        is_plam = False
        for arr in arrays:
            if arr['plam']:
                ps, pe = arr['plam']['start'], arr['plam']['end']
                if abs(start - ps) < 500 and abs(end - pe) < 500:
                    classification['type'] = 'plam'
                    classification['array_idx'] = arr['index']
                    classification['array_type'] = arr['type']
                    is_plam = True
                    break

        if is_plam:
            classifications.append(classification)
            continue

        # Check if cluster is WITHIN any array's boundaries
        in_array = False
        for arr in arrays:
            arr_start = arr['start']
            arr_end = arr['end']

            # Check if cluster is within this array's boundaries
            if start >= arr_start - 100 and end <= arr_end + 100:
                in_array = True
                classification['array_idx'] = arr['index']
                classification['array_type'] = arr['type']
                classification['type'] = 'array_element'

                # Find which element number (approximate by position)
                for elem_idx, (es, ee) in enumerate(arr['elements']):
                    if start >= es - 100 and end <= ee + 100:
                        classification['element_num'] = elem_idx + 1
                        break
                break

        # If not in any array, classify by exon content
        if not in_array:
            if fractions['exon'] >= TRANSITION_EXON_THRESHOLD:
                # High exon content = array element (just not in a formal array)
                classification['type'] = 'array_element'
                classification['array_type'] = 'fragment'  # Mark as orphan/fragment
            else:
                # Low exon content = true transition
                classification['type'] = 'transition'

        classifications.append(classification)

    return classifications


# Keep old function name as alias for compatibility
def get_bsr_clusters(bsat_regions, merge_threshold=BSR_MERGE_THRESHOLD):
    """Compatibility alias for get_bsr_units."""
    return get_bsr_units(bsat_regions, merge_threshold)


# Multi-tract detection constants
TRACT_GAP_THRESHOLD = 10000  # Gap > 10kb suggests separate tracts (potential segmental duplication)


def detect_multiple_tracts(bsat_regions, gap_threshold=TRACT_GAP_THRESHOLD):
    """
    Detect multiple separate D4Z4 tracts (potential segmental duplications).

    Some haplotypes have two or more distinct D4Z4 arrays separated by large gaps,
    suggesting segmental duplication or complex rearrangement.

    Args:
        bsat_regions: List of (start, end) tuples for bsat_specific regions
        gap_threshold: Minimum gap size to consider as tract boundary (default 10kb)

    Returns:
        List of dicts with tract info: {start, end, length, bsat_count, est_repeats}
    """
    if not bsat_regions:
        return []

    sorted_regions = sorted(bsat_regions, key=lambda x: x[0])

    tracts = []
    current_tract_start = sorted_regions[0][0]
    current_tract_regions = [sorted_regions[0]]

    for i, (start, end) in enumerate(sorted_regions[1:], 1):
        gap = start - sorted_regions[i-1][1]

        if gap > gap_threshold:
            # Large gap - end current tract and start new one
            tract_end = sorted_regions[i-1][1]
            tract_length = tract_end - current_tract_start
            tracts.append({
                'start': current_tract_start,
                'end': tract_end,
                'length': tract_length,
                'bsat_count': len(current_tract_regions),
                'est_repeats': round(tract_length / REPEAT_UNIT_SIZE),
                'gap_after': gap
            })
            current_tract_start = start
            current_tract_regions = [(start, end)]
        else:
            current_tract_regions.append((start, end))

    # Add final tract
    tract_end = sorted_regions[-1][1]
    tract_length = tract_end - current_tract_start
    tracts.append({
        'start': current_tract_start,
        'end': tract_end,
        'length': tract_length,
        'bsat_count': len(current_tract_regions),
        'est_repeats': round(tract_length / REPEAT_UNIT_SIZE),
        'gap_after': None
    })

    return tracts


def analyze_tract(tract, bsat_regions, gene_data, repeat_data, region_data, hap_name,
                  global_main_array_start=None):
    """
    Analyze a single tract to determine if it's a functional array or degenerate tract.

    Functional Array criteria:
    - >= MIN_ARRAY_CLUSTERS (3) exon-dense BSR clusters
    - Terminal element: pLAM-like (large terminal bsat) OR ct_region

    Degenerate Tract:
    - < MIN_ARRAY_CLUSTERS exon-dense clusters, OR
    - No proper terminal element (neither pLAM nor ct_region)

    Transition Zone:
    - Clusters before the first exon-dense cluster in the first tract
    - These are NOT counted as part of the tract

    Args:
        tract: Dict with 'start', 'end', 'length', etc.
        bsat_regions: List of (start, end) tuples for ALL bsat_specific regions
        gene_data: List of (chrom, start, end, annotation) for gene annotations
        repeat_data: List of (chrom, start, end, annotation) for repeat annotations
        region_data: List of (chrom, start, end, annotation) for region annotations
        hap_name: Haplotype name to filter annotations
        global_main_array_start: The main array start for the entire haplotype (for transition detection)

    Returns:
        Dict with tract analysis results:
        {
            'start': tract start,
            'end': tract end,
            'array_start': start of exon-dense array portion (or None),
            'transition_end': end of transition zone (clusters before this are transition),
            'is_functional': True if meets all criteria,
            'exon_dense_clusters': count of exon-dense clusters,
            'total_clusters': total BSR clusters in tract,
            'has_plam': True if terminal pLAM detected,
            'has_ct_terminal': True if ct_specific terminal,
            'terminal_feature': 'pLAM-like', 'ct_specific', or 'none',
            'repeat_count': estimated repeats in functional portion,
            'degenerate_reason': reason if degenerate (or None)
        }
    """
    tract_start = tract['start']
    tract_end = tract['end']

    # Get bsat regions within this tract
    tract_bsat = [(s, e) for s, e in bsat_regions if s >= tract_start and e <= tract_end + 1000]
    if not tract_bsat:
        return {
            'start': tract_start,
            'end': tract_end,
            'array_start': None,
            'transition_end': None,
            'is_functional': False,
            'exon_dense_clusters': 0,
            'total_clusters': 0,
            'has_plam': False,
            'has_ct_terminal': False,
            'terminal_feature': 'none',
            'repeat_count': 0,
            'degenerate_reason': 'no_bsat_regions'
        }

    # Get BSR clusters within this tract
    clusters = get_bsr_clusters(tract_bsat)
    if not clusters:
        return {
            'start': tract_start,
            'end': tract_end,
            'array_start': None,
            'transition_end': None,
            'is_functional': False,
            'exon_dense_clusters': 0,
            'total_clusters': 0,
            'has_plam': False,
            'has_ct_terminal': False,
            'terminal_feature': 'none',
            'repeat_count': 0,
            'degenerate_reason': 'no_clusters'
        }

    # Filter annotations for this haplotype
    hap_genes = [(s, e, a) for (c, s, e, a) in gene_data if c == hap_name]
    hap_repeats = [(s, e, a) for (c, s, e, a) in repeat_data if c == hap_name]
    hap_regions = [(s, e, a) for (c, s, e, a) in region_data if c == hap_name]

    # Get exon_specific regions
    exon_regions = [(s, e) for (s, e, a) in hap_genes if a == 'exon_specific']

    # Get ct_specific regions
    ct_regions = [(s, e) for (s, e, a) in hap_regions if a == 'ct_specific']

    def calc_exon_density(window_start, window_end):
        """Calculate fraction of window covered by exons"""
        window_size = window_end - window_start
        if window_size <= 0:
            return 0
        exon_bp = sum(min(e, window_end) - max(s, window_start)
                      for s, e in exon_regions
                      if s < window_end and e > window_start)
        return exon_bp / window_size

    def get_gene_annotation(pos_start, pos_end):
        """Get predominant gene annotation for a region"""
        overlaps = [(s, e, a) for (s, e, a) in hap_genes if s < pos_end and e > pos_start]
        if not overlaps:
            return 'none'
        annot_bp = {}
        for s, e, a in overlaps:
            overlap = min(e, pos_end) - max(s, pos_start)
            annot_bp[a] = annot_bp.get(a, 0) + overlap
        return max(annot_bp, key=annot_bp.get)

    def has_ct_overlap(pos_start, pos_end):
        """Check if position overlaps with ct_specific regions"""
        for s, e in ct_regions:
            if s < pos_end and e > pos_start:
                return True
        return False

    def is_plam_like(cluster_start, cluster_end):
        """Check if a cluster is pLAM-like (large + intron-rich)"""
        size = cluster_end - cluster_start
        if size < PLAM_MIN_SIZE:
            return False
        density = calc_exon_density(cluster_start, cluster_end)
        gene_annot = get_gene_annotation(cluster_start, cluster_end)
        # pLAM-like if large and low exon density or intron/intergenic annotation
        return density < 0.5 or gene_annot in ('intron_specific', 'intergenic_specific', 'none')

    def is_exon_dense(cluster_start, cluster_end):
        """Check if a cluster is exon-dense"""
        size = cluster_end - cluster_start
        if size < MIN_CLUSTER_SIZE:
            return False
        density = calc_exon_density(cluster_start, cluster_end)
        if density < EXON_DENSITY_THRESHOLD:
            return False
        gene_annot = get_gene_annotation(cluster_start, cluster_end)
        return gene_annot == 'exon_specific'

    # Classify each cluster
    cluster_info = []
    for cluster_start, cluster_end in clusters:
        info = {
            'start': cluster_start,
            'end': cluster_end,
            'size': cluster_end - cluster_start,
            'is_plam_like': is_plam_like(cluster_start, cluster_end),
            'is_exon_dense': is_exon_dense(cluster_start, cluster_end)
        }
        cluster_info.append(info)

    # Identify pLAM-like clusters (potential segment boundaries)
    plam_indices = [i for i, c in enumerate(cluster_info) if c['is_plam_like']]

    # Segment the tract by pLAM-like boundaries
    # Each segment ends with a pLAM-like cluster
    segments = []
    seg_start_idx = 0
    for plam_idx in plam_indices:
        segment = {
            'start_idx': seg_start_idx,
            'end_idx': plam_idx,
            'clusters': cluster_info[seg_start_idx:plam_idx + 1],
            'exon_dense_count': sum(1 for c in cluster_info[seg_start_idx:plam_idx] if c['is_exon_dense']),
            'has_plam': True,
            'plam_cluster': cluster_info[plam_idx]
        }
        segments.append(segment)
        seg_start_idx = plam_idx + 1

    # Handle remaining clusters after last pLAM (if any)
    if seg_start_idx < len(cluster_info):
        remaining = cluster_info[seg_start_idx:]
        # Check for ct_specific terminal
        last_cluster = remaining[-1]
        CT_TERMINAL_DISTANCE = 5000
        has_ct = has_ct_overlap(last_cluster['start'], last_cluster['end'] + CT_TERMINAL_DISTANCE)
        segment = {
            'start_idx': seg_start_idx,
            'end_idx': len(cluster_info) - 1,
            'clusters': remaining,
            'exon_dense_count': sum(1 for c in remaining if c['is_exon_dense']),
            'has_plam': False,
            'has_ct': has_ct
        }
        segments.append(segment)

    # Find the best functional segment (most exon-dense clusters with terminal)
    best_segment = None
    best_count = 0
    is_mosaic = len(plam_indices) > 1  # Multiple pLAM-like = mosaic

    for seg in segments:
        has_terminal = seg.get('has_plam', False) or seg.get('has_ct', False)
        if seg['exon_dense_count'] >= MIN_ARRAY_CLUSTERS and has_terminal:
            if seg['exon_dense_count'] > best_count:
                best_count = seg['exon_dense_count']
                best_segment = seg

    # Determine overall tract classification
    total_exon_dense = sum(1 for c in cluster_info if c['is_exon_dense'])

    # Find the first exon-dense cluster (start of array portion)
    first_exon_dense_cluster = next((c for c in cluster_info if c['is_exon_dense']), None)

    # Determine transition zone: clusters before the first exon-dense cluster
    # Only applies if this tract contains the global main array start
    transition_end = None
    if first_exon_dense_cluster and global_main_array_start:
        # If this tract contains clusters before the global main array start,
        # those are transition clusters
        if tract_start <= global_main_array_start <= tract_end:
            transition_end = first_exon_dense_cluster['start']

    # Check terminal of entire tract
    last_cluster = clusters[-1]
    last_size = last_cluster[1] - last_cluster[0]
    last_gene_annot = get_gene_annotation(last_cluster[0], last_cluster[1])
    has_plam = last_size >= PLAM_MIN_SIZE
    if not has_plam and last_gene_annot in ('intergenic_specific', 'intron_specific', 'none'):
        if last_size > EXPECTED_BSR_SIZE * 1.3:
            has_plam = True
    CT_TERMINAL_DISTANCE = 5000
    has_ct_terminal = has_ct_overlap(last_cluster[0], last_cluster[1] + CT_TERMINAL_DISTANCE)

    if has_plam:
        terminal_feature = 'β-satellite terminal'
    elif has_ct_terminal:
        terminal_feature = 'No β-satellite terminal'
    else:
        terminal_feature = 'none'

    # Determine if functional (either as whole or has functional sub-array)
    is_functional = False
    degenerate_reason = None
    array_start = None
    repeat_count = 0

    if best_segment is not None:
        # Has a functional sub-array (mosaic case)
        is_functional = True
        # Find the first exon-dense cluster in the best segment
        first_exon_dense = next((c for c in best_segment['clusters'] if c['is_exon_dense']), None)
        if first_exon_dense:
            array_start = first_exon_dense['start']
        repeat_count = best_segment['exon_dense_count']
    elif total_exon_dense >= MIN_ARRAY_CLUSTERS and terminal_feature != 'none':
        # Whole tract is functional (non-mosaic)
        is_functional = True
        if first_exon_dense_cluster:
            array_start = first_exon_dense_cluster['start']
        # Count exon-dense clusters, excluding the pLAM terminal if present
        repeat_count = total_exon_dense
        # Check if last cluster is pLAM-like - don't count it as a repeat
        if cluster_info and cluster_info[-1]['is_plam_like']:
            repeat_count = total_exon_dense  # pLAM is already not counted as exon-dense
    else:
        # Degenerate - set array_start to first exon-dense cluster for visualization
        if first_exon_dense_cluster:
            array_start = first_exon_dense_cluster['start']
        if total_exon_dense < MIN_ARRAY_CLUSTERS:
            degenerate_reason = f'insufficient_exon_dense_clusters ({total_exon_dense}<{MIN_ARRAY_CLUSTERS})'
        elif terminal_feature == 'none':
            degenerate_reason = 'no_terminal_element (neither pLAM nor ct_region)'

    # Build sub-array info for mosaic tracts
    sub_arrays = []
    if is_mosaic:
        for i, seg in enumerate(segments):
            has_terminal = seg.get('has_plam', False) or seg.get('has_ct', False)
            seg_functional = seg['exon_dense_count'] >= MIN_ARRAY_CLUSTERS and has_terminal
            # Find first exon-dense in this segment for proper start
            seg_first_exon = next((c for c in seg['clusters'] if c['is_exon_dense']), None)
            sub_arrays.append({
                'index': i + 1,
                'exon_dense': seg['exon_dense_count'],
                'has_terminal': has_terminal,
                'is_functional': seg_functional,
                'start': seg['clusters'][0]['start'] if seg['clusters'] else None,
                'end': seg['clusters'][-1]['end'] if seg['clusters'] else None,
                'array_start': seg_first_exon['start'] if seg_first_exon else None
            })

    return {
        'start': tract_start,
        'end': tract_end,
        'array_start': array_start,
        'transition_end': transition_end,
        'is_functional': is_functional,
        'exon_dense_clusters': total_exon_dense,
        'total_clusters': len(clusters),
        'has_plam': has_plam,
        'has_ct_terminal': has_ct_terminal,
        'terminal_feature': terminal_feature,
        'repeat_count': repeat_count,
        'degenerate_reason': degenerate_reason,
        'is_mosaic': is_mosaic,
        'num_plam_like': len(plam_indices),
        'sub_arrays': sub_arrays if is_mosaic else None,
        'best_sub_array': best_segment['exon_dense_count'] if best_segment else None
    }


def filter_plam_clusters(clusters, gene_data, hap_name, tar1_position=None):
    """
    Filter out terminal pLAM-like clusters from the repeat count.

    pLAM elements are identified by:
    1. Being the last cluster (closest to telomere)
    2. Being unusually large (>3500bp) OR
    3. Being in a non-exon region (intergenic/intron) AND close to TAR1

    Args:
        clusters: List of (start, end) tuples
        gene_data: Gene annotation data
        hap_name: Haplotype name
        tar1_position: Position of TAR1 if known

    Returns:
        Tuple of (filtered_clusters, plam_cluster, plam_info)
    """
    if not clusters:
        return clusters, None, None

    # Get gene annotations for this haplotype
    hap_genes = [(s, e, a) for (c, s, e, a) in gene_data if c == hap_name]

    def get_gene_annotation(pos_start, pos_end):
        """Get predominant gene annotation for a region"""
        overlaps = [(s, e, a) for (s, e, a) in hap_genes if s < pos_end and e > pos_start]
        if not overlaps:
            return 'none'
        annot_bp = {}
        for s, e, a in overlaps:
            overlap = min(e, pos_end) - max(s, pos_start)
            annot_bp[a] = annot_bp.get(a, 0) + overlap
        return max(annot_bp, key=annot_bp.get)

    # Check the last cluster
    last_cluster = clusters[-1]
    last_size = last_cluster[1] - last_cluster[0]
    last_gene_annot = get_gene_annotation(last_cluster[0], last_cluster[1])

    # Determine if last cluster is pLAM
    is_plam = False
    plam_reason = None

    # Criterion 1: Unusually large (>3500bp)
    if last_size >= PLAM_MIN_SIZE:
        is_plam = True
        plam_reason = f"large_terminal ({last_size}bp)"

    # Criterion 2: In non-exon region and close to TAR1
    elif last_gene_annot in ('intergenic_specific', 'intron_specific', 'none'):
        if tar1_position and (tar1_position - last_cluster[1]) < 5000:
            is_plam = True
            plam_reason = f"non_exon_terminal ({last_gene_annot})"
        elif last_size > EXPECTED_BSR_SIZE * 1.3:  # Also flag if moderately large + non-exon
            is_plam = True
            plam_reason = f"large_non_exon ({last_size}bp, {last_gene_annot})"

    if is_plam:
        return clusters[:-1], last_cluster, plam_reason
    else:
        return clusters, None, None


def detect_novel_sequence(region_data, gene_data, hap_name, main_array_start, main_array_end):
    """
    Detect novel/variant sequence within the D4Z4 array.

    Novel sequence is identified by KaryoScope as regions that don't match
    any known reference patterns. These could be interesting variants.

    Args:
        region_data: List of (chrom, start, end, annotation) for region annotations
        gene_data: List of (chrom, start, end, annotation) for gene annotations
        hap_name: Haplotype name
        main_array_start: Start of main D4Z4 array
        main_array_end: End of main D4Z4 array

    Returns:
        Tuple of (novel_bp, novel_regions) where novel_regions is a list of (start, end, size)
    """
    if main_array_start is None or main_array_end is None:
        return 0, []

    # Find "novel" annotations within the main array
    novel_regions = []

    # Check region annotations
    for (c, s, e, a) in region_data:
        if c == hap_name and a == 'novel':
            # Check if within main array
            if s >= main_array_start and e <= main_array_end:
                novel_regions.append((s, e, e - s))

    # Check gene annotations for novel
    for (c, s, e, a) in gene_data:
        if c == hap_name and a == 'novel':
            if s >= main_array_start and e <= main_array_end:
                # Avoid duplicates
                if not any(nr[0] == s and nr[1] == e for nr in novel_regions):
                    novel_regions.append((s, e, e - s))

    # Merge overlapping novel regions
    if novel_regions:
        novel_regions.sort()
        merged = []
        curr_start, curr_end, _ = novel_regions[0]
        for s, e, _ in novel_regions[1:]:
            if s <= curr_end:
                curr_end = max(curr_end, e)
            else:
                merged.append((curr_start, curr_end, curr_end - curr_start))
                curr_start, curr_end = s, e
        merged.append((curr_start, curr_end, curr_end - curr_start))
        novel_regions = merged

    novel_bp = sum(size for _, _, size in novel_regions)
    return novel_bp, novel_regions


def find_main_array_boundary(bsat_regions, gene_data, repeat_data, region_data, hap_name):
    """
    Find the start of the main contiguous D4Z4 array using gene and repeat annotations.

    The main array is identified by ALL of (v3 - stricter criteria):
    1. High exon density (>80%) - indicates DUX4 genes present
    2. Large cluster size (>1500bp) - full BSR region
    3. Gene annotation is exon_specific (not intergenic_specific)
    4. Absence of centromeric_multigroup1 in cluster region

    Args:
        bsat_regions: List of (start, end) tuples for bsat_specific
        gene_data: List of (chrom, start, end, annotation) for gene annotations
        repeat_data: List of (chrom, start, end, annotation) for repeat annotations
        region_data: List of (chrom, start, end, annotation) for region annotations
        hap_name: Haplotype name to filter annotations

    Returns:
        Tuple of (main_array_start, boundary_type, boundary_position)
        boundary_type: 'LINE', 'LTR', 'SINE', 'cluster_size', 'exon_specific', or 'none'
    """
    if not bsat_regions:
        return None, 'none', None

    # Get BSR clusters
    clusters = get_bsr_clusters(bsat_regions)
    if not clusters:
        return None, 'none', None

    bsat_start = min(s for s, e in bsat_regions)
    bsat_end = max(e for s, e in bsat_regions)

    # Filter annotations for this haplotype
    hap_genes = [(s, e, a) for (c, s, e, a) in gene_data if c == hap_name]
    hap_repeats = [(s, e, a) for (c, s, e, a) in repeat_data if c == hap_name]
    hap_regions = [(s, e, a) for (c, s, e, a) in region_data if c == hap_name]

    # Find LINE/LTR/SINE elements in the bsat region
    interspersed_repeats = [(s, e, a) for (s, e, a) in hap_repeats
                            if a in ('LINE_specific', 'LTR_specific', 'SINE_specific')
                            and s >= bsat_start - 50000  # Look 50kb before first bsat
                            and s <= bsat_end]

    # Get exon_specific regions for density calculation
    exon_regions = [(s, e) for (s, e, a) in hap_genes if a == 'exon_specific']

    # Get centromeric regions to exclude from main array
    centromeric_regions = [(s, e) for (s, e, a) in hap_regions
                           if 'centromeric' in a.lower()]

    def calc_exon_density(window_start, window_end):
        """Calculate fraction of window covered by exons"""
        window_size = window_end - window_start
        if window_size <= 0:
            return 0
        exon_bp = sum(min(e, window_end) - max(s, window_start)
                      for s, e in exon_regions
                      if s < window_end and e > window_start)
        return exon_bp / window_size

    def get_gene_annotation_at(pos_start, pos_end):
        """Get the predominant gene annotation type for a position range"""
        overlaps = [(s, e, a) for (s, e, a) in hap_genes
                    if s < pos_end and e > pos_start]
        if not overlaps:
            return 'none'
        # Find annotation with most overlap
        annot_bp = {}
        for s, e, a in overlaps:
            overlap = min(e, pos_end) - max(s, pos_start)
            annot_bp[a] = annot_bp.get(a, 0) + overlap
        return max(annot_bp, key=annot_bp.get)

    def has_centromeric_overlap(pos_start, pos_end):
        """Check if position overlaps with centromeric regions"""
        for s, e in centromeric_regions:
            if s < pos_end and e > pos_start:
                return True
        return False

    def has_repeat_interspersed(pos_start, pos_end):
        """Check if position has LINE/LTR/SINE interspersed"""
        for s, e, a in interspersed_repeats:
            if s < pos_end and e > pos_start:
                return True
        return False

    # Find where main array starts using stricter criteria
    main_array_start = None
    boundary_type = 'none'
    boundary_position = None

    # Check each BSR cluster against all criteria
    for i, (cluster_start, cluster_end) in enumerate(clusters):
        cluster_size = cluster_end - cluster_start

        # Criterion 1: Cluster size must be >= MIN_CLUSTER_SIZE
        if cluster_size < MIN_CLUSTER_SIZE:
            continue

        # Criterion 2: Calculate exon density around this cluster
        window_start = cluster_start
        window_end = min(cluster_start + EXON_WINDOW_SIZE, bsat_end)
        density = calc_exon_density(window_start, window_end)

        if density < EXON_DENSITY_THRESHOLD:
            continue

        # Criterion 3: Gene annotation should be exon_specific
        gene_annot = get_gene_annotation_at(cluster_start, cluster_end)
        if gene_annot != 'exon_specific':
            continue

        # Criterion 4: No centromeric regions overlapping
        if has_centromeric_overlap(cluster_start, cluster_end):
            continue

        # All criteria passed - this is the start of main array
        # Determine boundary type based on what precedes this cluster
        repeats_before = [(s, e, a) for (s, e, a) in interspersed_repeats
                         if e <= cluster_start and e > cluster_start - 10000]

        if repeats_before:
            # Use the last interspersed repeat as boundary
            last_repeat = max(repeats_before, key=lambda x: x[1])
            main_array_start = cluster_start
            boundary_type = last_repeat[2].replace('_specific', '')
            boundary_position = last_repeat[1]
        else:
            # No repeat boundary, classify by what criteria first passed
            main_array_start = cluster_start
            # Check what changed from previous cluster
            if i > 0:
                prev_cluster = clusters[i-1]
                prev_size = prev_cluster[1] - prev_cluster[0]
                prev_density = calc_exon_density(prev_cluster[0],
                                                 min(prev_cluster[0] + EXON_WINDOW_SIZE, bsat_end))
                prev_annot = get_gene_annotation_at(prev_cluster[0], prev_cluster[1])

                if prev_size < MIN_CLUSTER_SIZE:
                    boundary_type = 'cluster_size'
                elif prev_density < EXON_DENSITY_THRESHOLD:
                    boundary_type = 'exon_density'
                elif prev_annot != 'exon_specific':
                    boundary_type = 'exon_specific'
                else:
                    boundary_type = 'criteria_change'
            else:
                boundary_type = 'first_cluster'
            boundary_position = cluster_start
        break

    # Fallback: if no cluster passes all criteria, try relaxed approach
    if main_array_start is None:
        # Look for the last major LINE/LTR before bsat
        if interspersed_repeats:
            # Find the last LINE/LTR that's followed by significant bsat
            for s, e, a in sorted(interspersed_repeats, key=lambda x: x[0], reverse=True):
                # Check if there's substantial bsat after this repeat
                bsat_after = sum(be - bs for bs, be in bsat_regions if bs >= e)
                if bsat_after > 50000:  # At least 50kb of bsat after
                    main_array_start = e
                    boundary_type = a.replace('_specific', '') + '_fallback'
                    boundary_position = e
                    break

        # If still no boundary found, use the first large cluster
        if main_array_start is None:
            for cluster_start, cluster_end in clusters:
                if cluster_end - cluster_start >= MIN_CLUSTER_SIZE:
                    main_array_start = cluster_start
                    boundary_type = 'size_fallback'
                    break

        # Last resort: use first cluster
        if main_array_start is None:
            main_array_start = clusters[0][0]
            boundary_type = 'none'

    return main_array_start, boundary_type, boundary_position

def count_bsr_clusters(bsat_regions, merge_threshold=BSR_MERGE_THRESHOLD):
    """
    Count repeat units by identifying distinct BSR (beta satellite) clusters.

    Each D4Z4 repeat unit contains one BSR region. The KaryoScope bsat_specific
    annotations may be fragmented, so we merge nearby blocks (within merge_threshold)
    and count the resulting clusters as repeat units.

    Args:
        bsat_regions: List of (start, end) tuples for bsat_specific regions
        merge_threshold: Maximum gap between blocks to merge (default 500bp)

    Returns:
        Integer count of repeat units
    """
    if not bsat_regions:
        return 0

    # Sort by start position
    sorted_regions = sorted(bsat_regions, key=lambda x: x[0])

    cluster_count = 0
    cluster_end = -1

    for start, end in sorted_regions:
        if cluster_end == -1 or start - cluster_end > merge_threshold:
            # New cluster (new repeat unit)
            cluster_count += 1
        cluster_end = end

    return cluster_count


# =============================================================================
# v4 Haplotype Analysis (Array-First approach)
# =============================================================================

def analyze_haplotype_v4(hap_name, region_data, subtelomeric_data, gene_data=None, repeat_data=None):
    """
    Analyze a single haplotype using the v4 Array-First approach.

    This simplified analysis:
    1. Identifies BSR clusters using conservative merging
    2. Identifies arrays as contiguous exon-dense stretches
    3. Classifies arrays as canonical/degraded based on terminals
    4. Classifies each cluster as array_element, plam, or transition

    Returns dict with:
    - Core metrics: bsat_start, bsat_end, array_length, repeat counts
    - Array info: arrays_info JSON with per-array details
    - Cluster info: clusters_info JSON with per-cluster classification
    """
    # Filter data for this haplotype
    hap_regions = [(s, e, a) for (c, s, e, a) in region_data if c == hap_name]
    hap_subtelo = [(s, e, a) for (c, s, e, a) in subtelomeric_data if c == hap_name]
    hap_genes = [(s, e, a) for (c, s, e, a) in gene_data if c == hap_name] if gene_data else []

    if not hap_regions:
        return None

    # Find bsat_specific regions
    bsat_regions = [(s, e) for (s, e, a) in hap_regions if a == 'bsat_specific']

    if not bsat_regions:
        return {
            'haplotype': hap_name,
            'bsat_start': None,
            'bsat_end': None,
            'array_length': 0,
            'repeat_count_total': 0,
            'repeat_count_main': 0,
            'num_arrays': 0,
            'num_canonical': 0,
            'num_degraded': 0,
            'terminal_feature': 'no_bsat',
            'is_truncated': False,
            'arrays_info': None,
            'clusters_info': None,
            'notes': 'No bsat_specific regions found'
        }

    bsat_start = min(s for s, e in bsat_regions)
    bsat_end = max(e for s, e in bsat_regions)
    bsat_total_bp = sum(e - s for s, e in bsat_regions)

    # Note: is_truncated will be set after arrays are identified

    # Find TAR1 position
    tar1_regions = [(s, e) for (s, e, a) in hap_subtelo if a == 'TAR1_specific']
    tar1_after_bsat = [s for s, e in tar1_regions if s > bsat_end]
    tar1_position = min(tar1_after_bsat) if tar1_after_bsat else None

    # Calculate array length
    if tar1_position:
        array_length = tar1_position - bsat_start
    else:
        array_length = bsat_end - bsat_start

    # v4 Pipeline: Array-First approach
    clusters = get_bsr_units(bsat_regions)
    arrays, cluster_info = identify_arrays(clusters, hap_genes)

    # Merge adjacent arrays that were fragmented by small non-exon-dense clusters
    # or separated by large novel sequence insertions
    arrays = merge_adjacent_arrays(arrays, cluster_info, hap_genes)

    # Extend arrays to include leading/trailing exon-dense fragments
    arrays = extend_arrays_with_fragments(arrays, cluster_info, hap_genes, direction='leading')
    arrays = extend_arrays_with_fragments(arrays, cluster_info, hap_genes, direction='trailing')

    # Within-array merging: create consistent element sizes (includes collected fragments)
    arrays = merge_elements_within_arrays(arrays, cluster_info)

    arrays = classify_arrays(arrays, clusters, hap_genes, hap_regions)

    # Check if array is truncated (first array starts near position 0)
    # Use the first identified array's start, not scattered bsat fragments
    first_array_start = None
    if arrays:
        first_array_start = min(a['start'] for a in arrays)
    is_truncated = first_array_start is not None and first_array_start < TRUNCATION_THRESHOLD

    # Build clusters_info from MERGED elements (not original clusters)
    # This ensures visualization shows merged elements, not pre-merge fragments
    cluster_classifications = build_clusters_info_from_merged(arrays, cluster_info, hap_genes)

    # Count arrays by type (excluding fragments)
    real_arrays = [a for a in arrays if a['type'] in ('canonical', 'degraded')]
    num_canonical = sum(1 for a in real_arrays if a['type'] == 'canonical')
    num_degraded = sum(1 for a in real_arrays if a['type'] == 'degraded')

    # Calculate repeat counts
    repeat_count_total = len(clusters)
    repeat_count_main = sum(a['element_count'] for a in real_arrays)

    # Determine terminal feature from the last canonical array
    terminal_feature = 'none'
    for arr in reversed(arrays):
        if arr['type'] == 'canonical':
            if arr['plam']:
                terminal_feature = 'β-satellite terminal'
            elif arr['ct_terminal']:
                terminal_feature = 'No β-satellite terminal'
            break
    if terminal_feature == 'none' and num_degraded > 0:
        terminal_feature = 'degraded (no terminal)'

    # Build arrays_info for JSON output
    # Include canonical (1+ elem with pLAM, or 3+ with ct) and degraded (3+ no terminal)
    arrays_info = []
    for arr in arrays:
        if arr['type'] in ('canonical', 'degraded'):
            info = {
                'index': arr['index'],
                'start': arr['start'],
                'end': arr['end'],
                'element_count': arr['element_count'],
                'type': arr['type'],
                'has_plam': arr['plam'] is not None,
                'has_ct_terminal': arr['ct_terminal'] is not None,
                'plam_size': arr['plam']['size'] if arr['plam'] else None,
                'plam_end': arr['plam']['end'] if arr['plam'] else None,
                'median_element_size': arr.get('median_element_size'),
                'mean_element_size': round(arr.get('mean_element_size', 0), 0) if arr.get('mean_element_size') else None,
                'element_sizes': arr.get('element_sizes', [])
            }
            arrays_info.append(info)

    # Build clusters_info for JSON output
    clusters_info = []
    for cls in cluster_classifications:
        info = {
            'start': cls['start'],
            'end': cls['end'],
            'size': cls['size'],
            'type': cls['type'],
            'array_idx': cls['array_idx'],
            'element_num': cls['element_num'],
            'array_type': cls['array_type'],
            'terminal_type': cls.get('terminal_type'),
            'exon_fraction': round(cls['exon_fraction'], 2)
        }
        clusters_info.append(info)

    # Build notes
    notes_list = []
    if len(real_arrays) > 1:
        notes_list.append(f'MULTI-ARRAY ({len(real_arrays)} arrays)')
    if is_truncated:
        notes_list.append('Array extends beyond 250kb window')
    notes = '; '.join(notes_list) if notes_list else ''

    return {
        'haplotype': hap_name,
        'bsat_start': bsat_start,
        'bsat_end': bsat_end,
        'tar1_position': tar1_position,
        'array_length': array_length,
        'bsat_total_bp': bsat_total_bp,
        'repeat_count_total': repeat_count_total,
        'repeat_count_main': repeat_count_main,
        'num_arrays': len(real_arrays),
        'num_canonical': num_canonical,
        'num_degraded': num_degraded,
        'terminal_feature': terminal_feature,
        'is_truncated': is_truncated,
        'arrays_info': arrays_info,
        'clusters_info': clusters_info,
        'notes': notes
    }


def analyze_haplotype(hap_name, region_data, subtelomeric_data, gene_data=None, repeat_data=None, full_region_data=None):
    """
    Analyze a single haplotype to determine D4Z4 array characteristics

    Returns dict with:
    - bsat_start: First position of bsat_specific
    - bsat_end: Last position of bsat_specific
    - tar1_position: Position of TAR1 (telomeric boundary)
    - array_length: Calculated array length
    - repeat_count_total: Total BSR clusters (all bsat regions)
    - repeat_count_main: BSR clusters in main contiguous array only
    - terminal_feature: 'β-satellite terminal', 'No β-satellite terminal', or 'unknown'
    - is_truncated: Whether array extends beyond 250kb window
    """

    # Filter for this haplotype
    hap_regions = [(s, e, a) for (c, s, e, a) in region_data if c == hap_name]
    hap_subtelo = [(s, e, a) for (c, s, e, a) in subtelomeric_data if c == hap_name]

    if not hap_regions:
        return None

    # Find bsat_specific regions
    bsat_regions = [(s, e) for (s, e, a) in hap_regions if a == 'bsat_specific']

    if not bsat_regions:
        return {
            'haplotype': hap_name,
            'bsat_start': None,
            'bsat_end': None,
            'main_array_start': None,
            'main_array_boundary': 'no_bsat',
            'tar1_position': None,
            'array_length': 0,
            'main_array_length': 0,
            'repeat_count_total': 0,
            'repeat_count_main': 0,
            'repeat_count_by_length': 0,
            'terminal_feature': 'no_bsat',
            'is_truncated': False,
            'bsat_total_bp': 0,
            'num_tracts': 0,
            'num_functional_arrays': 0,
            'num_degenerate_tracts': 0,
            'tract_info': None,
            'arrays_info': None,
            'notes': 'No bsat_specific regions found'
        }

    bsat_start = min(s for s, e in bsat_regions)
    bsat_end = max(e for s, e in bsat_regions)
    bsat_total_bp = sum(e - s for s, e in bsat_regions)

    # Detect multiple tracts (potential segmental duplications)
    tracts = detect_multiple_tracts(bsat_regions)
    num_tracts = len(tracts)
    tract_info = None
    arrays_info = None  # New: detailed per-array analysis
    num_functional_arrays = 0
    num_degenerate_tracts = 0

    # Pre-calculate the global main array start for transition detection
    global_main_array_start = None
    if gene_data is not None and repeat_data is not None and full_region_data is not None:
        global_main_array_start, _, _ = find_main_array_boundary(
            bsat_regions, gene_data, repeat_data, full_region_data, hap_name
        )

    # Analyze each tract to determine if functional array or degenerate
    if gene_data is not None and repeat_data is not None and full_region_data is not None:
        analyzed_tracts = []
        for t in tracts:
            analysis = analyze_tract(t, bsat_regions, gene_data, repeat_data, full_region_data, hap_name,
                                     global_main_array_start=global_main_array_start)
            analyzed_tracts.append(analysis)
            if analysis['is_functional']:
                num_functional_arrays += 1
            else:
                num_degenerate_tracts += 1

        # Build tract_info string with functional/degenerate status
        tract_details = []
        array_details = []  # Detailed info for each array/tract
        for i, (t, a) in enumerate(zip(tracts, analyzed_tracts)):
            status = 'A' if a['is_functional'] else 'D'  # A=Array, D=Degenerate
            tract_details.append(f"{status}{i+1}:{t['start']//1000}-{t['end']//1000}kb(~{a['repeat_count']}rpt)")

            # Build detailed array info
            detail = {
                'index': i + 1,
                'start': t['start'],
                'end': t['end'],
                'array_start': a['array_start'],
                'transition_end': a.get('transition_end'),
                'is_functional': a['is_functional'],
                'terminal_feature': a['terminal_feature'],
                'repeat_count': a['repeat_count'],
                'exon_dense_clusters': a['exon_dense_clusters'],
                'has_plam': a['has_plam'],
                'degenerate_reason': a['degenerate_reason'],
                'is_mosaic': a.get('is_mosaic', False),
                'num_plam_like': a.get('num_plam_like', 0),
                'sub_arrays': a.get('sub_arrays'),
                'best_sub_array': a.get('best_sub_array')
            }
            array_details.append(detail)

        tract_info = "; ".join(tract_details)
        arrays_info = array_details
    elif num_tracts > 1:
        # Fallback if no annotation data available
        tract_details = []
        for i, t in enumerate(tracts):
            tract_details.append(f"T{i+1}:{t['start']//1000}-{t['end']//1000}kb(~{t['est_repeats']}rpt)")
        tract_info = "; ".join(tract_details)

    # Check if array is truncated (first array starts near position 0)
    # Use the first identified array's start, not scattered bsat fragments
    first_array_start = None
    if arrays:
        first_array_start = min(a['start'] for a in arrays)
    is_truncated = first_array_start is not None and first_array_start < TRUNCATION_THRESHOLD

    # Find TAR1 positions (telomeric boundary)
    tar1_regions = [(s, e) for (s, e, a) in hap_subtelo if a == 'TAR1_specific']

    # Find the TAR1 that's after the bsat (telomeric side)
    tar1_after_bsat = [s for s, e in tar1_regions if s > bsat_end]
    tar1_position = min(tar1_after_bsat) if tar1_after_bsat else None

    # Calculate total array length (all bsat)
    if tar1_position:
        array_length = tar1_position - bsat_start
    else:
        array_length = bsat_end - bsat_start

    # Count TOTAL repeat units (all BSR clusters)
    repeat_count_total = count_bsr_clusters(bsat_regions)

    # Find main array boundary using gene and repeat annotations
    main_array_start = None
    main_array_boundary = 'none'
    boundary_position = None
    repeat_count_main = repeat_count_total  # Default to total
    plam_cluster = None
    plam_reason = None

    if gene_data is not None and repeat_data is not None and full_region_data is not None:
        main_array_start, main_array_boundary, boundary_position = find_main_array_boundary(
            bsat_regions, gene_data, repeat_data, full_region_data, hap_name
        )

        if main_array_start is not None:
            # Get bsat clusters in main array (not just count)
            main_array_bsat = [(s, e) for s, e in bsat_regions if s >= main_array_start]
            main_array_clusters = get_bsr_clusters(main_array_bsat)

            # Filter out terminal pLAM clusters before counting
            filtered_clusters, plam_cluster, plam_reason = filter_plam_clusters(
                main_array_clusters, gene_data, hap_name, tar1_position
            )
            repeat_count_main = len(filtered_clusters)
        else:
            filtered_clusters = []
            plam_cluster = None
            plam_reason = None

    # Store pLAM detection info
    plam_detected = plam_cluster is not None
    plam_size = plam_cluster[1] - plam_cluster[0] if plam_cluster else None
    plam_info = plam_reason

    # Calculate main array length
    if main_array_start is not None:
        if tar1_position:
            main_array_length = tar1_position - main_array_start
        else:
            main_array_end = max(e for s, e in bsat_regions if s >= main_array_start)
            main_array_length = main_array_end - main_array_start
    else:
        main_array_length = array_length

    # Also calculate estimate based on main array length for comparison
    repeat_count_by_length = main_array_length / REPEAT_UNIT_SIZE

    # Detect novel/variant sequence within the main array
    novel_bp = 0
    novel_regions = []
    novel_info = None
    if gene_data is not None and full_region_data is not None and main_array_start is not None:
        main_end = tar1_position if tar1_position else bsat_end
        novel_bp, novel_regions = detect_novel_sequence(
            full_region_data, gene_data, hap_name, main_array_start, main_end
        )
        if novel_bp > 0:
            novel_info = f"{novel_bp}bp in {len(novel_regions)} region(s)"

    # Analyze terminal feature (what's between last bsat and TAR1)
    terminal_feature = 'unknown'
    terminal_gap = None
    largest_terminal_bsat = 0
    ct_bp_terminal = 0

    if tar1_position:
        terminal_gap = tar1_position - bsat_end
        # Look at annotations between last bsat and TAR1
        terminal_regions = [(s, e, a) for (s, e, a) in hap_regions
                           if s >= bsat_end and e <= tar1_position]

        # Count ct_specific in terminal region
        ct_bp_terminal = sum(e - s for s, e, a in terminal_regions if a == 'ct_specific')

        # Check for contiguous bsat at the end (potential pLAM)
        # Look at the last few bsat regions
        last_bsat_regions = sorted(bsat_regions, key=lambda x: x[0], reverse=True)[:5]
        if last_bsat_regions:
            largest_terminal_bsat = max(e - s for s, e in last_bsat_regions)

            # Terminal feature classification:
            # 1. If pLAM was detected by filter_plam_clusters() -> β-satellite terminal
            # 2. Short gap to TAR1 (<5kb) -> β-satellite terminal
            # 3. Large terminal bsat (>3kb) -> β-satellite terminal
            # 4. Long gap with ct_specific -> No β-satellite terminal
            if plam_detected:
                terminal_feature = 'β-satellite terminal'
            elif terminal_gap <= 5000:
                terminal_feature = 'β-satellite terminal'
            elif largest_terminal_bsat > 3000:
                terminal_feature = 'β-satellite terminal'
            elif terminal_gap > 5000 and ct_bp_terminal > 1000:
                terminal_feature = 'No β-satellite terminal'
            else:
                terminal_feature = 'standard'

    # Format repeat count string (using main array count)
    if is_truncated:
        repeat_count_str = f"> {repeat_count_main}"
    else:
        repeat_count_str = str(repeat_count_main)

    # FSHD relevance annotation (using main array count)
    # FSHD1: ≤10 repeats + β-satellite terminal (chr4 only)
    # FSHD2: 11-20 repeats + β-satellite terminal + SMCHD1 mutation
    # NOTE: β-satellite terminal is a PROXY for pLAM. The actual pLAM sequence
    # (ATTAAA = functional, ATCAAA = non-functional) cannot be determined from
    # k-mer annotations. Clinical FSHD diagnosis requires sequence confirmation.
    fshd_note = ''
    if not is_truncated and repeat_count_main <= 10 and terminal_feature == 'β-satellite terminal':
        fshd_note = 'FSHD1-like (short array + β-satellite)'
    elif not is_truncated and 11 <= repeat_count_main <= 20 and terminal_feature == 'β-satellite terminal':
        fshd_note = 'FSHD2-range (requires SMCHD1 check)'

    # Build notes with multiple annotations
    notes_list = []
    if num_tracts > 1:
        notes_list.append(f'MULTI-TRACT ({num_tracts} separate arrays - potential segmental duplication)')
    if is_truncated:
        notes_list.append('Array extends beyond 250kb window')
    if fshd_note:
        notes_list.append(fshd_note)
    if novel_bp > 0:
        notes_list.append(f'Novel sequence: {novel_info}')
    notes = '; '.join(notes_list) if notes_list else ''

    return {
        'haplotype': hap_name,
        'bsat_start': bsat_start,
        'bsat_end': bsat_end,
        'main_array_start': main_array_start,
        'main_array_boundary': main_array_boundary,
        'tar1_position': tar1_position,
        'array_length': array_length,
        'main_array_length': main_array_length,
        'bsat_total_bp': bsat_total_bp,
        'repeat_count_total': repeat_count_total,
        'repeat_count_main': repeat_count_main,
        'repeat_count_by_length': repeat_count_by_length,
        'repeat_count_str': repeat_count_str,
        'terminal_feature': terminal_feature,
        'terminal_gap': terminal_gap,
        'largest_terminal_bsat': largest_terminal_bsat,
        'ct_bp_terminal': ct_bp_terminal,
        'plam_detected': plam_detected,
        'plam_size': plam_size,
        'plam_info': plam_info,
        'novel_bp': novel_bp,
        'novel_info': novel_info,
        'is_truncated': is_truncated,
        'num_tracts': num_tracts,
        'num_functional_arrays': num_functional_arrays,
        'num_degenerate_tracts': num_degenerate_tracts,
        'tract_info': tract_info,
        'arrays_info': arrays_info,
        'notes': notes
    }


def generate_array_summary(df, output_file, logger):
    """
    Generate summary table from analysis results.

    Creates a per-haplotype summary with:
    - Number of arrays
    - Canonical vs degraded counts
    - Array types in order (C+=canonical with β-satellite, C-=canonical without, D=degraded)
    - Element counts for each array
    """
    rows = []

    for _, row in df.iterrows():
        chromosome = row['chromosome']
        sample = row['sample']
        hap_num = row['hap_num']
        haplotype = row['haplotype']

        # Extract contig from haplotype
        contig = haplotype.split('#')[-1] if '#' in haplotype else haplotype

        # Parse arrays_info JSON
        arrays_info_str = row.get('arrays_info')
        if pd.isna(arrays_info_str) or not arrays_info_str or arrays_info_str == '[]':
            # No valid arrays
            rows.append({
                'chromosome': chromosome,
                'sample': sample,
                'hap': hap_num,
                'contig': contig,
                'num_arrays': 0,
                'num_canonical': 0,
                'num_degraded': 0,
                'array_types': '-',
                'array_sizes': '-'
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
            element_count = arr.get('element_count', 0)
            has_plam = arr.get('has_plam', False)

            if arr_type_str == 'canonical':
                # C+ = has β-satellite terminal, C- = no β-satellite terminal
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

    # Log summary statistics
    logger.info(f"Total haplotypes: {len(summary_df)}")
    logger.info(f"Total canonical arrays: {summary_df['num_canonical'].sum()}")
    logger.info(f"Total degraded arrays: {summary_df['num_degraded'].sum()}")

    # Array count distribution
    logger.info("Array count distribution:")
    for count, num in summary_df['num_arrays'].value_counts().sort_index().items():
        logger.info(f"  {count} arrays: {num} haplotypes")

    return summary_df


def main():
    base_dir = os.environ.get('D4Z4_BASE_DIR', '.')
    data_dir = os.path.join(base_dir, 'data')
    results_dir = os.path.join(base_dir, 'results')

    # Create results directory if it doesn't exist
    os.makedirs(results_dir, exist_ok=True)

    # Set up logging to both file and console
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(results_dir, f'd4z4_analysis_{timestamp}.log')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("Starting D4Z4 array analysis")
    logger.info(f"Data directory: {data_dir}")
    logger.info(f"Results directory: {results_dir}")

    results = []

    for chrom in ['chr4', 'chr10']:
        logger.info(f"Processing {chrom}...")

        # Load annotation files
        region_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.region.bed')
        subtelo_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.subtelomeric.bed')
        gene_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.gene.bed')
        repeat_file = os.path.join(data_dir, f'pangenome.{chrom}.subtelo_q.KS_human_CHM13.presmoothed.repeat.bed')

        if not os.path.exists(region_file):
            logger.warning(f"File not found: {region_file}")
            continue
        if not os.path.exists(subtelo_file):
            logger.warning(f"File not found: {subtelo_file}")
            continue

        region_data = load_bed_file(region_file)
        subtelo_data = load_bed_file(subtelo_file)

        # Load gene and repeat data for main array boundary detection
        gene_data = None
        repeat_data = None
        if os.path.exists(gene_file):
            gene_data = load_bed_file(gene_file)
            logger.info(f"  Loaded gene annotations")
        else:
            logger.warning(f"  Gene file not found: {gene_file}")

        if os.path.exists(repeat_file):
            repeat_data = load_bed_file(repeat_file)
            logger.info(f"  Loaded repeat annotations")
        else:
            logger.warning(f"  Repeat file not found: {repeat_file}")

        haplotypes = get_haplotypes(region_data)
        logger.info(f"  Found {len(haplotypes)} haplotypes")

        for hap in haplotypes:
            # Use v4 Array-First analysis
            result = analyze_haplotype_v4(hap, region_data, subtelo_data, gene_data, repeat_data)
            if result:
                result['chromosome'] = chrom
                results.append(result)

    # Create DataFrame and save results
    df = pd.DataFrame(results)

    # Parse sample info from haplotype name
    df['sample'] = df['haplotype'].apply(lambda x: x.split('_')[2].split('#')[0] if '_' in x else x)
    df['hap_num'] = df['haplotype'].apply(lambda x: x.split('_')[1] if '_' in x else x)

    # Serialize arrays_info and clusters_info to JSON
    df['arrays_info'] = df['arrays_info'].apply(lambda x: json.dumps(x) if x is not None else None)
    df['clusters_info'] = df['clusters_info'].apply(lambda x: json.dumps(x) if x is not None else None)

    # Reorder columns (v4 format)
    cols = ['chromosome', 'sample', 'hap_num', 'haplotype',
            'bsat_start', 'bsat_end', 'tar1_position', 'array_length', 'bsat_total_bp',
            'repeat_count_main', 'repeat_count_total',
            'num_arrays', 'num_canonical', 'num_degraded',
            'terminal_feature', 'is_truncated',
            'arrays_info', 'clusters_info', 'notes']
    df = df[[c for c in cols if c in df.columns]]

    # Sort by chromosome, then sample, then haplotype
    df = df.sort_values(['chromosome', 'sample', 'hap_num'])

    # Save to TSV
    output_file = os.path.join(results_dir, 'd4z4_array_analysis.tsv')
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Results saved to {output_file}")
    logger.info(f"Log file saved to {log_file}")

    # Print summary statistics
    logger.info("=" * 60)
    logger.info("SUMMARY STATISTICS (v4 Array-First)")
    logger.info("=" * 60)

    for chrom in ['chr4', 'chr10']:
        chrom_df = df[df['chromosome'] == chrom]
        if len(chrom_df) == 0:
            continue

        logger.info(f"{chrom.upper()} Subtelomere:")
        logger.info(f"  Total haplotypes: {len(chrom_df)}")
        logger.info(f"  Haplotypes with bsat: {len(chrom_df[chrom_df['bsat_total_bp'] > 0])}")
        logger.info(f"  Truncated arrays (> 250kb): {len(chrom_df[chrom_df['is_truncated']])}")

        valid = chrom_df[chrom_df['repeat_count_main'] > 0]
        if len(valid) > 0:
            logger.info(f"  Array length (bp): Mean={valid['array_length'].mean():.0f}, Median={valid['array_length'].median():.0f}")
            logger.info(f"  Repeat count (array elements): Mean={valid['repeat_count_main'].mean():.1f}, Median={valid['repeat_count_main'].median():.0f}, Min={valid['repeat_count_main'].min()}, Max={valid['repeat_count_main'].max()}")
            logger.info(f"  Repeat count (total clusters): Mean={valid['repeat_count_total'].mean():.1f}, Median={valid['repeat_count_total'].median():.0f}")

            # Array type distribution
            logger.info(f"  Canonical arrays: {chrom_df['num_canonical'].sum()}")
            logger.info(f"  Degraded arrays: {chrom_df['num_degraded'].sum()}")

            # Terminal feature distribution
            terminal_counts = ", ".join([f"{feat}={len(chrom_df[chrom_df['terminal_feature'] == feat])}" for feat in chrom_df['terminal_feature'].unique()])
            logger.info(f"  Terminal features: {terminal_counts}")

    # Print notable samples
    logger.info("=" * 60)
    logger.info("NOTABLE SAMPLES")
    logger.info("=" * 60)

    # FSHD1-like: Short canonical arrays with pLAM (≤10 repeats on chr4)
    fshd1_like = df[(df['chromosome'] == 'chr4') &
                    (df['repeat_count_main'] > 0) &
                    (df['repeat_count_main'] <= 10) &
                    (df['num_canonical'] > 0) &
                    (~df['is_truncated'])]
    if len(fshd1_like) > 0:
        logger.info(f"FSHD1-like arrays (chr4, ≤10 elements, canonical): {len(fshd1_like)}")
        for _, row in fshd1_like.sort_values('repeat_count_main').iterrows():
            logger.info(f"  {row['haplotype']}: {int(row['repeat_count_main'])} elements, {row['terminal_feature']}")

    # Samples without any arrays
    no_arrays = df[(df['num_arrays'] == 0) | (df['arrays_info'].isna()) | (df['arrays_info'] == '[]')]
    if len(no_arrays) > 0:
        logger.info(f"Samples without arrays: {len(no_arrays)}")
        for _, row in no_arrays.iterrows():
            logger.info(f"  {row['haplotype']}: {int(row['repeat_count_total'])} BSR clusters (no valid arrays)")

    # Multi-array samples with array makeup
    multi_array = df[df['num_arrays'] > 1]
    if len(multi_array) > 0:
        logger.info(f"Multi-array samples: {len(multi_array)}")
        for _, row in multi_array.sort_values(['chromosome', 'num_arrays'], ascending=[True, False]).iterrows():
            try:
                arrays = json.loads(row['arrays_info'])
                # C+ = canonical with β-satellite, C- = canonical without, D = degraded
                types = []
                for a in arrays:
                    if a['type'] == 'canonical':
                        types.append('C+' if a.get('has_plam', False) else 'C-')
                    else:
                        types.append('D')
                sizes = [str(a['element_count']) for a in arrays]
                logger.info(f"  {row['haplotype']}: {','.join(types)} ({','.join(sizes)} elements)")
            except (json.JSONDecodeError, TypeError, KeyError):
                logger.info(f"  {row['haplotype']}: {row['num_arrays']} arrays")

    # Mixed arrays (both canonical + degraded on same haplotype)
    def has_mixed_arrays(row):
        if pd.isna(row.get('arrays_info')) or row['num_arrays'] <= 1:
            return False
        try:
            arrays = json.loads(row['arrays_info'])
            types = set(a.get('type') for a in arrays)
            return 'canonical' in types and 'degraded' in types
        except:
            return False

    mixed = df[df.apply(has_mixed_arrays, axis=1)]
    if len(mixed) > 0:
        logger.info(f"Mixed arrays (canonical + degraded): {len(mixed)}")
        for _, row in mixed.sort_values('chromosome').iterrows():
            try:
                arrays = json.loads(row['arrays_info'])
                types = []
                for a in arrays:
                    if a['type'] == 'canonical':
                        types.append('C+' if a.get('has_plam', False) else 'C-')
                    else:
                        types.append('D')
                sizes = [str(a['element_count']) for a in arrays]
                logger.info(f"  {row['haplotype']}: {','.join(types)} ({','.join(sizes)} elements)")
            except (json.JSONDecodeError, TypeError, KeyError):
                logger.info(f"  {row['haplotype']}: mixed arrays")

    # Only degraded arrays (no canonical - cannot express DUX4)
    def has_only_degraded(row):
        if pd.isna(row.get('arrays_info')) or row['num_arrays'] == 0:
            return False
        try:
            arrays = json.loads(row['arrays_info'])
            types = set(a.get('type') for a in arrays)
            return types == {'degraded'}
        except:
            return False

    only_degraded = df[df.apply(has_only_degraded, axis=1)]
    if len(only_degraded) > 0:
        logger.info(f"Only degraded arrays (no DUX4 expression): {len(only_degraded)}")
        for _, row in only_degraded.sort_values('chromosome').iterrows():
            try:
                arrays = json.loads(row['arrays_info'])
                sizes = [f"D({a['element_count']})" for a in arrays]
                logger.info(f"  {row['haplotype']}: {','.join(sizes)}")
            except (json.JSONDecodeError, TypeError, KeyError):
                logger.info(f"  {row['haplotype']}: degraded only")

    # Truncated arrays
    truncated = df[df['is_truncated']]
    if len(truncated) > 0:
        logger.info(f"Truncated arrays (extend beyond 250kb window): {len(truncated)}")
        for _, row in truncated.iterrows():
            logger.info(f"  {row['haplotype']}: >{int(row['repeat_count_main'])} elements, {row['terminal_feature']}")

    # Functional pLAM on chr10 (rare - potential 4q/10q translocation)
    plam_file = os.path.join(results_dir, 'd4z4_array_coordinates.plam.tsv')
    if os.path.exists(plam_file):
        plam_df = pd.read_csv(plam_file, sep='\t')
        chr10_func = plam_df[(plam_df['chromosome'] == 'chr10') & (plam_df['plam_functional'] == 'func_Plam')]
        if len(chr10_func) > 0:
            logger.info(f"Functional pLAM on chr10 (rare): {len(chr10_func)}")
            for _, row in chr10_func.iterrows():
                hap_match = df[(df['sample'] == row['sample']) & (df['chromosome'] == 'chr10') & (df['hap_num'] == row['hap'])]
                if len(hap_match) > 0:
                    logger.info(f"  {hap_match.iloc[0]['haplotype']}: {row['element_count']} elements, {row['plam_functional']}")

    # Generate and save summary table
    logger.info("=" * 60)
    logger.info("GENERATING SUMMARY TABLE")
    logger.info("=" * 60)
    summary_file = os.path.join(results_dir, 'd4z4_array_summary.tsv')
    summary_df = generate_array_summary(df, summary_file, logger)
    logger.info(f"Summary saved to {summary_file}")

    logger.info("Analysis complete (v4)")

if __name__ == '__main__':
    main()
