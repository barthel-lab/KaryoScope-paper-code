#!/usr/bin/env python3
# Vendored from barthel-lab/KaryoScope-ISCN at commit 4e42a98d6526ecb31ca6fb5d6c8910e4e84c605a, file KaryoScope_assembly_contig_zoom_plot.py.
# Copied unmodified. Authors (git log): Floris Barthel <fbarthel@tgen.org>; Rhyker Ranallo-Benavidez <rranallo-benavidez@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""
KaryoScope Assembly Contig Zoom Plot

Generates multi-track zoom plots for dicentric centromeres and translocation
breakpoints identified by KaryoScope_ISCN_report.py.

Features:
- 5 feature tracks side-by-side (chromosome, subtelomeric, repeat, gene, region)
- Multiple zoom levels (10Mb, 5Mb, 2Mb, 500kb)
- Zoom-adaptive binning for visualization clarity
- Scale bars and position axes
- Region indicators for centromeres/breakpoints

Usage:
    python KaryoScope_assembly_contig_zoom_plot.py \\
        --regions sample_ISCN_regions.tsv \\
        --bed-dir results/sample/assembly/ \\
        --colors-dir resources/databases/KS_human_CHM13/ \\
        --output-dir output/
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from math import floor

import drawsvg as draw


# --- CLI Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multi-track zoom plots for cytogenetic regions",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--regions', required=True, action='append',
                        help='Regions TSV from KaryoScope_ISCN_report.py. '
                             'May be repeated (once per sample) in stacked mode.')
    parser.add_argument('--bed-dir', required=True, action='append',
                        help='Base directory containing assembly BED files. '
                             'May be repeated (matched to --sample and --regions).')
    parser.add_argument('--colors-dir', required=True,
                        help='Directory with color files')
    parser.add_argument('--database', default='KS_human_CHM13',
                        help='Database name (default: KS_human_CHM13)')
    parser.add_argument('--sample', required=True, action='append',
                        help='Sample name for BED file pattern matching. '
                             'May be repeated (matched to --bed-dir and --regions).')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for per-scaffold SVG files')
    parser.add_argument('--stacked-output', default=None,
                        help='When set (and multiple samples/regions given), '
                             'produce one combined SVG stacking all contigs '
                             'vertically with a single shared legend.')
    parser.add_argument('--featuresets',
                        default='chromosome,subtelomeric,repeat,gene,region',
                        help='Comma-separated featuresets to display')
    parser.add_argument('--smoothness', default='presmoothed',
                        help='BED smoothness level (default: presmoothed)')
    parser.add_argument('--zoom-levels', default='10000000,5000000,2000000,500000',
                        help='Comma-separated zoom paddings in bp')
    parser.add_argument('--background', default='black', choices=['white', 'black'],
                        help='Background color')
    parser.add_argument('--bar-width', type=int, default=20,
                        help='Width of each feature track in pixels')
    parser.add_argument('--bar-spacing', type=int, default=2,
                        help='Spacing between tracks in pixels')
    parser.add_argument('--plot-height', type=int, default=500,
                        help='Height of feature tracks in pixels')
    parser.add_argument('--comparison', default=None, action='append',
                        help='KromaTid comparison TSV for alternate naming. '
                             'May be repeated to match --sample.')
    parser.add_argument('--show-arms', action='store_true',
                        help='Show p/q arm track next to chromosome (contig view only)')
    parser.add_argument('--no-refine', action='store_true', default=False,
                        help='Skip breakpoint refinement; use coordinates from '
                             'the regions TSV as-is')
    parser.add_argument('--hide-scale-bars', action='store_true', default=False,
                        help='Hide scale bars and CI lines in zoom panels')
    parser.add_argument('--panel-letters', action='store_true', default=False,
                        help='Label each stacked row A, B, C ... at its top-left '
                             '(stacked mode only)')
    parser.add_argument('--show-contig-breakpoints', action='store_true', default=False,
                        help='Show dashed breakpoint lines on the full contig view '
                             '(stacked mode only; integrated mode always shows them)')

    return parser.parse_args()


# --- Scaffold Name Handling ---
def get_haplotype_from_scaffold(scaffold_name):
    """Extract haplotype from scaffold name."""
    if '_hap1_' in scaffold_name or scaffold_name.endswith('_hap1'):
        return 'hap1'
    elif '_hap2_' in scaffold_name or scaffold_name.endswith('_hap2'):
        return 'hap2'
    return None


def get_contig_name(scaffold_name):
    """
    Extract contig name from full scaffold name.

    Full scaffold name format: {assigned_chrom}_{hap}_{contig_name}
    Example: chr20_hap1_chr18_chr20_hap1 -> chr18_chr20_hap1

    Returns the contig name used in feature BED files.
    """
    # Split on haplotype marker
    if '_hap1_' in scaffold_name:
        parts = scaffold_name.split('_hap1_', 1)
        if len(parts) == 2:
            return parts[1]
    elif '_hap2_' in scaffold_name:
        parts = scaffold_name.split('_hap2_', 1)
        if len(parts) == 2:
            return parts[1]
    # Fallback: return original name
    return scaffold_name


#: Fill opacity of the grey connector wedges between zoom levels. 0.2 matches
#: the published Figure 3; the value was 0.5, which the figure's hand-finishing
#: pass lightened. Over a white background 0.2 renders #606060 as #DFDFDF.
FAN_FILL_OPACITY = 0.2


# --- Binning Configuration ---
def get_bin_size_for_zoom(padding):
    """Choose bin size based on zoom level."""
    if padding >= 10_000_000:
        return 100_000   # 100kb bins for 10Mb+ zoom
    elif padding >= 5_000_000:
        return 50_000    # 50kb bins for 5Mb zoom
    elif padding >= 2_000_000:
        return 10_000    # 10kb bins for 2Mb zoom
    else:
        return 1_000     # 1kb bins for <2Mb zoom

def discover_prebinned_resolutions(bed_dir, sample, database, featuresets, smoothness):
    """
    Discover available pre-binned resolutions by scanning BED filenames.

    Looks for files matching:
        {sample}.{database}.{fs}.{smoothness}.scaffolded.binned{N}.features.bed[.gz]

    Returns sorted list of resolution ints (descending), e.g. [1000000, 100000, 10000, 1000, 100].
    Only includes resolutions that exist for ALL requested featuresets.
    """
    import re

    pattern = re.compile(
        rf'^{re.escape(sample)}\.{re.escape(database)}\.'
        rf'(\w+)\.{re.escape(smoothness)}'
        rf'\.scaffolded\.binned(\d+)\.features\.bed(\.gz)?$'
    )

    # Track which resolutions exist for each featureset
    resolutions_by_fs = {fs: set() for fs in featuresets}

    for filename in os.listdir(bed_dir):
        m = pattern.match(filename)
        if m:
            fs_name = m.group(1)
            resolution = int(m.group(2))
            if fs_name in resolutions_by_fs:
                resolutions_by_fs[fs_name].add(resolution)

    # Intersect: only keep resolutions available for ALL featuresets
    if not resolutions_by_fs:
        return []

    common = set.intersection(*resolutions_by_fs.values()) if resolutions_by_fs else set()

    if not common:
        print("Warning: No pre-binned resolutions found common to all featuresets",
              file=sys.stderr)
        # Fall back to union so we at least have something
        common = set.union(*resolutions_by_fs.values())

    result = sorted(common, reverse=True)
    return result


def get_prebinned_resolution(target_bin_size, available_prebinned):
    """
    Select the best pre-binned file resolution for a target bin size.

    Args:
        target_bin_size: Desired bin size in bp
        available_prebinned: List of available resolutions (sorted descending)

    Returns (prebinned_size, needs_rebinning) tuple.

    Strategy: Use the largest pre-binned resolution that:
      1. Is <= the target bin size, AND
      2. Evenly divides the target (no partial bins at re-bin boundaries)

    This guarantees that every re-bin boundary aligns with a pre-bin
    boundary, so max-overlap voting sees only whole pre-bins per target
    bin—eliminating tie artifacts from partial overlaps.
    """
    for prebinned in available_prebinned:
        if prebinned <= target_bin_size and target_bin_size % prebinned == 0:
            return prebinned, (prebinned != target_bin_size)

    # Fallback: smallest available
    if available_prebinned:
        smallest = available_prebinned[-1]
        if target_bin_size % smallest != 0:
            print(f"Warning: target bin size {target_bin_size} is not evenly divisible "
                  f"by smallest pre-binned resolution {smallest}bp", file=sys.stderr)
        return smallest, True

    # Nothing found at all
    print(f"Error: No pre-binned files found", file=sys.stderr)
    return target_bin_size, False

def get_scale_bar_config(padding):
    """Get scale bar length and label for zoom level."""
    if padding >= 10_000_000:
        return 2_000_000, "2 Mb"
    elif padding >= 5_000_000:
        return 1_000_000, "1 Mb"
    elif padding >= 2_000_000:
        return 500_000, "500 kb"
    elif padding >= 500_000:
        return 100_000, "100 kb"
    else:
        return 50_000, "50 kb"


def get_tick_interval(padding):
    """Get position axis tick interval for zoom level."""
    if padding >= 10_000_000:
        return 5_000_000   # 5Mb ticks
    elif padding >= 5_000_000:
        return 2_000_000   # 2Mb ticks
    elif padding >= 2_000_000:
        return 1_000_000   # 1Mb ticks
    elif padding >= 500_000:
        return 200_000     # 200kb ticks
    else:
        return 50_000      # 50kb ticks


def compute_view_window(center, padding, scaffold_length=None):
    """Compute a view window of size 2*padding centered on *center*.

    If the window would extend past 0 or *scaffold_length*, it is shifted
    (not truncated) so the full 2*padding span is shown whenever possible.
    """
    view_start = center - padding
    view_end = center + padding

    if view_start < 0:
        # Shift window right so it starts at 0
        view_end += -view_start
        view_start = 0

    if scaffold_length is not None and view_end > scaffold_length:
        # Shift window left so it ends at scaffold_length
        overshoot = view_end - scaffold_length
        view_start = max(0, view_start - overshoot)
        view_end = scaffold_length

    return view_start, view_end


def get_scaffold_length(bed_path, scaffold):
    """Get total length of a scaffold from BED data."""
    load_bed_file_to_cache(bed_path)
    cached = BED_CACHE.get(bed_path, {})
    scaffold_data = cached.get(scaffold, [])
    if not scaffold_data:
        return 0
    return max(end for _, end, _ in scaffold_data)


# --- File Loading ---
def load_colors(colors_dir, database, featuresets, background='black'):
    """Load color mappings for featuresets via karyoplot's library."""
    from karyoplot.core.colors import load_featureset_palettes
    return load_featureset_palettes(
        colors_dir, database, featuresets,
        on_missing="warn", value_format="tuple",
        background=background,
    )


def parse_regions_file(regions_file):
    """Parse regions TSV from ISCN report.

    Uses the header row for column lookup so that both the original format
    (without ``original_scaffold``) and the updated format (with it) are
    handled transparently.
    """
    regions = []
    with open(regions_file, 'r') as f:
        header = f.readline().rstrip('\n\r').split('\t')
        col = {name: idx for idx, name in enumerate(header)}

        # Required columns — fall back to positional indices if headers
        # are missing (backward compat with hand-edited files).
        def _idx(name, fallback):
            return col.get(name, fallback)

        i_scaffold    = _idx('scaffold', 0)
        i_region_type = _idx('region_type', 1)
        i_center      = _idx('center', 2)
        i_start       = _idx('start', 3)
        i_end         = _idx('end', 4)
        i_chrom1      = _idx('chrom1', 5)
        i_chrom2      = _idx('chrom2', 6)
        i_arm1        = _idx('arm1', 7)
        i_arm2        = _idx('arm2', 8)
        i_iscn        = _idx('iscn_notation', 9)

        min_cols = max(i_scaffold, i_region_type, i_center, i_start, i_end,
                       i_chrom1, i_chrom2, i_arm1, i_arm2, i_iscn) + 1

        for line in f:
            parts = line.rstrip('\n\r').split('\t')
            if len(parts) < min_cols:
                continue
            region = {
                'scaffold':      parts[i_scaffold],
                'region_type':   parts[i_region_type],
                'center':        int(parts[i_center]) if parts[i_center] else 0,
                'start':         int(parts[i_start]) if parts[i_start] else None,
                'end':           int(parts[i_end]) if parts[i_end] else None,
                'chrom1':        parts[i_chrom1] if parts[i_chrom1] else None,
                'chrom2':        parts[i_chrom2] if parts[i_chrom2] else None,
                'arm1':          parts[i_arm1] if parts[i_arm1] else None,
                'arm2':          parts[i_arm2] if parts[i_arm2] else None,
                'iscn_notation': parts[i_iscn] if parts[i_iscn] else ''
            }
            regions.append(region)
    return regions


def _extract_chrom_from_feature(feature_name):
    """Extract chromosome from feature like 'chr3_specific' -> 'chr3'."""
    if feature_name.endswith('_specific') and feature_name.startswith('chr'):
        return feature_name.replace('_specific', '')
    if feature_name.startswith('chr') and '_' not in feature_name:
        return feature_name
    return None


def refine_breakpoint_center(region, scaffold, bed_dir, sample, database,
                             smoothness, available_prebinned, search_radius=None):
    """
    Refine a breakpoint center coordinate using the finest-resolution BED data.

    Uses a content-dominance approach: aggregates chromosome-specific content
    into windows, identifies the zone where one chromosome's dominance ends
    and the other's begins, and places the breakpoint at the midpoint of the
    intervening non-dominant (multigroup/ambiguous) zone.

    This correctly handles real-world breakpoints where the transition between
    two chromosomes passes through a large block of shared/ambiguous sequence
    (e.g. autosome_multigroup) that may contain tiny specks of chromosome-
    specific features without being truly dominated by either chromosome.

    Args:
        region:             Region dict with 'center', 'start', 'end', 'chrom1', 'chrom2'
        scaffold:           Scaffold name
        bed_dir:            Directory with pre-binned BED files
        sample:             Sample name
        database:           Database name
        smoothness:         Smoothness level string
        available_prebinned: Sorted list of available pre-binned resolutions (desc)
        search_radius:      Override search radius in bp (default: use interval or 5 Mb)

    Returns:
        Refined center coordinate (int).  Falls back to the original center if
        no clear chromosome transition is found.
    """
    original_center = region['center']

    # Scan the entire scaffold rather than a small window.  The coarse
    # breakpoint estimate can be several Mb (even 10+ Mb) off when large
    # multigroup / ambiguous zones separate the two chromosomes, so a
    # narrow window around the estimate may miss the real transition
    # entirely.  With 100 kb aggregation windows this is still trivially
    # fast even for 200 Mb scaffolds (~2000 windows).
    if search_radius is not None:
        window_start = max(0, original_center - search_radius)
        window_end = original_center + search_radius
    else:
        # Full scaffold scan
        window_start = 0
        window_end = 500_000_000  # generous upper bound; load_bed_features clips

    # Pick the finest available pre-binned resolution
    if not available_prebinned:
        return original_center
    finest_resolution = available_prebinned[-1]  # list is sorted descending

    # Load chromosome featureset at finest resolution
    bed_path = os.path.join(
        bed_dir,
        f"{sample}.{database}.chromosome.{smoothness}"
        f".scaffolded.binned{finest_resolution}.features.bed.gz"
    )
    if not os.path.exists(bed_path):
        bed_path = bed_path.replace('.bed.gz', '.bed')
    if not os.path.exists(bed_path):
        return original_center

    features = load_bed_features(bed_path, scaffold, window_start, window_end)
    if not features:
        return original_center

    # Identify the two chromosomes flanking this breakpoint
    chrom1 = region.get('chrom1') or ''
    chrom2 = region.get('chrom2') or ''

    if not chrom1 or not chrom2:
        return original_center

    feature_list = sorted(features, key=lambda f: f['start'])

    # === Content-dominance approach ===
    # Aggregate chromosome-specific content into fixed-size windows and
    # classify each window as dominated by chrom1, chrom2, or neither.
    # Then find the transition zone between the last chrom-A-dominant
    # window and the first chrom-B-dominant window.
    #
    # Use 100kb windows — large enough to smooth over tiny specks but
    # small enough to resolve the transition zone.
    dominance_window = 100_000

    # Aggregate content per window
    window_content = defaultdict(lambda: defaultdict(int))  # pos -> chrom -> bp
    for feat in feature_list:
        chrom = _extract_chrom_from_feature(feat['feature'])
        if chrom is None:
            continue
        feat_start = feat['start']
        feat_end = feat['end']
        # Distribute across windows
        w_start = (feat_start // dominance_window) * dominance_window
        while w_start < feat_end:
            w_end = w_start + dominance_window
            overlap_start = max(feat_start, w_start)
            overlap_end = min(feat_end, w_end)
            overlap = overlap_end - overlap_start
            if overlap > 0:
                window_content[w_start][chrom] += overlap
            w_start += dominance_window

    if not window_content:
        return original_center

    # Classify each window by which of our two target chromosomes dominates.
    # "dominant" means the chromosome covers a substantial fraction of the
    # TOTAL window size — not just of chromosome-specific content.  This
    # prevents sparse multigroup zones with tiny chromosome specks from
    # being misclassified as chromosome-dominant.
    #
    # Thresholds:
    #   - The dominant chromosome must cover >= 25% of the window size
    #   - AND it must be >= 50% of the chromosome-specific content
    min_window_fraction = 0.25   # 25kb out of 100kb window
    specificity_threshold = 0.50  # >=50% of chr-specific content

    window_dominance = {}  # pos -> 'chrom1' | 'chrom2' | None
    for w_pos in sorted(window_content.keys()):
        content = window_content[w_pos]
        c1_bp = content.get(chrom1, 0)
        c2_bp = content.get(chrom2, 0)
        total_chr_bp = c1_bp + c2_bp

        c1_dominant = (c1_bp >= dominance_window * min_window_fraction and
                       total_chr_bp > 0 and
                       c1_bp / total_chr_bp >= specificity_threshold)
        c2_dominant = (c2_bp >= dominance_window * min_window_fraction and
                       total_chr_bp > 0 and
                       c2_bp / total_chr_bp >= specificity_threshold)

        if c1_dominant and not c2_dominant:
            window_dominance[w_pos] = chrom1
        elif c2_dominant and not c1_dominant:
            window_dominance[w_pos] = chrom2
        else:
            window_dominance[w_pos] = None  # Ambiguous / multigroup zone

    # Find the transition: scan windows in order, looking for the boundary
    # between chrom1-dominant and chrom2-dominant zones (in either direction).
    sorted_windows = sorted(window_dominance.keys())
    if not sorted_windows:
        return original_center

    # Find the last window dominated by each chromosome and the first window
    # dominated by the other, scanning outward from the original center.
    # This handles both orientations (chrom1 before chrom2 or vice versa).

    # Determine which chromosome is "upstream" (before center) and which
    # is "downstream" by looking at overall content distribution
    upstream_chrom = None
    downstream_chrom = None
    c1_before = sum(window_content[w].get(chrom1, 0)
                    for w in sorted_windows if w < original_center)
    c2_before = sum(window_content[w].get(chrom2, 0)
                    for w in sorted_windows if w < original_center)
    if c1_before >= c2_before:
        upstream_chrom = chrom1
        downstream_chrom = chrom2
    else:
        upstream_chrom = chrom2
        downstream_chrom = chrom1

    # Find last upstream-dominant window and first downstream-dominant window
    last_upstream_end = None
    first_downstream_start = None

    for w in sorted_windows:
        if window_dominance[w] == upstream_chrom:
            last_upstream_end = w + dominance_window

    for w in sorted_windows:
        if window_dominance[w] == downstream_chrom:
            first_downstream_start = w
            break

    if last_upstream_end is not None and first_downstream_start is not None:
        gap = first_downstream_start - last_upstream_end

        if gap > dominance_window:
            # Wide multigroup zone — use midpoint
            return (last_upstream_end + first_downstream_start) // 2

        # Sharp or narrow transition — drill down at display resolution
        # (1 kb bins) to find where each chromosome actually wins the
        # max-overlap vote in binned plots.  Using raw 100 bp bins would
        # place the breakpoint at tiny chromosome specks that are invisible
        # at any display zoom level.
        fine_bin = 1_000  # 1 kb — matches finest zoom level
        # Use a wider window: the downstream chromosome may be sparse for
        # several Mb before it starts winning display bins consistently
        fine_start = max(0, last_upstream_end - dominance_window)
        fine_end = first_downstream_start + 20 * dominance_window  # ±2 Mb downstream
        fine_features = load_bed_features(bed_path, scaffold, fine_start, fine_end)
        if fine_features:
            fine_list = sorted(fine_features, key=lambda f: f['start'])
            # Aggregate ALL features into 1 kb bins, then find which
            # chromosome wins each bin via max-overlap voting (same logic
            # as the display binning).
            fine_bins = defaultdict(lambda: defaultdict(int))
            for feat in fine_list:
                sb = feat['start'] // fine_bin
                eb = (feat['end'] - 1) // fine_bin
                for b_idx in range(sb, eb + 1):
                    bs = b_idx * fine_bin
                    be = bs + fine_bin
                    ovlp = min(feat['end'], be) - max(feat['start'], bs)
                    if ovlp > 0:
                        fine_bins[b_idx][feat['feature']] += ovlp

            last_up_pos = None
            first_down_pos = None
            up_feat = f"{upstream_chrom}_specific"
            down_feat = f"{downstream_chrom}_specific"
            prev_was_down = False
            for b_idx in sorted(fine_bins.keys()):
                winner = max(fine_bins[b_idx].items(), key=lambda x: x[1])[0]
                pos = b_idx * fine_bin
                if winner == up_feat:
                    last_up_pos = pos + fine_bin
                elif winner == down_feat:
                    # Only count the downstream start if it's a solid
                    # transition: the win is either consecutive with the
                    # previous bin (forming a real block) or it's the very
                    # first bin after the upstream chromosome (no gap for
                    # scattered specks to appear in).
                    if first_down_pos is None:
                        if prev_was_down:
                            # Consecutive: this is real.  Record the
                            # position of the FIRST of the consecutive pair.
                            first_down_pos = pos - fine_bin
                        # else: isolated win — don't record yet
                prev_was_down = (winner == down_feat)

            # If no consecutive downstream wins found, check if the first
            # downstream win is directly adjacent to the last upstream win
            # (gap <= 1 bin).  That's a genuinely sharp transition.
            if first_down_pos is None and last_up_pos is not None:
                for b_idx in sorted(fine_bins.keys()):
                    winner = max(fine_bins[b_idx].items(), key=lambda x: x[1])[0]
                    pos = b_idx * fine_bin
                    if winner == down_feat and pos <= last_up_pos + fine_bin:
                        first_down_pos = pos
                        break

            if last_up_pos is not None and first_down_pos is not None:
                return (last_up_pos + first_down_pos) // 2

        # Fallback to dominance-window midpoint
        return (last_upstream_end + first_downstream_start) // 2

    # Fallback: return original
    return original_center


def refine_breakpoint_interval(region, refined_center, scaffold, bed_dir,
                                sample, database, smoothness, available_prebinned):
    """
    Refine the breakpoint interval (start/end) using finest-resolution data.

    Uses a content-dominance approach matching refine_breakpoint_center():
    finds the last window dominated by the upstream chromosome and the first
    window dominated by the downstream chromosome, and returns those as the
    interval bounds.  This gives a tight interval around the ambiguous zone.

    Returns (refined_start, refined_end) tuple.
    """
    if not available_prebinned:
        return region.get('start'), region.get('end')

    chrom1 = region.get('chrom1') or ''
    chrom2 = region.get('chrom2') or ''
    if not chrom1 or not chrom2:
        return region.get('start'), region.get('end')

    finest_resolution = available_prebinned[-1]

    bed_path = os.path.join(
        bed_dir,
        f"{sample}.{database}.chromosome.{smoothness}"
        f".scaffolded.binned{finest_resolution}.features.bed.gz"
    )
    if not os.path.exists(bed_path):
        bed_path = bed_path.replace('.bed.gz', '.bed')
    if not os.path.exists(bed_path):
        return region.get('start'), region.get('end')

    # Scan full scaffold (same rationale as refine_breakpoint_center)
    features = load_bed_features(bed_path, scaffold, 0, 500_000_000)
    if not features:
        return region.get('start'), region.get('end')

    feature_list = sorted(features, key=lambda f: f['start'])

    # Aggregate into 100kb dominance windows (same approach as refine_breakpoint_center)
    dominance_window = 100_000
    window_content = defaultdict(lambda: defaultdict(int))
    for feat in feature_list:
        chrom = _extract_chrom_from_feature(feat['feature'])
        if chrom is None:
            continue
        w_start = (feat['start'] // dominance_window) * dominance_window
        while w_start < feat['end']:
            w_end = w_start + dominance_window
            overlap = min(feat['end'], w_end) - max(feat['start'], w_start)
            if overlap > 0:
                window_content[w_start][chrom] += overlap
            w_start += dominance_window

    if not window_content:
        return region.get('start'), region.get('end')

    min_window_fraction = 0.25
    specificity_threshold = 0.50

    sorted_windows = sorted(window_content.keys())

    # Determine upstream/downstream orientation
    c1_before = sum(window_content[w].get(chrom1, 0)
                    for w in sorted_windows if w < refined_center)
    c2_before = sum(window_content[w].get(chrom2, 0)
                    for w in sorted_windows if w < refined_center)
    upstream_chrom = chrom1 if c1_before >= c2_before else chrom2
    downstream_chrom = chrom2 if upstream_chrom == chrom1 else chrom1

    last_upstream_end = None
    first_downstream_start = None

    for w in sorted_windows:
        content = window_content[w]
        c_up = content.get(upstream_chrom, 0)
        c_down = content.get(downstream_chrom, 0)
        total = c_up + c_down

        up_dominant = (c_up >= dominance_window * min_window_fraction and
                       total > 0 and c_up / total >= specificity_threshold)
        down_dominant = (c_down >= dominance_window * min_window_fraction and
                         total > 0 and c_down / total >= specificity_threshold)

        if up_dominant:
            last_upstream_end = w + dominance_window
        if down_dominant and first_downstream_start is None:
            first_downstream_start = w

    if last_upstream_end is not None and first_downstream_start is not None:
        return last_upstream_end, first_downstream_start

    # Fallback: expand slightly around refined center for visibility
    return max(0, refined_center - 500_000), refined_center + 500_000


def parse_comparison_file(comparison_file):
    """
    Parse KromaTid comparison TSV to create KaryoScope -> KromaTid notation mapping.
    Returns dict mapping KaryoScope notation to KromaTid notation.
    """
    if not comparison_file or not os.path.exists(comparison_file):
        return {}

    mapping = {}
    with open(comparison_file, 'r') as f:
        header = f.readline().strip().split('\t')
        # Find column indices
        ks_col = None
        kt_col = None
        for i, col in enumerate(header):
            if col in ('KaryoScope', 'KaryoScope_retained', 'KaryoScope_breakpoint'):
                if ks_col is None:  # Use first KS column found
                    ks_col = i
            if col in ('KromaTid', 'KromaTid_call'):
                kt_col = i

        if ks_col is None or kt_col is None:
            return {}

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(ks_col, kt_col):
                ks_notation = parts[ks_col]
                kt_notation = parts[kt_col]
                if ks_notation and kt_notation and ks_notation != '.' and kt_notation != '.':
                    mapping[ks_notation] = kt_notation

    return mapping


# Global cache for BED data: {bed_path: {scaffold: [(start, end, feature), ...]}}
BED_CACHE = {}


def load_bed_file_to_cache(bed_path):
    """Load entire BED file into cache, indexed by scaffold."""
    if bed_path in BED_CACHE:
        return

    if not os.path.exists(bed_path):
        BED_CACHE[bed_path] = {}
        return

    scaffold_features = defaultdict(list)
    open_func = gzip.open if bed_path.endswith(".gz") else open
    mode = "rt" if bed_path.endswith(".gz") else "r"

    with open_func(bed_path, mode) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            scaf, start, end, feature = parts[:4]
            try:
                start, end = int(start), int(end)
            except ValueError:
                continue
            scaffold_features[scaf].append((start, end, feature))

    # Sort each scaffold's features by start position
    for scaf in scaffold_features:
        scaffold_features[scaf].sort(key=lambda x: x[0])

    BED_CACHE[bed_path] = dict(scaffold_features)


def load_bed_features(bed_path, scaffold, view_start, view_end):
    """Load features for a scaffold within view window from cache."""
    # Ensure file is cached
    load_bed_file_to_cache(bed_path)

    features = []
    cached = BED_CACHE.get(bed_path, {})
    scaffold_data = cached.get(scaffold, [])

    for start, end, feature in scaffold_data:
        # Check overlap with view window
        if end > view_start and start < view_end:
            # Clip to view window
            clipped_start = max(start, view_start)
            clipped_end = min(end, view_end)
            features.append({
                'start': clipped_start,
                'end': clipped_end,
                'feature': feature
            })

    return features


def bin_features(features, view_start, view_end, bin_size):
    """
    Bin features using max-overlap logic from bin_features.py.
    Returns list of binned features with merged adjacent same-feature bins.
    """
    if not features:
        return []

    # Calculate bin assignments
    bins = defaultdict(lambda: defaultdict(int))  # bin_idx -> feature -> overlap

    for feat in features:
        start_bin = feat['start'] // bin_size
        end_bin = (feat['end'] - 1) // bin_size

        for b_idx in range(start_bin, end_bin + 1):
            bin_start = b_idx * bin_size
            bin_end = bin_start + bin_size
            overlap_start = max(feat['start'], bin_start)
            overlap_end = min(feat['end'], bin_end)
            overlap_len = overlap_end - overlap_start
            if overlap_len > 0:
                bins[b_idx][feat['feature']] += overlap_len

    # Get best feature for each bin
    binned = []
    for b_idx in sorted(bins.keys()):
        bin_start = b_idx * bin_size
        bin_end = bin_start + bin_size

        # Clip to view window
        if bin_end <= view_start or bin_start >= view_end:
            continue
        clipped_start = max(bin_start, view_start)
        clipped_end = min(bin_end, view_end)

        # Get feature with max overlap
        # Novel only wins if it covers >50% of the bin
        best_feature = max(bins[b_idx], key=lambda f: (bins[b_idx][f], f != 'novel'))
        if best_feature in ('novel'):
            if bins[b_idx][best_feature] <= bin_size / 2:
                non_depr = {f: v for f, v in bins[b_idx].items() if f not in ('novel')}
                if non_depr:
                    best_feature = max(non_depr, key=non_depr.get)
        binned.append({
            'start': clipped_start,
            'end': clipped_end,
            'feature': best_feature
        })

    # Merge adjacent bins with same feature
    merged = []
    for feat in binned:
        if merged and merged[-1]['end'] == feat['start'] and merged[-1]['feature'] == feat['feature']:
            merged[-1]['end'] = feat['end']
        else:
            merged.append(feat.copy())

    return merged


# --- Drawing Functions ---
def draw_scale_bar(d, x, y, scale_length_bp, pixels_per_bp, text_color):
    """Draw horizontal scale bar with label."""
    bar_width = scale_length_bp * pixels_per_bp
    bar_height = 3

    # Draw bar
    d.append(draw.Rectangle(x, y, bar_width, bar_height, fill=text_color))

    # Draw end ticks
    tick_height = 8
    d.append(draw.Line(x, y - tick_height/2, x, y + tick_height/2,
                       stroke=text_color, stroke_width=1))
    d.append(draw.Line(x + bar_width, y - tick_height/2, x + bar_width, y + tick_height/2,
                       stroke=text_color, stroke_width=1))

    # Draw label
    scale_length, label = get_scale_bar_config(int(scale_length_bp * 2))  # Approximate padding
    d.append(draw.Text(label, font_size=13, x=x + bar_width/2, y=y + 15,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle'))


def draw_vertical_position_axis(d, x, y_start, y_end, view_start, view_end, pixels_per_bp, text_color, padding):
    """Draw vertical position axis on the left side with tick marks."""
    tick_interval = get_tick_interval(padding)

    # Draw main vertical line
    d.append(draw.Line(x, y_start, x, y_end, stroke=text_color, stroke_width=1))

    # Draw ticks and labels
    first_tick = ((view_start // tick_interval) + 1) * tick_interval
    tick_pos = first_tick
    while tick_pos < view_end:
        y = y_start + (tick_pos - view_start) * pixels_per_bp
        if y_start <= y <= y_end:
            # Tick mark (pointing left)
            d.append(draw.Line(x - 5, y, x, y, stroke=text_color, stroke_width=1))
            # Label
            if tick_pos >= 1_000_000:
                label = f"{tick_pos / 1_000_000:.0f}Mb"
            else:
                label = f"{tick_pos / 1_000:.0f}kb"
            d.append(draw.Text(label, font_size=11, x=x - 8, y=y + 3,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='end'))
        tick_pos += tick_interval


def draw_vertical_scale_bar(d, x, y, scale_length_bp, pixels_per_bp, text_color):
    """Draw vertical scale bar showing genomic distance."""
    bar_height = scale_length_bp * pixels_per_bp
    bar_width = 3

    # Draw vertical bar
    d.append(draw.Rectangle(x, y, bar_width, bar_height, fill=text_color))

    # Draw end ticks
    tick_width = 8
    d.append(draw.Line(x - tick_width/2 + bar_width/2, y,
                       x + tick_width/2 + bar_width/2, y,
                       stroke=text_color, stroke_width=1))
    d.append(draw.Line(x - tick_width/2 + bar_width/2, y + bar_height,
                       x + tick_width/2 + bar_width/2, y + bar_height,
                       stroke=text_color, stroke_width=1))

    # Label (rotated vertically)
    if scale_length_bp >= 1_000_000:
        label = f"{scale_length_bp / 1_000_000:.0f} Mb"
    elif scale_length_bp >= 1_000:
        label = f"{scale_length_bp / 1_000:.0f} kb"
    else:
        label = f"{scale_length_bp} bp"

    label_x = x + bar_width + 5
    label_y = y + bar_height / 2
    d.append(draw.Text(label, font_size=11, x=label_x, y=label_y,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle', dominant_baseline='middle',
                       transform=f"rotate(-90 {label_x} {label_y})"))


def _format_genomic_pos(bp_pos):
    """Format a genomic position (in bp) as a human-readable string.

    Two decimals, i.e. 10 kb resolution. e920c19 reduced this to one decimal,
    which is too coarse for the tightest zoom panels: a 200 kb window then reads
    "94.8 Mb" to "95.0 Mb", labelling a 200 kb span at 100 kb precision. Restored
    2026-08-28 so the panels read as they do in the published Figure 3.
    """
    if bp_pos >= 1_000_000:
        val = bp_pos / 1_000_000
        # Use decimal only when it adds information
        if val == int(val):
            return f"{int(val)} Mb"
        return f"{val:.2f} Mb"
    elif bp_pos >= 1_000:
        val = bp_pos / 1_000
        if val == int(val):
            return f"{int(val)} kb"
        return f"{val:.1f} kb"
    return f"{bp_pos} bp"


def draw_coordinate_labels(d, x, y_start, y_end, view_start, view_end,
                           text_color, font_size=9):
    """
    Draw genomic coordinate labels at the top-left and bottom-left of a
    zoom panel so the user can see the genomic span at a glance.

    Labels are placed just to the left of the panel at the top and bottom.
    """
    top_label = _format_genomic_pos(view_start)
    bot_label = _format_genomic_pos(view_end)

    d.append(draw.Text(top_label, font_size=font_size, x=x, y=y_start + 3,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='end', dominant_baseline='hanging'))
    d.append(draw.Text(bot_label, font_size=font_size, x=x, y=y_end - 3,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='end', dominant_baseline='auto'))


def draw_breakpoint_marker(d, region, x_start, x_end, y_start, pixels_per_bp,
                           view_start, label_side='right', text_color='white'):
    """
    Draw a yellow dashed horizontal line across all feature tracks at the
    region center position, with a coordinate label.

    Drawn for breakpoint regions (at center) and dicentric centromere
    regions (at midpoint of start/end).

    Args:
        label_side: 'right' to place label right of tracks,
                    'left' to place label left of tracks,
                    None to suppress label.
        text_color: Color for the position label text.
    """
    if region['region_type'] == 'breakpoint':
        bp_center = region['center']
    elif region['region_type'].startswith('dicentric_cen'):
        # Centromere center is midpoint of the refined start/end
        cen_start = region.get('start')
        cen_end = region.get('end')
        if cen_start is not None and cen_end is not None:
            bp_center = (cen_start + cen_end) // 2
        else:
            bp_center = region['center']
    else:
        return

    # Only draw if the center falls within the current view
    if bp_center < view_start:
        return

    y = y_start + (bp_center - view_start) * pixels_per_bp

    d.append(draw.Line(
        x_start, y, x_end, y,
        stroke='#000000',
        stroke_width=1.0,
        stroke_dasharray='4,2'
    ))

    # Breakpoint position label
    if label_side is not None:
        bp_label = _format_genomic_pos(bp_center)
        if label_side == 'right':
            d.append(draw.Text(bp_label, font_size=10, x=x_end + 4, y=y + 3,
                               fill='#000000', font_family='sans-serif',
                               text_anchor='start'))
        else:  # left
            d.append(draw.Text(bp_label, font_size=10, x=x_start - 4, y=y + 3,
                               fill='#000000', font_family='sans-serif',
                               text_anchor='end'))


def draw_region_indicator(d, region, x_start, x_end, y_start, pixels_per_bp,
                          view_start, text_color, bar_width, num_tracks, plot_height):
    """Draw indicator showing the region of interest."""
    total_width = x_end - x_start
    indicator_x = x_end + 10  # Right of tracks

    if region['region_type'].startswith('dicentric'):
        # Bracket around centromere region
        if region.get('start') and region.get('end'):
            cen_start_y = y_start + (region['start'] - view_start) * pixels_per_bp
            cen_end_y = y_start + (region['end'] - view_start) * pixels_per_bp

            # Bracket lines
            d.append(draw.Line(indicator_x, cen_start_y, indicator_x + 5, cen_start_y,
                               stroke=text_color, stroke_width=1))
            d.append(draw.Line(indicator_x + 5, cen_start_y, indicator_x + 5, cen_end_y,
                               stroke=text_color, stroke_width=1))
            d.append(draw.Line(indicator_x, cen_end_y, indicator_x + 5, cen_end_y,
                               stroke=text_color, stroke_width=1))

            # Label
            label = region['region_type'].replace('dicentric_', '').upper()
            d.append(draw.Text(label, font_size=12, x=indicator_x + 10,
                               y=(cen_start_y + cen_end_y) / 2 + 3,
                               fill=text_color, font_family='sans-serif'))

    elif region['region_type'] == 'breakpoint':
        # Get breakpoint interval bounds (or use ±500kb fallback)
        bp_interval_start = region.get('start')
        bp_interval_end = region.get('end')

        # Fallback to ±500kb if interval not defined
        if bp_interval_start is None or bp_interval_end is None:
            bp_center = region['center']
            bp_interval_start = max(0, bp_center - 500_000)
            bp_interval_end = bp_center + 500_000

        # Draw shaded interval region across all tracks
        interval_y_start = y_start + (bp_interval_start - view_start) * pixels_per_bp
        interval_y_end = y_start + (bp_interval_end - view_start) * pixels_per_bp
        interval_height = interval_y_end - interval_y_start

        # Only draw interval if meaningful: >1px and <80% of plot height
        # (At 500kb zoom, the ±500kb interval covers entire view - not useful)
        if interval_height > 1 and interval_height < plot_height * 0.8:
            d.append(draw.Rectangle(
                x_start - 5, interval_y_start,
                total_width + 10, interval_height,
                fill='#FF6666', fill_opacity=0.2
            ))


def draw_feature_track(d, features, colors, x, y_start, bar_width,
                       view_start, pixels_per_bp, background, track_height=None):
    """Draw vertical track of colored rectangles.

    If background is 'white' and track_height is provided, draws a single
    black outline around the entire track for visual separation.
    """
    for feat in features:
        y = y_start + (feat['start'] - view_start) * pixels_per_bp
        height = (feat['end'] - feat['start']) * pixels_per_bp

        color, opacity = colors.get(feat['feature'], ('#808080', 1.0))

        d.append(draw.Rectangle(x, y, bar_width, height,
                                fill=color, fill_opacity=opacity))

    # Outer track outline (white background only)
    if background == 'white' and track_height is not None:
        d.append(draw.Rectangle(x, y_start, bar_width, track_height,
                                fill='none', stroke='#000000',
                                stroke_width=0.5))


def draw_full_contig(d, features_by_track, colors, featuresets, scaffold_length,
                     x_start, y_start, total_height, feature_positions, args,
                     arm_features=None, regions=None,
                     show_coord_labels=True, show_breakpoint_markers=True,
                     show_scale_bar=True, show_track_labels=True):
    """
    Draw full contig view with 5 feature tracks.

    Args:
        d: Drawing object
        features_by_track: Dict of featureset -> features list
        colors: Feature colors
        featuresets: List of featureset names
        scaffold_length: Total scaffold length in bp
        x_start: X position for left edge
        y_start: Y position for top edge
        total_height: Height of the contig visualization
        feature_positions: List of (position, region_type) for markers
        args: Command line args
        arm_features: Optional list of arm features for p/q indicator line
        regions: Optional list of region dicts (for breakpoint lines + labels)

    Returns:
        x_end: X position of right edge of contig tracks
    """
    background = args.background
    text_color = 'white' if background == 'black' else 'black'
    bar_width = 10  # Match zoom track width
    bar_spacing = 2
    num_tracks = len(featuresets)

    # Check if arm indicator line should be shown
    show_arm_line = arm_features is not None and len(arm_features) > 0
    arm_line_width = 3  # Thin line width
    arm_line_gap = 2    # Gap between arm line and chromosome track

    pixels_per_bp = total_height / scaffold_length

    # Calculate x offset for tracks (shifted right if arm line shown)
    track_x_offset = (arm_line_width + arm_line_gap) if show_arm_line else 0

    # Draw feature tracks
    for i, fs in enumerate(featuresets):
        x = x_start + track_x_offset + i * (bar_width + bar_spacing)
        features = features_by_track.get(fs, [])
        fs_colors = colors.get(fs, {})
        draw_feature_track(d, features, fs_colors, x, y_start, bar_width,
                           0, pixels_per_bp, background,
                           track_height=total_height)

    # Draw p/q arm indicator line to the left of chromosome track
    if show_arm_line:
        arm_line_x = x_start
        arm_colors = {
            'p_arm_specific': '#A0A0A0',  # Light gray for p
            'q_arm_specific': '#505050',  # Dark gray for q
            'p_arm': '#A0A0A0',
            'q_arm': '#505050',
        }
        # Draw colored segments for arm annotations (line interrupted where no annotation)
        for feat in arm_features:
            # Handle both dict and tuple formats
            if isinstance(feat, dict):
                start, end, feature = feat['start'], feat['end'], feat['feature']
            else:
                start, end, feature = feat

            # Only draw p_arm or q_arm annotations
            if 'p_arm' not in feature and 'q_arm' not in feature:
                continue

            color = arm_colors.get(feature)
            if not color:
                continue

            y = y_start + start * pixels_per_bp
            height = (end - start) * pixels_per_bp
            d.append(draw.Rectangle(arm_line_x, y, arm_line_width, height, fill=color))

    total_track_width = track_x_offset + num_tracks * bar_width + (num_tracks - 1) * bar_spacing

    # Scale bar on left (shifted further left if arm line is shown)
    if show_scale_bar:
        scale_bar_offset = 10 + (arm_line_width + arm_line_gap if show_arm_line else 0)
        if scaffold_length >= 100_000_000:
            scale_length = 50_000_000
            scale_label = "50 Mb"
        elif scaffold_length >= 50_000_000:
            scale_length = 20_000_000
            scale_label = "20 Mb"
        elif scaffold_length >= 10_000_000:
            scale_length = 10_000_000
            scale_label = "10 Mb"
        else:
            scale_length = scaffold_length // 5
            scale_label = f"{scale_length // 1_000_000} Mb" if scale_length >= 1_000_000 else f"{scale_length // 1_000} kb"

        scale_bar_height = scale_length * pixels_per_bp
        scale_bar_x = x_start - scale_bar_offset
        scale_bar_y_top = y_start + 15
        d.append(draw.Rectangle(scale_bar_x, scale_bar_y_top, 2, scale_bar_height, fill=text_color))
        d.append(draw.Line(scale_bar_x - 2, scale_bar_y_top,
                           scale_bar_x + 4, scale_bar_y_top,
                           stroke=text_color, stroke_width=1))
        d.append(draw.Line(scale_bar_x - 2, scale_bar_y_top + scale_bar_height,
                           scale_bar_x + 4, scale_bar_y_top + scale_bar_height,
                           stroke=text_color, stroke_width=1))
        label_x = scale_bar_x - 5
        label_y = scale_bar_y_top + scale_bar_height / 2
        d.append(draw.Text(scale_label, font_size=9, x=label_x, y=label_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', dominant_baseline='middle',
                           transform=f"rotate(-90 {label_x} {label_y})"))

    # Calculate x_end for return value
    x_end = x_start + total_track_width

    # Track labels below (vertical) - matching KaryoScope_cluster_plot.py
    featureset_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomere',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }
    label_y = y_start + total_height + 6

    if show_track_labels:
        for i, fs in enumerate(featuresets):
            x = x_start + track_x_offset + i * (bar_width + bar_spacing) + bar_width / 2
            label = featureset_names.get(fs, fs[:6])
            d.append(draw.Text(label, font_size=10, x=x, y=label_y,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='end', dominant_baseline='middle',
                               transform=f"rotate(-90 {x} {label_y})"))

    # Add arm line label if shown (positioned at the arm line)
    if show_arm_line and show_track_labels:
        arm_label_x = x_start + arm_line_width / 2
        d.append(draw.Text("Arm", font_size=10, x=arm_label_x, y=label_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='end', dominant_baseline='middle',
                           transform=f"rotate(-90 {arm_label_x} {label_y})"))

    # Genomic coordinate labels at top-left and bottom-left
    arm_offset = (arm_line_width + arm_line_gap + 2) if show_arm_line else 0
    coord_label_x = x_start - 2 - arm_offset
    if show_coord_labels:
        draw_coordinate_labels(d, coord_label_x, y_start, y_start + total_height,
                               0, scaffold_length, text_color, font_size=9)

    # Breakpoint and centromere markers on the full contig.
    # Draw yellow dashed lines now; collect text labels to return for
    # deferred drawing (so they appear on top of fan polygons).
    deferred_labels = []
    if regions and show_breakpoint_markers:
        track_x_left = x_start + track_x_offset - 1
        track_x_right = x_start + total_track_width + 1

        # Zone occupied by scale bar + coordinate labels on the left side.
        # If a label falls in this zone, place it on the right instead.
        scale_bar_zone_top = y_start
        if show_scale_bar:
            scale_bar_zone_bottom = scale_bar_y_top + scale_bar_height + 10
        else:
            scale_bar_zone_bottom = y_start
        bottom_zone_top = y_start + total_height - 10

        for region in regions:
            rt = region.get('region_type', '')
            if rt == 'breakpoint':
                marker_center = region['center']
            elif rt.startswith('dicentric_cen'):
                cen_start = region.get('start')
                cen_end = region.get('end')
                if cen_start is not None and cen_end is not None:
                    marker_center = (cen_start + cen_end) // 2
                else:
                    marker_center = region['center']
            else:
                continue

            if 0 <= marker_center <= scaffold_length:
                marker_y = y_start + marker_center * pixels_per_bp
                # Yellow dashed line across all tracks
                d.append(draw.Line(
                    track_x_left, marker_y, track_x_right, marker_y,
                    stroke='#000000', stroke_width=1.0,
                    stroke_dasharray='4,2'
                ))
                # Position label — deferred so it draws on top of fans
                marker_label = _format_genomic_pos(marker_center)
                if marker_y < scale_bar_zone_bottom or marker_y > bottom_zone_top:
                    # Place on right side of tracks, tucked left to avoid
                    # overlapping the zoom panel's scale bar
                    deferred_labels.append((marker_label, track_x_right + 25,
                                            marker_y + 2, 'end', text_color))
                else:
                    # Place on left side, aligned with coordinate labels
                    deferred_labels.append((marker_label, coord_label_x,
                                            marker_y + 2, 'end', text_color))

    return x_end, deferred_labels


def draw_region_plot(region, features_by_track, colors, featuresets,
                     view_start, view_end, args, padding):
    """Generate SVG for one region at one zoom level."""
    background = args.background
    text_color = 'white' if background == 'black' else 'black'
    bar_width = args.bar_width
    bar_spacing = args.bar_spacing
    plot_height = args.plot_height

    num_tracks = len(featuresets)
    total_track_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing

    # Calculate pixels per bp
    view_span = view_end - view_start
    pixels_per_bp = plot_height / view_span

    # Layout dimensions
    margin_left = 80
    margin_top = 80
    margin_right = 100
    margin_bottom = 50

    image_width = margin_left + total_track_width + margin_right
    image_height = margin_top + plot_height + margin_bottom

    # Create drawing
    d = draw.Drawing(image_width, image_height)
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background))

    # Title
    title = region.get('iscn_notation', region['scaffold'])
    d.append(draw.Text(title, font_size=15, x=image_width/2, y=20,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle', font_weight='bold'))

    # Subtitle
    zoom_label = f"±{padding/1_000_000:.0f}Mb" if padding >= 1_000_000 else f"±{padding/1_000:.0f}kb"
    original_name = get_contig_name(region['scaffold'])
    subtitle = f"{original_name[:40]}... | {region['region_type']} | {zoom_label}"
    if len(original_name) <= 40:
        subtitle = f"{original_name} | {region['region_type']} | {zoom_label}"
    d.append(draw.Text(subtitle, font_size=12, x=image_width/2, y=35,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle'))

    # Track labels will be drawn below tracks (after position axis)

    # Draw feature tracks
    y_start = margin_top
    for i, fs in enumerate(featuresets):
        x = margin_left + i * (bar_width + bar_spacing)
        features = features_by_track.get(fs, [])
        fs_colors = colors.get(fs, {})
        draw_feature_track(d, features, fs_colors, x, y_start, bar_width,
                           view_start, pixels_per_bp, background,
                           track_height=plot_height)

    # Track borders
    for i in range(num_tracks + 1):
        x = margin_left + i * (bar_width + bar_spacing) - bar_spacing/2
        if i > 0:
            d.append(draw.Line(x, y_start, x, y_start + plot_height,
                               stroke=text_color, stroke_width=0.5, stroke_opacity=0.5))

    # Outer border
    x_start = margin_left - 1
    x_end = margin_left + total_track_width + 1
    d.append(draw.Rectangle(x_start, y_start, x_end - x_start, plot_height,
                            fill='none', stroke=text_color, stroke_width=1))

    # Vertical position axis on the left
    axis_x = margin_left - 10
    draw_vertical_position_axis(d, axis_x, y_start, y_start + plot_height,
                                view_start, view_end, pixels_per_bp, text_color, padding)

    # Vertical scale bar on the far left
    scale_length, _ = get_scale_bar_config(padding)
    scale_bar_x = 15
    scale_bar_y = y_start + (plot_height - scale_length * pixels_per_bp) / 2  # Centered vertically
    draw_vertical_scale_bar(d, scale_bar_x, scale_bar_y, scale_length, pixels_per_bp, text_color)

    # Vertical track labels below tracks
    label_y = y_start + plot_height + 10
    featureset_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomeric',
        'region': 'Region',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }
    for i, fs in enumerate(featuresets):
        x = margin_left + i * (bar_width + bar_spacing) + bar_width / 2
        label = featureset_names.get(fs, fs.title())
        d.append(draw.Text(label, font_size=10, x=x, y=label_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='end', dominant_baseline='middle',
                           transform=f"rotate(-90 {x} {label_y})"))

    # Genomic coordinate labels at top-left and bottom-left
    draw_coordinate_labels(d, x_start - 12, y_start, y_start + plot_height,
                           view_start, view_end, text_color)

    # Region indicator (shaded interval for breakpoints/centromeres)
    draw_region_indicator(d, region, x_start, x_end, y_start, pixels_per_bp,
                          view_start, text_color, bar_width, num_tracks, plot_height)

    # Breakpoint center marker (yellow dashed line across all tracks)
    draw_breakpoint_marker(d, region, x_start, x_end, y_start, pixels_per_bp,
                           view_start, label_side='right', text_color=text_color)

    return d


def draw_combined_zoom_plot(region, all_features_by_zoom, colors, featuresets, zoom_levels, args):
    """
    Generate SVG with all zoom levels side by side, connected by shaded fans.

    Args:
        region: Region dict with center, scaffold, etc.
        all_features_by_zoom: Dict of padding -> features_by_track
        colors: Feature colors
        featuresets: List of featureset names
        zoom_levels: List of padding values (sorted largest to smallest)
        args: Command line args
    """
    background = args.background
    text_color = 'white' if background == 'black' else 'black'
    bar_width = 12  # Narrower bars for combined view
    bar_spacing = 2
    plot_height = args.plot_height

    num_tracks = len(featuresets)
    num_zooms = len(zoom_levels)
    track_block_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing

    # Featureset labels
    featureset_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomeric',
        'region': 'Region',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }

    # Layout
    margin_left = 45  # Space for scale bars on left
    margin_top = 60
    margin_bottom = 55  # Space for track labels only
    gap_between_zooms = 30  # Space for shaded fans

    total_width = margin_left + num_zooms * track_block_width + (num_zooms - 1) * gap_between_zooms + 20
    total_height = margin_top + plot_height + margin_bottom

    # Create drawing
    d = draw.Drawing(total_width, total_height)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # Title
    title = region.get('iscn_notation', region['scaffold'])
    d.append(draw.Text(title, font_size=15, x=total_width/2, y=18,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle', font_weight='bold'))

    # Subtitle with scaffold name (use original contig name)
    original_name = get_contig_name(region['scaffold'])
    subtitle = original_name[:50] + ('...' if len(original_name) > 50 else '')
    d.append(draw.Text(subtitle, font_size=11, x=total_width/2, y=32,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle'))

    y_start = margin_top
    zoom_positions = []  # Store position info for each zoom level

    # Draw each zoom level
    for zoom_idx, padding in enumerate(sorted(zoom_levels, reverse=True)):  # Largest first
        view_start, view_end = compute_view_window(region['center'], padding)
        view_span = view_end - view_start
        pixels_per_bp = plot_height / view_span

        x_offset = margin_left + zoom_idx * (track_block_width + gap_between_zooms)

        features_by_track = all_features_by_zoom.get(padding, {})

        # Draw feature tracks
        for i, fs in enumerate(featuresets):
            x = x_offset + i * (bar_width + bar_spacing)
            features = features_by_track.get(fs, [])
            fs_colors = colors.get(fs, {})
            draw_feature_track(d, features, fs_colors, x, y_start, bar_width,
                               view_start, pixels_per_bp, background,
                               track_height=plot_height)

        # Outer border for this zoom block
        block_x_start = x_offset - 1
        block_x_end = x_offset + track_block_width + 1
        d.append(draw.Rectangle(block_x_start, y_start, block_x_end - block_x_start, plot_height,
                                fill='none', stroke=text_color, stroke_width=1))

        # Breakpoint center marker (yellow dashed line across all tracks)
        # (label drawn later, on top of fans)
        draw_breakpoint_marker(d, region, block_x_start, block_x_end,
                               y_start, pixels_per_bp, view_start,
                               label_side=None, text_color=text_color)

        # Store position info for fans
        zoom_positions.append({
            'x_start': block_x_start,
            'x_end': block_x_end,
            'x_offset': x_offset,
            'view_start': view_start,
            'view_end': view_end,
            'padding': padding,
            'pixels_per_bp': pixels_per_bp
        })

        # Zoom label at TOP of each subplot (show total window size)
        window_size = view_end - view_start
        if window_size >= 1_000_000:
            zoom_label = f"{window_size / 1_000_000:.0f} Mb"
        else:
            zoom_label = f"{window_size / 1_000:.0f} kb"
        d.append(draw.Text(zoom_label, font_size=13, x=(block_x_start + block_x_end)/2,
                           y=y_start - 6,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle'))

        # Vertical track labels BELOW each subplot
        label_y = y_start + plot_height + 8
        for i, fs in enumerate(featuresets):
            x = x_offset + i * (bar_width + bar_spacing) + bar_width / 2
            label = featureset_names.get(fs, fs.title())
            d.append(draw.Text(label, font_size=9, x=x, y=label_y,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='end', dominant_baseline='middle',
                               transform=f"rotate(-90 {x} {label_y})"))

        # Store scale bar info for drawing later (on top of fans)
        scale_length, _ = get_scale_bar_config(padding)
        scale_bar_height = scale_length * pixels_per_bp
        zoom_positions[-1]['scale_bar'] = {
            'x': block_x_start - 8,
            'y': y_start + 5,
            'height': scale_bar_height,
            'label': f"{scale_length // 1_000_000} Mb" if scale_length >= 1_000_000 else f"{scale_length // 1_000} kb"
        }

        # Store interval line info for drawing later (on top of fans)
        # Shows the region of interest: breakpoint interval or centromere extent
        interval_start = region.get('start')
        interval_end = region.get('end')
        if interval_start is not None and interval_end is not None:
            interval_bp = interval_end - interval_start

            # Skip error bar entirely when the interval is essentially zero
            if interval_bp < 1000:
                pass  # No CI line for point estimates
            else:
                if interval_bp >= 1_000_000:
                    interval_label = f"{interval_bp / 1_000_000:.1f} Mb"
                else:
                    interval_label = f"{interval_bp / 1_000:.0f} kb"

                interval_y_start = y_start + (interval_start - view_start) * pixels_per_bp
                interval_y_end = y_start + (interval_end - view_start) * pixels_per_bp
                interval_y_start = max(y_start, min(y_start + plot_height, interval_y_start))
                interval_y_end = max(y_start, min(y_start + plot_height, interval_y_end))

                if interval_y_end > interval_y_start and (interval_y_end - interval_y_start) < plot_height * 0.8:
                    zoom_positions[-1]['ci_line'] = {
                        'x': block_x_start - 4,
                        'y_start': interval_y_start,
                        'y_end': interval_y_end,
                        'label': interval_label
                    }

    # Draw shaded fans between zoom levels
    for i in range(len(zoom_positions) - 1):
        curr = zoom_positions[i]
        next_zoom = zoom_positions[i + 1]

        # Find where the next zoom's view falls within current zoom
        next_view_start = next_zoom['view_start']
        next_view_end = next_zoom['view_end']

        # Calculate Y positions in current zoom for next zoom's view bounds
        curr_y_for_next_start = y_start + (next_view_start - curr['view_start']) * curr['pixels_per_bp']
        curr_y_for_next_end = y_start + (next_view_end - curr['view_start']) * curr['pixels_per_bp']

        # Clamp to plot bounds
        curr_y_for_next_start = max(y_start, min(y_start + plot_height, curr_y_for_next_start))
        curr_y_for_next_end = max(y_start, min(y_start + plot_height, curr_y_for_next_end))

        # Draw shaded fan (polygon from narrow region to full height of next)
        fan_x1 = curr['x_end'] + 1  # Right edge of current
        fan_x2 = next_zoom['x_start'] - 1  # Left edge of next

        # Fan shape: narrow on left (current zoom), full height on right (next zoom)
        fan_points = [
            (fan_x1, curr_y_for_next_start),  # Top-left (narrow)
            (fan_x2, y_start),                 # Top-right (full height)
            (fan_x2, y_start + plot_height),   # Bottom-right (full height)
            (fan_x1, curr_y_for_next_end),     # Bottom-left (narrow)
        ]

        path_d = f"M {fan_points[0][0]},{fan_points[0][1]} "
        for px, py in fan_points[1:]:
            path_d += f"L {px},{py} "
        path_d += "Z"

        d.append(draw.Path(d=path_d, fill='#606060', fill_opacity=FAN_FILL_OPACITY,
                           stroke='#888888', stroke_width=0.5))

    # Draw scale bars, CI lines, coordinate labels, and breakpoint labels
    # ON TOP of fans
    for zp in zoom_positions:
        # Scale bar
        sb = zp.get('scale_bar')
        if sb and show_scale_bars:
            d.append(draw.Rectangle(sb['x'], sb['y'], 1.5, sb['height'], fill=text_color))
            d.append(draw.Line(sb['x'] - 1, sb['y'],
                               sb['x'] + 2.5, sb['y'],
                               stroke=text_color, stroke_width=1))
            d.append(draw.Line(sb['x'] - 1, sb['y'] + sb['height'],
                               sb['x'] + 2.5, sb['y'] + sb['height'],
                               stroke=text_color, stroke_width=1))
            label_y = sb['y'] + sb['height'] / 2
            d.append(draw.Text(sb['label'], font_size=9, x=sb['x'] - 2, y=label_y,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='end', dominant_baseline='middle'))

        # NOTE: not gated on show_scale_bars, unlike the stacked path's
        # equivalent. `show_scale_bars` is UNDEFINED in this function -- line
        # ~1623 above references it too and would raise NameError if a scale
        # bar were ever present here. Gating this line would turn that dormant
        # bug into a live one, so it is left alone; see the stacked path in
        # draw_scale_bars_and_ci, which has the parameter and is what the
        # Robertsonian figure uses.
        ci = zp.get('ci_line')
        if ci:
            d.append(draw.Line(ci['x'], ci['y_start'], ci['x'], ci['y_end'],
                               stroke='#000000', stroke_width=1.5, stroke_dasharray='4,2'))

        # Genomic coordinate labels (on top of fans)
        draw_coordinate_labels(d, zp['x_start'] - 2, y_start,
                               y_start + plot_height,
                               zp['view_start'], zp['view_end'],
                               text_color, font_size=9)

    return d


def draw_zoom_panel_inline(d, region, all_features_by_zoom, colors, featuresets,
                           zoom_levels, x_offset, y_offset, panel_height, args,
                           is_last_panel=False, scaffold_length=None,
                           show_breakpoint_lines=True, fan_stroke=True,
                           show_zoom_labels=True, gap_between_zooms=45):
    """
    Draw a multi-zoom panel inline within a larger drawing.

    This is a modified version of draw_combined_zoom_plot that draws into
    an existing drawing at a specified position.

    Args:
        is_last_panel: If True, draw track labels under all zoom levels

    Returns:
        dict with 'x_start', 'x_end', 'zoom_positions', 'text_color'
    """
    background = args.background
    text_color = 'white' if background == 'black' else 'black'
    bar_width = 10  # Narrower for integrated view
    bar_spacing = 2
    plot_height = panel_height

    num_tracks = len(featuresets)
    num_zooms = len(zoom_levels)
    track_block_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing

    # gap_between_zooms is a parameter: draw_stacked_plot derives it from the
    # coordinate labels, which hang into these gaps. It used to be hardcoded to 45
    # here, so widening the caller's copy moved the canvas edge but not the panels.

    # Labels matching KaryoScope_cluster_plot.py
    featureset_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomere',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }

    y_start = y_offset
    zoom_positions = []

    # Draw each zoom level
    for zoom_idx, padding in enumerate(sorted(zoom_levels, reverse=True)):
        view_start, view_end = compute_view_window(
            region['center'], padding, scaffold_length)
        view_span = view_end - view_start
        pixels_per_bp = plot_height / view_span

        block_x = x_offset + zoom_idx * (track_block_width + gap_between_zooms)

        features_by_track = all_features_by_zoom.get(padding, {})

        # Draw feature tracks
        for i, fs in enumerate(featuresets):
            x = block_x + i * (bar_width + bar_spacing)
            features = features_by_track.get(fs, [])
            fs_colors = colors.get(fs, {})
            draw_feature_track(d, features, fs_colors, x, y_start, bar_width,
                               view_start, pixels_per_bp, background,
                               track_height=plot_height)

        # Track block bounds (no border)
        block_x_start = block_x - 1
        block_x_end = block_x + track_block_width + 1

        # Breakpoint/centromere center marker (dashed line across all tracks)
        # (label and coordinate labels drawn later, on top of fans)
        if show_breakpoint_lines:
            draw_breakpoint_marker(d, region, block_x_start, block_x_end,
                                   y_start, pixels_per_bp, view_start,
                                   label_side=None, text_color=text_color)

        # Store position info
        zoom_positions.append({
            'x_start': block_x_start,
            'x_end': block_x_end,
            'x_offset': block_x,
            'view_start': view_start,
            'view_end': view_end,
            'padding': padding,
            'pixels_per_bp': pixels_per_bp
        })

        # Zoom label at top
        if show_zoom_labels:
            window_size = view_end - view_start
            if window_size >= 1_000_000:
                zoom_label = f"{window_size / 1_000_000:.0f} Mb"
            else:
                zoom_label = f"{window_size / 1_000:.0f} kb"
            d.append(draw.Text(zoom_label, font_size=13, x=(block_x_start + block_x_end)/2,
                               y=y_start - 6,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='middle'))

        # Store scale bar info (positioned 10% down from top to stay over fan area)
        scale_length, _ = get_scale_bar_config(padding)
        scale_bar_height = scale_length * pixels_per_bp
        scale_bar_y_offset = plot_height * 0.10  # 10% down from top
        zoom_positions[-1]['scale_bar'] = {
            'x': block_x_start - 6,
            'y': y_start + scale_bar_y_offset,
            'height': scale_bar_height,
            'label': f"{scale_length // 1_000_000} Mb" if scale_length >= 1_000_000 else f"{scale_length // 1_000} kb"
        }

        # Store interval line info (breakpoints and centromeres)
        interval_start = region.get('start')
        interval_end = region.get('end')
        if interval_start is not None and interval_end is not None:
            interval_bp = interval_end - interval_start

            # Skip error bar entirely when the interval is essentially zero
            if interval_bp >= 1000:
                if interval_bp >= 1_000_000:
                    interval_label = f"{interval_bp / 1_000_000:.1f} Mb"
                else:
                    interval_label = f"{interval_bp / 1_000:.0f} kb"

                interval_y_start = y_start + (interval_start - view_start) * pixels_per_bp
                interval_y_end = y_start + (interval_end - view_start) * pixels_per_bp
                interval_y_start = max(y_start, min(y_start + plot_height, interval_y_start))
                interval_y_end = max(y_start, min(y_start + plot_height, interval_y_end))

                if interval_y_end > interval_y_start and (interval_y_end - interval_y_start) < plot_height * 0.8:
                    zoom_positions[-1]['ci_line'] = {
                        'x': block_x_start - 6,  # Same x as scale bar
                        'y_start': interval_y_start,
                        'y_end': interval_y_end,
                        'label': interval_label
                    }

        # Store coordinate label info for drawing on top of fans
        zoom_positions[-1]['coord_labels'] = {
            'x': block_x_start - 2,
            'y_start': y_start,
            'y_end': y_start + plot_height,
            'view_start': view_start,
            'view_end': view_end,
        }

    # Draw shaded fans between zoom levels (with rectangular extensions)
    rect_width = 8  # Horizontal section width for inter-zoom fans
    for i in range(len(zoom_positions) - 1):
        curr = zoom_positions[i]
        next_zoom = zoom_positions[i + 1]

        next_view_start = next_zoom['view_start']
        next_view_end = next_zoom['view_end']

        curr_y_for_next_start = y_start + (next_view_start - curr['view_start']) * curr['pixels_per_bp']
        curr_y_for_next_end = y_start + (next_view_end - curr['view_start']) * curr['pixels_per_bp']

        curr_y_for_next_start = max(y_start, min(y_start + plot_height, curr_y_for_next_start))
        curr_y_for_next_end = max(y_start, min(y_start + plot_height, curr_y_for_next_end))

        fan_x1 = curr['x_end'] + 1
        fan_x2 = next_zoom['x_start'] - 1

        # 6-point polygon with rectangular extensions on both ends
        fan_points = [
            (fan_x1, curr_y_for_next_start),                      # Top-left
            (fan_x1 + rect_width, curr_y_for_next_start),         # Top of left rect
            (fan_x2 - rect_width, y_start),                       # Top of right rect
            (fan_x2, y_start),                                     # Top-right
            (fan_x2, y_start + plot_height),                       # Bottom-right
            (fan_x2 - rect_width, y_start + plot_height),         # Bottom of right rect
            (fan_x1 + rect_width, curr_y_for_next_end),           # Bottom of left rect
            (fan_x1, curr_y_for_next_end),                        # Bottom-left
        ]

        path_d = f"M {fan_points[0][0]},{fan_points[0][1]} "
        for px, py in fan_points[1:]:
            path_d += f"L {px},{py} "
        path_d += "Z"

        d.append(draw.Path(d=path_d, fill='#606060', fill_opacity=FAN_FILL_OPACITY,
                           stroke='#888888' if fan_stroke else 'none',
                           stroke_width=0.5 if fan_stroke else 0))

    # Scale bars and CI lines are returned for drawing later (on top of all fans)
    # They are stored in zoom_positions and drawn by the caller

    # Track labels under ALL zoom levels, but only for the bottom panel
    if is_last_panel:
        label_y = y_start + plot_height + 6
        for zoom_info in zoom_positions:
            for i, fs in enumerate(featuresets):
                x = zoom_info['x_offset'] + i * (bar_width + bar_spacing) + bar_width / 2
                label = featureset_names.get(fs, fs[:6])
                d.append(draw.Text(label, font_size=10, x=x, y=label_y,
                                   fill=text_color, font_family='sans-serif',
                                   text_anchor='end', dominant_baseline='middle',
                           transform=f"rotate(-90 {x} {label_y})"))

    return {
        'x_start': zoom_positions[0]['x_start'],
        'x_end': zoom_positions[-1]['x_end'],
        'zoom_positions': zoom_positions,
        'text_color': text_color
    }


def draw_scale_bars_and_ci(d, zoom_positions, text_color,
                           show_scale_bars=True):
    """Draw scale bars, CI lines, coordinate labels, and breakpoint labels on top of everything."""
    for zp in zoom_positions:
        sb = zp.get('scale_bar')
        if sb and show_scale_bars:
            # Scale bar with notches
            d.append(draw.Rectangle(sb['x'], sb['y'], 1.5, sb['height'], fill=text_color))
            d.append(draw.Line(sb['x'] - 1, sb['y'],
                               sb['x'] + 2.5, sb['y'],
                               stroke=text_color, stroke_width=1))
            d.append(draw.Line(sb['x'] - 1, sb['y'] + sb['height'],
                               sb['x'] + 2.5, sb['y'] + sb['height'],
                               stroke=text_color, stroke_width=1))
            # Vertical label (rotated)
            label_x = sb['x'] - 4
            label_y = sb['y'] + sb['height'] / 2
            d.append(draw.Text(sb['label'], font_size=9, x=label_x, y=label_y,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='middle', dominant_baseline='middle',
                               transform=f"rotate(-90 {label_x} {label_y})"))

        # Interval line: solid with notches and size label. Gated on
        # show_scale_bars -- see the note at the other draw site.
        ci = zp.get('ci_line')
        if ci and show_scale_bars:
            ci_color = '#000000'
            # Solid vertical line
            d.append(draw.Rectangle(ci['x'], ci['y_start'], 1.5, ci['y_end'] - ci['y_start'],
                                    fill=ci_color))
            # Top notch
            d.append(draw.Line(ci['x'] - 1, ci['y_start'],
                               ci['x'] + 2.5, ci['y_start'],
                               stroke=ci_color, stroke_width=1))
            # Bottom notch
            d.append(draw.Line(ci['x'] - 1, ci['y_end'],
                               ci['x'] + 2.5, ci['y_end'],
                               stroke=ci_color, stroke_width=1))
            # Vertical size label
            ci_label_x = ci['x'] - 4
            ci_label_y = (ci['y_start'] + ci['y_end']) / 2
            ci_label_text = ci.get('label', '')
            d.append(draw.Text(ci_label_text, font_size=9, x=ci_label_x, y=ci_label_y,
                               fill=ci_color, font_family='sans-serif',
                               text_anchor='middle', dominant_baseline='middle',
                               transform=f"rotate(-90 {ci_label_x} {ci_label_y})"))

        # Genomic coordinate labels at top-left and bottom-left
        cl = zp.get('coord_labels')
        if cl:
            draw_coordinate_labels(d, cl['x'], cl['y_start'], cl['y_end'],
                                   cl['view_start'], cl['view_end'],
                                   text_color, font_size=9)


def _find_common_label(feature_names):
    """
    Find a common/generic label for multiple features that share a color.

    Examples:
        - ['multigroup1', 'multigroup2', 'multigroup3'] → 'multigroup'
        - ['centromeric_multigroup1', 'hsat_multigroup1'] → 'multigroup'
        - ['p_arm_specific', 'q_arm_specific'] → 'arm'
        - ['chr1_specific', 'chr2_specific'] → 'chromosome'
        - ['autosome_multigroup1', 'categorized'] → 'autosome multigroup, categorized' (condensed)
    """
    import re

    if not feature_names:
        return 'unknown'

    # Clean up names: remove _specific suffix and replace _ with space
    def clean_name(f):
        name = f.replace('_specific', '').replace('_', ' ')
        # Strip trailing numbers from multigroup names
        if 'multigroup' in name.lower():
            name = re.sub(r'\s*\d+$', '', name).rstrip()
        return name

    cleaned = [clean_name(f) for f in feature_names]

    # Remove duplicates while preserving order
    seen = set()
    unique_cleaned = []
    for name in cleaned:
        if name not in seen:
            seen.add(name)
            unique_cleaned.append(name)

    # Check if all unique features are just "multigroup" variants
    if all('multigroup' in name.lower() for name in unique_cleaned):
        return 'multigroup'

    # Check if all features contain "arm" - use generic "arm"
    if all('arm' in name.lower() for name in unique_cleaned):
        return 'arm'

    # If only one unique name after cleaning, return it
    if len(unique_cleaned) == 1:
        return unique_cleaned[0]

    # Check for numbered suffixes with same base
    base_names = set()
    for name in unique_cleaned:
        base = re.sub(r'\d+$', '', name).strip()
        base_names.add(base)

    if len(base_names) == 1:
        return base_names.pop()

    # Find longest common prefix
    if unique_cleaned:
        prefix = unique_cleaned[0]
        for name in unique_cleaned[1:]:
            while not name.startswith(prefix) and prefix:
                prefix = prefix[:-1]
        if len(prefix) >= 3:  # Minimum meaningful prefix
            return prefix.rstrip()

    # Fallback: use first cleaned name
    return unique_cleaned[0] if unique_cleaned else 'unknown'


def _estimate_legend_width(colors, featuresets, displayed_features, show_arm_line):
    """Estimate the total width needed for the legend (border included).

    Returns the required margin_right value (legend content + padding + gap).
    """
    import re as _re

    featureset_display_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomere',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }

    _DISPLAY_NAME_OVERRIDES = {'noncentromeric_specific': 'rDNA'}

    max_label_len = 0
    if show_arm_line:
        max_label_len = max(max_label_len, len("p arm"), len("q arm"))

    for fs in featuresets:
        fs_colors = colors.get(fs, {})
        fs_displayed = displayed_features.get(fs, set())
        if not fs_colors or not fs_displayed:
            continue
        track_name = featureset_display_names.get(fs, fs.title())
        max_label_len = max(max_label_len, len(track_name))

        # Build display names using same logic as _draw_legend
        color_to_features = {}
        for feature_name in fs_displayed:
            color_info = fs_colors.get(feature_name)
            if color_info is None:
                color_hex = '#808080'
            elif isinstance(color_info, tuple):
                color_hex = color_info[0]
            else:
                color_hex = color_info
            if color_hex not in color_to_features:
                color_to_features[color_hex] = []
            color_to_features[color_hex].append(feature_name)

        def _is_multigroup_or_categorized(f):
            return 'multigroup' in f or f == 'categorized'
        def _is_novel(f):
            return f == 'novel'
        def _is_chr_specific(f):
            return _re.match(r'^chr[0-9XY]+_specific$', f) is not None
        def _is_non_category(f):
            if f in ('arm_multigroup1', 'p_arm_specific', 'q_arm_specific'):
                return True
            if f == 'noncentromeric_specific':
                return False
            if f.startswith('non') and 'canonical' not in f:
                return True
            return False

        for color_hex, feature_names in color_to_features.items():
            specifics = [f for f in feature_names
                         if f.endswith('_specific') and not _is_non_category(f)]
            multigroups = [f for f in feature_names
                           if _is_multigroup_or_categorized(f) and not _is_non_category(f)]
            non_cats = [f for f in feature_names if _is_non_category(f)]
            novels = [f for f in feature_names if _is_novel(f)]

            if novels:
                display_name = 'novel'
            elif non_cats and not specifics and not multigroups:
                cleaned = []
                for f in non_cats:
                    name = f.replace('_specific', '').replace('_multigroup1', '').replace('_', ' ')
                    if name not in cleaned:
                        cleaned.append(name)
                display_name = cleaned[0] if cleaned else 'unknown'
            elif specifics and multigroups:
                if len(specifics) == 1:
                    display_name = specifics[0].replace('_specific', '').replace('_', ' ')
                else:
                    # _draw_legend collapses to the shared multigroup root when
                    # there is exactly one; mirror that or we under-measure again.
                    mg_roots = set()
                    for mg in multigroups:
                        root = _re.sub(r'\s*multigroup\s*\d*', '', mg.replace('_', ' ')).strip()
                        if root:
                            mg_roots.add(root)
                    if len(mg_roots) == 1:
                        display_name = mg_roots.pop()
                    else:
                        display_name = ', '.join(sorted(
                            f.replace('_specific', '').replace('_', ' ') for f in specifics))
            elif multigroups:
                # Must mirror _draw_legend exactly. This branch used to return the
                # literal 'multigroup' (10 chars) while the drawing code renders a
                # single multigroup as its full root -- 'centromeric_multigroup1'
                # becomes "centromeric multigroup" (22). Estimating 12 characters
                # short put that label outside the legend border in the raster.
                if len(multigroups) == 1:
                    display_name = multigroups[0].replace('_', ' ').replace(' specific', '')
                    display_name = _re.sub(r'\d+$', '', display_name).rstrip()
                else:
                    display_name = 'multigroup'
            elif specifics:
                if len(specifics) == 1:
                    display_name = specifics[0].replace('_specific', '').replace('_', ' ')
                else:
                    display_name = ', '.join(sorted(
                        f.replace('_specific', '').replace('_', ' ') for f in specifics))
            else:
                display_name = feature_names[0].replace('_', ' ')

            if len(feature_names) == 1 and feature_names[0] in _DISPLAY_NAME_OVERRIDES:
                display_name = _DISPLAY_NAME_OVERRIDES[feature_names[0]]

            max_label_len = max(max_label_len, len(display_name))

    swatch_size = 8
    char_width = 5.5
    padding = 5
    legend_content_width = swatch_size + 6 + max_label_len * char_width + padding
    # margin_right = gap from zoom panels to legend_x (10) + legend content + padding + outer margin
    return int(legend_content_width + 2 * padding + 10 + 15)


def _draw_legend(d, legend_x, legend_y, colors, featuresets, displayed_features,
                 show_arm_line, arm_features, text_color, margin_right=None,
                 draw_border=True):
    """Draw the color legend into an existing drawing at (legend_x, legend_y).

    Shared by draw_integrated_plot and draw_stacked_plot.
    Returns the bottom Y of the drawn legend.
    """
    featureset_display_names = {
        'chromosome': 'Chromosome',
        'subtelomeric': 'Subtelomere',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
        'acrocentric': 'Acrocentric',
    }

    item_height = 11
    swatch_size = 8
    section_gap = 15
    current_y = legend_y
    first_section = True
    max_label_len = 0  # Track longest label for dynamic border width

    # Arm indicator legend
    if show_arm_line and arm_features:
        max_label_len = max(max_label_len, len("p arm"), len("q arm"))
        d.append(draw.Text("Arm", font_size=10, x=legend_x, y=current_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='start', font_weight='bold'))
        current_y += 4
        d.append(draw.Rectangle(legend_x + 3, current_y, swatch_size, swatch_size, fill='#A0A0A0', stroke=text_color, stroke_width=0.5))
        d.append(draw.Text("p arm", font_size=9, x=legend_x + swatch_size + 6,
                           y=current_y + swatch_size / 2 + 1,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='start', dominant_baseline='middle'))
        current_y += item_height
        d.append(draw.Rectangle(legend_x + 3, current_y, swatch_size, swatch_size, fill='#505050', stroke=text_color, stroke_width=0.5))
        d.append(draw.Text("q arm", font_size=9, x=legend_x + swatch_size + 6,
                           y=current_y + swatch_size / 2 + 1,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='start', dominant_baseline='middle'))
        current_y += item_height
        first_section = False

    import re as _re

    def _is_multigroup_or_categorized(f):
        return 'multigroup' in f or f == 'categorized'

    def _is_novel(f):
        return f == 'novel'

    _DISPLAY_NAME_OVERRIDES = {'noncentromeric_specific': 'rDNA'}

    def _is_non_category(f):
        if f in ('arm_multigroup1', 'p_arm_specific', 'q_arm_specific'):
            return True
        if f == 'noncentromeric_specific':
            return False
        if f.startswith('non') and 'canonical' not in f:
            return True
        return False

    def _chrom_sort_num(f):
        m = _re.match(r'chr(\d+)', f)
        if m:
            return int(m.group(1))
        if f.startswith('chrX'):
            return 23
        if f.startswith('chrY'):
            return 24
        return 99

    def _is_chr_specific(f):
        return _re.match(r'^chr[0-9XY]+_specific$', f) is not None

    for fs in featuresets:
        fs_colors = colors.get(fs, {})
        fs_displayed = displayed_features.get(fs, set())
        if not fs_colors or not fs_displayed:
            continue

        if not first_section:
            current_y += section_gap
        first_section = False

        track_name = featureset_display_names.get(fs, fs.title())
        max_label_len = max(max_label_len, len(track_name))
        d.append(draw.Text(track_name, font_size=10, x=legend_x, y=current_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='start', font_weight='bold'))
        current_y += 4

        color_to_features = {}
        for feature_name in fs_displayed:
            color_info = fs_colors.get(feature_name)
            if color_info is None:
                color_hex = '#808080'
            elif isinstance(color_info, tuple):
                color_hex = color_info[0]
            else:
                color_hex = color_info
            if color_hex not in color_to_features:
                color_to_features[color_hex] = []
            color_to_features[color_hex].append(feature_name)

        legend_entries = []
        for color_hex, feature_names in color_to_features.items():
            specifics = [f for f in feature_names
                         if f.endswith('_specific') and not _is_non_category(f)]
            multigroups = [f for f in feature_names
                           if _is_multigroup_or_categorized(f) and not _is_non_category(f)]
            non_cats = [f for f in feature_names if _is_non_category(f)]
            novels = [f for f in feature_names if _is_novel(f)]
            chr_specifics = [f for f in specifics if _is_chr_specific(f)]

            if novels:
                display_name = 'novel'
            elif non_cats and not specifics and not multigroups:
                cleaned = []
                for f in non_cats:
                    name = f.replace('_specific', '').replace('_multigroup1', '').replace('_', ' ')
                    if name not in cleaned:
                        cleaned.append(name)
                if len(cleaned) > 1:
                    words = set()
                    for name in cleaned:
                        words.update(name.split())
                    common = [w for w in words if all(w in n.split() for n in cleaned)]
                    display_name = ' '.join(common) if common else cleaned[0]
                else:
                    display_name = cleaned[0] if cleaned else 'unknown'
            elif specifics and multigroups:
                if len(specifics) == 1:
                    display_name = specifics[0].replace('_specific', '').replace('_', ' ')
                else:
                    mg_roots = set()
                    for mg in multigroups:
                        root = mg.replace('_', ' ')
                        root = _re.sub(r'\s*multigroup\s*\d*', '', root).strip()
                        if root:
                            mg_roots.add(root)
                    if len(mg_roots) == 1:
                        display_name = mg_roots.pop()
                    else:
                        display_name = ', '.join(sorted(
                            f.replace('_specific', '').replace('_', ' ')
                            for f in specifics
                        ))
            elif multigroups:
                if len(multigroups) == 1:
                    display_name = multigroups[0].replace('_', ' ').replace(' specific', '')
                    display_name = _re.sub(r'\d+$', '', display_name).rstrip()
                else:
                    display_name = 'multigroup'
            elif specifics:
                if len(specifics) == 1:
                    display_name = specifics[0].replace('_specific', '').replace('_', ' ')
                else:
                    display_name = ', '.join(sorted(
                        f.replace('_specific', '').replace('_', ' ')
                        for f in specifics
                    ))
            else:
                display_name = feature_names[0].replace('_', ' ')

            if len(feature_names) == 1 and feature_names[0] in _DISPLAY_NAME_OVERRIDES:
                display_name = _DISPLAY_NAME_OVERRIDES[feature_names[0]]

            if novels:
                sort_key = (4, 0, display_name)
            elif non_cats and not specifics:
                sort_key = (3, 0, display_name)
            elif multigroups and not specifics:
                sort_key = (2, 0, display_name)
            elif chr_specifics:
                sort_key = (0, _chrom_sort_num(chr_specifics[0]), display_name)
            else:
                sort_key = (1, 0, display_name)

            legend_entries.append((sort_key, color_hex, display_name))

        legend_entries.sort()
        for sort_key, color_hex, display_name in legend_entries:
            max_label_len = max(max_label_len, len(display_name))
            d.append(draw.Rectangle(legend_x + 3, current_y, swatch_size, swatch_size, fill=color_hex, stroke=text_color, stroke_width=0.5))
            d.append(draw.Text(display_name, font_size=9, x=legend_x + swatch_size + 6,
                               y=current_y + swatch_size / 2 + 1,
                               fill=text_color, font_family='sans-serif',
                               text_anchor='start', dominant_baseline='middle'))
            current_y += item_height

    # Draw border around the legend. The published Figure 3 has none -- its SVG
    # contains no rect wider than the 10 px track outlines -- so the stacked path
    # passes draw_border=False.
    padding = 5
    if not draw_border:
        return current_y

    if margin_right is not None:
        # Use margin_right to ensure border fits within the plot area
        # legend_x starts at total_width - margin_right + 10
        # so available width from legend_x to right edge is margin_right - 10
        legend_width = margin_right - 10 - 2 * padding - 5
    else:
        char_width = 5.5
        legend_width = swatch_size + 6 + max_label_len * char_width + padding
    legend_height = current_y - legend_y + padding
    d.append(draw.Rectangle(legend_x - padding, legend_y - padding - 4,
                            legend_width + 2 * padding, legend_height + padding,
                            fill='none', stroke=text_color, stroke_width=0.5,
                            rx=3, ry=3))

    return current_y


def draw_integrated_plot(scaffold, regions, scaffold_length, contig_features,
                         all_region_features, colors, featuresets, zoom_levels, args,
                         kromatid_mapping=None, arm_features=None):
    """
    Generate integrated visualization: full contig + stacked multi-zoom panels.

    Args:
        scaffold: Scaffold name
        regions: List of regions on this scaffold (sorted by position)
        scaffold_length: Total scaffold length in bp
        contig_features: Dict of featureset -> features for full contig
        all_region_features: Dict of region_idx -> padding -> features_by_track
        colors: Feature colors
        featuresets: List of featureset names
        zoom_levels: List of zoom paddings
        args: Command line args
        kromatid_mapping: Optional dict mapping KaryoScope notation to KromaTid notation
        arm_features: Optional list of arm features for p/q track (contig view only)
    """
    background = args.background
    text_color = 'white' if background == 'black' else 'black'

    num_features = len(regions)

    # Check if arm indicator line should be shown
    show_arm_line = arm_features is not None and len(arm_features) > 0

    # Pre-compute displayed features for legend width estimation
    displayed_features = {fs: set() for fs in featuresets}
    for fs, features in contig_features.items():
        for f in features:
            if isinstance(f, dict):
                displayed_features[fs].add(f['feature'])
            else:
                displayed_features[fs].add(f[2])
    for region_idx, zoom_data in all_region_features.items():
        for padding_val, features_by_track in zoom_data.items():
            for fs, features in features_by_track.items():
                for f in features:
                    if isinstance(f, dict):
                        displayed_features[fs].add(f['feature'])
                    else:
                        displayed_features[fs].add(f[2])

    # Layout constants
    margin_left = 50
    margin_top = 60 if kromatid_mapping else 50  # Extra space for dual titles
    margin_bottom = 80  # Extra space for track labels
    margin_right = _estimate_legend_width(colors, featuresets, displayed_features, show_arm_line)

    # 5 tracks × 10px + spacing (matching zoom tracks), plus thin arm line if shown (3px + 2px gap)
    num_tracks = len(featuresets)
    contig_width = num_tracks * 10 + (num_tracks - 1) * 1 + (5 if show_arm_line else 0)
    fan_gap = 45  # Gap for contig→zoom fan
    contig_height = 350  # Fixed contig height regardless of number of regions
    panel_gap = 15 if num_features > 1 else 40  # Tighter gap for multi-region

    # Zoom panels share the fixed contig height
    if num_features <= 1:
        zoom_panel_height = contig_height
    else:
        zoom_panel_height = (contig_height - (num_features - 1) * panel_gap) // num_features
    total_panel_height = contig_height

    # Estimate legend height (items per track × item_height + section gaps)
    # Roughly: 5 tracks × (header + 8 items) × 11px + gaps
    estimated_legend_height = len(featuresets) * (1 + 8) * 11 + len(featuresets) * 8

    # Calculate zoom panel width (estimate)
    num_tracks = len(featuresets)
    bar_width = 10
    bar_spacing = 2
    gap_between_zooms = 45  # Match contig-to-zoom fan gap
    track_block_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing
    num_zooms = len(zoom_levels)
    zoom_panel_width = num_zooms * track_block_width + (num_zooms - 1) * gap_between_zooms + 20

    total_width = margin_left + contig_width + fan_gap + zoom_panel_width + margin_right
    # Ensure height fits both panels and legend
    min_height_for_legend = margin_top + estimated_legend_height + 20
    total_height = max(margin_top + contig_height + margin_bottom, min_height_for_legend)

    # Create drawing
    d = draw.Drawing(total_width, total_height)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # Title - ISCN notation from first region
    iscn = regions[0].get('iscn_notation', scaffold)

    # Check for KromaTid alternate naming
    kromatid_name = None
    if kromatid_mapping:
        kromatid_name = kromatid_mapping.get(iscn)

    # Format scaffold size in Mb (1 decimal)
    scaffold_size_mb = scaffold_length / 1_000_000

    if kromatid_name:
        # Two-line title: KaryoScope on top, KromaTid below
        d.append(draw.Text(f"KaryoScope: {iscn}", font_size=14, x=total_width/2, y=12,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', font_weight='bold'))
        d.append(draw.Text(f"KromaTid: {kromatid_name}", font_size=14, x=total_width/2, y=26,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', font_weight='bold'))
        # Subtitle - scaffold name with size (use original contig name)
        original_name = get_contig_name(scaffold)
        subtitle = original_name[:50] + ('...' if len(original_name) > 50 else '')
        d.append(draw.Text(f"{subtitle} ({scaffold_size_mb:.1f} Mb)", font_size=11, x=total_width/2, y=40,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle'))
    else:
        # Single title (no KromaTid mapping)
        d.append(draw.Text(iscn, font_size=15, x=total_width/2, y=12,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', font_weight='bold'))
        # Subtitle - scaffold name with size (use original contig name)
        original_name = get_contig_name(scaffold)
        subtitle = original_name[:50] + ('...' if len(original_name) > 50 else '')
        d.append(draw.Text(f"{subtitle} ({scaffold_size_mb:.1f} Mb)", font_size=11, x=total_width/2, y=26,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle'))

    # Draw full contig
    contig_x = margin_left
    contig_y = margin_top
    feature_positions = [(r['center'], r['region_type']) for r in regions]

    # Add ISCN classification label above contig view (unbolded)
    contig_center_x = contig_x + contig_width / 2
    d.append(draw.Text(iscn, font_size=13, x=contig_center_x, y=contig_y - 6,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle'))

    contig_x_end, contig_deferred_labels = draw_full_contig(
        d, contig_features, colors, featuresets, scaffold_length,
        contig_x, contig_y, contig_height, feature_positions, args,
        arm_features=arm_features, regions=regions,
        show_scale_bar=not args.hide_scale_bars,
    )

    # Draw stacked zoom panels and connecting fans
    zoom_panel_x = contig_x_end + fan_gap

    # Collect all panel info for drawing scale bars/CI at the end
    all_panel_infos = []

    for i, region in enumerate(regions):
        panel_y = margin_top + i * (zoom_panel_height + panel_gap)

        # Get features for this region's zoom levels
        region_features = all_region_features.get(i, {})

        # Draw zoom panel (track labels only on bottom panel)
        is_last_panel = (i == len(regions) - 1)
        panel_info = draw_zoom_panel_inline(
            d, region, region_features, colors, featuresets,
            zoom_levels, zoom_panel_x, panel_y, zoom_panel_height, args,
            is_last_panel=is_last_panel, scaffold_length=scaffold_length,
            fan_stroke=False,
            show_zoom_labels=(i == 0),
        )
        all_panel_infos.append(panel_info)

        # Draw connecting fan from contig to first zoom
        # First zoom's view bounds mapped to contig
        first_zoom = panel_info['zoom_positions'][0]
        first_zoom_view_start = first_zoom['view_start']
        first_zoom_view_end = first_zoom['view_end']

        contig_y_for_zoom_start = contig_y + (first_zoom_view_start / scaffold_length) * contig_height
        contig_y_for_zoom_end = contig_y + (first_zoom_view_end / scaffold_length) * contig_height

        # Clamp to contig bounds
        contig_y_for_zoom_start = max(contig_y, min(contig_y + contig_height, contig_y_for_zoom_start))
        contig_y_for_zoom_end = max(contig_y, min(contig_y + contig_height, contig_y_for_zoom_end))

        # Fan with rectangular extensions on BOTH sides
        rect_width = 12  # Horizontal section width
        fan_x1 = contig_x_end  # At contig edge
        fan_x2 = panel_info['x_start'] - 1  # Before first zoom

        # 8-point polygon: rect on left (contig) + fan + rect on right (zoom)
        fan_points = [
            (fan_x1, contig_y_for_zoom_start),                    # Top-left at contig
            (fan_x1 + rect_width, contig_y_for_zoom_start),       # Top of left rect
            (fan_x2 - rect_width, panel_y),                       # Top of right rect
            (fan_x2, panel_y),                                     # Top-right at zoom
            (fan_x2, panel_y + zoom_panel_height),                 # Bottom-right at zoom
            (fan_x2 - rect_width, panel_y + zoom_panel_height),   # Bottom of right rect
            (fan_x1 + rect_width, contig_y_for_zoom_end),         # Bottom of left rect
            (fan_x1, contig_y_for_zoom_end),                      # Bottom-left at contig
        ]

        path_d = f"M {fan_points[0][0]},{fan_points[0][1]} "
        for px, py in fan_points[1:]:
            path_d += f"L {px},{py} "
        path_d += "Z"

        d.append(draw.Path(d=path_d, fill='#606060', fill_opacity=FAN_FILL_OPACITY,
                           stroke='none'))

    # Draw all scale bars and CI lines on top of everything (after all fans)
    for panel_info in all_panel_infos:
        draw_scale_bars_and_ci(d, panel_info['zoom_positions'], panel_info['text_color'],
                               show_scale_bars=not args.hide_scale_bars)

    # Draw deferred contig marker labels on top of fans
    for label_text, lx, ly, anchor, lcolor in contig_deferred_labels:
        d.append(draw.Text(label_text, font_size=9,
                           x=lx, y=ly,
                           fill=lcolor, font_family='sans-serif',
                           text_anchor=anchor))

    # Draw vertical color legend on right side - one section per track
    legend_x = total_width - margin_right + 10
    legend_y = margin_top
    _draw_legend(d, legend_x, legend_y, colors, featuresets,
                 displayed_features, show_arm_line, arm_features, text_color,
                 margin_right=margin_right)

    return d


def format_zoom_label(padding):
    """Format zoom level for filename."""
    if padding >= 1_000_000:
        return f"{padding // 1_000_000}Mb"
    else:
        return f"{padding // 1_000}kb"


def parse_size(size_str: str) -> int:
    """Parse human-readable size string (e.g., '20Mb', '100kb', '1000') to int."""
    size_str = size_str.strip().lower()
    if size_str.endswith('mb'):
        return int(float(size_str[:-2]) * 1_000_000)
    elif size_str.endswith('kb'):
        return int(float(size_str[:-2]) * 1_000)
    elif size_str.endswith('bp'):
        return int(size_str[:-2])
    else:
        return int(size_str)


def draw_stacked_plot(all_scaffolds, colors, featuresets, zoom_levels, args):
    """Render multiple scaffolds stacked vertically with one shared legend.

    all_scaffolds: list of dicts returned by _process_sample, each containing
        scaffold, regions, scaffold_length, contig_features,
        all_region_features, arm_features, kromatid_mapping, sample.

    Returns a drawsvg.Drawing.
    """
    background = args.background
    text_color = 'white' if background == 'black' else 'black'

    # Per-scaffold layout constants (same as draw_integrated_plot)
    margin_left = 50
    top_margin_per_scaffold = 15
    # The bottom margin is the rotated track-label band. Only the last row draws
    # labels (see is_last_panel below), so the other rows do not reserve it.
    bottom_margin_per_scaffold = 80
    bottom_margin_no_labels = 20
    inter_scaffold_gap = 0
    margin_bottom_figure = 40

    # Every feature that will be drawn, so the legend can be sized from its real
    # content. The fixed `margin_right = 120` this replaces was too narrow for the
    # longest label ("centromeric multigroup"), which ran through the legend border
    # and off the canvas. The integrated path already sized itself this way; only
    # the stacked path was hardcoded.
    displayed_features = {fs: set() for fs in featuresets}
    for _sd in all_scaffolds:
        for _fs, _feats in _sd['contig_features'].items():
            for _f in _feats:
                displayed_features[_fs].add(_f['feature'] if isinstance(_f, dict) else _f[2])
        for _zoom_data in _sd['all_region_features'].values():
            for _features_by_track in _zoom_data.values():
                for _fs, _feats in _features_by_track.items():
                    for _f in _feats:
                        displayed_features[_fs].add(_f['feature'] if isinstance(_f, dict) else _f[2])

    num_tracks = len(featuresets)
    bar_width = 10
    bar_spacing = 2
    zoom_panel_height = 200
    panel_gap = 40

    # Both inter-panel gaps have to clear the coordinate labels, which are drawn
    # text-anchor="end" just left of each zoom panel and so hang into the gap on
    # their left. `fan_gap` is the contig-to-first-panel gap and `gap_between_zooms`
    # the rest; at two decimals the widest ("100.57 Mb") measures 47.8 px in DejaVu
    # Sans at font-size 9, which overflowed the old 45 on both. Derive it from the
    # labels that will actually be drawn so a future change to _format_genomic_pos
    # cannot silently reintroduce the overlap.
    _coord_chars = 0
    for _sd in all_scaffolds:
        for _region in _sd['regions']:
            for _padding in zoom_levels:
                _vs, _ve = compute_view_window(
                    _region['center'], _padding, _sd['scaffold_length'])
                _coord_chars = max(_coord_chars,
                                   len(_format_genomic_pos(_vs)),
                                   len(_format_genomic_pos(_ve)))
    # 2 px label offset + text + 6 px clearance, at the estimator's 5.5 px/char
    # (measured 5.31 for digit-heavy labels, so this is the conservative side).
    _coord_gap = int(2 + _coord_chars * 5.5 + 6) + 1
    gap_between_zooms = max(45, _coord_gap)
    fan_gap = max(45, _coord_gap)

    track_block_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing
    num_zooms = len(zoom_levels)
    zoom_panel_width = num_zooms * track_block_width + (num_zooms - 1) * gap_between_zooms + 20

    # Check if any scaffold needs the arm track
    any_arm = any(
        sd['arm_features'] is not None and len(sd['arm_features']) > 0
        for sd in all_scaffolds
    )
    contig_width = num_tracks * 10 + (num_tracks - 1) * 1 + (5 if any_arm else 0)

    margin_right = _estimate_legend_width(colors, featuresets, displayed_features, any_arm)

    # Compute per-scaffold heights based on number of regions (zoom panels)
    scaffold_heights = []
    for _h_idx, sd in enumerate(all_scaffolds):
        num_features = len(sd['regions'])
        if num_features == 0:
            num_features = 1
        total_panel_h = num_features * zoom_panel_height + (num_features - 1) * panel_gap
        _bm = (bottom_margin_per_scaffold
               if _h_idx == len(all_scaffolds) - 1 else bottom_margin_no_labels)
        # top_margin_per_scaffold for titles, then contig, then bottom margin
        body_h = top_margin_per_scaffold + total_panel_h + _bm
        scaffold_heights.append(body_h)

    total_body_height = sum(scaffold_heights) + inter_scaffold_gap * (len(all_scaffolds) - 1)

    # Legend: estimate vertical space
    estimated_legend_height = len(featuresets) * (1 + 8) * 11 + len(featuresets) * 8 + 30

    # Figure dimensions
    total_width = margin_left + contig_width + fan_gap + zoom_panel_width + margin_right
    total_height = max(total_body_height + margin_bottom_figure,
                       top_margin_per_scaffold + estimated_legend_height + 20)

    d = draw.Drawing(total_width, total_height)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # Accumulate displayed features across all scaffolds (for shared legend)
    merged_displayed = {fs: set() for fs in featuresets}

    # Draw each scaffold at its y-offset
    current_y = 0
    for sd_idx, sd in enumerate(all_scaffolds):
        scaffold = sd['scaffold']
        regions = sd['regions']
        scaffold_length = sd['scaffold_length']
        contig_features = sd['contig_features']
        all_region_features = sd['all_region_features']
        arm_features = sd['arm_features']
        kromatid_mapping = sd.get('kromatid_mapping') or {}
        sample = sd['sample']

        show_arm_line = arm_features is not None and len(arm_features) > 0

        # Title band for this scaffold
        iscn = regions[0].get('iscn_notation', scaffold)
        kromatid_name = kromatid_mapping.get(iscn) if kromatid_mapping else None
        scaffold_size_mb = scaffold_length / 1_000_000

        # Draw full contig
        contig_x = margin_left
        contig_y = current_y + top_margin_per_scaffold
        num_features = len(regions)
        if num_features == 0:
            num_features = 1
        contig_height = num_features * zoom_panel_height + (num_features - 1) * panel_gap

        feature_positions = [(r['center'], r['region_type']) for r in regions]
        d.append(draw.Text(f"{scaffold_size_mb:.1f} Mb", font_size=13, x=contig_x + contig_width / 2,
                           y=contig_y - 6,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle'))

        contig_x_end, contig_deferred_labels = draw_full_contig(
            d, contig_features, colors, featuresets, scaffold_length,
            contig_x, contig_y, contig_height, feature_positions, args,
            arm_features=arm_features, regions=regions,
            show_coord_labels=args.show_contig_breakpoints,
            show_breakpoint_markers=args.show_contig_breakpoints,
            show_scale_bar=not args.hide_scale_bars,
            show_track_labels=(sd_idx == len(all_scaffolds) - 1),
        )

        # Panel letter (A, B, C ...) at the row's top-left corner.
        if getattr(args, 'panel_letters', False):
            d.append(draw.Text(chr(ord('A') + sd_idx), font_size=16,
                               x=8, y=contig_y - 4,
                               fill=text_color, font_family='sans-serif',
                               font_weight='bold', text_anchor='start'))

        # Vertical side label: ISCN notation, rotated -90°
        side_label = iscn
        if kromatid_name:
            side_label = f"{iscn}  /  {kromatid_name}"
        label_x = margin_left - 15
        label_y = contig_y + contig_height / 2
        d.append(draw.Text(side_label, font_size=12, x=label_x, y=label_y,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', dominant_baseline='middle',
                           transform=f'rotate(-90 {label_x} {label_y})'))

        # Zoom panels + fans
        zoom_panel_x = contig_x_end + fan_gap
        all_panel_infos = []

        for i, region in enumerate(regions):
            panel_y = contig_y + i * (zoom_panel_height + panel_gap)
            region_features = all_region_features.get(i, {})
            # Labels go under the bottom row of the whole stack, not under the
            # last region of every sample. Each sample here has one region, so
            # the old test was always true and every row carried its own band.
            is_last_panel = (i == len(regions) - 1) and (sd_idx == len(all_scaffolds) - 1)
            panel_info = draw_zoom_panel_inline(
                d, region, region_features, colors, featuresets,
                zoom_levels, zoom_panel_x, panel_y, zoom_panel_height, args,
                is_last_panel=is_last_panel, scaffold_length=scaffold_length,
                show_breakpoint_lines=args.show_contig_breakpoints,
                fan_stroke=False, gap_between_zooms=gap_between_zooms,
            )
            all_panel_infos.append(panel_info)

            first_zoom = panel_info['zoom_positions'][0]
            first_zoom_view_start = first_zoom['view_start']
            first_zoom_view_end = first_zoom['view_end']

            contig_y_for_zoom_start = contig_y + (first_zoom_view_start / scaffold_length) * contig_height
            contig_y_for_zoom_end = contig_y + (first_zoom_view_end / scaffold_length) * contig_height
            contig_y_for_zoom_start = max(contig_y, min(contig_y + contig_height, contig_y_for_zoom_start))
            contig_y_for_zoom_end = max(contig_y, min(contig_y + contig_height, contig_y_for_zoom_end))

            rect_width = 12
            fan_x1 = contig_x_end
            fan_x2 = panel_info['x_start'] - 1

            fan_points = [
                (fan_x1, contig_y_for_zoom_start),
                (fan_x1 + rect_width, contig_y_for_zoom_start),
                (fan_x2 - rect_width, panel_y),
                (fan_x2, panel_y),
                (fan_x2, panel_y + zoom_panel_height),
                (fan_x2 - rect_width, panel_y + zoom_panel_height),
                (fan_x1 + rect_width, contig_y_for_zoom_end),
                (fan_x1, contig_y_for_zoom_end),
            ]
            path_d = f"M {fan_points[0][0]},{fan_points[0][1]} "
            for px, py in fan_points[1:]:
                path_d += f"L {px},{py} "
            path_d += "Z"
            d.append(draw.Path(d=path_d, fill='#606060', fill_opacity=FAN_FILL_OPACITY,
                               stroke='none'))

        for panel_info in all_panel_infos:
            draw_scale_bars_and_ci(d, panel_info['zoom_positions'],
                                   panel_info['text_color'],
                                   show_scale_bars=not args.hide_scale_bars)

        for label_text, lx, ly, anchor, lcolor in contig_deferred_labels:
            d.append(draw.Text(label_text, font_size=9,
                               x=lx, y=ly,
                               fill=lcolor, font_family='sans-serif',
                               text_anchor=anchor))

        # Accumulate displayed features
        for fs, features in contig_features.items():
            for f in features:
                if isinstance(f, dict):
                    merged_displayed[fs].add(f['feature'])
                else:
                    merged_displayed[fs].add(f[2])
        for region_idx, zoom_data in all_region_features.items():
            for padding, features_by_track in zoom_data.items():
                for fs, features in features_by_track.items():
                    for f in features:
                        if isinstance(f, dict):
                            merged_displayed[fs].add(f['feature'])
                        else:
                            merged_displayed[fs].add(f[2])

        # Advance to the next scaffold's y position
        current_y += scaffold_heights[sd_idx] + inter_scaffold_gap

    # Draw ONE shared legend at the right (spanning all scaffolds)
    legend_x = total_width - margin_right + 10
    legend_y = top_margin_per_scaffold
    show_arm_line = any_arm
    arm_features_for_legend = None
    if show_arm_line:
        # Find any arm_features from any scaffold (for the legend entry)
        for sd in all_scaffolds:
            if sd['arm_features']:
                arm_features_for_legend = sd['arm_features']
                break
    _draw_legend(d, legend_x, legend_y, colors, featuresets,
                 merged_displayed, show_arm_line, arm_features_for_legend,
                 text_color, margin_right=margin_right,
                 draw_border=False)

    return d


def _prepare_scaffold_data(sample, bed_dir, scaffold, scaffold_regions,
                           featuresets, zoom_levels, available_prebinned, args):
    """Load all data needed to draw one scaffold's integrated plot.

    Returns dict with: scaffold_length, contig_features, all_region_features,
    arm_features, kromatid_mapping.
    Returns None if the scaffold cannot be processed.
    """
    # Determine scaffold length
    length_bed_path = os.path.join(
        bed_dir,
        f"{sample}.{args.database}.chromosome.{args.smoothness}.scaffolded.binned1000000.features.bed.gz"
    )
    if not os.path.exists(length_bed_path):
        length_bed_path = length_bed_path.replace('.bed.gz', '.bed')
    scaffold_length = get_scaffold_length(length_bed_path, scaffold)
    if scaffold_length == 0:
        return None

    # Load full contig features (1Mb binning)
    contig_features = {}
    for fs in featuresets:
        bed_path = os.path.join(
            bed_dir,
            f"{sample}.{args.database}.{fs}.{args.smoothness}.scaffolded.binned1000000.features.bed.gz"
        )
        if not os.path.exists(bed_path):
            bed_path = bed_path.replace('.bed.gz', '.bed')
        contig_features[fs] = load_bed_features(bed_path, scaffold, 0, scaffold_length)

    # Arm features (optional)
    arm_features = None
    if args.show_arms:
        region_bed_path = os.path.join(
            bed_dir,
            f"{sample}.{args.database}.region.{args.smoothness}.scaffolded.binned1000000.features.bed.gz"
        )
        if not os.path.exists(region_bed_path):
            region_bed_path = region_bed_path.replace('.bed.gz', '.bed')
        if os.path.exists(region_bed_path):
            all_region_feats = load_bed_features(region_bed_path, scaffold, 0, scaffold_length)
            arm_features = [f for f in all_region_feats
                           if 'p_arm' in f['feature'] or 'q_arm' in f['feature']]

    # Refine breakpoint centers
    for region in scaffold_regions:
        if args.no_refine:
            break
        if region['region_type'] == 'breakpoint':
            original = region['center']
            refined = refine_breakpoint_center(
                region, scaffold, bed_dir, sample,
                args.database, args.smoothness, available_prebinned
            )
            if refined != original:
                refined_start, refined_end = refine_breakpoint_interval(
                    region, refined, scaffold, bed_dir, sample,
                    args.database, args.smoothness, available_prebinned
                )
                print(f"    Refined breakpoint: {original:,} → {refined:,} "
                      f"(shifted {abs(refined - original):,} bp)")
                region['center'] = refined
                region['original_center'] = original
                if refined_start is not None:
                    region['start'] = refined_start
                if refined_end is not None:
                    region['end'] = refined_end
    scaffold_regions.sort(key=lambda r: r['center'])

    # Zoom features per region
    all_region_features = {}
    for region_idx, region in enumerate(scaffold_regions):
        all_features_by_zoom = {}
        for padding in zoom_levels:
            view_start, view_end = compute_view_window(
                region['center'], padding, scaffold_length)
            bin_size = get_bin_size_for_zoom(padding)
            prebinned_size, needs_rebinning = get_prebinned_resolution(bin_size, available_prebinned)
            features_by_track = {}
            for fs in featuresets:
                bed_path = os.path.join(
                    bed_dir,
                    f"{sample}.{args.database}.{fs}.{args.smoothness}.scaffolded.binned{prebinned_size}.features.bed.gz"
                )
                if not os.path.exists(bed_path):
                    bed_path = bed_path.replace('.bed.gz', '.bed')
                raw_features = load_bed_features(bed_path, scaffold, view_start, view_end)
                if needs_rebinning:
                    features_by_track[fs] = bin_features(raw_features, view_start, view_end, bin_size)
                else:
                    features_by_track[fs] = raw_features
            all_features_by_zoom[padding] = features_by_track
        all_region_features[region_idx] = all_features_by_zoom

    return {
        'scaffold_length': scaffold_length,
        'contig_features': contig_features,
        'all_region_features': all_region_features,
        'arm_features': arm_features,
    }


def _process_sample(sample, bed_dir, regions_path, comparison_path,
                    featuresets, zoom_levels, args):
    """Load regions + per-scaffold data for one sample. Returns list of
    dicts (one per scaffold) with all the data needed to draw each scaffold."""
    available_prebinned = discover_prebinned_resolutions(
        bed_dir, sample, args.database, featuresets, args.smoothness
    )

    regions = parse_regions_file(regions_path)
    print(f"  [{sample}] Found {len(regions)} regions in {regions_path}")

    kromatid_mapping = {}
    if comparison_path:
        kromatid_mapping = parse_comparison_file(comparison_path)

    scaffolds = defaultdict(list)
    for region in regions:
        scaffolds[region['scaffold']].append(region)

    out = []
    for scaffold, scaffold_regions in scaffolds.items():
        hap = get_haplotype_from_scaffold(scaffold)
        if hap is None:
            print(f"  Warning: Could not determine haplotype for {scaffold}", file=sys.stderr)
            continue
        scaffold_regions.sort(key=lambda r: r['center'])
        print(f"  [{sample}] Processing {scaffold} ({len(scaffold_regions)} features)...")
        data = _prepare_scaffold_data(
            sample, bed_dir, scaffold, scaffold_regions,
            featuresets, zoom_levels, available_prebinned, args
        )
        if data is None:
            print(f"    Warning: Could not determine length for {scaffold}", file=sys.stderr)
            continue
        out.append({
            'sample': sample,
            'bed_dir': bed_dir,
            'scaffold': scaffold,
            'regions': scaffold_regions,
            'kromatid_mapping': kromatid_mapping,
            **data,
        })
    return out


def main():
    args = parse_args()

    # Parse inputs
    featuresets = [fs.strip() for fs in args.featuresets.split(',')]
    zoom_levels = [parse_size(z) for z in args.zoom_levels.split(',')]

    # Validate parallel input lists
    samples = args.sample
    bed_dirs = args.bed_dir
    regions_files = args.regions
    comparisons = args.comparison if args.comparison else [None] * len(samples)
    if len(bed_dirs) != len(samples) or len(regions_files) != len(samples):
        print("Error: --sample, --bed-dir and --regions must be given the same "
              "number of times.", file=sys.stderr)
        sys.exit(1)
    if args.comparison and len(comparisons) != len(samples):
        print("Error: --comparison (if used) must match the number of samples.",
              file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load colors (shared across all samples — same database)
    print("Loading color files...")
    colors = load_colors(args.colors_dir, args.database, featuresets, background=args.background)

    # Process each sample
    print(f"\nLoading per-sample scaffold data...")
    all_scaffolds = []  # list of scaffold-data dicts
    for sample, bed_dir, regions_path, comparison_path in zip(
            samples, bed_dirs, regions_files, comparisons):
        all_scaffolds.extend(_process_sample(
            sample, bed_dir, regions_path, comparison_path,
            featuresets, zoom_levels, args
        ))

    print(f"\nTotal scaffolds to render: {len(all_scaffolds)}")

    # -----------------------------------------------------------------------
    # Stacked mode: render all scaffolds into one combined SVG
    # -----------------------------------------------------------------------
    if args.stacked_output and len(all_scaffolds) > 0:
        print(f"\nGenerating stacked integrated plot: {args.stacked_output}")
        stacked = draw_stacked_plot(
            all_scaffolds, colors, featuresets, zoom_levels, args
        )
        os.makedirs(os.path.dirname(os.path.abspath(args.stacked_output)) or '.',
                    exist_ok=True)
        stacked.save_svg(args.stacked_output)
        print(f"Saved: {args.stacked_output}")
        # PNG conversion
        png_file = args.stacked_output[:-4] + ".png" if args.stacked_output.endswith('.svg') \
                   else args.stacked_output + ".png"
        from karyoplot.svg.export import svg_to_png
        if svg_to_png(args.stacked_output, png_file, dpi=300, quiet=True):
            print(f"Saved: {png_file}")
        return

    # -----------------------------------------------------------------------
    # Default mode: one SVG per scaffold (backward compatible)
    # -----------------------------------------------------------------------
    print(f"\nGenerating integrated plots...")
    integrated_count = 0
    for sd in all_scaffolds:
        scaffold = sd['scaffold']
        safe_scaffold = scaffold.replace('/', '_')[:50]

        integrated_drawing = draw_integrated_plot(
            scaffold, sd['regions'], sd['scaffold_length'], sd['contig_features'],
            sd['all_region_features'], colors, featuresets, zoom_levels, args,
            kromatid_mapping=sd['kromatid_mapping'],
            arm_features=sd['arm_features'],
        )

        integrated_file = os.path.join(args.output_dir, f"{safe_scaffold}.integrated.svg")
        integrated_drawing.save_svg(integrated_file)
        integrated_count += 1

        # Try post-hoc PNG conversion via rsvg-convert (300 dpi)
        png_file = integrated_file[:-4] + ".png"
        from karyoplot.svg.export import svg_to_png
        svg_to_png(integrated_file, png_file, dpi=300, quiet=True)

    print(f"\nGenerated {integrated_count} integrated plots in {args.output_dir}")


if __name__ == "__main__":
    main()
