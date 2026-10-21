#!/usr/bin/env python3
"""
Stacked breakpoint density plot with KaryoScope zoom panel.

Produces a single SVG/PNG combining all three ROB samples. For each sample:
  - Left: KaryoScope feature tracks (chromosome, acrocentric, region, repeat, gene)
    at a single zoom level, rendered using drawsvg to match the integrated zoom
    plot style from KaryoScope_assembly_contig_zoom_plot.py.
  - Right: Density analysis panels (rotated 90° so genomic position runs
    vertically, aligned with the zoom panel):
      1. Raw chr-specific coverage
      2. Cumulative fraction curves
      3. Combined score with CI bands
      4. SST1 density overlay

Usage:
    python scripts/breakpoint_density_stacked.py [--resolution 100] [--output-dir results/breakpoint_density]
"""

import argparse
import gzip
import os
import re
import subprocess
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import drawsvg as draw
import numpy as np
import pandas as pd
import yaml

# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Color palette (loaded from KaryoScope color files)
# ---------------------------------------------------------------------------
def load_colors(colors_dir, database, featuresets, background='white'):
    """Load color mappings for featuresets (matches zoom plot code)."""
    featureset_colors = {}
    novel_color = "#ffffff" if background == 'white' else "#000000"

    for fs in featuresets:
        colors_path = os.path.join(colors_dir, f"{database}.{fs}.colors.txt")
        featureset_colors[fs] = {
            "novel": (novel_color, 1.0),
            "unknown": ("#808080", 1.0),
        }
        if not os.path.exists(colors_path):
            continue
        with open(colors_path) as f:
            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) >= 2:
                    feature, color = parts[0], parts[1]
                    if i == 0 and feature.lower() == 'feature':
                        continue
                    featureset_colors[fs][feature] = (color, 1.0)
                    if feature.endswith('_specific'):
                        featureset_colors[fs][feature[:-9]] = (color, 1.0)
    return featureset_colors


# ---------------------------------------------------------------------------
# BED loading (simplified from zoom plot code)
# ---------------------------------------------------------------------------
BED_CACHE = {}

def load_bed_file_to_cache(bed_path):
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
    for scaf in scaffold_features:
        scaffold_features[scaf].sort(key=lambda x: x[0])
    BED_CACHE[bed_path] = dict(scaffold_features)


def load_bed_features(bed_path, scaffold, view_start, view_end):
    load_bed_file_to_cache(bed_path)
    features = []
    cached = BED_CACHE.get(bed_path, {})
    for start, end, feature in cached.get(scaffold, []):
        if end > view_start and start < view_end:
            features.append({
                'start': max(start, view_start),
                'end': min(end, view_end),
                'feature': feature,
            })
    return features


def bin_features(features, view_start, view_end, bin_size):
    """Bin features using max-overlap logic (matches zoom plot code)."""
    if not features:
        return []
    bins = defaultdict(lambda: defaultdict(int))
    for feat in features:
        start_bin = feat['start'] // bin_size
        end_bin = (feat['end'] - 1) // bin_size
        for b_idx in range(start_bin, end_bin + 1):
            bs = b_idx * bin_size
            be = bs + bin_size
            overlap = min(feat['end'], be) - max(feat['start'], bs)
            if overlap > 0:
                bins[b_idx][feat['feature']] += overlap
    binned = []
    prev_feat = None
    prev_start = None
    prev_end = None
    for b_idx in sorted(bins.keys()):
        bs = b_idx * bin_size
        be = bs + bin_size
        if be <= view_start or bs >= view_end:
            continue
        cs = max(bs, view_start)
        ce = min(be, view_end)
        best = max(bins[b_idx], key=lambda f: (bins[b_idx][f], f != 'novel'))
        if best == 'novel' and bins[b_idx][best] <= bin_size / 2:
            non_novel = {f: v for f, v in bins[b_idx].items() if f != 'novel'}
            if non_novel:
                best = max(non_novel, key=non_novel.get)
        if best == prev_feat and cs == prev_end:
            binned[-1]['end'] = ce
            prev_end = ce
        else:
            binned.append({'start': cs, 'end': ce, 'feature': best})
            prev_feat = best
            prev_start = cs
            prev_end = ce
    return binned


def get_bin_size_for_zoom(padding):
    if padding >= 10_000_000:
        return 100_000
    elif padding >= 5_000_000:
        return 50_000
    elif padding >= 2_000_000:
        return 10_000
    else:
        return 1_000


def get_scale_bar_config(padding):
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


def _format_genomic_pos(bp_pos):
    if bp_pos >= 1_000_000:
        val = bp_pos / 1_000_000
        if val == int(val):
            return f"{int(val)} Mb"
        return f"{val:.1f} Mb"
    elif bp_pos >= 1_000:
        val = bp_pos / 1_000
        if val == int(val):
            return f"{int(val)} kb"
        return f"{val:.1f} kb"
    return f"{bp_pos} bp"


def format_width(val):
    if val >= 1_000_000:
        return f"{val / 1_000_000:.2f} Mb"
    elif val >= 1_000:
        return f"{val / 1_000:.1f} kb"
    return f"{val} bp"


# ---------------------------------------------------------------------------
# ROB definitions
# ---------------------------------------------------------------------------
def load_rob_definitions():
    config_path = Path(os.environ.get("ROBS_CONFIG", PROJECT_ROOT / "conf" / "robs.yaml"))
    with open(config_path) as f:
        config = yaml.safe_load(f)

    robs = []

    bp = config["breakpoint_groups"]["gm03417_robertsonian"][0]
    robs.append({
        "sample": "GM03417",
        "iscn": "rob(14;21)(q10;q10)",
        "scaffold": bp["scaffold"],
        "chrom_left": "chr21",
        "chrom_right": "chr14",
        "label_left": "chr21_specific",
        "label_right": "chr14_specific",
        "iscn_center": bp["center"],
        "iscn_start": bp["start"],
        "iscn_end": bp["end"],
        "sst1_center": bp["center"],
    })

    bp = config["breakpoint_groups"]["gm03786_robertsonian"][0]
    robs.append({
        "sample": "GM03786",
        "iscn": "rob(13;14)(q10;q10)",
        "scaffold": bp["scaffold"],
        "chrom_left": "chr14",
        "chrom_right": "chr13",
        "label_left": "chr14_specific",
        "label_right": "chr13_specific",
        "iscn_center": bp["center"],
        "iscn_start": bp["start"],
        "iscn_end": bp["end"],
        "sst1_center": bp["center"],
    })

    bp = config["breakpoint_groups"]["gm04890_robertsonian"][0]
    robs.append({
        "sample": "GM04890",
        "iscn": "rob(13;14)(q10;q10)",
        "scaffold": bp["scaffold"],
        "chrom_left": "chr14",
        "chrom_right": "chr13",
        "label_left": "chr14_specific",
        "label_right": "chr13_specific",
        "iscn_center": bp["center"],
        "iscn_start": bp["start"],
        "iscn_end": bp["end"],
        "sst1_center": bp["center"],
    })

    return robs, config


# ---------------------------------------------------------------------------
# Cumulative breakpoint analysis
# ---------------------------------------------------------------------------
def cumulative_breakpoint_analysis(bed_df, scaffold, label_left, label_right,
                                   focus_start, focus_end, resolution=100):
    scaffold_df = bed_df[bed_df["chrom"] == scaffold]
    left_df = scaffold_df[scaffold_df["label"] == label_left]
    right_df = scaffold_df[scaffold_df["label"] == label_right]

    n_bins = int((focus_end - focus_start) / resolution)
    positions = np.arange(focus_start, focus_start + n_bins * resolution, resolution)
    left_cov = np.zeros(n_bins)
    right_cov = np.zeros(n_bins)

    for _, row in left_df.iterrows():
        if row["end"] <= focus_start or row["start"] >= focus_end:
            continue
        s = max(0, int((row["start"] - focus_start) / resolution))
        e = min(n_bins, int(np.ceil((row["end"] - focus_start) / resolution)))
        for i in range(s, e):
            bin_start = focus_start + i * resolution
            bin_end = bin_start + resolution
            overlap = min(row["end"], bin_end) - max(row["start"], bin_start)
            left_cov[i] += max(0, overlap)

    for _, row in right_df.iterrows():
        if row["end"] <= focus_start or row["start"] >= focus_end:
            continue
        s = max(0, int((row["start"] - focus_start) / resolution))
        e = min(n_bins, int(np.ceil((row["end"] - focus_start) / resolution)))
        for i in range(s, e):
            bin_start = focus_start + i * resolution
            bin_end = bin_start + resolution
            overlap = min(row["end"], bin_end) - max(row["start"], bin_start)
            right_cov[i] += max(0, overlap)

    left_cumsum = np.cumsum(left_cov)
    right_cumsum = np.cumsum(right_cov[::-1])[::-1]
    total_left = left_cumsum[-1] if left_cumsum[-1] > 0 else 1
    total_right = right_cumsum[0] if right_cumsum[0] > 0 else 1
    left_frac = left_cumsum / total_left
    right_frac = right_cumsum / total_right
    combined = left_frac + right_frac

    max_score = combined.max()
    at_max = np.isclose(combined, max_score, rtol=1e-10)
    idx_max = np.where(at_max)[0]
    plateau_start = int(positions[idx_max[0]])
    plateau_end = int(positions[idx_max[-1]])
    plateau_center = (plateau_start + plateau_end) // 2

    def ci(threshold_frac):
        """Interval where the combined score stays within (1 - frac) of its max.

        NOTE ON THE NAME: `ci(0.99)` is what the manuscript calls the "99%
        plateau width" -- Methods define it as "the genomic interval over which
        the combined score remains within 1% of its maximum", which is exactly
        this. It is NOT a confidence interval, and it is NOT the argmax plateau
        (`plateau_start`/`plateau_end`), which is three orders of magnitude
        narrower. The `ci_` names are kept only because the emitted TSV columns
        already carry them.
        """
        mask = combined >= max_score * threshold_frac
        idx = np.where(mask)[0]
        return int(positions[idx[0]]), int(positions[idx[-1]])

    return {
        "positions": positions,
        "left_cov": left_cov,
        "right_cov": right_cov,
        "left_frac": left_frac,
        "right_frac": right_frac,
        "combined": combined,
        "max_score": max_score,
        "plateau_start": plateau_start,
        "plateau_end": plateau_end,
        "plateau_center": plateau_center,
        "ci_999": ci(0.999),
        "ci_99": ci(0.99),
        "ci_95": ci(0.95),
        "total_left_bp": int(total_left),
        "total_right_bp": int(total_right),
    }


def get_sst1_density(bed_df, scaffold, focus_start, focus_end, resolution=100):
    scaffold_df = bed_df[(bed_df["chrom"] == scaffold) & (bed_df["label"] == "SST1_specific")]
    n_bins = int((focus_end - focus_start) / resolution)
    positions = np.arange(focus_start, focus_start + n_bins * resolution, resolution)
    sst1_cov = np.zeros(n_bins)
    for _, row in scaffold_df.iterrows():
        if row["end"] <= focus_start or row["start"] >= focus_end:
            continue
        s = max(0, int((row["start"] - focus_start) / resolution))
        e = min(n_bins, int(np.ceil((row["end"] - focus_start) / resolution)))
        for i in range(s, e):
            bin_start = focus_start + i * resolution
            bin_end = bin_start + resolution
            overlap = min(row["end"], bin_end) - max(row["start"], bin_start)
            sst1_cov[i] += max(0, overlap)
    return positions, sst1_cov


def _read_scaffold_window(path, scaffold, focus_start, focus_end):
    """Stream a per-base BED, keeping only `scaffold` rows overlapping the window.

    These files are the per-base `presmoothed.scaffolded` BEDs -- 60-70 M rows,
    330-380 MB gzipped. Reading one whole into pandas materialises well over
    10 GB and gets the process OOM-killed, which is why the loaders below used
    to read the 100 bp binned copies instead. The binning was winner-take-all,
    so a 100 bp bin was assigned to whichever label covered most of it and the
    rest of the bin was discarded; that quantisation is what widened the
    breakpoint interval and made it unstable across k.

    Filtering to one scaffold and one focus window before pandas sees anything
    is exact rather than an approximation: KaryoScope annotations tile the
    input and are written in genome order, so every base in the window is
    represented in exactly one retained record. What reaches the caller is a
    few 10^5 rows.

    The prefilter runs in awk because a Python-level loop over 60 M lines costs
    minutes per file and there are six such reads per figure.
    """
    if not Path(path).exists():
        raise FileNotFoundError(
            f"per-base BED not found: {path}\n"
            "Fig. S5 reads the per-base presmoothed.scaffolded BEDs, not the "
            "binned100 copies. Stage them first (see run_stage_inputs.sh)."
        )
    prog = '$1==s && $3>a && $2<b'
    p1 = subprocess.Popen(["zcat", str(path)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(
        ["awk", "-F", "\t", "-v", f"s={scaffold}", "-v", f"a={focus_start}",
         "-v", f"b={focus_end}", prog],
        stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    df = pd.read_csv(p2.stdout, sep="\t", header=None,
                     names=["chrom", "start", "end", "label"],
                     dtype={"chrom": str, "start": int, "end": int, "label": str})
    p2.wait()
    p1.wait()
    if p2.returncode not in (0, None):
        raise RuntimeError(f"prefilter failed on {path} (awk rc={p2.returncode})")
    if df.empty:
        raise RuntimeError(
            f"no records for scaffold {scaffold!r} in {focus_start}-{focus_end} "
            f"from {path}. An empty window here yields a blank panel rather than "
            "an error, so this is raised deliberately -- check the scaffold name "
            "matches the scaffolded frame."
        )
    return df


def _per_base_path(sample, featureset, config):
    db = config["database"]
    assembly = config["samples"][sample]["assembly_name"]
    bed_dir = config["samples"][sample]["bed_dir"]
    fname = f"{assembly}.{db}.{featureset}.presmoothed.scaffolded.bed.gz"
    return PROJECT_ROOT / bed_dir / fname


def load_chrom_bed(sample, resolution, config, scaffold, focus_start, focus_end):
    return _read_scaffold_window(_per_base_path(sample, "chromosome", config),
                                 scaffold, focus_start, focus_end)


def load_acro_bed(sample, resolution, config, scaffold, focus_start, focus_end):
    return _read_scaffold_window(_per_base_path(sample, "acrocentric", config),
                                 scaffold, focus_start, focus_end)


# ---------------------------------------------------------------------------
# Drawing helpers (matching zoom plot style exactly)
# ---------------------------------------------------------------------------
def draw_feature_track(d, features, colors, x, y_start, bar_width,
                       view_start, pixels_per_bp, track_height):
    """Draw vertical track of colored rectangles (matches zoom plot code)."""
    for feat in features:
        y = y_start + (feat['start'] - view_start) * pixels_per_bp
        height = (feat['end'] - feat['start']) * pixels_per_bp
        color, opacity = colors.get(feat['feature'], ('#808080', 1.0))
        d.append(draw.Rectangle(x, y, bar_width, height,
                                fill=color, fill_opacity=opacity))
    # Outer track outline (white background)
    d.append(draw.Rectangle(x, y_start, bar_width, track_height,
                            fill='none', stroke='#000000', stroke_width=0.5))


def _clean_label(label):
    """Display form of a feature label.

    Drops the `_specific` suffix and any `multigroupN` tail. Those encode where
    a label sits in the hierarchy -- `_specific` marks a leaf, `multigroupN`
    marks an internal node that several leaves resolve to -- which is a fact
    about the database, not something a reader of this figure needs.
    """
    lab = label.replace('_specific', '')
    lab = re.sub(r'[ _]?multigroup\d+$', '', lab)
    return lab.replace('_', ' ').strip()


def draw_coordinate_labels(d, x, y_start, y_end, view_start, view_end,
                           text_color='black', font_size=11):
    top_label = _format_genomic_pos(view_start)
    bot_label = _format_genomic_pos(view_end)
    d.append(draw.Text(top_label, font_size=font_size, x=x, y=y_start + 3,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='end', dominant_baseline='hanging'))
    d.append(draw.Text(bot_label, font_size=font_size, x=x, y=y_end - 3,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='end', dominant_baseline='auto'))


def draw_scale_bar(d, x, y, height, label, text_color='black'):
    """Draw a scale bar with notches and rotated label."""
    d.append(draw.Rectangle(x, y, 1.5, height, fill=text_color))
    d.append(draw.Line(x - 1, y, x + 2.5, y, stroke=text_color, stroke_width=1))
    d.append(draw.Line(x - 1, y + height, x + 2.5, y + height,
                       stroke=text_color, stroke_width=1))
    label_x = x - 4
    label_y = y + height / 2
    d.append(draw.Text(label, font_size=11, x=label_x, y=label_y,
                       fill=text_color, font_family='sans-serif',
                       text_anchor='middle', dominant_baseline='middle',
                       transform=f"rotate(-90 {label_x} {label_y})"))


# ---------------------------------------------------------------------------
# Density panel drawing (rotated 90° — genomic position vertical)
# ---------------------------------------------------------------------------
def draw_density_panel_coverage(d, x_start, y_start, panel_width, panel_height,
                                positions, left_cov, right_cov,
                                view_start, view_end,
                                left_color, right_color,
                                left_label, right_label,
                                cov_frac=0.5):
    """Chromosome-specific k-mer panel: two half-width columns, one per partner.

    Bars are drawn at CONSTANT full sub-panel width rather than scaled by
    coverage. In the published figure every mark spans its whole half of the
    rectangle and the information is carried by presence/absence, not by bar
    length. White gaps are bins at or below the threshold.

    `cov_frac` is why this does not come out denser than the published figure.
    That figure was drawn from 100 bp WINNER-TAKE-ALL bins, where a bin showed a
    chromosome only if that chromosome won the bin -- an implicit majority
    filter. Per-base input has no such filter, so lighting every bin with any
    coverage lights 3-7x more of them (measured on the chromosome tracks of the
    three samples). Requiring more than
    `cov_frac` of the bin restores the published density. The ANALYSIS still
    uses true fractional coverage -- that is the point of the per-base switch,
    and this display threshold does not touch it.
    """
    n = len(positions)
    resolution = positions[1] - positions[0] if n > 1 else 100
    ppbp = panel_height / (view_end - view_start)

    chr_half = panel_width / 2
    chr_center_x = x_start + chr_half
    thresh = cov_frac * resolution

    for i in range(n):
        y = y_start + (positions[i] - view_start) * ppbp
        h = max(resolution * ppbp, 0.3)
        # Full-width marks: the left partner fills the left half, the right
        # partner the right half. No proportional scaling.
        if left_cov[i] > thresh:
            d.append(draw.Rectangle(chr_center_x - chr_half, y, chr_half, h,
                                    fill=left_color, fill_opacity=0.85))
        if right_cov[i] > thresh:
            d.append(draw.Rectangle(chr_center_x, y, chr_half, h,
                                    fill=right_color, fill_opacity=0.85))

    d.append(draw.Line(chr_center_x, y_start, chr_center_x, y_start + panel_height,
                       stroke='#000000', stroke_width=0.5, stroke_opacity=0.3))
    d.append(draw.Rectangle(x_start, y_start, panel_width, panel_height,
                            fill='none', stroke='#000000', stroke_width=0.5))

    # Labels sit INSIDE the rectangle, centred over their own half. The left
    # partner is labelled at the bottom and the right partner at the top, which
    # is how the published figure places them.
    lx = x_start + chr_half / 2
    rx = chr_center_x + chr_half / 2
    ll, rl = _clean_label(left_label), _clean_label(right_label)

    def _labelled(text, cx, cy, colour, baseline):
        # Opaque backing: these labels sit ON the track, and without it a dense
        # run of marks behind them makes them unreadable.
        w = len(text) * 6.1 + 5      # backing width tracks the font size
        top = cy - 8 if baseline == 'auto' else cy - 1
        d.append(draw.Rectangle(cx - w / 2, top, w, 9,
                                fill='#FFFFFF', fill_opacity=0.85))
        d.append(draw.Text(text, font_size=11, x=cx, y=cy, fill=colour,
                           font_family='sans-serif', text_anchor='middle',
                           dominant_baseline=baseline))

    _labelled(ll, lx, y_start + panel_height - 4, left_color, 'auto')
    _labelled(rl, rx, y_start + 4, right_color, 'hanging')


def draw_sst1_panel(d, x_start, y_start, panel_width, panel_height,
                    positions, sst1_cov, view_start, view_end,
                    cov_frac=0.5):
    """Standalone SST1 coverage strip, drawn full-width like the k-mer panel.

    Same majority threshold, for a slightly different reason. SST1 bins are
    mostly near-full already (most are over half covered), so this barely
    changes which bins survive. What it removes is the handful of 1-5 bp bins,
    which at full width and full opacity carried the same visual weight as a
    whole array and made this strip read as busier than the published one.
    """
    n = len(positions)
    resolution = positions[1] - positions[0] if n > 1 else 100
    ppbp = panel_height / (view_end - view_start)
    thresh = cov_frac * resolution
    for i in range(n):
        if sst1_cov[i] > thresh:
            y = y_start + (positions[i] - view_start) * ppbp
            h = max(resolution * ppbp, 0.3)
            d.append(draw.Rectangle(x_start, y, panel_width, h,
                                    fill='#00FF00', fill_opacity=0.9))
    d.append(draw.Rectangle(x_start, y_start, panel_width, panel_height,
                            fill='none', stroke='#000000', stroke_width=0.5))


def draw_density_panel_combined(d, x_start, y_start, panel_width, panel_height,
                                positions, combined, max_score,
                                view_start, view_end,
                                plateau_start, plateau_end,
                                ci_999, ci_99,
                                sst1_center,
                                left_frac=None, right_frac=None,
                                left_color='#000000', right_color='#000000'):
    """Draw combined score curve with CI annotations."""
    ppbp = panel_height / (view_end - view_start)

    # Score range for x-axis mapping
    score_min = combined.min()
    score_range = max_score - score_min
    if score_range == 0:
        score_range = 1

    # The yellow band: ci_99, which IS the manuscript's "99% plateau width"
    # (score within 1% of maximum). See the note on ci() above.
    ci99_y_start = y_start + (ci_99[0] - view_start) * ppbp
    ci99_y_end = y_start + (ci_99[1] - view_start) * ppbp
    ci99_y_start = max(y_start, min(y_start + panel_height, ci99_y_start))
    ci99_y_end = max(y_start, min(y_start + panel_height, ci99_y_end))
    d.append(draw.Rectangle(x_start, ci99_y_start, panel_width,
                            ci99_y_end - ci99_y_start,
                            fill='#FBBF24', fill_opacity=0.25))

    # No lines are drawn at conf/robs.yaml's start/end. Those are hand-set
    # display bounds (see the note in that file): the focus window is built by
    # padding them 1 Mb, so lines at them land exactly 1 Mb inside each panel
    # edge by construction and say nothing about the data. Removed 2026-08-28
    # on the PI's decision, together with the width bracket and the
    # "Transition interval" legend entry that named them.

    # SST1 center line
    sst1_y = y_start + (sst1_center - view_start) * ppbp
    if y_start <= sst1_y <= y_start + panel_height:
        d.append(draw.Line(x_start, sst1_y, x_start + panel_width, sst1_y,
                           stroke='#C4A9E8', stroke_width=2, stroke_dasharray='4,2'))

    # Two cumulative curves, one per fusion partner, coloured to match the
    # chromosome they belong to. `combined` (their sum) is what the breakpoint
    # is called from, but plotting the sum as a single black line hides the
    # crossing that IS the breakpoint -- the published figure draws both.
    def _draw_curve(values, colour):
        max_path_points = 500
        vals = values
        if len(vals) > max_path_points:
            from scipy.ndimage import gaussian_filter1d
            vals = gaussian_filter1d(vals, sigma=len(vals) / max_path_points)
            step = max(1, len(positions) // max_path_points)
            idx = np.arange(0, len(positions), step)
            if idx[-1] != len(positions) - 1:
                idx = np.append(idx, len(positions) - 1)
            pp, pv = positions[idx], vals[idx]
        else:
            pp, pv = positions, vals
        pts = []
        for i, pos in enumerate(pp):
            yy = y_start + (pos - view_start) * ppbp
            # Both curves are fractions on [0, 1], so they share one axis and
            # need no per-curve rescaling.
            xx = x_start + float(np.clip(pv[i], 0.0, 1.0)) * panel_width
            pts.append((xx, yy))
        if len(pts) >= 2:
            path = draw.Path(fill='none', stroke=colour, stroke_width=1.5)
            path.M(pts[0][0], pts[0][1])
            for px, py in pts[1:]:
                path.L(px, py)
            d.append(path)

    if left_frac is not None and right_frac is not None:
        _draw_curve(left_frac, left_color)
        _draw_curve(right_frac, right_color)
    else:
        _draw_curve((combined - combined.min()) / max(max_score - combined.min(), 1e-9),
                    '#000000')

    # Panel outline
    d.append(draw.Rectangle(x_start, y_start, panel_width, panel_height,
                            fill='none', stroke='#000000', stroke_width=0.5))


def draw_density_panel_sst1(d, x_start, y_start, panel_width, panel_height,
                            sst1_positions, sst1_cov, view_start, view_end):
    """Draw SST1 density as horizontal bars (standalone panel, kept for compatibility)."""
    ppbp = panel_height / (view_end - view_start)
    max_cov = sst1_cov.max() if sst1_cov.max() > 0 else 1
    resolution = sst1_positions[1] - sst1_positions[0] if len(sst1_positions) > 1 else 100

    for i, pos in enumerate(sst1_positions):
        if sst1_cov[i] <= 0:
            continue
        y = y_start + (pos - view_start) * ppbp
        h = resolution * ppbp
        if h < 0.3:
            h = 0.3
        bar_w = (sst1_cov[i] / max_cov) * panel_width
        d.append(draw.Rectangle(x_start, y, bar_w, h,
                                fill='#C4A9E8', fill_opacity=0.8))

    # Panel outline
    d.append(draw.Rectangle(x_start, y_start, panel_width, panel_height,
                            fill='none', stroke='#000000', stroke_width=0.5))


# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
# Depth of each node within its feature set, used to pick the legend
# representative when several features share one colour. Built lazily from the
# database's hierarchy.tsv; missing entries fall back to a large depth.
_HIER_DEPTH = {}


def _load_hier_depth(hier_path):
    parent = {}
    try:
        with open(hier_path) as f:
            f.readline()
            for line in f:
                fs, child, par = line.rstrip("\n").split("\t")[:3]
                parent[(fs, child)] = (fs, par)
    except OSError:
        return
    for key in parent:
        depth, cur, seen = 0, key, set()
        while cur in parent and cur not in seen:
            seen.add(cur)
            cur = parent[cur]
            depth += 1
        _HIER_DEPTH[key] = depth


# ---------------------------------------------------------------------------
# Text metrics
# ---------------------------------------------------------------------------
# The legend box used to be sized against a fixed column width, i.e. against
# ASSUMED glyph widths. Both this SVG and the assembled figure declare only the
# generic `sans-serif` family and embed no font, so each renderer is free to
# resolve it differently: rsvg-convert picks DejaVu Sans here, which is wider
# than the face the 118 px column was sized against, and the two longest
# Density Analysis labels crossed the border in the PNG while the SVG looked
# fine. Measuring the strings and sizing the box from the result survives a
# renderer swap; widening the box by a fixed margin does not, quite.

#: Deliberately WIDE per-character advances (em), used only when the real font
#: cannot be read. Sizing too generously leaves white space; sizing too tightly
#: puts glyphs through the border, so the fallback errs generous.
_FALLBACK_EM = {"upper": 0.73, "lower": 0.64, "digit": 0.70, "space": 0.35, "other": 0.42}


@lru_cache(maxsize=1)
def _sans_advances():
    """Advance widths (in em) keyed by character for the resolved `sans-serif`.

    Returns None if fontconfig or fontTools is unavailable, in which case
    :func:`text_width` falls back to `_FALLBACK_EM`.
    """
    try:
        from fontTools.ttLib import TTFont

        path = subprocess.run(["fc-match", "-f", "%{file}", "sans-serif"],
                              capture_output=True, text=True, check=True).stdout.strip()
        if not path:
            return None
        font = TTFont(path, fontNumber=0)
        upm = font["head"].unitsPerEm
        cmap = font.getBestCmap()
        hmtx = font["hmtx"]
        return {chr(cp): hmtx[g][0] / upm for cp, g in cmap.items() if g in hmtx.metrics}
    except Exception:
        return None


def text_width(text, font_size):
    """Rendered width of `text` at `font_size`, in user units."""
    adv = _sans_advances()
    total = 0.0
    for ch in text:
        if adv is not None:
            total += adv.get(ch, _FALLBACK_EM["other"])
        elif ch.isupper():
            total += _FALLBACK_EM["upper"]
        elif ch.islower():
            total += _FALLBACK_EM["lower"]
        elif ch.isdigit():
            total += _FALLBACK_EM["digit"]
        elif ch == " ":
            total += _FALLBACK_EM["space"]
        else:
            total += _FALLBACK_EM["other"]
    return total * font_size


def draw_legend(d, x, y, colors, featuresets, displayed_features,
                text_color='black'):
    """Two-column legend.

    The published figure runs the legend in two columns so it fits beside the
    other panels; a single column is too tall. Column 1 carries Chromosome,
    Acrocentric and Satellite, column 2 Repeat, Gene and Density Analysis.

    Returns ``(height, width)``. Both are measured from the rendered text, so a
    caller placing the legend must use the returned width rather than assume one.
    """
    item_height = 12
    swatch_size = 8
    # 12 -> 8. The legend renders at ~0.92 of natural size to fit panel C's
    # box, so every point of natural size costs ~8% on screen; recovering 12 px
    # of dead gap is what lets the item font go to 11 (10.1 px effective)
    # instead of 10 (9.2 px) inside the same rectangle.
    section_gap = 8
    header_size = 12
    item_size = 11
    item_indent = swatch_size + 8   # text sits this far right of the column edge
    col_gap = 8                     # between the two columns
    box_pad = 5                     # between the content and the drawn border

    featureset_names = {
        'chromosome': 'Chromosome',
        'acrocentric': 'Acrocentric',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
    }

    def collect(fs):
        """(label, colour) rows for one feature set.

        Sorted alphabetically by the DISPLAYED label. The previous sort used the
        raw feature name, which is why `active` appeared after `divergent`.
        Background-style catch-alls (`nonrepeat`, `nonacrocentric`) are pinned
        last, as in the published legend -- they are the leftover class, not a
        peer of the named ones.

        Representative choice matters here. In this database all six internal
        chromosome nodes -- autosome, acrocentric, metacentric, submetacentric,
        sex and categorized -- share ONE colour (#B0C4DE), so they collapse to a
        single legend row and only one name can be shown. Picking it
        alphabetically gave `acrocentric`; the published figure shows
        `autosome`. Neither is more correct on its own: the swatch means "shared
        by more than one chromosome". We take the SHALLOWEST node in the
        hierarchy, since that is the one whose meaning covers every other member
        of the group, excluding the artificial root `categorized`.
        """
        fs_colors = colors.get(fs, {})
        fs_displayed = displayed_features.get(fs, set())
        color_groups = defaultdict(list)
        for feat in fs_displayed:
            c, _ = fs_colors.get(feat, ('#808080', 1.0))
            color_groups[c].append(feat)
        rows = []
        for c, feats in color_groups.items():
            specifics = [f for f in feats if f.endswith('_specific')]
            if specifics:
                rep = _clean_label(sorted(specifics)[0])
            else:
                cands = [f for f in feats if f != 'categorized'] or list(feats)
                cands.sort(key=lambda f: (_HIER_DEPTH.get((fs, f), 99), f))
                rep = _clean_label(cands[0])
            rows.append((rep, c))
        # Only the gap-fill classes are pinned last. `unknown` is NOT one of
        # them -- the published legend sorts it alphabetically with the named
        # repeat classes (… SINE, Unknown, nonrepeat).
        BACKGROUNDISH = ('nonrepeat', 'nonacrocentric', 'novel')
        rows.sort(key=lambda r: (r[0].lower() in BACKGROUNDISH, r[0].lower()))
        return rows

    density_items = [
        ('#00FF00', 0.9, 'SST1 coverage'),
        ('#FBBF24', 0.25, '99% plateau width'),
        (None, None, 'SST1 center'),
    ]

    def draw_swatch_row(label, colour, opacity, cx, cy):
        if colour:
            d.append(draw.Rectangle(cx + 3, cy, swatch_size, swatch_size,
                                    fill=colour, fill_opacity=opacity if opacity else 1.0,
                                    stroke=text_color, stroke_width=0.5))
        elif label == 'SST1 center':
            my = cy + swatch_size / 2
            d.append(draw.Line(cx + 3, my, cx + 3 + swatch_size, my,
                               stroke='#C4A9E8', stroke_width=2, stroke_dasharray='4,2'))
        d.append(draw.Text(label, font_size=item_size,
                           x=cx + item_indent, y=cy + swatch_size / 2,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='start', dominant_baseline='middle'))

    col1 = [fs for fs in ('chromosome', 'acrocentric', 'region') if fs in featuresets]
    col2 = [fs for fs in ('repeat', 'gene') if fs in featuresets]
    rows_by_fs = {fs: collect(fs) for fs in col1 + col2}

    def column_width(fs_list, with_density):
        """Width the widest rendered string in this column actually needs.

        The two columns are sized independently, so column 1 (short satellite
        and chromosome names) does not hold width that column 2 (`99% plateau
        width`, and formerly `Transition interval`) needs. An equal 118 px split
        is what put those labels through the border.
        """
        w = 0.0
        for fs in fs_list:
            if not rows_by_fs.get(fs):
                continue
            w = max(w, text_width(featureset_names.get(fs, fs), header_size))
            for label, _colour in rows_by_fs[fs]:
                w = max(w, item_indent + text_width(label, item_size))
        if with_density:
            w = max(w, text_width("Density Analysis", header_size))
            for _colour, _opacity, label in density_items:
                w = max(w, item_indent + text_width(label, item_size))
        return w

    col1_w = column_width(col1, with_density=False)
    col2_w = column_width(col2, with_density=True)

    def draw_column(fs_list, cx, with_density):
        cy = y
        for fs in fs_list:
            rows = rows_by_fs.get(fs)
            if not rows:
                continue
            d.append(draw.Text(featureset_names.get(fs, fs), font_size=header_size, x=cx, y=cy,
                               fill=text_color, font_family='sans-serif',
                               font_weight='bold', dominant_baseline='hanging'))
            cy += item_height + 2
            for label, colour in rows:
                draw_swatch_row(label, colour, 1.0, cx, cy)
                cy += item_height
            cy += section_gap
        if with_density:
            d.append(draw.Text("Density Analysis", font_size=header_size, x=cx, y=cy,
                               fill=text_color, font_family='sans-serif',
                               font_weight='bold', dominant_baseline='hanging'))
            cy += item_height + 2
            for colour, opacity, label in density_items:
                draw_swatch_row(label, colour, opacity, cx, cy)
                cy += item_height
        return cy

    y1 = draw_column(col1, x, with_density=False)
    y2 = draw_column(col2, x + col1_w + col_gap, with_density=True)

    legend_height = max(y1, y2) - y + 10
    legend_width = col1_w + col_gap + col2_w + 2 * box_pad
    d.append(draw.Rectangle(x - box_pad, y - box_pad, legend_width, legend_height,
                            fill='none', stroke=text_color, stroke_width=0.5,
                            rx=3, ry=3))
    return legend_height, legend_width


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=100,
                        help="Density analysis bin resolution in bp (default: 100)")
    parser.add_argument("--output-dir", type=str,
                        default="results/breakpoint_density")
    parser.add_argument("--padding", type=float, default=1.0,
                        help="Extra padding beyond ISCN CI in Mb (default: 1.0)")
    parser.add_argument("--no-legend", action="store_true",
                        help="Omit the legend. The assembled Figure 5 positions "
                             "it separately, under panel D.")
    parser.add_argument("--emit-anchors", type=str, default=None,
                        help="Write panel geometry as JSON: the SST1 array "
                             "y-centres the assembler draws leader lines to, and "
                             "each sample's panel box.")
    parser.add_argument("--cov-frac", type=float, default=0.5,
                        help="Display threshold for the k-mer and SST1 tracks: a "
                             "bin is drawn only if the feature covers more than "
                             "this fraction of it (default 0.5, i.e. a majority). "
                             "This reproduces the winner-take-all binning the "
                             "published figure was drawn from. Display only -- it "
                             "does not affect the breakpoint analysis, which uses "
                             "true fractional coverage. Use 0 for every bin with "
                             "any coverage.")
    parser.add_argument("--db-root", type=str,
                        default="/tgen_labs/barthel/software/github/barthel/"
                                "KaryoScope/resources/registry",
                        help="Registry root, used only to find the database's "
                             "hierarchy.tsv for legend grouping.")
    parser.add_argument("--colors-dir", type=str, default=None,
                        help="Override conf/robs.yaml colors_dir. That default is "
                             "a macOS path; load_colors() skips a missing colors "
                             "file silently, so without a valid directory every "
                             "feature falls back to grey 'unknown' and the figure "
                             "renders wrong at exit 0.")
    args = parser.parse_args()

    output_dir = PROJECT_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    robs, config = load_rob_definitions()

    # Load KaryoScope color palette
    colors_dir = os.path.expanduser(args.colors_dir or config["colors_dir"])
    # Fail loudly. load_colors() `continue`s past a missing file, so a bad
    # colours path produces an all-grey figure and exit 0 rather than an error.
    if not os.path.isdir(colors_dir):
        raise SystemExit(
            f"colors_dir does not exist: {colors_dir}\n"
            "Pass --colors-dir pointing at the KS_human_CHM13 database directory "
            "holding KS_human_CHM13.<featureset>.colors.txt")
    database = config["database"]
    featuresets = ['chromosome', 'acrocentric', 'region', 'repeat', 'gene']
    ks_colors = load_colors(colors_dir, database, featuresets, background='white')
    # Node depths for the legend representative rule. Prefer a hierarchy.tsv
    # beside the colours; otherwise use the installed database's.
    _load_hier_depth(os.path.join(colors_dir, "hierarchy.tsv"))
    if not _HIER_DEPTH:
        _load_hier_depth(os.path.join(args.db_root, database, "hierarchy.tsv"))
    if not _HIER_DEPTH:
        print("  note: no hierarchy.tsv found; legend representatives fall back "
              "to alphabetical order")

    # --- Layout constants (matching zoom plot code) ---
    background = 'white'
    text_color = 'black'
    bar_width = 10
    bar_spacing = 2
    num_tracks = len(featuresets)
    track_block_width = num_tracks * bar_width + (num_tracks - 1) * bar_spacing

    margin_left = 62
    margin_right = 265 if not args.no_legend else 12
    margin_top = 15
    # Vertical gap between the sample panels. 20 left a 38 px band of white
    # between A/B and B/C; 8 brings it to 26, which is what the published figure
    # shows and still clears the coordinate labels.
    margin_bottom = 8
    inter_scaffold_gap = 0       # match stacked zoom plot
    # Panel height drives the whole figure's aspect. At 250 the assembled
    # figure came out square; the manuscript figure is portrait (w/h 0.77),
    # with the A-C column at 0.33 aspect. 280 puts A-C at 338x1023 = 0.33.
    zoom_panel_height = 280

    # Density panel dimensions
    density_gap = 12             # gap between zoom tracks and density panels
    density_panel_width = 80     # width of each WIDE density sub-panel
    sst1_panel_width = 14        # the SST1 strip is narrow, as published
    density_inter_gap = 10       # gap between density sub-panels

    # Three panels, in published order: chromosome-specific k-mers, cumulative
    # density, then the SST1 strip. SST1 used to be a sub-track inside the
    # k-mer panel; it is its own rectangle now.
    total_density_width = (2 * density_panel_width + sst1_panel_width +
                           2 * density_inter_gap)

    # Total figure width
    total_width = (margin_left + track_block_width + density_gap +
                   total_density_width + margin_right)

    # Total figure height (18px extra per sample for title row)
    # Add extra bottom margin on the last panel for rotated track labels
    # Clears whatever hangs below the last panel. The feature-track labels are
    # still rotated ("Chromosome" ~60 px); the density labels are now horizontal
    # two-liners (~24 px) plus the rotated "SST1" (~30 px).
    label_margin = 90
    total_height = (margin_top +
                    len(robs) * (zoom_panel_height + margin_bottom + 18) +
                    (len(robs) - 1) * inter_scaffold_gap +
                    label_margin)

    d = draw.Drawing(total_width, total_height)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=background))

    # --- Process each ROB ---
    results_rows = []
    anchor_rows = []
    current_y = margin_top
    merged_displayed = {fs: set() for fs in featuresets}

    featureset_names = {
        'chromosome': 'Chromosome',
        'acrocentric': 'Acrocentric',
        'region': 'Satellite',
        'repeat': 'Repeat',
        'gene': 'Gene',
    }

    for rob_idx, rob in enumerate(robs):
        print(f"Processing {rob['sample']} — {rob['iscn']}...")

        # Focus region
        focus_start = int(rob["iscn_start"] - args.padding * 1e6)
        focus_end = int(rob["iscn_end"] + args.padding * 1e6)
        focus_start = max(0, focus_start)
        view_span = focus_end - focus_start
        padding_bp = view_span // 2

        # Determine bin size for zoom track rendering
        bin_size = get_bin_size_for_zoom(padding_bp)

        bed_dir = str(PROJECT_ROOT / config["samples"][rob["sample"]]["bed_dir"])
        assembly_name = config["samples"][rob["sample"]]["assembly_name"]

        # Load BED features for zoom panel
        features_by_track = {}
        for fs in featuresets:
            for try_bin in [100, 100000, 1000000]:
                bed_path = os.path.join(
                    bed_dir,
                    f"{assembly_name}.{database}.{fs}.presmoothed.scaffolded.binned{try_bin}.features.bed.gz"
                )
                if os.path.exists(bed_path):
                    break
            raw_features = load_bed_features(bed_path, rob["scaffold"],
                                             focus_start, focus_end)
            if try_bin < bin_size:
                features_by_track[fs] = bin_features(raw_features, focus_start,
                                                     focus_end, bin_size)
            else:
                features_by_track[fs] = raw_features
            # Track displayed features for legend
            for f in features_by_track[fs]:
                merged_displayed[fs].add(f['feature'])

        pixels_per_bp = zoom_panel_height / view_span
        panel_y = current_y + 18  # extra space for title

        # --- Side label (rotated ISCN notation, like stacked plot) ---
        side_label = f"{rob['sample']} — {rob['iscn']}"
        label_x = margin_left - 34
        label_y_mid = panel_y + zoom_panel_height / 2
        d.append(draw.Text(side_label, font_size=11, x=label_x, y=label_y_mid,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle', dominant_baseline='middle',
                           font_weight='bold',
                           transform=f'rotate(-90 {label_x} {label_y_mid})'))

        # --- Zoom window label at top ---
        if view_span >= 1_000_000:
            zoom_label = f"{view_span / 1_000_000:.1f} Mb"
        else:
            zoom_label = f"{view_span / 1_000:.0f} kb"
        d.append(draw.Text(zoom_label, font_size=11,
                           x=margin_left + track_block_width / 2,
                           y=panel_y - 6,
                           fill=text_color, font_family='sans-serif',
                           text_anchor='middle'))

        # --- Draw feature tracks ---
        track_x = margin_left
        for i, fs in enumerate(featuresets):
            x = track_x + i * (bar_width + bar_spacing)
            features = features_by_track.get(fs, [])
            fs_colors = ks_colors.get(fs, {})
            draw_feature_track(d, features, fs_colors, x, panel_y, bar_width,
                               focus_start, pixels_per_bp, zoom_panel_height)

        block_x_start = track_x - 1
        block_x_end = track_x + track_block_width + 1

        density_x_start = track_x + track_block_width + density_gap

        # No width bracket is drawn. It measured conf/robs.yaml's hand-set
        # start/end (the display window), which read as a measurement of the
        # data and is not one. The panel-span scale bar above it is a real
        # scale annotation and stays. Removed 2026-08-28 with the dashed lines.

        # Coordinate labels
        draw_coordinate_labels(d, block_x_start - 2, panel_y,
                               panel_y + zoom_panel_height,
                               focus_start, focus_end, text_color)

        # Track labels below (rotated) — only for last sample
        if rob_idx == len(robs) - 1:
            label_y_bottom = panel_y + zoom_panel_height + 6
            for i, fs in enumerate(featuresets):
                x = track_x + i * (bar_width + bar_spacing) + bar_width / 2
                label = featureset_names.get(fs, fs[:6])
                d.append(draw.Text(label, font_size=11, x=x, y=label_y_bottom,
                                   fill=text_color, font_family='sans-serif',
                                   text_anchor='end', dominant_baseline='middle',
                                   transform=f"rotate(-90 {x} {label_y_bottom})"))

        # --- Density analysis ---
        print(f"  Loading density data...")
        chrom_bed = load_chrom_bed(rob["sample"], args.resolution, config,
                                   rob["scaffold"], focus_start, focus_end)
        result = cumulative_breakpoint_analysis(
            chrom_bed, rob["scaffold"],
            rob["label_left"], rob["label_right"],
            focus_start, focus_end, resolution=args.resolution
        )

        acro_bed = load_acro_bed(rob["sample"], args.resolution, config,
                                 rob["scaffold"], focus_start, focus_end)
        sst1_positions, sst1_cov = get_sst1_density(
            acro_bed, rob["scaffold"], focus_start, focus_end,
            resolution=args.resolution
        )

        # Get chromosome colors from palette
        left_hex = ks_colors.get('chromosome', {}).get(
            rob["label_left"], ('#60A5FA', 1.0))[0]
        right_hex = ks_colors.get('chromosome', {}).get(
            rob["label_right"], ('#F07167', 1.0))[0]

        # Panel 1: chromosome-specific k-mers
        draw_density_panel_coverage(
            d, density_x_start, panel_y,
            density_panel_width, zoom_panel_height,
            result["positions"], result["left_cov"], result["right_cov"],
            focus_start, focus_end, left_hex, right_hex,
            rob["label_left"], rob["label_right"],
            cov_frac=args.cov_frac
        )

        # Panel 2: cumulative density (two curves, one per fusion partner)
        comb_x = density_x_start + density_panel_width + density_inter_gap
        draw_density_panel_combined(
            d, comb_x, panel_y,
            density_panel_width, zoom_panel_height,
            result["positions"], result["combined"], result["max_score"],
            focus_start, focus_end,
            result["plateau_start"], result["plateau_end"],
            result["ci_999"], result["ci_99"],
            rob["sst1_center"],
            left_frac=result["left_frac"], right_frac=result["right_frac"],
            left_color=left_hex, right_color=right_hex
        )

        # Panel 3: SST1 strip
        sst1_x = comb_x + density_panel_width + density_inter_gap
        draw_sst1_panel(
            d, sst1_x, panel_y, sst1_panel_width, zoom_panel_height,
            sst1_positions, sst1_cov, focus_start, focus_end,
            cov_frac=args.cov_frac
        )

        # Anchor points for the assembled figure's leader lines: the y-centre of
        # each contiguous run of drawn SST1 bins, in this canvas's coordinates.
        # Recomputed with the same threshold the panel drew with, so the anchors
        # cannot drift away from the marks.
        ppbp_a = zoom_panel_height / (focus_end - focus_start)
        res_a = (sst1_positions[1] - sst1_positions[0]
                 if len(sst1_positions) > 1 else 100)
        lit = [i for i in range(len(sst1_positions))
               if sst1_cov[i] > args.cov_frac * res_a]
        runs, run = [], []
        for i in lit:
            if run and i != run[-1] + 1:
                runs.append(run); run = []
            run.append(i)
        if run:
            runs.append(run)
        sst1_anchors = []
        for r in runs:
            p0 = float(sst1_positions[r[0]])
            p1 = float(sst1_positions[r[-1]]) + res_a
            sst1_anchors.append({
                "bp_start": p0, "bp_end": p1, "n_bins": len(r),
                "y_center": panel_y + ((p0 + p1) / 2 - focus_start) * ppbp_a,
            })
        anchor_rows.append({
            "sample": rob["sample"],
            "panel_y": panel_y,
            "panel_height": zoom_panel_height,
            "sst1_x": sst1_x,
            "sst1_width": sst1_panel_width,
            "focus_start": focus_start,
            "focus_end": focus_end,
            "sst1_runs": sst1_anchors,
        })

        # Axis labels below the last panel only. The two wide panels get
        # horizontal labels; SST1 stays vertical because its strip is 14 px
        # wide, and it is drawn ONCE here rather than above every panel.
        if rob_idx == len(robs) - 1:
            label_y_bottom = panel_y + zoom_panel_height + 14
            wide_labels = [
                (density_x_start + density_panel_width / 2,
                 ['Chromosome-', 'specific k-mers']),
                (comb_x + density_panel_width / 2, ['Cumulative', 'density']),
            ]
            # The cumulative curves are fractions, so the panel spans 0-1.
            for tick_val, tick_x in ((0, comb_x), (1, comb_x + density_panel_width)):
                d.append(draw.Text(str(tick_val), font_size=9, x=tick_x,
                                   y=panel_y + zoom_panel_height + 2,
                                   fill=text_color, font_family='sans-serif',
                                   text_anchor='middle',
                                   dominant_baseline='hanging'))
            for px, lines in wide_labels:
                for li, line in enumerate(lines):
                    d.append(draw.Text(line, font_size=11, x=px,
                                       y=label_y_bottom + li * 12,
                                       fill=text_color, font_family='sans-serif',
                                       text_anchor='middle',
                                       dominant_baseline='hanging'))
            sx = sst1_x + sst1_panel_width / 2
            d.append(draw.Text('SST1', font_size=11, x=sx, y=label_y_bottom,
                               fill='#00AA00', font_family='sans-serif',
                               text_anchor='end', dominant_baseline='middle',
                               transform=f"rotate(-90 {sx} {label_y_bottom})"))

        # Print results
        print(f"  Plateau: {_format_genomic_pos(result['plateau_start'])} – "
              f"{_format_genomic_pos(result['plateau_end'])} "
              f"({format_width(result['plateau_end'] - result['plateau_start'])})")
        print(f"  99.9% CI: {_format_genomic_pos(result['ci_999'][0])} – "
              f"{_format_genomic_pos(result['ci_999'][1])} "
              f"({format_width(result['ci_999'][1] - result['ci_999'][0])})")
        # The band drawn in column 3 is ci_99, NOT the argmax plateau and NOT
        # ci_999 -- see the '99% CI' legend entry. The caption has called it a
        # "99% plateau width", which conflates two different quantities the
        # analysis reports separately. Record all three so the caption can name
        # the one actually drawn.
        sst1 = rob["sst1_center"]
        bp_called = result["plateau_center"]
        results_rows.append({
            "sample": rob["sample"],
            "iscn": rob["iscn"],
            "scaffold": rob["scaffold"],
            "sst1_center": sst1,
            "breakpoint": bp_called,
            "delta_sst1_bp": bp_called - sst1,
            "plateau_start": result["plateau_start"],
            "plateau_end": result["plateau_end"],
            "plateau_width_bp": result["plateau_end"] - result["plateau_start"],
            "ci99_start": result["ci_99"][0],
            "ci99_end": result["ci_99"][1],
            "ci99_width_bp": result["ci_99"][1] - result["ci_99"][0],
            "ci999_width_bp": result["ci_999"][1] - result["ci_999"][0],
        })
        print(f"  99% CI (THE DRAWN BAND): "
              f"{_format_genomic_pos(result['ci_99'][0])} – "
              f"{_format_genomic_pos(result['ci_99'][1])} "
              f"({format_width(result['ci_99'][1] - result['ci_99'][0])})")
        print(f"  breakpoint {bp_called:,}  delta-SST1 {(bp_called - sst1)/1000:+.1f} kb")

        current_y += zoom_panel_height + margin_bottom + inter_scaffold_gap + 18

    # --- Legend ---
    if not args.no_legend:
        legend_x = (margin_left + track_block_width + density_gap +
                    total_density_width + 15)
        legend_y = margin_top + 20
        draw_legend(d, legend_x, legend_y, ks_colors, featuresets,
                    merged_displayed, text_color)

    if args.emit_anchors:
        import json
        with open(args.emit_anchors, "w") as fh:
            json.dump({"canvas_width": total_width,
                       "canvas_height": total_height,
                       "samples": anchor_rows,
                       # Which features actually appeared, so the assembler can
                       # rebuild the identical legend elsewhere on the canvas.
                       "displayed": {k: sorted(v)
                                     for k, v in merged_displayed.items()}},
                      fh, indent=2)
        print(f"Saved {args.emit_anchors}")

    # Save
    svg_path = output_dir / "breakpoint_density_stacked.svg"
    png_path = output_dir / "breakpoint_density_stacked.png"

    res_path = output_dir / "breakpoint_results.tsv"
    pd.DataFrame(results_rows).to_csv(res_path, sep="\t", index=False)
    print(f"Saved {res_path}")
    d.save_svg(str(svg_path))
    print(f"\nSaved {svg_path}")

    # Try to convert to PNG
    try:
        import subprocess
        subprocess.run(
            ['rsvg-convert', '-d', '300', '-p', '300',
             str(svg_path), '-o', str(png_path)],
            check=True, capture_output=True
        )
        print(f"Saved {png_path}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Note: rsvg-convert not available; PNG not generated")


if __name__ == "__main__":
    main()
