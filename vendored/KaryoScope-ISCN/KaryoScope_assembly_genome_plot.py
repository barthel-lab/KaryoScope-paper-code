#!/usr/bin/env python3
# Vendored from barthel-lab/KaryoScope-ISCN at commit 4e42a98d6526ecb31ca6fb5d6c8910e4e84c605a, file KaryoScope_assembly_genome_plot.py.
# Copied unmodified except: one docstring example path replaced with a placeholder. Authors (git log): Floris Barthel <fbarthel@tgen.org>.
# Part of the KaryoScope-paper-code snapshot; see vendored/README.md.
"""
KaryoScope Assembly Genome Plot - Mirrored Haplotype Layout

Generates a genome-wide visualization showing each chromosome with both
haplotypes in a mirrored TACCAT layout, arranged horizontally:

    chr1                          chr2
    h1          h2                h1          h2
    T  A  C  |  C  A  T          T  A  C  |  C  A  T

Each panel is a single KaryoScope featureset track:
  - Genome (A): chromosome featureset at 1Mb bins
  - Subtelomere (T): subtelomeric featureset, first/last 50kb
  - Centromere (C): region featureset, centromeric region

Usage:
    python KaryoScope_assembly_genome_plot.py \\
        --bed-dir data/ \\
        --colors-dir <colors_dir>/KS_human_CHM13/ \\
        --centromeres HG002.assembly.KS_human_CHM13.region.presmoothed.scaffolded.centromeres.bed \\
        --sample HG002.assembly \\
        --sex XY \\
        --output HG002_genome_plot.svg
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from math import floor

import drawsvg as draw


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="KaryoScope mirrored haplotype genome plot",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed-dir", required=True,
                        help="Directory containing pre-binned BED files")
    parser.add_argument("--colors-dir", required=True,
                        help="Directory with color files")
    parser.add_argument("--centromeres", required=True,
                        help="Centromeres BED file")
    parser.add_argument("--output", "-o", required=True,
                        help="Output SVG file")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--sample", default="HG002.assembly",
                        help="Sample name (default: HG002.assembly)")
    parser.add_argument("--smoothness", default="presmoothed",
                        help="Smoothness level (default: presmoothed)")
    parser.add_argument("--background", default="white",
                        choices=["white", "black"],
                        help="Background color")
    parser.add_argument("--sex", default=None,
                        help="Sex: XX, XY, male, female, reference")
    parser.add_argument("--subtelo-span", type=int, default=25_000,
                        help="Subtelomere zoom span in bp (default: 25000)")
    parser.add_argument("--chromosomes", default=None,
                        help="Comma-separated chromosomes to plot (default: all)")
    parser.add_argument("--scale", type=float, default=1.4,
                        help="Pixels per Mb for genome track (default: 1.4)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# BED loading
# ---------------------------------------------------------------------------
BED_CACHE = {}


def open_bed(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "r")


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
    for start, end, feature in BED_CACHE.get(bed_path, {}).get(scaffold, []):
        if end > view_start and start < view_end:
            features.append({
                "start": max(start, view_start),
                "end": min(end, view_end),
                "feature": feature,
            })
    return features


# ---------------------------------------------------------------------------
# Binning
# ---------------------------------------------------------------------------
def bin_features(features, view_start, view_end, bin_size):
    if not features:
        return []
    bins = defaultdict(lambda: defaultdict(int))
    for feat in features:
        start_bin = feat["start"] // bin_size
        end_bin = (feat["end"] - 1) // bin_size
        for b_idx in range(start_bin, end_bin + 1):
            bs = b_idx * bin_size
            be = bs + bin_size
            overlap = min(feat["end"], be) - max(feat["start"], bs)
            if overlap > 0:
                bins[b_idx][feat["feature"]] += overlap
    binned = []
    for b_idx in sorted(bins.keys()):
        bs = b_idx * bin_size
        be = bs + bin_size
        if be <= view_start or bs >= view_end:
            continue
        best = max(bins[b_idx], key=lambda f: (bins[b_idx][f], f != "novel"))
        if best == "novel" and bins[b_idx][best] <= bin_size / 2:
            non = {f: v for f, v in bins[b_idx].items() if f != "novel"}
            if non:
                best = max(non, key=non.get)
        binned.append({
            "start": max(bs, view_start),
            "end": min(be, view_end),
            "feature": best,
        })
    merged = []
    for feat in binned:
        if merged and merged[-1]["end"] == feat["start"] and merged[-1]["feature"] == feat["feature"]:
            merged[-1]["end"] = feat["end"]
        else:
            merged.append(feat.copy())
    return merged


def get_bin_size_for_zoom(span):
    padding = span // 2
    if padding >= 10_000_000:
        return 100_000
    elif padding >= 5_000_000:
        return 50_000
    elif padding >= 2_000_000:
        return 10_000
    elif padding >= 500_000:
        return 1_000
    else:
        return 500


# ---------------------------------------------------------------------------
# Color loading
# ---------------------------------------------------------------------------
def load_colors(colors_dir, database, featuresets, background="white"):
    featureset_colors = {}
    novel_color = "#000000" if background == "white" else "#ffffff"
    for fs in featuresets:
        colors_path = os.path.join(colors_dir, f"{database}.{fs}.colors.txt")
        featureset_colors[fs] = {
            "novel": (novel_color, 1.0),
            "unknown": ("#808080", 1.0),
        }
        if not os.path.exists(colors_path):
            print(f"Warning: Colors file not found: {colors_path}", file=sys.stderr)
            continue
        with open(colors_path, "r") as f:
            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) >= 2:
                    feature, color = parts[0], parts[1]
                    if i == 0 and feature.lower() == "feature":
                        continue
                    featureset_colors[fs][feature] = (color, 1.0)
                    if feature.endswith("_specific"):
                        featureset_colors[fs][feature[:-9]] = (color, 1.0)
    return featureset_colors


# ---------------------------------------------------------------------------
# Chromosome helpers
# ---------------------------------------------------------------------------
def chromosome_sort_key(chrom):
    core = chrom.replace("chr", "")
    if core.isdigit():
        return (0, int(core))
    order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    return (1, order.get(core, 99), core)


_SEX_KARYOTYPE_MAP = {
    "XY": {"XX": "female", "XY": "male"},
    "X0": {"XX": "female", "X0": "male"},
    "ZW": {"ZZ": "male", "ZW": "female"},
    "ZO": {"ZZ": "male", "ZO": "female"},
}


def normalize_sex(raw_sex, sex_determination_system="XY"):
    if raw_sex is None:
        return None
    if raw_sex.lower() in ("male", "female", "reference"):
        return raw_sex.lower()
    kmap = _SEX_KARYOTYPE_MAP.get(sex_determination_system.upper(), {})
    result = kmap.get(raw_sex.upper())
    if result is None:
        print(f"Warning: unrecognized --sex '{raw_sex}'. Treating as unknown.",
              file=sys.stderr)
    return result


def get_expected_haps(chromosome, sex, haplotypes):
    system = {
        "heterogametic_sex": "male",
        "homogametic_chromosome": "chrX",
        "sex_chromosomes": ["chrY", "chrX"],
    }
    if sex is None:
        if chromosome in system["sex_chromosomes"]:
            return []
        return haplotypes
    if sex == "reference":
        return [haplotypes[0]]
    if chromosome not in system["sex_chromosomes"]:
        return haplotypes
    is_heterogametic = sex == system["heterogametic_sex"]
    if is_heterogametic:
        try:
            idx = system["sex_chromosomes"].index(chromosome)
            return [haplotypes[idx]] if idx < len(haplotypes) else []
        except ValueError:
            return []
    return haplotypes if chromosome == system["homogametic_chromosome"] else []


# ---------------------------------------------------------------------------
# Scaffold / centromere discovery
# ---------------------------------------------------------------------------
def discover_scaffolds(bed_path):
    """Discover scaffolds. For multi-scaffold (chrom, hap) groups, assign
    synthetic hap names: hap1, hap1b, hap1c, ... or hap2, hap2b, ..."""
    load_bed_file_to_cache(bed_path)
    cached = BED_CACHE.get(bed_path, {})
    scaffolds = {}
    # First pass: group scaffolds by (chrom, hap_base)
    grouped = defaultdict(list)  # (chrom, hap) -> [(length, scaf_name), ...]
    for scaf, entries in cached.items():
        if not entries:
            continue
        length = max(end for _, end, _ in entries)
        scaffolds[scaf] = length
        parts = scaf.split("_")
        if len(parts) >= 2:
            chrom = parts[0]
            hap = parts[1]
            grouped[(chrom, hap)].append((length, scaf))

    # Second pass: assign synthetic hap names; sort by length desc so the
    # longest gets 'a' (hap1a / hap2a). Single-scaffold groups use the
    # base name (hap1 / hap2) with no suffix.
    chrom_map = defaultdict(dict)
    for (chrom, hap_base), scaf_list in grouped.items():
        scaf_list.sort(key=lambda x: -x[0])
        if len(scaf_list) == 1:
            chrom_map[chrom][hap_base] = scaf_list[0][1]
        else:
            for i, (_, scaf_name) in enumerate(scaf_list):
                # hap1a, hap1b, hap1c, ...; hap2a, hap2b, hap2c, ...
                syn_hap = f"{hap_base}{chr(ord('a') + i)}"
                chrom_map[chrom][syn_hap] = scaf_name
    return scaffolds, dict(chrom_map)


def load_centromeres(centromere_file):
    coords = {}
    with open_bed(centromere_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            scaf, start, end = parts[0], int(parts[1]), int(parts[2])
            if scaf in coords:
                coords[scaf] = (min(coords[scaf][0], start),
                                max(coords[scaf][1], end))
            else:
                coords[scaf] = (start, end)
    return coords


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def draw_feature_track(d, features, colors, x, y_start, track_width,
                       view_start, pixels_per_bp):
    for feat in features:
        y = y_start + (feat["start"] - view_start) * pixels_per_bp
        height = (feat["end"] - feat["start"]) * pixels_per_bp
        color, opacity = colors.get(feat["feature"], ("#808080", 1.0))
        if height < 0.5:
            height = 0.5
        d.append(draw.Rectangle(x, y, track_width, height,
                                fill=color, fill_opacity=opacity))


def draw_zoom_fan(d, src_y_top, src_y_bot, tgt_y_top, tgt_y_bot,
                  x_src_edge, x_tgt_edge, direction="right", rect_width=6):
    if direction == "right":
        fan_points = [
            (x_src_edge, src_y_top),
            (x_src_edge + rect_width, src_y_top),
            (x_tgt_edge - rect_width, tgt_y_top),
            (x_tgt_edge, tgt_y_top),
            (x_tgt_edge, tgt_y_bot),
            (x_tgt_edge - rect_width, tgt_y_bot),
            (x_src_edge + rect_width, src_y_bot),
            (x_src_edge, src_y_bot),
        ]
    else:
        fan_points = [
            (x_src_edge, src_y_top),
            (x_src_edge - rect_width, src_y_top),
            (x_tgt_edge + rect_width, tgt_y_top),
            (x_tgt_edge, tgt_y_top),
            (x_tgt_edge, tgt_y_bot),
            (x_tgt_edge + rect_width, tgt_y_bot),
            (x_src_edge - rect_width, src_y_bot),
            (x_src_edge, src_y_bot),
        ]
    path_d = f"M {fan_points[0][0]},{fan_points[0][1]} "
    for px, py in fan_points[1:]:
        path_d += f"L {px},{py} "
    path_d += "Z"
    d.append(draw.Path(d=path_d, fill="#B0B0B0", fill_opacity=0.35,
                       stroke="#C0C0C0", stroke_width=0.5))


def draw_track_outline(d, x, y, width, height, stroke_color):
    d.append(draw.Rectangle(x, y, width, height,
                            fill="none", stroke=stroke_color,
                            stroke_width=0.5, stroke_opacity=0.6))


def format_size_label(bp):
    """Format a bp length as a compact size string."""
    if bp >= 1_000_000:
        val = bp / 1_000_000
        return f"{val:.1f} Mb" if val != int(val) else f"{int(val)} Mb"
    elif bp >= 1_000:
        val = bp / 1_000
        return f"{val:.1f} kb" if val != int(val) else f"{int(val)} kb"
    return f"{bp} bp"


# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------
TRACK_W = 10         # subtelomere/centromere tracks (wider)
GENOME_TRACK_W = 6   # genome track
FAN_GAP = 12         # gap for zoom fan polygons
HAP_GAP = 8          # gap between h1 and h2 within a chromosome
CHROM_GAP = 47       # gap between chromosomes
FAN_RECT_W = 3       # rectangular extension on fan ends
SUBTELO_OFFSET_GAP = 3  # extra gap when p-arm telomere is offset outward

# Header geometry (matching KaryoScope_assembly_CLI.py)
CHROM_LABEL_Y = 23
CHROM_LINE_Y = 28
HAP_LABEL_Y = 43
HAP_LINE_Y = 48
INITIAL_Y = 60       # Y where tracks start

# Subtelo/centromere panel sizing: fraction of genome track height
SUBTELO_FRAC = 0.214  # each subtelo panel covers ~21% of genome height
SUBTELO_GAP_FRAC = 0.04  # gap between p and q subtelo panels
CEN_FRAC = 0.70      # centromere panel covers 70% of genome height

MARGIN_LEFT = 50
MARGIN_RIGHT = 25
MARGIN_BOTTOM = 55


def hap_width():
    """Width of one haplotype: T + fan + A + fan + C."""
    return TRACK_W + FAN_GAP + GENOME_TRACK_W + FAN_GAP + TRACK_W


def chrom_width():
    """Width of one chromosome: h1 + gap + h2."""
    return 2 * hap_width() + HAP_GAP


# ---------------------------------------------------------------------------
# Feature loading per scaffold
# ---------------------------------------------------------------------------
def load_scaffold_features(hap_info, bed_dir, sample, database, smoothness,
                           subtelo_span):
    scaf = hap_info["scaffold"]
    scaf_len = hap_info["scaffold_length"]

    genome_bed = os.path.join(
        bed_dir,
        f"{sample}.{database}.chromosome.{smoothness}.scaffolded.binned1000000.features.bed"
    )
    if not os.path.exists(genome_bed):
        genome_bed += ".gz"
    genome_features = load_bed_features(genome_bed, scaf, 0, scaf_len)

    subtelo_bed = os.path.join(
        bed_dir,
        f"{sample}.{database}.subtelomeric.{smoothness}.scaffolded.binned100.features.bed"
    )
    if not os.path.exists(subtelo_bed):
        subtelo_bed += ".gz"

    # Subtelomere: use raw 100bp bins (no rebinning)
    p_arm_features = load_bed_features(subtelo_bed, scaf, 0, subtelo_span)

    q_start = max(0, scaf_len - subtelo_span)
    q_arm_features = load_bed_features(subtelo_bed, scaf, q_start, scaf_len)

    # Centromere: use raw 100kb bins (no rebinning)
    cen_start = hap_info["cen_start"]
    cen_end = hap_info["cen_end"]
    cen_features = []
    if cen_start is not None and cen_end is not None:
        cen_bed = os.path.join(
            bed_dir,
            f"{sample}.{database}.region.{smoothness}.scaffolded.binned100000.features.bed"
        )
        if not os.path.exists(cen_bed):
            cen_bed += ".gz"
        cen_features = load_bed_features(cen_bed, scaf, cen_start, cen_end)

    return {
        "genome": genome_features,
        "p_arm": p_arm_features,
        "q_arm": q_arm_features,
        "centromere": cen_features,
    }


# ---------------------------------------------------------------------------
# Scale bar
# ---------------------------------------------------------------------------
def draw_vertical_scale_bar(d, x, y, length_bp, pixels_per_bp, text_color):
    bar_h = length_bp * pixels_per_bp
    bar_w = 2
    d.append(draw.Rectangle(x, y, bar_w, bar_h, fill=text_color))
    tick_w = 5
    cx = x + bar_w / 2
    d.append(draw.Line(cx - tick_w / 2, y, cx + tick_w / 2, y,
                       stroke=text_color, stroke_width=0.8))
    d.append(draw.Line(cx - tick_w / 2, y + bar_h, cx + tick_w / 2, y + bar_h,
                       stroke=text_color, stroke_width=0.8))
    if length_bp >= 1_000_000:
        label = f"{length_bp / 1_000_000:.0f} Mb"
    elif length_bp >= 1_000:
        label = f"{length_bp / 1_000:.0f} kb"
    else:
        label = f"{length_bp} bp"
    ly = y + bar_h / 2
    lx = cx - 5
    d.append(draw.Text(label, font_size=8, x=lx, y=ly,
                       fill=text_color, font_family="sans-serif",
                       text_anchor="middle", dominant_baseline="middle",
                       transform=f"rotate(-90 {lx} {ly})"))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    pixels_per_bp_genome = args.scale / 1_000_000

    featuresets = ["chromosome", "subtelomeric", "region"]
    print("Loading colors...")
    colors = load_colors(args.colors_dir, args.database, featuresets,
                         background=args.background)

    print("Loading centromeres...")
    centromere_coords = load_centromeres(args.centromeres)

    genome_bed = os.path.join(
        args.bed_dir,
        f"{args.sample}.{args.database}.chromosome.{args.smoothness}"
        f".scaffolded.binned1000000.features.bed"
    )
    if not os.path.exists(genome_bed):
        genome_bed += ".gz"
    print(f"Discovering scaffolds...")
    scaffold_lengths, chrom_map = discover_scaffolds(genome_bed)
    print(f"  Found {len(scaffold_lengths)} scaffolds across "
          f"{len(chrom_map)} chromosomes")

    all_haps = set()
    for haps in chrom_map.values():
        all_haps.update(haps.keys())
    haplotypes = sorted(all_haps) or ["hap1", "hap2"]

    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")]
    else:
        chromosomes = sorted(chrom_map.keys(), key=chromosome_sort_key)

    sex = normalize_sex(args.sex)

    # -----------------------------------------------------------------------
    # Build per-chromosome layout (horizontal arrangement)
    # -----------------------------------------------------------------------
    # Determine which haplotypes to draw per chromosome
    chrom_draw_haps = {}
    for chrom in chromosomes:
        hap_data = chrom_map.get(chrom, {})
        expected = get_expected_haps(chrom, sex, haplotypes)
        actual = [h for h in haplotypes if h in hap_data]
        # Always draw actual; only fall back to expected for sex chromosomes
        # when no data exists at all.
        if actual:
            draw_haps = actual
        elif chrom in ("chrX", "chrY"):
            draw_haps = [h for h in haplotypes if h in expected]
        else:
            draw_haps = []
        if not draw_haps and not actual:
            continue
        chrom_draw_haps[chrom] = draw_haps

    # Global max length -> image height
    max_genome_len = 0
    for chrom, draw_haps in chrom_draw_haps.items():
        hap_data = chrom_map.get(chrom, {})
        for hap in draw_haps:
            scaf = hap_data.get(hap)
            if scaf and scaf in scaffold_lengths:
                max_genome_len = max(max_genome_len, scaffold_lengths[scaf])

    genome_height = max_genome_len * pixels_per_bp_genome
    subtelo_h = genome_height * SUBTELO_FRAC

    # Find max centromere span to set centromere scale
    # The largest centromere gets CEN_FRAC * genome_height; others scale proportionally
    max_cen_span = 0
    for scaf, (cs, ce) in centromere_coords.items():
        max_cen_span = max(max_cen_span, ce - cs)
    max_cen_panel_h = genome_height * CEN_FRAC
    cen_ppb_fixed = max_cen_panel_h / max_cen_span if max_cen_span > 0 else 0

    print(f"  Max genome: {max_genome_len / 1e6:.1f} Mb, "
          f"genome track height: {genome_height:.0f} px")

    # For male XY assemblies, flip chrX/chrY so centromeres face each other
    def _cen_on_right(chrom, hap):
        """Return True if centromere should be on the right (T-A-C layout)."""
        if sex == "male" and chrom == "chrX" and hap.startswith("hap2"):
            return True
        if sex == "male" and chrom == "chrY" and hap.startswith("hap1"):
            return False
        return hap.startswith("hap1")

    # Pre-compute which haplotypes need p-arm offset (telomere overlap)
    offset_extra = TRACK_W + SUBTELO_OFFSET_GAP
    hap_has_offset = {}  # (chrom, hap) -> bool
    for chrom, draw_haps in chrom_draw_haps.items():
        hap_data = chrom_map.get(chrom, {})
        for hap in draw_haps:
            scaf = hap_data.get(hap)
            if not scaf or scaf not in scaffold_lengths:
                hap_has_offset[(chrom, hap)] = False
                continue
            scaf_len = scaffold_lengths[scaf]
            scaf_h = scaf_len * pixels_per_bp_genome
            # Check centromere height for this scaffold
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = (cen_e - cen_s) * cen_ppb_fixed if cen_s is not None else 0
            eff_h = max(scaf_h, this_cen_h)
            # Overlap if p + q subtelo panels don't fit
            hap_has_offset[(chrom, hap)] = (subtelo_h * 2 > eff_h)

    # X positions per chromosome (accounting for p-arm offset)
    current_x = MARGIN_LEFT
    chrom_x_map = {}  # chrom -> x_start
    chrom_hap_x_map = {}  # chrom -> {hap -> x_start_of_hap}
    chrom_hap_label_info = {}  # (chrom, hap) -> (label_x, label_w)

    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        draw_haps = chrom_draw_haps[chrom]
        chrom_x_map[chrom] = current_x

        hap_x = {}
        hx = current_x
        for i, hap in enumerate(draw_haps):
            is_h1 = _cen_on_right(chrom, hap)
            has_off = hap_has_offset.get((chrom, hap), False)

            if has_off and is_h1:
                # h1 offset: p-arm shifts left, so hap block starts earlier
                hap_x[hap] = hx
                label_x = hx - offset_extra
                label_w = hap_width() + offset_extra
                hx += hap_width()
            elif has_off and not is_h1:
                # h2 offset: p-arm shifts right, hap block extends further right
                hap_x[hap] = hx
                label_x = hx
                label_w = hap_width() + offset_extra
                hx += hap_width() + offset_extra
            else:
                hap_x[hap] = hx
                label_x = hx
                label_w = hap_width()
                hx += hap_width()

            chrom_hap_label_info[(chrom, hap)] = (label_x, label_w)

            if i < len(draw_haps) - 1:
                hx += HAP_GAP
        chrom_hap_x_map[chrom] = hap_x
        current_x = hx + CHROM_GAP

    image_width = current_x - CHROM_GAP + MARGIN_RIGHT
    image_height = INITIAL_Y + genome_height + MARGIN_BOTTOM

    background = args.background
    text_color = "black" if background == "white" else "white"

    d = draw.Drawing(image_width, image_height, id_prefix="kg")
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background))

    # -----------------------------------------------------------------------
    # Draw header labels (chromosome + haplotype brackets)
    # -----------------------------------------------------------------------
    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        draw_haps = chrom_draw_haps[chrom]

        # Chromosome line spans from leftmost to rightmost label edge
        first_lx, first_lw = chrom_hap_label_info[(chrom, draw_haps[0])]
        last_lx, last_lw = chrom_hap_label_info[(chrom, draw_haps[-1])]
        chrom_line_x = first_lx
        chrom_line_w = (last_lx + last_lw) - first_lx

        d.append(draw.Rectangle(chrom_line_x, CHROM_LINE_Y, chrom_line_w, 1,
                                fill=text_color))
        d.append(draw.Text(chrom, font_size=14,
                           x=chrom_line_x + chrom_line_w / 2,
                           y=CHROM_LABEL_Y, fill=text_color,
                           font_family="sans-serif", text_anchor="middle"))

        # Haplotype lines + labels (using per-hap label bounds)
        if len(haplotypes) > 1:
            for hap in draw_haps:
                lx, lw = chrom_hap_label_info[(chrom, hap)]
                d.append(draw.Rectangle(lx, HAP_LINE_Y, lw, 1, fill=text_color))
                hap_text = f"h{hap[3:]}" if hap.startswith("hap") else hap
                d.append(draw.Text(hap_text, font_size=14,
                                   x=lx + lw / 2, y=HAP_LABEL_Y,
                                   fill=text_color, font_family="sans-serif",
                                   text_anchor="middle"))

    # -----------------------------------------------------------------------
    # Load all features
    # -----------------------------------------------------------------------
    print("Loading features...")
    all_chrom_features = {}  # chrom -> {hap -> features_dict}

    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        hap_data = chrom_map.get(chrom, {})
        chrom_feats = {}
        for hap in chrom_draw_haps[chrom]:
            scaf = hap_data.get(hap)
            if not scaf or scaf not in scaffold_lengths:
                continue
            scaf_len = scaffold_lengths[scaf]
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            hap_info = {
                "scaffold": scaf,
                "scaffold_length": scaf_len,
                "cen_start": cen_s,
                "cen_end": cen_e,
            }
            chrom_feats[hap] = load_scaffold_features(
                hap_info, args.bed_dir, args.sample, args.database,
                args.smoothness, args.subtelo_span,
            )
        all_chrom_features[chrom] = chrom_feats

    # -----------------------------------------------------------------------
    # Compute per-scaffold subtelo panel positions
    # q-arm tethered to scaffold end; if p/q overlap, offset p outward
    def compute_positions(scaf, scaf_len, is_h1, hx):
        """Compute all Y/X positions for one haplotype.

        Centromere panel height is proportional to centromere span (to scale).
        When the centromere panel is taller than the genome track, the genome
        is vertically centered and telomeres anchor to the centromere panel.

        Returns dict with all positioning info including per-scaffold cen_panel_h.
        """
        scaf_h = scaf_len * pixels_per_bp_genome

        if is_h1:
            base_subtelo_x = hx
            genome_x = hx + TRACK_W + FAN_GAP
            cen_x = genome_x + GENOME_TRACK_W + FAN_GAP
        else:
            cen_x = hx
            genome_x = hx + TRACK_W + FAN_GAP
            base_subtelo_x = genome_x + GENOME_TRACK_W + FAN_GAP

        # Per-scaffold centromere panel height (proportional to span)
        cen_s, cen_e = centromere_coords.get(scaf, (None, None))
        if cen_s is not None:
            this_cen_span = cen_e - cen_s
            this_cen_h = this_cen_span * cen_ppb_fixed
        else:
            this_cen_h = 0

        # If centromere panel is taller than genome, center genome vertically
        if this_cen_h > scaf_h:
            genome_y_offset = (this_cen_h - scaf_h) / 2
        else:
            genome_y_offset = 0
        genome_y0 = y0 + genome_y_offset

        # Centromere panel: centered on centromere position on genome track
        cen_panel_y = y0  # default
        if cen_s is not None:
            cen_mid_y = genome_y0 + ((cen_s + cen_e) / 2) * pixels_per_bp_genome
            cen_panel_y = cen_mid_y - this_cen_h / 2
            cen_panel_y = max(y0, min(y0 + genome_height - this_cen_h,
                                      cen_panel_y))

        # Effective vertical span: the taller of genome track or centromere panel
        if this_cen_h > scaf_h and cen_s is not None:
            eff_top = cen_panel_y
            eff_bot = cen_panel_y + this_cen_h
        else:
            eff_top = genome_y0
            eff_bot = genome_y0 + scaf_h

        p_panel_y = eff_top
        q_panel_y = eff_bot - subtelo_h
        q_panel_y = max(eff_top, q_panel_y)

        # Check overlap: if p bottom > q top, offset p further outward
        p_subtelo_x = base_subtelo_x
        if p_panel_y + subtelo_h > q_panel_y:
            if is_h1:
                p_subtelo_x = base_subtelo_x - TRACK_W - SUBTELO_OFFSET_GAP
            else:
                p_subtelo_x = base_subtelo_x + TRACK_W + SUBTELO_OFFSET_GAP

        return {
            "genome_y0": genome_y0,
            "p_panel_y": p_panel_y,
            "q_panel_y": q_panel_y,
            "p_stx": p_subtelo_x,
            "q_stx": base_subtelo_x,
            "genome_x": genome_x,
            "cen_x": cen_x,
            "cen_panel_y": cen_panel_y,
            "cen_panel_h": this_cen_h,
        }

    # -----------------------------------------------------------------------
    # Pass 1: Draw fans (background layer)
    # -----------------------------------------------------------------------
    print("Drawing fans...")
    y0 = INITIAL_Y

    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        hap_data = chrom_map.get(chrom, {})

        for hap in chrom_draw_haps[chrom]:
            scaf = hap_data.get(hap)
            if not scaf or scaf not in scaffold_lengths:
                continue
            scaf_len = scaffold_lengths[scaf]
            scaf_h = scaf_len * pixels_per_bp_genome

            hx = chrom_hap_x_map[chrom][hap]
            is_h1 = _cen_on_right(chrom, hap)
            pos = compute_positions(scaf, scaf_len, is_h1, hx)
            gy0 = pos["genome_y0"]

            # --- Subtelomere p-arm fan ---
            genome_p_src_top = gy0
            genome_p_src_bot = gy0 + args.subtelo_span * pixels_per_bp_genome
            if genome_p_src_bot - genome_p_src_top < 2:
                genome_p_src_bot = genome_p_src_top + 2

            if is_h1:
                draw_zoom_fan(d, genome_p_src_top, genome_p_src_bot,
                              pos["p_panel_y"], pos["p_panel_y"] + subtelo_h,
                              pos["genome_x"], pos["p_stx"] + TRACK_W,
                              direction="left")
            else:
                draw_zoom_fan(d, genome_p_src_top, genome_p_src_bot,
                              pos["p_panel_y"], pos["p_panel_y"] + subtelo_h,
                              pos["genome_x"] + GENOME_TRACK_W, pos["p_stx"],
                              direction="right")

            # --- Subtelomere q-arm fan ---
            q_start = max(0, scaf_len - args.subtelo_span)
            genome_q_src_top = gy0 + q_start * pixels_per_bp_genome
            genome_q_src_bot = gy0 + scaf_h
            if genome_q_src_bot - genome_q_src_top < 2:
                genome_q_src_top = genome_q_src_bot - 2

            if is_h1:
                draw_zoom_fan(d, genome_q_src_top, genome_q_src_bot,
                              pos["q_panel_y"], pos["q_panel_y"] + subtelo_h,
                              pos["genome_x"], pos["q_stx"] + TRACK_W,
                              direction="left")
            else:
                draw_zoom_fan(d, genome_q_src_top, genome_q_src_bot,
                              pos["q_panel_y"], pos["q_panel_y"] + subtelo_h,
                              pos["genome_x"] + GENOME_TRACK_W, pos["q_stx"],
                              direction="right")

            # --- Centromere fan ---
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = pos["cen_panel_h"]
            if cen_s is not None and this_cen_h > 0:
                genome_cen_src_top = gy0 + cen_s * pixels_per_bp_genome
                genome_cen_src_bot = min(gy0 + cen_e * pixels_per_bp_genome,
                                         gy0 + scaf_h)
                cpy = pos["cen_panel_y"]

                if is_h1:
                    draw_zoom_fan(d, genome_cen_src_top, genome_cen_src_bot,
                                  cpy, cpy + this_cen_h,
                                  pos["genome_x"] + GENOME_TRACK_W, pos["cen_x"],
                                  direction="right")
                else:
                    draw_zoom_fan(d, genome_cen_src_top, genome_cen_src_bot,
                                  cpy, cpy + this_cen_h,
                                  pos["genome_x"], pos["cen_x"] + TRACK_W,
                                  direction="left")

    # -----------------------------------------------------------------------
    # Pass 2: Draw tracks (foreground layer)
    # -----------------------------------------------------------------------
    print("Drawing tracks...")

    for chrom in chromosomes:
        if chrom not in chrom_draw_haps:
            continue
        hap_data = chrom_map.get(chrom, {})

        for hap in chrom_draw_haps[chrom]:
            scaf = hap_data.get(hap)
            if not scaf or scaf not in scaffold_lengths:
                continue
            feats = all_chrom_features.get(chrom, {}).get(hap)
            if feats is None:
                continue

            scaf_len = scaffold_lengths[scaf]
            scaf_h = scaf_len * pixels_per_bp_genome

            hx = chrom_hap_x_map[chrom][hap]
            is_h1 = _cen_on_right(chrom, hap)
            pos = compute_positions(scaf, scaf_len, is_h1, hx)
            gy0 = pos["genome_y0"]

            # Genome track (vertically centered when centromere is taller)
            draw_feature_track(d, feats["genome"], colors["chromosome"],
                               pos["genome_x"], gy0, GENOME_TRACK_W,
                               0, pixels_per_bp_genome)
            draw_track_outline(d, pos["genome_x"], gy0, GENOME_TRACK_W,
                               scaf_h, text_color)

            # Subtelomere p-arm
            p_ppb = subtelo_h / args.subtelo_span
            draw_feature_track(d, feats["p_arm"], colors["subtelomeric"],
                               pos["p_stx"], pos["p_panel_y"], TRACK_W,
                               0, p_ppb)
            draw_track_outline(d, pos["p_stx"], pos["p_panel_y"], TRACK_W,
                               subtelo_h, text_color)

            # Subtelomere q-arm
            q_start = max(0, scaf_len - args.subtelo_span)
            q_ppb = subtelo_h / max(1, scaf_len - q_start)
            draw_feature_track(d, feats["q_arm"], colors["subtelomeric"],
                               pos["q_stx"], pos["q_panel_y"], TRACK_W,
                               q_start, q_ppb)
            draw_track_outline(d, pos["q_stx"], pos["q_panel_y"], TRACK_W,
                               subtelo_h, text_color)

            # Centromere
            cen_s, cen_e = centromere_coords.get(scaf, (None, None))
            this_cen_h = pos["cen_panel_h"]
            if cen_s is not None and feats["centromere"] and this_cen_h > 0:
                cen_span = cen_e - cen_s
                cpy = pos["cen_panel_y"]
                cen_ppb = this_cen_h / cen_span

                draw_feature_track(d, feats["centromere"], colors["region"],
                                   pos["cen_x"], cpy, TRACK_W,
                                   cen_s, cen_ppb)
                draw_track_outline(d, pos["cen_x"], cpy, TRACK_W,
                                   this_cen_h, text_color)

            # --- Size labels anchored to bottom of each track, reading bottom-to-top ---
            label_font = 10
            label_pad = 4  # px below the track bottom

            # Genome size label
            genome_label = format_size_label(scaf_len)
            glx = pos["genome_x"] + GENOME_TRACK_W / 2
            gly = gy0 + scaf_h + label_pad
            d.append(draw.Text(genome_label, font_size=label_font,
                               x=glx, y=gly,
                               fill=text_color, font_family="sans-serif",
                               text_anchor="end", dominant_baseline="middle",
                               transform=f"rotate(-90 {glx} {gly})"))

            # Centromere size label
            if cen_s is not None and this_cen_h > 0:
                cen_span = cen_e - cen_s
                cen_label = format_size_label(cen_span)
                clx = pos["cen_x"] + TRACK_W / 2
                cly = pos["cen_panel_y"] + this_cen_h + label_pad
                d.append(draw.Text(cen_label, font_size=label_font,
                                   x=clx, y=cly,
                                   fill=text_color, font_family="sans-serif",
                                   text_anchor="end", dominant_baseline="middle",
                                   transform=f"rotate(-90 {clx} {cly})"))

        print(f"  {chrom}")

    # -----------------------------------------------------------------------
    # Scale bar (left margin)
    # -----------------------------------------------------------------------
    # Scale bar removed per user request

    # Save
    d.save_svg(args.output)
    print(f"\nSaved: {args.output}")
    print(f"Image size: {image_width:.0f} x {image_height:.0f} px")


if __name__ == "__main__":
    main()
