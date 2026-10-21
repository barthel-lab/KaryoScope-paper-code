#!/usr/bin/env python3
"""
Plot CHM13 chromosomes (region featureset) with lollipop markers showing
where SST1-derived chromosome-specific k-mers land.
"""

import gzip
import drawsvg as draw
import numpy as np
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CHM13_REGION_BED = (Path.home() / "Documents/software/KaryoScope/results/CHM13/"
                    "assembly/hap1/KaryoScope/KS_human_CHM13/"
                    "CHM13.hap1.KS_human_CHM13.region.presmoothed.features.bed.gz")

# Chromosome sizes from CHM13 FAI
CHROM_SIZES = {
    "chr7": 160_567_428,
    "chr8": 146_259_331,
    "chr9": 150_617_247,
    "chr17": 84_276_897,
    "chr20": 66_210_255,
}

# K-mer positions on CHM13 (from search results)
KMER_POSITIONS = {
    "chr7": [57_890_497, 57_890_498, 57_890_499, 57_890_500],
    "chr8": [39_858_704],
    "chr9": [
        41_644_530, 41_644_531, 41_646_516,
        41_647_701, 41_647_702, 41_647_703, 41_647_704, 41_647_705, 41_647_706,
        41_648_112, 41_648_113, 41_648_114, 41_648_115, 41_648_116, 41_648_117,
        41_648_118, 41_648_119, 41_648_120, 41_648_121, 41_648_122, 41_648_123,
        41_648_124, 41_648_125, 41_648_126, 41_648_127, 41_648_128, 41_648_129,
        41_648_130, 41_648_131, 41_648_132, 41_648_133, 41_648_134,
        41_650_715, 41_650_716, 41_650_717, 41_650_718, 41_650_719, 41_650_720,
        41_650_721, 41_650_722, 41_650_723, 41_650_724, 41_650_725, 41_650_726,
        41_650_727, 41_650_728, 41_650_729, 41_650_730, 41_650_731,
        43_196_630, 43_196_631, 43_196_632, 43_196_633, 43_196_634, 43_196_635,
        43_196_636, 43_196_637, 43_196_638, 43_196_639, 43_196_640, 43_196_641,
        43_196_642, 43_196_643, 43_196_644, 43_196_645, 43_196_646,
        43_196_693, 43_196_694,
        48_143_239, 48_143_240, 48_143_241, 48_143_242, 48_143_243, 48_143_244,
        48_143_245, 48_143_246, 48_143_247, 48_143_248, 48_143_249, 48_143_250,
        48_143_251, 48_143_252, 48_143_253, 48_143_254, 48_143_255, 48_143_256,
        48_143_257, 48_143_258, 48_143_259, 48_143_260, 48_143_261,
        48_144_058, 48_144_059,
        48_144_462, 48_144_463, 48_144_464, 48_144_465, 48_144_466, 48_144_472,
        48_144_473, 48_144_474, 48_144_475, 48_144_476,
        48_145_666,
        48_359_227,
        48_359_345, 48_359_346, 48_359_347, 48_359_348,
        48_359_449, 48_359_450, 48_359_451, 48_359_452, 48_359_453, 48_359_454,
        48_359_455, 48_359_456, 48_359_457, 48_359_458, 48_359_459,
        48_359_614, 48_359_615, 48_359_616, 48_359_617, 48_359_618, 48_359_619,
        79_502_506, 79_502_507, 79_502_508, 79_502_509, 79_502_510,
        79_503_125, 79_503_126,
        79_505_108,
        79_506_300, 79_506_301, 79_506_302, 79_506_303, 79_506_304, 79_506_305,
        79_506_723, 79_506_724, 79_506_725, 79_506_726, 79_506_727, 79_506_728,
        79_506_729, 79_506_730, 79_506_731, 79_506_732, 79_506_733, 79_506_734,
        79_506_735, 79_506_736, 79_506_737, 79_506_738, 79_506_739, 79_506_740,
        79_506_741, 79_506_742, 79_506_743, 79_506_744, 79_506_745, 79_506_746,
        79_506_747, 79_506_748, 79_506_749, 79_506_750, 79_506_751, 79_506_752,
        79_507_132, 79_507_133, 79_507_134, 79_507_135, 79_507_136, 79_507_137,
        79_509_338, 79_509_339, 79_509_340, 79_509_341, 79_509_342, 79_509_343,
        79_509_344, 79_509_345, 79_509_346, 79_509_347, 79_509_348, 79_509_349,
        79_509_350, 79_509_351, 79_509_352, 79_509_353, 79_509_354,
    ],
    "chr17": [
        21_628_481, 21_628_777, 21_628_778, 21_629_423, 21_631_105,
        22_684_585, 22_684_586, 22_685_041,
        22_686_886, 22_686_887, 22_688_050, 22_688_051, 22_689_213,
        22_689_454, 22_689_455, 22_689_456, 22_689_457, 22_689_458,
        22_689_459, 22_689_460, 22_689_461, 22_689_462, 22_689_463,
        22_689_464, 22_689_465, 22_689_466, 22_689_467, 22_689_468,
        22_689_469, 22_689_470, 22_689_471, 22_689_472,
    ],
    "chr20": [
        25_863_549, 25_863_796, 25_863_797, 25_863_798, 25_863_799,
        25_864_266, 25_865_927,
        25_870_081, 25_870_082, 25_870_083, 25_870_084, 25_870_085,
        25_870_086, 25_870_087, 25_870_088, 25_870_089, 25_870_090,
        25_870_091, 25_870_092, 25_870_093, 25_870_094, 25_870_095,
        26_062_488, 26_062_489, 26_062_490, 26_062_491, 26_062_492,
        26_062_493, 26_062_494,
        26_280_991, 26_280_992, 26_280_993, 26_280_994, 26_280_995,
        26_280_996, 26_280_997, 26_280_998, 26_280_999,
        31_177_493, 31_177_494, 31_177_495, 31_177_496, 31_177_497,
        31_177_498, 31_177_499, 31_177_500, 31_177_501, 31_177_502,
        31_177_503, 31_177_504, 31_177_505,
    ],
}

# Region featureset colors
REGION_COLORS = {
    "active_specific": "#fa0000",
    "inactive_specific": "#ff9200",
    "divergent_specific": "#990000",
    "monomeric_specific": "#ffcc99",
    "hsat1A_specific": "#91ff00",
    "hsat1B_specific": "#00994d",
    "hsat2_specific": "#0081fa",
    "hsat3_specific": "#0000fa",
    "bsat_specific": "#fa99ff",
    "gsat_specific": "#ac33c7",
    "censat_specific": "#00cccc",
    "ct_specific": "#B8AFA2",
    "centromeric_multigroup1": "#B0C4DE",
    "asat_multigroup1": "#fa0000",
    "hor_multigroup1": "#fa0000",
    "hsat_multigroup1": "#B0C4DE",
    "hsat1_multigroup1": "#B0C4DE",
    "arm_multigroup1": "#808080",
    "p_arm_specific": "#808080",
    "q_arm_specific": "#808080",
    "noncentromeric_specific": "#964B00",
}

# Chromosome-specific colors from KaryoScope
CHR_COLORS = {
    "chr7": "#a40000",
    "chr8": "#00d3e0",
    "chr9": "#ff0000",
    "chr17": "#003a6a",
    "chr20": "#a4d1ff",
}

# ---------------------------------------------------------------------------
# Load region features
# ---------------------------------------------------------------------------
_FEATURE_CACHE = {}


def _load_all_features(wanted):
    """One pass over the region BED, bucketing the chromosomes we need.

    The per-chromosome version re-read the whole file for every chromosome --
    five full passes over a multi-GB gzip. The BED is far too large for that to
    be free, and nothing about it changes between calls.
    """
    global _FEATURE_CACHE
    missing = [c for c in wanted if c not in _FEATURE_CACHE]
    if not missing:
        return
    for c in missing:
        _FEATURE_CACHE[c] = []
    want = set(missing)
    opener = gzip.open if str(CHM13_REGION_BED).endswith(".gz") else open
    with opener(CHM13_REGION_BED, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] not in want:
                continue
            _FEATURE_CACHE[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
    for c in missing:
        if not _FEATURE_CACHE[c]:
            raise SystemExit(
                f"no region records for {c} in {CHM13_REGION_BED}. "
                "An empty chromosome draws a blank track rather than failing, "
                "so this is raised deliberately.")


def load_region_features(chrom, bin_size=100_000):
    """Load and bin region features for a chromosome."""
    chrom_size = CHROM_SIZES[chrom]
    _load_all_features([chrom])
    features = _FEATURE_CACHE[chrom]

    # Bin using max-overlap logic
    n_bins = (chrom_size + bin_size - 1) // bin_size
    bins = defaultdict(lambda: defaultdict(int))
    for start, end, feat in features:
        start_bin = start // bin_size
        end_bin = (end - 1) // bin_size
        for b in range(start_bin, min(end_bin + 1, n_bins)):
            bs = b * bin_size
            be = min(bs + bin_size, chrom_size)
            overlap = min(end, be) - max(start, bs)
            if overlap > 0:
                bins[b][feat] += overlap

    binned = []
    for b in range(n_bins):
        bs = b * bin_size
        be = min(bs + bin_size, chrom_size)
        if b in bins:
            best = max(bins[b], key=lambda f: bins[b][f])
        else:
            best = "novel"
        binned.append((bs, be, best))

    # Merge consecutive same-feature bins
    merged = []
    for bs, be, feat in binned:
        if merged and merged[-1][2] == feat:
            merged[-1] = (merged[-1][0], be, feat)
        else:
            merged.append((bs, be, feat))
    return merged


# ---------------------------------------------------------------------------
# Cluster k-mer positions for lollipop markers
# ---------------------------------------------------------------------------
def cluster_positions(positions, merge_dist=500_000):
    """Cluster nearby positions into groups for lollipop display."""
    if not positions:
        return []
    positions = sorted(set(positions))
    clusters = []
    cluster_start = positions[0]
    cluster_end = positions[0]
    cluster_count = 1
    for p in positions[1:]:
        if p - cluster_end <= merge_dist:
            cluster_end = p
            cluster_count += 1
        else:
            clusters.append((cluster_start, cluster_end, cluster_count))
            cluster_start = p
            cluster_end = p
            cluster_count = 1
    clusters.append((cluster_start, cluster_end, cluster_count))
    return clusters


# ---------------------------------------------------------------------------
# Build the panel
# ---------------------------------------------------------------------------
def build(bed_path=None, chroms=None, draw_title=False, track_height=None):
    """Return (group, width, height, anchors) for panel E.

    `chroms` defaults to the five loaded here. Note the PUBLISHED figure shows
    four -- chr8 is dropped, and it is the one with a single k-mer position.
    Pass an explicit list to match it.
    """
    global CHM13_REGION_BED
    if bed_path is not None:
        CHM13_REGION_BED = Path(bed_path)
    if not Path(CHM13_REGION_BED).exists():
        raise SystemExit(
            f"CHM13 region BED not found: {CHM13_REGION_BED}\n"
            "Panel E needs CHM13 annotated against the v1 database "
            "(KS_human_CHM13); the v1 and v2 region vocabularies differ.")
    g = draw.Group()
    _chroms_override = chroms
    _track_height_override = track_height
    _load_all_features(chroms or ["chr7", "chr8", "chr9", "chr17", "chr20"])

    # ---------------------------------------------------------------------------
    # Layout (rotated 90°: chromosomes run vertically, lollipops extend right)
    # ---------------------------------------------------------------------------
    font_family = "sans-serif"
    text_color = "#000000"
    axis_color = "#545454"

    chroms = _chroms_override or ["chr7", "chr8", "chr9", "chr17", "chr20"]
    max_size = max(CHROM_SIZES.values())

    # Proportions follow the published panel: narrow tracks, long lollipop
    # sticks, and the count labels doing the work of filling the gap between
    # chromosomes. Tracks stay thin so they read as long rather than stubby --
    # the published aspect is roughly 24:1 height to width.
    track_width = 12          # horizontal width of each chromosome track
    lollipop_length = 24      # how far lollipop sticks extend to the right
    label_space = 17          # space for count labels (1-2 digits)
    col_width = track_width + lollipop_length + label_space + 5
    margin_left = 15
    margin_right = 8
    # Only the chromosome labels need headroom when the title is suppressed,
    # which is the case in the assembled figure.
    margin_top = 50 if draw_title else 22
    margin_bottom = 16
    # Overridable so the assembler can size this panel to its slot. Shrinking
    # the whole panel with a transform instead would shrink the text with it.
    track_height = _track_height_override or 500

    total_w = margin_left + len(chroms) * col_width + margin_right
    total_h = margin_top + track_height + margin_bottom

    pixels_per_bp = track_height / max_size

    if draw_title:
        g.append(draw.Text(
            "SST1-derived k-mer positions on non-acrocentric CHM13 chromosomes",
            font_size=11, x=total_w / 2, y=15,
            fill=text_color, font_family=font_family,
            text_anchor="middle", font_weight="bold"))

    # Scale bar (vertical)
    scale_mb = 50
    scale_px = scale_mb * 1_000_000 * pixels_per_bp
    sb_x = margin_left - 5
    sb_y = margin_top + track_height - scale_px
    g.append(draw.Line(sb_x, sb_y, sb_x, sb_y + scale_px,
                       stroke=text_color, stroke_width=1))
    g.append(draw.Line(sb_x - 2, sb_y, sb_x + 2, sb_y,
                       stroke=text_color, stroke_width=1))
    g.append(draw.Line(sb_x - 2, sb_y + scale_px, sb_x + 2, sb_y + scale_px,
                       stroke=text_color, stroke_width=1))
    # Set at the chromosome-name size: at 8 px this was the smallest text in the
    # assembled figure by some way, and panel E is placed at scale 1.0 so the
    # authored size is the printed size. The offset tracks the font size so the
    # rotated label keeps clear of the bar as the text grows.
    scale_label_size = 11
    lbl_x = sb_x - scale_label_size / 2 - 1
    lbl_y = sb_y + scale_px / 2
    g.append(draw.Text(f"{scale_mb} Mb", font_size=scale_label_size,
                       x=lbl_x, y=lbl_y,
                       fill=text_color, font_family=font_family,
                       text_anchor="middle", dominant_baseline="middle",
                       transform=f"rotate(-90 {lbl_x} {lbl_y})"))

    print("Loading region features...")
    for ci, chrom in enumerate(chroms):
        x_track = margin_left + ci * col_width
        chrom_size = CHROM_SIZES[chrom]
        chrom_height = chrom_size * pixels_per_bp

        # Chromosome label at top
        g.append(draw.Text(
            chrom, font_size=11,
            x=x_track + track_width / 2, y=margin_top - 8,
            fill=text_color, font_family=font_family,
            text_anchor="middle",
            font_weight="bold",
        ))

        # Load and draw region features (vertical track)
        features = load_region_features(chrom)
        for start, end, feat in features:
            y = margin_top + start * pixels_per_bp
            h = (end - start) * pixels_per_bp
            if h < 0.3:
                h = 0.3
            color = REGION_COLORS.get(feat, "#808080")
            g.append(draw.Rectangle(x_track, y, track_width, h, fill=color))

        # Track outline
        g.append(draw.Rectangle(x_track, margin_top, track_width, chrom_height,
                                fill="none", stroke="#000000", stroke_width=0.5))

        # Lollipop markers (extending to the right)
        positions = KMER_POSITIONS.get(chrom, [])
        clusters = cluster_positions(positions, merge_dist=500_000)
        chr_color = CHR_COLORS.get(chrom, "#000000")

        # Compute lollipop y-positions and spread labels to avoid overlap
        lollipops = []
        for cl_start, cl_end, count in clusters:
            mid = (cl_start + cl_end) / 2
            y = margin_top + mid * pixels_per_bp
            lollipops.append((y, count))

        # Spread label positions to avoid overlap (minimum 12px apart)
        min_gap = 13
        label_ys = [lp[0] for lp in lollipops]
        # Iterative relaxation
        for _ in range(50):
            moved = False
            for i in range(len(label_ys) - 1):
                gap = label_ys[i + 1] - label_ys[i]
                if gap < min_gap:
                    shift = (min_gap - gap) / 2
                    label_ys[i] -= shift
                    label_ys[i + 1] += shift
                    moved = True
            if not moved:
                break

        x_lollipop_base = x_track + track_width
        for li, (orig_y, count) in enumerate(lollipops):
            label_y = label_ys[li]
            r = min(3 + np.log2(count + 1), 6)
            circle_x = x_lollipop_base + lollipop_length - r

            # Stick from track to circle
            g.append(draw.Line(x_lollipop_base, orig_y,
                               circle_x - r, label_y,
                               stroke=chr_color, stroke_width=0.8))
            # Circle
            g.append(draw.Circle(circle_x, label_y, r,
                                 fill=chr_color, fill_opacity=0.85,
                                 stroke="white", stroke_width=0.5))
            # Count label
            g.append(draw.Text(
                str(count), font_size=10,
                x=circle_x + r + 3, y=label_y,
                fill=chr_color, font_family=font_family,
                text_anchor="start", dominant_baseline="middle",
            ))

        print(f"  {chrom}: {len(features)} regions, {len(clusters)} k-mer clusters "
              f"({len(positions)} k-mers)")

    # Save

    anchors = {"chroms": chroms, "col_width": col_width,
               "margin_left": margin_left, "track_height": track_height}
    return g, total_w, total_h, anchors


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--region-bed", required=True,
                    help="CHM13 region BED annotated against KS_human_CHM13 (v1).")
    ap.add_argument("--chroms", default=None,
                    help="Comma-separated chromosomes (default: all five). "
                         "The published panel E uses chr7,chr9,chr17,chr20.")
    ap.add_argument("--track-height", type=float, default=None)
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    chroms = args.chroms.split(",") if args.chroms else None
    g, w, h, anchors = build(args.region_bed, chroms=chroms, draw_title=True,
                             track_height=args.track_height)
    d = draw.Drawing(w, h)
    d.append(draw.Rectangle(0, 0, w, h, fill="white"))
    d.append(g)
    out_dir = Path(args.output_dir) if args.output_dir else (
        PROJECT_ROOT / "results" / "breakpoint_density")
    out_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / "chm13_kmer_lollipop.svg"
    d.save_svg(str(svg_path))
    print(f"Saved {svg_path}  ({w:.0f} x {h:.0f})  chroms={anchors['chroms']}")
    import subprocess
    png_path = out_dir / "chm13_kmer_lollipop.png"
    try:
        subprocess.run(["rsvg-convert", "-d", "300", "-p", "300",
                        str(svg_path), "-o", str(png_path)],
                       check=True, capture_output=True)
        print(f"Saved {png_path}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Note: rsvg-convert not available; PNG not generated")


if __name__ == "__main__":
    main()
